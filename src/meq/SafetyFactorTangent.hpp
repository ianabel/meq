#ifndef MEQ_SAFETY_FACTOR_TANGENT_HPP
#define MEQ_SAFETY_FACTOR_TANGENT_HPP

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/ContinuousField.hpp"
#include "meq/FluxExtraction.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/Source.hpp"

/**
 * @file SafetyFactorTangent.hpp
 * THE MFEM HALF OF THE q-DRIVEN SYSTEM: an equilibrium solved at a spline of
 * `g^2`, its family's `V' < R^-2 >` at the knots, and their derivatives with
 * respect to the knot values ALONG THE EXACT TANGENT.
 *
 * meq::SafetyFactorSystem is the rows and meq::solveLeastSquares the
 * iteration, both MFEM-free; this is what they call.
 *
 *
 * 1. THE JACOBIAN IS THE LINEARISED SOLVE'S, NOT A DIFFERENCE OF SOLVES
 *
 * Knot `j`'s column: `dR/du_j` is the source's derivative in `g g'` along the
 * knot's cardinal function (NormalisedMHDSource::ggPrimeDerivative()), and
 * GradShafranovSolver::linearisedTangent() returns `[ dx_j ; dpsi_ax,j ]` at
 * the retained, bordered linearisation -- one back-solve -- with `dpsi_bnd,j`,
 * the X-point and the exterior modes beside it where the edge is free. Then
 * `psi*`, `q*` and `psihat` of the tangent (reconstructFields(), linear), and
 * their continuous pair `( dpsi~, dq~ )` by the solver's kept interpolant
 * (GradShafranovSolver::continuousInterpolant()).
 *
 * THE EXTRACTION IS DIFFERENTIATED, NOT DIFFERENCED (SafetyFactorShape.cpp,
 * MEASUREMENTS.md M-408). Its quadrature is `( 1/N ) sum_j a_j/( | u_j . Q |
 * R_j^2 )` over the rays about the located axis, `a_j` each ray's root of
 * `psi~ = c_i`. Each node's change along a tangent is the root's motion,
 * which comes from the label `c_i` (psi_ax and psi_bnd move), from `dpsi~` and
 * from the origin. Then the metric's change follows from `dq~` and `q~`'s
 * gradient at the moved node. The ray origin moves by `-( dq_h/dx )^-1 dq` in
 * the axis's own element, the derivative of the root the map takes for its
 * axis. Holding it fixed misses 1e-5 to 4e-5 of the Jacobian (M-391).
 *
 * Each interior node's values are gathered from its element's polynomial by
 * meq::LinearBlocks, on the device under one. A band node's are sampled from
 * the tracer and differenced. It is the derivative of the map as
 * implemented, so it matches central differences of the extraction to their
 * truncation: 2e-10 to 4e-9 of the Jacobian at the best step. It costs one
 * extraction and m gathers, in place of 2m extractions. In single there is
 * no differencing noise to amplify (M-406). Options::differenced keeps the
 * differenced Jacobian, which is what the shape derivative is held to.
 *
 * The pins: R0B0's row reads the spline alone; the flux and the current are
 * linearised as the quadratures they are read with, through the solver's own
 * assemblers -- assembleSourceDerivative() and the bordered current solve's
 * row and corner (plasmaCurrentLinearisation()) -- contracted in double.
 *
 * 2. EVERY SOLVE IS WARM-STARTED, EXPLICITLY
 *
 * GradShafranovSolver::solve() starts from the stored initial guess, not
 * from the last answer: on `examples/q-driven.toml` every solve of the loop
 * took 4 Newton iterations, the differencing columns 1e-4 away included. So
 * each solve here is handed the last admitted state by setInitialState()
 * (MEASUREMENTS.md M-395).
 */

namespace meq
{

	/// The mesh's area in cells: a quadrilateral one, a triangle one half --
	/// toroidalKnotCount()'s argument on a fixed boundary, where the plasma
	/// is the whole of Omega_h.
	double plasmaCells( mfem::Mesh const &mesh );

	/// The plasma's area in cells: the support's where the solver's source is
	/// confined to one it has decided (GradShafranovSolver::elementInPlasma()),
	/// else, on a free edge, the elements whose centres the solved state puts
	/// at Psi > 0 -- a limited or diverted plasma inside a region it does not
	/// fill -- else the mesh's. After a solve.
	double plasmaCells( GradShafranovSolver const &solver );

	/**
	 * THE DELIVERED TOROIDAL FLUX Phi = sign int_plasma g( Psi )/R dA of the
	 * solver's state, g^2 the spline @a field and Psi = ( psi - psi_b )/span
	 * the total flux's label: on the source term's own points, support and
	 * conductor shift (assembleSourceDerivative()), so it is the quadrature
	 * the flux pin's row linearises, summed in double. What [source]
	 * ToroidalFlux pins and the .nc reports.
	 */
	double deliveredToroidalFlux( GradShafranovSolver &solver, ToroidalFieldSpline const &field,
	                              double sign );

	class SafetyFactorEquilibrium
	{
		public:
			struct Options
			{
				/// The knots: toroidalKnotCut()'s surfaces and cut, spaced in rho.
				FluxFamilyOptions family;

				ScalePin pin = ScalePin::EdgeField;

				/// The target q's sign; the flux is read signed as it.
				double sign = 1.0;

				/// The band, on a curved boundary: the extraction continues
				/// the field by the transfer lift through it, as the write
				/// phase does. Both null on a fitted one.
				mfem::Array<int> const *gammaHMarker = nullptr;
				mfem::TransferPath const *path = nullptr;

				/**
				 * The differenced extraction's step (@ref differenced only), as
				 * a fraction of the field: `s max|dpsi~| = step max|psi~|`,
				 * central unless @ref central is cleared. At 1e-5 in double the
				 * central difference is 1.2e-7 of the Jacobian from the shape
				 * derivative, and its truncation optimum is near 1e-6 to 1e-7
				 * (M-408). In single no step does better than 4e-5 on the box
				 * and 1e-3 on a free edge (M-406).
				 */
				double step = std::numeric_limits<mfem::real_t>::epsilon() < 1.0e-10 ? 1.0e-5 : 1.0e-4;
				bool central = true;
				/**
				 * THE DIFFERENCED EXTRACTION IN PLACE OF THE SHAPE DERIVATIVE:
				 * 2m extractions of the perturbed families at @ref step, the
				 * Jacobian the shape derivative is held to (M-408). Off by
				 * default; differentiate() takes the derivative of the
				 * extraction's own quadrature.
				 */
				bool differenced = false;
				/// Where the shape derivative gathers its node values: the
				/// device under any Device that allows one (Automatic), or the
				/// host loop -- the same loop, so `debug` is the host to the bit.
				ContinuousInterpolant::Route gatherRoute = ContinuousInterpolant::Route::Automatic;

				/// The solve, which may throw; solver.solve() when empty.
				std::function<void()> solve;
			};

			/// Borrows both; @a source is the solver's plasma source.
			SafetyFactorEquilibrium( GradShafranovSolver &solver, NormalisedMHDSource &source,
			                         Options options );

			/// The knots' labels, `Psi_N` ascending.
			std::vector<double> labels() const;

			/// Solve at @a spline, warm-started, and extract. False if the
			/// solve, the post-processing or the extraction fails; the solver
			/// then holds no admitted state until the next success.
			bool evaluate( ToroidalFieldSpline const &spline, SafetyFactorEvaluation &out );

			/// The derivative at the state the last evaluate() admitted,
			/// which must have been at @a spline.
			bool differentiate( ToroidalFieldSpline const &spline, SafetyFactorDerivative &out );

			/// The pinned scalar of the solver's current state at @a spline:
			/// the delivered flux or current; `g_edge` under EdgeField.
			double measure( ToroidalFieldSpline const &spline ) const;

			/// d psi_ax/d u_j from the last differentiate(): the border's part
			/// of each tangent -- what says which direction is the scale.
			std::vector<double> const &axisDerivative() const { return lastAxisDerivative; }

			/// Warm-start the next evaluate() from the state the solver holds
			/// now -- a solve made before the loop, at any g.
			void warmFromSolver();

			/// Why the last refused evaluate() was refused: the exception the
			/// solve or the extraction threw. Empty before any refusal.
			std::string const &refusal() const { return lastRefusal; }

			/// d psi_bnd/d u_j likewise; zeros where psi_bnd is not free.
			std::vector<double> const &boundaryDerivative() const { return lastBoundaryDerivative; }

			/// Inner Newton iterations summed over every solve, and the
			/// seconds spent in solves, tangents and extractions.
			int innerIterations() const { return newtonIterations; }
			double solveSeconds() const { return secondsSolving; }
			double tangentSeconds() const { return secondsTangents; }
			double extractionSeconds() const { return secondsExtracting; }
			int extractions() const { return extractionCount; }
			/// The last shape derivative's nodes, and how many of them sat in
			/// the band, where the fields are sampled and differenced rather
			/// than gathered from an element's polynomial.
			int shapeNodes() const { return lastShapeNodes; }
			int shapeBandNodes() const { return lastShapeBandNodes; }
			/// Whether the last shape derivative's gathers ran on the device.
			bool shapeGatheredOnDevice() const { return lastShapeOnDevice; }
			/// Options::gatherRoute, for the next differentiate(): a parity
			/// check differentiates one admitted state on each route.
			void setGatherRoute( ContinuousInterpolant::Route route ) { options.gatherRoute = route; }

		private:
			struct Fields
			{
				std::unique_ptr<mfem::GridFunction> psiStar, qStar, trace, psi, flux;
			};

			Fields solvedFields() const;
			Fields along( mfem::Vector const &dx ) const;
			std::vector<double> qOverG( Fields const &fields, double psiAxis ) const;
			std::vector<double> qOverG( Fields const &fields, CriticalPoint const &axis ) const;
			std::vector<double> qOverG( mfem::GridFunction const &potential, mfem::GridFunction const &flux,
			                            CriticalPoint const &axis, double boundary ) const;
			CriticalPoint locatedAxis( Fields const &fields ) const;
			FluxSurfaceFamily extractFamily( mfem::GridFunction const &potential, mfem::GridFunction const &flux,
			                                 CriticalPoint const &axis, double boundary ) const;
			/**
			 * d( q/g )_i/du_j of the extraction's own quadrature, by
			 * differentiating each node: the ray root's motion, the origin's,
			 * and the integrand's (SafetyFactorShape.cpp). @a originR and
			 * @a originZ are the axis root's motion along each tangent.
			 */
			bool shapeDerivative( ContinuousFields const &base, CriticalPoint const &axis, double boundary,
			                      std::vector<mfem::Vector> const &dx, std::vector<double> const &dAxis,
			                      std::vector<double> const &dBoundary, std::vector<double> const &originR,
			                      std::vector<double> const &originZ, std::vector<double> &out );

			GradShafranovSolver &solver;
			NormalisedMHDSource &source;
			Options options;
			mfem::Vector warm;
			std::vector<double> warmExterior;
			double warmBoundary = 0.0;
			double warmAxis = 0.0;
			mfem::Vector admitted;
			std::vector<double> admittedValues;
			std::vector<double> lastAxisDerivative;
			std::vector<double> lastBoundaryDerivative;
			std::string lastRefusal;
			int newtonIterations = 0;
			double secondsSolving = 0.0;
			double secondsTangents = 0.0;
			mutable double secondsExtracting = 0.0;
			mutable int extractionCount = 0;
			int lastShapeNodes = 0;
			int lastShapeBandNodes = 0;
			bool lastShapeOnDevice = false;
	};

}

#endif
