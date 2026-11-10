#ifndef MEQ_SAFETY_FACTOR_TANGENT_HPP
#define MEQ_SAFETY_FACTOR_TANGENT_HPP

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

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
 * the family extracted on the continuous fields of `x +- s dx_j` with the
 * labels `psi_ax +- s dpsi_ax,j` and `psi_bnd +- s dpsi_bnd,j` and the ray
 * origin MOVED WITH THE AXIS: the map
 * roots q_h for its axis, so the origin moves by `-( dq_h/dx )^-1 dq_j` in the
 * axis's own element, the derivative of that root. Holding it fixed misses
 * the extraction's dependence on where its rays start, 1e-5 to 4e-5 of the
 * Jacobian (MEASUREMENTS.md M-QR-1); re-locating it on each perturbed field
 * can pick the other element's root where the axis sits on a face, and
 * difference the O( h^( k + 1 ) ) jump between them. psi~ and q~ are linear
 * in the state, so the perturbed pair is the base's plus `s` times the
 * tangent's -- the solver's kept interpolant applied m + 1 times (GradShafranovSolver::continuousInterpolant()) -- and the extraction is differenced
 * centrally, 2m extractions; no equilibrium is re-solved (M-365, M-367,
 * M-QR-9). The extractions are most of the cost, and
 * the analytic shape derivative would remove them (TODO, *The q-driven
 * Jacobian's cost is the extraction*).
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
 * (MEASUREMENTS.md M-QR-5).
 */

namespace meq
{

	/// The mesh's area in cells: a quadrilateral one, a triangle one half --
	/// toroidalKnotCount()'s argument on a fixed boundary, where the plasma
	/// is the whole of Omega_h.
	double plasmaCells( mfem::Mesh const &mesh );

	/// The plasma's area in cells where the solver's source is confined to a
	/// support it has decided (GradShafranovSolver::elementInPlasma()) -- a
	/// limited or diverted plasma inside a vacuum region -- else the mesh's.
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
				 * The extraction's difference step, as a fraction of the
				 * field: `s max|dpsi~| = step max|psi~|`, central unless
				 * @ref central is cleared. CENTRAL, AT 1e-5 IN DOUBLE AND 1e-4
				 * IN SINGLE. One-sided at 1e-9 the columns are closer to
				 * re-solved differences (3.1e-8 against 1.2e-7) for m + 1
				 * extractions in place of 2m, but on a target with a misfit
				 * the error's projection through the near-null scale direction
				 * left the removable residual stalled at 1e-9 to 3e-9, above
				 * the tolerance, on QDrivenScalePin's base (M-QR-9). In single
				 * a step under the format's epsilon adds nothing to the field:
				 * at 1e-4 SafetyFactorSolver's round trip converges in 4
				 * steps, at 1e-3 its line search finds no decrease and at
				 * 1e-2 a perturbed family cannot be traced (M-QR-10).
				 */
				double step = std::numeric_limits<mfem::real_t>::epsilon() < 1.0e-10 ? 1.0e-5 : 1.0e-4;
				bool central = true;

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

			GradShafranovSolver &solver;
			NormalisedMHDSource &source;
			Options options;
			mfem::Vector warm;
			std::vector<double> warmExterior;
			double warmBoundary = 0.0;
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
	};

}

#endif
