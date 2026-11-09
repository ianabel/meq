#ifndef MEQ_SAFETY_FACTOR_TANGENT_HPP
#define MEQ_SAFETY_FACTOR_TANGENT_HPP

#include <functional>
#include <memory>
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
 * the retained, bordered linearisation -- one back-solve. Then `psi*`, `q*`
 * and `psihat` of the tangent (reconstructFields(), linear), and the family
 * extracted on the continuous fields of `x +- s dx_j` with the label
 * `psi_ax +- s dpsi_ax,j` and the ray origin MOVED WITH THE AXIS: the map
 * roots q_h for its axis, so the origin moves by `-( dq_h/dx )^-1 dq_j` in the
 * axis's own element, the derivative of that root. Holding it fixed misses
 * the extraction's dependence on where its rays start, 1e-5 to 4e-5 of the
 * Jacobian (MEASUREMENTS.md M-QR-1); re-locating it on each perturbed field
 * can pick the other element's root where the axis sits on a face, and
 * difference the O( h^( k + 1 ) ) jump between them. The extraction is differenced, centrally:
 * no equilibrium is re-solved, and the continuous fields make the family
 * smooth in `s` (M-365, M-367). The analytic shape derivative would remove
 * the `2m` extractions per Jacobian that are most of the cost (TODO, *The
 * q-driven Jacobian's cost is the extraction*).
 *
 * The pins: R0B0's row reads the spline alone; the flux and the current are
 * linearised as the quadratures they are read with
 * (GradShafranovSolver::integrateOverDomain() with a direction,
 * linearisedPlasmaCurrent()).
 *
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

				/// The extraction's difference step, as a fraction of the
				/// field: `s max|dpsi*| = step max|psi*|`.
				double step = 1.0e-5;

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
			static Fields combine( Fields const &base, Fields const &direction, double s );
			std::vector<double> qOverG( Fields const &fields, double psiAxis ) const;
			std::vector<double> qOverG( Fields const &fields, CriticalPoint const &axis ) const;
			CriticalPoint locatedAxis( Fields const &fields ) const;

			GradShafranovSolver &solver;
			NormalisedMHDSource &source;
			Options options;
			mfem::Vector warm;
			mfem::Vector admitted;
			std::vector<double> admittedValues;
			std::vector<double> lastAxisDerivative;
			int newtonIterations = 0;
			double secondsSolving = 0.0;
			double secondsTangents = 0.0;
			mutable double secondsExtracting = 0.0;
			mutable int extractionCount = 0;
	};

}

#endif
