#ifndef MEQ_ESTIMATOR_HPP
#define MEQ_ESTIMATOR_HPP

#include <array>
#include <memory>

#include "mfem.hpp"

#include "GradShafranov.hpp"
#include "Source.hpp"

/*
 * The residual error estimator of refs/HDG-GradShafranov-Adaptive.pdf eq (20),
 * and the two marking strategies of its section 3.2.
 *
 * The estimator is five terms. For an element K of the triangulation, with h_K
 * its diameter, h_e the length of an edge e, E the whole skeleton and E_I its
 * interior part:
 *
 *   eta_K^2 = h_K^2 || F( psi* )/r + div_bar q_h ||_K^2                 (eta_1)
 *           + || q_h - ( 1/r ) grad_bar psi* ||_K^2                     (eta_2)
 *           + ( 1/2 ) sum_{e in E_I, e in dK} h_e    || [[ q_h ]]  ||_e^2   (eta_3)
 *           + ( 1/2 ) sum_{e in E_I, e in dK} h_e^-1 || [[ psi* ]] ||_e^2   (eta_4)
 *           +         sum_{e in E,   e in dK} h_e^-1 || psihat_h - psi* ||_e^2  (eta_5)
 *
 * and eta^2 = sum_K eta_K^2. The jumps are the paper's own: [[ q ]] :=
 * q+ . n+ + q- . n- for the vector, [[ a ]] := a+ - a- for the scalar, so
 * [[ q ]] is a scalar and both are independent of which side is called plus.
 *
 * eta_1 is the residual of the strong divergence equation -div_bar q = F/r,
 * eta_2 that of the constitutive law q = ( 1/r ) grad_bar psi, eta_3 and eta_4
 * measure the loss of conformity of the flux and of the potential, and eta_5 the
 * disagreement between the hybrid unknown and the potential on element
 * boundaries. Each vanishes identically on the exact solution, which is the
 * property the rate table in tests/convergence/EstimatorConvergence.cpp checks:
 * a term whose sign or weight is wrong does not converge at all.
 *
 * PSI-STAR, NOT PSI-H, AND WHY THAT IS THE WHOLE POINT.
 *
 * The potential in eta_1, eta_2, eta_4 and eta_5 above is psi*, the
 * post-processed solution in P_(k+1) -- GradShafranovSolver::postProcess(). The
 * paper is explicit about eta_2: built on the raw psi_h it converges at REDUCED
 * order, because it differentiates the approximation, and substituting psi*
 * is what preserves k+1. meq's pre-modernisation estimator used raw psi_h
 * throughout and was a degraded copy of the published one.
 *
 * That is a claim, and this class is arranged so that it is a measurement:
 * setPotential( Potential::Raw ) rebuilds every psi-dependent term on psi_h, so
 * the two can be run side by side on the same solution and on the same
 * quadrature rule. Measured on the fitted Solov'ev benchmark, eta_2 converges at
 *
 *       k     on psi*   on psi_h
 *       1       2.002      0.998
 *       2       3.000      2.002
 *       3       3.992      2.981
 *
 * -- exactly the one order the paper predicts, at every k, and a factor of 124
 * to 407 in absolute size on the finest mesh. So the claim holds, and
 * Potential::Raw exists to keep on holding it rather than as an option anybody
 * should use.
 *
 * ETA_5 AS PRINTED DOES NOT VANISH ON THE EXACT SOLUTION.
 *
 * This is the one place meq departs from eq (20) as written, and the reason is
 * provable rather than a matter of taste. psihat_h lives in M_h = P_k( e );
 * psi*'s trace on e is a polynomial of degree k+1, which no element of M_h can
 * represent. Writing P_M for the L2( e ) projection onto M_h, orthogonality
 * splits the printed term as
 *
 *     || psihat_h - psi* ||_e^2 = || psihat_h - P_M psi* ||_e^2
 *                              + || ( I - P_M ) psi*   ||_e^2
 *
 * and the second piece is not an error at all: put the exact psi in for psi* and
 * its best trace P_M psi in for psihat_h, and it is still O( h^(k+3/2) ), which
 * h_e^-1 and O( h^-2 ) edges turn into an O( h^k ) floor. Measured, the printed
 * term does not merely sit above that floor, it IS the floor, agreeing with it
 * to between 0.05 and 1.8 per cent at every k and every mesh -- so it converges
 * at k, not k+1, and drags the total down with it.
 *
 * Taking the difference inside M_h removes exactly the piece that carries no
 * information, restores k+1, and reproduces the rates GS-2 Table 1 reports for
 * eta_5. It is also what the paper's own prose asks for: eta_5 compares the
 * hybrid variable and the post-processed solution "as approximations to the local
 * trace", and comparing a P_k( e ) function with a degree-(k+1) trace as
 * approximations to the same thing means comparing them in P_k( e ). So
 * TraceComparison::Projected is the default, and Literal is kept only so that
 * tests/convergence/EstimatorConvergence.cpp can keep measuring the difference.
 * That test carries the numbers and the argument in full.
 *
 * ON THE EXTENSION PATH, CALL setTransferredBoundary(). Four of the five terms
 * are fine there -- psi* still converges at k+2 even though
 * DarcyForm::Reconstruct() drops the boundary-face integrator carrying the
 * transferred datum, which was measured rather than assumed and came out the
 * opposite way round from the guess. eta_5 is not fine, for a different reason:
 * on Gamma_h psihat_h is pinned to zero and the datum actually imposed is phi_h,
 * which is never stored. eta_5 then compares psi* against zero on those faces
 * and the difference is O( dist( Gamma_h, Gamma ) ) = O( h ). Measured at k = 2,
 * eta becomes 4.09e-1 where eta_1 is 2.12e-3 and converges at about a half, so
 * the estimator is nothing but that one term. setTransferredBoundary() leaves
 * those faces out, which restores k+1 for eta and every component, and is an
 * omission rather than a repair -- see that method.
 *
 * AND NOT AT ALL THROUGH NEWTON. GradShafranovSolver::postProcess() refuses the
 * semi-linear path, because DarcyForm::Reconstruct() silently returns ~1e15
 * there. So this estimator is available on the linear path only, on either
 * boundary regime. The eta_1 constructor still takes a meq::Source, and should:
 * eta_1 is the residual of the semi-linear equation whatever solved it, and the
 * day the reconstruction works through Newton nothing here changes.
 */

namespace meq
{

	/**
	 * eq (20), as an mfem::ErrorEstimator, with the five terms kept apart.
	 *
	 * Borrows the solver and the source; both must outlive it. The solver must
	 * have been solved, and postProcess() must have been called on it, before
	 * GetLocalErrors() -- which is checked rather than assumed, because a
	 * post-processed potential left at zero would give a plausible-looking
	 * estimator that was really just measuring psi_h.
	 *
	 * A single total can hide one term being wrong; five rates cannot, which is
	 * why component() and localSquares() are part of the interface rather than
	 * diagnostics bolted on.
	 *
	 * The estimator caches. It recomputes when the mesh sequence changes, or when
	 * one of the setters or Reset() is called -- and NOT when the borrowed solver
	 * is solved again on the same mesh, which the mesh sequence cannot see. Call
	 * Reset() after a second solve, or build a fresh estimator, which is what the
	 * adaptive loop in tests/convergence/AdaptiveRefinement.cpp does.
	 */
	class ResidualEstimator : public mfem::ErrorEstimator
	{
		public:
			/// The five terms of eq (20), in the order the paper numbers them.
			enum class Term
			{
				Divergence = 0,   ///< eta_1, the residual of -div_bar q = F/r
				Constitutive,     ///< eta_2, the residual of q = grad_bar psi / r
				FluxJump,         ///< eta_3, [[ q_h ]] across interior edges
				PotentialJump,    ///< eta_4, [[ psi* ]] across interior edges
				TraceMismatch     ///< eta_5, psihat_h against psi* on dK
			};

			/// Which potential the psi-dependent terms are built on. See the file
			/// comment: PostProcessed is the published estimator, Raw is the
			/// degraded one, kept so that the difference stays measured.
			enum class Potential
			{
				PostProcessed,
				Raw
			};

			/// How eta_5 compares the hybrid unknown with the potential. See
			/// ETA_5 AS PRINTED DOES NOT VANISH ON THE EXACT SOLUTION in the file
			/// comment; Projected is the default and Literal is there to keep that
			/// finding measured.
			enum class TraceComparison
			{
				Projected,   ///< || psihat_h - P_M( psi* ) ||, inside M_h
				Literal      ///< || psihat_h - psi* ||, eq (20) as printed
			};

			static constexpr int termCount = 5;

			/// @param solverIn  solved, and post-processed if the potential is
			///                  PostProcessed. Borrowed.
			/// @param sourceIn  F( r, z, psi ). Borrowed. eta_1 evaluates it at
			///                  the potential in use, which is what makes the
			///                  term the residual of the semi-linear equation
			///                  rather than of a frozen one.
			ResidualEstimator( GradShafranovSolver &solverIn, Source const &sourceIn );

			/// The same, for a source that does not depend on psi. F is then
			/// evaluated as the coefficient F( r, z ), so eta_1 is the residual of
			/// the linear equation -- which is the same thing when F is the same
			/// thing, and is what the Solov'ev benchmark needs.
			ResidualEstimator( GradShafranovSolver &solverIn, mfem::Coefficient &sourceIn );

			ResidualEstimator( ResidualEstimator const & ) = delete;
			ResidualEstimator &operator=( ResidualEstimator const & ) = delete;

			void setPotential( Potential potentialIn );
			Potential potential() const;

			void setTraceComparison( TraceComparison comparisonIn );
			TraceComparison traceComparison() const;

			/**
			 * Boundary attributes whose Dirichlet datum is TRANSFERRED rather than
			 * stored, so that eta_5 leaves their faces out.
			 *
			 * Pass GradShafranovSolver::setExtension()'s own Gamma_h marker. On
			 * such a face psihat_h is not the condition that was imposed: phi_h is,
			 * and phi_h is a line integral of the extended flux that nothing keeps
			 * -- the trace dofs there are pinned to zero because nothing references
			 * them. eta_5 on those faces therefore compares psi* against zero, and
			 * that difference is O( dist( Gamma_h, Gamma ) ) = O( h ), not
			 * O( h^(k+2) ).
			 *
			 * MEASURED, on the stage-5 benchmark at k = 2, h = 0.213 / 0.106 /
			 * 0.053: with those faces in, eta_5 is 4.09e-1, 3.21e-1, 2.24e-1 -- a
			 * rate of about a half -- against eta_1 at 2.12e-3, 3.39e-4, 4.74e-5.
			 * So the term is four orders larger than the rest of the estimator and
			 * converging at h^(1/2), and eta is nothing but that. The other four
			 * terms are unaffected and keep k+1.
			 *
			 * Leaving those faces out is an OMISSION, not a fix: the paper's eta_5
			 * does sum over boundary edges, with b_h = phi_h there. The proper
			 * repair is to evaluate phi_h -- mfem::PathTraceCoefficient and
			 * mfem::ExtensionRegionQuadrature are the pieces -- and it is written up
			 * in the stage-6 report rather than guessed at here. Until then the
			 * elements along Gamma_h are still covered by eta_1 to eta_4.
			 *
			 * @param markerIn  sized by the largest boundary attribute, 1 where the
			 *                  datum is transferred. Copied. An empty array, the
			 *                  default, excludes nothing.
			 */
			void setTransferredBoundary( mfem::Array<int> const &markerIn );

			/// Quadrature order added to twice the degree of the potential in use.
			/// The default of 4 matches the rule the convergence tests measure the
			/// L2 error on, for the same reason: neither F/r nor 1/r is a
			/// polynomial, and a rule chosen for the polynomial part alone would
			/// be what limited the measured rate.
			void setExtraQuadratureOrder( int extraIn );

			/// eta_K, per element. MFEM's spelling, from mfem::ErrorEstimator.
			mfem::Vector const &GetLocalErrors() override; // NOLINT(readability-identifier-naming)

			/// eta, the square root of the sum over the mesh. MFEM's spelling,
			/// from mfem::ErrorEstimator. Const, and the base class says so, so
			/// the cache is mutable; it will compute if it has to.
			mfem::real_t GetTotalError() const override; // NOLINT(readability-identifier-naming)

			/// MFEM's spelling, from mfem::ErrorEstimator.
			void Reset() override; // NOLINT(readability-identifier-naming)

			/// eta_i, one term, summed over the mesh and square rooted.
			double component( Term term ) const;

			/// The per-element contribution of one term, SQUARED -- which is the
			/// additive quantity, and the one Doerfler marking needs.
			mfem::Vector const &localSquares( Term term ) const;

			/// "eta_1" ... "eta_5", for a table heading.
			static char const *name( Term term );

		private:
			void compute() const;

			/// F at a point, from whichever of the two constructors was used.
			double sourceValue( mfem::ElementTransformation &tr,
			                    mfem::IntegrationPoint const &ip,
			                    double r, double z, double psi ) const;

			GradShafranovSolver *solver;
			Source const *source;
			mfem::Coefficient *sourceCoeff;
			Potential potentialChoice;
			TraceComparison comparisonChoice;
			mfem::Array<int> transferredBoundary;
			int extraQuadratureOrder;

			/// Mesh::GetSequence() at the last compute, so that a refinement
			/// invalidates the cache without anyone having to remember to.
			mutable long sequence;
			mutable mfem::Vector errors;
			mutable std::array<mfem::Vector, termCount> squares;
			mutable std::array<double, termCount> sums;
	};

	/**
	 * Doerfler (bulk) marking, GS-2 section 3.2 criterion A.
	 *
	 * Marks a minimal set M with sum_{K in M} eta_K^2 >= gamma sum_K eta_K^2,
	 * obtained by sorting the elements by eta_K and taking the largest until the
	 * bulk is reached. gamma near zero marks few elements and refines very
	 * locally; gamma near one marks nearly all and approaches uniform refinement.
	 *
	 * This is the criterion the convergence proof assumes -- Cockburn, Nochetto
	 * and Zhang, refs/Refs.md -- which is the reason to have it even though the
	 * paper's own experiments used the maximum criterion below.
	 *
	 * @param localErrors  eta_K, as ResidualEstimator::GetLocalErrors() returns
	 *                     it. NOT the squares: this squares them itself, so that
	 *                     both marking functions take the same argument and
	 *                     neither can be fed the wrong one.
	 * @param gamma        the marking parameter, in ( 0, 1 ].
	 * @param marked       filled with the indices of the marked elements, in
	 *                     decreasing order of eta_K.
	 */
	void markDoerfler( mfem::Vector const &localErrors, double gamma,
	                   mfem::Array<int> &marked );

	/**
	 * Maximum marking, GS-2 section 3.2 criterion B.
	 *
	 * Marks every K with eta_K >= gamma max_K eta_K. The behaviour with gamma is
	 * the reverse of Doerfler's: large gamma refines locally, small gamma
	 * approaches uniform. GS-2's own experiments used this with gamma = 0.3, and
	 * note that its convergence analysis for HDG is still open.
	 *
	 * @param localErrors  eta_K, as GetLocalErrors() returns it.
	 * @param gamma        the marking parameter, in [ 0, 1 ].
	 * @param marked       filled with the indices of the marked elements, in
	 *                     increasing order of element index.
	 */
	void markMaximum( mfem::Vector const &localErrors, double gamma,
	                  mfem::Array<int> &marked );

	/**
	 * The computational and companion meshes of GS-2 section 3.3, and the local
	 * refinement that keeps dist( Gamma_h, Gamma ) = O( h_loc ).
	 *
	 * WHY THIS IS NOT STANDARD ADAPTIVE MESH REFINEMENT, which is the whole
	 * content of section 3.3. On the extension path Omega_h is a STRICT
	 * subdomain of Omega -- the background elements lying entirely inside it --
	 * and the transfer of the Dirichlet datum is only optimal while the gap
	 * between Gamma_h and Gamma stays of the order of the LOCAL element
	 * diameter. Refine an element of T_h and its children are still inside
	 * Omega, so Gamma_h does not move at all: the gap stays where it was while
	 * h_loc halves, and dist/h_loc doubles with every cycle. That is the paper's
	 * Figure 3, and tests/convergence/AdaptiveRefinement.cpp measures it: the
	 * ratio dist/h_loc goes 0.98, 1.97, 3.94, 7.88 over three cycles with the
	 * domain held fixed -- a clean factor of two each time -- against 0.98, 1.02,
	 * 1.39, 1.36 with the update in place.
	 *
	 * The fix needs a second mesh. The COMPANION mesh T_c^h is the minimal cover
	 * of Omega-bar: every background element that meets Omega, so T_h together
	 * with the band of elements that Gamma cuts through. Nothing is ever computed
	 * on it. Its only job is to be the thing that gets refined, so that the band
	 * subdivides and some of the children fall inside Omega -- which is how T_h
	 * grows towards Gamma instead of standing still.
	 *
	 * The five steps of section 3.3, and where each one is:
	 *
	 *   1) build T_h and T_c^h from the background mesh        the constructor
	 *   2) solve and estimate on T_h, giving a marked set M    the caller
	 *   3) mark the companion elements matching M, PLUS every  refine()
	 *      companion element that Gamma cuts and that shares
	 *      an edge with something in M -- giving M_c
	 *   4) refine M_c, with no hanging nodes                   refine()
	 *   5) re-select T_h and T_c^h from the refined mesh       refine()
	 *
	 * Step 3's second half is the part that does the work: it is what pushes the
	 * refinement out into the band, and without it steps 4 and 5 are just AMR on
	 * T_h.
	 *
	 * READ THE WARNING AT THE END OF SECTION 3.3 BEFORE READING A PICTURE OF THE
	 * OUTPUT. The refinement will look as though it is crowding the boundary.
	 * That is step 3's second half doing its job -- the local proximity condition
	 * forces new, necessarily smaller, elements to appear at Gamma whatever the
	 * indicator says -- and it is NOT the error estimator concentrating there.
	 * Reading it as the latter sends you looking for a bug in the estimator that
	 * is not there.
	 *
	 * ONE DELIBERATE SIMPLIFICATION. Step 5 in the paper discards the background
	 * elements that end up in neither mesh; this keeps the whole background box.
	 * That costs memory and nothing else, and it buys something: conforming
	 * refinement of triangles PROPAGATES, and a propagation that runs off the
	 * edge of a trimmed mesh has nowhere to go. Keeping the box means it always
	 * has somewhere.
	 */
	class AdaptiveDomain
	{
		public:
			/// @param backgroundIn  a uniform, shape-regular triangulation of a box
			///                      strictly containing Omega. Copied, because
			///                      refine() mutates it.
			/// @param levelSetIn    negative inside Omega. Copied.
			/// @param extraRefineIn how deeply to sample each element when
			///                      deciding inside/meets, as
			///                      mfem::MarkLevelSetSubdomain's own argument:
			///                      0 tests the vertices, which is exact only for
			///                      a convex Omega. The default of 1 is the cheap
			///                      insurance ExtensionConvergence.cpp also takes.
			AdaptiveDomain( mfem::Mesh const &backgroundIn,
			                mfem::PositionFunction levelSetIn,
			                int extraRefineIn = 1 );

			AdaptiveDomain( AdaptiveDomain const & ) = delete;
			AdaptiveDomain &operator=( AdaptiveDomain const & ) = delete;

			/// T_h, as a SubMesh of the background. Invalidated by refine(),
			/// which builds a new one -- so a solver, a transfer path or a
			/// GridFunction built on it must be rebuilt too.
			mfem::SubMesh &computational();

			/// The boundary attribute of Gamma_h on that SubMesh: the part SubMesh
			/// had to generate, which is what
			/// GradShafranovSolver::setExtension() wants marked. Every other
			/// attribute is inherited from the background box and is fitted.
			int gammaHAttribute() const;

			/// A marker over the largest boundary attribute, with Gamma_h set,
			/// ready for setExtension().
			mfem::Array<int> const &gammaHMarker() const;

			int numComputational() const;
			int numCompanion() const;
			int numBackground() const;

			/// How many times refine() has been called.
			int refinements() const;

			/// The largest and smallest element diameter of T_h. On a locally
			/// refined mesh these are the two ends of a spread, and it is the
			/// small one that a transfer path's search length has to respect.
			double largestElement() const;
			double smallestElement() const;

			/// Steps 3 to 5. @a marked indexes elements of computational().
			void refine( mfem::Array<int> const &marked );

			/// Steps 4 and 5 with step 3's second half left out: refine only the
			/// background elements that carry the marked elements of T_h. This is
			/// plain AMR on the computational mesh, it is what section 3.3 says
			/// does not work, and it is here to be measured failing rather than
			/// argued about.
			void refineWithoutCompanion( mfem::Array<int> const &marked );

			/// How many elements the last refine() added to the marked set through
			/// step 3's second half -- the band Gamma cuts -- rather than through
			/// the indicator. This is the warning at the end of section 3.3 as a
			/// number: refinement that appears to crowd Gamma is this count at
			/// work, and it is not the estimator concentrating there.
			int lastProximityAdditions() const;

		private:
			void select();
			void refineBackground( mfem::Array<int> const &marked, bool useCompanion );

			mfem::Mesh backgroundMesh;
			mfem::PositionFunction levelSet;
			int extraRefine;

			/// Over the background elements: 1 for T_h, 1 for T_c^h.
			mfem::Array<int> insideMarker;
			mfem::Array<int> companionMarker;

			std::unique_ptr<mfem::SubMesh> computationalMesh;
			int gammaH;
			mfem::Array<int> gammaHMarkerValue;
			int refinementCount;
			int proximityAdditions;
	};

	/// The diameter of a mesh element: the largest distance between two of its
	/// vertices, which for a simplex and for any convex element is the diameter
	/// itself. NOT mfem::Mesh::GetElementSize(), which is a Jacobian-derived
	/// length scale and is a different number on an anisotropic triangle. This is
	/// the h_K and h_e of eq (20) and the h_loc of section 3.3, so it is shared
	/// rather than written twice.
	double elementDiameter( mfem::Mesh &mesh, int element );

	/// The same for a face, which in two dimensions is its length.
	double faceDiameter( mfem::Mesh &mesh, int face );

}

#endif // MEQ_ESTIMATOR_HPP
