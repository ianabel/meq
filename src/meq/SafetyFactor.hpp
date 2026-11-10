#ifndef MEQ_SAFETY_FACTOR_HPP
#define MEQ_SAFETY_FACTOR_HPP

#include <functional>
#include <memory>
#include <vector>

#include "meq/FluxFamily.hpp"
#include "meq/Profiles.hpp"

/**
 * @file SafetyFactor.hpp
 * DRIVING AN EQUILIBRIUM BY `q( psi )` INSTEAD OF BY `g( psi )`.
 *
 * `ROADMAP.md` item 10. Every source in this tree takes the toroidal field
 * function as input -- `[source] GGPrimeFile` is `g dg/dPsi` -- and reports the
 * safety factor as an output. A transport code hands an equilibrium code the
 * other way round: `q( psi )` is the target and `g( psi )` is what has to be
 * found. This file is that inversion.
 *
 *
 * 1. THE ALGEBRA IS A DIVISION AND THE DIFFICULTY IS ALL IN THE GEOMETRY
 *
 * `meq::SurfaceAverages::safetyFactor` is
 *
 *     q = V' g < R^-2 > / 4 pi^2
 *
 * so
 *
 *     g = 4 pi^2 q / ( V' < R^-2 > )
 *
 * and that is the whole of the inversion at FIXED geometry. It is not a
 * quadrature, it is not an integral equation, and nothing here converges to
 * anything: given a family it is one division per surface.
 *
 * WHAT MAKES IT A SOLVER IS THAT `V'` AND `< R^-2 >` ARE FUNCTIONALS OF THE
 * SOLUTION, and the solution depends on `g`. So `g` is an UNKNOWN: it is
 * carried as its values `u_j = g^2` at a set of surfaces -- the KNOTS, which
 * refine with the mesh -- interpolated by a not-a-knot cubic in the source's
 * `Psi` (ToroidalFieldSpline), and the conditions are `q` at every knot and
 * one scale pin (SafetyFactorSystem). They are solved together by
 * Gauss-Newton (SafetyFactorSolve.hpp) on the linearised solve's EXACT
 * tangent: one bordered back-solve per knot, and the extraction's response
 * to each tangent field (SafetyFactorTangent.hpp). No equilibrium is
 * re-solved to form a Jacobian.
 *
 * THE DELIVERED q CONVERGES TO THE REQUESTED q AS h AND k REFINE, which is
 * the contract (MEASUREMENTS.md M-393): the knots follow the mesh
 * (toroidalKnotCount()), and so does the part of the request near the axis
 * and the edge that no knot sees (toroidalKnotFamily()).
 *
 *
 * 2. THE TWO NORMALISED FLUXES ARE OPPOSITE, AND THIS IS THE TRAP
 *
 * **`meq::FluxSurfaceFamily` and `meq::NormalisedSource` BOTH call their
 * argument the normalised flux and they run in opposite directions.**
 *
 *     family    Psi_N = ( psi_ax - psi )/span     0 on the axis, 1 at the edge
 *     source    Psi   = ( psi - psi_bnd )/span    1 on the axis, 0 at the edge
 *
 * so `Psi = 1 - Psi_N` and, which is the half that bites,
 *
 *     d/dPsi = - d/dPsi_N
 *
 * A `gg'` built by differentiating `g^2` against the family's label and handed
 * to the source **without that sign** is not a small error and it does not
 * fail: it is a plausible equilibrium with the shear reversed, which is a real
 * thing a real machine can have, so nothing downstream looks wrong. The
 * conversion happens once, here, in ggPrimeKnots(), and
 * `the_two_normalised_fluxes_run_in_opposite_directions` is the assertion that
 * keeps it.
 *
 *
 * 3. THE SIGN OF `g` IS A CONVENTION AND IT IS NOT RECOVERABLE FROM `q`
 *
 * `V'` and `< R^-2 >` are positive by construction, so `g` and `q` share a
 * sign -- and the source consumes `g dg/dPsi`, which is `(g^2)'/2` and cannot
 * see it at all. Reversing the toroidal field reverses `q` and leaves the
 * equilibrium alone. So this file takes `q` at the sign it is given, returns
 * the `g` that matches, and records that `g^2` is what actually reaches the
 * solve. A caller wanting the other field direction negates its own `q` and
 * gets the same equilibrium.
 *
 *
 * 4. THE CUT IS THE MODELLING DECISION, AND IT IS THE CALLER'S
 *
 * A family covers `Psi_N` in `[ innerLabel, outerLabel ]^2` -- `[ 0.05, 0.95 ]`
 * as shipped -- and **refuses to extrapolate**, which is IN-6's decision and
 * the right one: what changes across the cut is not the accuracy but what the
 * surface is MADE OF, and the band mask is the only thing that can say so.
 *
 * A source, though, is evaluated at every quadrature point of the domain,
 * including on the axis and outside the plasma. So a `gg'` built from a family
 * has to say something where the family does not. `meq::SplineProfile`
 * CONTINUES LINEARLY outside its knots -- the end knot's value and slope
 * carried on, so `dF/dpsi` is continuous there -- so what is carried past the
 * ends is the end knot's.
 *
 * **THAT IS A CHOICE ABOUT THE PHYSICS AND IT IS WRONG AT THE EDGE.** Outside
 * the plasma the toroidal field is the vacuum one, `g = const`, so `gg' = 0`;
 * carrying the outermost knot on instead keeps `g^2` growing. Whether
 * that matters depends on whether the source is confined -- `[source]
 * ConfineToPlasma` sets `F = 0` outside the plasma outright and the question
 * does not arise -- which is why this file offers the choice rather than making
 * it. See EdgeExtension.
 */

namespace meq
{

	/// What the edge knot says about `gg'` outside the family's own range.
	enum class EdgeExtension
	{
		/**
		 * Carry the outermost knot on -- its value and slope, SplineProfile's
		 * linear continuation -- which is what SplineProfile does with no help. Right when the source is confined to the plasma, where
		 * nothing outside is ever evaluated.
		 */
		Clamp,

		/**
		 * Add a knot at `Psi = 0` carrying `gg' = 0`, so the toroidal field goes
		 * over to a constant at the plasma edge and the source dies with it.
		 * The physical statement, and the one to take when the source is NOT
		 * confined.
		 */
		VacuumOutside
	};

	/// `g( psi )` recovered from a target `q`, against the FAMILY's label.
	///
	/// Everything here is in `Psi_N`, zero on the axis, because that is what the
	/// family that produced it is in. ggPrimeKnots() is the one place the
	/// reflection to the source's `Psi` happens.
	struct ToroidalField
	{
		/// The family's own labels, ascending, one entry per surface.
		std::vector<double> normalisedFlux;

		/// `q` as asked for, at those labels. Kept so that a caller can see
		/// what it requested beside what that implied.
		std::vector<double> safetyFactor;

		/// `g = 4 pi^2 q / ( V' < R^-2 > )`, and its square, which is the
		/// quantity the source is actually built from.
		std::vector<double> g;
		std::vector<double> gSquared;

		std::size_t size() const { return normalisedFlux.size(); }
		bool empty() const { return normalisedFlux.empty(); }
	};

	/**
	 * The `g` that would give @a family the safety factor @a target.
	 *
	 * One division per surface and no iteration: see section 1. @a target is
	 * called with the family's own `Psi_N`, zero on the axis.
	 *
	 * @throws std::invalid_argument if the family is empty, or if any surface
	 *         has a non-positive `V'` or `< R^-2 >` -- both are positive by
	 *         construction, so a non-positive one is a corrupt family rather
	 *         than a hard case, and dividing by it would return a signed
	 *         infinity that the spline below would carry into the source.
	 * @throws std::invalid_argument if @a target returns a non-finite value.
	 */
	ToroidalField invertSafetyFactor(
		FluxSurfaceFamily const &family,
		std::function<double( double normalisedFlux )> const &target );

	/**
	 * `gg'` for a meq::NormalisedMHDSource, in the SOURCE's `Psi`.
	 *
	 * Returns knots ascending in `Psi = 1 - Psi_N`, carrying `gg' = ( g^2 )'/2`
	 * and its own derivative, ready for `meq::SplineProfile`. The reflection of
	 * section 2 -- including the sign on the derivative -- happens here and
	 * nowhere else.
	 *
	 * THE SLOPES ARE FRITSCH-CARLSON LIMITED, which is the same rule
	 * meq::FluxSurfaceFamily interpolates with and is chosen for the same
	 * reason: `g^2` arrives as a table extracted from a solved field, so it
	 * carries that extraction's noise, and an unlimited cubic through noisy data
	 * OVERSHOOTS. An overshoot in `gg'` is a source term with a sign it should
	 * not have, fed back into the next solve of a fixed-point iteration -- which
	 * is how a loop that ought to converge instead oscillates.
	 *
	 * @throws std::invalid_argument if @a field has fewer than two surfaces, or
	 *         if its labels are not strictly ascending.
	 */
	std::vector<Knot> ggPrimeKnots( ToroidalField const &field,
	                                EdgeExtension extension
	                                    = EdgeExtension::Clamp );


	/**
	 * THE KNOT COUNT FOR A MESH: ABOUT TWO SURFACES PER ELEMENT ACROSS THE
	 * MINOR RADIUS.
	 *
	 * @a cells is the plasma's area in units of a cell -- a quadrilateral
	 * counts one and a triangle one half, so that a triangulated rectangle
	 * reads as its rectangle -- and `sqrt( cells/pi )` is then the number of
	 * elements across the minor radius of a disc of that area. The knots are
	 * twice that, so their spacing `Delta rho` follows `h`.
	 *
	 * WHY TWO AND WHY THE CUBIC IS ENOUGH TO k = 3. The delivered `g^2`'s
	 * error is the discretisation's, `O( h^( k + 1 ) )` in the averages, plus
	 * the spline's, `O( Delta^4 )`: with `Delta` proportional to `h` the
	 * spline keeps pace to `k = 3`. Above that the spline would cap the
	 * rate, so the count grows as `( a/h )^( ( k + 1 )/4 )` there -- measured
	 * at `k = 2` and 3 only (MEASUREMENTS.md M-393).
	 *
	 * At least six, at most 128.
	 *
	 * @throws std::invalid_argument for a non-positive cell count or order.
	 */
	std::size_t toroidalKnotCount( double cells, int order );

	/// The knots as a family's cut: how many surfaces, and the cut in `Psi_N`.
	struct ToroidalKnotCut
	{
		std::size_t surfaces = 0;
		double innerCut = 0.0;
		double outerCut = 0.0;
	};

	/**
	 * THE KNOTS AS A FAMILY: @a knots surfaces at `rho_j = j/( knots + 1 )`,
	 * `j = 1 .. knots`, so the cut is one knot spacing from the axis and one
	 * from the edge -- `meq::FluxFamilyOptions` with these, spaced in `rho`,
	 * labels them.
	 *
	 * **THE PART OF THE REQUEST NO KNOT SEES SHRINKS WITH THE MESH.** Inside
	 * the innermost knot and outside the outermost, `g^2` is the end cubic
	 * continued and the requested `q` is not imposed; that layer is one
	 * `Delta rho`, which the knot count makes proportional to `h`, so the
	 * unhonoured part of the request vanishes under refinement. A fixed cut
	 * does not: there the continuation is extrapolated over a width that the
	 * refining knots make ever more knot spacings wide (MEASUREMENTS.md
	 * M-392, M-394).
	 *
	 * @throws std::invalid_argument for fewer than four knots.
	 */
	ToroidalKnotCut toroidalKnotCut( std::size_t knots );

	/// The labels toroidalKnotCut() spaces, `Psi_N = ( j/( knots + 1 ) )^2`.
	std::vector<double> toroidalKnotLabels( std::size_t knots );

	/**
	 * `g^2` AS ITS VALUES AT THE KNOTS: A NOT-A-KNOT CUBIC IN THE SOURCE'S
	 * `Psi`, ITS END PIECES CONTINUED.
	 *
	 * The knots are given as the family's labels, `Psi_N` ascending, and the
	 * spline lives in `Psi = 1 - Psi_N` -- the reflection of section 2 happens
	 * here and in ggPrime(), nowhere else. Linear in the values, so the
	 * derivative with respect to value `j` is the cardinal function
	 * cardinal( j, . ).
	 *
	 * **IN `Psi`, NOT IN `rho`.** A smooth `g^2` is smooth in `Psi_N`; a spline
	 * in `rho = sqrt( Psi_N )` makes `gg' = d( g^2 )/dPsi` carry `1/rho` on
	 * the axis unless its slope there vanishes, which nothing imposes.
	 *
	 * **C^2, SO `gg'` IS C^1 AND `dF/dpsi` CONTINUOUS.** The source reads
	 * `gg' = ( g^2 )'/2`, the Newton matrix its derivative, and the bordered
	 * column one more (Profile::doublePrime()).
	 */
	class ToroidalFieldSpline
	{
		public:
			/// @throws std::invalid_argument for fewer than four labels, or
			///         labels that are not strictly ascending inside ( 0, 1 ).
			explicit ToroidalFieldSpline( std::vector<double> labels );

			std::size_t size() const { return knotLabels.size(); }

			/// The family's labels, `Psi_N`, ascending.
			std::vector<double> const &labels() const { return knotLabels; }

			/// `g^2` at each label, in the labels' order.
			/// @throws std::invalid_argument on a size mismatch or a
			///         non-finite value.
			void setValues( std::vector<double> values );
			std::vector<double> const &values() const { return knotValues; }

			/// Derivative @a derivative ( 0 to 3 ) of `g^2` with respect to the
			/// SOURCE's `Psi`, at @a psi; the end pieces continue outside the
			/// knots.
			double operator()( double psi, int derivative = 0 ) const;

			/// The same of the cardinal function of knot @a j: `d g^2/d u_j`.
			double cardinal( std::size_t j, double psi, int derivative = 0 ) const;

			/// `g g' = ( 1/2 ) d( g^2 )/dPsi` as a Profile for the source,
			/// exact on the spline, with its two derivatives.
			std::shared_ptr<Profile const> ggPrime() const;

			/// The same of cardinal function @a j: the direction in `gg'`
			/// that a change of value `j` makes.
			std::shared_ptr<Profile const> cardinalGGPrime( std::size_t j ) const;

			/// Whether `g^2 > 0` over all of `Psi` in [ 0, 1 ], sampled at
			/// @a samples points and at every knot: what the source will be
			/// evaluated over.
			bool positive( int samples = 256 ) const;

		private:
			std::vector<double> knotLabels, knotValues;
			std::vector<double> psi;               // source Psi, ascending
			std::vector<double> valuesToSeconds;   // ( n x n ), row-major
			std::vector<double> seconds;           // M_i at psi
			std::vector<double> sourceValues;      // values at psi

			double evaluate( std::vector<double> const &y, std::vector<double> const &m,
			                 double at, int derivative ) const;
	};

	/**
	 * THE ONE SCALAR BESIDE `p( Psi_N )` AND `q( Psi_N )` THAT SELECTS THE
	 * EQUILIBRIUM.
	 *
	 * Fixed-boundary equilibria sharing those two profiles form a ONE-PARAMETER
	 * FAMILY: `q` sets the field-line pitch and not the field's strength, which
	 * is VMEC's reason for wanting `PHIEDGE` beside `iota( s )`. The `q` rows
	 * leave that family the near-null direction of their Jacobian
	 * (MEASUREMENTS.md M-353), and each value here is a row that picks the
	 * member.
	 */
	enum class ScalePin
	{
		/// `g` on the boundary, `R0 B0`, read as the spline continued to
		/// `Psi = 0`: the row `( g^2( 0 ) - R0B0^2 )/R0B0^2`.
		EdgeField,

		/// `mu0 I_p` of the solved state, the row `( I - I* )/I*`. The inner
		/// solve's profile scale stays one.
		PlasmaCurrent,

		/// The DELIVERED toroidal flux, `Phi = int_plasma g( psi )/R dA` over
		/// the solved equilibrium with the iterate's own `g`, the row
		/// `( Phi - Phi* )/Phi*` -- the equilibrium's flux, not the target's
		/// (MEASUREMENTS.md M-362).
		ToroidalFlux
	};

	/// What one solve and extraction at a spline gives the system: per knot
	/// `q/g = V' < R^-2 >/4 pi^2`, and the pinned scalar as the solved state
	/// reads it (unused under ScalePin::EdgeField).
	struct SafetyFactorEvaluation
	{
		std::vector<double> qOverG;
		double measured = 0.0;
	};

	/// Their derivatives with respect to the knot values, THROUGH THE SOLVE:
	/// `d( q/g )_i/d u_j` row-major, knots by knots, at fixed `g` -- the
	/// geometry's response -- and the pinned scalar's total derivative.
	struct SafetyFactorDerivative
	{
		std::vector<double> qOverG;
		std::vector<double> measured;
	};

	/**
	 * THE q-DRIVEN SYSTEM: q AT EVERY KNOT AND ONE SCALE PIN, IN THE KNOT
	 * VALUES `u_j = g^2( knot j )`. MFEM-free: the solve and its derivative
	 * are callables (SafetyFactorTangent.hpp supplies them).
	 *
	 *     r_i     = q_i/q*_i - 1,   q_i = ( q/g )_i sqrt( u_i ),  i = 1 .. m
	 *     r_m+1   = the pin's row, ScalePin
	 *
	 * `m + 1` rows in `m` unknowns, solved by meq::solveLeastSquares with the
	 * pin row an EXACT constraint and the `q` rows in least squares within
	 * each level of the scale's coordinate (scaleGradient()). **NO q
	 * CONDITION GIVES WAY.** The `q` rows alone are nearly singular along the
	 * scale family, which the pin fixes; what is left is the discrete
	 * problem's inconsistency, the misfit, spread over the `q` rows and
	 * falling with `h` (MEASUREMENTS.md M-393).
	 * A square system -- one more value at the edge, held by the pin -- was
	 * measured and is ill-conditioned: the edge layer, which `q` does not
	 * see, takes the scale family's place and the pin barely reaches it
	 * (M-392).
	 *
	 * **THE PINS COMPOSE.** The answers under every pin lie on one curve --
	 * the q rows' least squares at each level of the scale coordinate -- and
	 * the pin, met exactly, picks the point; so an equilibrium solved under one
	 * pin is returned by another set to what it reads (QDrivenScalePin).
	 * Weighted into the least squares instead, a pin trades against the
	 * misfit, differently for each, and they returned answers 1.5e-6 of psi
	 * apart on an inconsistent target (M-396).
	 *
	 * **IT IS TOTAL.** A spline with `g^2 <= 0` anywhere on [ 0, 1 ], a solve
	 * that does not converge or an extraction that fails is a REFUSAL --
	 * residual() returns false -- and the least squares rejects the step.
	 * Nothing here throws for a reason the loop can meet.
	 *
	 * **THE TARGET IS ASKED IN THE SOURCE'S `Psi`**, one on the axis, and the
	 * knots are labelled by the family's `Psi_N`; the reflection is here.
	 */
	class SafetyFactorSystem
	{
		public:
			/// Solve and extract at this spline; false is a refusal.
			using Solve = std::function<bool( ToroidalFieldSpline const &, SafetyFactorEvaluation & )>;

			/// The derivative at the spline the last admissible Solve used;
			/// false if it cannot be formed.
			using Differentiate = std::function<bool( ToroidalFieldSpline const &, SafetyFactorDerivative & )>;

			struct Pin
			{
				ScalePin kind = ScalePin::EdgeField;

				/// `g_edge` (positive), `mu0 I_p` or `Phi`, by @ref kind.
				double value = 0.0;
			};

			/**
			 * @param labels the knots, `Psi_N` ascending.
			 * @param target `q` against the SOURCE's `Psi`.
			 * @param solveAccuracy the relative accuracy each equilibrium is
			 *        solved to, from which floor() is derived; zero claims an
			 *        exact solve.
			 * @throws std::invalid_argument on an unusable pin, a missing
			 *         callable, or a target that is zero or not finite at a
			 *         knot.
			 */
			SafetyFactorSystem( std::vector<double> labels,
			                    std::function<double( double psi )> const &target, Pin pin,
			                    Solve solve, Differentiate differentiate,
			                    double solveAccuracy );

			std::size_t unknowns() const { return spline.size(); }
			std::size_t rows() const { return spline.size() + 1; }

			/**
			 * THE SCALE'S COORDINATE, `sum_j ln u_j`, as its gradient `1/u_j`:
			 * the level functional meq::solveLeastSquares meets the q rows
			 * within, the pin row being its one constraint. The q rows' near-
			 * null direction is the scale family (M-353, M-396), and this
			 * moves along it; holding it, the q rows are a well-posed least
			 * squares, and every pin chooses from the one curve that traces.
			 */
			static std::vector<double> scaleGradient( std::vector<double> const &u );

			/// The residual at @a u, `rows()` entries; false is a refusal.
			bool residual( std::vector<double> const &u, std::vector<double> &r );

			/// The Jacobian, `rows() x unknowns()` row-major, at @a u, which
			/// must be the last point residual() admitted; false if it cannot
			/// be formed.
			bool jacobian( std::vector<double> const &u, std::vector<double> &j );

			/**
			 * THE ROUND-OFF FLOOR OF EACH ROW: what a residual
			 * evaluated twice from two warm starts differs by, modelled.
			 * `q = ( q/g ) sqrt( u )` carries `V'` and `< R^-2 >`, each in
			 * error by about the solve's relative accuracy, so `2 eps_solve`;
			 * four standard deviations, so `8 eps_solve`. The pin rows of the
			 * current and the flux carry one quadrature of the solved state,
			 * `4 eps_solve`; R0B0's row reads no solve. The per-average factor
			 * is the measured piece (CLAUDE_INVERSION.md, *Converged means*),
			 * and where it is under the noise the least squares measures the
			 * noise (M-363). One entry per row, the pin row's weighted.
			 */
			std::vector<double> floor() const;

			/// The spline at the last admissible point.
			ToroidalFieldSpline const &answer() const { return spline; }

			/// `q*` at the knots, the target as the rows read it.
			std::vector<double> const &targets() const { return target; }

			/// The pinned scalar at the last admissible point; `g_edge` under
			/// ScalePin::EdgeField, from the spline.
			double lastMeasured() const { return measuredValue; }

			/// Solves asked for, and refusals among them (including the
			/// splines refused before any solve).
			int solves() const { return solveCount; }
			int refusals() const { return refusalCount; }

		private:
			ToroidalFieldSpline spline;
			std::vector<double> target;
			Pin pin;
			Solve solve;
			Differentiate differentiate;
			double accuracy = 0.0;
			std::vector<double> lastGood;
			std::vector<double> qOverG;
			double measuredValue = 0.0;
			int solveCount = 0;
			int refusalCount = 0;
	};

	/**
	 * Fritsch-Carlson limited slopes for the data `( x, y )`.
	 *
	 * Exposed because it is the one piece of numerics here that is worth
	 * testing on its own: it is what stops the fed-back profile oscillating,
	 * and its defining property -- no overshoot, and a cubic reproduced where
	 * the data is smooth -- is a statement about the RULE rather than about any
	 * equilibrium.
	 *
	 * @throws std::invalid_argument unless the sizes match, there are at least
	 *         two points, and @a x is strictly ascending.
	 */
	std::vector<double> monotoneSlopes( std::vector<double> const &x,
	                                    std::vector<double> const &y );

}

#endif
