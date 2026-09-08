#ifndef MEQ_SAFETY_FACTOR_HPP
#define MEQ_SAFETY_FACTOR_HPP

#include <functional>
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
 * SOLUTION, and the solution depends on `g`. So the loop is
 *
 *     solve -> extract the surfaces -> invert -> rebuild gg' -> solve
 *
 * and the fixed point is an equilibrium whose own `q` is the one asked for.
 * That is a Picard iteration on the geometry, and it is deliberately the
 * halfway house rather than a Newton: making it quadratic needs
 * `d( V', < R^-2 > )/d psi`, which `INVERSION-PLAN.md` section 11.1 measures at
 * about **5.7 hours per Jacobian** by differencing, since it is `nFieldDOF`
 * complete re-extractions of the whole surface family. The shape derivative is
 * what would make that affordable and it is not built. See section 5.
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
 * has to say something where the family does not. `meq::SplineProfile` CLAMPS
 * outside its knots -- deliberately, since a linear extrapolation of a steep
 * edge profile turns an overshoot into a NaN -- so the value carried past the
 * ends is the end knot's.
 *
 * **THAT IS A CHOICE ABOUT THE PHYSICS AND IT IS WRONG AT THE EDGE.** Outside
 * the plasma the toroidal field is the vacuum one, `g = const`, so `gg' = 0`;
 * clamping at the outermost knot's value instead keeps `g^2` growing. Whether
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
		 * Carry the outermost knot's value, which is what SplineProfile does
		 * with no help. Right when the source is confined to the plasma, where
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
	 * `g^2` as a polynomial in the SOURCE's `Psi`, by least squares.
	 *
	 * Returns the coefficients `c` of `g^2 = sum_j c_j Psi^j`, ascending, so
	 * `gg' = ( 1/2 ) sum_j j c_j Psi^(j-1)` is exact and analytic rather than a
	 * slope rule applied to a point cloud.
	 *
	 * **THIS IS WHAT MAKES THE OUTER LOOP RUN, AND THE REASON IS THE
	 * CONDITIONING OF THE INVERSION AND NOT A TASTE FOR SMOOTH PROFILES.**
	 * `g = 4 pi^2 q/( V' < R^-2 > )` divides by extracted geometry, and near the
	 * axis a small change in the equilibrium moves `V'` by much more than it
	 * moves `g` -- measured on the loop of `SafetyFactorSolver`, the innermost
	 * surface's inverted `g` moved **13%** for a **1.3%** change in the profile
	 * that produced it. That is exactly where the extraction is least reliable,
	 * which is why the inner cut exists at all.
	 *
	 * Interpolating the cloud and differentiating it turns that into SHAPE: the
	 * loop's second iterate came back non-monotone in `g`, so `gg'` changed
	 * sign, and the next solve did not converge. A low-order fit cannot express
	 * that, which is the point -- it is the same reason a production code
	 * carries its profiles in a basis rather than as a table it differentiates.
	 *
	 * @param degree of the polynomial in `Psi`. Two or three is the usual; the
	 *        fit is refused if the family has fewer surfaces than coefficients.
	 * @throws std::invalid_argument if @a degree is zero, if there are fewer
	 *         surfaces than `degree + 1`, or if the labels are not ascending.
	 */
	std::vector<double> fitToroidalFieldSquared( ToroidalField const &field,
	                                             unsigned int degree );

	/// `gg'` knots from that fit, in the SOURCE's `Psi`, exact on the
	/// polynomial. Same reflection and the same edge extension as
	/// ggPrimeKnots(); @a samples knots are laid down uniformly over the
	/// family's own range.
	std::vector<Knot> ggPrimeKnotsFromFit( ToroidalField const &field,
	                                       unsigned int degree,
	                                       std::size_t samples = 17,
	                                       EdgeExtension extension
	                                           = EdgeExtension::Clamp );

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
