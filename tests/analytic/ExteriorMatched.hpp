#ifndef MEQ_TESTS_EXTERIORMATCHED_HPP
#define MEQ_TESTS_EXTERIORMATCHED_HPP

/*
 * A manufactured equilibrium that is EXACTLY ONE EXTERIOR MODE outside a
 * radius rho_0 and a polynomial inside it, for FB-1.
 *
 * FREE-BOUNDARY-PLAN.md section 8 makes FB-1 "the stage to protect" -- it
 * exercises ExteriorDtN, the transferred datum with a non-zero g, the
 * transmission condition, the augmented solve and the coil sources all at once
 * -- and gives its acceptance as "psi_h against that closed form at k+1", the
 * closed form being "the sum of the coils' loop fields".
 *
 *
 * WHY THE LOOP FIELDS CANNOT BE THAT CLOSED FORM
 * ----------------------------------------------
 *
 * They are the right check on ExteriorDtN alone -- tests/unit/ExteriorDtNTests
 * puts a current loop through the symbol and reads 6.9e-14, and section 3.3 is
 * explicit that a field the class knows nothing about is the only kind that can
 * catch a self-consistent misreading. What they cannot be is the exact
 * SOLUTION of an order study, and there are two independent reasons, matching
 * the two coil models section 5.4 offers.
 *
 * A FILAMENT IS A DELTA FUNCTION. In the ( R, z ) half-plane a circular loop is
 * a point, so its source is a Dirac mass and its psi has a logarithmic
 * singularity there -- CurrentLoop.hpp measures exactly that, psi reading
 * 3.4956, 5.7962 and 8.0828 at distances 1e-3, 1e-5 and 1e-7 against
 * -( 1/2 ) ln eps of 3.4539, 5.7565 and 8.0590, with the gradient going as
 * 1/eps. And ln( d ) IS NOT IN H^1 IN TWO DIMENSIONS: |grad psi| ~ 1/d makes
 * the energy integral 2 pi times the integral of dd/d, which diverges
 * logarithmically. A finite element solution cannot converge at k+1 to a
 * function whose energy is infinite; it converges at a reduced rate set by the
 * singularity, and a rate table would then be measuring the filament rather
 * than the coupling. That is the same species of obstruction CLAUDE.md records
 * for the re-entrant corner in *Testing stance* -- "no polynomial degree
 * recovers it" -- and the remedy is the same one: do not pose the study there.
 *
 * A FINITE CROSS-SECTION FIXES THE REGULARITY AND LOSES THE CLOSED FORM.
 * Section 5.4's other model, F_coil = mu0 R I_k / |Omega_ck| on a coil
 * subdomain, has psi in H^2 -- but its exact field is the two-dimensional
 * integral of loop fields over the cross-section, which is semi-analytic
 * (quadrature, not a formula) and whose integrand is weakly singular for a
 * field point inside the coil. A reference you have to converge is not a closed
 * form, and CLAUDE.md's testing stance is that a self-convergence study cannot
 * see a wrong sign convention. WORSE, THE TOP HAT CAPS THE RATE ANYWAY: a
 * source with a jump at the coil edge gives psi in C^1 there and no better, so
 * psi is in H^{5/2-eps} and the rate is capped near 2.5 whatever k is. That cap
 * is precisely what the contact order below exists to remove.
 *
 * So FB-1 needs a manufactured solution, and this is it. It is the coil
 * subdomain done properly: a source compactly supported strictly inside Gamma,
 * with a SMOOTH radial profile instead of a top hat, chosen so that the field
 * it produces is elementary.
 *
 *
 * THE CONSTRUCTION
 * ----------------
 *
 * In spherical coordinates about a centre on the axis -- R = rho sin theta,
 * z = zCentre + rho cos theta, mu = cos theta -- FREE-BOUNDARY-PLAN.md
 * section 3 separates the operator as
 *
 *     Delta*( F( rho ) C_n( mu ) ) = [ F'' - n( n - 1 ) F / rho^2 ] C_n( mu ),
 *
 * C_n being the Gegenbauer function of order -1/2 that src/meq/ExteriorDtN
 * already implements. Re-verified here symbolically in ( R, z ), exactly, for
 * n = 2..7 and exponents n, 1-n, n+2, n+4, n+6; the two homogeneous solutions
 * are rho^n (regular at the centre) and rho^{1-n} (decaying).
 *
 * TAKE psi TO BE THE DECAYING MODE OUTSIDE rho_0 AND A POLYNOMIAL INSIDE:
 *
 *     F( rho ) = rho^{1-n}                                    rho >= rho_0
 *     F( rho ) = A rho^n + sum_{i=0}^{p} c_i rho^{n+2i+2}     rho <  rho_0
 *
 * Every interior term is regular at the centre, because rho^{n+2j} C_n( mu ) is
 * a POLYNOMIAL in ( R, z ) -- R^2/2, R^2 z/2, R^2( 4z^2 - R^2 )/8 at
 * n = 2, 3, 4, which are Cerfon & Freidberg's own basis functions and are what
 * VacuumHarmonic.hpp carries. Every term also vanishes identically on the axis,
 * since each carries R^2, so the flat side of the half-disc is satisfied by
 * construction rather than imposed.
 *
 * Applying the operator to the interior form and using
 * Delta*( rho^{n+2j} C_n ) = 2j( 2n + 2j - 1 ) rho^{n+2j-2} C_n gives a source
 * that is a polynomial in rho^2 times rho^n, and the c_i are chosen to make it
 *
 *     Delta* psi = kappa rho^n ( rho_0^2 - rho^2 )^p C_n( mu )  rho < rho_0
 *     Delta* psi = 0                                        rho >= rho_0
 *
 * so that A and kappa are left to match VALUE and SLOPE at rho_0. Solving,
 * with binom the binomial coefficient and
 *
 *     D  = sum_{i=0}^{p} (-1)^i binom( p, i ) / ( 2n + 2i + 1 )
 *     S0 = sum_{i=0}^{p} (-1)^i binom( p, i )
 *                        / ( 2( i + 1 )( 2n + 2i + 1 ) ),
 *
 *     kappa = ( 1 - 2n ) rho_0^{-( 2n + 2p + 1 )} / D
 *     A     = rho_0^{1-2n} - kappa rho_0^{2p+2} S0
 *     c_i   = kappa (-1)^i binom( p, i ) rho_0^{2( p - i )}
 *             / ( 2( i + 1 )( 2n + 2i + 1 ) ).
 *
 * D IS NEVER ZERO, so kappa always exists: D is exactly the integral of
 * x^{2n}( 1 - x^2 )^p over [ 0, 1 ], which is positive for every admissible
 * n >= 2 and p >= 0. Checked symbolically as well as argued.
 *
 * THE SOURCE IS THEN IDENTICALLY ZERO OUTSIDE rho_0, EXACTLY, and that is the
 * whole point. psi decays at infinity, the source is compactly supported
 * strictly inside Gamma, and on ANY semicircle of radius rho_Gamma >= rho_0 the
 * exterior expansion is one mode with a known coefficient -- so the test can
 * pin the recovered coefficient VECTOR and not merely the field, with every
 * entry but one exactly zero.
 *
 *
 * THE CONTACT ORDER p, AND WHY IT IS NOT 0
 * ----------------------------------------
 *
 * p = 0 IS THE FAMILY THE BRIEF PROPOSED -- A rho^n + B rho^{n+2}, two
 * constants for two matching conditions -- and it does match, with a BOUNDED
 * source. Solving it by hand gives, and this file's own closed form reproduces
 * exactly at n = 2..7,
 *
 *     A = ( 2n + 1 )/2 rho_0^{1-2n},      B = -( 2n - 1 )/2 rho_0^{-( 2n+1 )},
 *     F_source = ( 4n^2 - 1 ) rho_0^{-( 2n+1 )} rho^n C_n( mu ),  rho < rho_0.
 *
 * BUT THE SOURCE IS THEN DISCONTINUOUS AT rho_0 -- inside it approaches
 * ( 4n^2 - 1 ) rho_0^{-( n+1 )} C_n and outside it is zero -- so psi is C^1
 * and no better, psi is in H^{5/2-eps}, AND THE RATE IS CAPPED NEAR 2.5
 * HOWEVER LARGE k IS. That is the top hat's defect over again, arrived at
 * from the other direction, and it would make a k = 3 study of the coupling
 * measure the matching radius instead.
 *
 * Requiring the source to vanish to order p at rho_0 buys smoothness back one
 * derivative at a time. Verified symbolically at n = 2, 3, 5 and p = 0..4: the
 * first jump is in d^{p+2}F, so psi is exactly C^{p+1} across the sphere and
 * lies in H^{p+5/2-eps}.
 *
 *     p     psi is       psi in         highest k with H^{k+2} available
 *     0     C^1          H^{5/2}        (0)
 *     1     C^2          H^{7/2}        1
 *     2     C^3          H^{9/2}        2
 *     3     C^4          H^{11/2}       3
 *     4     C^5          H^{13/2}       4
 *
 * THE DEFAULT IS p = 4, which covers the k = 1..4 that CLAUDE.md's testing
 * stance runs, with a margin of half a derivative rather than none. It also
 * makes deltaStarFD() clean across rho_0: a central second difference needs
 * four continuous derivatives for its O( h^2 ) truncation, and p = 4 supplies
 * five. p = 0 is kept reachable, as the control that shows what the smoothing
 * is worth, and measurement 1 below is sharper than expected: at p = 0 the
 * residual at the interface DOES NOT CONVERGE IN h AT ALL, while inside and
 * outside it are untouched.
 *
 *
 * THE COEFFICIENT SCALING, WHICH IS THE EASIEST THING HERE TO GET WRONG
 * --------------------------------------------------------------------
 *
 * ExteriorDtN::exterior() expands in ( rho/rho_Gamma )^{1-n}, NOT in rho^{1-n},
 * so the coefficient it works in is psi ON GAMMA and depends on where Gamma is.
 * With psi = sum_n alpha_n rho^{1-n} C_n( mu ) outside rho_0,
 *
 *     a_n = alpha_n rho_Gamma^{1-n},
 *
 * and exteriorCoefficients() returns exactly that. An amplitude read straight
 * into a_n instead would pass at rho_Gamma = 1 and fail everywhere else, which
 * is why both radii below are measured and why the test uses rho_Gamma = 1.5
 * rho_0 and 3 rho_0 rather than one of them.
 *
 *
 * WHAT WAS MEASURED
 * -----------------
 *
 * Everything below is on multiMode() -- degrees 2, 3 and 5 at amplitudes 1.0,
 * -0.6 and 0.35, zCentre = 0, rho_0 = 1, p = 4 -- unless another p is named,
 * with p = 0 kept as the control. DEGREE 4 IS DELIBERATELY ABSENT so that the
 * coefficient vector has an exact zero BETWEEN two nonzeros, where an
 * off-by-one in the degree indexing would show and where a single-mode case
 * could not see one. Reference scale |psi( 0.7, 0.2 )| = 4.129582e-01. The
 * grid is R in [ 0.05, 2.40 ], |z| <= 2.40 at spacing 0.05 -- 4656 points,
 * spanning both sides of rho_0 and out to 5.8 rho_0.
 *
 *  0. THE WHOLE CONSTRUCTION AGAINST AN INDEPENDENT ONE, WHICH IS THE CHECK
 *     THE OTHERS ARE ARRANGED AROUND. A sympy script solved the matching
 *     problem afresh from the two conditions -- not from the closed forms
 *     above -- differentiated the result in ( R, z ), and formed Delta*
 *     directly. It reports Delta*( psi ) + f = 0 SYMBOLICALLY, exactly zero
 *     rather than to a tolerance, for p = 0 and p = 4 alike, and against this
 *     header's evaluation at ten scattered points:
 *
 *         p = 0    psi 2.2e-15   d_r 2.9e-15   d_z 8.3e-16   f 1.2e-14
 *         p = 4    psi 1.2e-14   d_r 2.3e-13   d_z 5.8e-14   f 1.9e-13
 *
 *     The closed forms themselves were checked the same way, against a
 *     from-scratch symbolic solve at n = 2..8 and p = 0..5: A, kappa and every
 *     c_i EXACT. So is the separation identity, at n = 2..7 over exponents n,
 *     1-n, n+2, n+4 and n+6, and so is the order of contact -- the first jump
 *     is in d^{p+2}F at n = 2, 3 and 5 and p = 0..4, which is the C^{p+1}
 *     claim.
 *
 *     THE SMOOTHING COSTS ABOUT TWO DIGITS and that is worth knowing before
 *     reading any of the numbers below. p = 4 makes the interior a degree
 *     n + 2p + 2 = 12 polynomial whose coefficients cancel to an O( 1 ) answer,
 *     so p = 4 sits at 1e-13 where p = 0 sits at 1e-15. It buys four
 *     derivatives for two digits, which at these magnitudes is the right trade
 *     and is not free.
 *
 *  1. Delta* psi = -f, by deltaStarFD(), split by where the stencil sits
 *     relative to rho_0 -- worst |Delta* psi + f| over the grid:
 *
 *         p = 4    h        inside      band         outside
 *                  1e-2     1.694e-02   5.141e-03    1.889e-03
 *                  1e-3     1.695e-04   5.137e-05    1.888e-05
 *                  1e-4     1.717e-06   1.617e-06    2.345e-07
 *                  1e-5     7.173e-05   1.344e-04    4.553e-06
 *                  Rich     3.528e-08   4.745e-08    2.548e-09
 *
 *     Clean O( h^2 ) to h = 1e-4 -- three decades, three factors of a hundred
 *     -- then round-off, and Richardson from h = 1e-3 sees past the floor to
 *     3.5e-08, which is 8.5e-08 of the reference. THE FLOOR IS THE INSTRUMENT
 *     AND NOT THE FIELD, as it is everywhere else in this tree that a
 *     derivative meets a difference; item 0 is the statement without an
 *     instrument in it.
 *
 *     THE p = 0 CONTROL IS THE INTERESTING COLUMN, and it does not merely
 *     degrade -- IT DOES NOT CONVERGE AT ALL:
 *
 *         p = 0    h        inside      band         outside
 *                  1e-2     9.576e-04   4.710e+00    1.889e-03
 *                  1e-3     9.576e-06   4.805e+00    1.888e-05
 *                  1e-4     1.224e-07   4.814e+00    2.345e-07
 *
 *     The band sits at 4.81 whatever h is, because there Delta* psi is
 *     genuinely discontinuous and a stencil straddling it is differencing two
 *     different functions. Inside and outside are unaffected, and are in fact
 *     slightly BETTER than p = 4's for the cancellation reason in item 0. So
 *     the contact order is doing exactly one thing, in exactly one place, and
 *     the measurement can see which.
 *
 *  2. f is IDENTICALLY zero for rho >= rho_0 -- the literal 0.0, not a small
 *     number -- at 861 points per fixture running from rho_0 out to 6 rho_0,
 *     rho_0 itself included, for p = 0 and p = 4 alike. Nonzero count: 0 and 0.
 *     For contrast the same field reads -7.02e-01 at rho = 0.5 rho_0.
 *
 *     A TRAP FOUND WHILE MEASURING IT, AND IT IS THE READER'S TRAP AND NOT THE
 *     FIXTURE'S. A test point built as ( rho sin th, rho cos th ) does not
 *     round-trip: hypot of the pair can come back one ulp BELOW rho, which
 *     puts it inside the source. At p >= 1 that costs nothing, since the source
 *     vanishes there -- p = 4 leaked 1.1e-60 -- but AT p = 0 THE SAME POINT
 *     RETURNS -2.6e-01, the full discontinuity. One point in 861 hit it. Nudge
 *     the radius until sphericalRadius() agrees, which is what the figures above do.
 *
 *  3. psi is C^1 across rho_0, checked on the interior POLYNOMIAL rather than
 *     by stepping either side of the branch -- evaluating one ulp inside rho_0
 *     and comparing against rho_0^{1-n} and ( 1-n ) rho_0^{-n}, worst over the
 *     three modes:
 *
 *         p     F           F'
 *         0     2.220e-16   2.220e-15
 *         1     1.554e-15   3.664e-15
 *         2     7.661e-15   7.550e-15
 *         3     2.898e-14   7.083e-14
 *         4     8.460e-14   5.689e-13
 *
 *     -- item 0's two digits again, and nothing else. AND THE SECOND
 *     DERIVATIVE IS WHERE THE TWO CASES PART: differenced three steps either
 *     side at n = 2, p = 0 reads -12.8922 inside against 1.9821 outside, a
 *     jump of -14.87 against the -15/rho_0^3 the symbolic calculation predicts
 *     exactly, while p = 4 reads 2.0181 against 1.9821 -- continuous to the
 *     stencil's own truncation.
 *
 *  4. Analytic derivatives against central differences of psi(), worst
 *     relative to |d_r psi| + |d_z psi| over the grid:
 *
 *         h = 1e-2    d_r 2.480e-03    d_z 7.729e-04
 *         h = 1e-3        2.481e-05        7.732e-06
 *         h = 1e-4        2.481e-07        7.732e-08
 *         h = 1e-5        2.493e-09        9.240e-10
 *         h = 1e-6        8.490e-09        8.765e-09    round-off, rising
 *
 *     Richardson, ( 4 D( h/2 ) - D( h ) )/3 from h = 1e-3, reads 2.797e-11 and
 *     2.203e-11 -- below anything the plain column reaches, which is
 *     src/meq/Zernike.hpp's finding once more.
 *
 *  5. THE EXTERIOR COEFFICIENTS, which is the sharpest check here. psi on
 *     Gamma handed to ExteriorDtN::coefficients() with 12 modes, against
 *     exteriorCoefficients():
 *
 *         case          rho_Gamma    live entries     zero entries
 *         single n = 2  1.5 rho_0    0.000e+00 rel    11 of 12, worst 3.271e-15
 *         single n = 2  3.0 rho_0    0.000e+00 rel    11 of 12, worst 1.636e-15
 *         multi         1.5 rho_0    2.007e-16 rel     9 of 12, worst 3.353e-15
 *         multi         3.0 rho_0    6.825e-15 rel     9 of 12, worst 1.634e-15
 *
 *     The single-mode rows are EXACT -- the recovered coefficient is bit for
 *     bit the amplitude divided by rho_Gamma^{n-1} -- and degree 4 is among the
 *     zeros in the multi rows, sitting between two live entries. TWO RADII ARE
 *     RUN BECAUSE ONE CANNOT SEE THE SCALING: a_n = alpha_n rho_Gamma^{1-n}
 *     differs from alpha_n by 1.5^{-1} and 3^{-1} at n = 2 and by 5.06 and 81
 *     at n = 5, so a fixture returning the bare amplitude would fail both rows
 *     by a factor rather than passing one of them.
 *
 *     The zero entries grow with the degree -- 1e-17 at n = 3, 3e-15 at n = 12
 *     -- because the projection divides by h_n = 2/( n( n-1 )( 2n-1 ) ), which
 *     is 1e-3 at n = 12. That is the quadrature's round-off amplified by the
 *     mass and is a property of ExteriorDtN, not of this field.
 *
 *  6. THE DtN. ExteriorDtN::symbol() applied to those coefficients and summed
 *     against the basis, versus this fixture's analytic dPsiDrho(), at 39
 *     angles on Gamma:
 *
 *         single n = 2   4.905e-16 rel at both radii
 *         multi          6.341e-16 and 5.033e-16 rel
 *
 *     and dPsiDrho() against gradPsi() projected on the radial direction,
 *     worst 1.296e-15 -- so the two are independent paths and the DtN figure
 *     is not being compared against itself.
 *
 *  7. THE AXIS, WHICH FB-1'S GEOMETRY MAKES A BOUNDARY AND NOT AN INTERIOR
 *     LINE. psi( 0, z ) is 0.0 EXACTLY -- the literal, not a tolerance -- at
 *     121 values of z from -6 to +6, and gradPsi( 0, z ) is exactly ( 0, 0 ).
 *     That is inherited rather than arranged: every C_n carries the factor
 *     ( 1 - mu )( 1 + mu ), which is exactly zero at mu = +-1, so the flat
 *     side of the half-disc is satisfied by construction. It is the same
 *     property VacuumHarmonic.hpp was written for, arriving here for free.
 *
 *     The flux has a finite limit there and the expression does not reach it,
 *     exactly as CurrentLoop.hpp's does: psi( R, 0.3 )/R^2 settles to
 *     2.772735 as R falls through 1e-2, 1e-3, 1e-4, and q_r to 5.5454705 at
 *     R = 1e-3, 1e-5, 1e-7 while q_z goes as R. flux( 0, z ) is NaN in both
 *     components.
 *
 * WHAT IS NOT CHECKED HERE. Nothing in this file has been through a solver:
 * these are properties of the closed form, and whether MEQ's coupled system
 * reproduces it at k+1 is FB-1's own measurement and the reason the fixture
 * exists. The Sobolev exponents in the p table are the standard counting from
 * the measured order of contact, not themselves measured -- what is measured
 * is that the first discontinuous derivative is the ( p + 2 )-th. And the
 * claim that a filament caps the rate is an argument from |grad psi| ~ 1/d,
 * which CurrentLoop.hpp measures; no reduced-rate study was run to confirm it,
 * because running one is the thing this fixture exists to avoid.
 */

#include "meq/ExteriorDtN.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace meq
{
namespace analytic
{

/// A flux function that is exactly a sum of decaying exterior modes outside a
/// radius rho_0, polynomial inside it, and driven by a source compactly
/// supported in rho < rho_0.
///
/// MODES ARE INDEXED BY THEIR DEGREE n >= 2, as in meq::ExteriorDtN and for the
/// same reason: degrees 0 and 1 are inadmissible rather than omitted, and
/// getting it wrong is an off-by-two, which produces a plausible wrong answer
/// rather than an obvious one.
///
/// psi() and the gradient are valid everywhere, the centre included. flux()
/// divides by R and so is NaN on the axis, exactly as CurrentLoop.hpp's is;
/// the limit exists, since every mode carries R^2, but the expression reaching
/// it does not. A caller who needs the axis wants VacuumHarmonic.hpp.
class ExteriorMatched
{
	public:
		/// One exterior mode: its degree and its amplitude in the expansion
		/// psi = sum alpha_n rho^{1-n} C_n( mu ) valid for rho >= rho_0.
		///
		/// The amplitude is the coefficient of rho^{1-n}, NOT of
		/// ( rho/rho_Gamma )^{1-n}, which is what ExteriorDtN works in. See
		/// exteriorCoefficients() for the conversion and the file comment for
		/// why it is the easiest thing here to get wrong.
		struct Mode
		{
			int degree;
			double amplitude;
		};

		/// The order to which the source vanishes at rho_0, if none is given.
		/// Four, so that psi is C^5 and a k+1 study is unobstructed to k = 4.
		static constexpr int defaultContactOrder = 4;

		/// @param zCentreIn        the axial position of the centre, which is
		///                         the centre of FB-1's semicircle Gamma.
		/// @param rho0In           the source radius. Strictly positive, and
		///                         strictly less than any rho_Gamma the field
		///                         is later expanded on.
		/// @param modesIn          the modes, by degree and amplitude. At
		///                         least one, every degree at least 2, no
		///                         degree repeated.
		/// @param contactOrderIn   p, the order to which the source vanishes
		///                         at rho_0. Zero is allowed and is the brief's
		///                         own two-term family; see the file comment
		///                         for what it costs.
		///
		/// Throws std::invalid_argument rather than producing a field that is
		/// quietly not what was asked for.
		ExteriorMatched( double zCentreIn, double rho0In,
		                 std::vector<Mode> const &modesIn,
		                 int contactOrderIn = defaultContactOrder )
			: zCentreValue( zCentreIn ),
			  rho0Value( rho0In ),
			  contactOrderValue( contactOrderIn ),
			  angular( zCentreIn, rho0In,
			           checkedHighestDegree( rho0In, modesIn, contactOrderIn )
			           - ExteriorDtN::firstMode() + 1 )
		{
			buildRadialCoefficients( modesIn );
		}

		/// One mode, the dipole n = 2, at unit amplitude on the midplane with
		/// rho_0 = 1.
		///
		/// n = 2 IS THE FAR FIELD OF A CURRENT LOOP: rho^{-1} C_2( mu ) is
		/// R^2/( 2 rho^3 ), which is the dipole CurrentLoop.hpp's file comment
		/// records its field decaying to. So the simplest case here is the
		/// leading term of the field FB-1 was originally to be measured
		/// against, with the singularity replaced by a smooth source.
		static ExteriorMatched singleMode()
		{
			return ExteriorMatched( 0.0, 1.0, { { 2, 1.0 } } );
		}

		/// Three modes of both parities with degree 4 SKIPPED, which is the
		/// case to run first.
		///
		/// A single mode is a special case that can hide an indexing error --
		/// every entry of the coefficient vector but one is zero, and so is
		/// every entry of a vector that was built wrongly. Here degrees 2, 3
		/// and 5 are live and degree 4 is not, so a correct answer has a zero
		/// BETWEEN two nonzeros and an off-by-one cannot reproduce it. n = 3
		/// is odd in z, so the field is up-down asymmetric as well.
		static ExteriorMatched multiMode()
		{
			return ExteriorMatched( 0.0, 1.0,
			                        { { 2, 1.0 }, { 3, -0.6 }, { 5, 0.35 } } );
		}

		/// The poloidal flux. Valid everywhere, the centre and the axis
		/// included, where it is zero.
		double psi( double radius, double z ) const
		{
			double const rho = sphericalRadius( radius, z );

			double total = 0.0;
			for ( RadialMode const &mode : modeValues )
			{
				total += mode.amplitude*radialAt( mode, rho )
				         *angular.basis( mode.degree, radius, z );
			}
			return total;
		}

		/// d psi / d R. ANALYTIC, not a difference.
		///
		/// With psi = F( rho ) C_n( mu ), rho = hypot( R, z - zCentre ) and
		/// mu = ( z - zCentre )/rho, the chain rule needs d rho/d R = R/rho and
		/// d mu/d R = -mu R/rho^2, giving
		///
		///     d psi/d R = F'( rho )( R/rho ) C_n - F( rho ) C_n'( mu ) mu R/rho^2.
		double dPsiDr( double radius, double z ) const
		{
			double dR = 0.0;
			double dZ = 0.0;
			gradPsi( radius, z, dR, dZ );
			return dR;
		}

		/// d psi / d z. Analytic, by the same route, with d rho/d z = mu and
		/// d mu/d z = ( 1 - mu^2 )/rho -- WRITTEN AS R^2/rho^3, which is the
		/// same number with no cancellation in it. Near the axis mu is close to
		/// +-1 and 1 - mu*mu loses digits there; ExteriorDtN.cpp records eight
		/// orders on exactly that point, and R^2/rho^3 sidesteps it entirely.
		double dPsiDz( double radius, double z ) const
		{
			double dR = 0.0;
			double dZ = 0.0;
			gradPsi( radius, z, dR, dZ );
			return dZ;
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ). Not the HDG flux: that is
		/// this divided by R, see flux().
		///
		/// EXACT AT THE CENTRE, where the expression is 0/0. Every mode is
		/// A rho^n C_n + higher and rho^n C_n is a homogeneous polynomial of
		/// degree n >= 2, so the gradient there is zero and the branch returns
		/// it rather than a NaN.
		void gradPsi( double radius, double z, double &dR, double &dZ ) const
		{
			double const dz = z - zCentreValue;
			double const rho = std::hypot( radius, dz );

			if ( rho == 0.0 )
			{
				dR = 0.0;
				dZ = 0.0;
				return;
			}

			double const mu = dz/rho;
			double totalR = 0.0;
			double totalZ = 0.0;

			for ( RadialMode const &mode : modeValues )
			{
				double const f = radialAt( mode, rho );
				double const fPrime = radialDerivativeAt( mode, rho );
				double const c = angular.basis( mode.degree, radius, z );
				double const cPrime = angular.basisDerivative( mode.degree,
				                                               radius, z );

				totalR += mode.amplitude
				          *( fPrime*( radius/rho )*c - f*cPrime*mu*radius/( rho*rho ) );
				totalZ += mode.amplitude
				          *( fPrime*mu*c + f*cPrime*radius*radius/( rho*rho*rho ) );
			}

			dR = totalR;
			dZ = totalZ;
		}

		/// The HDG flux q = grad_bar( psi ) / R, in the convergence harness's
		/// own signature.
		///
		/// NaN at R = 0, like CurrentLoop.hpp and unlike VacuumHarmonic.hpp.
		/// The limit exists -- every mode carries R^2, so d_r psi ~ R -- but
		/// this expression does not reach it, and a fixture that quietly
		/// switched formula near the axis would be a fixture that can disagree
		/// with itself.
		void flux( double radius, double z, double &qR, double &qZ ) const
		{
			gradPsi( radius, z, qR, qZ );
			qR /= radius;
			qZ /= radius;
		}

		/// d psi / d rho at fixed direction: the OUTWARD radial derivative
		/// about the centre, which is what the Dirichlet-to-Neumann map
		/// returns and so is what FB-1's transmission condition is compared
		/// against.
		///
		/// Computed as sum alpha_n F_n'( rho ) C_n( mu ) rather than by
		/// projecting gradPsi() on the radial direction, so that the two are
		/// independent paths to the same number and a check between them means
		/// something. Measured, they agree to 1.3e-15.
		double dPsiDrho( double radius, double z ) const
		{
			double const rho = sphericalRadius( radius, z );
			if ( rho == 0.0 )
			{
				return 0.0;
			}

			double total = 0.0;
			for ( RadialMode const &mode : modeValues )
			{
				total += mode.amplitude*radialDerivativeAt( mode, rho )
				         *angular.basis( mode.degree, radius, z );
			}
			return total;
		}

		/// The source F, in MEQ's own convention: -Delta* psi = F, which is the
		/// F of eq (2) and NOT F/R. VacuumHarmonic.hpp and
		/// ManufacturedNonlinear.hpp document the same convention; the weak
		/// form's right-hand side is ( F/R, w ) and the division happens there.
		///
		/// IDENTICALLY ZERO FOR rho >= rho_0 -- a bare literal, not a small
		/// number -- which is the property the whole fixture exists for. At
		/// exactly rho = rho_0 the outside branch is taken; for every p >= 1
		/// the two agree there anyway, and at p = 0 the source is genuinely
		/// discontinuous and a value has to be chosen.
		double f( double radius, double z, double /*psi*/ ) const
		{
			double const rho = sphericalRadius( radius, z );
			if ( rho >= rho0Value )
			{
				return 0.0;
			}

			double total = 0.0;
			for ( RadialMode const &mode : modeValues )
			{
				total += mode.amplitude*radialSourceAt( mode, rho )
				         *angular.basis( mode.degree, radius, z );
			}
			return total;
		}

		/// And its derivative, which is zero: the problem is LINEAR in psi, so
		/// a Newton solve on it is affine and must finish in one step, exactly
		/// as Soloviev.hpp's does.
		double dFdPsi( double /*R*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// Delta*( psi ), by central differences of psi(), in exactly the
		/// arrangement Soloviev.hpp, CurrentLoop.hpp and VacuumHarmonic.hpp
		/// use.
		///
		/// It must come out at -f(), and that is what establishes the
		/// transcription rather than trusting the algebra -- see the file
		/// comment for the numbers, for the difference floor, and for the p = 0
		/// control that isolates the interface. Keep R well away from zero: the
		/// stencil divides by R.
		double deltaStarFD( double radius, double z, double h = 1.0e-4 ) const
		{
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( radius + h ) - innerR( radius - h ) )/( 2.0*h );
			double const dZZ = ( psi( radius, z + h ) - 2.0*psi( radius, z ) + psi( radius, z - h ) )
			                   /( h*h );

			return radius*dRInner + dZZ;
		}

		/// The coefficient vector meq::ExteriorDtN::coefficients() must return
		/// for a Gamma of this radius, one entry per mode in DEGREE ORDER from
		/// ExteriorDtN::firstMode().
		///
		/// THE SCALING IS THE POINT. ExteriorDtN expands in
		/// ( rho/rho_Gamma )^{1-n}, so its coefficient is psi on Gamma:
		///
		///     a_n = alpha_n rho_Gamma^{1-n},
		///
		/// with alpha_n this fixture's amplitude. Returning the amplitude
		/// itself would agree at rho_Gamma = 1 and nowhere else, which is why
		/// the acceptance uses two radii.
		///
		/// Every degree this fixture does not carry is an EXACT zero, and that
		/// is the discriminating part of the check rather than padding.
		///
		/// @param rhoGammaIn  the radius of Gamma. Must be at least rho_0: the
		///                    expansion is of the source-free exterior, and
		///                    inside rho_0 there is a source, so a smaller
		///                    radius is a different question and is refused.
		/// @param modeCountIn how many modes, i.e. degrees firstMode() ..
		///                    firstMode() + modeCountIn - 1. Must cover every
		///                    degree this fixture carries, or the vector cannot
		///                    represent the field and is refused rather than
		///                    silently truncated.
		std::vector<double> exteriorCoefficients( double rhoGammaIn,
		                                          int modeCountIn ) const
		{
			if ( !( rhoGammaIn >= rho0Value ) )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched::exteriorCoefficients: "
					"Gamma must enclose the source, so rhoGamma >= rho0 is "
					"required; inside rho0 the field is not a sum of exterior "
					"modes at all" );
			}
			if ( modeCountIn < 1 )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched::exteriorCoefficients: at "
					"least one mode is needed" );
			}

			int const highest = ExteriorDtN::firstMode() + modeCountIn - 1;
			for ( RadialMode const &mode : modeValues )
			{
				if ( mode.degree > highest )
				{
					throw std::invalid_argument(
						"meq::analytic::ExteriorMatched::exteriorCoefficients: "
						"degree " + std::to_string( mode.degree ) + " is "
						"carried by this field but is outside the requested "
						"range " + std::to_string( ExteriorDtN::firstMode() )
						+ " .. " + std::to_string( highest ) + ", so the vector "
						"could not represent it" );
				}
			}

			std::vector<double> out( static_cast<std::size_t>( modeCountIn ),
			                         0.0 );
			for ( RadialMode const &mode : modeValues )
			{
				std::size_t const i = static_cast<std::size_t>(
					mode.degree - ExteriorDtN::firstMode() );
				out[ i ] = mode.amplitude
				           *std::pow( rhoGammaIn,
				                      1.0 - static_cast<double>( mode.degree ) );
			}
			return out;
		}

		/// The same vector, sized and scaled from a live ExteriorDtN.
		///
		/// The safer call, and the one to prefer: it takes rho_Gamma and the
		/// mode count from the object the answer will be compared against, so
		/// the two cannot drift apart. It also REFUSES A CENTRE MISMATCH, which
		/// nothing else would catch -- two objects about different centres
		/// describe different fields, and the coefficients would be
		/// wrong rather than absent.
		std::vector<double> exteriorCoefficients( ExteriorDtN const &dtn ) const
		{
			if ( dtn.zCentre() != zCentreValue )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched::exteriorCoefficients: the "
					"ExteriorDtN is centred at a different z, so its modes are "
					"not this field's modes" );
			}
			return exteriorCoefficients( dtn.rhoGamma(), dtn.modeCount() );
		}

		/// The radial function F of one mode, at UNIT amplitude: rho^{1-n}
		/// outside rho_0 and the matched polynomial inside.
		///
		/// Exposed because it is the construction itself, and because the C^1
		/// match is checked on it directly rather than through psi().
		double radial( int degree, double rho ) const
		{
			return radialAt( modeValues[ indexOf( degree ) ], rho );
		}

		/// dF/d rho of one mode, at unit amplitude.
		double radialDerivative( int degree, double rho ) const
		{
			return radialDerivativeAt( modeValues[ indexOf( degree ) ], rho );
		}

		/// The radial factor of the SOURCE of one mode, at unit amplitude:
		/// -kappa rho^n ( rho_0^2 - rho^2 )^p inside rho_0 and a bare zero
		/// outside, so that f() is the sum of these times C_n.
		double radialSource( int degree, double rho ) const
		{
			if ( rho >= rho0Value )
			{
				return 0.0;
			}
			return radialSourceAt( modeValues[ indexOf( degree ) ], rho );
		}

		/// The amplitude of one degree in the exterior expansion, or zero if
		/// this field does not carry it. Zero rather than a throw, because
		/// "what is the coefficient of degree 4" has an answer here and it is
		/// none.
		double amplitude( int degree ) const
		{
			for ( RadialMode const &mode : modeValues )
			{
				if ( mode.degree == degree )
				{
					return mode.amplitude;
				}
			}
			return 0.0;
		}

		/// rho = |( R, z ) - centre|, the coordinate everything here is written
		/// in. One place, because the sign of z - zCentre is exactly what a
		/// second copy would get wrong.
		double sphericalRadius( double radius, double z ) const
		{
			return std::hypot( radius, z - zCentreValue );
		}

		double zCentre() const { return zCentreValue; }
		/// The support radius of the source: f is identically zero outside it.
		double rho0() const { return rho0Value; }
		/// p, the order to which the source vanishes at rho_0. psi is C^{p+1}.
		int contactOrder() const { return contactOrderValue; }
		/// How many modes this field carries, which is NOT the number of
		/// degrees spanned -- multiMode() carries three and spans four.
		int modeCount() const { return static_cast<int>( modeValues.size() ); }
		/// The degree of the i-th mode, in the order given to the constructor.
		int degree( int i ) const
		{
			return modeValues[ static_cast<std::size_t>( i ) ].degree;
		}
		/// The highest degree carried, which is what an ExteriorDtN compared
		/// against this field has to reach.
		int highestDegree() const
		{
			int highest = ExteriorDtN::firstMode();
			for ( RadialMode const &mode : modeValues )
			{
				highest = ( mode.degree > highest ) ? mode.degree : highest;
			}
			return highest;
		}

	private:
		/// One mode with its matched radial coefficients precomputed. A, kappa
		/// and the c_i are all at UNIT amplitude, so that radial() answers the
		/// construction and psi() does the scaling.
		struct RadialMode
		{
			int degree;
			double amplitude;
			double innerA;               ///< A, on rho^n
			double kappa;                ///< the source scale
			std::vector<double> inner;   ///< c_0 .. c_p, on rho^{n+2i+2}
		};

		/// Validate everything the constructor was given and return the
		/// highest degree, so that the angular object can be sized in the
		/// initialiser list before any of it has been stored.
		///
		/// STATIC AND CALLED FROM THE INITIALISER LIST on purpose: ExteriorDtN
		/// is constructed there and would otherwise report a bad radius or a
		/// bad mode count in its own words, which are about a boundary this
		/// caller has not mentioned.
		static int checkedHighestDegree( double rho0In,
		                                 std::vector<Mode> const &modesIn,
		                                 int contactOrderIn )
		{
			if ( !( rho0In > 0.0 ) )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched: the source radius rho0 "
					"must be positive" );
			}
			if ( contactOrderIn < 0 )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched: the contact order p must "
					"be at least zero; p = 0 is the two-term match, whose "
					"source is bounded but discontinuous at rho0" );
			}
			if ( modesIn.empty() )
			{
				throw std::invalid_argument(
					"meq::analytic::ExteriorMatched: at least one mode is "
					"needed" );
			}

			int highest = ExteriorDtN::firstMode();
			for ( std::size_t i = 0; i < modesIn.size(); ++i )
			{
				int const n = modesIn[ i ].degree;
				if ( n < ExteriorDtN::firstMode() )
				{
					throw std::invalid_argument(
						"meq::analytic::ExteriorMatched: mode degree "
						+ std::to_string( n ) + " is inadmissible; degrees "
						"start at 2, because C_0 and C_1 do not vanish on the "
						"axis and the mass divides by n( n - 1 )" );
				}
				for ( std::size_t j = 0; j < i; ++j )
				{
					if ( modesIn[ j ].degree == n )
					{
						throw std::invalid_argument(
							"meq::analytic::ExteriorMatched: degree "
							+ std::to_string( n ) + " is given twice; combine "
							"the amplitudes rather than relying on the order "
							"of the sum" );
					}
				}
				highest = ( n > highest ) ? n : highest;
			}
			return highest;
		}

		/// Solve the matching problem for every mode, once, in the constructor.
		///
		/// The closed forms are the file comment's, and they were checked
		/// against a from-scratch symbolic solve of the two matching conditions
		/// at n = 2..8 and p = 0..5 rather than transcribed. D is the integral
		/// of x^{2n}( 1 - x^2 )^p over [ 0, 1 ] and so is strictly positive,
		/// which is why no guard is needed on the division.
		void buildRadialCoefficients( std::vector<Mode> const &modesIn )
		{
			int const p = contactOrderValue;

			for ( Mode const &given : modesIn )
			{
				double const n = static_cast<double>( given.degree );

				double d = 0.0;
				double s0 = 0.0;
				for ( int i = 0; i <= p; ++i )
				{
					double const term = binomial( p, i )
					                    *( ( i % 2 == 0 ) ? 1.0 : -1.0 );
					double const denominator = 2.0*n + 2.0*i + 1.0;
					d += term/denominator;
					s0 += term/( 2.0*( i + 1.0 )*denominator );
				}

				RadialMode mode;
				mode.degree = given.degree;
				mode.amplitude = given.amplitude;
				mode.kappa = ( 1.0 - 2.0*n )
				             *std::pow( rho0Value, -( 2.0*n + 2.0*p + 1.0 ) )/d;
				mode.innerA = std::pow( rho0Value, 1.0 - 2.0*n )
				              - mode.kappa*std::pow( rho0Value, 2.0*p + 2.0 )*s0;

				mode.inner.resize( static_cast<std::size_t>( p ) + 1 );
				for ( int i = 0; i <= p; ++i )
				{
					double const term = binomial( p, i )
					                    *( ( i % 2 == 0 ) ? 1.0 : -1.0 );
					mode.inner[ static_cast<std::size_t>( i ) ] =
						mode.kappa*term
						*std::pow( rho0Value, 2.0*( p - i ) )
						/( 2.0*( i + 1.0 )*( 2.0*n + 2.0*i + 1.0 ) );
				}

				modeValues.push_back( mode );
			}
		}

		/// binom( p, i ), by the multiplicative recurrence. Exact in double for
		/// every p a fixture will ever be given, and it keeps this header free
		/// of a Boost include that only the .cpp of ExteriorDtN needs.
		static double binomial( int p, int i )
		{
			double result = 1.0;
			for ( int j = 0; j < i; ++j )
			{
				result = result*( p - j )/( j + 1 );
			}
			return result;
		}

		/// F( rho ) at unit amplitude. Horner in rho^2 inside, an integer power
		/// outside.
		double radialAt( RadialMode const &mode, double rho ) const
		{
			double const n = static_cast<double>( mode.degree );

			if ( rho >= rho0Value )
			{
				return std::pow( rho, 1.0 - n );
			}

			double const s = rho*rho;
			double sum = 0.0;
			for ( std::size_t i = mode.inner.size(); i > 0; --i )
			{
				sum = sum*s + mode.inner[ i - 1 ];
			}
			return std::pow( rho, n )*( mode.innerA + s*sum );
		}

		/// dF/d rho at unit amplitude.
		double radialDerivativeAt( RadialMode const &mode, double rho ) const
		{
			double const n = static_cast<double>( mode.degree );

			if ( rho >= rho0Value )
			{
				return ( 1.0 - n )*std::pow( rho, -n );
			}

			double const s = rho*rho;
			double sum = 0.0;
			for ( std::size_t i = mode.inner.size(); i > 0; --i )
			{
				double const power = n + 2.0*static_cast<double>( i - 1 ) + 2.0;
				sum = sum*s + mode.inner[ i - 1 ]*power;
			}
			return std::pow( rho, n - 1.0 )*( n*mode.innerA + s*sum );
		}

		/// The radial factor of the source at unit amplitude, INSIDE rho_0.
		/// The caller has already taken the outside branch; this one does not
		/// test, so that the zero is a literal in exactly one place.
		double radialSourceAt( RadialMode const &mode, double rho ) const
		{
			double const n = static_cast<double>( mode.degree );
			double const gap = rho0Value*rho0Value - rho*rho;

			return -mode.kappa*std::pow( rho, n )
			       *std::pow( gap, static_cast<double>( contactOrderValue ) );
		}

		/// Where a degree sits in modeValues, or a throw naming it.
		std::size_t indexOf( int degree ) const
		{
			for ( std::size_t i = 0; i < modeValues.size(); ++i )
			{
				if ( modeValues[ i ].degree == degree )
				{
					return i;
				}
			}
			throw std::invalid_argument(
				"meq::analytic::ExteriorMatched: degree "
				+ std::to_string( degree ) + " is not carried by this field. "
				"amplitude() answers zero for such a degree; the radial "
				"functions have nothing to answer" );
		}

		double zCentreValue;
		double rho0Value;
		int contactOrderValue;

		/// C_n and dC_n/dmu, from the unit under test rather than reimplemented
		/// here. Its OWN rhoGamma is set to rho0 and is never used: basis() and
		/// basisDerivative() depend on direction alone, and nothing in this
		/// fixture calls symbol(), mass() or exterior().
		ExteriorDtN angular;

		std::vector<RadialMode> modeValues;
};

}
}

#endif // MEQ_TESTS_EXTERIORMATCHED_HPP
