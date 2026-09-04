#ifndef MEQ_TESTS_ROTATINGSOLOVIEV_HPP
#define MEQ_TESTS_ROTATINGSOLOVIEV_HPP

/*
 * The ROTATING Solov'ev equilibrium: an exact solution of the generalised
 * Grad-Shafranov equation for a plasma in sonic toroidal rotation, with an
 * isothermal closure on each flux surface.
 *
 * It is the rotating counterpart of Soloviev.hpp, and it sits on the same rung
 * of the ladder tests/analytic keeps:
 *
 *   Soloviev.hpp               F constant in psi and in r beyond the r^2
 *   RotatingSoloviev.hpp       F constant in psi, EXPONENTIAL in r^2
 *   McCarthy.hpp               F linear in psi
 *   ManufacturedNonlinear.hpp  F nonlinear in psi
 *
 * WHAT IT TESTS, AND WHAT IT CANNOT. dF/dpsi is identically zero here, exactly
 * as it is for Soloviev.hpp and for the same reason: p1, F0, T0 and Omega0 are
 * all constants in the Solov'ev-with-rotation case, so the source does not
 * depend on psi at all. So this fixture measures the discretisation and the new
 * r-dependence -- which is the whole novelty of rotation, F acquiring a genuine
 * exponential in r^2 at fixed psi -- and says NOTHING WHATEVER about the Newton
 * Jacobian. FLOW-PLAN.md section 6.4 requires a manufactured nonlinear rotating
 * case for that, and this fixture is not a substitute for it.
 *
 * Source:
 *   refs/SpectralElementGSRotation.pdf  section 3.1, eqs (12), (14)-(16)
 *   refs/Refs.md   doi 10.1016/j.cpc.2020.107264 -- Li & Zhu, Comput. Phys.
 *                  Commun. 260 (2021) 107264
 *   FLOW-PLAN.md   section 6.2, which is why this fixture exists and what it
 *                  is allowed to be used for
 */

#include <array>
#include <cmath>
#include <stdexcept>

namespace meq
{
namespace analytic
{

/*
 * SIGN AND WEIGHT CONVENTION -- read this before using the class.
 *
 * Li & Zhu's eq (7) and MEQ agree on the operator:
 *
 *     Delta*( psi ) := d_rr psi - ( 1/r ) d_r psi + d_zz psi,
 *
 * which is the same thing as MEQ's r d_r( ( 1/r ) d_r psi ) + d_zz psi. Their
 * eq (12) is
 *
 *     Delta*( psi ) = -p1 r^2 exp[ M2 ( r^2/R0^2 - 1 ) ] - F0,
 *
 * and MEQ writes the equation as -Delta*( psi ) = F, so
 *
 *     F( r, z, psi ) = p1 r^2 exp[ M2 ( r^2/R0^2 - 1 ) ] + F0.            (*)
 *
 * POSITIVE, and F is the full right hand side numerator with no 1/r applied,
 * as everywhere else in this directory.
 *
 * (*) is also literally MEQ's own F = mu0 r^2 dp/dpsi + g dg/dpsi, read with
 *
 *     mu0 dp/dpsi|_r = p1 exp[ M2 ( r^2/R0^2 - 1 ) ],     g dg/dpsi = F0,
 *
 * i.e. with a pressure gradient that is a function of r as well as of psi.
 * That is the ONE structural change rotation makes to the source, and it is
 * what this fixture is for. See FLOW-PLAN.md section 2.2.
 *
 * deltaStarFD() below recomputes Delta*( psi ) by central differences and
 * RotatingSolovievConvergence.cpp asserts it against -f() over the benchmark
 * box, so the transcription of the twelve-odd terms below is checked rather
 * than trusted -- and so a closure that differs from MEQ's goes red at once
 * rather than converging beautifully to somebody else's equilibrium.
 */

/**
 * An exact rotating Solov'ev equilibrium, refs/SpectralElementGSRotation.pdf
 * eq (15):
 *
 *   psi = c1 + c2 r^2 + c3 ( r^4 - 4 r^2 z^2 ) + c4 [ r^2 ln r - z^2 ]
 *         - p1 ( R0^2/( 2 M2 ) )^2 { exp[ M2 ( r^2/R0^2 - 1 ) ]
 *                                    - ( M2/R0^2 )( r^2 - R0^2 ) - 1 }
 *         - ( F0/2 ) z^2
 *
 * The first four terms are Delta*-harmonic -- the same 1, r^2, r^4 - 4r^2z^2
 * and r^2 ln r - z^2 that Soloviev.hpp carries as psi_1, psi_2, psi_4 and (up
 * to sign) psi_3 -- so they fix the geometry and contribute nothing to F. The
 * last two are the particular solution.
 *
 * Lengths are arbitrary but consistent; r must be strictly positive, since the
 * expansion contains ln r.
 *
 * THE EXPONENT GROUP IS NAMED machSquared HERE, AND THAT IS DELIBERATE,
 * BECAUSE THE PAPER'S OWN SYMBOL IS INCONSISTENT WITH ITS OWN PROSE.
 *
 * Li & Zhu write M_0^2 in the exponent of eqs (12), (14) and (15), but their
 * prose defines
 *
 *     M_0 = m_i R_0^2 Omega_0^2 / ( 2 T_0 )                 -- as printed
 *
 * with no square root. Those two cannot both be right. Their eq (8) -- the
 * isothermal closure the whole construction rests on -- has exponent
 * coefficient m_i Omega^2 r^2 / ( 2 T ), so the group that multiplies
 * ( r^2/R_0^2 - 1 ) is
 *
 *     m_i Omega_0^2 R_0^2 / ( 2 T_0 ),
 *
 * which is what the equations call M_0^2 and what the prose calls M_0. It is
 * also dimensionally an ENERGY RATIO -- rotation energy over thermal energy --
 * that is to say a Mach number SQUARED, so the equations are right and the
 * prose is short a square root. Verified independently of the paper, both from
 * eq (8) and by requiring eq (16) to be the M -> 0 limit of eq (15).
 *
 * So this class never uses the symbol M_0 at all. machSquared is defined
 * unambiguously as
 *
 *     machSquared := m_i Omega_0^2 R_0^2 / ( 2 T_0 ),
 *
 * and anyone comparing a number against the paper must decide for themselves
 * which of its two definitions the number came from. This is the same class of
 * trap as the two coefficient errors Soloviev.hpp records, and it is written
 * out here for the same reason: it changes no convergence rate, it changes only
 * which equilibrium is being solved, and nothing in a rate table can see it.
 *
 * NUMERICAL STABILITY AS machSquared -> 0. The prefactor is 1/M2^2 and the
 * brace vanishes like M2^2, so eq (15) as printed is 0/0 and loses every digit
 * it has for small M2. It is not evaluated as printed. Writing
 *
 *     v := r^2/R0^2 - 1,     u := M2 v,
 *
 * the whole rotating term is exactly
 *
 *     -p1 ( R0^2/( 2 M2 ) )^2 ( e^u - u - 1 ) = -( p1 R0^4/4 ) v^2 G2( u ),
 *     G2( u ) := ( e^u - u - 1 )/u^2,
 *
 * in which M2 has cancelled algebraically and nothing is divided by it. G2 is
 * an entire function with G2( 0 ) = 1/2, so the M2 = 0 case needs no special
 * handling at all -- and eq (16), the static Solov'ev particular solution
 * -p1( r^2 - R0^2 )^2/8, falls out of G2( 0 ) = 1/2 rather than being a
 * separate branch. The radial derivative goes the same way, through
 * G1( u ) := ( e^u - 1 )/u.
 *
 * What DOES need care is small u, which happens for any M2 whenever r is near
 * R0: e^u - u - 1 loses all significance there. G1 and G2 therefore switch to
 * their Taylor series below |u| = seriesThreshold(). See the note on that
 * function for why both switch at the same point even though only G2 needs to.
 */
class RotatingSolovievEquilibrium
{
	public:
		/// @param majorRadiusIn  R0, the reference radius the rotation profile
		///                       is normalised to. Strictly positive.
		/// @param machSquaredIn  m_i Omega_0^2 R_0^2/( 2 T_0 ), the exponent
		///                       group. See the note above on the paper's
		///                       symbol. Zero is legal and gives a static
		///                       Solov'ev equilibrium.
		/// @param p1In           mu0 dp/dpsi at r = R0, Li & Zhu's p1.
		/// @param f0In           g dg/dpsi, Li & Zhu's F0. Constant.
		/// @param cIn            c1 ... c4, the four Delta*-harmonic
		///                       coefficients that fix the geometry.
		RotatingSolovievEquilibrium( double majorRadiusIn, double machSquaredIn,
		                             double p1In, double f0In,
		                             std::array<double, 4> const &cIn )
			: majorRadius( majorRadiusIn ), machSquared( machSquaredIn ),
			  p1( p1In ), f0( f0In ), c( cIn )
		{
			if ( !( majorRadius > 0.0 ) )
			{
				throw std::invalid_argument( "the major radius must be positive" );
			}
		}

		/*
		 * THE THREE CONFIGURATIONS BELOW, AND HOW THEIR COEFFICIENTS WERE FIXED.
		 *
		 * Li & Zhu print no c1 ... c4: their section 3.1 gives the closed form
		 * and leaves the geometry to the reader, exactly as HDG-GradShafranov.pdf
		 * does for its Examples 1-3. So each set below was SOLVED, from four
		 * conditions chosen to put a closed psi = 0 contour well inside MEQ's
		 * standard benchmark box [0.6,1.4] x [-0.6,0.6]:
		 *
		 *     psi( 1.3, 0    ) = 0        the outer equatorial point
		 *     psi( 0.7, 0    ) = 0        the inner equatorial point
		 *     psi( 1.0, 0.45 ) = 0        the high point
		 *     d_r psi( 1.0, 0.45 ) = 0    the high point is a high point
		 *
		 * i.e. a geometric centre at R0 = 1, minor radius 0.3, elongation 1.5
		 * and zero triangularity. Every term of eq (15) is even in z, so the
		 * configuration is up-down symmetric for free and four conditions are
		 * enough; there is no X-point and no need for C&F's twelve.
		 *
		 * WHY c3 AND c4 ARE NOT ZERO, WHICH IS THE ONLY REALLY LOAD-BEARING
		 * CHOICE HERE. c1 and c2 alone would have given a perfectly good closed
		 * contour. But r^4 - 4 r^2 z^2 and r^2 ln r - z^2 are Delta*-harmonic,
		 * so a TYPO IN EITHER OF THEM IS INVISIBLE TO deltaStarFD() WHEN ITS
		 * COEFFICIENT IS ZERO -- the check would pass on a fixture whose psi is
		 * wrong. Non-zero c3 and c4 are what put those two terms, and their
		 * hand-derived gradients, under the source check and the gradient check.
		 *
		 * MEASURED, at the coefficients as they are written below and evaluated
		 * in double precision, which is what the printed digits actually buy:
		 *
		 *                             stationary     rotating    fastRotating
		 *   psi( 1.3, 0 )               -6.9e-18      0.0e+00        0.0e+00
		 *   psi( 0.7, 0 )               -2.8e-17     -3.1e-17        0.0e+00
		 *   psi( 1.0, 0.45 )             0.0e+00      0.0e+00        6.9e-17
		 *   d_r psi( 1.0, 0.45 )         0.0e+00      0.0e+00        0.0e+00
		 *   max psi on the box           0.0620       0.0663         0.1139
		 *   min psi on the box          -0.2093      -0.2643        -0.7922
		 *   max psi on its boundary     -0.0346      -0.0327        -0.0368
		 *
		 * The last row is the property that matters and is the one worth
		 * re-measuring after any edit: psi is strictly negative on the whole
		 * boundary of the box and positive at the axis, so the psi = 0 level set
		 * is a closed curve strictly inside it. psi itself is O( 0.1 ) there,
		 * and analytic, r staying between 0.6 and 1.4.
		 *
		 * AS WITH Soloviev.hpp, NO CHOICE OF c CHANGES F, Delta*( psi ) OR ANY
		 * CONVERGENCE RATE. The absolute-error ceilings in
		 * RotatingSolovievConvergence.cpp are the only thing in the suite that
		 * can see a wrong one, which is why they are recorded there beside the
		 * measurement they were set from.
		 */

		/// machSquared = 0: no rotation, so this is a plain static Solov'ev
		/// equilibrium and eq (15) collapses to eq (16),
		/// psi_h - p1( r^2 - R0^2 )^2/8 - ( F0/2 ) z^2.
		///
		/// It is the control for the whole fixture. Anything measured on
		/// rotating() that is not also visible here is a property of the
		/// rotation rather than of the geometry or of the solver.
		static RotatingSolovievEquilibrium stationary()
		{
			return RotatingSolovievEquilibrium( 1.0, 0.0, 1.0, 1.0, {
				-0.086884123488195214,
				 0.12078984601508387,
				 0.02741928690676745,
				-0.30683759486827428
			} );
		}

		/// machSquared = 1: the representative rotating case, and the one the
		/// convergence study is run on.
		///
		/// The pressure gradient varies by exp( 0.96 ) = 2.6 across the box,
		/// which is a real exponential and not a perturbation, while leaving
		/// deltaStarFD() comfortably inside its truncation floor -- see
		/// fastRotating() for where that stops being true.
		static RotatingSolovievEquilibrium rotating()
		{
			return RotatingSolovievEquilibrium( 1.0, 1.0, 1.0, 1.0, {
				-0.11096425540112245,
				 0.12890481706657952,
				 0.046287651495763429,
				-0.36797424469307599
			} );
		}

		/// machSquared = 4, i.e. Mach 2, which is FLOW-PLAN.md section 8's
		/// "a factor of 7 at M = 2": exp( machSquared/2 ) = 7.4 across a flux
		/// surface, and exp( 3.84 ) = 46 across the whole box.
		///
		/// NOT USED IN THE deltaStarFD SCAN, and the reason is the check rather
		/// than the fixture. A central difference of a function whose fourth
		/// derivative carries ( 2 M2 r/R0^2 )^4 e^u has a truncation floor of
		/// about h^2 f''''/12, and h = 1e-4 is already the value that balances
		/// truncation against round-off, so there is no h that rescues it.
		/// MEASURED, max | Delta*psi + F | over the box at h = 1e-4:
		///
		///     stationary   1.6e-08      rotating   2.9e-07
		///     fastRotating 4.8e-05
		///
		/// so the scan's 1e-5 -- which is Soloviev.hpp's own figure, and which
		/// the first two clear by three and two orders -- would have to be
		/// loosened by an order for this one alone. A looser tolerance on the
		/// single guard standing between a mistyped term and a solver that
		/// converges beautifully to the wrong equilibrium is not worth an extra
		/// Mach number, so the scan runs at machSquared = 0 and 1 and this
		/// configuration is kept for the solve.
		static RotatingSolovievEquilibrium fastRotating()
		{
			return RotatingSolovievEquilibrium( 1.0, 4.0, 1.0, 1.0, {
				-0.20159393373560924,
				 0.14242973485054433,
				 0.15288363162526158,
				-0.64872251296921124
			} );
		}

		/// The poloidal flux function, eq (15), evaluated in the form that is
		/// stable as machSquared -> 0. See the class comment.
		double psi( double r, double z ) const
		{
			double const r2 = r*r;
			double const z2 = z*z;
			double const logR = std::log( r );

			double value = c[ 0 ]
			             + c[ 1 ]*r2
			             + c[ 2 ]*( r2*r2 - 4.0*r2*z2 )
			             + c[ 3 ]*( r2*logR - z2 );

			value += rotatingTerm( r );
			value -= 0.5*f0*z2;
			return value;
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ), differentiated by hand rather
		/// than differenced. Not the HDG flux: that is this divided by r, see
		/// flux().
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			double const r2 = r*r;
			double const z2 = z*z;
			double const logR = std::log( r );

			dPsiDr = 2.0*c[ 1 ]*r
			       + c[ 2 ]*( 4.0*r2*r - 8.0*r*z2 )
			       + c[ 3 ]*( 2.0*r*logR + r )
			       + rotatingTermPrime( r );

			dPsiDz = -8.0*c[ 2 ]*r2*z
			       - 2.0*c[ 3 ]*z
			       - f0*z;
		}

		/// The HDG flux q = grad_bar( psi )/r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The Grad-Shafranov source,
		/// F = p1 r^2 exp[ M2 ( r^2/R0^2 - 1 ) ] + F0, from eq (12) and MEQ's
		/// -Delta*( psi ) = F. Returns F, not F/r.
		///
		/// Independent of psi, which is what makes this equilibrium linear --
		/// and EXPONENTIAL in r^2, which is what makes it a new test rather
		/// than a restatement of Soloviev.hpp.
		double f( double r, double /*z*/, double /*psiValue*/ ) const
		{
			double const v = r*r/( majorRadius*majorRadius ) - 1.0;
			return p1*r*r*std::exp( machSquared*v ) + f0;
		}

		/// dF/dpsi. Identically zero: p1, F0, T0 and Omega0 are all constants
		/// in the Solov'ev-with-rotation case, so the source does not depend on
		/// psi at all and a Newton solve on this problem must converge in a
		/// single step.
		///
		/// It is also why this fixture cannot test a Jacobian, and why
		/// FLOW-PLAN.md section 6.4 asks for a manufactured nonlinear rotating
		/// case as well.
		double dFdPsi( double, double, double ) const
		{
			return 0.0;
		}

		/// Delta*( psi ), by central differences of psi(). Used to verify the
		/// sign convention documented above, and to catch a mistyped term in
		/// the expansion, without depending on the hand-derived gradients.
		///
		/// The implementation McCarthy.hpp carries, and Soloviev.hpp's but for
		/// its two missing consts, deliberately: this is the one check every
		/// fixture in this directory shares, and it is worth being able to diff
		/// them and see nothing.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			// r d_r( ( 1/r ) d_r psi ) as a second difference of the inner
			// quantity, plus d_zz psi.
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) ) / ( 2.0 * h ) / rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) ) / ( 2.0 * h );
			double const dZZ = ( psi( r, z + h ) - 2.0 * psi( r, z ) + psi( r, z - h ) )
			                 / ( h * h );

			return r * dRInner + dZZ;
		}

		/// R0, the radius the rotation profile is referenced to.
		double getMajorRadius() const
		{
			return majorRadius;
		}

		/// m_i Omega_0^2 R_0^2/( 2 T_0 ). See the class comment on why it is
		/// not called M_0 or M_0^2.
		double getMachSquared() const
		{
			return machSquared;
		}

		/// Li & Zhu's p1 = mu0 dp/dpsi at r = R0.
		double getP1() const
		{
			return p1;
		}

		/// Li & Zhu's F0 = g dg/dpsi, constant.
		double getF0() const
		{
			return f0;
		}

		/// c1 ... c4, the Delta*-harmonic coefficients.
		std::array<double, 4> const &getCoefficients() const
		{
			return c;
		}

		/// |u| below which G1 and G2 are summed rather than evaluated, with u
		/// the exponent machSquared( r^2/R0^2 - 1 ). Public because it is the
		/// one number a test of the crossover has to know: a sweep that never
		/// straddles it measures one branch twice.
		///
		/// 0.5 costs G2 about one significant digit to cancellation on the
		/// closed-form side ( e^0.5 = 1.649 against e^0.5 - 1.5 = 0.1487 ), and
		/// the series needs a dozen terms to reach round-off on the other, so
		/// there is a wide plateau either side and the exact value is not
		/// critical.
		///
		/// BOTH FUNCTIONS SWITCH HERE EVEN THOUGH ONLY G2 NEEDS TO. G1 could
		/// use std::expm1 for every non-zero u and lose nothing, since expm1
		/// has no cancellation to suffer. Sharing the threshold means the value
		/// and the radial derivative change branch at the same place, so ONE
		/// continuity test covers both -- and there is one place to move rather
		/// than two if it ever needs moving.
		static double seriesThreshold()
		{
			return 0.5;
		}

	private:
		double majorRadius;
		double machSquared;
		double p1;
		double f0;
		std::array<double, 4> c;

		/// ( e^u - 1 )/u, with its removable singularity filled: G1( 0 ) = 1.
		static double g1( double u )
		{
			if ( std::abs( u ) < seriesThreshold() )
			{
				// sum_{k >= 0} u^k/( k + 1 )!, by the ratio u/( k + 1 ). At
				// |u| <= 0.5 the twenty-fifth term is below 1e-40, so the loop
				// bound is a fixed count rather than a tolerance -- a
				// deterministic number of flops, and no branch that could
				// behave differently on two builds.
				double term = 1.0;
				double sum = 1.0;
				for ( int k = 1; k < 25; ++k )
				{
					term *= u/( k + 1 );
					sum += term;
				}
				return sum;
			}
			return std::expm1( u )/u;
		}

		/// ( e^u - u - 1 )/u^2, with its removable singularity filled:
		/// G2( 0 ) = 1/2.
		static double g2( double u )
		{
			if ( std::abs( u ) < seriesThreshold() )
			{
				// sum_{k >= 0} u^k/( k + 2 )!, by the ratio u/( k + 2 ).
				double term = 0.5;
				double sum = 0.5;
				for ( int k = 1; k < 25; ++k )
				{
					term *= u/( k + 2 );
					sum += term;
				}
				return sum;
			}
			return ( std::exp( u ) - u - 1.0 )/( u*u );
		}

		/// The particular solution's rotating part,
		///
		///     -p1 ( R0^2/( 2 M2 ) )^2 { e^u - ( M2/R0^2 )( r^2 - R0^2 ) - 1 }
		///       = -( p1 R0^4/4 ) v^2 G2( u ),
		///
		/// with v = r^2/R0^2 - 1 and u = M2 v. M2 has cancelled: this is finite
		/// and accurate at machSquared = 0, where G2( 0 ) = 1/2 recovers
		/// eq (16)'s -p1( r^2 - R0^2 )^2/8.
		double rotatingTerm( double r ) const
		{
			double const r0 = majorRadius;
			double const v = r*r/( r0*r0 ) - 1.0;
			double const u = machSquared*v;
			double const r0Fourth = r0*r0*r0*r0;
			return -0.25*p1*r0Fourth*v*v*g2( u );
		}

		/// d/dr of rotatingTerm(). Differentiating the printed form gives
		/// -( p1 R0^2 r/( 2 M2 ) )( e^u - 1 ), and the same cancellation of M2
		/// applies: e^u - 1 = u G1( u ) = M2 v G1( u ), leaving
		/// -( p1 R0^2 r v/2 ) G1( u ).
		double rotatingTermPrime( double r ) const
		{
			double const r0 = majorRadius;
			double const v = r*r/( r0*r0 ) - 1.0;
			double const u = machSquared*v;
			return -0.5*p1*r0*r0*r*v*g1( u );
		}
};

}
}

#endif // MEQ_TESTS_ROTATINGSOLOVIEV_HPP
