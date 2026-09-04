#ifndef MEQ_TESTS_HIGHBETAPOLOIDAL_HPP
#define MEQ_TESTS_HIGHBETAPOLOIDAL_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

/*
 * HIGH POLOIDAL BETA, WITH THE PROFILES SPECIFIED AS POLYNOMIALS IN NORMALISED
 * FLUX -- which is how equilibrium codes actually pose the problem, and which
 * every other source in tests/analytic/ avoids.
 *
 * refs/GourdainContour.pdf section V, eq (39):
 *
 *     p( Psi )   = sum_{i>0} a_i Psi^i
 *     F^2( Psi ) = F_axis^2 - ( F_axis^2 - F_bnd^2 ) sum_{i>0} b_i Psi^i
 *     Psi        = ( psi - psi_bnd ) / ( psi_axis - psi_bnd )
 *
 * and the paper's reason for that form is worth keeping in view: "this concise
 * decomposition can match a lot of real plasma configurations". It is the shape
 * meq::Profile has, and the shape a free-boundary solver has to live with.
 *
 * WHY IT IS HERE AT ALL: it is a candidate STIFF case, and unlike GS-2's
 * sections 4.2 to 4.5 it is stiff for a physical reason rather than a chosen
 * one. refs/CowleyHighBeta.pdf and refs/HsuCowleyHighBeta.pdf construct these
 * equilibria asymptotically in the limit beta_p >> 1, where the core is pushed
 * outward into a thin boundary layer against the separatrix -- a genuinely
 * multi-scale problem -- and refs/GourdainContour.pdf demonstrates that they can
 * be computed numerically, reporting that its own contour-dynamics scheme
 * "converges for extreme high beta configurations" where flux-conserving methods
 * do not.
 *
 * STIFFNESS IS A HYPOTHESIS UNTIL IT IS MEASURED, and this project has been
 * wrong about it before: sections 4.2, 4.3 and 4.5 were called stiff for months
 * and were under-resolved. HighBetaConvergence.cpp measures
 * max| dF/dpsi | / lambda_1 over the range the CONVERGED SOLUTION occupies --
 * which is the only range worth sampling, and which needs the normalisation to
 * be solved for rather than assumed -- so the question is answered rather than
 * assumed.
 *
 * psi_axis IS SETTABLE, AND THE SOLVER IS WHAT SETS IT. psi_bnd is zero, because
 * MEQ solves the fixed-boundary problem with psi = 0 on Gamma, so
 * Psi = psi / psi_axis throughout and psi_axis is the whole of the
 * normalisation. It arrives as a constructor argument and setPsiAxis() moves it,
 * which is what meq::NormalisedSource requires: psi_axis is an UNKNOWN of the
 * non-linear system, carried in the bordered Newton of
 * GradShafranovSolver::setSource( NormalisedSource &, double ), and the solver
 * writes it here before every residual evaluation.
 *
 * WHY IT CANNOT BE FIXED, which is what measuring it turned up and what the
 * solver work exists for. Hand the solver a psi_axis the solution does not reach
 * and the profile is never sampled: with psi_axis = 1 and a peaked pressure the
 * computed solution peaks at Psi = 0.0013, so Psi^(nu-1) is around 1e-9 and the
 * pressure term contributes NOTHING. Measured, the answers at amplitude 1 and
 * amplitude 512 were identical to every digit -- psi in
 * [ 2.570e-07, 1.259e-03 ] for both -- because the solve was driven entirely by
 * the constant g g' term and the pressure profile might as well not have been
 * there. Every case converged in one or two Newton steps, at reaction ratios up
 * to 2523, for the same reason: there was no non-linearity switched on.
 *
 * AND THE SELF-CONSISTENT SOLUTION IS NOT THE ONE NEWTON FINDS FROM ZERO. At a
 * FIXED psi_axis the equation has two positive solutions: a small one, which is
 * what an unguided Newton lands on and where the profile is inert, and a large
 * one with Psi of order 1, which is the equilibrium. Measured at nu = 4,
 * amplitude 1, psi_axis = 0.42: 3.00e-3 from zero, 6.23e-1 from a bump. Only the
 * large one can satisfy max psi = psi_axis, so the constraint removes the small
 * branch from the solution set -- but not from the iteration's reach, which is
 * why a run of this fixture needs an initial guess of about the right height.
 *
 * Closing the loop OUTSIDE the solver does not rescue it either: the map
 * psi_axis -> max psi has a pole beside its own fixed point (at nu = 2,
 * amplitude 1 the fixed point is 0.3059 and the pole is at 0.2996), so an outer
 * iteration falls off the branch. That is the measurement that put psi_axis
 * inside the residual.
 */

namespace meq
{
namespace analytic
{

class HighBetaPoloidal
{
	public:
		/**
		 * @param pressureIn    the a_i of eq (39), starting at a_1. p( 0 ) = 0 at
		 *                      the boundary by construction, which is why there is
		 *                      no a_0.
		 * @param toroidalIn    the b_i, starting at b_1, with the same convention.
		 * @param fAxisSquared  F^2 on the magnetic axis.
		 * @param fBndSquared   F^2 on the boundary. The vacuum toroidal field is
		 *                      F = R B_T, so this is the one that is really known.
		 * @param psiAxisIn     the flux on the axis, as a starting value:
		 *                      setPsiAxis() moves it and the solver owns it. See
		 *                      the file comment.
		 * @param mu0In         1 for a problem in normalised units.
		 */
		HighBetaPoloidal( std::vector<double> pressureIn,
		                  std::vector<double> toroidalIn,
		                  double fAxisSquared, double fBndSquared,
		                  double psiAxisIn, double mu0In )
			: a( std::move( pressureIn ) ), b( std::move( toroidalIn ) ),
			  fAxisSq( fAxisSquared ), fBndSq( fBndSquared ),
			  psiAxisValue( psiAxisIn ), mu0Value( mu0In )
		{
			if ( !( psiAxisValue > 0.0 ) )
				throw std::invalid_argument( "HighBetaPoloidal: psi_axis must be positive" );
		}

		/**
		 * A moderate-beta case, and the reference point everything else is scaled
		 * against.
		 *
		 * p = a1 Psi + a2 Psi^2 with a1 = a2 = 0.5, so p( axis ) = 1 and dp/dPsi
		 * is 0.5 at the boundary and 1.5 on the axis -- a peaked but unremarkable
		 * profile. F^2 falls by 4% from axis to boundary, which for F_bnd = 1 is
		 * a diamagnetic shift of the size an ordinary tokamak has.
		 */
		static HighBetaPoloidal moderate()
		{
			return HighBetaPoloidal( { 0.5, 0.5 }, { 1.0 }, 1.04, 1.0, 1.0, 1.0 );
		}

		/**
		 * A PEAKED profile: p proportional to Psi^nu, so the pressure gradient is
		 * concentrated near the axis and F is genuinely non-linear in psi.
		 *
		 * moderate() above is NOT, and measuring it is what showed the difference
		 * matters. With p = a1 Psi + a2 Psi^2 the pressure derivative is linear in
		 * psi and, with a single b_1, g g' is constant -- so F is LINEAR in psi
		 * and Newton finishes in one step whatever the amplitude. Measured: at a
		 * reaction ratio of 22.5, comparable to the GS-2 current hole's 26, it
		 * still took one iteration. A large dF/dpsi is not by itself a hard
		 * problem; a large VARIATION in dF/dpsi is.
		 *
		 * nu = 1 or 2 reproduces the mild case. nu >= 3 gives a dF/dpsi that
		 * varies over the domain, which is what makes the Jacobian move between
		 * Newton steps.
		 *
		 * @param nu         the peaking exponent, at least 1.
		 * @param amplitude  p on the axis, since Psi = 1 there.
		 */
		static HighBetaPoloidal peaked( int nu, double amplitude, double psiAxisIn = 1.0 )
		{
			if ( nu < 1 )
				throw std::invalid_argument( "HighBetaPoloidal::peaked: nu must be at least 1" );

			std::vector<double> pressure( static_cast<std::size_t>( nu ), 0.0 );
			pressure.back() = amplitude;
			return HighBetaPoloidal( pressure, { 1.0 }, 1.04, 1.0, psiAxisIn, 1.0 );
		}

		/// The same shape with the pressure scaled by @a factor and nothing else
		/// touched, which is how the beta sweep is built. Scaling the pressure at
		/// fixed F is what raises beta_p.
		HighBetaPoloidal scaledPressure( double factor ) const
		{
			std::vector<double> scaled( a );
			for ( double &value : scaled )
				value *= factor;
			return HighBetaPoloidal( scaled, b, fAxisSq, fBndSq, psiAxisValue, mu0Value );
		}

		/// Move the normalisation. This is what meq::NormalisedSource's
		/// setNormalisation() forwards to, and the solver calls it once per
		/// residual evaluation, so it does no work beyond the assignment.
		void setPsiAxis( double psiAxisIn )
		{
			if ( !( psiAxisIn > 0.0 ) )
				throw std::invalid_argument( "HighBetaPoloidal: psi_axis must be positive" );
			psiAxisValue = psiAxisIn;
		}

		/// Normalised flux. psi_bnd = 0, so this is psi / psi_axis.
		double normalised( double psi ) const
		{
			return psi/psiAxisValue;
		}

		/// p( psi ), for reporting beta rather than for the solve.
		double p( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = s;
			for ( double coefficient : a )
			{
				sum += coefficient*power;
				power *= s;
			}
			return sum;
		}

		/// dp/dpsi. The chain rule contributes one factor of 1/psi_axis, and that
		/// factor is the whole reason a normalised profile is not the same object
		/// as an unnormalised one.
		double pPrime( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = 1.0;
			for ( std::size_t i = 0; i < a.size(); ++i )
			{
				sum += static_cast<double>( i + 1 )*a[ i ]*power;
				power *= s;
			}
			return sum/psiAxisValue;
		}

		double pDoublePrime( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = 1.0;
			for ( std::size_t i = 1; i < a.size(); ++i )
			{
				sum += static_cast<double>( i + 1 )*static_cast<double>( i )*a[ i ]*power;
				power *= s;
			}
			return sum/( psiAxisValue*psiAxisValue );
		}

		/// F^2( psi ), eq (39). Falls from F_axis^2 to F_bnd^2 as Psi goes 1 -> 0
		/// when the b_i sum to one, which is the normalisation the paper intends.
		double fSquared( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = s;
			for ( double coefficient : b )
			{
				sum += coefficient*power;
				power *= s;
			}
			return fAxisSq - ( fAxisSq - fBndSq )*( 1.0 - sum );
		}

		/// g dg/dpsi, which is half of d( F^2 )/dpsi. That identity is the reason
		/// eq (2)'s second term is usually written as F F' and MEQ writes it as
		/// g g': they are the same quantity and this fixture never has to form g
		/// itself, which would need a square root and a sign convention.
		double ggPrime( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = 1.0;
			for ( std::size_t i = 0; i < b.size(); ++i )
			{
				sum += static_cast<double>( i + 1 )*b[ i ]*power;
				power *= s;
			}
			return 0.5*( fAxisSq - fBndSq )*sum/psiAxisValue;
		}

		double ggDoublePrime( double psi ) const
		{
			double const s = normalised( psi );
			double sum = 0.0;
			double power = 1.0;
			for ( std::size_t i = 1; i < b.size(); ++i )
			{
				sum += static_cast<double>( i + 1 )*static_cast<double>( i )*b[ i ]*power;
				power *= s;
			}
			return 0.5*( fAxisSq - fBndSq )*sum/( psiAxisValue*psiAxisValue );
		}

		/// F( r, z, psi ) = mu0 r^2 dp/dpsi + g dg/dpsi, eq (2) exactly. No z
		/// dependence: a profile equilibrium has none.
		double f( double r, double, double psi ) const
		{
			return mu0Value*r*r*pPrime( psi ) + ggPrime( psi );
		}

		double dFdPsi( double r, double, double psi ) const
		{
			return mu0Value*r*r*pDoublePrime( psi ) + ggDoublePrime( psi );
		}

		double psiAxis() const { return psiAxisValue; }
		double mu0() const { return mu0Value; }

	private:
		std::vector<double> a;
		std::vector<double> b;
		double fAxisSq;
		double fBndSq;
		double psiAxisValue;
		double mu0Value;
};

}
}

#endif // MEQ_TESTS_HIGHBETAPOLOIDAL_HPP
