#ifndef MEQ_TESTS_MANUFACTUREDANISOTROPIC_HPP
#define MEQ_TESTS_MANUFACTUREDANISOTROPIC_HPP

/*
 * A manufactured solution of the ANISOTROPIC Grad-Shafranov equation:
 * Example 5's flux (ManufacturedNonlinear.hpp) under a prescribed anisotropy.
 *
 * Source:
 *   refs/AlphaGK.pdf   eq (6.18); refs/Refs.md records its two typos, neither
 *                      of which touches the operator used here
 *   ANISOTROPY-PLAN.md AN-0, which this fixture is, and AN-1 and AN-2, which
 *                      solve it
 *   tools/anisotropic_equation.py  which checks the closed forms below
 *                      symbolically and by 50-digit central differences, and
 *                      prints the reference values the convergence test pins
 *
 * THE EQUATION, in MEQ's convention:
 *
 *     -div_bar( ( 1 - Delta )/R grad_bar psi ) = F( R, z ) / R,
 *
 * whose left-hand side times -R is (6.18)'s, ( 1 - Delta ) Delta* psi -
 * grad psi . grad Delta. The diffusion coefficient 1/R becomes
 * ( 1 - Delta )/R, and with it the HDG flux: q = ( 1 - Delta )/R grad_bar psi.
 * flux() below returns THAT q, the one an anisotropic solve carries, and not
 * grad_bar( psi )/R.
 *
 * THE FLUX is Example 5's, unchanged, so the domain, the meshes and the
 * harness of NewtonConvergence.cpp carry over:
 *
 *     psi( R, z ) = sin( kr ( R + R_0 ) ) cos( kz z ),
 *     R_0 = -0.5,   kr = 1.15 pi,   kz = 1.15,
 *
 * and F_0 := -Delta*( psi ) = ( kr^2 + kz^2 ) psi
 *                             + ( kr/R ) cos( kr ( R + R_0 ) ) cos( kz z ),
 * which is Example 5's F at its own solution.
 *
 * TWO MODELS OF Delta, because AN-1 and AN-2 are different solves:
 *
 *   Position (AN-1)   Delta = a cos( alpha ( R - 1 ) ) cos( beta z ), a FIELD.
 *                     The flux mass stays linear. Default a = -0.3,
 *                     alpha = 2.5, beta = 2, so Delta takes both signs on
 *                     the box and a sign error in the grad Delta term cannot
 *                     hide in a one-signed coefficient.
 *
 *                     F = ( 1 - Delta ) F_0 + grad Delta . grad psi
 *
 *   Flux (AN-2)       Delta = D( psi ) = d psi^2, a function of the
 *                     SOLUTION. The flux mass is non-linear and couples to
 *                     the potential. Default d = -0.5, so 1 - Delta lies in
 *                     [ 1, 1.5 ] wherever |psi| <= 1.
 *
 *                     F = ( 1 - D( psi ) ) F_0 + D'( psi ) |grad psi|^2
 *
 * IN BOTH, F IS A FUNCTION OF POSITION ONLY: the manufactured right-hand side
 * is evaluated at the exact flux, and every non-linearity of the flux model is
 * in the operator, where an anisotropic solve puts it. So dF/dpsi = 0, and a
 * solve of the flux model that converges in one Newton step has lost its
 * operator's dependence on psi.
 *
 * f() RETURNS F, NOT F/R, as ManufacturedNonlinear.hpp's does and for the
 * reason it gives at length: meq::GradShafranovSolver applies the 1/R.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

class ManufacturedAnisotropic
{
	public:
		enum class DeltaModel
		{
			Position,
			Flux
		};

		/// Delta = amplitude cos( alpha ( R - 1 ) ) cos( beta z ).
		static ManufacturedAnisotropic position( double amplitude = -0.3, double alpha = 2.5,
		                                         double beta = 2.0 )
		{
			return ManufacturedAnisotropic( DeltaModel::Position, amplitude, alpha, beta );
		}

		/// Delta = amplitude psi^2.
		static ManufacturedAnisotropic flux( double amplitude = -0.5 )
		{
			return ManufacturedAnisotropic( DeltaModel::Flux, amplitude, 0.0, 0.0 );
		}

		DeltaModel model() const
		{
			return modelValue;
		}

		/// The poloidal flux function, Example 5's.
		double psi( double radius, double z ) const
		{
			return std::sin( kr*( radius + r0 ) )*std::cos( kz*z );
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ).
		void gradPsi( double radius, double z, double &dPsiDr, double &dPsiDz ) const
		{
			dPsiDr =  kr*std::cos( kr*( radius + r0 ) )*std::cos( kz*z );
			dPsiDz = -kz*std::sin( kr*( radius + r0 ) )*std::sin( kz*z );
		}

		/// Delta at the exact solution, under either model.
		double delta( double radius, double z ) const
		{
			if ( modelValue == DeltaModel::Flux )
				return deltaOfPsi( psi( radius, z ) );
			return amplitudeValue*std::cos( alphaValue*( radius - 1.0 ) )*std::cos( betaValue*z );
		}

		/// grad_bar( Delta ) at the exact solution, under either model.
		void gradDelta( double radius, double z, double &dDeltaDr, double &dDeltaDz ) const
		{
			if ( modelValue == DeltaModel::Flux )
			{
				double dPsiDr = 0.0;
				double dPsiDz = 0.0;
				gradPsi( radius, z, dPsiDr, dPsiDz );
				double const slope = dDeltaDPsi( psi( radius, z ) );
				dDeltaDr = slope*dPsiDr;
				dDeltaDz = slope*dPsiDz;
				return;
			}
			double const cr = std::cos( alphaValue*( radius - 1.0 ) );
			double const sr = std::sin( alphaValue*( radius - 1.0 ) );
			double const cz = std::cos( betaValue*z );
			double const sz = std::sin( betaValue*z );
			dDeltaDr = -amplitudeValue*alphaValue*sr*cz;
			dDeltaDz = -amplitudeValue*betaValue*cr*sz;
		}

		/// D( psi ), the flux model's Delta as a function of the solution. Only
		/// meaningful under DeltaModel::Flux.
		double deltaOfPsi( double psiValue ) const
		{
			return amplitudeValue*psiValue*psiValue;
		}

		/// D'( psi ). Only meaningful under DeltaModel::Flux.
		double dDeltaDPsi( double psiValue ) const
		{
			return 2.0*amplitudeValue*psiValue;
		}

		/// The anisotropic HDG flux q = ( 1 - Delta )/R grad_bar( psi ).
		void flux( double radius, double z, double &qR, double &qZ ) const
		{
			gradPsi( radius, z, qR, qZ );
			double const weight = ( 1.0 - delta( radius, z ) )/radius;
			qR *= weight;
			qZ *= weight;
		}

		/// F_0 = -Delta*( psi ), Example 5's source at its own solution.
		double isotropicSource( double radius, double z ) const
		{
			return ( kr*kr + kz*kz )*psi( radius, z )
			       + ( kr/radius )*std::cos( kr*( radius + r0 ) )*std::cos( kz*z );
		}

		/// The source F( R, z ), NOT F/R. Independent of psi under both models.
		double f( double radius, double z ) const
		{
			double dPsiDr = 0.0;
			double dPsiDz = 0.0;
			gradPsi( radius, z, dPsiDr, dPsiDz );
			double dDeltaDr = 0.0;
			double dDeltaDz = 0.0;
			gradDelta( radius, z, dDeltaDr, dDeltaDz );
			return ( 1.0 - delta( radius, z ) )*isotropicSource( radius, z )
			       + dDeltaDr*dPsiDr + dDeltaDz*dPsiDz;
		}

		/// R div_bar( ( 1 - Delta )/R grad_bar psi ), by central differences of
		/// psi() and delta() alone. Must equal -f(): the check that catches a
		/// closed form that is consistent with itself and wrong, which no
		/// convergence rate can see.
		double anisotropicOperatorFD( double radius, double z, double h = 1.0e-4 ) const
		{
			auto fluxR = [ & ]( double rr, double zz )
			{
				return ( 1.0 - delta( rr, zz ) )/rr*( psi( rr + h, zz ) - psi( rr - h, zz ) )/( 2.0*h );
			};
			auto fluxZ = [ & ]( double rr, double zz )
			{
				return ( 1.0 - delta( rr, zz ) )/rr*( psi( rr, zz + h ) - psi( rr, zz - h ) )/( 2.0*h );
			};

			double const divergence = ( fluxR( radius + h, z ) - fluxR( radius - h, z ) )/( 2.0*h )
			                          + ( fluxZ( radius, z + h ) - fluxZ( radius, z - h ) )/( 2.0*h );
			return radius*divergence;
		}

	private:
		ManufacturedAnisotropic( DeltaModel modelIn, double amplitudeIn, double alphaIn,
		                         double betaIn )
			: modelValue( modelIn ), amplitudeValue( amplitudeIn ), alphaValue( alphaIn ),
			  betaValue( betaIn )
		{
		}

		static constexpr double pi = 3.14159265358979323846;
		static constexpr double r0 = -0.5;
		static constexpr double kr = 1.15*pi;
		static constexpr double kz = 1.15;

		DeltaModel modelValue;
		double amplitudeValue;
		double alphaValue;
		double betaValue;
};

}
}

#endif
