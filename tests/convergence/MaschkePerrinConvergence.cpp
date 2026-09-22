#define BOOST_TEST_MODULE MaschkePerrinConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"

#include "analytic/MaschkePerrin.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * The HDG Grad-Shafranov operator, and meq::RotatingSource, measured against
 * Maschke & Perrin's exact rotating equilibrium -- refs/MaschkePerrin.pdf
 * section 4, Plasma Physics 22 (1980) 579.
 *
 * WHAT IS NEW HERE, AND IT IS NOT THE DISCRETISATION. Maschke & Perrin's (4.10)
 * is Li & Zhu's (12) renamed, and their (4.17) is a member of the same solution
 * family, so a study driven from the fixture's own f() measures exactly what
 * RotatingSolovievConvergence.cpp measures. That study is here as the CONTROL
 * and is asserted with the linear benchmark's slack, but it is not the point.
 *
 * The point is the SOURCE. Every closed-form check of meq::RotatingSource in
 * this tree runs at CONSTANT temperature and CONSTANT omega:
 * RotatingSourceConvergence.cpp builds its species on ConstantMassProfile for
 * both, because Li & Zhu's rotating Solov'ev case has T_0 and Omega_0 constant
 * and no closed form survives otherwise. Maschke & Perrin's (4.7) constrains
 * only the RATIO omega^2/( Rbar T ), so T( psi ) stays an arbitrary surface
 * function and omega( psi ) follows it. That is the one exact solution MEQ has
 * that admits psi-DEPENDENT profiles, and
 * theRotatingSourceReproducesTheClosedForm is what it is for: five profiles
 * varying in psi, whose variations must cancel to a source independent of psi.
 *
 * The alternative check on that region is a central difference of MEQ's own
 * f(), which is what RotatingSourceTests.cpp does. CLAUDE.md records why it is
 * not the same thing: f() and dFdPsi() are evaluated at the same profiles and
 * agree with each other however wrong both are. A closed form does not.
 *
 * WHAT IS STILL NOT COVERED. (4.7) makes the shared exponent coefficient C a
 * CONSTANT -- that is precisely what collapses (4.6) to (4.8) -- so C'( psi ) is
 * zero here as it is for Li & Zhu, and the C' term remains touched by nothing
 * but RotatingSourceTests' dFdPsi sweep. No exact rotating solution can close
 * that, because a varying C is what makes the equation unsolvable in closed
 * form. Do not read this file as having closed it.
 *
 * ORDER OF THE TEST CASES IS LOAD BEARING. The Delta* scan and the two
 * geometric checks come first, because everything below them is measured
 * against a closed form that would otherwise be taken on trust; a mistyped term
 * in (3.16) or (4.16), or section 3's polytrope closure mistaken for section
 * 4's isothermal one, converges at the right rate to the wrong function and no
 * rate table in this file could tell.
 *
 * The domain is meq::tests::standardBox(), the same rectangle
 * [0.6,1.4] x [-0.6,0.6] SolovievConvergence.cpp and
 * RotatingSolovievConvergence.cpp use, so the three studies differ in the source
 * and in nothing else.
 */

namespace
{

	using Equilibrium = meq::analytic::MaschkePerrinEquilibrium;

	/// The scan SolovievConvergence.cpp uses, at the same step: 0.2 in R and
	/// 0.3 in z over the benchmark box.
	template<typename Check>
	void overTheBox( Check check )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		for ( double radius = box.minRadius; radius <= box.maxRadius + 1.0e-12; radius += 0.2 )
		{
			for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.3 )
			{
				check( radius, z );
			}
		}
	}

	/// (4.16)'s g_T EXACTLY AS THE PAPER PRINTS IT, prefactor 1/( gamma^2
	/// Omega^4 ) and all. Not what the fixture evaluates, and that is the point:
	/// it is the control for theSmallMachLimitIsContinuous, which asserts that
	/// this form loses everything for small m while the fixture does not.
	double naiveParticularAsPrinted( double radius, double machSquared, double pressure,
	                                 double majorRadius )
	{
		double const w = machSquared*radius*radius/( 2.0*majorRadius*majorRadius );
		return pressure/( machSquared*machSquared )*( 1.0 + w - std::exp( w ) );
	}

	/*
	 * ------------------------------------------------------------------
	 * The section 4 plasma, as meq::RotatingSource takes it
	 * ------------------------------------------------------------------
	 *
	 * MEQ's rotating source is RoPP (96) closed by (97): two species whose
	 * densities carry exp[ m_s omega^2( R^2 - R_ref^2 )/2T_s - Z_s e phi_0/T_s ],
	 * with phi_0 fixed by quasineutrality. At two species the exponents come out
	 * equal and shared,
	 *
	 *     C( psi ) = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ),
	 *
	 * so p = Sum_s n_s T_s is p_0( psi ) exp[ C( psi )( R^2 - R_ref^2 )/2 ] with
	 * p_0 = Sum_s n_s0 T_s the pressure on the reference curve.
	 *
	 * MASCHKE & PERRIN'S SECTION 4 IS EXACTLY TWO CONDITIONS ON THAT:
	 *
	 *   (4.7)  omega^2/( Rbar T ) constant  <=>  C constant,
	 *   (4.9)  p_T linear in F              <=>  p_0 linear in psi,
	 *
	 * and then F = mu0 R^2 p_0' exp[ C( R^2 - R_ref^2 )/2 ] + g g' IS (4.10).
	 * Setting R_ref = R0 and C = m/R0^2 makes the exponent m R^2/2R0^2 - m/2, so
	 * the amplitude match is p_0' = ( P/mu0 R0^4 ) exp( m/2 ).
	 *
	 * THE CONSTRUCTION BELOW IS THE CHEAPEST WAY TO SATISFY BOTH WHILE LEAVING
	 * EVERY PROFILE VARYING, and each piece is there for a reason:
	 *
	 *     theta( psi ) = 1 + sigma psi      one shape function
	 *     T_s( psi )   = tau_s theta        so Z_1 T_2 - Z_2 T_1 carries theta
	 *     omega( psi ) = omega_0 sqrt theta so omega^2 carries theta, and the
	 *                                       two thetas cancel in C
	 *     n_s0( psi )  = c_s N( psi )/theta with N linear and Sum_s Z_s c_s = 0,
	 *                                       so each n_s0 T_s = c_s tau_s N is
	 *                                       linear and neutrality is exact
	 *
	 * omega( psi ) = omega_0 sqrt( theta ) IS (4.13) ITSELF, which prints
	 * omega^2 = gamma( Omega^2/R0^2 ) Rbar T -- the rotation frequency is
	 * whatever the temperature makes it, and that is the physical content of the
	 * whole section.
	 *
	 * WHAT THE CANCELLATION COSTS, WHICH IS WHAT MAKES THIS A SHARP TEST. In
	 * d2p/dpsi2 the surviving terms are n_s0'' T_s + 2 n_s0' T_s' + n_s0 T_s''
	 * per species, and the first two are individually O( 1 ) and cancel exactly:
	 *
	 *     n_s0'' T_s = -2 sigma c_s tau_s( beta - sigma nbar )/theta^2
	 *     2 n_s0' T_s' = +2 sigma c_s tau_s( beta - sigma nbar )/theta^2
	 *
	 * so beta != sigma nbar is a requirement rather than a convenience -- at
	 * equality both terms vanish and the test is vacuous.
	 * theCancellationIsNotVacuous measures the magnitude and asserts it.
	 */

	double const speciesTauOne = 0.2;    // T_i = tau_1 theta( psi )
	double const speciesTauTwo = 0.1;    // T_e = tau_2 theta( psi )
	double const temperatureSlope = 0.5; // sigma
	double const densityOffset = 12.0;   // nbar
	double const ionMass = 1.0;
	double const electronMass = 1.0/1836.0;
	double const normalisedMu0 = 1.0;

	/// theta( psi ) = 1 + sigma psi, and tau times it: T_s( psi ).
	class AffineProfile : public meq::Profile
	{
		public:
			AffineProfile( double value, double slope ) : v( value ), s( slope ) {}

			double operator()( double psi ) const override { return v + s*psi; }
			double prime( double ) const override { return s; }
			double doublePrime( double ) const override { return 0.0; }

		private:
			double v, s;
	};

	class ConstantValueProfile : public meq::Profile
	{
		public:
			explicit ConstantValueProfile( double value ) : v( value ) {}

			double operator()( double ) const override { return v; }
			double prime( double ) const override { return 0.0; }
			double doublePrime( double ) const override { return 0.0; }

		private:
			double v;
	};

	/// n_s0( psi ) = c ( nbar + beta psi )/( 1 + sigma psi ), exact at both
	/// derivative levels. Not a ratio of two Profiles, because a general
	/// quotient rule is a second thing to get wrong and this one is three lines.
	class RatioProfile : public meq::Profile
	{
		public:
			RatioProfile( double scale, double offset, double slope, double sigma )
				: c( scale ), n0( offset ), b( slope ), s( sigma ) {}

			double operator()( double psi ) const override
			{
				return c*( n0 + b*psi )/( 1.0 + s*psi );
			}

			double prime( double psi ) const override
			{
				double const d = 1.0 + s*psi;
				return c*( b - s*n0 )/( d*d );
			}

			double doublePrime( double psi ) const override
			{
				double const d = 1.0 + s*psi;
				return -2.0*c*s*( b - s*n0 )/( d*d*d );
			}

		private:
			double c, n0, b, s;
	};

	/// omega( psi ) = omega_0 sqrt( 1 + sigma psi ), which is (4.13).
	class RootProfile : public meq::Profile
	{
		public:
			RootProfile( double amplitude, double sigma ) : a( amplitude ), s( sigma ) {}

			double operator()( double psi ) const override
			{
				return a*std::sqrt( 1.0 + s*psi );
			}

			double prime( double psi ) const override
			{
				return 0.5*a*s/std::sqrt( 1.0 + s*psi );
			}

			double doublePrime( double psi ) const override
			{
				double const d = 1.0 + s*psi;
				return -0.25*a*s*s/( d*std::sqrt( d ) );
			}

		private:
			double a, s;
	};

	/// The exponent coefficient (4.7) fixes, C = m/R0^2, and the reference
	/// pressure gradient that matches (4.10)'s amplitude.
	double sectionFourC( Equilibrium const &eq )
	{
		return eq.getMachSquared()/( eq.getMajorRadius()*eq.getMajorRadius() );
	}

	double sectionFourPressureSlope( Equilibrium const &eq )
	{
		double const r0Sq = eq.getMajorRadius()*eq.getMajorRadius();
		return eq.getPressure()/( normalisedMu0*r0Sq*r0Sq )
		       *std::exp( 0.5*eq.getMachSquared() );
	}

	/// meq::RotatingSource configured to be the fixture, with R_ref = R0 and
	/// mu0 = 1. Every constant is derived from the equilibrium rather than
	/// tabulated, so the mapping above is stated as code and cannot drift from
	/// the fixture it claims to reproduce.
	std::shared_ptr<meq::RotatingSource const> matching( Equilibrium const &eq,
	                                                     bool breakTheRatio = false )
	{
		double const tauSum = speciesTauOne + speciesTauTwo;
		double const massSum = ionMass + electronMass;

		// C = omega_0^2 ( m_1 + m_2 )/( tau_1 + tau_2 ), the thetas cancelling.
		double const omegaAmplitude = std::sqrt( sectionFourC( eq )*tauSum/massSum );

		// p_0( psi ) = ( c tau_1 + c tau_2 ) N( psi ) with c = 1, so
		// p_0' = tauSum * beta.
		double const beta = sectionFourPressureSlope( eq )/tauSum;

		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = ionMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature
			= std::make_shared<AffineProfile const>( speciesTauOne,
			                                         speciesTauOne*temperatureSlope );
		species[ 0 ].density
			= std::make_shared<RatioProfile const>( 1.0, densityOffset, beta,
			                                        temperatureSlope );

		species[ 1 ].mass = electronMass;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature
			= std::make_shared<AffineProfile const>( speciesTauTwo,
			                                         speciesTauTwo*temperatureSlope );
		// Z_1 c_1 + Z_2 c_2 = 0 with Z = +-1 means c_2 = c_1, so the two
		// reference densities are the same function. neutralisingDensity() would
		// build it too; writing it out keeps the parallel with species 0 visible.
		species[ 1 ].density
			= std::make_shared<RatioProfile const>( 1.0, densityOffset, beta,
			                                        temperatureSlope );

		// THE CONTROL. Breaking (4.7) means letting omega^2/T vary, which is done
		// here by giving omega a psi-dependence theta does not cancel. C is then
		// a function of psi, the source acquires a psi-dependence, and the closed
		// form stops describing it -- which is what says the agreement below is
		// (4.7) doing work rather than an identity.
		std::shared_ptr<meq::Profile const> omega;
		if ( eq.getMachSquared() == 0.0 && !breakTheRatio )
			omega = nullptr;
		else if ( breakTheRatio )
			omega = std::make_shared<AffineProfile const>( omegaAmplitude,
			                                              0.5*omegaAmplitude );
		else
			omega = std::make_shared<RootProfile const>( omegaAmplitude,
			                                             temperatureSlope );

		double const r0Sq = eq.getMajorRadius()*eq.getMajorRadius();

		return std::make_shared<meq::RotatingSource const>(
			species, omega,
			std::make_shared<ConstantValueProfile const>( eq.getCurrent()/r0Sq ),
			eq.getMajorRadius(), normalisedMu0 );
	}

	/// An Equilibrium for meq::tests::checkOrder whose psi and flux are the
	/// fixture's closed form and whose F comes from meq::RotatingSource. It is
	/// what makes the rate study a statement about the PRODUCTION source rather
	/// than about the fixture's own f().
	class SourceDrivenEquilibrium
	{
		public:
			explicit SourceDrivenEquilibrium( Equilibrium const &eqIn )
				: eq( eqIn ), source( matching( eqIn ) )
			{
			}

			double psi( double radius, double z ) const { return eq.psi( radius, z ); }

			void flux( double radius, double z, double &qR, double &qZ ) const
			{
				eq.flux( radius, z, qR, qZ );
			}

			double f( double radius, double z, double psiValue ) const
			{
				return source->f( radius, z, psiValue );
			}

			double dFdPsi( double radius, double z, double psiValue ) const
			{
				return source->dFdPsi( radius, z, psiValue );
			}

		private:
			Equilibrium eq;
			std::shared_ptr<meq::RotatingSource const> source;
	};

}

/// The benchmark before the solver: Delta*( psi ) must equal -F, or everything
/// measured against it is measured against the wrong thing.
///
/// THIS IS THE GUARD THE WHOLE FILE RESTS ON. It catches a mistyped term in
/// (3.16) or (4.16), a sign slip between the paper's L and MEQ's Delta*, and --
/// the trap the fixture's own header opens with -- section 3's POLYTROPE
/// closure transcribed in place of section 4's isothermal one, whose g_S is a
/// power law where g_T is an exponential. Without it a wrong fixture would
/// produce a perfect k+1 table for somebody else's equilibrium.
BOOST_AUTO_TEST_CASE( theMaschkePerrinSourceMatchesTheOperator )
{
	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::withPoloidalCurrent() } )
	{
		double worst = 0.0;
		overTheBox( [ &eq, &worst ]( double radius, double z )
		{
			double const deltaStar = eq.deltaStarFD( radius, z );
			double const minusF = -eq.f( radius, z, 0.0 );
			worst = std::max( worst, std::abs( deltaStar - minusF ) );

			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at m = " << eq.getMachSquared() << ", M = " << eq.getCurrent()
			            << ", ( " << radius << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		} );

		std::printf( "  Delta* check at m = %.1f, M = %.1f: worst %10.3e\n",
		             eq.getMachSquared(), eq.getCurrent(), worst );
	}
	std::fflush( stdout );
}

/// gradPsi() is differentiated by hand and the solver's flux error is measured
/// against it, so a slip there would move every L2( q ) in this file without
/// moving a single rate. deltaStarFD() cannot see it: that is built on psi()
/// alone, deliberately.
BOOST_AUTO_TEST_CASE( theGradientsMatchFiniteDifferences )
{
	double const h = 1.0e-6;

	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::withPoloidalCurrent() } )
	{
		overTheBox( [ &eq, h ]( double radius, double z )
		{
			double analyticR, analyticZ;
			eq.gradPsi( radius, z, analyticR, analyticZ );

			double const differencedR = ( eq.psi( radius + h, z ) - eq.psi( radius - h, z ) )
			                            /( 2.0*h );
			double const differencedZ = ( eq.psi( radius, z + h ) - eq.psi( radius, z - h ) )
			                            /( 2.0*h );

			BOOST_TEST( std::abs( analyticR - differencedR ) < 1.0e-8,
			            "at m = " << eq.getMachSquared() << ", ( " << radius << ", " << z
			            << " ): d_r psi = " << analyticR
			            << " but the difference is " << differencedR );
			BOOST_TEST( std::abs( analyticZ - differencedZ ) < 1.0e-8,
			            "at m = " << eq.getMachSquared() << ", ( " << radius << ", " << z
			            << " ): d_z psi = " << analyticZ
			            << " but the difference is " << differencedZ );
		} );

		// And the flux really is grad_bar( psi )/R, since that is what the
		// solver's fluxError() is handed.
		double qR, qZ, gR, gZ;
		eq.flux( 1.1, 0.2, qR, qZ );
		eq.gradPsi( 1.1, 0.2, gR, gZ );
		BOOST_TEST( std::abs( qR - gR/1.1 ) < 1.0e-15 );
		BOOST_TEST( std::abs( qZ - gZ/1.1 ) < 1.0e-15 );
	}
}

/// (4.18) and (4.19), the paper's own statements about the geometry, checked
/// against psi's own derivatives.
///
/// THIS IS THE ONLY THING IN THE FILE THAT CAN SEE A WRONG C OR A WRONG eps_a,
/// and it is here for the reason Soloviev.hpp records about its own
/// coefficients: those constants multiply Delta*-HARMONIC terms, so a wrong one
/// leaves F, Delta*( psi ) and every convergence rate exact and changes only
/// which equilibrium is being solved. Verifying the solution by substitution
/// into (4.10) -- which is what the Delta* scan above does -- cannot reach them
/// either, for the same reason.
///
/// (4.19) is the sharper of the two because it is not the condition C was
/// derived from: it is a separately published closed form for a quantity
/// recoverable from psi's Hessian, so agreement is two routes to one number
/// rather than a check of a formula against itself.
BOOST_AUTO_TEST_CASE( theAxisConditionsAreThePapersOwn )
{
	double const h = 1.0e-4;

	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::withPoloidalCurrent() } )
	{
		double const ra = eq.getAxisRadius()*eq.getMajorRadius();

		// (4.18): C is fixed by dpsi/dR = 0 at ( X = 0, R = R_a ). dpsi/dX
		// vanishes there for free, psi being even in z.
		double dPsiDr, dPsiDz;
		eq.gradPsi( ra, 0.0, dPsiDr, dPsiDz );

		BOOST_TEST( std::abs( dPsiDr ) < 1.0e-14,
		            "at m = " << eq.getMachSquared() << ": (4.18)'s C leaves "
		            "d_r psi = " << dPsiDr << " at the nominated axis R = " << ra
		            << ", so R_a is not an extremum and it is not the magnetic axis" );
		BOOST_TEST( std::abs( dPsiDz ) < 1.0e-15 );

		// (4.19): ( b/a )^2 = d2psi/dR2 / d2psi/dX2 at the axis, both taken by
		// central differences of psi() so that the check does not go through the
		// hand-derived gradients.
		double const psiRR = ( eq.psi( ra + h, 0.0 ) - 2.0*eq.psi( ra, 0.0 )
		                       + eq.psi( ra - h, 0.0 ) )/( h*h );
		double const psiZZ = ( eq.psi( ra, h ) - 2.0*eq.psi( ra, 0.0 )
		                       + eq.psi( ra, -h ) )/( h*h );
		double const measured = psiRR/psiZZ;
		double const printed = eq.axisEllipticitySquared();

		std::printf( "  m = %.1f, M = %.1f: (b/a)^2 measured %12.9f, (4.19) %12.9f\n",
		             eq.getMachSquared(), eq.getCurrent(), measured, printed );

		// Both second derivatives are negative -- the axis is a MAXIMUM, F being
		// positive here where Soloviev.hpp's is negative.
		BOOST_TEST( psiRR < 0.0 );
		BOOST_TEST( psiZZ < 0.0 );

		if ( eq.getCurrent() == 0.0 )
		{
			BOOST_TEST( std::abs( measured - printed ) < 1.0e-7*std::max( 1.0, printed ),
			            "at m = " << eq.getMachSquared() << ": the axis ellipticity "
			            "is " << measured << " where (4.19) prints " << printed );
		}
		else
		{
			// (4.19) IS DERIVED UNDER M = 0 AND IS WRONG WITHOUT IT, and this
			// branch is what keeps that from being a footnote. The M term is
			// -M X^2/( 2 R0^2 ), so it moves d2psi/dX2 and leaves d2psi/dR2
			// alone; measured, 1.767 against (4.19)'s 2.297, 23% apart. A test
			// that applied (4.19) to every configuration would fail here, and
			// one that quietly skipped M != 0 would not record why.
			BOOST_TEST( std::abs( measured - printed ) > 0.1*printed,
			            "at M = " << eq.getCurrent() << " the axis ellipticity is "
			            << measured << " and (4.19) prints " << printed
			            << "; they agree, so either M has stopped reaching "
			               "d2psi/dX2 or (4.19) has been generalised and this "
			               "control has stopped being a control" );
		}
	}
	std::fflush( stdout );
}

/// (4.16) as printed is 0/0 as m -> 0 -- a 1/m^2 prefactor on a brace that
/// vanishes like m^2 -- and the fixture evaluates an algebraically equivalent
/// form in which m has cancelled. This is the test of that, and the middle part
/// is the control that makes the other two mean anything.
BOOST_AUTO_TEST_CASE( theSmallMachLimitIsContinuous )
{
	Equilibrium const stationary = Equilibrium::stationary();
	double const p = stationary.getPressure();
	double const radius0 = stationary.getMajorRadius();

	// PART 1: the paper's own limit under (4.16), g_T -> -( P/8 )( R/R0 )^4 as
	// Omega -> 0. The fixture must BE it at m = 0, not merely approach it.
	// psi at m = 0 minus its harmonic part is exactly that, and the harmonic
	// part is what a second fixture with the same C, eps_a and F_0 but P = 0
	// would give -- so it is read off by subtraction rather than reimplemented.
	{
		double worst = 0.0;
		overTheBox( [ & ]( double radius, double z )
		{
			// The harmonic part written out from (3.16), so that what is
			// compared is psi minus (3.16) against the paper's printed limit
			// rather than the fixture against itself.
			double const r0Sq = radius0*radius0;
			double const c = stationary.coefficientC();
			double const harmonic = c*p*radius*radius/r0Sq
			                        + ( stationary.getEllipticity() - 1.0 )*p*radius*radius
			                          *( z*z - 0.25*radius*radius )/( 4.0*r0Sq*r0Sq )
			                        + stationary.getFluxOffset();
			double const limit = -p*std::pow( radius/radius0, 4 )/8.0;
			worst = std::max( worst,
			                  std::abs( stationary.psi( radius, z ) - harmonic - limit ) );
		} );
		BOOST_TEST( worst < 1.0e-15,
		            "at m = 0 the fixture differs from the paper's own -( P/8 )"
		            "( R/R0 )^4 limit by " << worst );
		std::printf( "  (4.16) at m = 0 against -(P/8)(R/R0)^4: %10.3e\n", worst );
	}

	// PART 2, THE CONTROL. (4.16) as printed is not merely less accurate for
	// small m, it is useless -- and a test that only checked the fixture against
	// itself would pass with the printed form in place. MEASURED at R = 1.3,
	// comparing the printed form against the fixture's AT THE SAME m:
	//
	//     m        printed form, relative error against the fixture
	//     1e-2            1.1e-12
	//     1e-4            3.1e-08
	//     1e-6            9.8e-05
	//     1e-8            1.0e+00     <-- returns exactly zero
	//     1e-10           1.0e+00
	//
	// Asserted at 1e-2 and 1e-8 together, so the test states both halves: the
	// printed form is fine where there is no cancellation and useless where
	// there is. The first half matters as much as the second -- it is what says
	// the two expressions really are the same algebra, so that the second half
	// is measuring precision loss and not a transcription error in either.
	{
		double const radius = 1.3;

		auto particularOf = [ & ]( Equilibrium const &eq )
		{
			// psi minus the harmonic part, at z = 0.
			double const r0Sq = radius0*radius0;
			double const harmonic = eq.coefficientC()*p*radius*radius/r0Sq
			                        + ( eq.getEllipticity() - 1.0 )*p*radius*radius
			                          *( -0.25*radius*radius )/( 4.0*r0Sq*r0Sq )
			                        + eq.getFluxOffset();
			return eq.psi( radius, 0.0 ) - harmonic;
		};

		double const limit = -p*std::pow( radius/radius0, 4 )/8.0;

		double const mildMach = 1.0e-2;
		Equilibrium const mild( radius0, mildMach, p, 0.0, 0.0, 1.0, 0.0 );
		double const stableMild = particularOf( mild );
		double const printedMild = naiveParticularAsPrinted( radius, mildMach, p, radius0 );
		BOOST_TEST( std::abs( printedMild - stableMild )/std::abs( stableMild ) < 1.0e-9,
		            "the two forms disagree by "
		            << std::abs( printedMild - stableMild )/std::abs( stableMild )
		            << " at m = " << mildMach
		            << ", where there is no cancellation to blame -- so they are "
		               "not the same algebra and one of them is mistyped" );

		double const smallMach = 1.0e-8;
		Equilibrium const small( radius0, smallMach, p, 0.0, 0.0, 1.0, 0.0 );
		double const stableSmall = particularOf( small );
		double const printedSmall = naiveParticularAsPrinted( radius, smallMach, p, radius0 );
		BOOST_TEST( std::abs( printedSmall - stableSmall )/std::abs( stableSmall ) > 1.0e-2,
		            "(4.16) as printed is accurate to "
		            << std::abs( printedSmall - stableSmall )/std::abs( stableSmall )
		            << " at m = " << smallMach
		            << ", so the stable form is not buying anything and this "
		               "control has stopped being a control" );

		// And it is the printed form that is wrong, not the fixture.
		BOOST_TEST( std::abs( stableSmall - limit )/std::abs( limit ) < 1.0e-7,
		            "the fixture is out by "
		            << std::abs( stableSmall - limit )/std::abs( limit )
		            << " at m = " << smallMach );
	}

	// PART 3: the crossover. g1() and g2() change branch at
	// |w| = seriesThreshold(), and a series that disagreed with the closed form
	// there would put a step into psi and into the flux -- small enough to
	// survive every other assertion in this file and large enough to spoil a
	// convergence rate. Straddle it in m at fixed R and require the jump to
	// scale like the perturbation, which is what "no step" means. A genuine
	// branch mismatch does not scale with delta at all, so it fails this at the
	// smallest one however loose the constant is.
	{
		double const radius = 1.3;
		double const crossover = 2.0*Equilibrium::seriesThreshold()*radius0*radius0/( radius*radius );

		for ( double delta : { 1.0e-6, 1.0e-9, 1.0e-12 } )
		{
			Equilibrium const below( radius0, crossover - delta, p, 0.0, 0.0, 1.0, 0.0 );
			Equilibrium const above( radius0, crossover + delta, p, 0.0, 0.0, 1.0, 0.0 );

			double belowR, belowZ, aboveR, aboveZ;
			below.gradPsi( radius, 0.3, belowR, belowZ );
			above.gradPsi( radius, 0.3, aboveR, aboveZ );

			double const jumpPsi = std::abs( above.psi( radius, 0.3 ) - below.psi( radius, 0.3 ) );
			double const jumpFlux = std::abs( aboveR - belowR );

			BOOST_TEST( jumpPsi < 10.0*delta,
			            "psi jumps by " << jumpPsi << " across the series crossover "
			            "at m = " << crossover << " for a perturbation of " << delta );
			BOOST_TEST( jumpFlux < 10.0*delta,
			            "d_r psi jumps by " << jumpFlux << " across the series "
			            "crossover at m = " << crossover << " for a perturbation of "
			            << delta );
		}
	}
	std::fflush( stdout );
}

/// psi = 0 is a closed curve strictly inside the benchmark box, which is what
/// makes this a physically sensible equilibrium rather than merely a solution of
/// the equation -- and what makes the Dirichlet data non-homogeneous, so that
/// the study below exercises the boundary path a homogeneous one would not.
///
/// It is also the only assertion that can see a wrong F_0, which changes no
/// rate and no source and only which surface is psi = 0.
BOOST_AUTO_TEST_CASE( theZeroContourIsClosedInsideTheBox )
{
	meq::tests::Rectangle const box = meq::tests::standardBox();

	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::withPoloidalCurrent() } )
	{
		// F_0 is pinned by psi( 1.2, 0 ) = 0, the outer equatorial point of the
		// plasma. Only ONE geometric condition is available here, unlike
		// RotatingSoloviev.hpp's four: eps_a, R_a and m fix the shape outright
		// and F_0 chooses only which level set is psi = 0.
		BOOST_TEST( std::abs( eq.psi( 1.2, 0.0 ) ) < 1.0e-15,
		            "m = " << eq.getMachSquared() << ": psi at the outer equatorial "
		            "point is " << eq.psi( 1.2, 0.0 ) );

		double worstOnBoundary = -1.0e300;
		for ( int i = 0; i <= 400; ++i )
		{
			double const s = static_cast<double>( i )/400.0;
			double const radius = box.minRadius + s*box.width();
			double const z = box.zMin + s*box.height();
			for ( double value : { eq.psi( radius, box.zMin ), eq.psi( radius, box.zMax ),
			                       eq.psi( box.minRadius, z ), eq.psi( box.maxRadius, z ) } )
				worstOnBoundary = std::max( worstOnBoundary, value );
		}

		std::printf( "  m = %.1f, M = %.1f: psi at the axis %10.6f, worst on the "
		             "box boundary %10.6f\n", eq.getMachSquared(), eq.getCurrent(),
		             eq.psi( eq.getAxisRadius()*eq.getMajorRadius(), 0.0 ),
		             worstOnBoundary );

		// Measured: -0.0270 for stationary and -0.0392 for the other two.
		BOOST_TEST( worstOnBoundary < -1.0e-2,
		            "m = " << eq.getMachSquared() << ": psi reaches " << worstOnBoundary
		            << " on the boundary of the benchmark box, so the psi = 0 "
		               "contour is not closed inside it" );
		BOOST_TEST( eq.psi( eq.getAxisRadius()*eq.getMajorRadius(), 0.0 ) > 1.0e-2,
		            "m = " << eq.getMachSquared() << ": psi at the axis is only "
		            << eq.psi( eq.getAxisRadius()*eq.getMajorRadius(), 0.0 ) );
	}
	std::fflush( stdout );
}

/*
 * ------------------------------------------------------------------------
 * WHAT THIS FIXTURE IS FOR: meq::RotatingSource with psi-DEPENDENT PROFILES
 * ------------------------------------------------------------------------
 */

/// meq::RotatingSource, built on section 4's profiles, reproduces the closed
/// form -- and its dF/dpsi is zero, which is where five varying profiles have to
/// cancel.
///
/// EVERY OTHER CLOSED-FORM CHECK OF meq::RotatingSource IN THIS TREE RUNS AT
/// CONSTANT T AND CONSTANT omega. RotatingSourceConvergence.cpp's species carry
/// ConstantMassProfile for both, because Li & Zhu's rotating Solov'ev case has
/// T_0 and Omega_0 constant and no closed form survives otherwise. What checks
/// the varying case today is a central difference of MEQ's own f(), in
/// RotatingSourceTests.cpp -- and CLAUDE.md records that this cannot see a term
/// missing from both f() and dFdPsi(), because both are evaluated at the same
/// profiles and agree with each other however wrong they are.
///
/// (4.7) is what makes a varying case reachable at all: it constrains
/// omega^2/( Rbar T ) and nothing else, so T( psi ) stays free and omega( psi )
/// follows it through (4.13).
BOOST_AUTO_TEST_CASE( theRotatingSourceReproducesTheClosedForm )
{
	meq::tests::Rectangle const box = meq::tests::standardBox();

	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::withPoloidalCurrent() } )
	{
		std::shared_ptr<meq::RotatingSource const> const source = matching( eq );

		double worstF = 0.0;
		double worstJacobian = 0.0;
		double worstExponentDrift = 0.0;

		for ( double radius = box.minRadius; radius <= box.maxRadius + 1.0e-12; radius += 0.05 )
		{
			for ( double z : { -0.4, 0.0, 0.4 } )
			{
				// Deliberately NOT evaluated at psi( R, z ). F is independent of
				// psi here, so sweeping psi across the whole range the solve
				// visits is what makes the check say so.
				for ( double psiValue : { -0.7, -0.3, 0.0, 0.03, 0.06 } )
				{
					double const expected = eq.f( radius, z, psiValue );
					double const actual = source->f( radius, z, psiValue );
					double const scale = std::max( 1.0, std::abs( expected ) );

					worstF = std::max( worstF, std::abs( actual - expected )/scale );
					worstJacobian = std::max( worstJacobian,
					                          std::abs( source->dFdPsi( radius, z, psiValue ) ) );

					BOOST_TEST( std::abs( actual - expected ) <= 1.0e-12*scale,
					            "m = " << eq.getMachSquared() << " at ( " << radius << ", "
					            << z << " ), psi = " << psiValue
					            << ": meq::RotatingSource gives F = " << actual
					            << " where (4.10) gives " << expected );
				}

				// (4.7) itself, as the source sees it: the shared exponent must
				// not move with psi. It is the hypothesis the whole section
				// rests on, and it is checkable directly rather than only
				// through its consequences.
				double const reference = source->densityExponent( 0, radius, 0.0 );
				for ( double psiValue : { -0.7, 0.06 } )
					worstExponentDrift
						= std::max( worstExponentDrift,
						            std::abs( source->densityExponent( 0, radius, psiValue )
						                      - reference ) );
			}
		}

		std::printf( "  m = %.1f, M = %.1f: worst |F - (4.10)|/scale %10.3e, "
		             "worst |dF/dpsi| %10.3e, worst exponent drift %10.3e\n",
		             eq.getMachSquared(), eq.getCurrent(), worstF, worstJacobian,
		             worstExponentDrift );

		// dF/dpsi is mu0 R^2 d2p/dpsi2 + ( g g' )'. Every term is individually
		// O( 1 ) -- see theCancellationIsNotVacuous -- and they cancel, so this
		// is a round-off assertion on an O( 1 ) cancellation and not a triviality.
		BOOST_TEST( worstJacobian < 1.0e-12,
		            "m = " << eq.getMachSquared() << ": dF/dpsi reaches "
		            << worstJacobian << " where (4.9)'s linear p_T makes it zero" );

		// (4.7) as the source sees it.
		BOOST_TEST( worstExponentDrift < 1.0e-14,
		            "m = " << eq.getMachSquared() << ": the shared density exponent "
		            "moves by " << worstExponentDrift << " with psi, so C is not "
		            "constant and (4.7) is not being satisfied" );
	}
	std::fflush( stdout );
}

/// The profiles really do vary, and the terms that cancel in dF/dpsi really are
/// O( 1 ).
///
/// WITHOUT THIS THE TEST ABOVE IS SATISFIED BY A PLASMA WITH CONSTANT PROFILES,
/// which is exactly what the rest of the suite already covers. It asserts the
/// four derivatives that must be non-zero for the fixture to be doing what it
/// claims, and it measures the size of the cancellation that makes dF/dpsi
/// vanish -- n_s0'' T_s + 2 n_s0' T_s' summing to zero from two terms of
/// magnitude 2 sigma c tau ( beta - sigma nbar )/theta^2 each.
BOOST_AUTO_TEST_CASE( theCancellationIsNotVacuous )
{
	Equilibrium const eq = Equilibrium::rotating();
	std::shared_ptr<meq::RotatingSource const> const source = matching( eq );

	std::vector<meq::Species> const &species = source->species();
	BOOST_TEST_REQUIRE( species.size() == 2u );

	double const psiValue = 0.0;

	double const tPrime = species[ 0 ].temperature->prime( psiValue );
	double const nPrime = species[ 0 ].density->prime( psiValue );
	double const nDoublePrime = species[ 0 ].density->doublePrime( psiValue );
	double const omegaPrime = source->omega()->prime( psiValue );

	std::printf( "  T'(0) = %.6f, n'(0) = %.6f, n''(0) = %.6f, omega'(0) = %.6f\n",
	             tPrime, nPrime, nDoublePrime, omegaPrime );

	// Measured: T' = 0.1, n' = 4.991475, n'' = -4.991475, omega' = 0.136893.
	BOOST_TEST( std::abs( tPrime ) > 1.0e-2,
	            "the temperature is constant in psi, so this fixture is testing "
	            "nothing RotatingSourceConvergence.cpp does not already test" );
	BOOST_TEST( std::abs( omegaPrime ) > 1.0e-2,
	            "omega is constant in psi, so (4.13) is not being exercised" );
	BOOST_TEST( std::abs( nPrime ) > 1.0e-2 );
	BOOST_TEST( std::abs( nDoublePrime ) > 1.0e-2,
	            "the reference density is affine, so n'' T + 2 n' T' is a cancellation "
	            "of one term against zero rather than of two O( 1 ) terms" );

	// The two terms that cancel, per species, and the largest F they sit inside.
	double largestPair = 0.0;
	for ( std::size_t s = 0; s < species.size(); ++s )
	{
		double const a = std::abs( species[ s ].density->doublePrime( psiValue )
		                           *( *species[ s ].temperature )( psiValue ) );
		double const b = std::abs( 2.0*species[ s ].density->prime( psiValue )
		                           *species[ s ].temperature->prime( psiValue ) );
		largestPair = std::max( largestPair, std::max( a, b ) );
	}
	double const scale = 1.4*1.4*largestPair;

	std::printf( "  the cancelling terms in mu0 R^2 d2p/dpsi2 are %.6f at R = 1.4, "
	             "against dF/dpsi = %.3e\n", scale,
	             std::abs( source->dFdPsi( 1.4, 0.0, psiValue ) ) );

	// Measured: 1.956658 against a dF/dpsi of 0.000e+00.
	BOOST_TEST( scale > 0.1,
	            "the terms cancelling in d2p/dpsi2 are only " << scale
	            << ", so beta is too close to sigma nbar and the zero Jacobian "
	               "asserted above is vacuous" );
	BOOST_TEST( std::abs( source->dFdPsi( 1.4, 0.0, psiValue ) ) < 1.0e-12*scale );
}

/// The control: break (4.7) and the agreement goes away.
///
/// Without it, theRotatingSourceReproducesTheClosedForm establishes only that
/// two expressions agree, not that (4.7) is what makes them agree -- and a
/// source that ignored omega's psi-dependence entirely would pass it.
BOOST_AUTO_TEST_CASE( breakingTheConstantRatioBreaksTheAgreement )
{
	Equilibrium const eq = Equilibrium::rotating();
	std::shared_ptr<meq::RotatingSource const> const conforming = matching( eq );
	std::shared_ptr<meq::RotatingSource const> const broken = matching( eq, true );

	double worst = 0.0;
	double worstJacobian = 0.0;
	for ( double radius : { 0.6, 1.0, 1.4 } )
	{
		for ( double psiValue : { -0.5, 0.0, 0.05 } )
		{
			double const expected = eq.f( radius, 0.0, psiValue );
			worst = std::max( worst, std::abs( broken->f( radius, 0.0, psiValue ) - expected )
			                         /std::max( 1.0, std::abs( expected ) ) );
			worstJacobian = std::max( worstJacobian,
			                          std::abs( broken->dFdPsi( radius, 0.0, psiValue ) ) );
		}
	}

	std::printf( "  (4.7) broken: worst |F - (4.10)|/scale %10.3e, "
	             "worst |dF/dpsi| %10.3e\n", worst, worstJacobian );

	// Measured: 2.894e-01 and 5.769e+00, against the conforming source's
	// 4.255e-16 and 2.611e-15.
	BOOST_TEST( worst > 1.0e-3,
	            "omega^2/T was made psi-dependent and F still matches (4.10) to "
	            << worst << ", so (4.7) is not what the agreement rests on" );
	BOOST_TEST( worstJacobian > 1.0e-3,
	            "omega^2/T was made psi-dependent and dF/dpsi is still " << worstJacobian
	            << ", so the source is not seeing omega's psi-dependence at all" );

	// And the conforming one is unmoved, so the two differ in (4.7) and in
	// nothing else.
	BOOST_TEST( std::abs( conforming->dFdPsi( 1.4, 0.0, 0.0 ) ) < 1.0e-12 );
}

/*
 * THE RATES. meq::tests::checkOrder asserts k+1 in psi AND in q over the four
 * dyadic meshes, plus an absolute ceiling on the finest -- which is the only
 * assertion here that can see a solution converging beautifully to the wrong
 * function.
 *
 * THE SOURCE IS meq::RotatingSource, NOT THE FIXTURE'S OWN f(), and that is the
 * whole reason these three cases exist rather than being a restatement of
 * RotatingSolovievConvergence.cpp. What is being measured is the PRODUCTION
 * source, on five psi-dependent profiles, driving the solver to a published
 * exact solution. A source that agreed pointwise but could not be assembled --
 * or one whose Jacobian mass term appeared from somewhere it should not have --
 * would pass everything above and fail here.
 *
 * The slack is 0.15 rather than the harness default of 0.2, matching
 * SolovievConvergence.cpp: this problem is linear and has an exact solution, so
 * there is no reason to allow the width the non-linear studies need.
 *
 * The Newton column of the printed table must read 1 at every mesh and every
 * order. dF/dpsi is identically zero, so the system is affine and step one must
 * be exact; anything else means the Jacobian has acquired a mass term.
 */

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Ceilings at roughly 3x the measured finest-mesh error:
	// psi 1.516748e-04, q 2.543225e-04 at h = 0.025.
	meq::tests::checkOrder( SourceDrivenEquilibrium( Equilibrium::withPoloidalCurrent() ),
	                        "Maschke & Perrin section 4, m = 1, M = 0.3", 1,
	                        4.6e-4, 7.7e-4, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured: psi 9.150755e-07, q 1.262753e-06.
	meq::tests::checkOrder( SourceDrivenEquilibrium( Equilibrium::withPoloidalCurrent() ),
	                        "Maschke & Perrin section 4, m = 1, M = 0.3", 2,
	                        2.8e-6, 3.8e-6, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured: psi 3.561544e-09, q 5.318085e-09.
	meq::tests::checkOrder( SourceDrivenEquilibrium( Equilibrium::withPoloidalCurrent() ),
	                        "Maschke & Perrin section 4, m = 1, M = 0.3", 3,
	                        1.1e-8, 1.6e-8, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}

/// The m = 0 control, all the way through the solver. meq::RotatingSource is
/// handed a NULL omega there, which is a different branch of the closure, and
/// the fixture's g_T is its G2( 0 ) = 1/2 series branch -- so this is not merely
/// a weaker version of the rotating study.
BOOST_AUTO_TEST_CASE( theStationaryControlConvergesToo )
{
	// Measured: psi 4.318208e-07, q 5.362945e-07.
	meq::tests::checkOrder( SourceDrivenEquilibrium( Equilibrium::stationary() ),
	                        "Maschke & Perrin section 4, m = 0", 2,
	                        1.3e-6, 1.7e-6, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}
