#define BOOST_TEST_MODULE VaryingCentrifugalConvergence

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"

#include "analytic/VaryingCentrifugal.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * meq::RotatingSource AT C'( psi ) != 0, AGAINST AN INDEPENDENT CLOSED FORM.
 *
 * C( psi ) = omega^2( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ) is the exponent
 * both species of a two-species rotating plasma share, and C' is non-zero
 * exactly when omega^2/T varies from flux surface to flux surface. It is the
 * term Li & Zhu print with the wrong sign in their (9) -- twice, on the
 * dOmega/dpsi and the dT/dpsi corrections -- and neither of their benchmarks nor
 * Maschke & Perrin's can see it, both holding C constant. tests/analytic/
 * VaryingCentrifugal.hpp is the fixture and its header carries the derivation,
 * the construction and the reason no exact solution can reach this.
 *
 * WHAT THIS FILE ADDS, STATED AGAINST WHAT IS ALREADY THERE, because the gap is
 * narrower than "nothing touches C'" and overstating it would be the same
 * species of error the fixture is guarding against.
 *
 *   RotatingSourceTests.cpp runs at C' != 0 and already checks pressure() and
 *   densityExponent() against an independent closed form there, to 1e-14. So
 *   C ITSELF is pinned. C' and C'' are not: they appear in nothing but dp/dpsi
 *   and d2p/dpsi2, which is to say in nothing but f() and dFdPsi(), and those
 *   are checked today only by CENTRAL DIFFERENCES of quantities MEQ computes,
 *   at 1e-6 to 1e-7. A difference of f() cannot see a term missing from both
 *   f() and dFdPsi().
 *
 *   RotatingNewtonConvergence.cpp also runs at C' != 0 -- its profiles are
 *   unrelated polynomials, so C varies -- and it drives a full manufactured
 *   solve there. But its manufactured remainder subtracts meq::RotatingSource's
 *   OWN f() at the exact solution, so psiExact solves the equation whatever
 *   f() computes. It measures dFdPsi against f, by Newton's order and by
 *   differencing the assembled residual; it cannot measure either against the
 *   equations.
 *
 * So what is new here is one thing and it is the thing the campaign was
 * missing: a SECOND IMPLEMENTATION of (96), (97) and (136), written from the
 * paper in tests/analytic/ and calling no MEQ code, agreeing with
 * meq::RotatingSource pointwise at C' != 0. Measured below at 7.1e-16 in F and
 * 5.3e-16 in dF/dpsi.
 *
 * THE FIXTURE CARRIES TWO ROUTES OF ITS OWN AND THAT IS NOT DECORATION. Two
 * transcriptions of one piece of algebra agreeing says only that the algebra was
 * copied consistently. So VaryingCentrifugal.hpp also solves (97) by BISECTION,
 * with C appearing nowhere, and theFixtureAlgebraFollowsFromQuasineutrality
 * differences that numerically to reach the same dp/dpsi and d2p/dpsi2. A slip
 * shared between MEQ's derivation and the fixture's would survive the pointwise
 * comparison and would not survive this.
 *
 * ORDER OF THE TEST CASES IS LOAD BEARING, as in MaschkePerrinConvergence.cpp.
 * The fixture's own guards come first -- positivity, quasineutrality, the
 * bisection anchor -- because everything after them is measured against a closed
 * form that would otherwise be taken on trust.
 *
 * THE SOLVE AT THE END IS THE SECONDARY ACCEPTANCE AND ITS VALUE IS NARROW.
 * A manufactured solve IS reachable -- the remainder trick works here exactly as
 * it does next door -- but the pointwise comparison above is nine orders sharper
 * than any rate table could be, so the solve is not what closes the gap. What it
 * buys is that the source and its Jacobian are ASSEMBLED and driven through
 * Newton at a C'-dominated Jacobian and a charge-asymmetric species set, neither
 * of which the pointwise sweep touches. Its remainder is built from the
 * FIXTURE's f rather than from MEQ's, which is the one substantive difference
 * from RotatingNewtonConvergence.cpp: here psiExact stops being the exact
 * solution if the two implementations disagree.
 *
 * The domain is meq::tests::standardBox(), and psiExact is
 * ManufacturedNonlinear's shape, so this study sits on the same box and the same
 * flux as NewtonConvergence.cpp's Example 5 and RotatingNewtonConvergence.cpp's
 * rotating one and differs from them in the source alone.
 */

namespace
{

	using Plasma = meq::analytic::VaryingCentrifugalPlasma;
	using meq::tests::Measurement;
	using meq::tests::bestNewtonOrder;
	using meq::tests::newtonOrder;
	using meq::tests::standardBox;

	/// psiExact = sin( kr( R + R_0 ) )cos( kz z ), ManufacturedNonlinear's
	/// parameters. R_0 is a radial offset that places the arch of the sine inside
	/// the domain and is NOT a major radius.
	double const pi = 3.14159265358979323846;
	double const psiR0 = -0.5;
	double const psiKr = 1.15*pi;
	double const psiKz = 1.15;

	/*
	 * ------------------------------------------------------------------
	 * The fixture, as meq::RotatingSource takes it
	 * ------------------------------------------------------------------
	 *
	 * ONE TRANSCRIPTION OF THE PROFILES, TWO CONSUMERS. The adapter below wraps
	 * the fixture's own accessors as meq::Profile, so the production source and
	 * the closed form are handed the SAME six flux functions. That is deliberate
	 * and it is where the independence is and is not: the profiles are shared,
	 * the equation solved on them is not. Writing the profiles out twice -- once
	 * in the fixture and once here -- would buy nothing and would let the two
	 * drift, which is a failure mode with no upside.
	 */
	class AdaptedProfile : public meq::Profile
	{
		public:
			AdaptedProfile( std::function<double( double )> valueIn,
			                std::function<double( double )> primeIn,
			                std::function<double( double )> doublePrimeIn )
				: valueFunction( std::move( valueIn ) ),
				  primeFunction( std::move( primeIn ) ),
				  doublePrimeFunction( std::move( doublePrimeIn ) )
			{
			}

			double operator()( double psi ) const override
			{
				return valueFunction( psi );
			}

			double prime( double psi ) const override
			{
				return primeFunction( psi );
			}

			double doublePrime( double psi ) const override
			{
				return doublePrimeFunction( psi );
			}

		private:
			std::function<double( double )> valueFunction;
			std::function<double( double )> primeFunction;
			std::function<double( double )> doublePrimeFunction;
	};

	std::shared_ptr<meq::Profile const> adapt( std::function<double( double )> value,
	                                           std::function<double( double )> prime,
	                                           std::function<double( double )> doublePrime )
	{
		return std::make_shared<AdaptedProfile const>( std::move( value ),
		                                              std::move( prime ),
		                                              std::move( doublePrime ) );
	}

	/// meq::RotatingSource configured to be the fixture. The plasma is captured
	/// BY VALUE in every lambda, so the returned source owns everything it reads
	/// and outliving the caller's fixture is not a hazard.
	std::shared_ptr<meq::RotatingSource const>
	matching( Plasma const &plasma,
	          meq::RotatingSource::Closure closure = meq::RotatingSource::Closure::Automatic )
	{
		std::vector<meq::Species> species( Plasma::speciesCount() );

		for ( std::size_t s = 0; s < Plasma::speciesCount(); ++s )
		{
			species[ s ].mass = plasma.mass( s );
			species[ s ].charge = plasma.charge( s );
			species[ s ].temperature = adapt(
				[ plasma, s ]( double psi ) { return plasma.temperature( s, psi ); },
				[ plasma, s ]( double psi ) { return plasma.temperaturePrime( s, psi ); },
				[ plasma, s ]( double psi ) { return plasma.temperatureDoublePrime( s, psi ); } );
			species[ s ].density = adapt(
				[ plasma, s ]( double psi ) { return plasma.referenceDensity( s, psi ); },
				[ plasma, s ]( double psi ) { return plasma.referenceDensityPrime( s, psi ); },
				[ plasma, s ]( double psi ) { return plasma.referenceDensityDoublePrime( s, psi ); } );
		}

		std::shared_ptr<meq::Profile const> const omega = adapt(
			[ plasma ]( double psi ) { return plasma.omega( psi ); },
			[ plasma ]( double psi ) { return plasma.omegaPrime( psi ); },
			[ plasma ]( double psi ) { return plasma.omegaDoublePrime( psi ); } );

		// g g' is affine, so its second derivative is zero -- and it is a genuine
		// zero rather than an unimplemented one, which is why it is written out
		// rather than left to a default.
		std::shared_ptr<meq::Profile const> const ggPrime = adapt(
			[]( double psi ) { return Plasma::ggPrime( psi ); },
			[]( double psi ) { return Plasma::ggPrimeDerivative( psi ); },
			[]( double ) { return 0.0; } );

		return std::make_shared<meq::RotatingSource const>(
			species, omega, ggPrime, Plasma::referenceRadius(), Plasma::mu0(), closure );
	}

	/// The sweep every pointwise comparison below runs: the standard box in R,
	/// three heights in z -- F does not depend on z and the sweep is what says so
	/// -- and psi over [ -0.2, 1.0 ], which covers the range the manufactured
	/// solve visits with margin at both ends.
	///
	/// R RUNS TO BOTH EDGES OF THE BOX ON PURPOSE. Every C' term carries
	/// h = ( R^2 - R_ref^2 )/2, which vanishes at R_ref = 1 and changes sign
	/// there; a sweep confined to one side would measure a term of one sign only,
	/// and a sweep placed at R_ref would measure nothing at all. That is the trap
	/// CLAUDE_FLOW.md records from FL-5.
	template<typename Check>
	void overTheSweep( Check check )
	{
		meq::tests::Rectangle const box = standardBox();
		for ( double radius = box.minRadius; radius <= box.maxRadius + 1.0e-12; radius += 0.05 )
		{
			for ( double z : { -0.4, 0.0, 0.4 } )
			{
				for ( double psi = -0.2; psi <= 1.0 + 1.0e-12; psi += 0.05 )
				{
					check( radius, z, psi );
				}
			}
		}
	}

	/// A relative deviation against a floor of one, which is the scaling the rest
	/// of this directory uses. F passes through zero inside the sweep -- it is a
	/// pressure gradient times R^2 plus g g' and there is no reason it should not
	/// -- so a bare relative error is not defined everywhere and a floored one is.
	double deviation( double actual, double expected )
	{
		return std::abs( actual - expected )/std::max( 1.0, std::abs( expected ) );
	}

	/*
	 * ------------------------------------------------------------------
	 * The manufactured solve
	 * ------------------------------------------------------------------
	 *
	 * F_total( R, z, psi ) = rot.f( R, z, psi ) + remainder( R, z ),
	 * remainder( R, z ) = -Delta*( psiExact ) - FIXTURE.f( R, z, psiExact ).
	 *
	 * THE REMAINDER IS THE FIXTURE'S AND NOT THE SOURCE'S, which is the whole
	 * difference from RotatingNewtonConvergence.cpp. There the two rotating terms
	 * cancel identically at psiExact and the exact solution is exact whatever
	 * f() computes; here they cancel only insofar as the two implementations
	 * agree, so a disagreement stops psiExact being a solution and the rate table
	 * degrades. It is a far blunter instrument than the pointwise sweep above --
	 * it could not see a disagreement below the discretisation error -- but it is
	 * the half that goes through the assembly.
	 *
	 * The remainder is a function of ( R, z ) ALONE either way, so
	 * dF_total/dpsi IS meq::RotatingSource::dFdPsi exactly, with nothing added to
	 * it. That is what makes Newton's observed order below a statement about the
	 * production Jacobian.
	 */
	class ManufacturedEquilibrium
	{
		public:
			/// @param jacobianScaleIn  a deliberate error in dF/dpsi, for the
			///        control that calibrates the Newton-order assertion. 1.0 is
			///        the only value a measurement may use; 1.05 is CLAUDE_FLOW.md's
			///        +5% experiment, reproduced here on this fixture rather than
			///        quoted from another one. It scales ONLY the Jacobian, so
			///        the residual -- and therefore the converged answer and every
			///        rate -- is untouched, which is the whole point of the
			///        experiment.
			explicit ManufacturedEquilibrium( Plasma const &plasmaIn = Plasma::standard(),
			                                  double jacobianScaleIn = 1.0 )
				: plasma( plasmaIn ), rot( matching( plasmaIn ) ),
				  jacobianScale( jacobianScaleIn )
			{
			}

			double psi( double radius, double z ) const
			{
				return std::sin( psiKr*( radius + psiR0 ) )*std::cos( psiKz*z );
			}

			void gradPsi( double radius, double z, double &dPsiDr, double &dPsiDz ) const
			{
				dPsiDr =  psiKr*std::cos( psiKr*( radius + psiR0 ) )*std::cos( psiKz*z );
				dPsiDz = -psiKz*std::sin( psiKr*( radius + psiR0 ) )*std::sin( psiKz*z );
			}

			/// The HDG flux q = grad_bar( psi )/R.
			void flux( double radius, double z, double &qR, double &qZ ) const
			{
				gradPsi( radius, z, qR, qZ );
				qR /= radius;
				qZ /= radius;
			}

			/// Delta*( psiExact ) in closed form, from
			/// Delta* := d_rr - ( 1/R )d_r + d_zz.
			double deltaStar( double radius, double z ) const
			{
				return -( psiKr*psiKr + psiKz*psiKz )*psi( radius, z )
				       - ( psiKr/radius )*std::cos( psiKr*( radius + psiR0 ) )*std::cos( psiKz*z );
			}

			/// Delta*( psiExact ) by central differences, arranged as
			/// R d_r( ( 1/R )d_r psi ) + d_zz psi -- the form every fixture in
			/// tests/analytic uses, so that it is independent of the closed form
			/// above rather than a rearrangement of it.
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

			double remainder( double radius, double z ) const
			{
				return -deltaStar( radius, z ) - plasma.f( radius, z, psi( radius, z ) );
			}

			double f( double radius, double z, double psiValue ) const
			{
				return rot->f( radius, z, psiValue ) + remainder( radius, z );
			}

			double dFdPsi( double radius, double z, double psiValue ) const
			{
				return jacobianScale*rot->dFdPsi( radius, z, psiValue );
			}

		private:
			Plasma plasma;
			std::shared_ptr<meq::RotatingSource const> rot;
			double jacobianScale;
	};

}

/*
 * ------------------------------------------------------------------------
 * The fixture's own guards, first
 * ------------------------------------------------------------------------
 */

/// Every shape function stays positive, and omega stays real, over a range far
/// wider than any iterate reaches.
///
/// THIS IS A PRECONDITION AND NOT TIDINESS. The two-species closure divides by
/// Z_1 T_2 - Z_2 T_1 and the fixture takes a square root to get omega from the
/// prescribed C, so a shape function that turns over gives a throw or a NaN --
/// from inside a quadrature loop, which is a poor place to find out. Each
/// quadratic in the fixture has a negative discriminant; this is the assertion
/// that says so, over psi in [ -3, 3 ] rather than over the visited range, so
/// that a later change to the coefficients fails with a sentence.
BOOST_AUTO_TEST_CASE( theShapeFunctionsStayPositiveEverywhere )
{
	Plasma const plasma = Plasma::standard();

	double worstTemperature = 1.0e300;
	double worstDenominator = 1.0e300;
	double worstDensity = 1.0e300;
	double worstCoefficient = 1.0e300;

	for ( int i = -300; i <= 300; ++i )
	{
		double const psi = i/100.0;

		// T_1 is affine and DOES turn over, at psi = -2.857; it is the one shape
		// function that is not positive everywhere, and the margin is what is
		// asserted rather than positivity on the whole line.
		worstTemperature = std::min( worstTemperature, plasma.temperature( 1, psi ) );
		worstDenominator = std::min( worstDenominator, plasma.closureDenominator( psi ) );
		worstCoefficient = std::min( worstCoefficient, plasma.exponentCoefficient( psi ) );
		for ( std::size_t s = 0; s < Plasma::speciesCount(); ++s )
		{
			worstDensity = std::min( worstDensity, plasma.referenceDensity( s, psi ) );
		}

		BOOST_TEST( std::isfinite( plasma.omega( psi ) ),
		            "omega is not finite at psi = " << psi
		            << ", so C D/K has gone negative and the square root has "
		               "produced a NaN that a quadrature loop would carry silently" );
	}

	// Measured over [ -3, 3 ]: T_2 >= 0.4125, its own minimum, at psi = 1.25 --
	// the quadratic's vertex, which is inside the range; D >= 2.0154 at
	// psi = 0.594; n_s0 >= 0.0247 at psi = -0.667; C >= 1.5100 at psi = -0.7.
	// Every one of those minima is interior, so the sweep is measuring the
	// quadratics' vertices rather than their endpoints.
	std::printf( "\n  over psi in [ -3, 3 ]: min T_2 %.4f, min D %.4f, "
	             "min n_s0 %.4f, min C %.4f\n",
	             worstTemperature, worstDenominator, worstDensity, worstCoefficient );
	std::fflush( stdout );

	BOOST_TEST( worstTemperature > 0.0 );
	BOOST_TEST( worstDenominator > 0.0,
	            "Z_1 T_2 - Z_2 T_1 reaches " << worstDenominator
	            << ", so the two-species closure is singular somewhere a Newton "
	               "iterate could reach" );
	BOOST_TEST( worstDensity > 0.0 );
	BOOST_TEST( worstCoefficient > 0.0,
	            "C reaches " << worstCoefficient << ", so omega is imaginary there" );

	// T_1 turns over outside the box but inside the real line, and the margin is
	// the thing to know rather than the fact.
	BOOST_TEST( plasma.temperature( 0, -2.0 ) > 0.0 );
	BOOST_TEST( plasma.temperature( 0, -3.0 ) < 0.0,
	            "T_1 has stopped turning over at psi = -3, so the affine shape has "
	               "been changed and the margin recorded here is stale" );
}

/// Quasineutrality on R = R_ref, which every closed form in the fixture and in
/// meq::RotatingSource is derived from, holds identically -- and holds outside
/// [ 0, 1 ], which is the only range meq::RotatingSource's constructor samples.
BOOST_AUTO_TEST_CASE( theSpeciesAreNeutralOnTheReferenceCurve )
{
	Plasma const plasma = Plasma::standard();

	double worst = 0.0;
	for ( int i = -200; i <= 200; ++i )
	{
		double const psi = i/100.0;
		double sum = 0.0;
		double scale = 0.0;
		for ( std::size_t s = 0; s < Plasma::speciesCount(); ++s )
		{
			sum += plasma.charge( s )*plasma.referenceDensity( s, psi );
			scale += std::abs( plasma.charge( s )*plasma.referenceDensity( s, psi ) );
		}
		worst = std::max( worst, std::abs( sum )/std::max( 1.0e-300, scale ) );
	}

	// Measured: 0.000e+00 -- EXACTLY zero, not merely small. n_20 is -( Z_1/Z_2 )
	// times n_10 with the same shape function, so the two terms of the sum are the
	// same product with opposite signs and cancel bit for bit.
	std::printf( "  worst |Sum Z_s n_s0|/scale over psi in [ -2, 2 ]: %.3e\n", worst );
	std::fflush( stdout );

	BOOST_TEST( worst < 1.0e-15,
	            "the reference densities violate neutrality by " << worst
	            << " relative; the closures are DERIVED from Sum Z_s n_s0 = 0 and "
	               "are simply wrong without it" );

}

/// The closure's charge-weighted combinations are not plain sums here, which is
/// a coverage gap this fixture closes incidentally and which is worth recording
/// separately from the C' work.
///
/// meq::RotatingSource's two-species closed form is built on
/// Z_1 T_2 - Z_2 T_1, Z_1 m_2 - Z_2 m_1 and m_1 T_2 - m_2 T_1. AT Z = +-1 THOSE
/// COLLAPSE TO T_1 + T_2, m_1 + m_2 AND m_1 T_2 - m_2 T_1, and every other
/// two-species rotating configuration in this tree -- RotatingSourceTests,
/// RotatingSourceConvergence, RotatingNewtonConvergence,
/// RotatingNormalisedConvergence, MaschkePerrinConvergence -- runs at Z = +-1.
/// Z = +6 appears only in three-species sets, which take Closure::RootFind and
/// never form the closed-form combinations at all. So a closed form that wrote
/// the plain sums would have passed everything in the tree.
///
/// This fixture is at Z_1 = +2, and the agreement measured below in
/// theRotatingSourceReproducesTheIndependentClosedForm therefore says the
/// weighting is right as well as the chain rule. The temperature combination is
/// the one that matters, being the one that carries psi-dependence into C.
BOOST_AUTO_TEST_CASE( theChargeWeightedCombinationsAreNotPlainSums )
{
	Plasma const plasma = Plasma::standard();

	double const psi = 0.0;
	double const weightedT = plasma.closureDenominator( psi );
	double const plainT = plasma.temperature( 0, psi ) + plasma.temperature( 1, psi );
	double const weightedM = plasma.closureMassFactor();
	double const plainM = plasma.mass( 0 ) + plasma.mass( 1 );

	std::printf( "  Z_1 T_2 - Z_2 T_1 = %.6f against T_1 + T_2 = %.6f  (%.1f%% apart)\n",
	             weightedT, plainT, 100.0*std::abs( weightedT - plainT )/plainT );
	std::printf( "  Z_1 m_2 - Z_2 m_1 = %.9f against m_1 + m_2 = %.9f  (%.3e apart)\n",
	             weightedM, plainM, std::abs( weightedM - plainM )/plainM );
	std::fflush( stdout );

	// Measured: 2.100000 against 1.500000, 40.0% apart; and 4.001089 against
	// 4.000545, 1.36e-04 apart. The mass pair separates by only a part in ten
	// thousand -- m_2 is an electron and Z_1 m_2 is small either way -- but that
	// is still nine orders above the 1e-13 the comparison below is asserted at,
	// so both combinations are pinned rather than only the temperature one.
	BOOST_TEST( std::abs( weightedT - plainT )/plainT > 0.2,
	            "Z_1 T_2 - Z_2 T_1 is within " << 100.0*std::abs( weightedT - plainT )/plainT
	            << "% of T_1 + T_2, so the charges have gone back to +-1 and the "
	               "closure's weighting is untested again" );
	BOOST_TEST( std::abs( weightedM - plainM )/plainM > 1.0e-10 );
}

/*
 * ------------------------------------------------------------------------
 * The fixture's algebra IS (96) and (97)
 * ------------------------------------------------------------------------
 */

/// The bisection really solves (97), and the closed forms above it really are
/// what (96) and (97) give.
///
/// THIS IS THE CASE THAT MAKES THE COMPARISON WITH MEQ MEAN SOMETHING. Two
/// transcriptions of one derivation agreeing says the copying was consistent,
/// not that the derivation is right. So the fixture carries a second route --
/// bisection on Sum_s Z_s n_s = 0, a different algorithm from
/// meq::RotatingSource's safeguarded Newton, with C appearing nowhere in it --
/// and the closed-form p, dp/dpsi and d2p/dpsi2 are measured against it, the
/// derivatives by a fourth-order central difference of the bisected pressure.
///
/// A slip shared between MEQ's C-chain and the fixture's would pass
/// theRotatingSourceReproducesTheIndependentClosedForm and would fail here.
BOOST_AUTO_TEST_CASE( theFixtureAlgebraFollowsFromQuasineutrality )
{
	Plasma const plasma = Plasma::standard();

	// h = 1e-3 balances the difference formula's O( h^4 ) truncation against the
	// round-off it amplifies by 1/h^2: measured, the second derivative lands at
	// 4.5e-10 and neither a larger nor a smaller step improves it.
	double const step = 1.0e-3;

	double worstResidual = 0.0;
	double worstPressure = 0.0;
	double worstFirst = 0.0;
	double worstSecond = 0.0;

	overTheSweep( [ & ]( double radius, double, double psi )
	{
		double const y = plasma.potentialByBisection( radius, psi );
		worstResidual = std::max( worstResidual,
		                          std::abs( plasma.neutralityResidual( y, radius, psi ) ) );

		double const bisected = plasma.pressureByBisection( radius, psi );
		worstPressure = std::max( worstPressure,
		                          deviation( bisected, plasma.pressure( radius, psi ) ) );

		auto sampled = [ & ]( double offset )
		{
			return plasma.pressureByBisection( radius, psi + offset );
		};

		double const first = ( sampled( -2.0*step ) - 8.0*sampled( -step )
		                       + 8.0*sampled( step ) - sampled( 2.0*step ) )/( 12.0*step );
		double const second = ( -sampled( -2.0*step ) + 16.0*sampled( -step )
		                        - 30.0*sampled( 0.0 ) + 16.0*sampled( step )
		                        - sampled( 2.0*step ) )/( 12.0*step*step );

		worstFirst = std::max( worstFirst, deviation( plasma.dPressureDPsi( radius, psi ), first ) );
		worstSecond = std::max( worstSecond,
		                        deviation( plasma.d2PressureDPsi2( radius, psi ), second ) );
	} );

	std::printf( "\n  the fixture's two routes, over the sweep\n"
	             "    (97) residual at the bisected root       %10.3e\n"
	             "    p:          closed form vs bisection     %10.3e\n"
	             "    dp/dpsi:    closed form vs a difference  %10.3e\n"
	             "    d2p/dpsi2:  closed form vs a difference  %10.3e\n",
	             worstResidual, worstPressure, worstFirst, worstSecond );
	std::fflush( stdout );

	// Measured: 6.94e-17, 3.33e-16, 4.04e-12, 4.53e-10. The last two are the
	// difference formula's own floor, not a disagreement -- and both are five
	// orders below the size of the terms being checked, which is what matters.
	BOOST_TEST( worstResidual < 1.0e-14 );
	BOOST_TEST( worstPressure < 1.0e-14,
	            "the closed-form pressure differs from the quasineutral one by "
	            << worstPressure << ", so p_0 exp( C h ) is not what (96) and (97) give" );
	BOOST_TEST( worstFirst < 1.0e-9,
	            "dp/dpsi differs from a difference of the bisected pressure by "
	            << worstFirst << "; the C' term in the fixture is wrong, and MEQ "
	               "agreeing with it would then mean nothing" );
	BOOST_TEST( worstSecond < 1.0e-7,
	            "d2p/dpsi2 differs from a difference of the bisected pressure by "
	            << worstSecond << "; the C'' and C'^2 terms in the fixture are wrong" );

	// The gauge itself, since the bisection is the one route that could violate
	// it: phi_0 must vanish on R = R_ref exactly, that being what makes n_s0 the
	// physical density there.
	for ( double psi : { -0.2, 0.0, 0.5, 1.0 } )
	{
		BOOST_TEST( std::abs( plasma.potentialByBisection( Plasma::referenceRadius(), psi ) )
		            < 1.0e-15,
		            "phi_0 does not vanish on the reference curve at psi = " << psi );
	}
}

/*
 * ------------------------------------------------------------------------
 * C' is real and it is most of the answer
 * ------------------------------------------------------------------------
 */

/// C really does drift with psi, and by a lot.
///
/// WITHOUT THIS THE FILE IS SATISFIED BY MASCHKE & PERRIN'S PLASMA, where C is
/// constant and every measurement below reduces to one the tree already has. It
/// reports the drift and asserts it, and it asserts the two derivatives
/// separately -- a linear C would leave C'' identically zero and one of the
/// three drift terms of d2p/dpsi2 untested.
BOOST_AUTO_TEST_CASE( theExponentCoefficientGenuinelyVaries )
{
	Plasma const plasma = Plasma::standard();

	double lowest = 1.0e300, highest = -1.0e300;
	double smallestPrime = 1.0e300, largestPrime = -1.0e300;
	double lowestExponent = 1.0e300, highestExponent = -1.0e300;

	overTheSweep( [ & ]( double radius, double, double psi )
	{
		lowest = std::min( lowest, plasma.exponentCoefficient( psi ) );
		highest = std::max( highest, plasma.exponentCoefficient( psi ) );
		smallestPrime = std::min( smallestPrime,
		                          std::abs( plasma.exponentCoefficientPrime( psi ) ) );
		largestPrime = std::max( largestPrime,
		                         std::abs( plasma.exponentCoefficientPrime( psi ) ) );
		lowestExponent = std::min( lowestExponent, plasma.densityExponent( radius, psi ) );
		highestExponent = std::max( highestExponent, plasma.densityExponent( radius, psi ) );
	} );

	std::printf( "\n  C  in [ %8.4f, %8.4f ]   drift %.2f x\n", lowest, highest,
	             highest/lowest );
	std::printf( "  |C'| in [ %8.4f, %8.4f ]   C'' = %.4f\n", smallestPrime,
	             largestPrime, plasma.exponentCoefficientDoublePrime( 0.0 ) );
	std::printf( "  C h in [ %8.4f, %8.4f ]   so n_s varies by %.1f x across the box\n",
	             lowestExponent, highestExponent,
	             std::exp( highestExponent - lowestExponent ) );
	std::fflush( stdout );

	// Measured: C over [ 1.7600, 4.4000 ], a drift of 2.50x; |C'| over
	// [ 1.0000, 3.4000 ] and never zero on the sweep; C'' = 2.0000; C h over
	// [ -1.4080, 2.1120 ], so the densities vary by a factor of 34 across the box.
	BOOST_TEST( highest/lowest > 2.0,
	            "C moves by only " << highest/lowest << " across the sweep; that is "
	               "Maschke & Perrin's constant-C case again and this file has "
	               "stopped testing what it exists for" );
	BOOST_TEST( smallestPrime > 0.5,
	            "|C'| falls to " << smallestPrime << " somewhere on the sweep, so part "
	               "of it is measuring the constant-C case" );
	BOOST_TEST( std::abs( plasma.exponentCoefficientDoublePrime( 0.0 ) ) > 0.5,
	            "C'' is zero, so C is affine and the C'' term of d2p/dpsi2 is "
	               "untested" );
	BOOST_TEST( highestExponent - lowestExponent > 2.0,
	            "the shared exponent spans only " << highestExponent - lowestExponent
	            << ", so the rotation is a perturbation rather than sonic" );
}

/*
 * ------------------------------------------------------------------------
 * THE DELIVERABLE: two implementations of (96), (97) and (136)
 * ------------------------------------------------------------------------
 */

/// meq::RotatingSource reproduces the independent closed form, at C' != 0.
///
/// THIS IS WHAT THE FILE IS FOR. Nothing in tests/analytic/VaryingCentrifugal.hpp
/// calls into src/meq: it derives phi_0, C, p, dp/dpsi and d2p/dpsi2 from RoPP
/// (96) and (97) itself, and f and dF/dpsi from (136). The two share the six
/// profiles and nothing else.
///
/// BOTH CLOSURES ARE DRIVEN. Closure::ClosedForm is the route two species
/// normally take; Closure::RootFind is the general path, a safeguarded Newton on
/// (97) with phi_0's derivatives by implicit differentiation, and at two species
/// it must give the same answer. Checking only the first would leave the code
/// that runs for three or more species measured against nothing but itself.
BOOST_AUTO_TEST_CASE( theRotatingSourceReproducesTheIndependentClosedForm )
{
	Plasma const plasma = Plasma::standard();

	for ( auto closure : { meq::RotatingSource::Closure::ClosedForm,
	                       meq::RotatingSource::Closure::RootFind } )
	{
		std::shared_ptr<meq::RotatingSource const> const source = matching( plasma, closure );

		double worstF = 0.0;
		double worstJacobian = 0.0;
		double worstPressure = 0.0;
		double worstPotential = 0.0;
		double worstPotentialPrime = 0.0;
		double worstDensity = 0.0;

		overTheSweep( [ & ]( double radius, double z, double psi )
		{
			double const expectedF = plasma.f( radius, z, psi );
			double const actualF = source->f( radius, z, psi );
			double const expectedJ = plasma.dFdPsi( radius, z, psi );
			double const actualJ = source->dFdPsi( radius, z, psi );

			worstF = std::max( worstF, deviation( actualF, expectedF ) );
			worstJacobian = std::max( worstJacobian, deviation( actualJ, expectedJ ) );
			worstPressure = std::max( worstPressure,
			                          deviation( source->pressure( radius, psi ),
			                                     plasma.pressure( radius, psi ) ) );
			worstPotential = std::max( worstPotential,
			                           deviation( source->potential( radius, psi ),
			                                      plasma.potential( radius, psi ) ) );
			worstPotentialPrime = std::max( worstPotentialPrime,
			                                deviation( source->dPotentialDPsi( radius, psi ),
			                                           plasma.potentialPrime( radius, psi ) ) );
			for ( std::size_t s = 0; s < Plasma::speciesCount(); ++s )
			{
				worstDensity = std::max( worstDensity,
				                         deviation( source->density( s, radius, psi ),
				                                    plasma.density( s, radius, psi ) ) );
			}

			BOOST_TEST( deviation( actualF, expectedF ) < 1.0e-13,
			            "at ( " << radius << ", " << z << " ), psi = " << psi
			            << ": meq::RotatingSource gives F = " << actualF
			            << " where (136) closed by (96) and (97) gives " << expectedF );
			BOOST_TEST( deviation( actualJ, expectedJ ) < 1.0e-13,
			            "at ( " << radius << ", " << z << " ), psi = " << psi
			            << ": meq::RotatingSource gives dF/dpsi = " << actualJ
			            << " where the closed form gives " << expectedJ );
		} );

		std::printf( "\n  meq::RotatingSource against the independent closed form, "
		             "closure %s\n",
		             closure == meq::RotatingSource::Closure::ClosedForm
		                 ? "ClosedForm" : "RootFind" );
		std::printf( "    F                %10.3e     dF/dpsi   %10.3e\n",
		             worstF, worstJacobian );
		std::printf( "    p                %10.3e     e phi_0   %10.3e\n",
		             worstPressure, worstPotential );
		std::printf( "    d( e phi_0 )/dpsi %9.3e     n_s       %10.3e\n",
		             worstPotentialPrime, worstDensity );
		std::fflush( stdout );

		// Measured, ClosedForm: F 7.13e-16, dF/dpsi 5.26e-16, p 4.44e-16,
		// e phi_0 2.22e-16, its derivative 1.11e-16, n_s 4.44e-16.
		// RootFind reaches 1.07e-15 in F, the inner iteration's own floor.
		BOOST_TEST( worstF < 1.0e-13 );
		BOOST_TEST( worstJacobian < 1.0e-13 );
		BOOST_TEST( worstPressure < 1.0e-13 );
		BOOST_TEST( worstPotential < 1.0e-13 );
		BOOST_TEST( worstPotentialPrime < 1.0e-13 );
		BOOST_TEST( worstDensity < 1.0e-13 );
	}
}

/// The C' terms are most of dp/dpsi and most of d2p/dpsi2, so the agreement
/// above is not a comparison of two expressions that both reduce to the static
/// answer.
///
/// THE PATTERN IS MaschkePerrinConvergence.cpp's theCancellationIsNotVacuous, in
/// the other direction: there the terms that must CANCEL are measured, here the
/// terms that must SURVIVE are. Both are the same guard -- a test whose subject
/// is numerically negligible cannot fail for the reason it was written.
///
/// THE SHARE IS TAKEN AGAINST THE SUM OF MAGNITUDES AND NOT AGAINST F ITSELF,
/// which is not squeamishness. F is mu0 R^2 dp/dpsi + g g' and it passes through
/// zero inside the sweep -- at R = 1.4 it does so just below psi = 0, where the
/// drift term and the rest very nearly cancel. A share taken against F reads
/// 394% there, which is a true statement about a small denominator and a useless
/// one about the term. It is also why deviation() above floors at one.
BOOST_AUTO_TEST_CASE( theDriftTermsAreMostOfTheAnswer )
{
	Plasma const plasma = Plasma::standard();

	// The outboard edge, where h and therefore every drift term is largest. The
	// share is reported at three psi so that it is clear it is not one lucky
	// point.
	double const radius = standardBox().maxRadius;
	double const h = Plasma::radialFactor( radius );

	std::printf( "\n  the drift terms at R = %.2f, h = %.3f\n", radius, h );
	std::printf( "  %6s %9s %9s %7s %9s %9s %9s %9s %7s\n",
	             "psi", "p0'", "p0C'h", "share", "p0''", "2p0'C'h", "p0C''h",
	             "p0C'^2h^2", "share" );

	double smallestFirstShare = 1.0;
	double smallestSecondShare = 1.0;

	for ( double psi : { 0.0, 0.5, 1.0 } )
	{
		double const p0 = plasma.referencePressure( psi );
		double const p0Prime = plasma.referencePressurePrime( psi );
		double const p0DoublePrime = plasma.referencePressureDoublePrime( psi );
		double const cPrime = plasma.exponentCoefficientPrime( psi );
		double const cDoublePrime = plasma.exponentCoefficientDoublePrime( psi );

		double const firstDrift = p0*cPrime*h;
		double const firstShare = std::abs( firstDrift )
		                          /( std::abs( p0Prime ) + std::abs( firstDrift ) );

		double const secondA = 2.0*p0Prime*cPrime*h;
		double const secondB = p0*cDoublePrime*h;
		double const secondC = p0*cPrime*cPrime*h*h;
		double const secondDrift = std::abs( secondA ) + std::abs( secondB )
		                           + std::abs( secondC );
		double const secondShare = secondDrift
		                           /( std::abs( p0DoublePrime ) + secondDrift );

		smallestFirstShare = std::min( smallestFirstShare, firstShare );
		smallestSecondShare = std::min( smallestSecondShare, secondShare );

		std::printf( "  %6.2f %9.4f %9.4f %6.1f%% %9.4f %9.4f %9.4f %9.4f %6.1f%%\n",
		             psi, p0Prime, firstDrift, 100.0*firstShare,
		             p0DoublePrime, secondA, secondB, secondC, 100.0*secondShare );
	}
	std::fflush( stdout );

	// Measured at R = 1.4: the drift is 62.8 to 67.5% of dp/dpsi and 68.8 to
	// 88.8% of d2p/dpsi2 over psi in [ 0, 1 ]. THREE of the four terms of
	// d2p/dpsi2 carry C' or C''.
	BOOST_TEST( smallestFirstShare > 0.4,
	            "the C' term is only " << 100.0*smallestFirstShare << "% of dp/dpsi, "
	               "so the agreement measured above is largely an agreement about "
	               "the static source and the tolerance is not buying what it looks "
	               "like it is buying" );
	BOOST_TEST( smallestSecondShare > 0.5,
	            "the C' and C'' terms are only " << 100.0*smallestSecondShare
	            << "% of d2p/dpsi2" );
}

/// Break the drift term and the agreement goes away.
///
/// A TEST THAT CANNOT FAIL IS WORSE THAN NO TEST, and this repository has been
/// bitten by exactly that -- see CLAUDE.md on
/// theTwoBorderSolveReportsATrueMagneticAxis asserting a tautology. Two
/// mutations, both of which a real code could commit:
///
///   DropExponentDrift   differentiate p_0 and forget that the exponent is a
///                       flux function too. This is the term absent from every
///                       published rotating benchmark, so it is the mistake
///                       those benchmarks cannot catch.
///   FlipExponentDrift   Li & Zhu (9)'s two reversed signs. The sharper of the
///                       two, since it moves the term by twice its size -- and
///                       the more interesting, since it is a mistake somebody
///                       actually made in print.
///
/// BOTH MUTATE THE FIXTURE AND NOT src/meq, so this runs on every build rather
/// than being an experiment somebody once did with a patched library. The
/// direction of the inference is the same either way: the two implementations
/// agree, and they stop agreeing when one of them drops the term.
BOOST_AUTO_TEST_CASE( mutatingTheExponentDriftBreaksTheAgreement )
{
	Plasma const plasma = Plasma::standard();
	std::shared_ptr<meq::RotatingSource const> const source = matching( plasma );

	for ( auto how : { Plasma::Mutation::DropExponentDrift,
	                   Plasma::Mutation::FlipExponentDrift } )
	{
		Plasma const broken = plasma.mutated( how );

		double worstF = 0.0;
		double worstJacobian = 0.0;
		double worstPressure = 0.0;

		overTheSweep( [ & ]( double radius, double z, double psi )
		{
			worstF = std::max( worstF, deviation( broken.f( radius, z, psi ),
			                                      source->f( radius, z, psi ) ) );
			worstJacobian = std::max( worstJacobian,
			                          deviation( broken.dFdPsi( radius, z, psi ),
			                                     source->dFdPsi( radius, z, psi ) ) );
			worstPressure = std::max( worstPressure,
			                          deviation( broken.pressure( radius, psi ),
			                                     source->pressure( radius, psi ) ) );
		} );

		std::printf( "\n  %s: worst F deviation %10.3e, dF/dpsi %10.3e, p %10.3e\n",
		             how == Plasma::Mutation::DropExponentDrift
		                 ? "C' dropped " : "C' sign flipped",
		             worstF, worstJacobian, worstPressure );
		std::fflush( stdout );

		// Measured: dropped 7.90e-01 in F and 8.72e-01 in dF/dpsi; flipped
		// 1.58e+00 and 1.02e+00. Against the conforming fixture's 7.13e-16, which
		// is fifteen orders away.
		BOOST_TEST( worstF > 1.0e-3,
		            "the C' term was removed from F and meq::RotatingSource still "
		               "agrees to " << worstF << ", so the agreement measured above "
		               "does not rest on that term and this file is not testing it" );
		BOOST_TEST( worstJacobian > 1.0e-3,
		            "the C' and C'' terms were removed from dF/dpsi and "
		               "meq::RotatingSource still agrees to " << worstJacobian );

		// AND p IS UNMOVED, which is what says the mutation is confined to the
		// derivatives. If it moved the pressure too, the failure above would be
		// caught by the pressure check and would say nothing about the drift.
		BOOST_TEST( worstPressure < 1.0e-13,
		            "mutating C' moved the pressure by " << worstPressure
		            << ", so it is not the isolated control it claims to be" );
	}
}

/// The omega = 0 control: with C_0 = 0 the source must be mu0 R^2 p_0' + g g',
/// which is what meq::MHDSource would give on the same profiles.
///
/// It is here for the same reason MaschkePerrin.hpp's stationary() is: anything
/// the rotating case shows that this one also shows belongs to the profiles
/// rather than to the rotation. It is also the one configuration in which the
/// fixture's omegaPrime() takes its 0/0 branch -- omega is identically zero, so
/// its derivative is zero rather than 1/( 2 sqrt 0 ) -- and a NaN there would be
/// carried silently into the Jacobian.
BOOST_AUTO_TEST_CASE( theStationaryControlAgreesToo )
{
	Plasma const plasma = Plasma::stationary();
	std::shared_ptr<meq::RotatingSource const> const source = matching( plasma );

	double worstF = 0.0;
	double worstJacobian = 0.0;
	double worstStatic = 0.0;

	overTheSweep( [ & ]( double radius, double z, double psi )
	{
		worstF = std::max( worstF, deviation( source->f( radius, z, psi ),
		                                      plasma.f( radius, z, psi ) ) );
		worstJacobian = std::max( worstJacobian,
		                          deviation( source->dFdPsi( radius, z, psi ),
		                                     plasma.dFdPsi( radius, z, psi ) ) );

		// And the source really has collapsed to the static one, which is the
		// half that says C_0 = 0 reached the exponent rather than merely the
		// tolerance.
		double const staticF = Plasma::mu0()*radius*radius*plasma.referencePressurePrime( psi )
		                       + Plasma::ggPrime( psi );
		worstStatic = std::max( worstStatic, deviation( source->f( radius, z, psi ), staticF ) );
	} );

	std::printf( "\n  the C_0 = 0 control: F %10.3e, dF/dpsi %10.3e, "
	             "against mu0 R^2 p_0' + g g' %10.3e\n",
	             worstF, worstJacobian, worstStatic );
	std::fflush( stdout );

	// Measured: 0.000e+00 in all three columns. Exactly zero rather than merely
	// small, because at C_0 = 0 every exponent is exp( 0 ) = 1 to the bit and the
	// two routes are then the same sum of the same products.
	BOOST_TEST( worstF < 1.0e-14 );
	BOOST_TEST( worstJacobian < 1.0e-14 );
	BOOST_TEST( worstStatic < 1.0e-14,
	            "at C_0 = 0 the source differs from the non-rotating one by "
	            << worstStatic << ", so the exponent has not gone away" );
	BOOST_TEST( std::abs( plasma.omega( 0.5 ) ) < 1.0e-300 );
	BOOST_TEST( std::abs( plasma.omegaPrime( 0.5 ) ) < 1.0e-300,
	            "omega' is " << plasma.omegaPrime( 0.5 ) << " where omega is "
	               "identically zero, so the 0/0 branch has produced a NaN" );
}

/*
 * ------------------------------------------------------------------------
 * The manufactured solve
 * ------------------------------------------------------------------------
 */

/// The guard, and it comes before every rate: psiExact really does solve the
/// equation the solver is about to be handed.
///
/// Two independent statements, in the order a failure should be read: the
/// hand-derived Delta* agrees with a central difference of psiExact, and
/// Delta*( psiExact ) = -F_total at psiExact. The second is the construction --
/// the fixture's remainder cancelling meq::RotatingSource's f at the exact
/// solution -- and here it is a REAL cancellation between two implementations
/// rather than an identity, which is why it is worth measuring rather than
/// arguing.
BOOST_AUTO_TEST_CASE( theManufacturedSolutionSatisfiesTheEquation )
{
	ManufacturedEquilibrium const eq;
	meq::tests::Rectangle const box = standardBox();

	double worstClosedForm = 0.0;
	double worstEquation = 0.0;

	for ( double radius = box.minRadius; radius <= box.maxRadius + 1.0e-12; radius += 0.05 )
	{
		for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.075 )
		{
			double const fd = eq.deltaStarFD( radius, z );
			double const closed = eq.deltaStar( radius, z );
			double const minusF = -eq.f( radius, z, eq.psi( radius, z ) );

			worstClosedForm = std::max( worstClosedForm, std::abs( closed - fd ) );
			worstEquation = std::max( worstEquation, std::abs( fd - minusF ) );

			BOOST_TEST( std::abs( closed - fd ) < 1.0e-5,
			            "at ( " << radius << ", " << z << " ): the closed-form Delta* gives "
			            << closed << " where a central difference gives " << fd );
			BOOST_TEST( std::abs( fd - minusF ) < 1.0e-5,
			            "at ( " << radius << ", " << z << " ): Delta*( psi ) = " << fd
			            << " but -F_total = " << minusF << ", so meq::RotatingSource "
			               "and the fixture disagree by enough to stop psiExact being "
			               "a solution" );
		}
	}

	std::printf( "\n  the manufactured varying-C solution\n"
	             "    closed-form Delta* vs central difference: worst %.3e\n"
	             "    Delta*( psi ) vs -F_total:                worst %.3e\n",
	             worstClosedForm, worstEquation );
	std::fflush( stdout );
}

/*
 * THE RATES. k+1 in psi and in q over four dyadic meshes, plus an absolute
 * ceiling on the finest -- which is the only assertion here that can see a
 * solution converging beautifully to the wrong function.
 *
 * WHAT THEY ADD OVER RotatingNewtonConvergence.cpp'S IDENTICAL-LOOKING TABLES,
 * which is a fair question and the answer is narrow. Those are driven by a
 * remainder built from meq::RotatingSource's own f, so their exact solution is
 * exact by construction; these are driven by the fixture's, so they degrade if
 * the two implementations part company. That is a blunt instrument beside the
 * pointwise sweep above -- it could not see a disagreement smaller than the
 * discretisation error -- and it is not what closes the C' gap. What it does
 * reach that the sweep cannot is the ASSEMBLY: meq::SourceIntegrator evaluating
 * this source and its Jacobian at quadrature points, on a species set with
 * Z_1 = +2 and a quadratic C.
 *
 * The ceilings sit at roughly three times the measured finest-mesh error.
 */

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Measured at h = 0.025: psi 4.292075e-04, q 7.473691e-04, rates 1.999 and
	// 1.993.
	meq::tests::checkOrder( ManufacturedEquilibrium(),
	                        "manufactured varying-C rotating, two species", 1,
	                        1.3e-3, 2.3e-3 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured at h = 0.025: psi 3.575050e-06, q 6.625524e-06, rates 2.997 and
	// 2.998.
	meq::tests::checkOrder( ManufacturedEquilibrium(),
	                        "manufactured varying-C rotating, two species", 2,
	                        1.1e-5, 2.0e-5 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured at h = 0.025: psi 2.239704e-08, q 4.061468e-08, rates 3.999 and
	// 3.998.
	meq::tests::checkOrder( ManufacturedEquilibrium(),
	                        "manufactured varying-C rotating, two species", 3,
	                        6.8e-8, 1.3e-7 );
}

/// Newton's observed order, which is the assertion with teeth on dF/dpsi once
/// the source is assembled.
///
/// CLAUDE_FLOW.md's calibration is the reason this case exists beside the rate
/// tables: a dF/dpsi 5% too large leaves every error and every rate unchanged to
/// all seven digits printed and drops the observed order from 1.980 to 1.055.
/// The whole psi-dependence of F_total is meq::RotatingSource's, the remainder
/// being a function of ( R, z ) alone, so there is nothing else in the Jacobian
/// to be right.
///
/// Per the testing stance the assertion is on the BEST triple above the
/// round-off floor -- Newton is not monotone on a source like this -- and it is
/// bounded BOTH sides, because a best-of reading much above 2 is a triple
/// straddling the floor rather than a better Jacobian.
BOOST_AUTO_TEST_CASE( newtonConvergesAtOrderTwo )
{
	std::vector<double> history;
	Measurement const point = meq::tests::measure( ManufacturedEquilibrium(),
	                                              standardBox(), 3, 8, &history );

	BOOST_TEST_REQUIRE( history.size() >= 3u,
	                    "Newton produced only " << history.size() << " residuals, so "
	                    "there is nothing to estimate an order from" );

	std::printf( "\n  Newton residual history, varying-C manufactured case, "
	             "k = 3, h = %.5f, %d trace dofs\n", point.h, point.traceDofs );
	std::printf( "  %5s %16s %16s %8s\n", "it", "||r||", "||r||/||r_0||", "order" );
	for ( std::size_t i = 0; i < history.size(); ++i )
	{
		if ( i >= 2 )
		{
			std::printf( "  %5zu %16.6e %16.6e %8.3f\n", i, history[ i ],
			             history[ i ]/history[ 0 ],
			             newtonOrder( history[ i - 2 ], history[ i - 1 ], history[ i ] ) );
		}
		else
		{
			std::printf( "  %5zu %16.6e %16.6e %8s\n", i, history[ i ],
			             history[ i ]/history[ 0 ], "-" );
		}
	}
	std::fflush( stdout );

	double const best = bestNewtonOrder( history );

	// THE THRESHOLD IS CALIBRATED AND NOT COPIED, and the case below is the
	// calibration: a consistent Jacobian reads 1.792 here and a 5% inconsistent
	// one reads 1.173, so 1.5 sits between them with room on both sides.
	//
	// 1.792 RATHER THAN ~2 IS THE HISTORY AND NOT THE JACOBIAN, which is worth
	// saying because the neighbouring file reads 1.980 on the same box, the same
	// degree and the same mesh. The whole run is 5.73e-01 -> 1.07e-03 -> 1.38e-08
	// -> 7.72e-16: the FIRST step already takes five hundred off the residual, so
	// the triple ( R_0, R_1, r2 ) straddles the pre-asymptotic step and the
	// three-point estimate log( r2/R_1 )/log( R_1/R_0 ) reads low; and the triple
	// ( R_1, r2, r3 ) is dropped by bestNewtonOrder's floor because r3 is at
	// 1.3e-15 of R_0, which is round-off. Checked directly, r2/R_1^2 = 1.2e-2 and
	// r3/r2^2 = 4.0, both O( 1 ) -- the iteration IS quadratic and there is
	// simply no clean tail to read an order off. That is the trap CLAUDE_FLOW.md
	// records as "the best observed Newton order is not the order", met from the
	// other side.
	BOOST_TEST( best >= 1.5,
	            "the best observed Newton order is " << best << ", below the 1.5 that "
	               "separates a consistent Jacobian from a 5% inconsistent one on "
	               "this fixture. Every psi-dependent term here comes from "
	               "meq::RotatingSource, and the C' terms are most of it" );
	BOOST_TEST( best <= 2.5,
	            "the best observed Newton order is " << best << ", above 2; that is a "
	               "triple straddling the round-off floor rather than a better "
	               "Jacobian" );

	// The problem must actually be non-linear here, or the order above measures
	// nothing. One step is what a Solov'ev, a Li & Zhu rotating or a Maschke &
	// Perrin source gives, and is exactly why none of them can stand in for this.
	BOOST_TEST( point.newtonIterations >= 2,
	            "Newton converged in " << point.newtonIterations << " iterations, so "
	               "the source is behaving linearly in psi and dFdPsi is untested" );
	BOOST_TEST( point.newtonIterations <= 8,
	            "Newton took " << point.newtonIterations << " iterations to reach the "
	               "tolerance" );
}

/// THE CALIBRATION, AND WITHOUT IT THE CASE ABOVE IS A NUMBER WITH NO SCALE ON
/// IT.
///
/// CLAUDE_FLOW.md records the +5% experiment on RotatingSource::dFdPsi -- every
/// L2 error and every rate unchanged to all seven digits, Newton 3 iterations to
/// 6, observed order 1.980 to 1.055 -- and this reproduces it HERE, on this
/// fixture, rather than borrowing the calibration from another one. The
/// perturbation is applied to the Jacobian alone, so the residual and therefore
/// the converged answer are untouched: that is what makes the L2 half of the
/// measurement meaningful.
///
/// It is a control and not a defect report. What it establishes is that
/// newtonConvergesAtOrderTwo's threshold sits between a consistent Jacobian and
/// a 5% inconsistent one, which is the only thing that makes the threshold worth
/// asserting.
BOOST_AUTO_TEST_CASE( aPerturbedJacobianCostsNewtonItsOrderAndNotItsAnswer )
{
	std::vector<double> exactHistory;
	Measurement const exact = meq::tests::measure( ManufacturedEquilibrium(),
	                                              standardBox(), 3, 8, &exactHistory );

	std::vector<double> perturbedHistory;
	Measurement const perturbed
		= meq::tests::measure( ManufacturedEquilibrium( Plasma::standard(), 1.05 ),
		                       standardBox(), 3, 8, &perturbedHistory );

	std::printf( "\n  the +5%% Jacobian experiment, k = 3, h = %.5f\n", exact.h );
	std::printf( "  %14s %14s %14s %10s %8s\n",
	             "dF/dpsi", "L2(psi)", "L2(q)", "order", "Newton" );
	std::printf( "  %14s %14.6e %14.6e %10.3f %8d\n", "exact",
	             exact.errorPsi, exact.errorFlux, bestNewtonOrder( exactHistory ),
	             exact.newtonIterations );
	std::printf( "  %14s %14.6e %14.6e %10.3f %8d\n", "x 1.05",
	             perturbed.errorPsi, perturbed.errorFlux,
	             bestNewtonOrder( perturbedHistory ), perturbed.newtonIterations );
	std::fflush( stdout );

	// Measured: order 1.792 -> 1.173 and Newton 3 -> 4 iterations, while
	// L2( psi ) and L2( q ) do not move in the seventh significant figure. The
	// gap is what newtonConvergesAtOrderTwo's 1.5 threshold sits in.
	BOOST_TEST( bestNewtonOrder( perturbedHistory ) < 1.5,
	            "a 5% error in dF/dpsi still gives an observed Newton order of "
	            << bestNewtonOrder( perturbedHistory )
	            << ", above the threshold the case above asserts -- so that "
	               "threshold no longer separates a consistent Jacobian from an "
	               "inconsistent one and it is not testing anything" );
	BOOST_TEST( bestNewtonOrder( exactHistory ) - bestNewtonOrder( perturbedHistory )
	            > 0.3,
	            "the perturbation moves the observed order by only "
	            << bestNewtonOrder( exactHistory ) - bestNewtonOrder( perturbedHistory ) );
	BOOST_TEST( perturbed.newtonIterations > exact.newtonIterations,
	            "a 5% error in dF/dpsi costs no extra Newton iterations" );

	// AND THE ANSWER DOES NOT MOVE, which is the half that says a rate table
	// cannot see this. Both errors agree to every digit printed.
	BOOST_TEST( std::abs( perturbed.errorPsi - exact.errorPsi ) < 1.0e-6*exact.errorPsi,
	            "the perturbed Jacobian changed L2( psi ) from " << exact.errorPsi
	            << " to " << perturbed.errorPsi << "; it should change the WORK and "
	               "not the ANSWER, and if it changes the answer the residual has "
	               "been perturbed too and the experiment is not the one described" );
	BOOST_TEST( std::abs( perturbed.errorFlux - exact.errorFlux ) < 1.0e-6*exact.errorFlux );
}
