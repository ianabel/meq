// Unit tests for the coil model in src/meq/Coils.
//
// FREE-BOUNDARY-PLAN.md stage FB-2 puts a prescribed current inside the domain
// and checks the exterior. meq::CoilSet is the current, and its acceptance --
// "ComputeOutwardFlux against the total current, which is the sharpest
// whole-assembly test available" -- rests on two things being right that no
// convergence rate could see:
//
//   * THE SOURCE FACTOR. F = mu0 r j_phi, derived in Coils.cpp from
//     Delta* psi = -mu0 r j_phi rather than transcribed. An extra or missing r
//     converges at the full rate to the wrong function, and CLAUDE.md records
//     that failure happening in this codebase already.
//   * THE SIGN of the outward flux. FREE-BOUNDARY-PLAN.md section 7 predicts
//     plainly that the coupling sign will be got wrong at least once. This file
//     pins it here, where the answer is mu0 I exactly and there is nothing else
//     in the way.
//
// MFEM-FREE, like the unit under test, so CI can run it -- which for a stage
// whose acceptance is these two numbers means CI can gate FB-2's foundation
// even though it cannot build the solver.

#define BOOST_TEST_MODULE CoilsTests
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <utility>

#include "meq/Coils.hpp"

// The fixture this file's CurrentFilament is checked against. It is MFEM-free,
// and tests/CMakeLists.txt already puts tests/ on the include path so that a
// unit test can reach the closed forms in analytic/ -- SourceFactoryTests does
// the same thing for the same reason.
#include "analytic/CurrentLoop.hpp"

using meq::Coil;
using meq::CoilSet;
using meq::CurrentFilament;

namespace
{
	/// One coil, well away from the axis, with a current large enough that the
	/// numbers are readable and small enough that nothing overflows.
	Coil standardCoil()
	{
		return Coil( 2.0, 0.0, 0.15, 0.20, 1.0e6 );
	}

	/// Delta* psi = d_rr psi - ( 1/r ) d_r psi + d_zz psi, by central
	/// differences, in EXACTLY the arrangement tests/analytic/VacuumHarmonic.hpp
	/// and ManufacturedNonlinear.hpp use. Copied rather than reinvented so that
	/// a disagreement between this file and those is a real one.
	template <typename Psi>
	double deltaStarFD( Psi psi, double r, double z, double h )
	{
		auto innerR = [ & ]( double rr )
		{
			return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
		};
		double const dRInner = ( innerR( r + h ) - innerR( r - h ) )/( 2.0*h );
		double const dZZ = ( psi( r, z + h ) - 2.0*psi( r, z ) + psi( r, z - h ) )
		                   /( h*h );
		return r*dRInner + dZZ;
	}
}

/*
 * THE SOURCE FACTOR, AND A CONTROL THAT SAYS AN r WOULD BE CAUGHT.
 *
 * f() and psi() are two independent statements of the same physics: one is the
 * source term the solver is fed, the other is the field that source produces.
 * Checking them against each other by recomputing Delta* is the only check here
 * that can see a missing or extra factor of r -- a convergence study cannot,
 * because a wrong constant converges at the full rate to the wrong function.
 *
 * The interior points are what pin the constant. OUTSIDE the coil Delta* psi is
 * zero and stays zero under ANY rescaling of psi, so an exterior check
 * establishes that the field is harmonic and nothing at all about its size.
 * That is worth stating because the exterior residual is the smaller and
 * prettier number, and it is the one that means less.
 */
BOOST_AUTO_TEST_CASE( the_source_is_mu0_r_jphi_and_an_extra_r_would_be_caught )
{
	CoilSet set;
	set.add( standardCoil() );

	auto psi = [ &set ]( double r, double z ) { return set.psi( r, z ); };

	double const centreR = 2.0;
	double const centreZ = 0.0;
	double const expected = set.f( centreR, centreZ );

	std::printf( "\n  the source at the coil centre: F = %.6e\n", expected );

	// Delta* psi must be -F, inside the coil.
	double const h = 5.0e-3;
	double const measured = deltaStarFD( psi, centreR, centreZ, h );

	std::printf( "    Delta*_FD psi = %.6e against -F = %.6e, relative %.3e\n",
	             measured, -expected,
	             std::fabs( measured + expected )/std::fabs( expected ) );

	BOOST_TEST( std::fabs( measured + expected )/std::fabs( expected ) < 1.0e-3,
	            "Delta* of the coil's own field is not -f(). One of the two is "
	            "wrong, and since they are independent statements of the same "
	            "physics that is a real defect rather than a tolerance" );

	/*
	 * THE CONTROL. The plausible wrong answers are F/r and F*r, and at r = 2 they
	 * differ from F by a factor of two -- far outside any finite-difference
	 * floor. Asserting that they are REJECTED is what turns the check above from
	 * "these two agree" into "these two agree and the neighbouring conventions
	 * do not".
	 */
	double const wrongDivided = expected/centreR;
	double const wrongMultiplied = expected*centreR;

	std::printf( "    the neighbouring conventions: F/r = %.6e, F*r = %.6e\n",
	             wrongDivided, wrongMultiplied );

	BOOST_TEST( std::fabs( measured + wrongDivided )/std::fabs( wrongDivided ) > 0.1,
	            "Delta* psi is as close to -F/r as it is to -F, so this test "
	            "cannot tell the two conventions apart and is not doing its job" );
	BOOST_TEST( std::fabs( measured + wrongMultiplied )/std::fabs( wrongMultiplied ) > 0.1,
	            "Delta* psi is as close to -F*r as it is to -F" );
}

/*
 * The source is uniform inside a coil and exactly zero outside it.
 *
 * A top-hat, and the header says what that costs: a discontinuous right-hand
 * side caps the achievable convergence rate unless the mesh resolves the coil
 * edges. That is precisely why FB-2's acceptance is the outward flux rather
 * than a rate -- and it is the same defect tests/analytic/ExteriorMatched.hpp
 * exists to avoid, met from the other side.
 */
BOOST_AUTO_TEST_CASE( the_source_is_a_top_hat )
{
	Coil const c = standardCoil();
	CoilSet set;
	set.add( c );

	double const inside = set.f( c.centreR(), c.centreZ() );
	BOOST_TEST( inside != 0.0 );

	// Uniform within: f depends on r through mu0 * r * j, so it is NOT constant
	// in r -- it is the CURRENT DENSITY that is uniform. Checking the ratio is
	// what states that correctly.
	double const atEdge = set.f( c.rMax() - 1.0e-9, c.centreZ() );
	double const ratio = atEdge/inside;
	double const expectedRatio = ( c.rMax() - 1.0e-9 )/c.centreR();

	std::printf( "\n  f scales as r within the coil: %.9f against %.9f\n",
	             ratio, expectedRatio );
	BOOST_TEST( std::fabs( ratio - expectedRatio ) < 1.0e-8,
	            "f is not proportional to r inside the coil, so it is not "
	            "mu0 r j_phi with j uniform" );

	// And exactly zero outside, at every side.
	for ( double r : { c.rMin() - 1.0e-6, c.rMax() + 1.0e-6 } )
		BOOST_TEST( set.f( r, c.centreZ() ) == 0.0 );
	for ( double z : { c.zMin() - 1.0e-6, c.zMax() + 1.0e-6 } )
		BOOST_TEST( set.f( c.centreR(), z ) == 0.0 );

	BOOST_TEST( set.indexContaining( c.centreR(), c.centreZ() ) == 0 );
	BOOST_TEST( set.indexContaining( c.rMax() + 1.0, c.centreZ() ) == -1 );
}

/*
 * FB-2's OWN ACCEPTANCE, AT UNIT SCALE: the outward flux is -mu0 I.
 *
 * Integrating the Grad-Shafranov equation over a region enclosing the coils and
 * applying the divergence theorem gives
 *
 *     oint ( 1/r ) dpsi/dn dl  =  -mu0 * totalCurrent()
 *
 * with no discretisation anywhere in it. FREE-BOUNDARY-PLAN.md calls the same
 * identity "the sharpest whole-assembly test available" for FB-2; here it is
 * checked on the exact field, so that when FB-2 checks it on a SOLVE any
 * discrepancy is the solve rather than the identity.
 *
 * THE SIGN IS NEGATIVE AND IS ASSERTED AS SUCH. Section 7 says the coupling
 * sign will be got wrong at least once; this is the cheapest place to pin it.
 * The same identity written as a counterclockwise circulation of B in ( r, z )
 * comes out POSITIVE, because phi-hat = z-hat x r-hat -- so both signs are
 * defensible statements about different quantities, which is exactly how a sign
 * error survives review.
 */
BOOST_AUTO_TEST_CASE( the_outward_flux_is_minus_mu0_times_the_total_current )
{
	CoilSet set;
	set.add( Coil( 1.6, -0.3, 0.10, 0.12, 6.0e5 ) );
	set.add( Coil( 2.2, 0.4, 0.12, 0.10, -2.0e5 ) );

	BOOST_TEST( set.totalCurrent() == 4.0e5 );

	// A rectangle enclosing both coils with room to spare.
	double const rA = 1.0, rB = 3.0, zA = -1.2, zB = 1.2;

	// The flux of ( 1/r ) grad-bar psi outward through it, by a midpoint rule
	// refined and Richardson-extrapolated -- the header records the same
	// construction giving 1.07e-08.
	auto flux = [ &set ]( int m, double rA_, double rB_, double zA_, double zB_ )
	{
		double const dr = ( rB_ - rA_ )/m;
		double const dz = ( zB_ - zA_ )/m;
		double const eps = 1.0e-6;
		double total = 0.0;

		// d psi/dn by a central difference, divided by r, times the length.
		auto normalDerivative = [ &set, eps ]( double r, double z,
		                                       double nr, double nz )
		{
			double const plus = set.psi( r + eps*nr, z + eps*nz );
			double const minus = set.psi( r - eps*nr, z - eps*nz );
			return ( plus - minus )/( 2.0*eps )/r;
		};

		for ( int i = 0; i < m; ++i )
		{
			double const r = rA_ + ( i + 0.5 )*dr;
			total += normalDerivative( r, zB_, 0.0, 1.0 )*dr;   // top,    n = +z
			total += normalDerivative( r, zA_, 0.0, -1.0 )*dr;  // bottom, n = -z
			double const z = zA_ + ( i + 0.5 )*dz;
			total += normalDerivative( rB_, z, 1.0, 0.0 )*dz;   // right,  n = +r
			total += normalDerivative( rA_, z, -1.0, 0.0 )*dz;  // left,   n = -r
		}
		return total;
	};

	double const coarse = flux( 200, rA, rB, zA, zB );
	double const fine = flux( 400, rA, rB, zA, zB );
	// Midpoint is second order, so Richardson is ( 4 fine - coarse )/3.
	double const extrapolated = ( 4.0*fine - coarse )/3.0;
	double const expected = -set.mu0()*set.totalCurrent();

	std::printf( "\n  the outward flux against -mu0 I\n" );
	std::printf( "    m = 200      : %.10e\n", coarse );
	std::printf( "    m = 400      : %.10e\n", fine );
	std::printf( "    Richardson   : %.10e\n", extrapolated );
	std::printf( "    -mu0 I       : %.10e\n", expected );
	std::printf( "    relative     : %.3e\n",
	             std::fabs( extrapolated - expected )/std::fabs( expected ) );
	std::fflush( stdout );

	BOOST_TEST( std::fabs( extrapolated - expected )/std::fabs( expected ) < 1.0e-6,
	            "the outward flux of ( 1/r ) grad psi is not -mu0 I. This is an "
	            "identity with no discretisation in it, so a failure is the coil "
	            "field or the source constant, not a resolution question -- and "
	            "it is what FB-2 will check a SOLVE against" );

	// The SIGN, separately and loudly.
	BOOST_TEST( extrapolated < 0.0,
	            "the outward flux came out positive for a positive total "
	            "current. The identity is oint ( 1/r ) dpsi/dn = -mu0 I; the "
	            "POSITIVE convention belongs to the circulation of B, which is a "
	            "different quantity. Both are defensible sentences, which is how "
	            "a sign error survives review" );
}

/*
 * A thin coil is a filament, and it approaches one at second order.
 *
 * This is what says the finite cross-section is a smoothing of the filament
 * rather than a different object -- which matters because FB-1's reference
 * field, tests/analytic/CurrentLoop.hpp, IS the filament, and because the plan
 * originally proposed filament loop fields as FB-1's exact answer.
 */
BOOST_AUTO_TEST_CASE( a_thin_coil_approaches_the_filament_at_second_order )
{
	double const current = 1.0e6;
	double const centreR = 1.5;
	double const fieldR = 2.4;
	double const fieldZ = 0.7;

	std::printf( "\n  a shrinking coil against the filament it becomes\n" );
	std::printf( "    %10s %18s %12s\n", "half-width", "psi", "ratio" );

	double previous = 0.0;
	double previousDifference = 0.0;
	bool first = true;

	// The filament limit, from a very thin coil rather than from a second
	// implementation: the point is the LIMIT, and taking it from the same code
	// keeps this a statement about the cross-section rather than about two
	// transcriptions of an elliptic integral.
	CoilSet reference;
	reference.add( Coil( centreR, 0.0, 1.0e-7, 1.0e-7, current ) );
	double const limit = reference.psi( fieldR, fieldZ );

	for ( double halfSize : { 0.16, 0.08, 0.04, 0.02, 0.01 } )
	{
		CoilSet set;
		set.add( Coil( centreR, 0.0, halfSize, halfSize, current ) );
		double const value = set.psi( fieldR, fieldZ );
		double const difference = std::fabs( value - limit );

		double ratio = 0.0;
		if ( !first && difference > 0.0 )
			ratio = previousDifference/difference;

		std::printf( "    %10.4f %18.10e %12.3f\n", halfSize, value, ratio );

		if ( !first && previousDifference > 0.0 && difference > 0.0 )
			BOOST_TEST( ratio > 3.0,
			            "halving the cross-section reduced the difference from "
			            "the filament by only " << ratio << "; second order "
			            "demands about 4" );

		previous = value;
		previousDifference = difference;
		first = false;
	}
	(void)previous;
	std::fflush( stdout );
}

/*
 * The refusals, which are the contract.
 *
 * A coil reaching the axis is the one worth singling out: the operator's 1/r is
 * not integrable through r = 0, so a coil straddling it is not merely unusual,
 * it is unsolvable -- the same refusal meq::BoundaryShape makes, for the same
 * reason, and stage FB-A measured what the axis costs even when it is only
 * touched by the mesh.
 */
BOOST_AUTO_TEST_CASE( the_refusals_are_the_contract )
{
	// Reaching or crossing the axis.
	BOOST_CHECK_THROW( Coil( 0.10, 0.0, 0.10, 0.1, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( Coil( 0.05, 0.0, 0.10, 0.1, 1.0 ), std::invalid_argument );

	// Degenerate cross-sections.
	BOOST_CHECK_THROW( Coil( 2.0, 0.0, 0.0, 0.1, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( Coil( 2.0, 0.0, 0.1, 0.0, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( Coil( 2.0, 0.0, -0.1, 0.1, 1.0 ), std::invalid_argument );

	// Zero current is legitimate: a coil that is present and carrying nothing
	// is an ordinary configuration, and refusing it would make a scan through
	// zero impossible.
	BOOST_CHECK_NO_THROW( Coil( 2.0, 0.0, 0.1, 0.1, 0.0 ) );

	CoilSet set;
	// An empty set is legitimate too, and its total current is zero.
	BOOST_TEST( set.size() == 0u );
	BOOST_TEST( set.totalCurrent() == 0.0 );
	BOOST_TEST( set.f( 2.0, 0.0 ) == 0.0 );

	BOOST_CHECK_THROW( set.coil( 0 ), std::out_of_range );

	set.add( standardCoil() );
	BOOST_CHECK_THROW( set.psiOf( 5, 2.0, 0.0 ), std::out_of_range );
	BOOST_CHECK_THROW( set.setQuadratureOrder( 0 ), std::invalid_argument );
}

/*
 * THE ADAPTERS, AND THE ONE THING ABOUT THEM THAT CAN BE SILENTLY WRONG.
 *
 * meq::CoilAugmentedSource and its normalised sibling exist so that a coil set
 * can be handed to the solver as a meq::Source, which is what the driver does.
 * The sum and the derivative pass-through are easy to get right and are checked
 * below for completeness. What is NOT easy is the plasma support: it is
 * consulted by whichever object evaluates the profiles, and that is the WRAPPED
 * source -- so a setPlasmaSupport() that set only the wrapper's flag would
 * compile, run, converge, and leave the plasma unconfined. The three cases at
 * the end are about that, and the second of them goes through a base reference
 * on purpose, because that is the call a non-virtual method would get wrong.
 */
namespace
{
	/// A plasma source that is a plain function of psi, so the sum below has a
	/// closed form and the coil half is what is being isolated.
	class LinearPlasma : public meq::NormalisedSource
	{
		public:
			double f( double r, double, double psi ) const override
			{
				if ( !insidePlasma( psi ) )
					return 0.0;
				return r*( psi - boundaryValue )/( axisValue - boundaryValue );
			}

			double dFdPsi( double r, double, double psi ) const override
			{
				if ( !insidePlasma( psi ) )
					return 0.0;
				return r/( axisValue - boundaryValue );
			}

			void setNormalisation( double psiAxis, double psiBoundary ) override
			{
				axisValue = psiAxis;
				boundaryValue = psiBoundary;
			}

			double normalisation() const override { return axisValue; }
			double boundaryNormalisation() const override { return boundaryValue; }

		private:
			double axisValue = 1.0;
			double boundaryValue = 0.0;
	};
}

BOOST_AUTO_TEST_CASE( the_augmented_source_is_the_sum_and_the_coils_are_not_in_the_jacobian )
{
	auto coils = std::make_shared<CoilSet>();
	coils->add( standardCoil() );

	auto plasma = std::make_shared<LinearPlasma>();
	plasma->setNormalisation( 2.0, 0.0 );

	meq::CoilAugmentedNormalisedSource sum( plasma, coils );

	// Inside the coil, where both terms are non-zero, and outside it, where
	// only the plasma is. The coil's own contribution is checked against
	// CoilSet::f rather than recomputed, since Coils.hpp's derivation is what
	// the earlier cases in this file pin.
	for ( auto point : { std::make_pair( 2.0, 0.0 ), std::make_pair( 1.0, 0.5 ) } )
	{
		double const r = point.first, z = point.second;
		for ( double psi : { -0.7, 0.0, 0.4, 1.9 } )
		{
			BOOST_TEST( sum.f( r, z, psi )
			            == plasma->f( r, z, psi ) + coils->f( r, z ),
			            boost::test_tools::tolerance( 1.0e-15 ) );

			// EXACT EQUALITY, not a tolerance: the coils contribute nothing at
			// all to the Jacobian, so this is the same double travelling
			// through one more function call. A tolerance here would accept a
			// coil term that had leaked into dF/dpsi and happened to be small.
			BOOST_TEST( sum.dFdPsi( r, z, psi ) == plasma->dFdPsi( r, z, psi ) );
		}
	}

	// The coil term really is doing something at the first point, or the check
	// above is satisfied by two zeros.
	BOOST_TEST( coils->f( 2.0, 0.0 ) != 0.0 );
	BOOST_TEST( coils->f( 1.0, 0.5 ) == 0.0 );

	// The plain adapter agrees with the normalised one on the same pair, which
	// is what says neither of them applies the normalisation twice.
	meq::CoilAugmentedSource plain( plasma, coils );
	BOOST_TEST( plain.f( 2.0, 0.0, 0.4 ) == sum.f( 2.0, 0.0, 0.4 ) );
}

BOOST_AUTO_TEST_CASE( the_normalisation_is_the_wrapped_sources_and_there_is_only_one_of_it )
{
	auto coils = std::make_shared<CoilSet>();
	coils->add( standardCoil() );
	auto plasma = std::make_shared<LinearPlasma>();

	meq::CoilAugmentedNormalisedSource sum( plasma, coils );

	// Set through the WRAPPER and read back through the plasma, then the other
	// way round. Two stored copies of the pair would pass one direction and
	// fail the other.
	sum.setNormalisation( 3.0, -0.5 );
	BOOST_TEST( plasma->normalisation() == 3.0 );
	BOOST_TEST( plasma->boundaryNormalisation() == -0.5 );
	BOOST_TEST( sum.normalisation() == 3.0 );
	BOOST_TEST( sum.boundaryNormalisation() == -0.5 );

	// The one-argument convenience of the base class must still be reachable
	// through the override; a plain `override` without the using-declaration
	// hides it and the fixed-boundary call stops compiling.
	sum.setNormalisation( 4.0 );
	BOOST_TEST( plasma->normalisation() == 4.0 );
	BOOST_TEST( plasma->boundaryNormalisation() == 0.0 );

	plasma->setNormalisation( 5.0, 1.0 );
	BOOST_TEST( sum.normalisation() == 5.0 );
	BOOST_TEST( sum.boundaryNormalisation() == 1.0 );
}

BOOST_AUTO_TEST_CASE( the_plasma_support_reaches_the_source_that_evaluates_the_profiles )
{
	auto coils = std::make_shared<CoilSet>();
	coils->add( standardCoil() );
	auto plasma = std::make_shared<LinearPlasma>();
	plasma->setNormalisation( 2.0, 0.0 );

	auto sum = std::make_shared<meq::CoilAugmentedNormalisedSource>( plasma, coils );

	// Off by default, and then on -- through a BASE REFERENCE, which is the
	// handle the driver holds and the call a non-virtual setPlasmaSupport()
	// would send to the wrong object.
	BOOST_TEST( sum->plasmaSupport() == false );
	meq::NormalisedSource &asBase = *sum;
	asBase.setPlasmaSupport( true );

	BOOST_TEST( sum->plasmaSupport() == true );
	BOOST_TEST( plasma->plasmaSupport() == true,
	            "setPlasmaSupport() did not reach the wrapped source, so the "
	            "moving plasma support would silently do nothing" );

	// AND THE COILS ARE OUTSIDE THE SUPPORT, WHICH IS THE POINT OF THE ORDER
	// THE SUM IS TAKEN IN. At psi below the boundary value the plasma term is
	// switched off and the coil term is not -- a coil sits in the vacuum
	// region by construction, so confining it to the plasma would switch off
	// every coil in the machine.
	double const outside = -1.0;
	BOOST_TEST( plasma->f( 2.0, 0.0, outside ) == 0.0 );
	BOOST_TEST( sum->f( 2.0, 0.0, outside ) == coils->f( 2.0, 0.0 ) );
	BOOST_TEST( sum->f( 2.0, 0.0, outside ) != 0.0 );

	// Inside, both terms are live.
	BOOST_TEST( sum->f( 2.0, 0.0, 1.0 )
	            == plasma->f( 2.0, 0.0, 1.0 ) + coils->f( 2.0, 0.0 ),
	            boost::test_tools::tolerance( 1.0e-15 ) );

	asBase.setPlasmaSupport( false );
	BOOST_TEST( plasma->plasmaSupport() == false );
	BOOST_TEST( sum->f( 2.0, 0.0, outside ) != coils->f( 2.0, 0.0 ) );
}

BOOST_AUTO_TEST_CASE( the_adapters_refuse_a_null_half )
{
	auto coils = std::make_shared<CoilSet>();
	auto plasma = std::make_shared<LinearPlasma>();

	BOOST_CHECK_THROW( meq::CoilAugmentedSource( nullptr, coils ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( meq::CoilAugmentedSource( plasma, nullptr ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( meq::CoilAugmentedNormalisedSource( nullptr, coils ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( meq::CoilAugmentedNormalisedSource( plasma, nullptr ),
	                   std::invalid_argument );
}

/*
 * THE GRADIENT, AGAINST A RICHARDSON-EXTRAPOLATED DIFFERENCE OF psi.
 *
 * A plain central difference cannot settle this and it is worth saying why
 * rather than discovering it: it carries its own O( h^2 ) truncation, so the
 * comparison FLOORS at the instrument rather than at the derivative, and no
 * choice of h gets past it -- CLAUDE.md records that in four separate places,
 * from Zernike's derivative to the flux-surface identity, and each time the cure
 * is ( 4 D( h/2 ) - D( h ) )/3.
 *
 * Both columns are printed, because the PLAIN one is the control: if the two
 * agreed there would be nothing to demonstrate, and its being three or four
 * orders worse is what says the extrapolation is doing the work.
 */
BOOST_AUTO_TEST_CASE( the_gradient_is_analytic_and_a_plain_difference_cannot_see_it )
{
	auto richardson = [ ]( auto value, double h )
	{
		double const coarse = ( value( h ) - value( -h ) )/( 2.0*h );
		double const fine = ( value( 0.5*h ) - value( -0.5*h ) )/h;
		return std::make_pair( coarse, ( 4.0*fine - coarse )/3.0 );
	};

	std::printf( "\n  the analytic gradient against a difference of psi\n" );
	std::printf( "    %-28s %12s %12s\n", "case", "plain", "Richardson" );

	double worstPlain = 0.0;
	double worstRich = 0.0;

	// A FILAMENT and a COIL, because they are different code paths: the
	// filament evaluates the derivative kernel once and the coil integrates it
	// over a cross-section, so a sign or a factor could be right in one and
	// wrong in the other.
	CurrentFilament const filament( 1.5, 0.1, 1.0e6 );
	Coil const coil = standardCoil();

	struct Probe { char const *name; double r; double z; bool isCoil; };
	std::vector<Probe> const probes = {
		{ "filament, outboard",       2.40, 0.70, false },
		{ "filament, above",          1.60, 1.30, false },
		{ "filament, inboard",        0.60, 0.20, false },
		{ "coil, 5 cm clear",         2.20, 0.00, true  },
		{ "coil, far field",          5.00, 3.00, true  },
		{ "coil, inboard",            1.00, 0.40, true  },
	};

	for ( Probe const &probe : probes )
	{
		double dR = 0.0;
		double dZ = 0.0;
		auto psiAt = [ & ]( double r, double z )
		{
			return probe.isCoil ? meq::coilPsi( coil, r, z )
			                    : meq::filamentPsi( filament, r, z );
		};
		if ( probe.isCoil )
			meq::coilGradPsi( coil, probe.r, probe.z, dR, dZ );
		else
			meq::filamentGradPsi( filament, probe.r, probe.z, dR, dZ );

		double const h = 1.0e-3;
		auto const inR = richardson(
			[ & ]( double d ) { return psiAt( probe.r + d, probe.z ); }, h );
		auto const inZ = richardson(
			[ & ]( double d ) { return psiAt( probe.r, probe.z + d ); }, h );

		// Relative to the size of the gradient rather than to either
		// component, so that a component which is small by symmetry does not
		// dominate the measure.
		double const scale = std::fabs( dR ) + std::fabs( dZ );
		double const plain =
			( std::fabs( inR.first - dR ) + std::fabs( inZ.first - dZ ) )/scale;
		double const rich =
			( std::fabs( inR.second - dR )
			  + std::fabs( inZ.second - dZ ) )/scale;

		std::printf( "    %-28s %12.3e %12.3e\n", probe.name, plain, rich );
		worstPlain = std::max( worstPlain, plain );
		worstRich = std::max( worstRich, rich );
	}
	std::fflush( stdout );

	BOOST_TEST( worstRich < 1.0e-9,
	            "the analytic gradient disagrees with a Richardson-extrapolated "
	            "difference of psi by " << worstRich << ". These are the same "
	            "function differentiated two ways, so a disagreement is a "
	            "defect in the chain rule and not a tolerance to widen" );

	// THE CONTROL. If the plain difference were as good, the extrapolation
	// would be decoration and this file would be claiming something it had not
	// shown.
	BOOST_TEST( worstPlain > 100.0*worstRich,
	            "the plain central difference is " << worstPlain
	            << " against the extrapolated " << worstRich
	            << ", so the two are comparable and the Richardson step is "
	            "buying nothing here. Either h has been chosen where the "
	            "truncation and the round-off happen to cross, or the "
	            "comparison is no longer measuring what it claims" );
}

/*
 * meq::CurrentFilament AGAINST tests/analytic/CurrentLoop.hpp.
 *
 * Two implementations of one closed form, sharing no code: the fixture writes
 * the textbook  psi = ( mu0 I/2 pi ) d [ ( 1 - k^2/2 )K - E ]  with
 * std::comp_ellint of the MODULUS, and the library writes Carlson's symmetric
 * forms of the COMPLEMENTARY modulus squared. So this is the check that the
 * promotion is faithful, and it is the same shape as
 * SolovievGeometryConvergence checking coefficients the solver also uses.
 *
 * NOT ROUND-OFF, AND THE HEADER ALREADY SAID SO. Coils.hpp records a worst
 * 1.4e-12 relative between these two over the benchmark box, and that is a
 * property of the two ALGEBRAIC arrangements rather than of either
 * implementation -- ( 1 - k^2/2 )K - E is formed by cancellation and Carlson's
 * k^2( R_D/3 - R_F/2 ) is not. Asserting round-off here would be asserting
 * something known to be false; 1e-11 is the bar, and where they genuinely
 * diverge -- at the conductor -- Coils.hpp measures that THIS file is the one
 * that is right.
 */
BOOST_AUTO_TEST_CASE( the_filament_agrees_with_the_analytic_fixture )
{
	double const radius = 1.5;
	double const height = 0.1;
	double const current = 1.0e6;

	CurrentFilament const filament( radius, height, current );
	meq::analytic::CurrentLoop const loop( radius, height, current );

	std::printf( "\n  meq::CurrentFilament against analytic::CurrentLoop\n" );
	std::printf( "    %8s %8s %14s %12s %12s\n",
	             "r", "z", "psi", "rel psi", "rel grad" );

	double worstPsi = 0.0;
	double worstGrad = 0.0;

	for ( double r : { 0.4, 0.9, 2.2, 3.5 } )
	{
		for ( double z : { -1.1, 0.35, 1.7 } )
		{
			double const mine = meq::filamentPsi( filament, r, z );
			double const theirs = loop.psi( r, z );
			double const scale = std::max( std::fabs( theirs ), 1.0e-300 );
			double const relPsi = std::fabs( mine - theirs )/scale;

			double dR = 0.0;
			double dZ = 0.0;
			meq::filamentGradPsi( filament, r, z, dR, dZ );
			double tR = 0.0;
			double tZ = 0.0;
			loop.gradPsi( r, z, tR, tZ );
			double const gradScale =
				std::max( std::fabs( tR ) + std::fabs( tZ ), 1.0e-300 );
			double const relGrad =
				( std::fabs( dR - tR ) + std::fabs( dZ - tZ ) )/gradScale;

			std::printf( "    %8.2f %8.2f %14.6e %12.3e %12.3e\n",
			             r, z, mine, relPsi, relGrad );
			worstPsi = std::max( worstPsi, relPsi );
			worstGrad = std::max( worstGrad, relGrad );
		}
	}
	std::fflush( stdout );

	BOOST_TEST( worstPsi < 1.0e-11,
	            "meq::filamentPsi and analytic::CurrentLoop::psi differ by "
	            << worstPsi << ". They are two arrangements of one closed form, "
	            "so this is a transcription error rather than a tolerance" );
	BOOST_TEST( worstGrad < 1.0e-11,
	            "the two gradients differ by " << worstGrad
	            << ". The library's is Carlson throughout and the fixture's "
	            "divides by 1 - k^2; away from the conductor they must still "
	            "agree, and near it Coils.hpp records which one to believe" );
}

/*
 * THE AXIS, WHERE THE LIBRARY REACHES A LIMIT THE FIXTURE DOES NOT.
 *
 * This is the one place the two DELIBERATELY differ, so it is asserted rather
 * than left to be discovered. CurrentLoop.hpp writes
 * dk/dr = k[ 1/( 2r ) - ( a + r )/d^2 ] literally and says of it: *"the
 * 1/( 2r ) is why this is NaN at r = 0. It is a real 1/r and not an artefact:
 * psi ~ r^2 there, so d psi/d r ~ r and the limit exists, but the expression as
 * written does not reach it."*
 *
 * The library's arrangement reaches it. Both brackets carry k^2 as an explicit
 * factor and k^2 = 4ar/d^2 is exactly zero on the axis, so grad_bar psi is
 * ( 0, 0 ) BIT EXACTLY -- not to round-off. That is worth having exactly for
 * the same reason psi( 0, z ) = 0 is: it is the boundary condition the
 * free-boundary problem imposes there.
 *
 * AND THE FLUX q IS STILL NaN, WHICH IS A SEPARATE STATEMENT. q divides that
 * exact zero by r, and 0/0 is not 0. The limit is finite and this does not
 * reach it. It matters because a semicircle centred on the axis MEETS the axis
 * at both ends, so a sweep of such a Gamma samples exactly the point where q is
 * unavailable -- a caller must handle its endpoints rather than discover NaN in
 * a quadrature sum.
 */
BOOST_AUTO_TEST_CASE( the_gradient_is_exactly_zero_on_the_axis_where_the_flux_is_not )
{
	CurrentFilament const filament( 1.5, 0.1, 1.0e6 );

	for ( double z : { -0.8, 0.0, 0.1, 2.0 } )
	{
		BOOST_TEST( meq::filamentPsi( filament, 0.0, z ) == 0.0,
		            "psi is not exactly zero on the axis at z = " << z );

		double dR = 1.0;
		double dZ = 1.0;
		meq::filamentGradPsi( filament, 0.0, z, dR, dZ );
		BOOST_TEST( dR == 0.0,
		            "d psi/d r on the axis at z = " << z << " is " << dR
		            << " and not exactly 0.0. Both brackets carry k^2, which is "
		            "exactly zero there, so anything else means the factored "
		            "form has been replaced by the literal chain rule -- which "
		            "is NaN here, as CurrentLoop.hpp records" );
		BOOST_TEST( dZ == 0.0,
		            "d psi/d z on the axis at z = " << z << " is " << dZ );

		// AND THE FLUX IS NOT. 0/0, and deliberately not special-cased.
		double qR = 0.0;
		double qZ = 0.0;
		meq::filamentFlux( filament, 0.0, z, qR, qZ );
		BOOST_TEST( std::isnan( qR ),
		            "q_r on the axis is " << qR << " rather than NaN. If this "
		            "has become finite, someone has substituted the limit -- "
		            "which is a kindness that hides from a caller sweeping a "
		            "Gamma that MEETS the axis that its endpoints are special" );
		BOOST_TEST( std::isnan( qZ ), "q_z on the axis is " << qZ );
	}

	// The fixture's behaviour, asserted so that the difference above is a
	// measured contrast and not an assumption about someone else's code.
	meq::analytic::CurrentLoop const loop( 1.5, 0.1, 1.0e6 );
	double tR = 0.0;
	double tZ = 0.0;
	loop.gradPsi( 0.0, 0.3, tR, tZ );
	BOOST_TEST( std::isnan( tR ),
	            "analytic::CurrentLoop::dPsiDr is no longer NaN on the axis, so "
	            "the contrast this case draws has gone away and its comment is "
	            "stale" );
}

/*
 * A FILAMENT AS THE LIMIT OF A SHRINKING COIL -- FOR THE GRADIENT THIS TIME.
 *
 * a_thin_coil_approaches_the_filament_at_second_order does this for psi, and
 * took its limit from a very thin COIL rather than from a second
 * implementation, deliberately: it wanted a statement about the cross-section
 * and not about two transcriptions.
 *
 * THIS ONE TAKES THE LIMIT FROM meq::CurrentFilament, and that is the point.
 * The two classes are separate code -- one integrates the derivative kernel
 * over a rectangle at a quadrature order, the other evaluates it once -- so a
 * convention mismatch between them is invisible to either alone. A factor of
 * the area, a sign, a missing mu0, a current density confused with a current:
 * every one of those leaves both classes internally consistent and makes them
 * disagree here. It is the same argument the psi case makes for its own
 * existence, applied across the seam instead of along it.
 */
BOOST_AUTO_TEST_CASE( a_shrinking_coils_gradient_approaches_the_filaments )
{
	double const current = 1.0e6;
	double const centreR = 1.5;
	double const centreZ = 0.1;
	double const fieldR = 2.4;
	double const fieldZ = 0.7;

	CurrentFilament const filament( centreR, centreZ, current );
	double limitR = 0.0;
	double limitZ = 0.0;
	meq::filamentGradPsi( filament, fieldR, fieldZ, limitR, limitZ );
	double const limitScale = std::fabs( limitR ) + std::fabs( limitZ );

	std::printf( "\n  a shrinking coil's GRADIENT against the filament's\n" );
	std::printf( "    %10s %16s %16s %12s %10s\n",
	             "half-size", "d psi/d r", "d psi/d z", "rel", "ratio" );

	double previousDifference = 0.0;
	bool first = true;
	double worst = 0.0;

	for ( double halfSize : { 0.16, 0.08, 0.04, 0.02, 0.01 } )
	{
		Coil const coil( centreR, centreZ, halfSize, halfSize, current );
		double dR = 0.0;
		double dZ = 0.0;
		meq::coilGradPsi( coil, fieldR, fieldZ, dR, dZ );

		double const difference =
			( std::fabs( dR - limitR ) + std::fabs( dZ - limitZ ) )/limitScale;
		double ratio = 0.0;
		if ( !first && difference > 0.0 )
			ratio = previousDifference/difference;

		std::printf( "    %10.4f %16.8e %16.8e %12.3e %10.3f\n",
		             halfSize, dR, dZ, difference, ratio );

		if ( !first && previousDifference > 0.0 && difference > 0.0 )
			BOOST_TEST( ratio > 3.0,
			            "halving the cross-section reduced the gradient's "
			            "difference from the filament by only " << ratio
			            << "; second order demands about 4. A ratio near 1 "
			            "means the two classes disagree by a CONSTANT, which "
			            "is a convention mismatch -- an area, a sign, a mu0 -- "
			            "and not a cross-section effect" );

		previousDifference = difference;
		worst = difference;
		first = false;
	}
	std::fflush( stdout );

	// AND IT ACTUALLY ARRIVES. A rate alone is satisfied by two quantities
	// converging to different limits at second order.
	BOOST_TEST( worst < 1.0e-4,
	            "the thinnest coil's gradient is still " << worst
	            << " from the filament's. The RATE above can be met by two "
	            "sequences approaching different limits, so this is the half "
	            "that says they approach the SAME one" );
}

/*
 * AMPERE'S LAW ON THE FILAMENT, which is FB-2's acceptance identity applied to
 * the class that has no cross-section.
 *
 *     oint ( 1/r ) dpsi/dn dl = -mu0 I
 *
 * Coils.hpp derives the sign and the_outward_flux_is_minus_mu0_times_the_total_current
 * pins it for a rectangle, at 1.07e-08 by central differences of psi. THIS one
 * uses the ANALYTIC gradient, so the only error left is the contour quadrature
 * -- which is why it reaches so much further, and which makes it a check on the
 * gradient rather than on the difference stencil.
 *
 * The contour encloses the filament, and the identity is exact for ANY contour
 * that does: it is Ampere's law, and what is enclosed is the whole current.
 */
BOOST_AUTO_TEST_CASE( amperes_law_holds_on_the_filament )
{
	double const current = 1.0e6;
	CurrentFilament const filament( 1.5, 0.1, current );
	double const mu0 = meq::vacuumPermeability;

	double const rLo = 0.7;
	double const rHi = 2.6;
	double const zLo = -0.9;
	double const zHi = 1.2;

	auto circulation = [ & ]( int panels )
	{
		double total = 0.0;
		double const dr = ( rHi - rLo )/panels;
		double const dz = ( zHi - zLo )/panels;
		for ( int i = 0; i < panels; ++i )
		{
			double const r = rLo + ( i + 0.5 )*dr;
			double const z = zLo + ( i + 0.5 )*dz;
			double dR = 0.0;
			double dZ = 0.0;

			// The two vertical sides: n = +/- r-hat, dl = dz.
			meq::filamentGradPsi( filament, rHi, z, dR, dZ );
			total += ( dR/rHi )*dz;
			meq::filamentGradPsi( filament, rLo, z, dR, dZ );
			total += ( -dR/rLo )*dz;

			// The two horizontal sides: n = +/- z-hat, dl = dr.
			meq::filamentGradPsi( filament, r, zHi, dR, dZ );
			total += ( dZ/r )*dr;
			meq::filamentGradPsi( filament, r, zLo, dR, dZ );
			total += ( -dZ/r )*dr;
		}
		return total;
	};

	double const coarse = circulation( 400 );
	double const fine = circulation( 800 );
	// Midpoint is second order in the panel count, so Richardson is
	// ( 4 fine - coarse )/3 -- the same step the gradient test needs and for
	// the same reason.
	double const extrapolated = ( 4.0*fine - coarse )/3.0;
	double const expected = -mu0*current;

	std::printf( "\n  Ampere's law on a filament: oint ( 1/r ) dpsi/dn dl\n" );
	std::printf( "    400 panels        %18.10e\n", coarse );
	std::printf( "    800 panels        %18.10e\n", fine );
	std::printf( "    Richardson        %18.10e\n", extrapolated );
	std::printf( "    -mu0 I            %18.10e\n", expected );
	std::printf( "    relative          %18.3e\n",
	             std::fabs( extrapolated - expected )/std::fabs( expected ) );
	std::fflush( stdout );

	BOOST_TEST( std::fabs( extrapolated - expected )
	            < 1.0e-9*std::fabs( expected ),
	            "the circulation is " << extrapolated << " against -mu0 I = "
	            << expected << ". THE SIGN IS THE THING TO CHECK FIRST: "
	            "Coils.hpp records that the same identity as a circulation of B "
	            "counterclockwise comes out POSITIVE, because phi-hat = "
	            "z-hat x r-hat, and FREE-BOUNDARY-PLAN.md section 7 predicts it "
	            "being got wrong at least once" );
}
