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

using meq::Coil;
using meq::CoilSet;

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
