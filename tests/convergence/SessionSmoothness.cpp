#define BOOST_TEST_MODULE SessionSmoothness
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem/config/config.hpp"

#include "meq/EquilibriumSnapshot.hpp"
#include "meq/Run.hpp"

#include "SessionCases.hpp"
#include "SourceTree.hpp"

/*
 * IS THE GEOMETRY A CONTINUOUS FUNCTION OF THE TABLES? -- MANTA-EQUILIBRIUM-
 * PLAN.md sections 4.5 and 4.6.
 *
 * MaNTA differences its equilibrium model at relative steps of 1e-8, so G
 * must be continuous there, and the noise of a resolve below 1e-12 of its
 * value. What breaks continuity is measured (M-363): a flux-surface node
 * whose ray root crosses an element face, onto a broken field's other
 * polynomial. A scan that happens to straddle no crossing passes whatever the
 * extraction does, so this test FINDS one: Omega is stepped from the base
 * until some node of some surface changes element (FluxSurface::elements),
 * the step is bisected on that node's element to 1e-9, and G on that surface
 * must not jump across it -- its change there, less the smooth change, below
 * 1e-11 of its value. And the same tables from two warm starts must agree to
 * 1e-12.
 *
 * On rotating-normalised's box alone, about a minute. fixed-d-tcv's MXH
 * boundary, its band and the long sweep M-SESSION-4 was taken from are
 * ExtractionJumps' cases (sweepBox, sweepTcv), built and not registered.
 */
namespace
{
	using namespace meq::test::session;

	using Family = std::shared_ptr<meq::FluxSurfaceFamily const>;

	std::vector<double> quantities( meq::FluxSurfaceFamily const &family, std::size_t i )
	{
		meq::FluxSurface const &s = family.surfaces[ i ];
		double const rho = std::sqrt( s.normalisedFlux );
		double const dRhoDPsi = 1.0/( 2.0*rho*( family.psiBoundary - family.psiAxis ) );
		return { s.vPrime/std::abs( dRhoDPsi ), s.gradPsiSquared*dRhoDPsi*dRhoDPsi, s.inverseRSquared,
		         s.gradPsiSquaredOverRSquared };
	}

	char const *quantityName( std::size_t q )
	{
		static char const *const names[] = { "V'_rho", "<|grad rho|^2>", "<R^-2>", "<|grad psi|^2/R^2>" };
		return names[ q ];
	}

	meq::ProfileTables tables( double h )
	{
		return meq::ProfileTables{ { "Omega", omega( 4.0e5*( 1.0 + h ) ) },
		                           { "species.D.Density", density( 0.0 ) } };
	}

	/// The first node, on any surface, whose element differs between @a a
	/// and @a b; false if none does.
	bool crossing( meq::FluxSurfaceFamily const &a, meq::FluxSurfaceFamily const &b,
	               std::size_t &surface, std::size_t &node )
	{
		for ( surface = 0; surface < a.surfaces.size(); ++surface )
			for ( node = 0; node < a.surfaces[ surface ].count(); ++node )
				if ( a.surfaces[ surface ].elements[ node ] != b.surfaces[ surface ].elements[ node ] )
					return true;
		return false;
	}
}

BOOST_AUTO_TEST_CASE( theGeometryIsContinuousWhereANodeCrossesAFace )
{
	std::string const path = meq::test::sourcePath( "examples/rotating-normalised.toml" );
	// TO ROUND-OFF, as a session serving MaNTA solves: at the shipped 1e-10 the
	// noise across warm starts reaches 2.8e-11 (M-SESSION-4).
	auto run = solvedRun( path, adapted( meq::test::shippedConfig( "examples/rotating-normalised.toml" ),
	                                     1.0e-16 ) );
	auto at = [ & ]( double h, std::shared_ptr<meq::EquilibriumSnapshot const> const &warm )
	{
		BOOST_TEST_REQUIRE( run->resolve( tables( h ), warm ) == meq::Solved );
		return run->snapshot();
	};
	std::shared_ptr<meq::EquilibriumSnapshot const> const base = at( 0.0, run->snapshot() );

	// THE NOISE: the base tables from the base and from a state 1e-4 away.
	Family const fromBase = at( 0.0, base )->family;
	Family const fromAway = at( 0.0, at( 1.0e-4, base ) )->family;
	double floor = 0.0;
	for ( std::size_t i = 0; i < fromBase->surfaces.size(); ++i )
	{
		std::vector<double> const a = quantities( *fromBase, i ), b = quantities( *fromAway, i );
		for ( std::size_t q = 0; q < a.size(); ++q )
			floor = std::max( floor, std::abs( b[ q ] - a[ q ] )/std::abs( a[ q ] ) );
	}
	std::printf( "\n  the same tables from two warm starts differ by %.2e of the value\n", floor );
	// REPEATABILITY, so max( 1e-12, 20 eps ): two warm starts reach the same
	// tables to the solve's own accuracy, ten eps, and the bound is twice that,
	// as QDrivenScalePin's is. In double 20 eps is 4.4e-15 and 1e-12 stands --
	// MaNTA's contract, which is double (MANTA-EQUILIBRIUM-PLAN.md, D5); in
	// single it is 2.4e-6. The jump below is not a repeatability bound.
	double const twentyEps = 20.0*static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() );
	BOOST_TEST( floor < std::max( 1.0e-12, twentyEps ),
		"the same tables solved from two warm starts differ by " << floor << " of the value, so "
		"a difference at MaNTA's 1e-8 step carries " << floor/1.0e-8 << " of noise. A session must "
		"converge to round-off, not to a relative tolerance (NewtonRelativeTolerance below ten eps "
		"stops at the round-off floor, M-306)." );

	// A CROSSING: Omega stepped by 1e-4 until a node changes element. The box
	// meets about 1.4 a step (M-363), so twenty steps is a wide margin.
	double const step = 1.0e-4;
	double lo = 0.0, hi = 0.0;
	Family left = at( 0.0, base )->family, right;
	std::size_t surface = 0, node = 0;
	bool found = false;
	for ( int j = 1; j <= 20 && !found; ++j )
	{
		right = at( j*step, base )->family;
		found = crossing( *left, *right, surface, node );
		if ( !found )
		{
			left = right;
			lo = j*step;
		}
		hi = j*step;
	}
	BOOST_TEST_REQUIRE( found, "no node of any surface changed element over 2e-3 of Omega, so the "
	                           "test has no crossing to look at; widen the scan" );

	// BISECTED to 1e-9 on that node's element.
	int const leftElement = left->surfaces[ surface ].elements[ node ];
	std::vector<double> const coarseLeft = quantities( *left, surface );
	std::vector<double> const coarseRight = quantities( *right, surface );
	double const coarseWidth = hi - lo;
	while ( hi - lo > 1.0e-9 )
	{
		double const mid = 0.5*( lo + hi );
		Family const f = at( mid, base )->family;
		if ( f->surfaces[ surface ].elements[ node ] == leftElement )
		{
			left = f;
			lo = mid;
		}
		else
		{
			right = f;
			hi = mid;
		}
	}
	int const rightElement = right->surfaces[ surface ].elements[ node ];
	std::vector<double> const a = quantities( *left, surface ), b = quantities( *right, surface );
	double jump = 0.0;
	std::size_t worst = 0;
	for ( std::size_t q = 0; q < a.size(); ++q )
	{
		// Less the smooth change over the interval, from the coarse step's slope.
		double const smooth = ( coarseRight[ q ] - coarseLeft[ q ] )*( hi - lo )/coarseWidth;
		double const change = std::abs( b[ q ] - a[ q ] - smooth )/std::abs( a[ q ] );
		if ( change > jump )
		{
			jump = change;
			worst = q;
		}
	}
	std::printf( "  surface %zu's node %zu crosses from element %d to %d between Omega x ( 1 %+.10e ) "
	             "and ( 1 %+.10e ); across it G moves by %.2e of %s beyond its smooth change\n",
	             surface, node, leftElement, rightElement, lo, hi, jump, quantityName( worst ) );
	BOOST_TEST( jump < 1.0e-11,
		"G jumps by " << jump << " of " << quantityName( worst ) << " on surface " << surface
		<< " where its node " << node << " crosses from element " << leftElement << " to "
		<< rightElement << ", within 1e-9 of Omega. A central difference at MaNTA's 1.5e-8 straddling "
		"it is wrong by the jump over the step (M-363). What removes it is an extraction "
		"continuous in the solved state -- MANTA-EQUILIBRIUM-PLAN.md section 4.5, M-364." );
}
