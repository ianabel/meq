#define BOOST_TEST_MODULE ExtractionContinuity
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

#include "ExtractionCrossing.hpp"

/*
 * THE FLUX-SURFACE GEOMETRY DOES NOT JUMP WHERE A FAMILY NODE CROSSES AN
 * ELEMENT FACE -- MEASUREMENTS.md M-364 to M-367, CLAUDE_INVERSION.md IN-0.
 *
 * On psi* and q*, broken across every face, a node whose ray root slides
 * across one lands on the neighbour's polynomial and the family's averages
 * jump (M-364): inside the discretisation's accuracy, but divided by the step
 * of anything that differences the geometry. The extraction's default pair,
 * psi~ and q~ (ContinuousField.hpp), is continuous, and the family on it must
 * not jump. A scan that happens to straddle no crossing passes whatever the
 * extraction does, so this test FINDS one: a Solov'ev-type source
 * F = -( ( 1 - A ) R^2 + A ) on the standard box, psi = 0 on its edge, the
 * mix A stepped until some node of some surface changes element
 * (FluxSurface::elements), the step bisected on that node's element to 1e-9,
 * and the family's V', < R^-2 >, < |grad psi|^2/R^2 > and < |grad psi|^2 >
 * on that surface compared across it, less their smooth change
 * (ExtractionCrossing.hpp has the fixture, ExtractionJumps the picture).
 *
 * The control is the same procedure on psi* and q*, where the jump must be
 * there: without it this case could pass on an instrument that sees nothing.
 * About 20 s; every state a linear solve on n = 12, k = 2.
 */
namespace
{
	std::vector<double> quantities( meq::FluxSurface const &s )
	{
		return { s.vPrime, s.inverseRSquared, s.gradPsiSquaredOverRSquared, s.gradPsiSquared };
	}

	struct Crossing
	{
		bool found = false;
		std::size_t surface = 0, node = 0;
		int from = 0, to = 0;
		double lo = 0.0, hi = 0.0;
		double jump = 0.0;   ///< relative, worst quantity, beyond the smooth change
	};

	/// Step A from @a start by @a step until a node changes element, bisect
	/// that node's element to an interval of 1e-9 in A, and measure the jump.
	Crossing measure( meq::Potential which, double start, double step )
	{
		Crossing out;
		auto familyAt = [ & ]( double mix )
		{
			return meq::tests::crossing::familyAt( mix, 12, 2, which );
		};
		meq::FluxSurfaceFamily left = familyAt( start ), right;
		double lo = start, hi = start;
		for ( int j = 1; j <= 40 && !out.found; ++j )
		{
			hi = start + j*step;
			right = familyAt( hi );
			for ( std::size_t i = 0; i < left.surfaces.size() && !out.found; ++i )
				for ( std::size_t n = 0; n < left.surfaces[ i ].count(); ++n )
					if ( left.surfaces[ i ].elements[ n ] != right.surfaces[ i ].elements[ n ] )
					{
						out.found = true;
						out.surface = i;
						out.node = n;
						break;
					}
			if ( !out.found )
			{
				left = right;
				lo = hi;
			}
		}
		if ( !out.found )
			return out;
		std::vector<double> const coarseLeft = quantities( left.surfaces[ out.surface ] );
		std::vector<double> const coarseRight = quantities( right.surfaces[ out.surface ] );
		double const coarseWidth = hi - lo;
		out.from = left.surfaces[ out.surface ].elements[ out.node ];
		while ( hi - lo > 1.0e-9 )
		{
			double const mid = 0.5*( lo + hi );
			meq::FluxSurfaceFamily f = familyAt( mid );
			if ( f.surfaces[ out.surface ].elements[ out.node ] == out.from )
			{
				left = std::move( f );
				lo = mid;
			}
			else
			{
				right = std::move( f );
				hi = mid;
			}
		}
		out.to = right.surfaces[ out.surface ].elements[ out.node ];
		out.lo = lo;
		out.hi = hi;
		std::vector<double> const a = quantities( left.surfaces[ out.surface ] );
		std::vector<double> const b = quantities( right.surfaces[ out.surface ] );
		for ( std::size_t q = 0; q < a.size(); ++q )
		{
			double const smooth = ( coarseRight[ q ] - coarseLeft[ q ] )*( hi - lo )/coarseWidth;
			out.jump = std::max( out.jump, std::abs( b[ q ] - a[ q ] - smooth )/std::abs( a[ q ] ) );
		}
		return out;
	}

	void report( char const *name, Crossing const &c )
	{
		std::printf( "  %s: surface %zu's node %zu crosses element %d -> %d within A in "
		             "[ %.10f, %.10f ]; across it the family moves %.2e beyond its smooth change\n",
		             name, c.surface, c.node, c.from, c.to, c.lo, c.hi, c.jump );
	}

	double const start = 0.30;
	double const step = 1.0e-3;
}

BOOST_AUTO_TEST_CASE( theFamilyIsContinuousWhereANodeCrossesAFace )
{
	Crossing const c = measure( meq::Potential::Continuous, start, step );
	BOOST_TEST_REQUIRE( c.found, "no node of any surface changed element over A in [ 0.30, 0.34 ]; "
	                             "the test has no crossing to look at, so widen the scan" );
	report( "psi~ and q~", c );
	// max( 1e-11, 20 eps ), the rule QDrivenScalePin's repeatability bounds
	// follow: in double 20 eps is 4.4e-15 and 1e-11 stands; in single it is
	// 2.4e-6. Two states 1e-9 apart differ by their round-off as well as by any
	// jump, and in single the family's own is about 1e-5 (M-SX-1), so this is
	// red there until that floor comes down.
	double const twentyEps = 20.0*static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() );
	BOOST_TEST( c.jump < std::max( 1.0e-11, twentyEps ),
		"the family jumps by " << c.jump << " where surface " << c.surface << "'s node " << c.node
		<< " crosses from element " << c.from << " to " << c.to << ". On the continuous pair it "
		"must not: a jump is divided by the step of anything that differences the geometry "
		"(MEASUREMENTS.md M-364). Suspect the fields' continuity first (ContinuousFieldConvergence)" );
}

/*
 * AND THE FAMILY IS LABELLED ON THE FIELD IT ROOTS -- meq::labelledAxis(),
 * MEASUREMENTS.md M-SX-1. With no psi_ax unknown in the solve, Psi_N must be
 * zero at the axis of psi~, the field the levels are in: the finder's psi_h
 * there differs by the discretisation's error at the axis, which entered every
 * inner label as an offset divided by Psi_N, and in single by the element-local
 * recovery's round-off, 27 times psi~'s. The normalised half -- labelled by the
 * solve's own psi_ax -- is SafetyFactorSolver's
 * theReferenceIsAFixedPointOfItsOwnTarget, which no other label passes.
 */
BOOST_AUTO_TEST_CASE( theFamilyIsLabelledOnTheFieldItRoots )
{
	meq::tests::crossing::Solved s( start, 12, 2 );
	meq::FluxSurfaceFamily const family = s.family( meq::Potential::Continuous );
	mfem::DenseMatrix point( 2, 1 );
	point( 0, 0 ) = family.axisR;
	point( 1, 0 ) = family.axisZ;
	mfem::Array<int> elements;
	mfem::Array<mfem::IntegrationPoint> ips;
	s.theMesh().FindPoints( point, elements, ips, false );
	BOOST_TEST_REQUIRE( elements[ 0 ] >= 0 );
	double const tilde = s.theSolver().continuousPotential().GetValue( elements[ 0 ], ips[ 0 ] );
	double const finder = s.theSolver().potential().GetValue( elements[ 0 ], ips[ 0 ] );
	double const span = std::abs( family.psiAxis - family.psiBoundary );
	std::printf( "  the family's psi_ax less psi~ at its axis: %.2e of the span; psi_h there differs from "
	             "psi~ by %.2e\n", ( family.psiAxis - tilde )/span, ( finder - tilde )/span );
	BOOST_TEST( std::abs( family.psiAxis - tilde ) < 64.0*std::numeric_limits<mfem::real_t>::epsilon()*span,
		"the family's psi_ax is " << family.psiAxis << " where psi~ at its axis is " << tilde
		<< " (psi_h there " << finder << "): Psi_N is not zero on the axis of the field its levels are in. "
		"meq::labelledAxis() is what reads it there" );
}

/*
 * THE CONTROL. The same procedure on psi* and q*, where M-364 measured the
 * jump: the instrument must see it, or the case above proves nothing.
 */
BOOST_AUTO_TEST_CASE( theBrokenPairJumpsThereAndTheTestSeesIt )
{
	Crossing const c = measure( meq::Potential::PostProcessed, start, step );
	BOOST_TEST_REQUIRE( c.found, "no node crossed a face on psi*; widen the scan" );
	report( "psi* and q*", c );
	BOOST_TEST( c.jump > 1.0e-10,
		"on psi* and q* the family moves only " << c.jump << " across a node's face crossing; "
		"the broken pair jumps there by the face jump (M-364), so an instrument that cannot see "
		"it cannot see one on the continuous pair either" );
}
