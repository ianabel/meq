/*
 * WHAT EACH ELEMENT IS MADE OF -- ELEMENT-REGIONS-PLAN.md ER-0.
 *
 * meq::GradShafranovSolver::elementMaterials() asks the solver's own support
 * decision at the solver's own quadrature points, and assigns each point to
 * vacuum, plasma or a conductor. Every case here asserts it against something
 * it cannot fake:
 *
 *   1. thePlasmaAreaConvergesOnTheExactSoloviev -- the partition sums to one
 *      in every element, the plasma area converges to the area of
 *      { psi < 0 } computed from the closed form by a quadrature of its own,
 *      and the plasma_boundary count grows like 1/h.
 *   2. theFillAndTheExclusionsAreTheIntegratorsOwn -- the private flux
 *      region below the X-point reads plasma_fraction 0 under the fill and 1
 *      without it; an excluded attribute reads 0.
 *   3. aPlasmaInsideACoilIsAFault -- posed deliberately, and its control.
 *   4. theConductorsAndTheCutoffOnTheMirror -- geq's machine: a meshed coil
 *      is pure; the same machine under a subtracting model cuts copper cells
 *      on every edge of every pack; AlignConductors makes them pure again;
 *      the aligned cutoff is straddled by no element and M-262's unaligned
 *      one is.
 *   5. thePredicateIsTheIntegratedOne -- `F/R` summed over the points called
 *      plasma is the assembled plasma current, under a conductor split.
 *   6. alignConductorsCost -- the element count of the key on the shipped
 *      mirror and limited-tokamak-filament, through the driver's own mesh
 *      command, and the purity it buys on the second machine.
 *
 * NO SOLVE RUNS. Each field is a closed form projected onto the potential
 * space after prepare(), or zero under a split whose psi_c is the vacuum
 * field, and the classifier reads it as it would read a converged answer.
 */

#define BOOST_TEST_MODULE ElementRegions
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "MirrorMachine.hpp"

#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Materials.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"

using meq::ElementMaterials;
using meq::ElementRegion;
using meq::GradShafranovSolver;
using meq::MaterialTable;
using meq::analytic::SolovievEquilibrium;

namespace
{
	/// The worst | sum of an element's fractions - 1 | over the mesh.
	double worstPartition( ElementMaterials const &record )
	{
		double worst = 0.0;
		for ( int e = 0; e < record.elementCount(); ++e )
		{
			double sum = 0.0;
			for ( int i = record.offsets[ static_cast< std::size_t >( e ) ];
			      i < record.offsets[ static_cast< std::size_t >( e ) + 1 ]; ++i )
				sum += record.entries[ static_cast< std::size_t >( i ) ].fraction;
			worst = std::max( worst, std::abs( sum - 1.0 ) );
		}
		return worst;
	}

	double totalArea( ElementMaterials const &record )
	{
		double sum = 0.0;
		for ( double a : record.area )
			sum += a;
		return sum;
	}

	void printRegions( char const *label, ElementMaterials const &record )
	{
		std::printf( "    %-34s %6d elements:", label, record.elementCount() );
		for ( int r = 0; r < meq::elementRegionCount; ++r )
		{
			int const n = record.regionCount( static_cast< ElementRegion >( r ) );
			if ( n > 0 )
				std::printf( " %s %d", meq::elementRegionName(
					static_cast< ElementRegion >( r ) ), n );
		}
		std::printf( "; faults %d\n", record.faults );
		std::fflush( stdout );
	}

	/// The vertices' extent in R and z.
	struct Extent
	{
		double rMin = 1.0e300, rMax = -1.0e300, zMin = 1.0e300, zMax = -1.0e300;
	};

	Extent extentOf( mfem::Mesh const &mesh, int e )
	{
		Extent x;
		mfem::Array< int > vertices;
		mesh.GetElementVertices( e, vertices );
		for ( int i = 0; i < vertices.Size(); ++i )
		{
			mfem::real_t const *v = mesh.GetVertex( vertices[ i ] );
			x.rMin = std::min< double >( x.rMin, v[ 0 ] );
			x.rMax = std::max< double >( x.rMax, v[ 0 ] );
			x.zMin = std::min< double >( x.zMin, v[ 1 ] );
			x.zMax = std::max< double >( x.zMax, v[ 1 ] );
		}
		return x;
	}

	/// Whether an element's vertices lie strictly on both sides of z = line.
	bool straddlesZ( mfem::Mesh const &mesh, int e, double line )
	{
		Extent const x = extentOf( mesh, e );
		return x.zMin < line && x.zMax > line;
	}

	/// Of a rectangle's four edges, how many have a conductor_boundary element
	/// of @a material on them: one whose vertices straddle the edge's line
	/// strictly and overlap the edge's extent along it.
	int edgesCut( mfem::Mesh const &mesh, ElementMaterials const &record,
	              int material, meq::Coil const &coil )
	{
		std::array< bool, 4 > cut{ false, false, false, false };
		for ( int e = 0; e < record.elementCount(); ++e )
		{
			if ( record.region[ static_cast< std::size_t >( e ) ]
			     != ElementRegion::ConductorBoundary
			     || record.fraction( e, material ) <= 0.0 )
				continue;
			Extent const x = extentOf( mesh, e );
			bool const alongR = x.rMax > coil.minRadius() && x.rMin < coil.maxRadius();
			bool const alongZ = x.zMax > coil.zMin() && x.zMin < coil.zMax();
			if ( alongZ && x.rMin < coil.minRadius() && x.rMax > coil.minRadius() )
				cut[ 0 ] = true;
			if ( alongZ && x.rMin < coil.maxRadius() && x.rMax > coil.maxRadius() )
				cut[ 1 ] = true;
			if ( alongR && x.zMin < coil.zMin() && x.zMax > coil.zMin() )
				cut[ 2 ] = true;
			if ( alongR && x.zMin < coil.zMax() && x.zMax > coil.zMax() )
				cut[ 3 ] = true;
		}
		return static_cast< int >( std::count( cut.begin(), cut.end(), true ) );
	}

	/// The full ( flux, potential, trace ) state with the potential the solver
	/// holds, which is what the border assemblers take.
	mfem::Vector stateOf( GradShafranovSolver &solver )
	{
		int const fluxSize = solver.fluxSpace().GetVSize();
		int const potentialSize = solver.potentialSpace().GetVSize();
		mfem::Vector state( fluxSize + potentialSize + solver.traceSpace().GetVSize() );
		state = 0.0;
		mfem::GridFunction const &psi = solver.potential();
		for ( int i = 0; i < potentialSize; ++i )
			state( fluxSize + i ) = psi( i );
		return state;
	}

	/// `int F_plasma/R` over the points elementMaterials() called plasma,
	/// summed in element order as assemblePlasmaCurrent() sums.
	double plasmaLoadOverClassifiedPoints(
		meq::NormalisedSource const &source,
		std::vector< GradShafranovSolver::ClassifiedPoint > const &points,
		int elements )
	{
		std::vector< double > perElement( static_cast< std::size_t >( elements ), 0.0 );
		for ( GradShafranovSolver::ClassifiedPoint const &p : points )
			if ( p.material == MaterialTable::plasma )
				perElement[ static_cast< std::size_t >( p.element ) ] +=
					p.weight*source.scaledF( p.radius, p.z, p.psi )/p.radius;
		double total = 0.0;
		for ( double v : perElement )
			total += v;
		return total;
	}

	// ---- the Solov'ev box ------------------------------------------------

	/// ITER Example 2: the separatrix is psi = 0, the axis an interior
	/// MINIMUM, and the X-point at ( 0.88384, -0.704 ).
	SolovievEquilibrium const &iter()
	{
		static SolovievEquilibrium const eq = SolovievEquilibrium::iterExample2();
		return eq;
	}

	/// A box of triangles holding ITER Example 2's psi, projected, with a
	/// confined one-sided source normalised to ( psi_ax, 0 ) and the fill
	/// refreshed at it -- the state a converged answer leaves.
	struct SolovievBox
	{
		mfem::Mesh mesh;
		meq::NormalisedMHDSource source;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr< GradShafranovSolver > solver;

		SolovievBox( int nR, int nZ, double r0, double r1, double z0, double z1,
		             int order,
		             std::vector< std::array< double, 4 > > const &excluded = {} )
			: mesh( mfem::Mesh::MakeCartesian2D( nR, nZ, mfem::Element::TRIANGLE,
			                                     false, r1 - r0, z1 - z0 ) ),
			  source( std::make_shared< meq::ConstantProfile const >( 1.0 ),
			          std::make_shared< meq::ConstantProfile const >( 0.0 ),
			          -0.04, 1.0 )
		{
			for ( int v = 0; v < mesh.GetNV(); ++v )
			{
				mesh.GetVertex( v )[ 0 ] += r0;
				mesh.GetVertex( v )[ 1 ] += z0;
			}

			// An excluded patch is attribute 2, by element centre.
			if ( !excluded.empty() )
			{
				for ( int e = 0; e < mesh.GetNE(); ++e )
				{
					mfem::Vector centre;
					mesh.GetElementCenter( e, centre );
					for ( std::array< double, 4 > const &box : excluded )
						if ( centre( 0 ) > box[ 0 ] && centre( 0 ) < box[ 1 ]
						     && centre( 1 ) > box[ 2 ] && centre( 1 ) < box[ 3 ] )
							mesh.SetAttribute( e, 2 );
				}
				mesh.SetAttributes();
			}

			solver = std::make_unique< GradShafranovSolver >( mesh, order );
			solver->setSource( source, -0.04 );
			solver->setBoundaryData( zero );
			if ( !excluded.empty() )
			{
				mfem::Array< int > attributes( 1 );
				attributes[ 0 ] = 2;
				solver->setPlasmaExclusion( attributes );
			}
			solver->prepare();

			mfem::FunctionCoefficient field( [ ]( mfem::Vector const &x )
			{
				return iter().psi( x( 0 ), x( 1 ) );
			} );
			solver->potential().ProjectCoefficient( field );

			// psi_ax: the core's minimum, which is all the one-sided test needs
			// the sign of.
			double axis = 0.0;
			for ( int i = 0; i <= 400; ++i )
				for ( int j = 0; j <= 400; ++j )
					axis = std::min( axis, iter().psi( 0.75 + 0.5*i/400.0,
					                                   -0.4 + 0.8*j/400.0 ) );
			source.setNormalisation( axis, 0.0 );
			source.setPlasmaSupport( true );
			solver->refreshPlasmaComponent( solver->potential() );
		}
	};

	/// The minimum of psi in R at fixed z, by golden section on the core's
	/// span. Unimodal there for every z the area integral visits.
	double coreMinimum( double z, double &where )
	{
		double a = 0.70;
		double b = 1.30;
		double const g = 0.5*( std::sqrt( 5.0 ) - 1.0 );
		double c = b - g*( b - a );
		double d = a + g*( b - a );
		while ( b - a > 1.0e-13 )
		{
			if ( iter().psi( c, z ) < iter().psi( d, z ) )
				b = d;
			else
				a = c;
			c = b - g*( b - a );
			d = a + g*( b - a );
		}
		where = 0.5*( a + b );
		return iter().psi( where, z );
	}

	/// The root of psi( ., z ) between @a a and @a b, by bisection.
	double rootInR( double z, double a, double b )
	{
		double fa = iter().psi( a, z );
		for ( int k = 0; k < 200 && b - a > 1.0e-15; ++k )
		{
			double const m = 0.5*( a + b );
			double const fm = iter().psi( m, z );
			if ( ( fm > 0.0 ) == ( fa > 0.0 ) )
			{
				a = m;
				fa = fm;
			}
			else
				b = m;
		}
		return 0.5*( a + b );
	}

	/**
	 * THE AREA OF { psi < 0 } ABOVE z = @a z0, FROM THE CLOSED FORM ALONE.
	 *
	 * The width w( z ) = R_out - R_in is a root pair per height, and the
	 * integral is taken in t with z = zTop - t^2: w vanishes like the square
	 * root of zTop - z at the high point, so in t the integrand is smooth and
	 * Gauss-Legendre converges spectrally. @a panels composite panels of
	 * order @a order.
	 */
	double exactArea( double z0, int panels, int order )
	{
		// zTop: where the core's minimum in R reaches zero.
		double lo = 0.5;
		double hi = 0.72;
		double where = 0.0;
		for ( int k = 0; k < 200 && hi - lo > 1.0e-15; ++k )
		{
			double const m = 0.5*( lo + hi );
			if ( coreMinimum( m, where ) < 0.0 )
				lo = m;
			else
				hi = m;
		}
		double const zTop = 0.5*( lo + hi );
		double const span = std::sqrt( zTop - z0 );

		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, order );
		double area = 0.0;
		for ( int p = 0; p < panels; ++p )
			for ( int q = 0; q < rule.GetNPoints(); ++q )
			{
				double const t = span*( p + rule.IntPoint( q ).x )/panels;
				double const z = zTop - t*t;
				double centre = 0.0;
				if ( coreMinimum( z, centre ) >= 0.0 )
					continue;
				double const width = rootInR( z, centre, 1.45 ) - rootInR( z, 0.55, centre );
				area += rule.IntPoint( q ).weight*span/panels*width*2.0*t;
			}
		return area;
	}

	// ---- the mirror ------------------------------------------------------

	namespace mirror = meqtest::mirror;

	/// The mirror's packs as the material table's coil rows.
	MaterialTable mirrorTable()
	{
		MaterialTable table;
		table.addCoil( mirror::pack( -1.0 ), "lower" );
		table.addCoil( mirror::pack( +1.0 ), "upper" );
		return table;
	}

	/// One of the mirror's meshes with a box-supported source on it, labelled
	/// at MR-4's two points on the VACUUM field. Under @a split the field is
	/// geq's filaments subtracted, and the potential -- the remainder -- is
	/// zero, which is the vacuum's remainder exactly; otherwise the same
	/// field is projected onto the potential, as a meshed solve holds it.
	struct MirrorCase
	{
		std::unique_ptr< mirror::Domain > domain;
		std::shared_ptr< meq::NormalisedMHDSource > plasma;
		std::unique_ptr< meq::NormalisedSource > wrapper;
		std::unique_ptr< meq::ConductorField > filaments;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr< GradShafranovSolver > solver;
		double axisLabel = 0.0;
		double boundaryLabel = 0.0;

		meq::NormalisedSource &source()
		{
			return wrapper ? *wrapper : *plasma;
		}

		MirrorCase( std::string const &stem, mirror::PlasmaMesh const &shape,
		            bool split )
			: domain( mirror::buildDomain( 0, stem, shape ) ),
			  plasma( std::make_shared< meq::NormalisedMHDSource >(
				  std::make_shared< meq::ConstantProfile const >( 1.0 ),
				  std::make_shared< meq::ConstantProfile const >( 0.0 ),
				  1.0e-3, mirror::machine().mu0 ) ),
			  filaments( std::make_unique< meq::ConductorField >(
				  mirror::geqFilaments() ) )
		{
			axisLabel = filaments->psi( 0.1, 0.0 );
			boundaryLabel = filaments->psi( 0.35, 0.0 );

			if ( !split )
				wrapper = std::make_unique< meq::CoilAugmentedNormalisedSource >(
					plasma, mirror::machineCoils() );

			solver = std::make_unique< GradShafranovSolver >( *domain->sub, 2 );
			solver->setSource( source(), axisLabel );
			solver->setBoundaryData( zero );
			if ( split )
				solver->setConductorField( *filaments );
			solver->prepare();

			if ( split )
				solver->potential() = 0.0;
			else
			{
				meq::ConductorField const *field = filaments.get();
				mfem::FunctionCoefficient total( [ field ]( mfem::Vector const &x )
				{
					return field->psi( x( 0 ), x( 1 ) );
				} );
				solver->potential().ProjectCoefficient( total );
			}

			source().setNormalisation( axisLabel, boundaryLabel );
			source().setPlasmaSupportShape( meq::PlasmaSupportShape::Box, 1.5 );
		}
	};

	/// The background's element attributes through the SubMesh's parent map,
	/// since the cut overwrote them with its own marker.
	std::vector< int > parentAttributes( mirror::Domain const &d,
	                                     std::string const &path )
	{
		mfem::Mesh original( path.c_str(), 1, 1 );
		mfem::Array< int > const &parent = d.sub->GetParentElementIDMap();
		std::vector< int > out( static_cast< std::size_t >( parent.Size() ) );
		for ( int e = 0; e < parent.Size(); ++e )
			out[ static_cast< std::size_t >( e ) ] = original.GetAttribute( parent[ e ] );
		return out;
	}
}

/*
 * THE PARTITION, AND THE PLASMA AREA AGAINST THE CLOSED FORM.
 *
 * The box stops at z = -0.70, just above the X-point at -0.704, so { psi < 0 }
 * inside it is the core alone and its area is a root pair per height. Each
 * mesh holds the closed form projected at k = 2, and the plasma area is
 * sum( plasma_fraction * |e| ). The error has two parts: psi_h's own, which
 * moves the edge by O( h^3 ), and the source rule integrating an indicator
 * on the cut elements, O( h^2 ) each over O( 1/h ) of them -- so the rate is
 * at least one, and the count of cut elements doubles per halving.
 */
BOOST_AUTO_TEST_CASE( thePlasmaAreaConvergesOnTheExactSoloviev )
{
	double const r0 = 0.6, r1 = 1.4, z0 = -0.70, z1 = 0.75;

	double const reference = exactArea( z0, 64, 20 );
	double const check = exactArea( z0, 32, 20 );
	std::printf( "\n  ER-0: THE PLASMA AREA ON ITER EXAMPLE 2, box R [%.2f, %.2f] "
	             "z [%.2f, %.2f], k = 2\n", r0, r1, z0, z1 );
	std::printf( "    exact area %.15f (64 panels), %.3e from 32 panels\n",
	             reference, std::abs( reference - check ) );
	BOOST_TEST_REQUIRE( std::abs( reference - check ) <= 1.0e-12*reference,
	                    "the reference area is not converged: " << reference
	                    << " against " << check );

	std::printf( "    %5s %8s %18s %12s %6s %10s %8s %12s\n", "n", "elements",
	             "plasma area", "error", "rate", "boundary", "ratio", "partition" );

	std::vector< double > errors;
	std::vector< int > boundaries;
	for ( int n : { 8, 16, 32, 64, 128, 256 } )
	{
		SolovievBox box( n, 2*n, r0, r1, z0, z1, 2 );
		ElementMaterials const record = box.solver->elementMaterials( MaterialTable() );

		double plasmaArea = 0.0;
		for ( int e = 0; e < record.elementCount(); ++e )
			plasmaArea += record.plasmaFraction[ static_cast< std::size_t >( e ) ]
			              *record.area[ static_cast< std::size_t >( e ) ];
		double const error = std::abs( plasmaArea - reference );
		int const boundary = record.regionCount( ElementRegion::PlasmaBoundary );
		double const partition = worstPartition( record );

		std::printf( "    %5d %8d %18.15f %12.4e", n, record.elementCount(),
		             plasmaArea, error );
		if ( !errors.empty() )
			std::printf( " %6.2f", std::log2( errors.back()/error ) );
		else
			std::printf( " %6s", "-" );
		std::printf( " %10d", boundary );
		if ( !boundaries.empty() )
			std::printf( " %8.3f", static_cast< double >( boundary )/boundaries.back() );
		else
			std::printf( " %8s", "-" );
		std::printf( " %12.3e\n", partition );
		std::fflush( stdout );

		BOOST_TEST( partition <= 1.0e-15, "n = " << n << ": the fractions sum to "
		            "1 within " << partition );
		// RELATIVE, and not round-off of one sum: MakeCartesian2D's vertices
		// are i*h + r0, rounded, so the elements tile a box within an ulp of
		// the one named.
		BOOST_TEST( std::abs( totalArea( record ) - ( r1 - r0 )*( z1 - z0 ) )
		            <= 1.0e-11*( r1 - r0 )*( z1 - z0 ), "the element areas are the box's: "
		            << totalArea( record ) - ( r1 - r0 )*( z1 - z0 ) );
		BOOST_TEST( record.faults == 0 );
		BOOST_TEST( record.regionCount( ElementRegion::Mixed ) == 0 );
		BOOST_TEST( record.regionCount( ElementRegion::Plasma ) > 0 );
		BOOST_TEST( std::abs( record.materialArea[ MaterialTable::plasma ] - plasmaArea )
		            <= 1.0e-12*plasmaArea );

		errors.push_back( error );
		boundaries.push_back( boundary );
	}

	/*
	 * THE RATE IS A FIT AND AN ENVELOPE, NOT A STEP-BY-STEP ORDER. The cut
	 * elements' indicator error is O( h^2 ) each with a sign the mesh decides,
	 * so over O( 1/h ) of them the sum is O( h ) at worst and O( h^1.5 ) when
	 * the signs are uncorrelated -- and single steps scatter about that. What
	 * is asserted is the worst case: a least-squares slope of at least one
	 * over the whole ladder, and every level under the coarsest one's error
	 * carried down at first order.
	 */
	double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
	for ( std::size_t i = 0; i < errors.size(); ++i )
	{
		double const x = static_cast< double >( i );
		double const y = -std::log2( errors[ i ] );
		sx += x;
		sy += y;
		sxx += x*x;
		sxy += x*y;
	}
	double const count = static_cast< double >( errors.size() );
	double const slope = ( count*sxy - sx*sy )/( count*sxx - sx*sx );
	std::printf( "    least-squares rate %.3f over %zu levels\n", slope, errors.size() );
	BOOST_TEST( slope >= 1.0, "the plasma area's error falls at a fitted rate "
	            << slope );
	for ( std::size_t i = 1; i < errors.size(); ++i )
	{
		BOOST_TEST( errors[ i ] <= errors[ 0 ]/std::pow( 2.0, static_cast< double >( i ) ),
		            "level " << i << " is above the first-order envelope" );
		double const ratio = static_cast< double >( boundaries[ i ] )/boundaries[ i - 1 ];
		BOOST_TEST( ( ratio > 1.8 && ratio < 2.2 ), "the plasma_boundary count "
		            "grew by " << ratio << " under a halving of h, against 2 "
		            "for a curve's length over h" );
	}
}

/*
 * THE FILL AND THE EXCLUSIONS ARE THE INTEGRATOR'S OWN.
 *
 * The box reaches z = -0.9, below the X-point, so { psi < 0 } has a second
 * lobe -- the private flux region -- whose Psi stays below 0.47 and so cannot
 * seed the fill. Under PlasmaConnectivity::Component it carries no source and
 * must read plasma_fraction 0; under Pointwise it does and must read 1. An
 * attribute named by setPlasmaExclusion(), a patch on the axis, reads 0 under
 * both.
 */
BOOST_AUTO_TEST_CASE( theFillAndTheExclusionsAreTheIntegratorsOwn )
{
	double const r0 = 0.55, r1 = 1.45, z0 = -0.9, z1 = 0.75;
	std::array< double, 4 > const patch{ 1.0, 1.1, -0.05, 0.1 };

	auto deepInTheLobe = [ ]( mfem::Mesh const &mesh, int e )
	{
		return extentOf( mesh, e ).zMax < -0.78;
	};

	struct Reading
	{
		int lobe = 0;
		int lobeCarrying = 0;
		int lobeWhole = 0;
		int patch = 0;
		int patchCarrying = 0;
		int patchFlagged = 0;
		int core = 0;
		int coreWhole = 0;
	};

	auto read = [ & ]( GradShafranovSolver::PlasmaConnectivity choice )
	{
		SolovievBox box( 30, 55, r0, r1, z0, z1, 2, { patch } );
		box.solver->setPlasmaConnectivity( choice );
		box.solver->refreshPlasmaComponent( box.solver->potential() );
		ElementMaterials const record = box.solver->elementMaterials( MaterialTable() );
		BOOST_TEST( worstPartition( record ) <= 1.0e-15 );

		Reading r;
		for ( int e = 0; e < record.elementCount(); ++e )
		{
			double const f = record.plasmaFraction[ static_cast< std::size_t >( e ) ];
			if ( deepInTheLobe( box.mesh, e ) && iter().psi(
				0.5*( extentOf( box.mesh, e ).rMin + extentOf( box.mesh, e ).rMax ),
				0.5*( extentOf( box.mesh, e ).zMin + extentOf( box.mesh, e ).zMax ) ) < 0.0 )
			{
				++r.lobe;
				r.lobeCarrying += f > 0.0 ? 1 : 0;
				r.lobeWhole += f == 1.0 ? 1 : 0;
			}
			if ( box.mesh.GetAttribute( e ) == 2 )
			{
				++r.patch;
				r.patchCarrying += f > 0.0 ? 1 : 0;
				r.patchFlagged += record.excluded[ static_cast< std::size_t >( e ) ];
			}
			else
			{
				mfem::Vector centre;
				box.mesh.GetElementCenter( e, centre );
				if ( std::abs( centre( 0 ) - 1.05 ) < 0.15 && std::abs( centre( 1 ) ) < 0.3 )
				{
					++r.core;
					r.coreWhole += f == 1.0 ? 1 : 0;
				}
			}
		}
		return r;
	};

	Reading const component = read( GradShafranovSolver::PlasmaConnectivity::Component );
	Reading const pointwise = read( GradShafranovSolver::PlasmaConnectivity::Pointwise );

	std::printf( "\n  ER-0: THE FILL AND THE EXCLUSION, ITER Example 2 below its "
	             "X-point\n" );
	std::printf( "    %-10s %6s %10s %8s %6s %10s %8s %6s %8s\n", "fill", "lobe",
	             "carrying", "whole", "patch", "carrying", "flagged", "core",
	             "whole" );
	for ( auto const &[ name, r ] : { std::make_pair( "component", component ),
	                                  std::make_pair( "pointwise", pointwise ) } )
		std::printf( "    %-10s %6d %10d %8d %6d %10d %8d %6d %8d\n", name, r.lobe,
		             r.lobeCarrying, r.lobeWhole, r.patch, r.patchCarrying,
		             r.patchFlagged, r.core, r.coreWhole );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( component.lobe > 10, "the box must reach the private "
	                    "flux region for this case to ask anything" );
	BOOST_TEST( component.lobeCarrying == 0, component.lobeCarrying
	            << " elements of the disconnected lobe read plasma under the fill" );
	BOOST_TEST( pointwise.lobeWhole > component.lobe/2, "THE CONTROL: without "
	            "the fill the lobe must read plasma, or the zero above is not the "
	            "fill's" );

	BOOST_TEST_REQUIRE( component.patch > 0 );
	for ( Reading const &r : { component, pointwise } )
	{
		BOOST_TEST( r.patchCarrying == 0, "an excluded element reads plasma" );
		BOOST_TEST( r.patchFlagged == r.patch, "an excluded element is not "
		            "flagged excluded" );
		BOOST_TEST( r.coreWhole == r.core, "the core beside the patch must still "
		            "be plasma, or the patch's zero is not the exclusion's" );
	}
}

/*
 * A PLASMA THAT REACHES A COIL IS A FAULT, POSED DELIBERATELY: a coil row
 * whose rectangle sits on the magnetic axis. Every point in it is plasma and
 * copper at once, so every element it touches is faulted and the copper keeps
 * the points -- its area is what it is with the plasma labelled away. The control
 * moves the same rectangle into the vacuum corner of the box: no fault, and
 * cut copper cells, since this mesh is not aligned to it.
 */
BOOST_AUTO_TEST_CASE( aPlasmaInsideACoilIsAFault )
{
	SolovievBox box( 16, 32, 0.6, 1.4, -0.70, 0.75, 2 );

	MaterialTable inside;
	int const insideRow = inside.addCoil( meq::Coil( 1.0513, 0.0237, 0.0517, 0.0411, 1.0e5 ),
	                                      "on-axis" );
	MaterialTable outside;
	int const outsideRow = outside.addCoil( meq::Coil( 1.3313, 0.6637, 0.0517, 0.0411, 1.0e5 ),
	                                        "corner" );

	ElementMaterials const posed = box.solver->elementMaterials( inside );
	ElementMaterials const control = box.solver->elementMaterials( outside );

	std::printf( "\n  ER-0: A FAULT, POSED\n" );
	printRegions( "coil on the axis", posed );
	printRegions( "coil in the vacuum corner", control );

	// The same rectangle with the plasma labelled away -- no point of the box
	// is plasma at ( psi_ax, psi_bnd ) = ( -2, -1 ) -- is the copper's area
	// with nothing to fault against.
	box.source.setNormalisation( -2.0, -1.0 );
	ElementMaterials const alone = box.solver->elementMaterials( inside );
	BOOST_TEST_REQUIRE( alone.regionCount( ElementRegion::Plasma ) == 0 );

	BOOST_TEST( posed.faults > 0 );
	BOOST_TEST( alone.faults == 0 );
	BOOST_TEST( posed.materialArea[ static_cast< std::size_t >( insideRow ) ]
	            == alone.materialArea[ static_cast< std::size_t >( insideRow ) ],
	            "the copper keeps its points under a fault, so its area cannot "
	            "depend on the plasma" );
	int faultedOutsideTheCoil = 0;
	for ( int e = 0; e < posed.elementCount(); ++e )
		if ( posed.faulted[ static_cast< std::size_t >( e ) ]
		     && posed.fraction( e, insideRow ) == 0.0 )
			++faultedOutsideTheCoil;
	BOOST_TEST( faultedOutsideTheCoil == 0 );
	BOOST_TEST( worstPartition( posed ) <= 1.0e-15 );

	BOOST_TEST( control.faults == 0 );
	BOOST_TEST( control.regionCount( ElementRegion::ConductorBoundary ) > 0 );
	BOOST_TEST( control.materialArea[ static_cast< std::size_t >( outsideRow ) ] > 0.0 );
	BOOST_TEST( control.regionCount( ElementRegion::Mixed ) == 0 );
}

/*
 * THE CONDUCTORS AND THE CUTOFF, ON geq'S MIRROR.
 *
 * Four meshes of one machine, level 0, cut at Gamma:
 *
 *   meshed      the vacuum's mesh, the packs a --coil each. No cell may be
 *               conductor_boundary, each pack's area is its rectangle's to
 *               round-off, and its cells are exactly attribute 10 + i.
 *   unaligned   M-262's unaligned plasma mesh: no --coil, as a subtracting
 *               model meshes it, and the box's edges cut elements. The
 *               copper is cut on every edge of every pack, and the cutoff
 *               |z| = 1.5 runs through plasma_boundary cells.
 *   aligned     M-262's aligned mesh: no element straddles the cutoff over
 *               the box's span R < 0.75, where it is a line of the geometry.
 *   conductors  the aligned mesh with AlignConductors' --align-coil: the
 *               copper pure again, at an element count printed beside the
 *               aligned mesh's, which is that key's cost.
 *
 * The three plasma meshes carry geq's filaments subtracted, as the shipped
 * mirror does; the meshed one projects the same field. No mesh may fault.
 */
BOOST_AUTO_TEST_CASE( theConductorsAndTheCutoffOnTheMirror )
{
	MaterialTable const table = mirrorTable();
	double const edge = 1.5;

	mirror::PlasmaMesh aligned;
	aligned.axialEdge = edge;
	aligned.aligned = true;
	mirror::PlasmaMesh unaligned = aligned;
	unaligned.aligned = false;
	mirror::PlasmaMesh conductors = aligned;
	conductors.alignConductors = true;

	MirrorCase meshed( "er-meshed", mirror::PlasmaMesh(), false );
	MirrorCase cut( "er-unaligned", unaligned, true );
	MirrorCase fitted( "er-aligned", aligned, true );
	MirrorCase pure( "er-align-conductors", conductors, true );

	ElementMaterials const meshedRecord = meshed.solver->elementMaterials( table );
	ElementMaterials const cutRecord = cut.solver->elementMaterials( table );
	ElementMaterials const fittedRecord = fitted.solver->elementMaterials( table );
	ElementMaterials const pureRecord = pure.solver->elementMaterials( table );

	std::printf( "\n  ER-0: geq'S MIRROR, labels %.6e at ( 0.1, 0 ) and %.6e at "
	             "( 0.35, 0 ), box support |z| <= %.1f\n", meshed.axisLabel,
	             meshed.boundaryLabel, edge );
	printRegions( "meshed packs (--coil)", meshedRecord );
	printRegions( "M-262 unaligned, no --coil", cutRecord );
	printRegions( "M-262 aligned, no --coil", fittedRecord );
	printRegions( "aligned + AlignConductors", pureRecord );

	for ( ElementMaterials const *r : { &meshedRecord, &cutRecord, &fittedRecord,
	                                    &pureRecord } )
	{
		BOOST_TEST( worstPartition( *r ) <= 1.0e-15 );
		BOOST_TEST( r->faults == 0 );
		BOOST_TEST( r->regionCount( ElementRegion::Plasma ) > 0 );
		BOOST_TEST( r->regionCount( ElementRegion::Mixed ) == 0 );
	}

	// MESHED: pure, the rectangle's area, and the generator's own attribute.
	std::vector< int > const attributes = parentAttributes(
		*meshed.domain, mirror::meshFile( "er-meshed" ) );
	BOOST_TEST( meshedRecord.regionCount( ElementRegion::ConductorBoundary ) == 0 );
	for ( std::size_t i = 0; i < table.coilCount(); ++i )
	{
		int const row = table.coilMaterial( i );
		double const rectangle = table.coil( i ).area();
		std::printf( "    coil %zu: area %.15e meshed, %.15e aligned, rectangle "
		             "%.15e; %d of 4 edges cut unaligned\n", i,
		             meshedRecord.materialArea[ static_cast< std::size_t >( row ) ],
		             pureRecord.materialArea[ static_cast< std::size_t >( row ) ],
		             rectangle,
		             edgesCut( *cut.domain->sub, cutRecord, row, table.coil( i ) ) );
		BOOST_TEST( std::abs( meshedRecord.materialArea[ static_cast< std::size_t >( row ) ]
		                      - rectangle ) <= 1.0e-12*rectangle );
		BOOST_TEST( std::abs( pureRecord.materialArea[ static_cast< std::size_t >( row ) ]
		                      - rectangle ) <= 1.0e-12*rectangle );

		int disagree = 0;
		for ( int e = 0; e < meshedRecord.elementCount(); ++e )
		{
			bool const copper = meshedRecord.fraction( e, row ) == 1.0;
			bool const tagged = attributes[ static_cast< std::size_t >( e ) ]
			                    == 10 + static_cast< int >( i );
			disagree += copper != tagged ? 1 : 0;
		}
		BOOST_TEST( disagree == 0, disagree << " elements disagree with "
		            "halfdisc.py's attribute " << 10 + i << " about pack " << i );

		// THE CONTROL: unaligned, every edge of every pack cuts cells.
		BOOST_TEST( edgesCut( *cut.domain->sub, cutRecord, row, table.coil( i ) ) == 4 );
	}

	// ALIGNCONDUCTORS: pure, and its cost against the same mesh without it.
	BOOST_TEST( cutRecord.regionCount( ElementRegion::ConductorBoundary ) > 0 );
	BOOST_TEST( fittedRecord.regionCount( ElementRegion::ConductorBoundary ) > 0 );
	BOOST_TEST( pureRecord.regionCount( ElementRegion::ConductorBoundary ) == 0 );
	std::printf( "    AlignConductors: %d elements inside Gamma against %d without "
	             "it (+%d)\n", pureRecord.elementCount(), fittedRecord.elementCount(),
	             pureRecord.elementCount() - fittedRecord.elementCount() );
	BOOST_TEST( pureRecord.elementCount() >= fittedRecord.elementCount() );

	// THE CUTOFF.
	auto straddling = [ & ]( MirrorCase const &c, ElementMaterials const &r,
	                         bool boundaryOnly )
	{
		int count = 0;
		mfem::Mesh const &mesh = *c.domain->sub;
		for ( int e = 0; e < r.elementCount(); ++e )
		{
			if ( boundaryOnly && r.region[ static_cast< std::size_t >( e ) ]
			                     != ElementRegion::PlasmaBoundary )
				continue;
			// Over the box's own span in R, R < 0.75, which is where the
			// cutoff is a line of the geometry; beyond it the line is nothing.
			if ( extentOf( mesh, e ).rMin >= 0.75 )
				continue;
			if ( straddlesZ( mesh, e, edge ) || straddlesZ( mesh, e, -edge ) )
				++count;
		}
		return count;
	};
	int const fittedStraddle = straddling( fitted, fittedRecord, false );
	int const cutStraddle = straddling( cut, cutRecord, true );
	std::printf( "    the cutoff |z| = %.1f: %d elements straddle it aligned, %d "
	             "plasma_boundary elements straddle it unaligned\n", edge,
	             fittedStraddle, cutStraddle );
	BOOST_TEST( fittedStraddle == 0 );
	BOOST_TEST( straddling( pure, pureRecord, false ) == 0 );
	BOOST_TEST( cutStraddle > 0, "THE CONTROL: M-262's unaligned mesh must put "
	            "plasma_boundary cells on the cutoff" );
}

/*
 * THE PREDICATE IS THE INTEGRATED ONE, UNDER A SPLIT.
 *
 * `F/R` summed over the points elementMaterials() called plasma, at the flux
 * it recorded for them, must reproduce assemblePlasmaCurrent() -- the plasma
 * load the bordered Newton constrains. Only the summation order may differ.
 * The field is geq's filaments subtracted with a zero remainder, so a
 * classifier that judged psi_h without psi_c would find no plasma at all.
 * A Solov'ev box, no split, is the second arm.
 *
 * The sum alone sees only one direction: the plasma term is zero outside the
 * support, so a classifier calling too MANY points plasma leaves it unchanged.
 * The Solov'ev arm's profiles are constant, so its plasma term is non-zero at
 * every point the support admits, and a point called plasma where it is zero
 * is a point the plasma current was not integrated over.
 */
BOOST_AUTO_TEST_CASE( thePredicateIsTheIntegratedOne )
{
	mirror::PlasmaMesh aligned;
	aligned.axialEdge = 1.5;
	aligned.aligned = true;
	MirrorCase m( "er-predicate", aligned, true );
	SolovievBox box( 16, 32, 0.6, 1.4, -0.70, 0.75, 2 );

	struct Arm
	{
		char const *name;
		GradShafranovSolver *solver;
		meq::NormalisedSource *source;
	};

	std::printf( "\n  ER-0: int F/R OVER THE POINTS CALLED PLASMA\n" );
	for ( Arm const &arm : { Arm{ "mirror, filaments subtracted", m.solver.get(),
	                              &m.source() },
	                         Arm{ "Solov'ev box", box.solver.get(), &box.source } } )
	{
		std::vector< GradShafranovSolver::ClassifiedPoint > points;
		ElementMaterials const record =
			arm.solver->elementMaterials( MaterialTable(), &points );
		double const classified = plasmaLoadOverClassifiedPoints(
			*arm.source, points, record.elementCount() );
		double const assembled = arm.solver->assemblePlasmaCurrent( stateOf( *arm.solver ) );
		double const relative = std::abs( classified - assembled )/std::abs( assembled );

		std::printf( "    %-30s classified %.15e assembled %.15e relative %.3e\n",
		             arm.name, classified, assembled, relative );
		std::fflush( stdout );

		int unloaded = 0;
		for ( GradShafranovSolver::ClassifiedPoint const &p : points )
			if ( p.material == MaterialTable::plasma
			     && arm.source->scaledF( p.radius, p.z, p.psi ) == 0.0 )
				++unloaded;
		std::printf( "    %-30s points called plasma with no plasma term: %d\n",
		             arm.name, unloaded );
		std::fflush( stdout );

		BOOST_TEST_REQUIRE( assembled != 0.0 );
		BOOST_TEST( relative <= 1.0e-13, arm.name << ": the points called plasma "
		            "are not the points the plasma current was integrated over" );
		if ( arm.source == &box.source )
			BOOST_TEST( unloaded == 0, arm.name << ": " << unloaded << " points are "
			            "called plasma where the plasma term is zero, so the classifier "
			            "admits points the support does not" );
	}
}

#ifdef MEQ_DRIVER_PATH
/*
 * WHAT AlignConductors COSTS, THROUGH THE DRIVER'S OWN COMMAND.
 *
 * The shipped mirror and limited-tokamak-filament, each written out with its
 * File redirected into the build tree and the key off and on, meshed by the
 * command `meq --mesh-command` prints -- so the count is of the key as a user
 * reaches it. On limited-tokamak-filament the copper is then classified on
 * both meshes: cut cells without the key, none with it.
 */
BOOST_AUTO_TEST_CASE( alignConductorsCost )
{
	struct Example
	{
		char const *name;
		double gamma;
	};

	auto run = [ ]( std::string const &command )
	{
		return std::system( ( command + " > /dev/null 2>&1" ).c_str() );
	};

	std::printf( "\n  ER-0: AlignConductors' COST, elements (background, and "
	             "centre inside Gamma)\n" );
	std::printf( "    %-26s %10s %10s %10s %10s\n", "example", "off", "on",
	             "off in G", "on in G" );

	for ( Example const &ex : { Example{ "mirror-rotating", 3.0 },
	                            Example{ "limited-tokamak-filament", 2.40 } } )
	{
		std::array< int, 2 > background{ 0, 0 };
		std::array< int, 2 > inside{ 0, 0 };
		std::array< ElementMaterials, 2 > records;
		MaterialTable table;

		for ( int on = 0; on < 2; ++on )
		{
			std::string const stem = std::string( MEQ_TEST_OUTPUT_DIR ) + "/er-cost-"
			                         + ex.name + ( on ? "-on" : "-off" );
			std::ifstream in( std::string( "examples/" ) + ex.name + ".toml" );
			BOOST_TEST_REQUIRE( in.good() );
			std::stringstream text;
			text << in.rdbuf();
			std::string config = text.str();
			std::string const file = std::string( "File = \"examples/" ) + ex.name + ".msh\"";
			std::size_t const at = config.find( file );
			BOOST_TEST_REQUIRE( at != std::string::npos );
			config.replace( at, file.size(), "File = \"" + stem + ".msh\"" );
			if ( on )
			{
				std::size_t const block = config.find( "[mesh.generate]\n" );
				BOOST_TEST_REQUIRE( block != std::string::npos );
				config.insert( block + 16, "AlignConductors = true\n" );
			}
			{
				std::ofstream out( stem + ".toml" );
				out << config;
			}

			std::string const commandFile = stem + ".command";
			BOOST_TEST_REQUIRE( std::system( ( std::string( MEQ_DRIVER_PATH )
				+ " --mesh-command " + stem + ".toml > " + commandFile
				+ " 2> /dev/null" ).c_str() ) == 0 );
			std::ifstream commandIn( commandFile );
			std::string command;
			std::getline( commandIn, command );
			BOOST_TEST_REQUIRE( command.rfind( "halfdisc ", 0 ) == 0, command );
			BOOST_TEST( ( command.find( "--align-coil" ) != std::string::npos ) == ( on == 1 ),
			            command );
			BOOST_TEST( command.find( "--coil " ) == std::string::npos, command );
			BOOST_TEST_REQUIRE( run( std::string( "python3 " ) + MEQ_HALFDISC_PATH + " "
			                         + command.substr( 9 ) ) == 0, command );

			mfem::Mesh mesh( ( stem + ".msh" ).c_str(), 1, 1 );
			background[ on ] = mesh.GetNE();
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				mfem::Vector centre;
				mesh.GetElementCenter( e, centre );
				inside[ on ] += std::hypot( centre( 0 ), centre( 1 ) ) < ex.gamma ? 1 : 0;
			}

			// The copper alone: a confined source labelled so that no point
			// is plasma, which leaves conductor membership the only question.
			if ( std::string( ex.name ) == "limited-tokamak-filament" )
			{
				if ( table.coilCount() == 0 )
				{
					table.addCoil( meq::Coil( 1.75, 0.90, 0.05, 0.05, 1.0 ), "P1U" );
					table.addCoil( meq::Coil( 1.75, -0.90, 0.05, 0.05, 1.0 ), "P1L" );
					table.addCoil( meq::Coil( 0.55, 1.10, 0.05, 0.05, 1.0 ), "P2U" );
					table.addCoil( meq::Coil( 0.55, -1.10, 0.05, 0.05, 1.0 ), "P2L" );
				}
				meq::NormalisedMHDSource source(
					std::make_shared< meq::ConstantProfile const >( 1.0 ),
					std::make_shared< meq::ConstantProfile const >( 0.0 ), 1.0, 1.0 );
				mfem::ConstantCoefficient zero( 0.0 );
				GradShafranovSolver solver( mesh, 1 );
				solver.setSource( source, 1.0 );
				solver.setBoundaryData( zero );
				solver.prepare();
				solver.potential() = 0.0;
				source.setNormalisation( 1.0, 0.0 );
				source.setPlasmaSupport( true );
				records[ static_cast< std::size_t >( on ) ] = solver.elementMaterials( table );
			}
		}

		std::printf( "    %-26s %10d %10d %10d %10d   +%.1f%% in Gamma\n", ex.name,
		             background[ 0 ], background[ 1 ], inside[ 0 ], inside[ 1 ],
		             100.0*( inside[ 1 ] - inside[ 0 ] )/inside[ 0 ] );
		std::fflush( stdout );
		BOOST_TEST( background[ 1 ] >= background[ 0 ] );

		if ( table.coilCount() > 0 )
		{
			printRegions( "  without AlignConductors", records[ 0 ] );
			printRegions( "  with AlignConductors", records[ 1 ] );
			BOOST_TEST( records[ 0 ].regionCount( ElementRegion::ConductorBoundary ) > 0 );
			BOOST_TEST( records[ 1 ].regionCount( ElementRegion::ConductorBoundary ) == 0 );
			BOOST_TEST( records[ 1 ].regionCount( ElementRegion::Plasma ) == 0 );
			for ( std::size_t i = 0; i < table.coilCount(); ++i )
			{
				double const rectangle = table.coil( i ).area();
				double const area = records[ 1 ].materialArea[
					static_cast< std::size_t >( table.coilMaterial( i ) ) ];
				BOOST_TEST( std::abs( area - rectangle ) <= 1.0e-12*rectangle,
				            "coil " << i << ": " << area << " against " << rectangle );
			}
		}
	}
}
#endif
