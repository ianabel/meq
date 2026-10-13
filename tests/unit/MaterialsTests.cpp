// Unit tests for the material table and the region rule in src/meq/Materials.
//
// ELEMENT-REGIONS-PLAN.md §2.1 and §2.1b. Which material a quadrature point
// belongs to is asked by meq::GradShafranovSolver::elementMaterials() and is
// tested on meshes in tests/convergence/ElementRegions.cpp; what is tested here
// is everything after that question -- the table's rows, the fractions an
// element's weights make, and the region those fractions derive -- which is
// MFEM-free so that CI runs it.

#define BOOST_TEST_MODULE MaterialsTests
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "meq/Config.hpp"
#include "meq/Materials.hpp"

using meq::ElementMaterials;
using meq::ElementRegion;
using meq::Material;
using meq::MaterialFraction;
using meq::MaterialKind;
using meq::MaterialTable;

namespace
{
	/// Two coils and one of each reserved kind, as ER-1's reader must hold.
	MaterialTable machine()
	{
		MaterialTable table;
		table.addCoil( meq::Coil( 1.0, 0.5, 0.1, 0.2, 1.0e5 ), "upper" );
		table.addCoil( meq::Coil( 1.0, -0.5, 0.1, 0.2, -2.0e5 ), "lower" );
		Material wall;
		wall.kind = MaterialKind::Passive;
		wall.name = "wall";
		wall.resistivity = 7.0e-7;
		table.addMaterial( wall );
		Material iron;
		iron.kind = MaterialKind::Ferromagnetic;
		iron.name = "core";
		iron.relativePermeability = 500.0;
		table.addMaterial( iron );
		return table;
	}

	std::vector< MaterialFraction > fractions( std::vector< int > const &materials )
	{
		std::vector< MaterialFraction > out;
		for ( int m : materials )
			out.push_back( MaterialFraction{ m, 1.0/static_cast< double >( materials.size() ) } );
		return out;
	}
}

BOOST_AUTO_TEST_CASE( theTableStartsWithVacuumAndPlasmaAndKeepsCoilOrder )
{
	MaterialTable const table = machine();
	BOOST_TEST( table.size() == 6 );
	BOOST_TEST( ( table.material( MaterialTable::vacuum ).kind == MaterialKind::Vacuum ) );
	BOOST_TEST( ( table.material( MaterialTable::plasma ).kind == MaterialKind::Plasma ) );
	BOOST_TEST( table.coilCount() == 2u );
	BOOST_TEST( table.coilMaterial( 0 ) == 2 );
	BOOST_TEST( table.coilMaterial( 1 ) == 3 );
	BOOST_TEST( table.material( 2 ).name == "upper" );
	BOOST_TEST( table.material( 3 ).current == -2.0e5 );

	// THE FILL VALUES, which is what lets a new kind add rows and no columns.
	BOOST_TEST( std::isnan( table.material( MaterialTable::plasma ).current ) );
	BOOST_TEST( std::isnan( table.material( 2 ).resistivity ) );
	BOOST_TEST( table.material( 2 ).relativePermeability == 1.0 );
	BOOST_TEST( table.material( 2 ).circuit == -1 );
	BOOST_TEST( table.material( 2 ).permeabilityCurve == -1 );
	BOOST_TEST( ( table.material( 4 ).kind == MaterialKind::Passive ) );
	BOOST_TEST( table.material( 5 ).relativePermeability == 500.0 );

	BOOST_CHECK_THROW( table.material( 6 ), std::out_of_range );
	BOOST_CHECK_THROW( table.material( -1 ), std::out_of_range );

	MaterialTable other;
	Material coilWithoutShape;
	coilWithoutShape.kind = MaterialKind::Coil;
	BOOST_CHECK_THROW( other.addMaterial( coilWithoutShape ), std::invalid_argument );
	Material secondPlasma;
	secondPlasma.kind = MaterialKind::Plasma;
	BOOST_CHECK_THROW( other.addMaterial( secondPlasma ), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( conductorMembershipIsTheRectangleAndOnlyTheRectangle )
{
	MaterialTable const table = machine();
	BOOST_TEST( table.conductorContaining( 1.0, 0.5 ) == 2 );
	BOOST_TEST( table.conductorContaining( 0.95, -0.4 ) == 3 );
	// Closed, as meq::Coil::contains is.
	BOOST_TEST( table.conductorContaining( 1.1, 0.7 ) == 2 );
	BOOST_TEST( table.conductorContaining( 1.1 + 1.0e-12, 0.5 ) == -1 );
	BOOST_TEST( table.conductorContaining( 1.0, 0.0 ) == -1 );

	// A reserved row has no shape and is never found.
	MaterialTable bare;
	Material wall;
	wall.kind = MaterialKind::Passive;
	bare.addMaterial( wall );
	BOOST_TEST( bare.conductorContaining( 0.0, 0.0 ) == -1 );
}

/// §2.1's table, every row of it, and §2.1b's three pairings.
BOOST_AUTO_TEST_CASE( theRegionIsDerivedByTheTable )
{
	MaterialTable const table = machine();
	auto region = [ & ]( std::vector< int > const &materials )
	{
		return meq::deriveRegion( fractions( materials ), table );
	};

	BOOST_TEST( ( region( { 0 } ) == ElementRegion::Vacuum ) );
	BOOST_TEST( ( region( { 1 } ) == ElementRegion::Plasma ) );
	BOOST_TEST( ( region( { 3 } ) == ElementRegion::Coil ) );
	BOOST_TEST( ( region( { 4 } ) == ElementRegion::Passive ) );
	BOOST_TEST( ( region( { 5 } ) == ElementRegion::Ferromagnetic ) );

	BOOST_TEST( ( region( { 0, 1 } ) == ElementRegion::PlasmaBoundary ) );
	BOOST_TEST( ( region( { 0, 2 } ) == ElementRegion::ConductorBoundary ) );
	BOOST_TEST( ( region( { 0, 4 } ) == ElementRegion::ConductorBoundary ) );
	BOOST_TEST( ( region( { 0, 2, 3 } ) == ElementRegion::ConductorBoundary ) );

	// All three, or plasma against copper: a fault or a cell to split.
	BOOST_TEST( ( region( { 0, 1, 2 } ) == ElementRegion::Mixed ) );
	BOOST_TEST( ( region( { 1, 2 } ) == ElementRegion::Mixed ) );
	// Two conductors with no vacuum between them.
	BOOST_TEST( ( region( { 2, 3 } ) == ElementRegion::Mixed ) );

	BOOST_CHECK_THROW( meq::deriveRegion( {}, table ), std::invalid_argument );

	BOOST_TEST( meq::elementRegionName( ElementRegion::PlasmaBoundary )
	            == std::string( "plasma_boundary" ) );
	BOOST_TEST( meq::materialKindName( MaterialKind::Ferromagnetic )
	            == std::string( "ferromagnetic" ) );
	BOOST_TEST( meq::isConductor( MaterialKind::Passive ) );
	BOOST_TEST( !meq::isConductor( MaterialKind::Plasma ) );
}

/// The fractions are each weight over the sum: one material reads exactly 1,
/// several sum to 1 to round-off, and repeated materials are merged.
BOOST_AUTO_TEST_CASE( weightsBecomeFractionsThatSumToOne )
{
	MaterialTable const table = machine();
	ElementMaterials record;

	record.addElement( table, { { 0, 0.3 }, { 0, 0.1 }, { 0, 0.7 } }, false, false, 1 );
	record.addElement( table, { { 1, 0.1 }, { 0, 0.2 }, { 1, 0.3 } }, false, true, 2 );
	record.addElement( table, { { 2, 0.25 }, { 1, 0.25 }, { 0, 0.5 } }, true, false, 10 );

	BOOST_TEST( record.elementCount() == 3 );
	BOOST_TEST( record.fraction( 0, 0 ) == 1.0 );
	BOOST_TEST( ( record.region[ 0 ] == ElementRegion::Vacuum ) );

	BOOST_TEST( record.offsets[ 2 ] - record.offsets[ 1 ] == 2 );
	BOOST_TEST( record.entries[ 1 ].material == 0 );
	BOOST_TEST( record.entries[ 2 ].material == 1 );
	BOOST_TEST( std::abs( record.fraction( 1, 1 ) - 2.0/3.0 ) <= 1.0e-15 );
	BOOST_TEST( record.plasmaFraction[ 1 ] == record.fraction( 1, 1 ) );
	BOOST_TEST( ( record.region[ 1 ] == ElementRegion::PlasmaBoundary ) );
	BOOST_TEST( record.excluded[ 1 ] == 1 );

	BOOST_TEST( ( record.region[ 2 ] == ElementRegion::Mixed ) );
	BOOST_TEST( record.faults == 1 );
	BOOST_TEST( record.attribute[ 2 ] == 10 );

	for ( int e = 0; e < record.elementCount(); ++e )
	{
		double sum = 0.0;
		for ( int i = record.offsets[ e ]; i < record.offsets[ e + 1 ]; ++i )
			sum += record.entries[ static_cast< std::size_t >( i ) ].fraction;
		BOOST_TEST( std::abs( sum - 1.0 ) <= 4.0e-16, "element " << e );
	}

	BOOST_TEST( std::abs( record.area[ 0 ] - 1.1 ) <= 1.0e-15 );
	BOOST_TEST( std::abs( record.materialArea[ 0 ] - ( 1.1 + 0.2 + 0.5 ) ) <= 1.0e-15 );
	BOOST_TEST( std::abs( record.materialArea[ 1 ] - 0.65 ) <= 1.0e-15 );
	BOOST_TEST( record.materialArea[ 2 ] == 0.25 );
	BOOST_TEST( record.regionCount( ElementRegion::Mixed ) == 1 );
	BOOST_TEST( record.fraction( 2, 4 ) == 0.0 );

	BOOST_CHECK_THROW( record.addElement( table, {}, false, false, 1 ), std::invalid_argument );
	BOOST_CHECK_THROW( record.addElement( table, { { 6, 1.0 } }, false, false, 1 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( record.addElement( table, { { 0, -1.0 } }, false, false, 1 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( record.addElement( table, { { 0, 0.0 } }, false, false, 1 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( record.fraction( 3, 0 ), std::out_of_range );
	BOOST_TEST( record.elementCount() == 3, "a refused element must not be half-appended" );
}

namespace
{
	meq::CoilParameters coilBlock( std::string const &name, double centreZ,
	                               double current )
	{
		meq::CoilParameters coil;
		coil.name = name;
		coil.centreR = 1.6;
		coil.centreZ = centreZ;
		coil.halfWidth = 0.1;
		coil.halfHeight = 0.05;
		coil.current = current;
		return coil;
	}

	/// An anti-series pair on one supply beside a coil on none, resolved as
	/// the parse resolves it.
	meq::CoilConfig pairAndSolo()
	{
		meq::CoilConfig config;
		config.coils = { coilBlock( "PU", 0.8, 2.0e5 ), coilBlock( "PL", -0.8, -2.0e5 ),
		                 coilBlock( "Solo", 0.0, 3.0e5 ) };
		meq::CircuitParameters pair;
		pair.name = "P";
		pair.current = 1.0e5;
		pair.members = { meq::CircuitMember{ "PU", 2.0, 0 },
		                 meq::CircuitMember{ "PL", -2.0, 1 } };
		config.coils[ 0 ].circuit = "P";
		config.coils[ 0 ].multiplier = 2.0;
		config.coils[ 1 ].circuit = "P";
		config.coils[ 1 ].multiplier = -2.0;
		config.circuits = { pair };
		return config;
	}
}

/// §2.1: every coil outside a [[circuits]] block is in ONE group named
/// `default`, a grouping and not a supply, present only when it has a member.
BOOST_AUTO_TEST_CASE( everyCoilIsInItsCircuitOrInTheDefaultGroup )
{
	MaterialTable const table = meq::buildMaterialTable( pairAndSolo() );
	BOOST_TEST_REQUIRE( table.circuits().size() == 2u );
	BOOST_TEST( table.circuits()[ 0 ].name == "P" );
	BOOST_TEST( table.circuits()[ 0 ].supplied );
	BOOST_TEST( table.circuits()[ 0 ].current == 1.0e5 );
	BOOST_TEST( table.circuits()[ 1 ].name == meq::defaultCircuitName );
	BOOST_TEST( !table.circuits()[ 1 ].supplied );
	BOOST_TEST( std::isnan( table.circuits()[ 1 ].current ) );

	BOOST_TEST_REQUIRE( table.coilCount() == 3u );
	BOOST_TEST( table.material( table.coilMaterial( 0 ) ).circuit == 0 );
	BOOST_TEST( table.material( table.coilMaterial( 1 ) ).circuit == 0 );
	BOOST_TEST( table.material( table.coilMaterial( 2 ) ).circuit == 1 );
	BOOST_TEST( table.material( table.coilMaterial( 1 ) ).name == "PL" );
	BOOST_TEST( table.material( table.coilMaterial( 1 ) ).current == -2.0e5 );

	// THE SOLVED CURRENTS, when an inverse solve moved them: each coil's own,
	// and the circuit's as its first member's over that member's multiplier.
	MaterialTable const solved =
		meq::buildMaterialTable( pairAndSolo(), { 3.0e5, -3.0e5, 1.0e5 } );
	BOOST_TEST( solved.circuits()[ 0 ].current == 1.5e5 );
	BOOST_TEST( solved.material( solved.coilMaterial( 1 ) ).current == -3.0e5 );
	BOOST_TEST( solved.material( solved.coilMaterial( 2 ) ).current == 1.0e5 );
	BOOST_CHECK_THROW( meq::buildMaterialTable( pairAndSolo(), { 1.0 } ),
	                   std::invalid_argument );

	// ALL WIRED: no empty group.
	meq::CoilConfig wired = pairAndSolo();
	wired.coils.pop_back();
	MaterialTable const allWired = meq::buildMaterialTable( wired );
	BOOST_TEST_REQUIRE( allWired.circuits().size() == 1u );
	BOOST_TEST( allWired.circuits()[ 0 ].name == "P" );

	// No coils, no circuits: a fixed-boundary run's table is vacuum and plasma.
	MaterialTable const none = meq::buildMaterialTable( meq::CoilConfig() );
	BOOST_TEST( none.size() == 2 );
	BOOST_TEST( none.circuits().empty() );

	MaterialTable named;
	meq::Circuit row;
	row.name = "P";
	named.addCircuit( row );
	BOOST_CHECK_THROW( named.addCircuit( row ), std::invalid_argument );
	BOOST_CHECK_THROW( named.addCircuit( meq::Circuit() ), std::invalid_argument );
	BOOST_CHECK_THROW( named.addCoil( meq::Coil( 1.0, 0.0, 0.1, 0.1, 1.0 ), "c", 1 ),
	                   std::invalid_argument );
}

/// §2.1: a ferromagnetic row holds a constant OR a curve, the curve names its
/// material, and neither reserved kind is on a supply.
BOOST_AUTO_TEST_CASE( aReservedRowHoldsOnePermeabilityAndNoCircuit )
{
	MaterialTable table;
	meq::PermeabilityCurve curve;
	curve.material = "yoke";
	curve.fieldStrength = { 0.0, 100.0, 1.0e4 };
	curve.fluxDensity = { 0.0, 1.0, 1.8 };
	int const index = table.addPermeabilityCurve( curve );
	BOOST_TEST( index == 0 );

	Material iron;
	iron.kind = MaterialKind::Ferromagnetic;
	iron.name = "yoke";
	iron.permeabilityCurve = index;

	// Both a constant and a curve: refused, whichever is meant.
	BOOST_CHECK_THROW( table.addMaterial( iron ), std::invalid_argument );
	iron.relativePermeability = std::nan( "" );
	BOOST_TEST( table.addMaterial( iron ) == 2 );

	// A curve that is another material's, or that is not there.
	Material other = iron;
	other.name = "pole";
	BOOST_CHECK_THROW( table.addMaterial( other ), std::invalid_argument );
	other.permeabilityCurve = 1;
	BOOST_CHECK_THROW( table.addMaterial( other ), std::invalid_argument );

	// Linear iron needs a usable constant.
	Material linear;
	linear.kind = MaterialKind::Ferromagnetic;
	linear.relativePermeability = 0.0;
	BOOST_CHECK_THROW( table.addMaterial( linear ), std::invalid_argument );

	// Neither kind is on a supply, and a passive row is not magnetic.
	Material wall;
	wall.kind = MaterialKind::Passive;
	wall.circuit = 0;
	BOOST_CHECK_THROW( table.addMaterial( wall ), std::invalid_argument );
	wall.circuit = -1;
	wall.relativePermeability = 2.0;
	BOOST_CHECK_THROW( table.addMaterial( wall ), std::invalid_argument );

	// A curve from ( 0, 0 ), increasing in H.
	meq::PermeabilityCurve bad = curve;
	bad.fieldStrength[ 0 ] = 1.0;
	BOOST_CHECK_THROW( table.addPermeabilityCurve( bad ), std::invalid_argument );
	bad = curve;
	bad.fieldStrength[ 2 ] = 50.0;
	BOOST_CHECK_THROW( table.addPermeabilityCurve( bad ), std::invalid_argument );
	bad = curve;
	bad.fluxDensity.pop_back();
	BOOST_CHECK_THROW( table.addPermeabilityCurve( bad ), std::invalid_argument );
	bad = curve;
	bad.material.clear();
	BOOST_CHECK_THROW( table.addPermeabilityCurve( bad ), std::invalid_argument );
}

/*
 * THE VTK'S CONDUCTOR COLUMN: the conductor with the largest share of an
 * element, of any conductor kind, with ties to the lower row -- and plasma
 * and vacuum are never it, however large their share.
 */
BOOST_AUTO_TEST_CASE( theDominantConductorIsTheLargestAndTiesGoLow )
{
	MaterialTable const table = machine();   // rows 2, 3 coils; 4 wall; 5 iron
	ElementMaterials record;
	record.addElement( table, { { MaterialTable::plasma, 0.7 }, { 3, 0.1 }, { 4, 0.2 } },
	                   true, false, 1 );
	record.addElement( table, { { 5, 0.25 }, { 2, 0.25 }, { MaterialTable::vacuum, 0.5 } },
	                   false, false, 1 );
	record.addElement( table, { { MaterialTable::plasma, 0.4 },
	                            { MaterialTable::vacuum, 0.6 } }, false, false, 1 );

	MaterialFraction const passive = record.dominantConductor( 0, table );
	BOOST_TEST( passive.material == 4 );
	BOOST_TEST( passive.fraction == 0.2 );
	BOOST_TEST( record.dominantConductor( 1, table ).material == 2 );
	BOOST_TEST( record.dominantConductor( 1, table ).fraction == 0.25 );
	BOOST_TEST( record.dominantConductor( 2, table ).material == -1 );
	BOOST_TEST( record.dominantConductor( 2, table ).fraction == 0.0 );
	BOOST_CHECK_THROW( record.dominantConductor( 3, table ), std::out_of_range );
}
