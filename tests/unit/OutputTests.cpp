// Unit tests for the material record in src/meq/Output: NetCDFWriter::materials(),
// NetCDFWriter::nodeMaterial(), and the VTK cell data of meq::ElementCellData.
//
// ELEMENT-REGIONS-PLAN.md §2.3 and ER-1. The reserved kinds -- passive and
// ferromagnetic -- must be READABLE BEFORE THEY ARE PRODUCED, so this writes a
// table holding one passive row and two ferromagnetic rows, one linear and one
// with a permeability curve, and reads it back. The read is through ncdump and
// not through MEQ's writer, so the check shares no code with the thing it
// checks; the driver's own files are read back in DriverAcceptance.
//
// The grid the writer is constructed on is a formality here: NetCDFWriter
// always carries one, and a two-by-two box is the cheapest that locates.

#define BOOST_TEST_MODULE OutputTests
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <unistd.h>
#include <vector>

#include "mfem.hpp"

#include "meq/Materials.hpp"
#include "meq/Output.hpp"
#include "meq/Sampler.hpp"

#include "../convergence/VtuReader.hpp"

using meq::ElementMaterials;
using meq::Material;
using meq::MaterialKind;
using meq::MaterialTable;

namespace
{
	/// A path in the system's temporary directory, removed when this goes out
	/// of scope, so a red run leaves nothing behind in the build tree.
	struct Scratch
	{
		std::string path;
		explicit Scratch( std::string const &name )
			: path( ( std::filesystem::temp_directory_path()
			          / ( name + "-" + std::to_string( ::getpid() ) ) ).string() ) { }
		~Scratch() { std::remove( path.c_str() ); }
	};

	std::string slurp( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( ( std::istreambuf_iterator<char>( stream ) ),
		                    std::istreambuf_iterator<char>() );
	}

	/// ncdump's whole output for the file at 17 significant digits, which is
	/// what makes a double read back bit for bit.
	std::string dump( std::string const &file )
	{
		Scratch const text( "meq-output-tests-dump.txt" );
		std::string const command = "ncdump -p 9,17 " + file + " > " + text.path + " 2>&1";
		if ( std::system( command.c_str() ) != 0 )
			return std::string();
		return slurp( text.path );
	}

	/// One variable's values as the tokens ncdump prints between `name = ` and
	/// `;` after `data:`: a quoted string loses its quotes, and a fill is `_`.
	std::vector<std::string> tokens( std::string const &text, std::string const &name )
	{
		std::vector<std::string> out;
		std::size_t at = text.find( "\ndata:" );
		if ( at == std::string::npos )
			return out;
		// " =" and not " = ": a two-dimensional variable's values start on
		// the next line, and the loop below skips the whitespace either way.
		at = text.find( "\n " + name + " =", at );
		if ( at == std::string::npos )
			return out;
		at += name.size() + 4;
		while ( at < text.size() && text[ at ] != ';' )
		{
			char const c = text[ at ];
			if ( c == '"' )
			{
				std::size_t const end = text.find( '"', at + 1 );
				out.push_back( text.substr( at + 1, end - at - 1 ) );
				at = end + 1;
			}
			else if ( c == ',' || c == ' ' || c == '\n' )
				++at;
			else
			{
				std::size_t end = at;
				while ( end < text.size() && text[ end ] != ',' && text[ end ] != ';'
				        && text[ end ] != ' ' && text[ end ] != '\n' )
					++end;
				out.push_back( text.substr( at, end - at ) );
				at = end;
			}
		}
		return out;
	}

	/// The same as numbers, NaN for a fill. A byte prints with no suffix in
	/// the data section, so this reads every numeric type.
	std::vector<double> numbers( std::string const &text, std::string const &name )
	{
		std::vector<double> out;
		for ( std::string const &one : tokens( text, name ) )
			out.push_back( one == "_" ? std::nan( "" ) : std::strtod( one.c_str(), nullptr ) );
		return out;
	}

	/// A global attribute, as a double, NaN if absent.
	double globalAttribute( std::string const &text, std::string const &name )
	{
		std::string const needle = "\t\t:" + name + " = ";
		std::size_t const at = text.find( needle );
		if ( at == std::string::npos )
			return std::nan( "" );
		return std::strtod( text.c_str() + at + needle.size(), nullptr );
	}

	/// A curve whose points do not survive a float, so that a writer storing
	/// them single precision is caught: 1/3 and its friends.
	meq::PermeabilityCurve steelCurve()
	{
		meq::PermeabilityCurve curve;
		curve.material = "yoke";
		curve.fieldStrength = { 0.0, 1.0/3.0, 100.0/7.0, 1.0e3*std::sqrt( 2.0 ), 1.0e5 };
		curve.fluxDensity = { 0.0, 0.1/3.0, 0.9/7.0, 1.4142135623730951, 2.15 };
		return curve;
	}

	/// The table ER-1 asks for: two wired coils and one unwired, a passive
	/// wall, a linear iron, and an iron with a curve -- in an order that puts
	/// the curve's row LAST, so an index-for-row confusion shows.
	MaterialTable reservedKinds()
	{
		MaterialTable table;
		meq::Circuit supply;
		supply.name = "P";
		supply.current = 1.0e5;
		supply.supplied = true;
		int const p = table.addCircuit( supply );
		meq::Circuit grouping;
		grouping.name = meq::defaultCircuitName;
		int const d = table.addCircuit( grouping );

		table.addCoil( meq::Coil( 1.0, 0.5, 0.1, 0.1, 2.0e5 ), "upper", p );
		table.addCoil( meq::Coil( 1.0, -0.5, 0.1, 0.1, -2.0e5 ), "lower", p );
		table.addCoil( meq::Coil( 0.4, 0.0, 0.05, 0.2, 3.0e4 ), "solenoid", d );

		Material wall;
		wall.kind = MaterialKind::Passive;
		wall.name = "wall";
		wall.resistivity = 6.9e-7;
		table.addMaterial( wall );

		Material linear;
		linear.kind = MaterialKind::Ferromagnetic;
		linear.name = "core";
		linear.relativePermeability = 750.0;
		table.addMaterial( linear );

		Material saturating;
		saturating.kind = MaterialKind::Ferromagnetic;
		saturating.name = "yoke";
		saturating.relativePermeability = std::nan( "" );
		saturating.permeabilityCurve = table.addPermeabilityCurve( steelCurve() );
		table.addMaterial( saturating );
		return table;
	}

	/// Three elements: plasma whole, a plasma boundary, and the curve iron
	/// beside the wall -- so the iron is in the ragged record as well as the
	/// table.
	ElementMaterials threeElements( MaterialTable const &table )
	{
		ElementMaterials record;
		record.addElement( table, { { MaterialTable::plasma, 0.5 } }, false, false, 1 );
		record.addElement( table, { { MaterialTable::plasma, 0.25 },
		                            { MaterialTable::vacuum, 0.75 } }, false, false, 1 );
		record.addElement( table, { { 7, 0.125 }, { 5, 0.375 } }, false, true, 30 );
		return record;
	}
}

BOOST_AUTO_TEST_CASE( theReservedKindsAreReadableBeforeTheyAreProduced )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	MaterialTable const table = reservedKinds();
	BOOST_TEST_REQUIRE( table.size() == 8 );
	ElementMaterials const record = threeElements( table );

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 2, 2, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	meq::GridSampler const sampler( mesh, 0.1, 0.9, 3, 0.1, 0.9, 3 );

	Scratch const file( "meq-output-tests-materials.nc" );
	{
		meq::NetCDFWriter writer( file.path, sampler );
		writer.materials( table, record );
		writer.close();
	}
	std::string const text = dump( file.path );
	BOOST_TEST_REQUIRE( !text.empty(), "ncdump could not read the file the writer made" );

	// THE TABLE, every column of §2.1, with each kind's fills.
	std::vector<std::string> const names = tokens( text, "material_name" );
	std::vector<std::string> const wanted = { "vacuum", "plasma", "upper", "lower",
	                                          "solenoid", "wall", "core", "yoke" };
	BOOST_TEST_REQUIRE( names.size() == wanted.size() );
	for ( std::size_t m = 0; m < wanted.size(); ++m )
		BOOST_TEST( names[ m ] == wanted[ m ] );

	std::vector<double> const kind = numbers( text, "material_kind" );
	std::vector<double> const circuit = numbers( text, "material_circuit" );
	std::vector<double> const current = numbers( text, "material_current" );
	std::vector<double> const rho = numbers( text, "material_resistivity" );
	std::vector<double> const mu = numbers( text, "material_relative_permeability" );
	std::vector<double> const curveOf = numbers( text, "material_permeability_curve" );
	BOOST_TEST_REQUIRE( kind.size() == 8u );
	BOOST_TEST_REQUIRE( circuit.size() == 8u );
	BOOST_TEST_REQUIRE( current.size() == 8u );
	BOOST_TEST_REQUIRE( rho.size() == 8u );
	BOOST_TEST_REQUIRE( mu.size() == 8u );
	BOOST_TEST_REQUIRE( curveOf.size() == 8u );

	BOOST_TEST( kind[ 5 ] == 3.0 );    // passive
	BOOST_TEST( kind[ 6 ] == 4.0 );    // ferromagnetic
	BOOST_TEST( kind[ 7 ] == 4.0 );
	BOOST_TEST( circuit[ 2 ] == 0.0 );
	BOOST_TEST( circuit[ 4 ] == 1.0 );
	BOOST_TEST( circuit[ 5 ] == -1.0 );
	BOOST_TEST( current[ 3 ] == -2.0e5 );
	BOOST_TEST( std::isnan( current[ 5 ] ) );
	BOOST_TEST( rho[ 5 ] == 6.9e-7 );
	BOOST_TEST( std::isnan( rho[ 6 ] ) );
	BOOST_TEST( mu[ 5 ] == 1.0 );
	BOOST_TEST( mu[ 6 ] == 750.0 );
	// A CURVE OR A CONSTANT, NEVER BOTH: the curve row's constant is the fill.
	BOOST_TEST( std::isnan( mu[ 7 ] ) );
	BOOST_TEST( curveOf[ 6 ] == -1.0 );
	BOOST_TEST( curveOf[ 7 ] == 0.0 );

	// THE CIRCUITS: one supplied, one grouping whose current is the fill.
	std::vector<std::string> const circuitNames = tokens( text, "circuit_name" );
	BOOST_TEST_REQUIRE( circuitNames.size() == 2u );
	BOOST_TEST( circuitNames[ 0 ] == "P" );
	BOOST_TEST( circuitNames[ 1 ] == meq::defaultCircuitName );
	std::vector<double> const circuitCurrent = numbers( text, "circuit_current" );
	std::vector<double> const supplied = numbers( text, "circuit_supplied" );
	BOOST_TEST_REQUIRE( circuitCurrent.size() == 2u );
	BOOST_TEST_REQUIRE( supplied.size() == 2u );
	BOOST_TEST( circuitCurrent[ 0 ] == 1.0e5 );
	BOOST_TEST( std::isnan( circuitCurrent[ 1 ] ) );
	BOOST_TEST( supplied[ 0 ] == 1.0 );
	BOOST_TEST( supplied[ 1 ] == 0.0 );

	// THE CHECK docs/output.rst TELLS A READER TO MAKE: the curve's own name is
	// the name of the row that points to it.
	std::vector<std::string> const owners = tokens( text, "permeability_curve_material" );
	BOOST_TEST_REQUIRE( owners.size() == 1u );
	int const pointing = 7;
	BOOST_TEST_REQUIRE( ( curveOf[ pointing ] >= 0.0
	                      && curveOf[ pointing ] < static_cast<double>( owners.size() ) ),
	                    "the iron's material_permeability_curve is " << curveOf[ pointing ]
	                    << " of " << owners.size() << " curves" );
	BOOST_TEST( owners[ static_cast<std::size_t>( curveOf[ pointing ] ) ]
	            == names[ pointing ],
	            "permeability_curve_material reads \"" << owners.front()
	            << "\" and the row pointing to it is \"" << names[ pointing ] << "\"" );

	// AND THE POINTS ROUND-TRIP EXACTLY, as a contiguous ragged array.
	std::vector<double> const pointCount = numbers( text, "permeability_curve_point_count" );
	std::vector<double> const h = numbers( text, "permeability_curve_h" );
	std::vector<double> const b = numbers( text, "permeability_curve_b" );
	meq::PermeabilityCurve const curve = steelCurve();
	BOOST_TEST_REQUIRE( pointCount.size() == 1u );
	BOOST_TEST( pointCount[ 0 ] == static_cast<double>( curve.fieldStrength.size() ) );
	BOOST_TEST_REQUIRE( h.size() == curve.fieldStrength.size() );
	BOOST_TEST_REQUIRE( b.size() == curve.fluxDensity.size() );
	for ( std::size_t i = 0; i < h.size(); ++i )
	{
		BOOST_TEST( h[ i ] == curve.fieldStrength[ i ], "H[" << i << "] did not round-trip" );
		BOOST_TEST( b[ i ] == curve.fluxDensity[ i ], "B[" << i << "] did not round-trip" );
	}
	for ( std::string const &needle :
	      { std::string( "permeability_curve_point_count:sample_dimension = \"permeability_curve_point\"" ),
	        std::string( "permeability_curve_b:interpolation = \"" )
	            + meq::permeabilityInterpolation + "\"" } )
		BOOST_TEST( text.find( needle ) != std::string::npos, "missing: " << needle );
}

BOOST_AUTO_TEST_CASE( theElementRecordIsAContiguousRaggedArray )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	MaterialTable const table = reservedKinds();
	ElementMaterials const record = threeElements( table );

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 2, 2, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	meq::GridSampler const sampler( mesh, 0.1, 0.9, 3, 0.1, 0.9, 3 );

	Scratch const file( "meq-output-tests-record.nc" );
	{
		meq::NetCDFWriter writer( file.path, sampler );
		writer.materials( table, record );
		writer.close();
	}
	std::string const text = dump( file.path );
	BOOST_TEST_REQUIRE( !text.empty() );

	std::vector<double> const count = numbers( text, "element_material_count" );
	std::vector<double> const material = numbers( text, "element_material" );
	std::vector<double> const fraction = numbers( text, "element_fraction" );
	BOOST_TEST_REQUIRE( count.size() == 3u );
	BOOST_TEST( count[ 0 ] == 1.0 );
	BOOST_TEST( count[ 1 ] == 2.0 );
	BOOST_TEST( count[ 2 ] == 2.0 );
	BOOST_TEST_REQUIRE( material.size() == 5u );
	BOOST_TEST_REQUIRE( fraction.size() == 5u );
	// Sorted by material within an element: the wall (5) before the yoke (7).
	BOOST_TEST( material[ 3 ] == 5.0 );
	BOOST_TEST( material[ 4 ] == 7.0 );
	BOOST_TEST( fraction[ 3 ] == 0.75 );
	BOOST_TEST( fraction[ 4 ] == 0.25 );
	BOOST_TEST( text.find( "element_material_count:sample_dimension = \"element_entry\"" )
	            != std::string::npos );

	std::vector<double> const region = numbers( text, "element_region" );
	std::vector<double> const plasma = numbers( text, "element_plasma_fraction" );
	std::vector<double> const excluded = numbers( text, "element_excluded" );
	std::vector<double> const attribute = numbers( text, "element_attribute" );
	BOOST_TEST_REQUIRE( region.size() == 3u );
	BOOST_TEST( region[ 0 ] == 1.0 );   // plasma
	BOOST_TEST( region[ 1 ] == 5.0 );   // plasma_boundary
	BOOST_TEST( region[ 2 ] == 7.0 );   // two conductors: mixed
	BOOST_TEST_REQUIRE( plasma.size() == 3u );
	BOOST_TEST( plasma[ 1 ] == 0.25 );
	BOOST_TEST_REQUIRE( excluded.size() == 3u );
	BOOST_TEST( excluded[ 2 ] == 1.0 );
	BOOST_TEST_REQUIRE( attribute.size() == 3u );
	BOOST_TEST( attribute[ 2 ] == 30.0 );
	BOOST_TEST( text.find( "element_region:flag_meanings = \"vacuum plasma coil passive "
	                       "ferromagnetic plasma_boundary conductor_boundary mixed\"" )
	            != std::string::npos );

	BOOST_TEST( globalAttribute( text, "region_elements_plasma" ) == 1.0 );
	BOOST_TEST( globalAttribute( text, "region_elements_plasma_boundary" ) == 1.0 );
	BOOST_TEST( globalAttribute( text, "region_elements_mixed" ) == 1.0 );
	BOOST_TEST( globalAttribute( text, "region_elements_coil" ) == 0.0 );
	BOOST_TEST( globalAttribute( text, "region_faults" ) == 0.0 );
	BOOST_TEST( globalAttribute( text, "material_area_plasma" ) == 0.75 );
	BOOST_TEST( globalAttribute( text, "material_area_yoke" ) == 0.125 );
	BOOST_TEST( globalAttribute( text, "material_area_upper" ) == 0.0 );
}

/// A coil may share its name with another coil or with a fixed row, and a
/// second attribute of one name would overwrite the first without a word.
BOOST_AUTO_TEST_CASE( twoRowsOfOneNameKeepTwoAreaAttributes )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	MaterialTable table;
	table.addCoil( meq::Coil( 1.0, 0.5, 0.1, 0.1, 1.0 ), "vacuum" );
	table.addCoil( meq::Coil( 1.0, -0.5, 0.1, 0.1, 1.0 ), "PF 1" );
	ElementMaterials record;
	record.addElement( table, { { MaterialTable::vacuum, 0.5 } }, false, false, 1 );
	record.addElement( table, { { 2, 0.25 } }, false, false, 1 );
	record.addElement( table, { { 3, 0.125 } }, false, false, 1 );

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 2, 2, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	meq::GridSampler const sampler( mesh, 0.1, 0.9, 3, 0.1, 0.9, 3 );
	Scratch const file( "meq-output-tests-names.nc" );
	{
		meq::NetCDFWriter writer( file.path, sampler );
		writer.materials( table, record );
		writer.close();
	}
	std::string const text = dump( file.path );
	BOOST_TEST( globalAttribute( text, "material_area_vacuum" ) == 0.5 );
	BOOST_TEST( globalAttribute( text, "material_area_vacuum_2" ) == 0.25 );
	BOOST_TEST( globalAttribute( text, "material_area_PF_1" ) == 0.125 );

	// And a record made with another table is refused rather than mislabelled.
	MaterialTable const shorter;
	meq::NetCDFWriter writer( file.path, sampler );
	BOOST_CHECK_THROW( writer.materials( shorter, record ), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( theNodeMaterialIsAnIntOnTheGridWithAFill )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 2, 2, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	meq::GridSampler const sampler( mesh, 0.1, 0.9, 3, 0.1, 0.9, 2 );
	Scratch const file( "meq-output-tests-node.nc" );
	{
		meq::NetCDFWriter writer( file.path, sampler );
		writer.nodeMaterial( { 0, 1, 2, -1, 7, 1 } );
		// One per node, or refused rather than written short.
		BOOST_CHECK_THROW( writer.nodeMaterial( { 0, 1 } ), std::runtime_error );
		writer.close();
	}
	std::string const text = dump( file.path );
	BOOST_TEST( text.find( "int node_material(Z, R)" ) != std::string::npos, text );
	BOOST_TEST( text.find( "node_material:_FillValue = -1" ) != std::string::npos, text );
	// ( Z, R ) with R fastest, and the -1 read back as the fill.
	std::vector<std::string> const values = tokens( text, "node_material" );
	BOOST_TEST_REQUIRE( values.size() == 6u );
	BOOST_TEST( values[ 0 ] == "0" );
	BOOST_TEST( values[ 2 ] == "2" );
	BOOST_TEST( values[ 3 ] == "_" );
	BOOST_TEST( values[ 4 ] == "7" );
}

namespace
{
	/// Every element of a record whose values differ from element to
	/// element, so a field written against the wrong element, or averaged
	/// with a neighbour, cannot pass by coincidence.
	ElementMaterials variedRecord( MaterialTable const &table, mfem::Mesh const &mesh )
	{
		ElementMaterials record;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			double const share = 0.05 + 0.9*( e + 1.0 )/( mesh.GetNE() + 1.0 );
			switch ( e % 5 )
			{
				case 0:
					record.addElement( table, { { MaterialTable::plasma, 1.0 } },
					                   false, false, mesh.GetAttribute( e ) );
					break;
				case 1:
					record.addElement( table, { { MaterialTable::plasma, share },
					                            { MaterialTable::vacuum, 1.0 - share } },
					                   false, false, mesh.GetAttribute( e ) );
					break;
				case 2:
					// Two conductors, the larger the HIGHER row.
					record.addElement( table, { { 2, 0.1*share }, { 6, 0.3*share },
					                            { MaterialTable::vacuum, 1.0 - 0.4*share } },
					                   false, false, mesh.GetAttribute( e ) );
					break;
				case 3:
					// A tie, which goes to the lower row.
					record.addElement( table, { { 4, 0.25 }, { 3, 0.25 },
					                            { MaterialTable::vacuum, 0.5 } },
					                   false, false, mesh.GetAttribute( e ) );
					break;
				default:
					record.addElement( table, { { MaterialTable::vacuum, 1.0 } },
					                   false, false, mesh.GetAttribute( e ) );
			}
		}
		return record;
	}

	/// What the cell data say, read back from @a vtu, against @a record:
	/// one value per cell at every point of it, no point shared, and the
	/// attribute as cell data. Returns the number of points.
	int checkCellData( std::string const &vtu, MaterialTable const &table,
	                   ElementMaterials const &record, mfem::Mesh &mesh )
	{
		std::vector<double> const connectivity = meq::tests::readVtuArray( vtu, "connectivity" );
		std::vector<double> const offsets = meq::tests::readVtuArray( vtu, "offsets" );
		std::vector<double> const points = meq::tests::readVtuArray( vtu, "Points" );
		BOOST_TEST_REQUIRE( offsets.size() == static_cast<std::size_t>( mesh.GetNE() ) );
		int const pointCount = static_cast<int>( points.size()/3 );

		// NO POINT BELONGS TO TWO CELLS, which is what makes a point array of
		// an order-0 field a cell value: nothing can be averaged across a face.
		std::vector<int> owners( static_cast<std::size_t>( pointCount ), 0 );
		for ( double const p : connectivity )
			++owners.at( static_cast<std::size_t>( p ) );
		int shared = 0;
		for ( int const n : owners )
			shared += n != 1 ? 1 : 0;
		BOOST_TEST( shared == 0, shared << " of " << pointCount << " points are not "
		            "in exactly one cell, so a per-element field is smeared across "
		            "faces" );

		std::vector<double> expectedRegion, expectedPlasma, expectedConductor, expectedShare;
		for ( int e = 0; e < record.elementCount(); ++e )
		{
			meq::MaterialFraction const conductor = record.dominantConductor( e, table );
			expectedRegion.push_back( static_cast<int>( record.region[ static_cast<std::size_t>( e ) ] ) );
			expectedPlasma.push_back( record.plasmaFraction[ static_cast<std::size_t>( e ) ] );
			expectedConductor.push_back( conductor.material );
			expectedShare.push_back( conductor.fraction );
		}
		struct Named { char const *name; std::vector<double> const *expected; };
		for ( Named const &field : { Named{ "region", &expectedRegion },
		                             Named{ "plasma_fraction", &expectedPlasma },
		                             Named{ "conductor_material", &expectedConductor },
		                             Named{ "conductor_fraction", &expectedShare } } )
		{
			std::vector<double> const values = meq::tests::readVtuArray( vtu, field.name );
			BOOST_TEST_REQUIRE( values.size() == static_cast<std::size_t>( pointCount ),
			                    field.name << " has " << values.size() << " values for "
			                    << pointCount << " points" );
			BOOST_TEST( meq::tests::vtuArrayIsPointData( vtu, field.name ) );
			// The cell data are GridFunctions, so they hold the record's double
			// in this build's real_t: in single a fraction arrives rounded to
			// float, exactly. Comparing against that rounding keeps the check
			// exact in both formats, and in double the cast is the identity.
			int wrong = 0;
			for ( std::size_t c = 0; c < offsets.size(); ++c )
			{
				std::size_t const begin = c == 0 ? 0 : static_cast<std::size_t>( offsets[ c - 1 ] );
				std::size_t const end = static_cast<std::size_t>( offsets[ c ] );
				double const want = static_cast<double>(
					static_cast<mfem::real_t>( ( *field.expected )[ c ] ) );
				for ( std::size_t k = begin; k < end; ++k )
					if ( values[ static_cast<std::size_t>( connectivity[ k ] ) ] != want )
					{
						++wrong;
						break;
					}
			}
			BOOST_TEST( wrong == 0, field.name << ": " << wrong << " of "
			            << offsets.size() << " cells do not carry their element's "
			            "value at every one of their points" );
		}

		// THE ATTRIBUTE IS THE COLLECTION'S OWN CELL DATA, NOT A FIELD OF OURS.
		std::vector<double> const attribute = meq::tests::readVtuArray( vtu, "attribute" );
		BOOST_TEST_REQUIRE( attribute.size() == offsets.size() );
		BOOST_TEST( !meq::tests::vtuArrayIsPointData( vtu, "attribute" ) );
		for ( std::size_t c = 0; c < attribute.size(); ++c )
			BOOST_TEST( attribute[ c ] == record.attribute[ c ] );
		std::printf( "  %s: %zu cells, %d points, %d shared\n", vtu.c_str(),
		             offsets.size(), pointCount, shared );
		return pointCount;
	}

	/// A potential and a field on @a mesh, to stand beside the cell data.
	struct Fields
	{
		mfem::L2_FECollection collection;
		mfem::FiniteElementSpace scalar, vector;
		mfem::GridFunction psi, field;
		explicit Fields( mfem::Mesh &mesh )
			: collection( 2, 2 ), scalar( &mesh, &collection ),
			  vector( &mesh, &collection, 2 ), psi( &scalar ), field( &vector )
		{
			mfem::FunctionCoefficient shape( []( mfem::Vector const &x )
			                                 { return x( 0 )*x( 1 ); } );
			psi.ProjectCoefficient( shape );
			field = 0.0;
		}
	};
}

/*
 * ER-2: THE VTK CARRIES WHAT EACH ELEMENT IS MADE OF, ONE VALUE PER CELL.
 *
 * On a CURVED mesh at the driver's high-order output, because that is where
 * an order-0 field could be smeared: if ParaViewDataCollection shared a
 * Lagrange point between neighbours, the point would carry one element's
 * value for both. Read back by hand (VtuReader.hpp), so the check shares no
 * code with the writer. And a frame series, whose second frame is a refined
 * mesh with a record of its own.
 */
BOOST_AUTO_TEST_CASE( theVtkCarriesOneMaterialValuePerCell )
{
	MaterialTable const table = reservedKinds();
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 4, 3, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	for ( int e = 0; e < mesh.GetNE(); ++e )
		mesh.SetAttribute( e, 1 + e % 3 );
	mesh.SetAttributes();
	mesh.SetCurvature( 3 );
	// Bend it, so the cells are genuinely curved Lagrange cells.
	mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
	{
		out( 0 ) = 1.0 + in( 0 ) + 0.05*std::sin( 3.0*in( 1 ) );
		out( 1 ) = in( 1 ) + 0.05*std::sin( 3.0*in( 0 ) );
	} );

	ElementMaterials const record = variedRecord( table, mesh );
	meq::ElementCellData const cells( mesh, table, record );
	Fields const fields( mesh );

	std::filesystem::path const directory = std::filesystem::temp_directory_path()
		/ ( "meq-output-tests-vtk-" + std::to_string( ::getpid() ) );
	std::filesystem::create_directories( directory );
	std::string const stem = ( directory / "cells" ).string();
	meq::writeVtu( stem, mesh, fields.psi, fields.field, 4, &cells );
	std::string const vtu = stem + "/Cycle000000/proc000000.vtu";
	BOOST_TEST_REQUIRE( std::filesystem::exists( vtu ) );
	int const points = checkCellData( vtu, table, record, mesh );
	BOOST_TEST( points > 3*mesh.GetNE(), "the output is not high order: "
	            << points << " points for " << mesh.GetNE() << " cells" );

	// THE CONDUCTOR COLUMNS, spelled out once: element 2 holds rows 2 and 6
	// and the higher row is the larger, element 3 ties rows 3 and 4, element
	// 4 holds none.
	BOOST_TEST( record.dominantConductor( 2, table ).material == 6 );
	BOOST_TEST( record.dominantConductor( 3, table ).material == 3 );
	BOOST_TEST( record.dominantConductor( 4, table ).material == -1 );
	BOOST_TEST( record.dominantConductor( 4, table ).fraction == 0.0 );

	// ---- the series: two frames, the second on a refined mesh ----
	{
		meq::VtuSeries series( ( directory / "frames" ).string(), 3 );
		series.append( mesh, fields.psi, fields.field, 0, 0.0, &cells );

		mfem::Mesh refined( mesh );
		refined.UniformRefinement();
		ElementMaterials const second = variedRecord( table, refined );
		meq::ElementCellData const secondCells( refined, table, second );
		Fields const secondFields( refined );
		series.append( refined, secondFields.psi, secondFields.field, 1, 1.0, &secondCells );

		checkCellData( ( directory / "frames_cycles/Cycle000000/proc000000.vtu" ).string(),
		               table, record, mesh );
		checkCellData( ( directory / "frames_cycles/Cycle000001/proc000000.vtu" ).string(),
		               table, second, refined );
	}

	// A RECORD OF ANOTHER MESH, OR ANOTHER TABLE, IS REFUSED.
	mfem::Mesh other = mfem::Mesh::MakeCartesian2D( 2, 2, mfem::Element::TRIANGLE,
	                                                false, 1.0, 1.0 );
	BOOST_CHECK_THROW( meq::ElementCellData( other, table, record ), std::invalid_argument );
	MaterialTable const plain;
	BOOST_CHECK_THROW( meq::ElementCellData( mesh, plain, record ), std::invalid_argument );

	std::filesystem::remove_all( directory );
}
