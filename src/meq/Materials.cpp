#include "Materials.hpp"

#include "Config.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace meq
{
	bool isConductor( MaterialKind kind )
	{
		return kind == MaterialKind::Coil || kind == MaterialKind::Passive
		       || kind == MaterialKind::Ferromagnetic;
	}

	char const *materialKindName( MaterialKind kind )
	{
		switch ( kind )
		{
			case MaterialKind::Vacuum:        return "vacuum";
			case MaterialKind::Plasma:        return "plasma";
			case MaterialKind::Coil:          return "coil";
			case MaterialKind::Passive:       return "passive";
			case MaterialKind::Ferromagnetic: return "ferromagnetic";
			default:                          break;
		}
		throw std::invalid_argument( "meq::materialKindName: not a MaterialKind" );
	}

	char const *elementRegionName( ElementRegion region )
	{
		switch ( region )
		{
			case ElementRegion::Vacuum:            return "vacuum";
			case ElementRegion::Plasma:            return "plasma";
			case ElementRegion::Coil:              return "coil";
			case ElementRegion::Passive:           return "passive";
			case ElementRegion::Ferromagnetic:     return "ferromagnetic";
			case ElementRegion::PlasmaBoundary:    return "plasma_boundary";
			case ElementRegion::ConductorBoundary: return "conductor_boundary";
			case ElementRegion::Mixed:             return "mixed";
			default:                               break;
		}
		throw std::invalid_argument( "meq::elementRegionName: not an ElementRegion" );
	}

	MaterialTable::MaterialTable()
	{
		Material vacuumRow;
		vacuumRow.kind = MaterialKind::Vacuum;
		vacuumRow.name = "vacuum";
		rows.push_back( vacuumRow );

		Material plasmaRow;
		plasmaRow.kind = MaterialKind::Plasma;
		plasmaRow.name = "plasma";
		rows.push_back( plasmaRow );
	}

	int MaterialTable::addCoil( Coil const &shape, std::string const &name,
	                            int circuit )
	{
		if ( circuit < -1 || circuit >= static_cast< int >( circuitRows.size() ) )
			throw std::invalid_argument( "meq::MaterialTable::addCoil: circuit "
			                             + std::to_string( circuit ) + " is not a row of the "
			                             + std::to_string( circuitRows.size() ) + "-row circuit table" );
		Material row;
		row.kind = MaterialKind::Coil;
		row.name = name;
		row.circuit = circuit;
		row.current = shape.current();
		rows.push_back( row );
		coilShapes.push_back( shape );
		coilRows.push_back( static_cast< int >( rows.size() ) - 1 );
		return coilRows.back();
	}

	int MaterialTable::addMaterial( Material const &row )
	{
		// A coil without its rectangle would be a conductor nothing can find,
		// and a second vacuum or plasma row would make the two fixed indices a
		// convention rather than a fact.
		if ( row.kind == MaterialKind::Vacuum || row.kind == MaterialKind::Plasma
		     || row.kind == MaterialKind::Coil )
			throw std::invalid_argument( std::string( "meq::MaterialTable::addMaterial: a " )
			                             + materialKindName( row.kind )
			                             + " row has its own constructor" );

		std::string const where = "meq::MaterialTable::addMaterial: \"" + row.name + "\" ";
		// NEITHER RESERVED KIND IS ON A SUPPLY. A passive structure's current
		// is induced, and an iron's is none.
		if ( row.circuit != -1 )
			throw std::invalid_argument( where + "is a " + materialKindName( row.kind )
			                             + " row and has no circuit" );
		if ( row.kind == MaterialKind::Passive )
		{
			if ( row.relativePermeability != 1.0 || row.permeabilityCurve != -1 )
				throw std::invalid_argument( where + "is passive, whose relative permeability is 1 and which has no curve; a magnetic structure is ferromagnetic" );
		}
		else if ( row.permeabilityCurve == -1 )
		{
			if ( !std::isfinite( row.relativePermeability ) || !( row.relativePermeability > 0.0 ) )
				throw std::invalid_argument( where + "is a linear iron and needs a finite, positive relative permeability" );
		}
		else
		{
			// A CURVE OR A CONSTANT, NEVER BOTH: the other holds its fill, so
			// that a reader cannot be handed two permeabilities to choose from.
			if ( row.permeabilityCurve < 0
			     || row.permeabilityCurve >= static_cast< int >( curves.size() ) )
				throw std::invalid_argument( where + "points to permeability curve "
				                             + std::to_string( row.permeabilityCurve ) + " of "
				                             + std::to_string( curves.size() ) );
			if ( !std::isnan( row.relativePermeability ) )
				throw std::invalid_argument( where + "has a permeability curve, so its constant relative permeability must be NaN, the fill; a row holds one or the other" );
			std::string const &owner =
				curves[ static_cast< std::size_t >( row.permeabilityCurve ) ].material;
			if ( owner != row.name )
				throw std::invalid_argument( where + "points to the curve of \"" + owner
				                             + "\"; a curve belongs to one material, and two regions of one steel hold two copies" );
		}
		rows.push_back( row );
		return static_cast< int >( rows.size() ) - 1;
	}

	int MaterialTable::addCircuit( Circuit const &row )
	{
		if ( row.name.empty() )
			throw std::invalid_argument( "meq::MaterialTable::addCircuit: a circuit needs a name" );
		for ( Circuit const &existing : circuitRows )
			if ( existing.name == row.name )
				throw std::invalid_argument( "meq::MaterialTable::addCircuit: \"" + row.name
				                             + "\" is already a circuit" );
		circuitRows.push_back( row );
		return static_cast< int >( circuitRows.size() ) - 1;
	}

	std::vector< Circuit > const &MaterialTable::circuits() const
	{
		return circuitRows;
	}

	int MaterialTable::addPermeabilityCurve( PermeabilityCurve const &curve )
	{
		std::string const where = "meq::MaterialTable::addPermeabilityCurve: \""
		                          + curve.material + "\" ";
		if ( curve.material.empty() )
			throw std::invalid_argument( "meq::MaterialTable::addPermeabilityCurve: a curve names its material" );
		if ( curve.fieldStrength.size() != curve.fluxDensity.size()
		     || curve.fieldStrength.size() < 2 )
			throw std::invalid_argument( where + "needs at least two ( H, B ) points, one B per H" );
		if ( curve.fieldStrength.front() != 0.0 || curve.fluxDensity.front() != 0.0 )
			throw std::invalid_argument( where + "starts at ( 0, 0 )" );
		for ( std::size_t i = 1; i < curve.fieldStrength.size(); ++i )
			if ( !std::isfinite( curve.fieldStrength[ i ] ) || !std::isfinite( curve.fluxDensity[ i ] )
			     || !( curve.fieldStrength[ i ] > curve.fieldStrength[ i - 1 ] )
			     || curve.fluxDensity[ i ] < curve.fluxDensity[ i - 1 ] )
				throw std::invalid_argument( where + "must be strictly increasing in H, with B finite and non-decreasing; point "
				                             + std::to_string( i ) + " is not" );
		curves.push_back( curve );
		return static_cast< int >( curves.size() ) - 1;
	}

	std::vector< PermeabilityCurve > const &MaterialTable::permeabilityCurves() const
	{
		return curves;
	}

	MaterialTable buildMaterialTable( CoilConfig const &coils,
	                                  std::vector< double > const &currents )
	{
		if ( !currents.empty() && currents.size() != coils.coils.size() )
			throw std::invalid_argument( "meq::buildMaterialTable: "
			                             + std::to_string( currents.size() ) + " currents for "
			                             + std::to_string( coils.coils.size() ) + " coils" );
		auto currentOf = [ & ]( std::size_t i )
		{
			return currents.empty() ? coils.coils[ i ].current : currents[ i ];
		};

		MaterialTable table;
		std::vector< int > circuitOf( coils.coils.size(), -1 );
		for ( CircuitParameters const &block : coils.circuits )
		{
			Circuit row;
			row.name = block.name;
			row.supplied = true;
			row.current = block.current;
			if ( !block.members.empty() )
				row.current = currentOf( block.members.front().index )
				              /block.members.front().multiplier;
			int const index = table.addCircuit( row );
			for ( CircuitMember const &member : block.members )
				circuitOf.at( member.index ) = index;
		}

		// THE GROUPING, ONLY WHEN IT HAS A MEMBER.
		if ( std::find( circuitOf.begin(), circuitOf.end(), -1 ) != circuitOf.end() )
		{
			Circuit grouping;
			grouping.name = defaultCircuitName;
			int const index = table.addCircuit( grouping );
			std::replace( circuitOf.begin(), circuitOf.end(), -1, index );
		}

		for ( std::size_t i = 0; i < coils.coils.size(); ++i )
		{
			CoilParameters const &block = coils.coils[ i ];
			table.addCoil( Coil( block.centreR, block.centreZ, block.halfWidth,
			                     block.halfHeight, currentOf( i ) ),
			               block.name, circuitOf[ i ] );
		}
		return table;
	}

	int MaterialTable::size() const
	{
		return static_cast< int >( rows.size() );
	}

	Material const &MaterialTable::material( int index ) const
	{
		if ( index < 0 || index >= size() )
			throw std::out_of_range( "meq::MaterialTable::material: index "
			                         + std::to_string( index ) + " of "
			                         + std::to_string( size() ) );
		return rows[ static_cast< std::size_t >( index ) ];
	}

	std::vector< Material > const &MaterialTable::materials() const
	{
		return rows;
	}

	std::size_t MaterialTable::coilCount() const
	{
		return coilShapes.size();
	}

	Coil const &MaterialTable::coil( std::size_t which ) const
	{
		return coilShapes.at( which );
	}

	int MaterialTable::coilMaterial( std::size_t which ) const
	{
		return coilRows.at( which );
	}

	int MaterialTable::conductorContaining( double radius, double z ) const
	{
		for ( std::size_t i = 0; i < coilShapes.size(); ++i )
			if ( coilShapes[ i ].contains( radius, z ) )
				return coilRows[ i ];
		return -1;
	}

	ElementRegion deriveRegion( std::vector< MaterialFraction > const &entries,
	                            MaterialTable const &table )
	{
		if ( entries.empty() )
			throw std::invalid_argument( "meq::deriveRegion: an element holds at least one material" );

		// ONE MATERIAL IS ITS OWN KIND. The ElementRegion numbers of the five
		// pure regions are the MaterialKind numbers on purpose.
		if ( entries.size() == 1 )
			return static_cast< ElementRegion >(
				static_cast< int >( table.material( entries.front().material ).kind ) );

		bool vacuum = false;
		bool plasma = false;
		bool conductor = false;
		for ( MaterialFraction const &entry : entries )
		{
			MaterialKind const kind = table.material( entry.material ).kind;
			vacuum = vacuum || kind == MaterialKind::Vacuum;
			plasma = plasma || kind == MaterialKind::Plasma;
			conductor = conductor || isConductor( kind );
		}

		// §2.1b: plasma beside vacuum is the boundary nobody can mesh to;
		// copper beside vacuum is a mesh not aligned to it. Plasma beside
		// copper, with or without vacuum, is a fault or a cell to split.
		if ( vacuum && plasma && !conductor )
			return ElementRegion::PlasmaBoundary;
		if ( vacuum && conductor && !plasma )
			return ElementRegion::ConductorBoundary;
		return ElementRegion::Mixed;
	}

	void ElementMaterials::addElement( MaterialTable const &table,
	                                   std::vector< Weight > weights,
	                                   bool faultedIn, bool excludedIn,
	                                   int attributeIn )
	{
		if ( weights.empty() )
			throw std::invalid_argument( "meq::ElementMaterials::addElement: an element has at least one quadrature point" );

		std::sort( weights.begin(), weights.end(),
		           []( Weight const &a, Weight const &b )
		           { return a.material < b.material; } );

		// MERGED BY MATERIAL, so a caller may hand one weight per point.
		std::vector< Weight > merged;
		double total = 0.0;
		for ( Weight const &w : weights )
		{
			if ( w.material < 0 || w.material >= table.size() )
				throw std::invalid_argument( "meq::ElementMaterials::addElement: material "
				                             + std::to_string( w.material )
				                             + " is not a row of the table" );
			if ( !std::isfinite( w.weight ) || w.weight < 0.0 )
				throw std::invalid_argument( "meq::ElementMaterials::addElement: a weight must be finite and non-negative" );
			if ( !merged.empty() && merged.back().material == w.material )
				merged.back().weight += w.weight;
			else
				merged.push_back( w );
			total += w.weight;
		}
		if ( !( total > 0.0 ) )
			throw std::invalid_argument( "meq::ElementMaterials::addElement: the element has no area" );

		if ( materialArea.empty() )
			materialArea.assign( static_cast< std::size_t >( table.size() ), 0.0 );

		std::vector< MaterialFraction > own;
		double plasmaShare = 0.0;
		for ( Weight const &w : merged )
		{
			if ( w.weight == 0.0 )
				continue;
			// w/total, so that one material reads EXACTLY 1: x/x is 1 in
			// floating point, where 1 - (the others) need not be.
			MaterialFraction entry;
			entry.material = w.material;
			entry.fraction = w.weight/total;
			own.push_back( entry );
			materialArea[ static_cast< std::size_t >( w.material ) ] += w.weight;
			if ( w.material == MaterialTable::plasma )
				plasmaShare = entry.fraction;
		}

		entries.insert( entries.end(), own.begin(), own.end() );
		offsets.push_back( static_cast< int >( entries.size() ) );
		region.push_back( deriveRegion( own, table ) );
		plasmaFraction.push_back( plasmaShare );
		area.push_back( total );
		faulted.push_back( faultedIn ? 1 : 0 );
		excluded.push_back( excludedIn ? 1 : 0 );
		attribute.push_back( attributeIn );
		if ( faultedIn )
			++faults;
	}

	int ElementMaterials::elementCount() const
	{
		return static_cast< int >( offsets.size() ) - 1;
	}

	double ElementMaterials::fraction( int element, int material ) const
	{
		if ( element < 0 || element >= elementCount() )
			throw std::out_of_range( "meq::ElementMaterials::fraction: element "
			                         + std::to_string( element ) + " of "
			                         + std::to_string( elementCount() ) );
		for ( int i = offsets[ static_cast< std::size_t >( element ) ];
		      i < offsets[ static_cast< std::size_t >( element ) + 1 ]; ++i )
			if ( entries[ static_cast< std::size_t >( i ) ].material == material )
				return entries[ static_cast< std::size_t >( i ) ].fraction;
		return 0.0;
	}

	MaterialFraction ElementMaterials::dominantConductor( int element,
	                                                     MaterialTable const &table ) const
	{
		if ( element < 0 || element >= elementCount() )
			throw std::out_of_range( "meq::ElementMaterials::dominantConductor: element "
			                         + std::to_string( element ) + " of "
			                         + std::to_string( elementCount() ) );
		MaterialFraction best{ -1, 0.0 };
		for ( int i = offsets[ static_cast< std::size_t >( element ) ];
		      i < offsets[ static_cast< std::size_t >( element ) + 1 ]; ++i )
		{
			MaterialFraction const &entry = entries[ static_cast< std::size_t >( i ) ];
			if ( !isConductor( table.material( entry.material ).kind ) )
				continue;
			if ( best.material < 0 || entry.fraction > best.fraction
			     || ( entry.fraction == best.fraction && entry.material < best.material ) )
				best = entry;
		}
		return best;
	}

	int ElementMaterials::regionCount( ElementRegion which ) const
	{
		return static_cast< int >( std::count( region.begin(), region.end(), which ) );
	}
}
