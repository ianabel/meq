#include "ConductorStore.hpp"

#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <netcdf>

namespace meq
{
	namespace
	{
		constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
		constexpr std::uint64_t fnvPrime = 1099511628211ULL;

		constexpr char const *formatName = "meq-conductor-cache";
		constexpr int formatVersion = 2;

		std::uint64_t mixBytes( std::uint64_t digest, std::uint64_t bits )
		{
			for ( int byte = 0; byte < 8; ++byte )
			{
				digest ^= ( bits >> ( 8*byte ) ) & 0xffULL;
				digest *= fnvPrime;
			}
			return digest;
		}

		/// A digest is only ever compared with another digest written by this
		/// same code, so the IEEE bytes are the right thing to mix -- and
		/// -0.0 and +0.0 hashing differently is CORRECT here, because a vertex
		/// at -0.0 and one at +0.0 are the same point and `psi_c` is even in
		/// `R`, so the values agree and a spurious refusal costs a recompute
		/// rather than an answer.
		std::string toHex( std::uint64_t value )
		{
			std::ostringstream out;
			out << std::hex << std::setw( 16 ) << std::setfill( '0' ) << value;
			return out.str();
		}

		std::uint64_t fromHex( std::string const &text, std::string const &path )
		{
			try
			{
				return std::stoull( text, nullptr, 16 );
			}
			catch ( std::exception const & )
			{
				throw std::runtime_error(
					"meq::readConductorCache: " + path + " carries a digest "
					"that is not 16 hexadecimal digits: \"" + text + "\"" );
			}
		}

		std::string textAttribute( netCDF::NcGroup const &group,
		                           std::string const &name,
		                           std::string const &path )
		{
			netCDF::NcGroupAtt attribute = group.getAtt( name );
			if ( attribute.isNull() )
				throw std::runtime_error(
					"meq::readConductorCache: " + path + " claims to be a "
					+ formatName + " but has no \"" + name + "\" attribute" );

			std::string value;
			attribute.getValues( value );
			return value;
		}

		int intAttribute( netCDF::NcGroup const &group,
		                  std::string const &name,
		                  std::string const &path )
		{
			netCDF::NcGroupAtt attribute = group.getAtt( name );
			if ( attribute.isNull() )
				throw std::runtime_error(
					"meq::readConductorCache: " + path + " claims to be a "
					+ formatName + " but has no \"" + name + "\" attribute" );

			int value = 0;
			attribute.getValues( &value );
			return value;
		}

		bool carriesTheFormat( netCDF::NcGroup const &group )
		{
			if ( group.isNull() )
				return false;

			netCDF::NcGroupAtt attribute = group.getAtt( "format" );
			if ( attribute.isNull() )
				return false;

			std::string value;
			attribute.getValues( value );
			return value == formatName;
		}
	}

	bool ConductorCacheSignature::operator==(
		ConductorCacheSignature const &other ) const
	{
		return meshDigest == other.meshDigest
		       && conductorDigest == other.conductorDigest
		       && polynomialDegree == other.polynomialDegree
		       && elements == other.elements
		       && potentialDofs == other.potentialDofs;
	}

	bool ConductorCacheSignature::operator!=(
		ConductorCacheSignature const &other ) const
	{
		return !( *this == other );
	}

	std::string ConductorCacheSignature::difference(
		ConductorCacheSignature const &other ) const
	{
		/*
		 * ORDERED SO THE FIRST LINE IS THE MOST LIKELY CAUSE, which is what a
		 * reader acts on. A coupled run changes the profiles and nothing else,
		 * so a conductor digest that moved means the machine was edited; a mesh
		 * digest that moved means the run re-meshed, which is the ordinary
		 * reason a cache goes stale and is not a mistake.
		 */
		if ( conductorDigest != other.conductorDigest )
			return "the conductors differ ( digest " + toHex( conductorDigest )
			       + " against " + toHex( other.conductorDigest )
			       + " ): a current, a position, a size, the permeability or "
			         "the cross-section quadrature order has changed";

		if ( meshDigest != other.meshDigest )
			return "the mesh differs ( digest " + toHex( meshDigest )
			       + " against " + toHex( other.meshDigest )
			       + " ): psi_c was evaluated at points this mesh does not have";

		if ( polynomialDegree != other.polynomialDegree )
			return "the polynomial degree differs ( "
			       + std::to_string( polynomialDegree ) + " against "
			       + std::to_string( other.polynomialDegree )
			       + " ): the dof positions are the element's nodes and those "
			         "move with the degree";

		if ( elements != other.elements )
			return "the element count differs ( "
			       + std::to_string( elements ) + " against "
			       + std::to_string( other.elements ) + " )";

		if ( potentialDofs != other.potentialDofs )
			return "the potential dof count differs ( "
			       + std::to_string( potentialDofs ) + " against "
			       + std::to_string( other.potentialDofs ) + " )";

		return std::string();
	}

	bool ConductorCache::hasCriticalPointTable() const
	{
		return !criticalFluxQ.empty() && !criticalFluxOffset.empty()
		       && !criticalPotentialOffset.empty();
	}

	bool ConductorCache::empty() const
	{
		return nodalPsi.empty() && quadraturePsi.empty()
		       && !hasCriticalPointTable();
	}

	bool ConductorCache::consistent() const
	{
		if ( signature.potentialDofs < 0 || signature.elements < 0 )
			return false;

		if ( !nodalPsi.empty()
		     && nodalPsi.size()
		        != static_cast< std::size_t >( signature.potentialDofs ) )
			return false;

		// The quadrature half is present as a whole or not at all: an offset
		// table without values indexes nothing, and values without offsets
		// cannot be indexed.
		if ( quadraturePsi.empty() != quadratureOffset.empty() )
			return false;

		if ( quadratureOffset.empty() )
			return criticalPointTableIsConsistent();

		if ( quadratureOffset.size()
		     != static_cast< std::size_t >( signature.elements ) + 1 )
			return false;

		if ( quadratureOffset.front() != 0 )
			return false;

		for ( std::size_t i = 1; i < quadratureOffset.size(); ++i )
			if ( quadratureOffset[ i ] < quadratureOffset[ i - 1 ] )
				return false;

		if ( static_cast< std::size_t >( quadratureOffset.back() )
		     != quadraturePsi.size() )
			return false;

		return criticalPointTableIsConsistent();
	}

	bool ConductorCache::criticalPointTableIsConsistent() const
	{
		// ABSENT IS FINE AND PARTIAL IS NOT. The offsets index the values and
		// the screen indexes the flux nodes, so any one of the four missing
		// leaves the others indexing something that is not there.
		bool const any = !criticalPotentialPsi.empty()
		                 || !criticalPotentialOffset.empty()
		                 || !criticalFluxQ.empty()
		                 || !criticalFluxUsable.empty()
		                 || !criticalFluxOffset.empty();
		if ( !any )
			return true;

		if ( criticalPotentialOffset.size()
		     != static_cast< std::size_t >( signature.elements ) + 1
		     || criticalFluxOffset.size()
		        != static_cast< std::size_t >( signature.elements ) + 1 )
			return false;

		for ( std::vector< int > const *offsets :
		      { &criticalPotentialOffset, &criticalFluxOffset } )
		{
			if ( offsets->front() != 0 )
				return false;
			for ( std::size_t i = 1; i < offsets->size(); ++i )
				if ( ( *offsets )[ i ] < ( *offsets )[ i - 1 ] )
					return false;
		}

		if ( static_cast< std::size_t >( criticalPotentialOffset.back() )
		     != criticalPotentialPsi.size() )
			return false;

		std::size_t const fluxNodes
			= static_cast< std::size_t >( criticalFluxOffset.back() );

		// TWO DOUBLES PER NODE, ONE BYTE PER NODE. Getting this pair wrong is
		// the one way a reader could walk off the end of q_c while every
		// offset still looked sane.
		return criticalFluxQ.size() == 2*fluxNodes
		       && criticalFluxUsable.size() == fluxNodes;
	}

	std::uint64_t digestSeed()
	{
		return fnvOffsetBasis;
	}

	std::uint64_t digestAppend( std::uint64_t digest, double value )
	{
		std::uint64_t bits = 0;
		std::memcpy( &bits, &value, sizeof( bits ) );
		return mixBytes( digest, bits );
	}

	std::uint64_t digestAppend( std::uint64_t digest, std::int64_t value )
	{
		return mixBytes( digest, static_cast< std::uint64_t >( value ) );
	}

	char const *conductorCacheGroupName()
	{
		return "psi_coil";
	}

	void writeConductorCache( std::string const &path,
	                          ConductorCache const &cache,
	                          std::string const &group )
	{
		if ( !cache.consistent() )
			throw std::invalid_argument(
				"meq::writeConductorCache: the cache is not self-consistent, "
				"so it would be refused on the way back in. Check that the "
				"offsets are a monotone prefix sum of length elements + 1 and "
				"that the nodal count is the dof count." );

		try
		{
			netCDF::NcFile file;

			// An empty group name means a file of its own, so it is created;
			// a named group is added to an equilibrium that already exists,
			// which is why the modes differ.
			if ( group.empty() )
				file.open( path, netCDF::NcFile::replace );
			else
				file.open( path, netCDF::NcFile::write );

			netCDF::NcGroup target
				= group.empty() ? static_cast< netCDF::NcGroup >( file )
				                : file.addGroup( group );

			target.putAtt( "format", formatName );
			target.putAtt( "format_version", netCDF::ncInt, formatVersion );
			target.putAtt( "mesh_digest", toHex( cache.signature.meshDigest ) );
			target.putAtt( "conductor_digest",
			               toHex( cache.signature.conductorDigest ) );
			target.putAtt( "polynomial_degree", netCDF::ncInt,
			               cache.signature.polynomialDegree );
			target.putAtt( "elements", netCDF::ncInt, cache.signature.elements );
			target.putAtt( "potential_dofs", netCDF::ncInt,
			               cache.signature.potentialDofs );
			target.putAtt( "note",
			               "psi_c at exactly the points MEQ evaluates it at. "
			               "Not a grid and not interpolable: see "
			               "src/meq/ConductorStore.hpp" );

			if ( !cache.nodalPsi.empty() )
			{
				netCDF::NcDim dofDim
					= target.addDim( "potential_dof", cache.nodalPsi.size() );
				netCDF::NcVar var
					= target.addVar( "psi_c_dof", netCDF::ncDouble, dofDim );
				var.putAtt( "long_name",
				            "psi_c at every potential degree of freedom" );
				var.putAtt( "units", "Wb/rad" );
				var.putVar( cache.nodalPsi.data() );
			}

			if ( !cache.quadraturePsi.empty() )
			{
				netCDF::NcDim pointDim = target.addDim(
					"quadrature_point", cache.quadraturePsi.size() );
				netCDF::NcDim offsetDim = target.addDim(
					"element_offset", cache.quadratureOffset.size() );

				netCDF::NcVar var = target.addVar(
					"psi_c_quadrature", netCDF::ncDouble, pointDim );
				var.putAtt( "long_name",
				            "psi_c at every source quadrature point, flattened" );
				var.putAtt( "units", "Wb/rad" );
				var.putVar( cache.quadraturePsi.data() );

				netCDF::NcVar offsets = target.addVar(
					"quadrature_offset", netCDF::ncInt, offsetDim );
				offsets.putAtt( "long_name",
				                "first psi_c_quadrature entry of each element" );
				offsets.putVar( cache.quadratureOffset.data() );
			}

			if ( cache.hasCriticalPointTable() )
			{
				// ONE SHARED element_node_offset DIMENSION would be wrong: the
				// two spaces have their own node counts and MEQ does not
				// require them to agree, so each offset table gets its own.
				netCDF::NcDim elementDim = target.addDim(
					"critical_element_offset",
					cache.criticalPotentialOffset.size() );

				netCDF::NcVar potentialOffsets = target.addVar(
					"critical_potential_offset", netCDF::ncInt, elementDim );
				potentialOffsets.putVar( cache.criticalPotentialOffset.data() );

				netCDF::NcVar fluxOffsets = target.addVar(
					"critical_flux_offset", netCDF::ncInt, elementDim );
				fluxOffsets.putVar( cache.criticalFluxOffset.data() );

				if ( !cache.criticalPotentialPsi.empty() )
				{
					netCDF::NcDim nodeDim = target.addDim(
						"critical_potential_node",
						cache.criticalPotentialPsi.size() );
					netCDF::NcVar var = target.addVar(
						"critical_psi_c", netCDF::ncDouble, nodeDim );
					var.putAtt( "long_name",
					            "psi_c at each potential-space element node" );
					var.putAtt( "units", "Wb/rad" );
					var.putVar( cache.criticalPotentialPsi.data() );
				}

				netCDF::NcDim fluxNodeDim = target.addDim(
					"critical_flux_node", cache.criticalFluxUsable.size() );
				netCDF::NcDim componentDim
					= target.getDim( "component" ).isNull()
					  ? target.addDim( "component", 2 )
					  : target.getDim( "component" );

				std::vector< netCDF::NcDim > const shape{ fluxNodeDim,
				                                          componentDim };
				netCDF::NcVar flux = target.addVar( "critical_q_c",
				                                    netCDF::ncDouble, shape );
				flux.putAtt( "long_name",
				             "q_c = ( 1/R ) grad_bar psi_c at each flux-space "
				             "element node, ( qR, qZ )" );
				flux.putVar( cache.criticalFluxQ.data() );

				netCDF::NcVar usable = target.addVar(
					"critical_q_c_usable", netCDF::ncByte, fluxNodeDim );
				usable.putAtt( "long_name",
				               "0 where q_c is NaN, which is the axis" );
				usable.putVar( cache.criticalFluxUsable.data() );

				target.putAtt( "critical_flux_scale", netCDF::ncDouble,
				               cache.criticalFluxScale );
			}
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error(
				"meq::writeConductorCache: could not write " + path + ": "
				+ error.what() );
		}
	}

	std::optional< ConductorCache >
	readConductorCache( std::string const &path, std::string const &group )
	{
		try
		{
			netCDF::NcFile file( path, netCDF::NcFile::read );

			netCDF::NcGroup target;
			if ( !group.empty() )
			{
				target = file.getGroup( group );
				if ( target.isNull() )
					return std::nullopt;
			}
			else if ( carriesTheFormat( file ) )
			{
				target = file;
			}
			else
			{
				target = file.getGroup( conductorCacheGroupName() );
				if ( target.isNull() )
					return std::nullopt;
			}

			if ( !carriesTheFormat( target ) )
				return std::nullopt;

			int const version = intAttribute( target, "format_version", path );
			if ( version != formatVersion )
				throw std::runtime_error(
					"meq::readConductorCache: " + path + " is format version "
					+ std::to_string( version ) + " and this MEQ writes "
					+ std::to_string( formatVersion )
					+ ". Delete it and let the run rebuild the cache." );

			ConductorCache cache;
			cache.signature.meshDigest
				= fromHex( textAttribute( target, "mesh_digest", path ), path );
			cache.signature.conductorDigest = fromHex(
				textAttribute( target, "conductor_digest", path ), path );
			cache.signature.polynomialDegree
				= intAttribute( target, "polynomial_degree", path );
			cache.signature.elements = intAttribute( target, "elements", path );
			cache.signature.potentialDofs
				= intAttribute( target, "potential_dofs", path );

			netCDF::NcVar nodal = target.getVar( "psi_c_dof" );
			if ( !nodal.isNull() )
			{
				cache.nodalPsi.resize( nodal.getDim( 0 ).getSize() );
				nodal.getVar( cache.nodalPsi.data() );
			}

			netCDF::NcVar quadrature = target.getVar( "psi_c_quadrature" );
			netCDF::NcVar offsets = target.getVar( "quadrature_offset" );
			if ( !quadrature.isNull() && !offsets.isNull() )
			{
				cache.quadraturePsi.resize( quadrature.getDim( 0 ).getSize() );
				quadrature.getVar( cache.quadraturePsi.data() );
				cache.quadratureOffset.resize( offsets.getDim( 0 ).getSize() );
				offsets.getVar( cache.quadratureOffset.data() );
			}

			netCDF::NcVar potentialOffsets
				= target.getVar( "critical_potential_offset" );
			netCDF::NcVar fluxOffsets = target.getVar( "critical_flux_offset" );
			netCDF::NcVar fluxValues = target.getVar( "critical_q_c" );
			netCDF::NcVar fluxScreen = target.getVar( "critical_q_c_usable" );
			if ( !potentialOffsets.isNull() && !fluxOffsets.isNull()
			     && !fluxValues.isNull() && !fluxScreen.isNull() )
			{
				cache.criticalPotentialOffset.resize(
					potentialOffsets.getDim( 0 ).getSize() );
				potentialOffsets.getVar( cache.criticalPotentialOffset.data() );

				cache.criticalFluxOffset.resize(
					fluxOffsets.getDim( 0 ).getSize() );
				fluxOffsets.getVar( cache.criticalFluxOffset.data() );

				netCDF::NcVar potentialValues = target.getVar( "critical_psi_c" );
				if ( !potentialValues.isNull() )
				{
					cache.criticalPotentialPsi.resize(
						potentialValues.getDim( 0 ).getSize() );
					potentialValues.getVar( cache.criticalPotentialPsi.data() );
				}

				std::size_t const fluxNodes
					= fluxValues.getDim( 0 ).getSize();
				cache.criticalFluxQ.resize( 2*fluxNodes );
				fluxValues.getVar( cache.criticalFluxQ.data() );

				cache.criticalFluxUsable.resize(
					fluxScreen.getDim( 0 ).getSize() );
				fluxScreen.getVar( cache.criticalFluxUsable.data() );

				netCDF::NcGroupAtt scale
					= target.getAtt( "critical_flux_scale" );
				if ( !scale.isNull() )
					scale.getValues( &cache.criticalFluxScale );
			}

			// A FILE CLAIMING TO BE ONE OF THESE AND FAILING ITS OWN INVARIANTS
			// IS AN ERROR RATHER THAN AN ABSENCE. Returning nullopt would send
			// the caller down the "no cache here, compute one" path and hide a
			// corrupt file behind a slow run.
			if ( !cache.consistent() )
				throw std::runtime_error(
					"meq::readConductorCache: " + path + " carries a "
					+ formatName + " whose arrays do not agree with its own "
					"attributes -- truncated or written by something else" );

			return cache;
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error(
				"meq::readConductorCache: could not read " + path + ": "
				+ error.what() );
		}
	}
}
