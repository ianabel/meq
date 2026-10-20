#include "StoredEquilibrium.hpp"

#include <sstream>
#include <stdexcept>

#ifdef MEQ_USE_NETCDF
#include <netcdf>
#endif

namespace meq
{
#ifdef MEQ_USE_NETCDF
	namespace
	{
		std::string textAttribute( netCDF::NcFile const &file, std::string const &name,
		                           std::string const &path )
		{
			netCDF::NcGroupAtt const attribute = file.getAtt( name );
			if ( attribute.isNull() )
				throw std::runtime_error( path + " has no " + name + " attribute; it was "
				                          "not written by this MEQ" );
			std::string value;
			attribute.getValues( value );
			return value;
		}

		int integerAttribute( netCDF::NcFile const &file, std::string const &name,
		                      std::string const &path )
		{
			netCDF::NcGroupAtt const attribute = file.getAtt( name );
			if ( attribute.isNull() )
				throw std::runtime_error( path + " has no " + name + " attribute; it was "
				                          "not written by this MEQ" );
			int value = 0;
			attribute.getValues( &value );
			return value;
		}

		double doubleAttribute( netCDF::NcFile const &file, std::string const &name,
		                        std::string const &path )
		{
			netCDF::NcGroupAtt const attribute = file.getAtt( name );
			if ( attribute.isNull() )
				throw std::runtime_error( path + " has no " + name + " attribute" );
			double value = 0.0;
			attribute.getValues( &value );
			return value;
		}

		/// A one-dimensional variable as doubles, or empty when it is absent.
		std::vector<double> column( netCDF::NcFile const &file, std::string const &name )
		{
			netCDF::NcVar const var = file.getVar( name );
			if ( var.isNull() )
				return {};
			std::size_t size = 1;
			for ( netCDF::NcDim const &dim : var.getDims() )
				size *= dim.getSize();
			std::vector<double> values( size );
			var.getVar( values.data() );
			return values;
		}

		/// Coefficients into a fresh GridFunction on @a space, refusing a count
		/// that does not fit it -- a file and a space that disagree are not a
		/// field.
		std::unique_ptr<mfem::GridFunction> field( mfem::FiniteElementSpace &space,
		                                           std::vector<double> const &values,
		                                           std::string const &name,
		                                           std::string const &path )
		{
			if ( values.size() != static_cast<std::size_t>( space.GetVSize() ) )
				throw std::runtime_error(
					path + ": " + name + " holds " + std::to_string( values.size() )
					+ " coefficients and its space " + std::to_string( space.GetVSize() ) );
			auto result = std::make_unique<mfem::GridFunction>( &space );
			mfem::real_t *data = result->HostWrite();
			std::copy( values.begin(), values.end(), data );
			return result;
		}
	}

	StoredEquilibrium::StoredEquilibrium( std::string const &path )
	{
		netCDF::NcFile file;
		try
		{
			file.open( path, netCDF::NcFile::read );
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error( "cannot read " + path + ": " + error.what() );
		}

		netCDF::NcVar const meshVar = file.getVar( "mesh" );
		if ( meshVar.isNull() )
			throw std::runtime_error(
				path + " carries no mesh: it was written before MEQ embedded the solve's "
				"mesh in its .nc. Restart from the .gf and .mesh that run wrote "
				"([initialguess] File and MeshFile), or re-run it" );
		{
			char *pointer = nullptr;
			meshVar.getVar( &pointer );
			text = pointer ? pointer : "";
			nc_free_string( 1, &pointer );
		}
		// generate_edges = 1, refine = 0, fix_orientation = false: see the header.
		std::istringstream meshStream( text );
		meshPointer = std::make_unique<mfem::Mesh>( meshStream, 1, 0, false );

		degree = integerAttribute( file, "polynomial_degree", path );
		axisGiven = !file.getAtt( "psi_axis" ).isNull();
		if ( axisGiven )
			axis = doubleAttribute( file, "psi_axis", path );
		contentText = textAttribute( file, "content", path );
		total = contentText.rfind( "total", 0 ) == 0;

		potentialCollection.reset( mfem::FiniteElementCollection::New(
			textAttribute( file, "fe_collection", path ).c_str() ) );
		potentialSpace = std::make_unique<mfem::FiniteElementSpace>(
			meshPointer.get(), potentialCollection.get(),
			integerAttribute( file, "fe_vdim", path ) );
		potentialField = field( *potentialSpace, column( file, "psi_coefficients" ),
		                        "psi_coefficients", path );

		fluxCollection.reset( mfem::FiniteElementCollection::New(
			textAttribute( file, "flux_collection", path ).c_str() ) );
		fluxSpace = std::make_unique<mfem::FiniteElementSpace>(
			meshPointer.get(), fluxCollection.get(),
			integerAttribute( file, "flux_vdim", path ),
			textAttribute( file, "flux_ordering", path ) == "byNODES"
				? mfem::Ordering::byNODES : mfem::Ordering::byVDIM );
		fluxField = field( *fluxSpace, column( file, "flux_coefficients" ),
		                   "flux_coefficients", path );

		postCollection.reset( mfem::FiniteElementCollection::New(
			textAttribute( file, "psistar_collection", path ).c_str() ) );
		postSpace = std::make_unique<mfem::FiniteElementSpace>(
			meshPointer.get(), postCollection.get() );
		postField = field( *postSpace, column( file, "psistar_coefficients" ),
		                   "psistar_coefficients", path );

		if ( !file.getAtt( "shape_minor_radius" ).isNull() )
			curve = std::make_unique<BoundaryShape>(
				doubleAttribute( file, "shape_major_radius", path ),
				doubleAttribute( file, "shape_centre_height", path ),
				doubleAttribute( file, "shape_minor_radius", path ),
				doubleAttribute( file, "shape_elongation", path ),
				column( file, "shape_cosines" ), column( file, "shape_sines" ) );

		if ( !file.getAtt( "xpoint_located" ).isNull()
		     && integerAttribute( file, "xpoint_located", path ) == 1 )
		{
			xPointGiven = true;
			xPointRadius = doubleAttribute( file, "xpoint_r", path );
			xPointHeight = doubleAttribute( file, "xpoint_z", path );
		}
		axisPointGiven = !file.getAtt( "axis_r" ).isNull() && !file.getAtt( "axis_z" ).isNull();
		if ( axisPointGiven )
		{
			axisPointRadius = doubleAttribute( file, "axis_r", path );
			axisPointHeight = doubleAttribute( file, "axis_z", path );
		}
		boundaryGiven = !file.getAtt( "psi_boundary" ).isNull();
		if ( boundaryGiven )
			boundary = doubleAttribute( file, "psi_boundary", path );
		scaleGiven = !file.getAtt( "profile_scale" ).isNull();
		if ( scaleGiven )
			scale = doubleAttribute( file, "profile_scale", path );
		if ( !file.getAtt( "inverse_currents" ).isNull() )
		{
			std::istringstream in( textAttribute( file, "inverse_currents", path ) );
			double value = 0.0;
			while ( in >> value )
				currents.push_back( value );
		}

		if ( !file.getAtt( "g_squared_coefficients" ).isNull() )
		{
			std::istringstream in( textAttribute( file, "g_squared_coefficients", path ) );
			double value = 0.0;
			while ( in >> value )
				gSquared.push_back( value );
		}

		trace = column( file, "trace_coefficients" );
		if ( !file.getAtt( "exterior_modes" ).isNull() )
		{
			int const modes = integerAttribute( file, "exterior_modes", path );
			for ( int m = 0; m < modes; ++m )
				exterior.push_back( doubleAttribute(
					file, "exterior_a" + std::to_string( m + 2 ), path ) );
		}

		region = column( file, "element_region" );
		plasmaFraction = column( file, "element_plasma_fraction" );
		conductorMaterial = column( file, "element_conductor_material" );
		conductorFraction = column( file, "element_conductor_fraction" );
	}
	void StoredEquilibrium::wholeState( mfem::Vector &state ) const
	{
		int const nFlux = fluxField->Size(), nPotential = potentialField->Size();
		int const nTrace = static_cast<int>( trace.size() );
		state.SetSize( nFlux + nPotential + nTrace );
		mfem::real_t *const to = state.HostWrite();
		mfem::real_t const *const q = fluxField->HostRead();
		mfem::real_t const *const psi = potentialField->HostRead();
		for ( int i = 0; i < nFlux; ++i )
			to[ i ] = -q[ i ];                                // DarcyForm's -q
		std::copy( psi, psi + nPotential, to + nFlux );
		std::copy( trace.begin(), trace.end(), to + nFlux + nPotential );
	}
#else
	StoredEquilibrium::StoredEquilibrium( std::string const &path )
	{
		throw std::runtime_error( "cannot read " + path
		                          + ": MEQ was built without netcdf-cxx4" );
	}
	void StoredEquilibrium::wholeState( mfem::Vector & ) const { }
#endif
}
