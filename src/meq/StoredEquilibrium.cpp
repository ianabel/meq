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

		region = column( file, "element_region" );
		plasmaFraction = column( file, "element_plasma_fraction" );
		conductorMaterial = column( file, "element_conductor_material" );
		conductorFraction = column( file, "element_conductor_fraction" );
	}
#else
	StoredEquilibrium::StoredEquilibrium( std::string const &path )
	{
		throw std::runtime_error( "cannot read " + path
		                          + ": MEQ was built without netcdf-cxx4" );
	}
#endif
}
