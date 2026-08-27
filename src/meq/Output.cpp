#include "Output.hpp"

#include <fstream>
#include <limits>
#include <stdexcept>

#ifdef MEQ_USE_NETCDF
#include <netcdf>
#endif

namespace meq
{
	bool hasNetCDF()
	{
#ifdef MEQ_USE_NETCDF
		return true;
#else
		return false;
#endif
	}

	void writeMfem( std::string const &stem, mfem::Mesh &mesh,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux )
	{
		// Precision 16, not MFEM's default 8. These files are read back for an
		// exact restart, and eight digits is not exact -- it would put a 1e-8
		// perturbation into a warm start whose whole point is to be the previous
		// answer.
		auto open = [ &stem ]( char const *suffix )
		{
			std::ofstream out( stem + suffix );
			if ( !out )
				throw std::runtime_error( "meq::writeMfem: cannot write " + stem + suffix );
			out.precision( 16 );
			return out;
		};

		{ std::ofstream out = open( ".mesh" );          mesh.Print( out ); }
		{ std::ofstream out = open( "_psi.gf" );        potential.Save( out ); }
		{ std::ofstream out = open( "_grad_psi.gf" );   flux.Save( out ); }
	}

#ifdef MEQ_USE_NETCDF

	struct NetCDFWriter::State
	{
		netCDF::NcFile file;
		netCDF::NcDim rDim;
		netCDF::NcDim zDim;
		int nR;
		int nZ;
		bool closed;
	};

	NetCDFWriter::NetCDFWriter( std::string const &path, GridSampler const &sampler )
		: state( new State{ {}, {}, {}, sampler.nodesR(), sampler.nodesZ(), false } )
	{
		try
		{
			state->file.open( path, netCDF::NcFile::replace );
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error( "meq::NetCDFWriter: cannot create " + path
			                          + ": " + error.what() );
		}

		state->rDim = state->file.addDim( "R", static_cast<std::size_t>( state->nR ) );
		state->zDim = state->file.addDim( "Z", static_cast<std::size_t>( state->nZ ) );

		// The coordinates first, so the file is self-describing even if a later
		// write fails.
		std::vector<double> coordinate;

		coordinate.resize( static_cast<std::size_t>( state->nR ) );
		for ( int i = 0; i < state->nR; ++i ) coordinate[ i ] = sampler.rAt( i );
		netCDF::NcVar rVar = state->file.addVar( "R", netCDF::ncDouble, state->rDim );
		rVar.putAtt( "long_name", "Major radius" );
		rVar.putAtt( "units", "m" );
		rVar.putVar( coordinate.data() );

		coordinate.resize( static_cast<std::size_t>( state->nZ ) );
		for ( int j = 0; j < state->nZ; ++j ) coordinate[ j ] = sampler.zAt( j );
		netCDF::NcVar zVar = state->file.addVar( "Z", netCDF::ncDouble, state->zDim );
		zVar.putAtt( "long_name", "Height" );
		zVar.putAtt( "units", "m" );
		zVar.putVar( coordinate.data() );

		// The mask, which is the thing a reader can always trust. Written as
		// signed char because NetCDF's byte is one, and 0/1 rather than a
		// boolean because the format has no boolean.
		std::vector<signed char> mask( static_cast<std::size_t>( state->nR )*state->nZ, 0 );
		for ( int j = 0; j < state->nZ; ++j )
			for ( int i = 0; i < state->nR; ++i )
				mask[ static_cast<std::size_t>( j )*state->nR + i ] =
					sampler.located( i, j ) ? 1 : 0;

		std::vector<netCDF::NcDim> const shape{ state->zDim, state->rDim };
		netCDF::NcVar maskVar = state->file.addVar( "inside", netCDF::ncByte, shape );
		maskVar.putAtt( "long_name",
		                "1 where the node lies inside the computational domain" );
		maskVar.putVar( mask.data() );
	}

	NetCDFWriter::~NetCDFWriter()
	{
		// Swallowed on purpose: a throw from a destructor terminates, and the
		// caller who wanted the error had close() available.
		try { close(); } catch ( ... ) { }
	}

	void NetCDFWriter::close()
	{
		if ( state && !state->closed )
		{
			state->file.close();
			state->closed = true;
		}
	}

	void NetCDFWriter::attribute( std::string const &name, std::string const &value )
	{
		state->file.putAtt( name, value );
	}

	void NetCDFWriter::attribute( std::string const &name, double value )
	{
		state->file.putAtt( name, netCDF::ncDouble, value );
	}

	void NetCDFWriter::attribute( std::string const &name, int value )
	{
		state->file.putAtt( name, netCDF::ncInt, value );
	}

	void NetCDFWriter::field( std::string const &name,
	                          std::vector<double> const &values,
	                          std::string const &longName,
	                          std::string const &units )
	{
		std::size_t const expected = static_cast<std::size_t>( state->nR )*state->nZ;
		if ( values.size() != expected )
			throw std::runtime_error( "meq::NetCDFWriter::field: " + name + " has "
			                          + std::to_string( values.size() )
			                          + " values where the grid has "
			                          + std::to_string( expected ) );

		std::vector<netCDF::NcDim> const shape{ state->zDim, state->rDim };
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncDouble, shape );
		var.putAtt( "long_name", longName );
		var.putAtt( "units", units );
		// Both the fill attribute and the `inside` mask describe the same
		// absence. See the header for why both.
		var.putAtt( "_FillValue", netCDF::ncDouble,
		            std::numeric_limits<double>::quiet_NaN() );
		var.putVar( values.data() );
	}

	void NetCDFWriter::boundary( std::vector<double> const &r,
	                             std::vector<double> const &z )
	{
		if ( r.size() != z.size() )
			throw std::runtime_error( "meq::NetCDFWriter::boundary: the two coordinate arrays differ in length" );
		if ( r.empty() )
			return;

		netCDF::NcDim dim = state->file.addDim( "boundary", r.size() );
		netCDF::NcVar rVar = state->file.addVar( "boundary_R", netCDF::ncDouble, dim );
		rVar.putAtt( "long_name", "Major radius of the prescribed boundary" );
		rVar.putAtt( "units", "m" );
		rVar.putVar( r.data() );

		netCDF::NcVar zVar = state->file.addVar( "boundary_Z", netCDF::ncDouble, dim );
		zVar.putAtt( "long_name", "Height of the prescribed boundary" );
		zVar.putAtt( "units", "m" );
		zVar.putVar( z.data() );
	}

#else   // MEQ_USE_NETCDF

	struct NetCDFWriter::State { };

	namespace
	{
		[[noreturn]] void unavailable()
		{
			throw std::runtime_error( "meq::NetCDFWriter: meq was built without netcdf-cxx4" );
		}
	}

	NetCDFWriter::NetCDFWriter( std::string const &, GridSampler const & ) { unavailable(); }
	NetCDFWriter::~NetCDFWriter() = default;
	void NetCDFWriter::close() { }
	void NetCDFWriter::attribute( std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, double ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, int ) { unavailable(); }
	void NetCDFWriter::field( std::string const &, std::vector<double> const &,
	                          std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::boundary( std::vector<double> const &,
	                             std::vector<double> const & ) { unavailable(); }

#endif  // MEQ_USE_NETCDF
}
