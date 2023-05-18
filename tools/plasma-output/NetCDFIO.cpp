
#include "NetCDFIO.hpp"
#include "MirrorPlasma.hpp"

// Code for NetCDF interface
//

using namespace netCDF;
#include <cstddef>

template <typename ProfileType> NetCDFPlasma<ProfileType>::NetCDFPlasma( const std::string &filename, std::vector<Species> PlasmaSpecies, size_t Ncells, std::vector<double> cellBoundaries, size_t polynomialOrder );
	: psiPoints( Ncells * ( polynomialOrder + 1 ) )
{
	filename = file;
	data_file.open( file, netCDF::NcFile::FileMode::replace );
	TimeDim = data_file.addDim( "t" );
	TimeVar = data_file.addVar( "t", netCDF::NcDouble(), TimeDim );
	TimeVar.putAtt( "description", "Time since start of simulation" );
	TimeVar.putAtt( "units", "s" );

	SpeciesNCType = data_file.addCompoundType( "Species", sizeof( NetCDFPlasma<ProfileType>::Species ) );
	SpeciesNCType.addMember( "Z", netCDF::NcInt(), offsetof( NetCDFPlasma<ProfileType>::Species, Z ) );
	SpeciesNCType.addMember( "Mass", netCDF::NcDouble(), offsetof( NetCDFPlasma<ProfileType>::Species, Mass ) );

	N_spec = PlasmaSpecies.size();
	SpeciesDim = data_file.addDim( "Species", N_spec );
	SpeciesVar = data_file.addVar( "Species", SpeciesNCType, SpeciesDim );
	SpeciesVar.putAtt( "description", "List of species in the plasma" );

	// Because std::vector is guaranteed contiguous in memory we can just blat into the NetCDF variable now
	// as the whole point of NetCDF compound types is that they have the same memory layout as the corresponding C struct

	SpeciesVar.putVar( PlasmaSpecies.data() );

	PsiDim = data_file.addDim( "Psi", psiPoints.size() );
	PsiVar = data_file.addVar( "Psi", netCDF::NcDouble(), PsiDim );

	for ( size_t i = 0; i < cellBoundaries.size()-1; ++i )
	{
		double x_l = cellBoundaries[ i ];
		double x_u = cellBoundaries[ i + 1 ];
		double h = x_u - x_l;
		for ( unsigned int j=0; j <= polynomialOrder; ++j )
		{
			psiPoints[ i*( polynomialOrder + 1 ) + j ] = x_l + ( h*j ) /  ( polynomialOrder + 1.0 );
		}
	}

	PsiVar.putAtt( "description", "Axial magnetic flux inside a flux surface" );
	PsiVar.putAtt( "units", "Wb" );
	PsiVar.putVar( psiPoints.data() );

	Density = data_file.addVar( "Density", netCDF::NcDouble(), { SpeciesDim, TimeDim, PsiDim } );
	Temperature = data_file.addVar( "Temperature", netCDF::NcDouble(), { SpeciesDim, TimeDim, PsiDim } );
	Omega = data_file.addVar( "Omega", netCDF::NcDouble(), { TimeDim, PsiDim } );
}

template <typename ProfileType> void NetCDFPlasma<ProfileType>::WriteTimeslice( double T, std::vector< ProfileType& > const& Densities, std::vector< ProfileType& > const& Temperatures, ProfileType & omega )
{
	size_t next = TimeDim.getSize();
	TimeVar.putVar( {next}, T );
	std::vector<double> data( N_psi ), data2( N_psi );

	for ( size_t i=0; i < N_psi; ++i ) {
		data[ i ] = omega( psiPoints[ i ] );
	}

	// std::vector<size_t> start = {next, 0}, count = {1, N_psi};

	Omega.putVar( {next, 0}, {1, N_psi}, data.data() );

	for ( size_t s=0; s < N_spec; ++s ) {
		for ( size_t i=0; i < N_psi; ++i ) {
			data[ i ]  = Densities[ s ]( psiPoints[ i ] );
			data2[ i ] = Temperatures[ s ]( psiPoints[ i ] );
		}
		Density.putVar( {s, next, 0}, {1, 1, N_psi}, data.data() );
		Temperature.putVar( {s, next, 0}, {1, 1, N_psi}, data2.data() )
	}
}

template <typename ProfileType> void NetCDFPlasma<ProfileType>::Close()
{
	filename = "";
	data_file.close();
}

template <typename ProfileType> NetCDFPlasma<ProfileType>::~NetCDFPlasma()
{
	if ( filename != "" )
		Close();
}


