#ifndef NETCDFIO_HPP
#define NETCDFIO_HPP

#include <netcdf>

/*
 * Class for storing plasma transport output in a NetCDF file.
 */

template< typename ProfileType > class NetCDFPlasma
{
	public:
		NetCDFPlasma( const std::string &filename, std::vector<Species> PlasmaSpecies, size_t Ncells, std::vector<double> cellBoundaries, size_t polynomialOrder );
		~NetCDFPlasma();

		WriteTimeslice( double T, std::vector< ProfileType& > const& Densities, std::vector< ProfileType& > const& Temperatures, 
		                 ProfileType & omega );

		struct Species {
			double Z; // Charge in units of the elemental charge (electrons have Z = -1)
			double Mass; // in Kilograms
		};


	private:
		std::vector<double> psiPoints;
		size_t N_spec;
		std::string filename;
		netCDF::NcFile data_file;
		netCDF::NcDim TimeDim,SpeciesDim,PsiDim;
		netCDF::NcVar TimeVar,SpeciesVar,PsiVar;

		netCDF::NcVar Density,Temperature,Omega;

		netCDF::NcCompoundType SpeciesNCType;

		void WriteHeader();
};


#endif // NETCDFIO_HPP
