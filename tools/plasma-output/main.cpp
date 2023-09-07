

#include <cmath>
#include <cstdint>
#include <vector>
#include <netcdf>
#include <functional>
#include <boost/math/tools/roots.hpp>

using RealFunc = std::function<double( double )>;

struct Species {
	RealFunc Temperature; // T( Psi )
	RealFunc MidplaneDensity; // n( Psi ) on midplane
	double Z; // Charge in units of the elemental charge (electrons have Z = -1)
	double mass; 
	Species( RealFunc T, RealFunc Dens, double Z_, double mass_ ) :
		Temperature( T ), MidplaneDensity( Dens ), Z( Z_ ), mass( mass_ )
	{
	}
};

std::size_t IonIndex,ElectronIndex;
std::vector<Species> PlasmaContent;
std::function<double( double )> Omega;

double Density( Species const& s, double psi, double R, double z, double phi )
{
	double T =  s.Temperature( psi );
	double Z = s.Z;
	double m = s.mass;
	double omega = Omega( psi );

	return s.MidplaneDensity( psi ) * std::exp( -Z * phi / T + m*omega*omega*R*R / ( 2.0 * T )  );

}


struct cd_functor {
	double R,Z,psi;
	cd_functor( double R_, double Z_, double psi_ ) : R( R_ ), Z( Z_ ), psi( psi_ ) {};

	double charge_density( double phi ) {
		double rho_c = 0;
		for ( auto const& s : PlasmaContent ) {
			rho_c += s.Z * Density( s, psi, R, Z, phi );
		}
		return rho_c;
	};
	double charge_density_prime( double phi ) {
		double rho_c = 0;
		for ( auto const& s : PlasmaContent ) {
			rho_c -= s.Z * ( s.Z / s.Temperature( psi ) ) * Density( s, psi, R, Z, phi );
		}
		return rho_c;
	};
	std::pair<double,double> operator()( double phi ) {
		return std::make_pair( charge_density( phi ), charge_density_prime( phi ) );
	};
};

double phi0( double psi, double R, double z )
{
	cd_functor cdf( R, z, psi );
	int digits = 6;
	std::uintmax_t it = 20;
	Species const& s = PlasmaContent[ IonIndex ];
	double T = s.Temperature( psi );
	double Z = s.Z;
	double m = s.mass;
	double omega = Omega( psi );
	double guess = ( m * omega * omega * R * R )/( 2.0 * T ); // Definitely wrong, but should be good enough for Newton

	double max = guess * 4;
	double ans = boost::math::tools::newton_raphson_iterate( cdf, guess, -max, max, digits, it );
	return ans;
}

double Density( Species const& s, double psi, double R, double z )
{
	double phi = phi0( psi, R, z );
	return Density( s, psi, R, z, phi );
}

int main( int, char ** )
{
	double *psiData;
	netCDF::NcFile data_file;
	std::string netcdf_file = "exact-vacuum.nc";
	data_file.open( netcdf_file, netCDF::NcFile::FileMode::read );
	netCDF::NcDim R_dim,Z_dim;
	R_dim = data_file.getDim( "R" );
	Z_dim = data_file.getDim( "Z" );

	size_t N_R = R_dim.getSize(), N_Z = Z_dim.getSize();

	std::vector<double> R_data,Z_data;

	R_data.resize( N_R );
	Z_data.resize( N_Z );

	netCDF::NcVar R_var,Z_var;

	R_var = data_file.getVar( "R" );
	Z_var = data_file.getVar( "Z" );

	R_var.getVar( R_data.data() );
	Z_var.getVar( Z_data.data() );

	psiData = new double[ N_R * N_Z ];

	netCDF::NcVar psi = data_file.getVar( "psi" );
	psi.getVar( psiData );

	// Given N / T / omega as functions of psi, project
	//
	
	constexpr double ElectronMass = 9.1094e-31;	//Electron Mass, kg
	constexpr double IonMass = 1.6726e-27;	//Ion Mass ( = proton mass) kg
	constexpr double e_charge = 1.60217663e-19; //Coulombs
	
	auto ElectronTemeprature = [=]( double Psi ) {
		return 100 * e_charge; // Isothermal with 100 eV
	};

	double psi_l,psi_u; // bounds in psi for the plasma
	double nMax = 1e20;

	auto ElectronMidplaneDensity = [ = ]( double psi ) {
		double psiL = psi_u - psi_l;
		return nMax*( psi - psi_l )*( psi_u - psi )/( psiL * psiL * 0.25 );
	};

	auto IonMidplaneDensity = [ = ]( double psi ) {
		return ElectronMidplaneDensity( psi );
	};

	auto IonTemperature = [=]( double Psi ) {
		return 150 * e_charge; // Isothermal with 150 eV
	};

	PlasmaContent.clear();
	PlasmaContent.emplace_back( ElectronTemeprature, ElectronMidplaneDensity, -1, ElectronMass );
	ElectronIndex = 0;
	PlasmaContent.emplace_back( IonTemperature, IonMidplaneDensity, -1, IonMass );
	IonIndex = 1;


	
	return 0;
}
