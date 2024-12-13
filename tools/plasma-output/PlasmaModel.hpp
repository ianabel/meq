#ifndef PLASMAMODEL_HPP
#define PLASMAMODEL_HPP

#include <vector>

class PlasmaModel {
	public:
		using Index = std::size_t;
		using RealFunc = std::function<double( double )>;
	private:
		PlasmaModel();
		~PlasmaModel();
		double Psi( double R, double Z );
		double phi0( double R, double Z );
		double N( Index s, double psi );
	public:
		double Density( Index species, double R,  double z );
		double Temperature( Index species, double R,  double z );
		double Omega( double R, double z );

		struct Species {
			RealFunc Temperature; // T( Psi )
			RealFunc MidplaneDensity; // n( Psi ) on midplane
			double Z; // Charge in units of the elemental charge (electrons have Z = -1)
			double mass; 
			Species( double T_in, double Dens_in, double Z_in, double mass_in ) :
				Temperature( T_in ), Density( Dens_in ), Z( Z_in ), mass( mass_in )
			{
			}
			// Copy constructor
			Species( Species const& other ) :
				Temperature( other.Temperature ), Density( other.Density ), Z( other.Z ), mass( other.mass ), fprim( other.fprim ) 
			{
			};
		};

	private:
		std::size_t N_spec;
		std::vector<Species> Content;
		RealFunc Omega_psi,RMid_psi;
		Index IonIdx,ElectronIdx;

}

#endif // PLASMAMODEL_HPP
