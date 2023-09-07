#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include <boost/math/interpolators/makima.hpp>
#include <vector>
#include <string>

class Equilibrium {
	public:
		Equilibrium( std::string const& netcdf_file );
		~Equilibrium() {
			if ( psiData != nullptr )
				delete psiData;
		};
		double Psi( double , double );
	private:
		std::vector<double> R_data,Z_data;
		double *psiData = nullptr;
};

#endif // EQUILIBRIUM_HPP

