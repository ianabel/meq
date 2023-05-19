#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include <boost/math/interpolators/makima.hpp>

class Equilibrium {
	public:
		Equilibrium( std::string const& netcdf_file );
		~Equilibrium();
		double Psi( double , double );
	private:
		double *pData;
}

#endif // EQUILIBRIUM_HPP

