#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include <boost/math/interpolators/makima.hpp>

class Equilibrium {
	public:
		Equilibrium( std::string const& netcdf_file );
		~Equilibrium();
		double Psi( double , double );
	private:
		std::vector< boost::math::interpolators::makima<std::vector<double>> > psiInterpolantZ;
}

#endif // EQUILIBRIUM_HPP

