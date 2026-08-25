#ifndef MODEL_HPP
#define MODEL_HPP

/*
 * Classes for different plasma models
 */

#include "SplineInterpolant.hpp"

namespace meq {

// Store as a spline interpolant internally, 
// which  has constructors to read from .dat / NetCDF / etc.
using UserSuppliedFunction = SplineInterpolant;

class ToroidalField {
	private:
		UserSuppliedFunction FFPrime;
	public:
		ToroidalField( std::string const & datafile, UserSuppliedFunction::DataFileType type ) 
			: FFPrime( datafile, type )
		{
		};
		~ToroidalField() {};
		double operator()( double psi ) {
			return FFPrime( psi );
		}

};

class PlasmaModel {
	private:
		double PlasmaPsiMin,PlasmaPsiMax;
		bool ToroidalField;
		bool Rotation;
	public:
		virtual ~PlasmaModel() {};
		virtual double operator()( const mfem::Vector & pt, double psi ) = 0;
};


class StaticMHDPlasma : PlasmaModel {
	private:
		UserSuppliedFunction PPrime;
	public:
		StaticMHDPlasma( std::string const& datafile, UserSuppliedFunction::DataFileType type )
			: PPrime( datafile, type )
		{
		};

		double operator()( const mfem::Vector& pt, double psi ) {
			return pt( 0 )*pt( 0 )*PPrime( psi );
		};
};

};

#endif // MODEL_HPP

