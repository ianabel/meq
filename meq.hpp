#ifndef MEQ_HPP
#define MEQ_HPP
/*
 * Utility classes for MEQ
 */

#include "mfem.hpp"

class Coil {
	private:
		double R_l,R_u,Z_l,Z_u;
	public:
		double R,z;
		double J;
		double h,w;
		Coil( double R_0, double z_0, double height, double width, double Current )
			: R( R_0 ), z( z_0 ), h( height ), w( width ), J( Current )
		{
			R_l = R - w/2;
			R_u = R + w/2;
			Z_l = z - h/2;
			Z_u = z + h/2;
		};

		bool inline Contains( mfem::Vector const& pt ) const
		{
			return ( ( pt( 0 ) >= R_l ) && ( pt( 0 ) <= R_u ) && ( pt( 1 ) >= Z_l ) && ( pt( 1 ) <= Z_u ) );
		}

		~Coil() {};
};

class Jtor {
	protected:
		std::vector<Coil> Coils;
	public:
		Jtor() {Coils.clear();};
		~Jtor() {};
		void AddCoil( double R, double z, double h, double w, double J ) { Coils.emplace_back( R, z, h, w, J );};
		void AddCoil( Coil const& othercoil )
		{
			Coils.emplace_back( othercoil ); // Copy construct
		};
		void AddCoils( std::vector<Coil> const& othercoils )
		{
			Coils.insert( Coils.end(), othercoils.begin(), othercoils.end() );
		};
		void AddCoils( std::initializer_list<Coil> il )
		{
			Coils.insert( Coils.end(), il );
		};

		// Returns R*j_tor
		double operator()( mfem::Vector const& pt ) {
			double JtorPt = 0.0;
			for ( auto const &coil : Coils )
			{
				if ( coil.Contains( pt ) )
					JtorPt += coil.J;
			}
			return pt( 0 ) * JtorPt;
		};
};

class BoundaryCondition {
	private:
	BoundaryConditionType BCType;
	public:
	
}

enum BoundaryConditionType { VonHagenow, Prescribed, Zero };

class Domain {
	double RMin,RMax,ZMin,ZMax;
	double psiResolution;
	BoundaryConditionType BCType;
	BoundaryCondition* BCs;
};


// Store as a spline interpolant internally, 
// which  has constructors to read from .dat / NetCDF / etc.
using UserSuppliedFunction = SplineInterpolant;

class ToridalField {
	private:
		UserSuppliedFunction FFprime;
	public:
		ToroidalField();
		ToroidalField( std::string const & datafile, UserSuppliedFunction::DataFileType type ) 
			: FFPrime( datafile, type )
		{
		};
		~ToroidalField() {};
		double operator( double psi ) {

		}

}
class PlasmaModel {
	private:
		double PlasmaPsiMin,PlasmaPsiMax;
		bool ToroidalField;
		bool Rotation;
	public:
		double operator( const mfem::Vector & pt, double psi ) = 0;
};

class StaticMHDPlasma : KineticPlasmaModel {
	private:
		UserSuppliedFunction PPrime;
	public:
		double operator( const mfem::Vector& pt, double psi ) {
				
		}

}


#endif // MEQ_HPP
