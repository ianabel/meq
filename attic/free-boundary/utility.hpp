#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "mfem.hpp"

/*
 * Utility classes for MEQ
 */

namespace meq {

class Coil {
	private:
		double R_l,R_u,Z_l,Z_u;
	public:
		double R,z;
		double h,w;
		double J;
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

		// Returns j_tor
		double operator()( mfem::Vector const& pt ) {
			double JtorPt = 0.0;
			for ( auto const &coil : Coils )
			{
				if ( coil.Contains( pt ) )
					JtorPt += coil.J;
			}
			return JtorPt;
		};
};

enum BoundaryConditionType { VonHagenow, Prescribed, Zero, Unknown };

class BoundaryCondition {
	private:
		BoundaryConditionType BCType;
	public:
		virtual ~BoundaryCondition() {};
		virtual double operator()( const mfem::Vector &boundary_point ) = 0;

		bool CheckBCType( BoundaryConditionType Type ) { return BCType == Type;};
		BoundaryConditionType getBCType() { return BCType; };
};


class Domain {
	private:
		BoundaryConditionType BCType;
		BoundaryCondition* BCs;
	public:
		Domain( double RMin_in, double RMax_in, double ZMin_in, double ZMax_in, double resolution, BoundaryConditionType type ) 
			: BCType( type ), BCs( nullptr ), RMin( RMin_in ), RMax( RMax_in ), ZMin( ZMin_in ), ZMax( ZMax_in ),CellSize( resolution )
		{

		};
		double RMin,RMax,ZMin,ZMax;
		double CellSize;
};

};

#endif // UTILITY_HPP
