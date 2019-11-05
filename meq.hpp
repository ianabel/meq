#ifndef MEQ_HPP
#define MEQ_HPP
/*
 * Main include file for MEQ
 */

#include <vector>
#include <memory>
#include <tuple>

#include "mfem.hpp"
#include "SplineInterpolant.hpp"
#include "GSSolver.hpp"


namespace meq {

	using RealScalarField = std::function<double( const mfem::Vector & )>;
	using RealVectorField = std::function<void( const mfem::Vector &, mfem::Vector & )>;

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

class Solution {
	private:
		mfem::Vector qu;
		mfem::Vector u_star;
		GSSolver &solver;
		mfem::GridFunction q_variable,u_variable,u_hat_variable;
		static const int dim = 2;
	public:
		Solution( GSSolver & solver_ref )
			: solver( solver_ref )
		{
			qu.SetSize( solver.Height() );

			q_variable.MakeRef( solver.QSpace().get(), qu, 0 );
			u_variable.MakeRef( solver.USpace().get(), qu, solver.QSpace()->GetVSize() );
			u_hat_variable.MakeRef( solver.MSpace().get(), qu, solver.QSpace()->GetVSize() + solver.USpace()->GetVSize() );

		}
		std::pair<double,double> l2_errors( RealScalarField uFun_ex, RealVectorField qFun_ex, int order ) {
			int order_quad = std::max(2, 2*order+4);
			const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
			for (int i=0; i < mfem::Geometry::NumGeom; ++i)
			{
				irs[i] = &(mfem::IntRules.Get(i, order_quad));
			}
			mfem::StdFunctionCoefficient ucoeff(uFun_ex);
			mfem::VectorStdFunctionCoefficient qcoeff(dim, qFun_ex);

			double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
			double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
			return std::make_pair( err_u, err_q );
		};

		std::pair<double,double> l2_errors( mfem::GridFunction *uFn, mfem::GridFunction *qFn, int order )
		{
			mfem::GridFunctionCoefficient ucoeff( uFn ),qcoeff( qFn );
int order_quad = std::max(2, 2*order+4);
			const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
			for (int i=0; i < mfem::Geometry::NumGeom; ++i)
			{
				irs[i] = &(mfem::IntRules.Get(i, order_quad));
			}

			double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
			double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
			return std::make_pair( err_u, err_q );
		};

		std::pair<double,double> lInf_errors( RealScalarField uFun_ex, RealVectorField qFun_ex, int order ) {
			int order_quad = std::max(2, 2*order+2);
			const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
			for (int i=0; i < mfem::Geometry::NumGeom; ++i)
			{
				irs[i] = &(mfem::IntRules.Get(i, order_quad));
			}
			mfem::StdFunctionCoefficient ucoeff(uFun_ex);
			mfem::VectorStdFunctionCoefficient qcoeff(dim, qFun_ex);

			double err_u    = u_variable.ComputeMaxError(ucoeff, irs);
			double err_q    = q_variable.ComputeMaxError(qcoeff, irs);
			return std::make_pair( err_u, err_q );
		};

		std::pair<double,double> lInf_errors( mfem::GridFunction *uFn, mfem::GridFunction *qFn, int order )
		{
			mfem::GridFunctionCoefficient ucoeff( uFn ),qcoeff( qFn );
int order_quad = std::max(2, 2*order+4);
			const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
			for (int i=0; i < mfem::Geometry::NumGeom; ++i)
			{
				irs[i] = &(mfem::IntRules.Get(i, order_quad));
			}

			double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
			double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
			return std::make_pair( err_u, err_q );
		};


};

}
#endif // MEQ_HPP
