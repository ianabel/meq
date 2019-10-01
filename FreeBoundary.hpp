#include <functional>
#include "mfem.hpp"

namespace meq {

double GreensFunction( mfem::Vector const& r, mfem::Vector const& r_star );
double SinhTanhQuad( double a, double b, std::function<double( double )> F, unsigned int N=20 ) ;
double BoundaryPsi( mfem::FiniteElementSpace *q_space, mfem::Vector & zero_soln, mfem::Vector const& r );
double GreensFunctionPsi( mfem::Mesh * mesh, mfem::Vector r, std::function<double( const mfem::Vector& )> const& j_coil );


class GreensFunctionBoundaryCoefficient : public mfem::Coefficient
{
	protected:
		mfem::Mesh const *mesh;
		mfem::FiniteElementSpace *Q_space;
		mfem::Vector & psi_hat;
	public:
		GreensFunctionBoundaryCoefficient( mfem::Mesh const* mesh_r, mfem::FiniteElementSpace *q_space, mfem::Vector & zero_soln )
			: mesh( mesh_r ), Q_space( q_space ), psi_hat( zero_soln )
		{
		}

		virtual double Eval( mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip ) override;
		
};


class FullGreensFunctionBoundaryCoefficient : public mfem::Coefficient
{
	protected:
		mfem::Mesh *mesh;
		std::function<double( const mfem::Vector & )> const& j_tor;

	public:
		FullGreensFunctionBoundaryCoefficient( mfem::Mesh * mesh_r, std::function<double( const mfem::Vector& )> const& j_phi )
			: mesh( mesh_r ),  j_tor( j_phi )
		{
		}

		virtual double Eval( mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip ) override;
		
		
};

}
