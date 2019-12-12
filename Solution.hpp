#ifndef SOLUTION_HPP
#define SOLUTION_HPP

/*
 * Solution class, wrapping a full solution ot the Grad-Shafronv equation
 */

#include <fstream>
#include <iostream>

#include "config.hpp"

namespace meq {

	using RealScalarField = std::function<double( const mfem::Vector & )>;
	using RealVectorField = std::function<void( const mfem::Vector &, mfem::Vector & )>;
	class DGSpace {
		public:
			DGSpace( std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order ) :
				Order( order ),
				Dim( 2 ),
				mesh( meshPtr )
		{

			dg_coll = new mfem::DG_FECollection( Order, Dim );
			face_coll = new mfem::DG_Interface_FECollection( Order, Dim );

			V_space = new mfem::FiniteElementSpace( mesh.get(), dg_coll, Dim );
			W_space = new mfem::FiniteElementSpace( mesh.get(), dg_coll );
			M_space = new mfem::FiniteElementSpace( mesh.get(), face_coll );

			postproc_coll = new mfem::DG_FECollection( Order + 1, Dim );
			postproc_space = new mfem::FiniteElementSpace( mesh.get(), postproc_coll );

			dimV = V_space->GetVSize();
			dimW = W_space->GetVSize();
			dimM = M_space->GetVSize();

			bOffsets.SetSize( 4 );
			bOffsets[ 0 ] = 0;
			bOffsets[ 1 ] = dimV;
			bOffsets[ 2 ] = dimV + dimW;
			bOffsets[ 3 ] = dimV + dimW + dimM;
		};

			~DGSpace() {
				delete V_space;
				delete W_space;
				delete M_space;
				delete dg_coll;
				delete face_coll;
				delete postproc_space;
				delete postproc_coll;
			};

			mfem::FiniteElementSpace const * QSpace() const { return V_space; };
			mfem::FiniteElementSpace const * USpace() const { return W_space; };
			mfem::FiniteElementSpace const * MSpace() const { return M_space; };
			mfem::FiniteElementSpace const * UStarSpace() const { return postproc_space; };

			mfem::FiniteElementSpace * QSpace()  { return V_space; };
			mfem::FiniteElementSpace * USpace()  { return W_space; };
			mfem::FiniteElementSpace * MSpace()  { return M_space; };
			mfem::FiniteElementSpace * UStarSpace() { return postproc_space; };

			std::shared_ptr<mfem::Mesh> Mesh() { return mesh; };

			void Update() {
				V_space->Update(true);
				W_space->Update(true);
				M_space->Update(false);
				postproc_space->Update( true );

				int dimV = V_space->GetVSize();
				int dimW = W_space->GetVSize();
				int dimM = M_space->GetVSize();

				bOffsets.SetSize( 4 );
				bOffsets[ 0 ] = 0;
				bOffsets[ 1 ] = dimV;
				bOffsets[ 2 ] = dimV + dimW;
				bOffsets[ 3 ] = dimV + dimW + dimM;
			};

			mfem::Array<int> const & GetOffsets() { return bOffsets; };
			int Order,Dim;
		protected:
			// Owned
			mfem::FiniteElementCollection *dg_coll, *face_coll;
			mfem::FiniteElementSpace *V_space,*W_space,*M_space;
			mfem::FiniteElementCollection *postproc_coll;
			mfem::FiniteElementSpace *postproc_space;

			// Unowned
			std::shared_ptr<mfem::Mesh> mesh; 
			int dimV,dimM,dimW;
			mfem::Array<int> bOffsets;

	};

	class Solution {
		private:
			std::shared_ptr<DGSpace> SolutionSpace;
			bool hasUStar;
			static const int dim = 2;
		public:
			mfem::Vector qu;
			mfem::Vector u_star;
			mfem::GridFunction q_variable,u_variable,u_hat_variable;
			mfem::GridFunction u_star_variable;

			Solution( std::shared_ptr<DGSpace> SolSpace )
				: SolutionSpace( SolSpace )
			{
				mfem::Array<int> const & offsets = SolutionSpace->GetOffsets();
				qu.SetSize( offsets[ 3 ] );

				q_variable.MakeRef( SolutionSpace->QSpace(), qu, offsets[ 0 ] );
				u_variable.MakeRef( SolutionSpace->USpace(), qu, offsets[ 1 ] );
				u_hat_variable.MakeRef( SolutionSpace->MSpace(), qu, offsets[ 2 ] );
				hasUStar = false;
			}

			Solution( std::shared_ptr<DGSpace> SolSpace, double *data )
				: SolutionSpace( SolSpace )
			{
				mfem::Array<int> const & offsets = SolutionSpace->GetOffsets();
				qu.SetDataAndSize( data, offsets[ 3 ] );

				q_variable.MakeRef( SolutionSpace->QSpace(), qu, offsets[ 0 ] );
				u_variable.MakeRef( SolutionSpace->USpace(), qu, offsets[ 1 ] );
				u_hat_variable.MakeRef( SolutionSpace->MSpace(), qu, offsets[ 2 ] );
				hasUStar = false;
			}

			void Zero()
			{
				qu = 0.0;
				if ( hasUStar )
					u_star = 0.0;
			}


			void AllocateUStar() 
			{
				u_star.SetSize( SolutionSpace->UStarSpace()->GetVSize() );
				u_star_variable.MakeRef( SolutionSpace->UStarSpace(), u_star, 0 );
				hasUStar = true;
			}

			void Reset()
			{
				mfem::Array<int> const & offsets = SolutionSpace->GetOffsets();
				qu.SetSize( offsets[ 3 ] );

				q_variable.MakeRef( SolutionSpace->QSpace(), qu, offsets[ 0 ] );
				u_variable.MakeRef( SolutionSpace->USpace(), qu, offsets[ 1 ] );
				u_hat_variable.MakeRef( SolutionSpace->MSpace(), qu, offsets[ 2 ] );
				u_star.Destroy();
				hasUStar = false;
			}


			/* 
			 * Prolongs a vector containing q & u from the old mesh to the 
			 * new. This will not handle increasing polynomial order.
			 * Should be called after the underlying SolutionSpace is updated
			 */

			void Prolong()
			{
				const mfem::Operator* U_update = SolutionSpace->USpace()->GetUpdateOperator();
				const mfem::Operator* Q_update = SolutionSpace->QSpace()->GetUpdateOperator();
				int U_old_dim = U_update->Width();
				int U_new_dim = U_update->Height();
				int Q_old_dim = Q_update->Width();
				int Q_new_dim = Q_update->Height();
				mfem::Vector qu_old( qu );

				qu.SetSize( U_new_dim + Q_new_dim + SolutionSpace->MSpace()->GetVSize() );

				// So the new Lambda variable is zero.
				qu = 0.;

				mfem::Array<int> oldOffsets; oldOffsets.SetSize( 4 );
				mfem::Array<int> const &bOffsets = SolutionSpace->GetOffsets();
				oldOffsets[ 0 ] = 0; oldOffsets[ 1 ] = Q_old_dim; oldOffsets[ 2 ] = U_old_dim + Q_old_dim;
				oldOffsets[ 3 ] = qu_old.Size();
				mfem::BlockVector QU_old_blk( qu_old.GetData(), oldOffsets );
				mfem::BlockVector QU_new_blk( qu.GetData(), bOffsets );

				Q_update->Mult( QU_old_blk.GetBlock( 0 ), QU_new_blk.GetBlock( 0 ) );
				U_update->Mult( QU_old_blk.GetBlock( 1 ), QU_new_blk.GetBlock( 1 ) );
				// perhaps also prolong u* if we have it?
			}

			std::tuple<double,double,double> l2_errors( RealScalarField uFun_ex, RealVectorField qFun_ex ) {
				int order = SolutionSpace->Order;
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
				double err_u_star;

				if ( hasUStar )
					err_u_star = u_star_variable.ComputeL2Error( ucoeff, irs );
				else
					err_u_star = std::nan("");

				return { err_u, err_q, err_u_star };
			};

			std::tuple<double,double,double> l2_diff( Solution & other ) {
				int order = std::max( SolutionSpace->Order, other.SolutionSpace->Order );
				int order_quad = std::max(2, 2*order+4);
				const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
				for (int i=0; i < mfem::Geometry::NumGeom; ++i)
				{
					irs[i] = &(mfem::IntRules.Get(i, order_quad));
				}

				mfem::GridFunctionCoefficient ucoeff( &other.u_variable );
				mfem::GridFunctionCoefficient qcoeff( &other.q_variable );
				// compare with the u* if they have one
				if ( other.hasUStar )
					ucoeff.SetGridFunction( &other.u_variable );

				double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
				double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
				double err_u_star;

				if ( hasUStar )
				{
					err_u_star = u_star_variable.ComputeL2Error( ucoeff, irs );
				}
				else
					err_u_star = std::nan("");

				return { err_u, err_q, err_u_star };
			};

			void WriteOutputMFEM( std::string prefix )
			{
				std::string meshFile( prefix );
				std::string GradPsiFile( prefix );
				std::string PsiFile( prefix );

				meshFile += ".mesh";
				GradPsiFile += "_grad_psi.gf";
				PsiFile += "_psi.gf";

				std::ofstream mesh_ofs( meshFile );
				SolutionSpace->Mesh()->Print( mesh_ofs );

				std::ofstream q_solution_ofs( GradPsiFile );
				q_solution_ofs.precision(8);
				q_variable.Save(q_solution_ofs);

				std::ofstream u_solution_ofs( PsiFile );
				u_solution_ofs.precision(8);
				if ( !hasUStar )
					u_variable.Save(u_solution_ofs);
				else
					u_star_variable.Save( u_solution_ofs );
			}

			void WriteOutputMFEM( meq::Configuration const& config )
			{
				std::string meshFile( config.GetFinalMeshFile() );
				std::string GradPsiFile( config.GetGradPsiFile() );
				std::string PsiFile( config.GetPsiFile() );

				std::ofstream mesh_ofs( meshFile );
				SolutionSpace->Mesh()->Print( mesh_ofs );

				std::ofstream q_solution_ofs( GradPsiFile );
				q_solution_ofs.precision(8);
				q_variable.Save(q_solution_ofs);

				std::ofstream u_solution_ofs( PsiFile );
				u_solution_ofs.precision(8);
				if ( !hasUStar )
					u_variable.Save(u_solution_ofs);
				else
					u_star_variable.Save( u_solution_ofs );
			}
	};


}
#endif // SOLUTION_HPP

