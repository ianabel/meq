#ifndef GSINVERTER_HPP
#define GSINVERTER_HPP

#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"

#include <memory>

class DGSpace {
	public:
		DGSpace( std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order ) :
			mesh( meshPtr ),
			Order( order ),
			Dim( 2 )
	{

		dg_coll = new DG_FECollection( Order, Dim );
		face_coll = new DG_Interface_FECollection( Order, Dim );

		V_space = new FiniteElementSpace( mesh.get(), dg_coll, Dim );
		W_space = new FiniteElementSpace( mesh.get(), dg_coll );
		M_space = new FiniteElementSpace( mesh.get(), face_coll );

		postproc_coll = new DG_FECollection( Order + 1, Dim );
		postproc_space = new FiniteElementSpace( mesh.get(), postproc_coll );

		dimV = V_space->GetVSize();
		dimW = W_space->GetVSize();
		dimM = M_space->GetVSize();

		bOffsets.SetSize( 4 );
		bOffsets[ 0 ] = 0;
		bOffsets[ 1 ] = dimV;
		bOffsets[ 2 ] = dimV + dimW;
		bOffsets[ 3 ] = dimV + dimW + dimM;


	}

		~DGSpace() {
			delete V_space;
			delete W_space;
			delete M_space;
			delete dg_coll;
			delete face_coll;
		};

		mfem::FiniteElementSpace const * QSpace() const { return V_space; };
		mfem::FiniteElementSpace const * USpace() const { return W_space; };
		mfem::FiniteElementSpace const * MSpace() const { return M_space; };
		mfem::FiniteElementSpace const * UStarSpace() const { return postproc_space; };

		mfem::FiniteElementSpace * QSpace()  { return V_space; };
		mfem::FiniteElementSpace * USpace()  { return W_space; };
		mfem::FiniteElementSpace * MSpace()  { return M_space; };
		mfem::FiniteElementSpace * UStarSpace() { return postproc_space; };


		mfem::Array<int> const & GetOffsets() { return bOffsets; };
	protected:
		int Order,Dim;
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

class GSInverter : public mfem::Operator
{
	public:
		using SharedFES = std::shared_ptr<mfem::FiniteElementSpace>;
		using SharedFEC = std::shared_ptr<mfem::FiniteElementCollection>;
	protected:
		DGSpace SolutionSpace;

		// Owned

		mfem::Coefficient *boundary_conditions;


		SharedFEC dg_coll, face_coll;
		SharedFES V_space, W_space, M_space;
		mfem::HDGBilinearForm *AVarf;

		const double tau_D;

		SharedFEC postproc_coll;
		SharedFES postproc_space;

	public:
		GSInverter(std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order);

		SharedFES QSpace() const { return V_space; };
		SharedFES USpace() const { return W_space; };
		SharedFES MSpace() const { return M_space; };

		SharedFES UStarSpace() const { return postproc_space; };

		std::shared_ptr<mfem::Mesh> Mesh() const {return mesh;};

		void Prolong( mfem::Vector const& soln_old, mfem::Vector &soln_new ) const;

		mfem::Array<int> const & GetOffsets() { return bOffsets; };
		void SetBCs( mfem::Coefficient& coeff );

		// Computes (-Delta^*)^{-1}( qu_in )
		// this class chould not be used directly, but is intended to be wrapped in 
		// a class that constructs the RHS vector.
		virtual void Mult( const mfem::Vector& u_in , mfem::Vector& qu_out ) const ;

		void Postprocess( mfem::GridFunction &u_out, mfem::Vector &qu_in ) const;

		void Update();

		~GSInverter()
		{
			delete AVarf;
		};
		
};

#endif // GSINVERTER_HPP
