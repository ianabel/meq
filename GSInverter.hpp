#ifndef GSINVERTER_HPP
#define GSINVERTER_HPP

#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"


class GSInverter : public mfem::Operator
{

	protected:
		mfem::Mesh *mesh; // Unowned
		unsigned int Order,Dim;
		int dimV,dimM,dimW;

		// Owned

		mfem::Coefficient *boundary_conditions;

		mfem::FiniteElementCollection *dg_coll, *face_coll;
		mfem::FiniteElementSpace *V_space,*W_space,*M_space;
		mfem::HDGBilinearForm *AVarf;


		const double tau_D;


		mfem::FiniteElementCollection *postproc_coll;
		mfem::FiniteElementSpace *postproc_space;

		mfem::Array<int> bOffsets;
	public:
		GSInverter(mfem::Mesh *meshPtr, unsigned int order);
		mfem::FiniteElementSpace const* GetQSpace() const { return V_space; };
		mfem::FiniteElementSpace const* GetUSpace() const { return W_space; };
		mfem::FiniteElementSpace const* GetUStarSpace() const { return postproc_space; };


		mfem::FiniteElementSpace * GetQSpace()  { return V_space; };
		mfem::FiniteElementSpace * GetUSpace()  { return W_space; };
		mfem::FiniteElementSpace * GetUStarSpace()  { return postproc_space; };

		mfem::FiniteElementSpace * GetMSpace() { return M_space;};
		mfem::FiniteElementSpace const* GetMSpace() const { return M_space;};

		mfem::Mesh *GetMesh() {return mesh;};
		mfem::Mesh const *GetMesh() const {return mesh;};

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
			delete V_space;
			delete W_space;
			delete M_space;
			delete AVarf;
			delete dg_coll;
			delete face_coll;
			delete postproc_space;
			delete postproc_coll;
		};
		
};

#endif // GSINVERTER_HPP
