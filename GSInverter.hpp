#ifndef GSINVERTER_HPP
#define GSINVERTER_HPP

#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"

#include <memory>


class GSInverter : public mfem::Operator
{

	protected:
		std::shared_ptr<mfem::Mesh> mesh; // Unowned
		unsigned int Order,Dim;
		int dimV,dimM,dimW;

		// Owned

		mfem::Coefficient *boundary_conditions;

		using SharedFES = std::shared_ptr<mfem::FiniteElementSpace>;
		using SharedFEC = std::shared_ptr<mfem::FiniteElementCollection>;

		SharedFEC dg_coll, face_coll;
		SharedFES V_space, W_space, M_space;
		mfem::HDGBilinearForm *AVarf;

		const double tau_D;

		SharedFEC postproc_coll;
		SharedFES postproc_space;

		mfem::Array<int> bOffsets;
	public:
		GSInverter(std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order);

		SharedFES QSpace() { return V_space; };
		SharedFES USpace() { return W_space; };
		SharedFES MSpace() { return M_space; };

		const SharedFES GetQSpace() { return V_space; } const;
		const SharedFES GetUSpace() { return W_space; } const;
		const SharedFES GetMSpace() { return M_space; } const;

		SharedFES GetUStarSpace()  { return postproc_space; };
		const SharedFES GetUStarSpace() { return postproc_space; } const;

		std::shared_ptr<mfem::Mesh> GetMesh() {return mesh;};
		const std::shared_ptr<mfem::Mesh> GetMesh() { return mesh; } const;

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
