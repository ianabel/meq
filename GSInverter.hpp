#ifndef GSINVERTER_HPP
#define GSINVERTER_HPP

#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"

#include <memory>


class GSInverter : public mfem::Operator
{
	public:
		using SharedFES = std::shared_ptr<mfem::FiniteElementSpace>;
		using SharedFEC = std::shared_ptr<mfem::FiniteElementCollection>;
	protected:
		std::shared_ptr<mfem::Mesh> mesh; // Unowned
		unsigned int Order,Dim;
		int dimV,dimM,dimW;

		// Owned

		mfem::Coefficient *boundary_conditions;


		SharedFEC dg_coll, face_coll;
		SharedFES V_space, W_space, M_space;
		mfem::HDGBilinearForm *AVarf;

		const double tau_D;

		SharedFEC postproc_coll;
		SharedFES postproc_space;

		mfem::Array<int> bOffsets;
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
