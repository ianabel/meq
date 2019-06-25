#ifndef GSINVERTER_HPP
#define GSINVERTER_HPP

#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"
#include "StdFnCoeffs.hpp"


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
		mfem::HDGBilinearForm3 *AVarf;


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

class GSSolver
{
	public:
		using RealFunc = std::function< double( const mfem::Vector & )>;
	protected:
		GSInverter solver;
		RealFunc RHS;
	public:
		GSSolver(mfem::Mesh *meshPtr, unsigned int order, RealFunc fRHS )
			: solver( meshPtr, order ), RHS( fRHS )
		{
		};

		void SetBCs( mfem::Coefficient& coeff ) 
		{
			solver.SetBCs( coeff );
		};

		void Solve( mfem::Vector& soln_out )
		{
			mfem::Vector rhs_F( solver.NumCols() );
			// Assemble the RHS and the Schur complement
			mfem::LinearForm *fform = new mfem::LinearForm;
			mfem::StdFunctionCoefficient fcoeff( RHS );
			fform->AddDomainIntegrator( new mfem::DomainLFIntegrator( fcoeff ) );
			fform->Update(solver.GetUSpace(), rhs_F, 0);
			fform->Assemble();

			solver.Mult( rhs_F, soln_out );
			delete fform;
		};

		void Update() { 
			solver.Update(); 
		};

		void ApplyAdaptiveRefinement(  mfem::Vector & soln_vector )
		{
			mfem::StdFunctionCoefficient fFunCoeff( RHS );
			auto kappaF = []( const mfem::Vector& pt ) { return 1.0 / pt( 0 ); };
			mfem::StdFunctionCoefficient kappa( kappaF );

			mfem::GridFunction q_variable,u_variable,u_hat_variable;
			q_variable.MakeRef( solver.GetQSpace(), soln_vector, 0 );
			u_hat_variable.MakeRef( solver.GetMSpace(), soln_vector, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );

			mfem::GridFunction u_star( solver.GetUStarSpace() );
			solver.Postprocess( u_star, soln_vector );

			mfem::CockburnZhangEstimator errorEstimator( q_variable, u_star, u_hat_variable, kappa, fFunCoeff );
			mfem::ThresholdRefiner refiner( errorEstimator );
			refiner.SetTotalErrorFraction( 0.1 );
			refiner.Apply( *solver.GetMesh() );
			solver.Update();
		}

		void Postprocess( mfem::GridFunction &u_out, mfem::Vector & qu_in )
		{
			solver.Postprocess( u_out, qu_in );
		};

		mfem::FiniteElementSpace const* GetQSpace() const { return solver.GetQSpace(); };
		mfem::FiniteElementSpace const* GetUSpace() const { return solver.GetUSpace(); };
		mfem::FiniteElementSpace const* GetUStarSpace() const { return solver.GetUStarSpace(); };

		mfem::FiniteElementSpace * GetQSpace()  { return solver.GetQSpace(); };
		mfem::FiniteElementSpace * GetUSpace()  { return solver.GetUSpace(); };
		mfem::FiniteElementSpace * GetUStarSpace()  { return solver.GetUStarSpace(); };

		mfem::FiniteElementSpace * GetMSpace() { return solver.GetMSpace();};
		mfem::FiniteElementSpace const* GetMSpace() const { return solver.GetMSpace();};

		mfem::Array<int> const& GetOffsets() { return solver.GetOffsets();};
};


#endif // GSINVERTER_HPP
