#include "HDGGSIntegrator.hpp"
#include "StdFnCoeffs.hpp"

class GSInverter : public mfem::Operator
{

	protected:
		mfem::Mesh *mesh; // Unowned
		unsigned int Order,Dim;
		int dimV,dimM,dimW;

		// Owned

		mfem::FiniteElementCollection *dg_coll, *face_coll;
		mfem::FiniteElementSpace *V_space,*W_space,*M_space;
		mfem::HDGBilinearForm3 *AVarf;


		const double tau_D;


		mfem::GridFunction BoundaryConditions;

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


		void QUUpdate( mfem::Vector const& qu_old, mfem::Vector &qu_new ) const;

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

		void Solve( mfem::Vector& qu_out )
		{
			mfem::Vector rhs_F( solver.NumCols() );
			qu_out.SetSize( solver.NumRows() );
			// Assemble the RHS and the Schur complement
			mfem::LinearForm *fform = new mfem::LinearForm;
			mfem::StdFunctionCoefficient fcoeff( RHS );
			fform->AddDomainIntegrator( new mfem::DomainLFIntegrator( fcoeff ) );
			fform->Update(solver.GetUSpace(), rhs_F, 0);
			fform->Assemble();

			solver.Mult( rhs_F, qu_out );
			delete fform;
		};

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


