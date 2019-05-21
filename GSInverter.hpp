#include "mfem.hpp"

class GSInverter : public mfem::Operator
{

	public:
		using RealFunc = double (*)( const mfem::Vector& );
	protected:
		mfem::Mesh *mesh; // Unowned
		unsigned int Order,Dim;
		int dimV,dimM,dimW;

		// Owned

		mfem::FiniteElementCollection *dg_coll, *face_coll;
		mfem::FiniteElementSpace *V_space,*W_space,*M_space;
		mfem::HDGBilinearForm3 *AVarf;


		mfem::FunctionCoefficient diffusion; // diffusion constant
		const double tau_D;

		RealFunc RHS_ref;

		mfem::GridFunction BoundaryConditions;

		mfem::Array<int> bOffsets;
		static double Rinv( const mfem::Vector& rz_pt ) { return 1.0/rz_pt( 0 );};
	public:
		GSInverter(mfem::Mesh *meshPtr, unsigned int order, RealFunc fRHS );
		mfem::FiniteElementSpace* GetQSpace() { return V_space; };
		mfem::FiniteElementSpace* GetUSpace() { return W_space; };

		void QUUpdate( mfem::Vector const& qu_old, mfem::Vector &qu_new ) const;

		mfem::Array<int> const & GetOffsets() { return bOffsets; };

		// Actually solve the problem:
		// which in this case doesn't depend on the input vector
		// and store in the Vector y
		virtual void Mult( const mfem::Vector& qu_in , mfem::Vector& qu_out ) const ;

		void Update();

		~GSInverter();
		{
			delete V_space;
			delete W_space;
			delete M_space;
			delete AVarf;
			delete dg_coll;
			delete face_coll;
		};
		
}
