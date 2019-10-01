
#include "meq.hpp"
#include "GSInverter.hpp"
#include "NLRHSIntegrator.hpp"
#include "CockburnEstimator.hpp"

namespace meq {

	class GSSolver : public mfem::Operator
	{
		public:
			using Func = std::function< double( const mfem::Vector&, double ) >;
			using SharedFES = std::shared_ptr<mfem::FiniteElementSpace>;
			using SharedFEC = std::shared_ptr<mfem::FiniteElementCollection>;
		protected:
			GSInverter solver;
			Func PlasmaRHS;
		public:
			GSSolver(std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order, Func J_RHS )
				: solver( meshPtr, order ), PlasmaRHS( J_RHS )
			{
				height = solver.NumRows();
				width = solver.NumRows();
			};

			void SetBCs( mfem::Coefficient& coeff ) {
				solver.SetBCs( coeff );
			};

			void Update() { 
				solver.Update(); 
				height = solver.NumRows();
				width = solver.NumRows();
			};

			virtual void Mult( mfem::Vector const& qu_in, mfem::Vector & qu_out ) const override
			{
				mfem::Vector rhs_F( solver.NumCols() );
				qu_out.SetSize( solver.NumRows() );
				// Assemble the RHS and the Schur complement
				mfem::LinearForm *fform = new mfem::LinearForm;

				mfem::GridFunction u;
				u.MakeRef( const_cast<mfem::FiniteElementSpace* >( solver.USpace().get() ), static_cast<double*>( qu_in.GetData() + solver.QSpace()->GetVSize() ) );



				fform->AddDomainIntegrator( new NonlinearDomainLFIntegrator( u, PlasmaRHS, 4, 2 ) );
				fform->Update(const_cast<mfem::FiniteElementSpace* >( solver.USpace().get() ), rhs_F, 0);
				fform->Assemble();

				solver.Mult( rhs_F, qu_out );
				delete fform;
			};

			void Postprocess( mfem::GridFunction &u_out, mfem::Vector & qu_in )
			{
				solver.Postprocess( u_out, qu_in );
			};

			SharedFES MSpace() { return solver.MSpace(); };
			SharedFES QSpace() { return solver.QSpace(); };
			SharedFES USpace() { return solver.USpace(); };
			SharedFES UStarSpace() { return solver.UStarSpace(); };

			void Prolong( mfem::Vector const& qu_old, mfem::Vector & qu_new ) const
			{
				solver.Prolong( qu_old, qu_new );
			};

			void ApplyAdaptiveRefinement(  mfem::Vector & soln_vector )
			{
				mfem::GridFunction q_variable,u_variable,u_hat_variable;

				q_variable.MakeRef( solver.QSpace().get(), soln_vector, 0 );
				u_variable.MakeRef( solver.USpace().get(), soln_vector, solver.QSpace()->GetVSize() );
				u_hat_variable.MakeRef( solver.MSpace().get(), soln_vector, solver.QSpace()->GetVSize() + solver.USpace()->GetVSize() );


				mfem::GridFunction u_star( solver.UStarSpace().get() );
				solver.Postprocess( u_star, soln_vector );

				mfem::GradShafranovEstimator errorEstimator( q_variable, u_star, u_hat_variable, PlasmaRHS );
				mfem::DoerflerMarkingRefiner refiner( errorEstimator );

				refiner.SetGamma( 0.5 );
				refiner.Apply( *( solver.Mesh().get() ) );
				Update();
			};
	};

}
