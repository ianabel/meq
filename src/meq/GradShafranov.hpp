#ifndef GSSOLVER_HPP
#define GSSOLVER_HPP

#include "mfem.hpp"

#include "Solution.hpp"
#include "HDGGSIntegrator.hpp"
#include "CockburnEstimator.hpp"


#include <memory>

namespace meq {

	class GSSolver : public mfem::Operator
	{
		public:
			using Func = std::function< double( const mfem::Vector&, double ) >;
		protected:
			int Order, Dim;
		public:
			std::shared_ptr<DGSpace> SolutionSpace;
		protected:
			// Owned
			mfem::Coefficient *boundary_conditions;
			mfem::HDGBilinearForm *AVarf;

			const double tau_D;

			Func PlasmaRHS;
		public:
			GSSolver(std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order, Func JPlasma);

			mfem::FiniteElementSpace *QSpace() const { return SolutionSpace->QSpace(); };
			mfem::FiniteElementSpace *USpace() const { return SolutionSpace->USpace(); };
			mfem::FiniteElementSpace *MSpace() const { return SolutionSpace->MSpace(); };

			mfem::FiniteElementSpace *UStarSpace() const { return SolutionSpace->UStarSpace(); };

			void SetBCs( mfem::Coefficient& coeff );

			virtual void Mult( const mfem::Vector& u_in, mfem::Vector& qu_out ) const override;
			void Solve( Solution &soln );

			void Postprocess( Solution &soln ) const;

			void Update();
			void ApplyAdaptiveRefinement( Solution& soln );

			~GSSolver()
			{
				delete AVarf;
			};
	};

	class NonlinearDomainLFIntegrator : public mfem::LinearFormIntegrator
	{
		using NLFunc = std::function< double( const mfem::Vector &, double ) >;

		mfem::Coefficient *Q;
		NLFunc F;
		mfem::GridFunction &u;

		int oa, ob;
		public:

		/// Constructs a domain integrator with a given Coefficient
		NonlinearDomainLFIntegrator(mfem::Coefficient &QF, mfem::GridFunction &u_fn, NLFunc f, int a = 2, int b = 0)
			: Q(&QF), F( f ), u( u_fn ), oa(a), ob(b) 
		{
		};

		NonlinearDomainLFIntegrator(mfem::GridFunction &u_fn, NLFunc f, int a = 2, int b = 0)
			: Q(nullptr), F( f ), u( u_fn ), oa(a), ob(b) 
		{
		};

		/** Given a particular Finite Element and a transformation (Tr)
		  computes the element right hand side element vector, elvect. */
		virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
				mfem::ElementTransformation &Tr,
				mfem::Vector &elvect);

		using mfem::LinearFormIntegrator::AssembleRHSElementVect;
	};
};

#endif // GSSOLVER_HPP
