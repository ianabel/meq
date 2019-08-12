
#include "mfem.hpp"
#include <functional>

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

void NonlinearDomainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                mfem::ElementTransformation &Tr,
                                                mfem::Vector &elvect)
{
   int dof = el.GetDof();
	mfem::Vector shape( dof );

   elvect.SetSize(dof);
   elvect = 0.0;

   const mfem::IntegrationRule *ir;
   ir = &mfem::IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);

	mfem::Vector u_vals( ir->GetNPoints() );

	u.GetValues( Tr.ElementNo, *ir, u_vals );

	mfem::Vector point( 2 );


   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);

      Tr.SetIntPoint (&ip);
		Tr.Transform( ip, point );

      double val = Tr.Weight() * F( point, u_vals( i ) );

		if ( Q != nullptr )
			val *= Q->Eval( Tr, ip );

      el.CalcShape(ip, shape);
		shape *= ( ip.weight * val );

		elvect += shape;
   }
}

