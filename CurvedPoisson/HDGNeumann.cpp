
#include "HDGNeumann.hpp"

using namespace mfem;

void HDGNeumannLFIntegrator::AssembleRHSElementVect(
   const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
{
   mfem_error("Not implemented \n");
}

void HDGNeumannLFIntegrator::AssembleRHSElementVect(
   const FiniteElement &face, FaceElementTransformations &Tr,
   Vector &favect)
{
   int dim, ndof_face;
   double w, uin;

   ndof_face = face.GetDof();

   Vector shape_f(ndof_face);
   favect.SetSize(ndof_face);
   favect = 0.0;

   if (Tr.Elem2No >= 0)
   {
      // Interior face, do nothing
   }
   else
   {
      // Boundary face
      const IntegrationRule *ir = IntRule;
      if (ir == NULL)
      {
         int order = 2 * face.GetOrder();
         if (face.GetMapType() == FiniteElement::VALUE)
         {
            order += Tr.Face->OrderW();
         }

         ir = &IntRules.Get(Tr.FaceGeom, order);
      }

      for (int p = 0; p < ir->GetNPoints(); p++)
      {
         const IntegrationPoint &ip = ir->IntPoint(p);

			Tr.SetAllIntPoints( &ip );

			const IntegrationPoint &eip = Tr.GetElement1IntPoint();

			double val = Tr.Face->Weight() * ip.weight * gN.Eval( *Tr.Elem1, eip );

			face.CalcShape( eip, shape_f );

			add( favect, val, shape_f, favect );
      }
   }
}
