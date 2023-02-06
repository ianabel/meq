#ifndef HDGGSINTEGRATOR_HPP
#define HDGGSINTEGRATOR_HPP
// Integrators for the HDG discretizations
// Contributed by: T. Horvath, S. Rhebergen, A. Sivas
//                 University of Waterloo


#include "mfem.hpp"

namespace mfem
{

//---------------------------------------------------------------------
//---------------------------------------------------------------------

// Diffusion integrator: to compute all the domain based integrals
//
// The output is
//
//         [elemmat1 elemmat2]
// elmat = [elemmat3   0.0   ]
//
// elemmat1 = -(\nu^{-1} q, v)
// elemmat2 = (u, div(v))
// elemmat3 = (div(q), w)
//
// elemmat3 = elemmat2^T
//
// \nu is the constant diffusion coefficient
class HDGDomainIntegratorGS: public BilinearFormIntegrator
{
private:
   Vector shape, divshape;
   DenseMatrix partelmat, dshape, gshape, Jadj;

public:
   HDGDomainIntegratorGS() {};

   using BilinearFormIntegrator::AssembleElementMatrix2FES;
   virtual void AssembleElementMatrix2FES(const FiniteElement &fe_q,
                                          const FiniteElement &fe_u,
                                          ElementTransformation &Trans,
                                          DenseMatrix &elmat);

};

// Diffusion integrator to compute all the face based integrals
//
// The output is
//
//          [ 0.0   0.0  ]
// elmat1 = [ 0.0 local2 ]  - the face based integral for matrix A
//
//          [ local1 ]
// elmat2 = [ local3 ]  - the face based integral for matrix B
//
// elmat3 = [ local4  local5 ]  - the face based integral for matrix C
//
// elmat4 = local6  - the face based integral for matrix D
//
// where
// local1 = < \lambda,v\cdot n>
// local2 = < \tau u, w>
// local3 = -< tau \lambda, w>
// local4 = < \lambda, v\cdot n>
// local5 = -< \tau \lambda, w>
// local6 = < \tau \lambda, \mu>
//
// q_diff_coeff is the constant diffusion coefficient
// local4 = local1^T
// local5 = local3^T
class HDGFaceIntegratorGS : public BilinearFormIntegrator
{
private:
   double tauD;

   Vector shapeu, shapeq, normal, shape_face;
   DenseMatrix shape_dot_n;

public:
   HDGFaceIntegratorGS(double a)
      { tauD = a; }

   using BilinearFormIntegrator::AssembleFaceMatrixOneElement2and1FES;
   virtual void AssembleFaceMatrixOneElement2and1FES(const FiniteElement &fe_q,
                                                     const FiniteElement &fe_u,
                                                     const FiniteElement &face_fe,
                                                     FaceElementTransformations &Trans,
                                                     const int elem1or2,
                                                     const bool onlyB,
                                                     DenseMatrix &elmat1,
                                                     DenseMatrix &elmat2,
                                                     DenseMatrix &elmat3,
                                                     DenseMatrix &elmat4) override;

};


}

// Boundary Integrator for free-boundary work
namespace mfem {
	class HDGBoundaryTraceIntegrator : public LinearFormIntegrator
	{
		protected:
			Coefficient &f;

		public:
			HDGBoundaryTraceIntegrator(Coefficient &f_ref )
				: f( f_ref ) {};

			using LinearFormIntegrator::AssembleRHSElementVect;
			virtual void AssembleRHSElementVect(const FiniteElement &el,
					ElementTransformation &Tr,
					Vector &elvect)
			{
				mfem_error("Not implemented -- This is a Boundary Integrator not a Domain Integrator \n");
			}

			virtual void AssembleRHSElementVect(const FiniteElement &el,
					FaceElementTransformations &Tr,
					Vector &elvect);
	};

	class HDGBoundaryNormalTraceIntegrator : public LinearFormIntegrator
	{
		protected:
			Coefficient &f;
			unsigned int f_order;

		public:
			HDGBoundaryNormalTraceIntegrator(Coefficient &f_ref, unsigned int f_order_in = 2 )
				: f( f_ref ), f_order( f_order_in ) {};

			using LinearFormIntegrator::AssembleRHSElementVect;
			virtual void AssembleRHSElementVect(const FiniteElement &el,
					ElementTransformation &Tr,
					Vector &elvect)
			{
				mfem_error("Not implemented -- This is a Boundary Integrator not a Domain Integrator \n");
			}

			virtual void AssembleRHSElementVect(const FiniteElement &el,
					FaceElementTransformations &Tr,
					Vector &elvect);
	};

};

#endif // HDGGSINTEGRATOR_HPP

