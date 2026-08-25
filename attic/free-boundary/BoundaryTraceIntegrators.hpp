#ifndef ATTIC_BOUNDARYTRACEINTEGRATORS_HPP
#define ATTIC_BOUNDARYTRACEINTEGRATORS_HPP

// Extracted from HDGGSIntegrator.hpp at tag v0-legacy. Free-boundary only:
// these assemble the Green's-function boundary data of the von Hagenow /
// Lackner scheme. Not ported to MFEM 4.9.1 -- see the README in this directory.

#include "mfem.hpp"

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

#endif // ATTIC_BOUNDARYTRACEINTEGRATORS_HPP
