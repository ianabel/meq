#ifndef HDGNEUMANN_HPP
#define HDGNEUMANN_HPP

#include "mfem.hpp"

namespace mfem {

class HDGNeumannLFIntegrator : public LinearFormIntegrator
{
protected:
   Coefficient &gN;

public:
   HDGNeumannLFIntegrator( Coefficient &_u )
		: gN( _u )
	{
	}
   

   using LinearFormIntegrator::AssembleRHSElementVect;
   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect);

   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       FaceElementTransformations &Tr,
                                       Vector &elvect);
};

}
#endif // HDGNEUMANN_HPP
