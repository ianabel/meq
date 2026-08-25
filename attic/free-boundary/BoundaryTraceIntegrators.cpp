// Extracted from HDGGSIntegrator.cpp at tag v0-legacy. See the header.

#include "BoundaryTraceIntegrators.hpp"

namespace mfem {

void HDGBoundaryTraceIntegrator::AssembleRHSElementVect(
		const FiniteElement &bdr_cell, FaceElementTransformations &Trans,
		Vector &favect)
{
		double w, f_val;

		if ( Trans.Elem2No >= 0 )
		{
			// Interior Face. Do Nothing.
			return;
		}

		int ndof_cell = bdr_cell.GetDof();

		Vector shape_f( ndof_cell );
		favect.SetSize( ndof_cell );
		favect = 0.0;


		// Boundary face
		const IntegrationRule *ir = IntRule;
		if (ir == NULL)
		{
			int order = 2 * bdr_cell.GetOrder();
			if (bdr_cell.GetMapType() == FiniteElement::VALUE)
			{
				order += Trans.Face->OrderW();
			}

			ir = &IntRules.Get(Trans.FaceGeom, order);
		}

		for (int p = 0; p < ir->GetNPoints(); p++)
		{
			const IntegrationPoint &ip = ir->IntPoint(p);

			// eip_L is inside the boundary cell, but is on the edge of it that corresponds
			// to the Boundary Face that we are integrating over
			IntegrationPoint eip_L;
			Trans.Loc1.Transform(ip, eip_L);
			bdr_cell.CalcShape(eip_L, shape_f);
			Trans.Face->SetIntPoint(&ip);

			f_val = f.Eval(*Trans.Elem1, eip_L);

			w = ip.weight * Trans.Face->Weight(); // The weight factor is exactly the normalization of normal

			for (int i = 0; i < ndof_cell; i++)
			{
				favect(i) += w * f_val * shape_f(i);
			}
		}
	};
	
	void HDGBoundaryNormalTraceIntegrator::AssembleRHSElementVect(
			const FiniteElement &bdr_cell, FaceElementTransformations &Trans,
			Vector &favect)
	{
		double w, f_val;

		if ( Trans.Elem2No >= 0 )
		{
			// Interior Face. Do Nothing.
			return;
		}

		int dim = Trans.Face->GetSpaceDim();
		int ndof_cell = bdr_cell.GetDof();

		Vector shape_f( ndof_cell );
		favect.SetSize( ndof_cell * dim );
		favect = 0.0;


		// Boundary face
		const IntegrationRule *ir = IntRule;
		if (ir == NULL)
		{
			int order = 2 * bdr_cell.GetOrder() + 2*f_order;
			if (bdr_cell.GetMapType() == FiniteElement::VALUE)
			{
				order += Trans.Face->OrderW();
			}

			ir = &IntRules.Get(Trans.FaceGeom, order);
		}

		for (int p = 0; p < ir->GetNPoints(); p++)
		{
			const IntegrationPoint &ip = ir->IntPoint(p);

			// eip_L is inside the boundary cell, but is on the edge of it that corresponds
			// to the Boundary Face that we are integrating over
			IntegrationPoint eip_L;
			Trans.Loc1.Transform(ip, eip_L);
			bdr_cell.CalcShape(eip_L, shape_f);
			Trans.Face->SetIntPoint(&ip);
			Vector normal( dim );
			CalcOrtho( Trans.Face->Jacobian(), normal );


			f_val = f.Eval(*Trans.Elem1, eip_L);

			w = ip.weight; // The weight factor is exactly the normalization of normal

			for ( int k = 0; k < dim; k++ )
			{
				double v_k = w * f_val * normal( k );
				for (int i = 0; i < ndof_cell; i++)
				{
					favect(ndof_cell * k + i) += v_k * shape_f(i);
				}
			}
		}
	};

}
