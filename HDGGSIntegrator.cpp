// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

// Implementation of Bilinear Form Integrators

#include "HDGGSIntegrator.hpp"

#include <cmath>
#include <algorithm>

namespace mfem
{

//---------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////
void HDGDomainIntegratorGS::AssembleElementMatrix2FES(
   const FiniteElement &fe_q,
   const FiniteElement &fe_u,
   ElementTransformation &Trans,
   DenseMatrix &elmat)
{
   // get the number of degrees of freedoms and the dimension of the problem
   int ndof_u = fe_u.GetDof();
   int ndof_q = fe_q.GetDof();
   int dim  = fe_q.GetDim();
   double norm;

   int vdim = dim ;

   // set the vector and matrix sizes
   dshape.SetSize (ndof_u, dim); // for nabla u_reference
   gshape.SetSize (ndof_u, dim); // for nabla u
   Jadj.SetSize (dim);  // the Jacobian
   divshape.SetSize (vdim*ndof_u);  // divergence of q
   shape.SetSize (ndof_q);  // shape of q (and u)

   // for vector diffusion the matrix is built up from partial matrices
   partelmat.SetSize(ndof_q);

   DenseMatrix elmat1, elmat2, elmat3;

   // setting the sizes of the local element matrices
   elmat1.SetSize(dim*ndof_q, dim*ndof_q);
   elmat2.SetSize(ndof_q, vdim*ndof_u);
   elmat3.SetSize(ndof_q, vdim*ndof_u);

   elmat.SetSize(dim*ndof_q + ndof_u);

   elmat1 = 0.0;
   elmat2 = 0.0;
   elmat3 = 0.0;
   elmat  = 0.0;

   // setting the order of integration
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order1 = 2 * fe_q.GetOrder();
      int order2 = 2 * fe_q.GetOrder() + Trans.OrderW();
      int order = std::max(order1, order2) + 1;
      ir = &IntRules.Get(fe_u.GetGeomType(), order);
   }

   for (int i = 0; i < ir -> GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      // compute the shape and the gradient values on the reference element
      fe_u.CalcDShape (ip, dshape);
      fe_q.CalcShape (ip, shape);

      // calculate the adjugate of the Jacobian
      Trans.SetIntPoint (&ip);
      CalcAdjugate(Trans.Jacobian(), Jadj);

      // Calculate the gradient of the function of the physical element
      Mult (dshape, Jadj, gshape);

      // the weight is the product of the integral weight and the
      // determinant of the Jacobian
      norm = ip.weight * Trans.Weight();
      MultVVt(shape, partelmat);

      double c = ip.weight;

      // transform the the matrix to divergence vector
      gshape.GradToDiv (divshape);

      // multiply by R
		Vector point( 3 );
		Trans.Transform( ip, point );
      partelmat *= point( 0 );

      shape *= c;
      // compute the (u, \div v) term
      AddMultVWt (shape, divshape, elmat3);

      // assemble -(q, v) from the partial matrices
      partelmat *= norm*(-1.0);
      for (int k = 0; k < vdim; k++)
      {
         elmat1.AddMatrix(partelmat, ndof_q*k, ndof_q*k);
      }

   }

   elmat2.Transpose(elmat3);

   int block_size1 = dim*ndof_q;
   int block_size2 = ndof_u;

   for (int i = 0; i<block_size1; i++)
   {
      for (int j = 0; j<block_size1; j++)
      {
         elmat(i,j) = elmat1(i,j);
      }
      for (int j = 0; j<block_size2; j++)
      {
         elmat(i,j+block_size1) = elmat2(i,j);
      }
   }

   for (int i = 0; i<block_size2; i++)
   {
      for (int j = 0; j<block_size1; j++)
      {
         elmat(i+block_size1,j) = elmat3(i,j);
      }
   }
}


void HDGFaceIntegratorGS::AssembleFaceMatrixOneElement2and1FES(
   const FiniteElement &fe_qL,
   const FiniteElement &fe_qR,
   const FiniteElement &fe_uL,
   const FiniteElement &fe_uR,
   const FiniteElement &face_fe,
   FaceElementTransformations &Trans,
   const int elem1or2,
   const bool onlyB,
   DenseMatrix &elmat1,
   DenseMatrix &elmat2,
   DenseMatrix &elmat3,
   DenseMatrix &elmat4)
{
   // Get DoF from faces and the dimension
   int ndof_face = face_fe.GetDof();
   int ndof_q, ndof_u;
   int dim = fe_qL.GetDim();
   int vdim = dim;
   int order, order1, order2, order3, order4;

   DenseMatrix shape1_n_mtx;

   // set the dofs for u and q according to the element
   if (elem1or2 == 1)
   {
      ndof_u = fe_uL.GetDof();
      ndof_q = fe_qL.GetDof();
   }
   else
   {
      ndof_u = fe_uR.GetDof();
      ndof_q = fe_qR.GetDof();
   }

   DenseMatrix local1, local2, local3, local4, local5, local6;

   // set the shape functions, the normal and the advection
   shapeu.SetSize(ndof_u);
   shapeq.SetSize(ndof_q);
   shape_face.SetSize(ndof_face);
   normal.SetSize(dim);

   // set the proper size for the matrices
   local1.SetSize(vdim*ndof_q, ndof_face);
   local1 = 0.0;
   local2.SetSize(ndof_u, ndof_u);
   local2 = 0.0;
   local3.SetSize(ndof_u, ndof_face);
   local3 = 0.0;
   local4.SetSize(vdim*ndof_q, ndof_face);
   local4 = 0.0;
   local5.SetSize(ndof_u, ndof_face);
   local5 = 0.0;
   local6.SetSize(ndof_face, ndof_face);
   local6 = 0.0;

   int sub_block_size1 = vdim*ndof_q;
   int sub_block_size2 = ndof_u;

   int block_size1 = sub_block_size1 + sub_block_size2;
   int block_size2 = ndof_face;

   elmat1.SetSize(block_size1);
   elmat1 = 0.0;
   elmat2.SetSize(block_size1, block_size2);
   elmat2 = 0.0;
   elmat3.SetSize(block_size2, block_size1);
   elmat3 = 0.0;
   elmat4.SetSize(block_size2);
   elmat4 = 0.0;


   shape1_n_mtx.SetSize(ndof_q,dim);
   shape_dot_n.SetSize(ndof_q,dim);

   // set the order of integration
   // using the fact that q and u has the same order!
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      if (Trans.Elem2No >= 0)
      {
         order1 = std::max(fe_qL.GetOrder(), fe_qR.GetOrder());
         order2 = 2*fe_qR.GetOrder();
      }
      else
      {
         order1 = fe_qL.GetOrder();
         order2 = 2*fe_qL.GetOrder();
      }
      order1 += face_fe.GetOrder();
      order1 += 2;
      order2 += 2;

      order3 = std::max(order1, order2);
      order4 = 2*face_fe.GetOrder();
      order4 += 2;

      order = std::max(order3, order4);

      // IntegrationRule depends on the Geometry of the face (pont, line, triangle, rectangular)
      ir = &IntRules.Get(Trans.FaceGeom, order);
   }

   for (int p = 0; p < ir->GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir->IntPoint(p);
      IntegrationPoint eip1, eip2; // integration points on the elements

      // Trace finite element shape function
      Trans.Face->SetIntPoint(&ip);
      face_fe.CalcShape(ip, shape_face);

      // calculate the normal at the integration point
      if (dim == 1)
      {
         normal(0) = 2*eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(Trans.Face->Jacobian(), normal);
      }

      if (elem1or2 == 1)
      {
         // Side 1 finite element shape function
         Trans.Loc1.Transform(ip, eip1);
         fe_uL.CalcShape(eip1, shapeu);
         fe_qL.CalcShape(eip1, shapeq);
         MultVWt(shapeq, normal, shape_dot_n) ;
      }
      else
      {

         // Side 2 finite element shape function
         Trans.Loc2.Transform(ip, eip2);
         fe_uR.CalcShape(eip2, shapeu);
         fe_qR.CalcShape(eip2, shapeq);
         MultVWt(shapeq, normal, shape_dot_n) ;
      }

      // set the coefficients for the different terms
      // if the normal is involved Trans.Face->Weight() is not required
      double w1 = ip.weight*(-1.0);

      if (elem1or2 == 2)
      {
         w1 *=-1.0;
      }

      double w2 = tauD*Trans.Face->Weight()* ip.weight;

      double w3 = -w2;

      // local1 = < \lambda,\nu v\cdot n>
      for (int i = 0; i < vdim; i++)
         for (int k = 0; k < ndof_q; k++)
            for (int j = 0; j < ndof_face; j++)
            {
               local1(i*ndof_q + k, j) +=
                  shape_face(j) * shape_dot_n(k,i) * w1;
            }

      // local2 =  < \tau u, w>
      // local3 = -< tau \lambda, w>
      // local5 = -< tau \lambda, w>
      for (int i = 0; i < ndof_u; i++)
      {
         for (int j = 0; j < ndof_u; j++)
         {
            local2(i, j) += w2 * shapeu(i) * shapeu(j);
         }

         for (int j = 0; j < ndof_face; j++)
         {
            local3(i, j) += w3 * shapeu(i) * shape_face(j);
         }
      }

      if (!onlyB)
      {
         // local6 = < \tau \lambda, \mu>

         AddMult_a_VVt(w2, shape_face, local6);
      }
   }

   local4.Transpose(local1);
   local5.Transpose(local3);

   for (int i = 0; i<sub_block_size2; i++)
   {
      for (int j = 0; j<sub_block_size2; j++)
      {
         elmat1(i+sub_block_size1,j+sub_block_size1) = local2(i,j);
      }
   }

   for (int i = 0; i<block_size2; i++)
   {
      for (int j = 0; j<sub_block_size1; j++)
      {
         elmat2(j,i) = local1(j,i);

         elmat3(i,j) = local4(i,j);
      }
      for (int j = 0; j<sub_block_size2; j++)
      {
         elmat2(j+sub_block_size1,i) = local3(j,i);

         elmat3(i,j+sub_block_size1) = local5(i,j);
      }
   }

   elmat4 = local6;
}

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
