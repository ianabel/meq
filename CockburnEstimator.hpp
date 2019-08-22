#ifndef COCKBURNESTIMATOR_HPP
#define COCKBURNESTIMATOR_HPP

#include "mfem.hpp"
#include <functional>


namespace mfem {



class CockburnZhangEstimator : public mfem::ErrorEstimator
{
	protected:
		mfem::Mesh *mesh;
		mfem::FiniteElementSpace *q_space;
		mfem::FiniteElementSpace *u_space;
		mfem::FiniteElementSpace *m_space;

		mfem::Vector localErrors;
		GridFunction &q_sol,&u_sol,&lambda;
		Coefficient &kappa;
		Coefficient &rhs;

		bool valid;
	protected:
		static double ElementDiameter( Mesh const *m, Element const* e )
		{
			unsigned int nVerts = e->GetNVertices();
			mfem::Array<int> VertIndices( nVerts );
			e->GetVertices( VertIndices );
			double const * pts[ 3 ];
			unsigned int N = m->SpaceDimension();
			double len[ 3 ];

			switch ( e->GetType() )
			{
				case mfem::Element::Type::SEGMENT:
					if ( nVerts != 2 )
						throw new std::logic_error( "Whargl, SEGMENT should have 2 vertices" );
					pts[ 0 ] = m->GetVertex( VertIndices[ 0 ] );
					pts[ 1 ] = m->GetVertex( VertIndices[ 1 ] );

					return mfem::Distance( pts[ 0 ], pts[ 1 ], N );
					break;

				case mfem::Element::Type::TRIANGLE:

					if ( nVerts != 3 )
						throw new std::logic_error( "Whargl, Triangle should have 3 vertices" );

					pts[ 0 ] = m->GetVertex( VertIndices[ 0 ] );
					pts[ 1 ] = m->GetVertex( VertIndices[ 1 ] );
					pts[ 2 ] = m->GetVertex( VertIndices[ 2 ] );

					len[ 0 ] = mfem::Distance( pts[ 0 ], pts[ 1 ], N );
					len[ 1 ] = mfem::Distance( pts[ 1 ], pts[ 2 ], N );
					len[ 2 ] = mfem::Distance( pts[ 2 ], pts[ 0 ], N );

					return std::max( std::max( len[ 0 ], len[ 1 ] ), len[ 2 ] );
					break;

				default:
					throw new std::logic_error( "Unsupported geometric element type." );
			}
			return -1.0;
		}
	public:
		CockburnZhangEstimator( GridFunction &q, GridFunction &u, GridFunction &lambda_ref, Coefficient &kappa_ref, Coefficient &rhs_ref )
			: q_sol( q ), u_sol( u ), lambda( lambda_ref ), kappa( kappa_ref ), rhs( rhs_ref ),valid( false )
		{
			q_space = q.FESpace();
			u_space = u.FESpace();
			m_space = lambda.FESpace();
			mesh = u_space->GetMesh();
		};


		double ComputeElementError( int i ) 
		{
			mfem::ElementTransformation *trans = mesh->GetElementTransformation( i );
			mfem::Element const *K = mesh->GetElement( i );

			double eta_1 = 0;
			double eta_2 = 0;

			unsigned int order = 2 * q_space->GetOrder(i) + 3;
			mfem::IntegrationRule const& CellIntegrator = mfem::IntRules.Get( K->GetType(), order );


			double h_K = ElementDiameter( mesh, K );

			// Sum over integration points
			for (int j=0; j < CellIntegrator.GetNPoints(); j++)
			{
				const IntegrationPoint &ip = CellIntegrator.IntPoint( j );
				trans->SetIntPoint( &ip );

				double eta_1_tmp;
				// eta_1 = h_K^2 * || F_RHS + div q ||^2
				eta_1_tmp = ( rhs.Eval( *trans, ip ) + q_sol.GetDivergence( *trans ) );
				eta_1 += h_K * h_K * eta_1_tmp * eta_1_tmp * ip.weight;

				Vector eta_2_tmp;
				// eta_2 = || q - kappa*(grad u) ||^2
				Vector q_val;
				Vector GradU;
				q_sol.GetVectorValue( i, ip, q_val );
				u_sol.GetGradient( *trans, GradU );
				double kappa_val = kappa.Eval( *trans, ip );
				GradU *= kappa_val;
				eta_2_tmp = q_val - GradU;	
				eta_2 += ( eta_2_tmp * eta_2_tmp ) * ip.weight;
			}

			double eta_3 = 0;
			double eta_4 = 0;
			double eta_5 = 0;

			mfem::Array<int> faces,orientations;
			mesh->GetElementEdges( i, faces, orientations );

			// Sum over all faces
			for (int j=0; j < faces.Size(); j++)
			{
				mfem::Element const *e = mesh->GetFace( faces[ j ] );
				double h_e = ElementDiameter( mesh, e );

				mfem::FaceElementTransformations *feTrans = mesh->GetFaceElementTransformations( faces[ j ] );


				// Integrate over e
				mfem::IntegrationRule const& EdgeIntegrator = mfem::IntRules.Get( e->GetType(), order );
				for (int k=0; k < EdgeIntegrator.GetNPoints(); k++)
				{
					const IntegrationPoint& ip = EdgeIntegrator.IntPoint( k );
					IntegrationPoint eip1,eip2; 
					feTrans->Loc1.Transform( ip, eip1 );
					feTrans->Elem1->SetIntPoint( &eip1 );
					feTrans->Face->SetIntPoint( &ip );

					// eta_5 = (.5/h_e) * || lambda - u ||^2
					double lambda_val,u_val;
					lambda_val = lambda.GetInterfaceValue( faces[ j ], ip );
					if ( feTrans->Elem1No == i )
						u_val = u_sol.GetValue( feTrans->Elem1No, eip1 );
					else if ( feTrans->Elem2No == i )
						u_val = u_sol.GetValue( feTrans->Elem2No, eip2 );
					else 
						throw new std::logic_error( "Element is neither of the ones attached to the face. Wat." );

					eta_5 += ( .5/h_e ) * ip.weight * ( lambda_val - u_val )*( lambda_val - u_val );

					if ( feTrans->Elem2No == -1 )
						continue;

					// eta_3 & eta_4 are only calculated for interior edges

					feTrans->Loc2.Transform( ip, eip2 );
					feTrans->Elem2->SetIntPoint( &eip2 );

					// eta_3 = (1/2) * h_e * || (q_plus - q_minus) . n ||^2
					Vector q_plus,q_minus;
					q_sol.GetVectorValue( feTrans->Elem1No, eip1, q_plus );
					q_sol.GetVectorValue( feTrans->Elem2No, eip2, q_minus );

					Vector normal( mesh->SpaceDimension() );	
					mfem::CalcOrtho( feTrans->Face->Jacobian(), normal );

					double q_jump = q_plus * normal - q_minus * normal;
					eta_3 += 0.5 * h_e * q_jump * q_jump;

					// eta_4 = .5 * h_e^{-1} || u_plus - u_minus ||^2
					double u_plus,u_minus;

					u_plus  = u_sol.GetValue( feTrans->Elem1No, eip1 );
					u_minus = u_sol.GetValue( feTrans->Elem2No, eip2 );

					eta_4 += ( .5 / h_e ) * ip.weight * ( u_plus - u_minus ) * ( u_plus - u_minus );

				}
			}
			return eta_1 + eta_2 + eta_3 + eta_4 + eta_5;
		}

		/// Get a Vector with all element errors.
		virtual const Vector &GetLocalErrors() 
		{
			if ( valid )
				return localErrors;

			unsigned int nElems = mesh->GetNE();
			localErrors.SetSize( nElems );
			for ( unsigned int j=0; j < nElems; j++ )
				localErrors( j ) = ComputeElementError( j );

			valid = true;
			return localErrors;
		}

		/// Force recomputation of the estimates on the next call to GetLocalErrors.
		virtual void Reset() 
		{
			valid = false;
		}
};

class GradShafranovEstimator : public mfem::ErrorEstimator
{
	protected:
		mfem::Mesh *mesh;
		mfem::FiniteElementSpace *q_space;
		mfem::FiniteElementSpace *u_space;
		mfem::FiniteElementSpace *m_space;

		mfem::Vector localErrors;
		GridFunction &q_sol,&u_sol,&lambda;



		bool valid;
	protected:
		static double ElementDiameter( Mesh const *m, Element const* e )
		{
			unsigned int nVerts = e->GetNVertices();
			mfem::Array<int> VertIndices( nVerts );
			e->GetVertices( VertIndices );
			double const * pts[ 3 ];
			unsigned int N = m->SpaceDimension();
			double len[ 3 ];

			switch ( e->GetType() )
			{
				case mfem::Element::Type::SEGMENT:
					if ( nVerts != 2 )
						throw new std::logic_error( "Whargl, SEGMENT should have 2 vertices" );
					pts[ 0 ] = m->GetVertex( VertIndices[ 0 ] );
					pts[ 1 ] = m->GetVertex( VertIndices[ 1 ] );

					return mfem::Distance( pts[ 0 ], pts[ 1 ], N );
					break;

				case mfem::Element::Type::TRIANGLE:

					if ( nVerts != 3 )
						throw new std::logic_error( "Whargl, Triangle should have 3 vertices" );

					pts[ 0 ] = m->GetVertex( VertIndices[ 0 ] );
					pts[ 1 ] = m->GetVertex( VertIndices[ 1 ] );
					pts[ 2 ] = m->GetVertex( VertIndices[ 2 ] );

					len[ 0 ] = mfem::Distance( pts[ 0 ], pts[ 1 ], N );
					len[ 1 ] = mfem::Distance( pts[ 1 ], pts[ 2 ], N );
					len[ 2 ] = mfem::Distance( pts[ 2 ], pts[ 0 ], N );

					return std::max( std::max( len[ 0 ], len[ 1 ] ), len[ 2 ] );
					break;

				default:
					throw new std::logic_error( "Unsupported geometric element type." );
			}
			return -1.0;
		}
	public:
		// Type for F((R,z),psi)
		using NLFunc = std::function<double( const mfem::Vector &, double )>;
		NLFunc  F;

		GradShafranovEstimator( GridFunction &q, GridFunction &u, GridFunction &lambda_ref,NLFunc F_ref )
			: q_sol( q ), u_sol( u ), lambda( lambda_ref ),valid( false ),F( F_ref )
		{
			q_space = q.FESpace();
			u_space = u.FESpace();
			m_space = lambda.FESpace();
			mesh = u_space->GetMesh();
		};


		double ComputeElementError( int i ) 
		{
			mfem::ElementTransformation *trans = mesh->GetElementTransformation( i );
			mfem::Element const *K = mesh->GetElement( i );

			double eta_1 = 0;
			double eta_2 = 0;

			unsigned int order = 2 * q_space->GetOrder(i) + 3;
			mfem::IntegrationRule const& CellIntegrator = mfem::IntRules.Get( K->GetType(), order );


			double h_K = ElementDiameter( mesh, K );

			// Sum over integration points
			for (int j=0; j < CellIntegrator.GetNPoints(); j++)
			{
				const IntegrationPoint &ip = CellIntegrator.IntPoint( j );
				trans->SetIntPoint( &ip );

				double eta_1_tmp;
				// eta_1 = h_K^2 * || F_RHS + div q ||^2
				Vector pt;
				trans->Transform( ip, pt );
				double R = pt( 0 );
				double psi_val = u_sol.GetValue( i, ip );
				eta_1_tmp = ( F( pt, psi_val )/R  + q_sol.GetDivergence( *trans ) );
				eta_1 += h_K * h_K * eta_1_tmp * eta_1_tmp * ip.weight;

				Vector eta_2_tmp;
				// eta_2 = || q - (grad u)/R ||^2
				Vector q_val;
				Vector GradU;
				q_sol.GetVectorValue( i, ip, q_val );
				u_sol.GetGradient( *trans, GradU );
				GradU /= R;
				eta_2_tmp = q_val - GradU;	
				eta_2 += ( eta_2_tmp * eta_2_tmp ) * ip.weight;
			}

			double eta_3 = 0;
			double eta_4 = 0;
			double eta_5 = 0;

			mfem::Array<int> faces,orientations;
			mesh->GetElementEdges( i, faces, orientations );

			// Sum over all faces
			for (int j=0; j < faces.Size(); j++)
			{
				mfem::Element const *e = mesh->GetFace( faces[ j ] );
				double h_e = ElementDiameter( mesh, e );

				mfem::FaceElementTransformations *feTrans = mesh->GetFaceElementTransformations( faces[ j ] );


				// Integrate over e
				mfem::IntegrationRule const& EdgeIntegrator = mfem::IntRules.Get( e->GetType(), order );
				for (int k=0; k < EdgeIntegrator.GetNPoints(); k++)
				{
					const IntegrationPoint& ip = EdgeIntegrator.IntPoint( k );
					IntegrationPoint eip1,eip2; 
					feTrans->Loc1.Transform( ip, eip1 );
					feTrans->Elem1->SetIntPoint( &eip1 );
					feTrans->Face->SetIntPoint( &ip );

					// eta_5 = (.5/h_e) * || lambda - u ||^2
					double lambda_val,u_val;
					lambda_val = lambda.GetInterfaceValue( faces[ j ], ip );
					if ( feTrans->Elem1No == i )
						u_val = u_sol.GetValue( feTrans->Elem1No, eip1 );
					else if ( feTrans->Elem2No == i )
						u_val = u_sol.GetValue( feTrans->Elem2No, eip2 );
					else 
						throw new std::logic_error( "Element is neither of the ones attached to the face. Wat." );

					eta_5 += ( .5/h_e ) * ip.weight * ( lambda_val - u_val )*( lambda_val - u_val );

					if ( feTrans->Elem2No == -1 )
						continue;

					// eta_3 & eta_4 are only calculated for interior edges

					feTrans->Loc2.Transform( ip, eip2 );
					feTrans->Elem2->SetIntPoint( &eip2 );

					// eta_3 = (1/2) * h_e * || (q_plus - q_minus) . n ||^2
					Vector q_plus,q_minus;
					q_sol.GetVectorValue( feTrans->Elem1No, eip1, q_plus );
					q_sol.GetVectorValue( feTrans->Elem2No, eip2, q_minus );

					Vector normal( mesh->SpaceDimension() );	
					mfem::CalcOrtho( feTrans->Face->Jacobian(), normal );

					double q_jump = q_plus * normal - q_minus * normal;
					eta_3 += 0.5 * h_e * q_jump * q_jump;

					// eta_4 = .5 * h_e^{-1} || u_plus - u_minus ||^2
					double u_plus,u_minus;

					u_plus  = u_sol.GetValue( feTrans->Elem1No, eip1 );
					u_minus = u_sol.GetValue( feTrans->Elem2No, eip2 );

					eta_4 += ( .5 / h_e ) * ip.weight * ( u_plus - u_minus ) * ( u_plus - u_minus );

				}
			}
			return eta_1 + eta_2 + eta_3 + eta_4 + eta_5;
		}

		/// Get a Vector with all element errors.
		virtual const Vector &GetLocalErrors() 
		{
			if ( valid )
				return localErrors;

			unsigned int nElems = mesh->GetNE();
			localErrors.SetSize( nElems );
			for ( unsigned int j=0; j < nElems; j++ )
				localErrors( j ) = ComputeElementError( j );

			valid = true;
			return localErrors;
		}

		/// Force recomputation of the estimates on the next call to GetLocalErrors.
		virtual void Reset() 
		{
			valid = false;
		}
};

}
#endif // COCKBURNESTIMATOR_HPP
