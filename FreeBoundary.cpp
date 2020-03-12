
#include "FreeBoundary.hpp"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace mfem;

namespace meq {

double GreensFunction( mfem::Vector const& r, mfem::Vector const& r_star )
{
	double R = r[ 0 ];
	double R_star = r_star[ 0 ];
	double Z = r[ 1 ];
	double Z_star = r_star[ 1 ];

	double pi = 3.14159265358979323844;
	double answer = 1.0 / ( 2.0 * pi );


	double k_squared = 4 * R * R_star / ( ( R + R_star )*( R + R_star ) + ( Z - Z_star )*( Z - Z_star ) );
	double k = ::sqrt( k_squared );

	double dR = ( R_star - R )/R;
	double dZ = ( Z_star - Z )/R;

	if ( ( ::fabs( dR ) > 1e-7 ) || ( ::fabs( dZ ) > 1e-7 ) )
	{
		answer *= ::sqrt( ( R + R_star )*( R + R_star ) + ( Z - Z_star )*( Z - Z_star ) ) * ( ( 1.0 - 0.5*k_squared )*boost::math::ellint_1( k ) - boost::math::ellint_2( k ) );
	}
	else
	{
		// k_squared = (1 + dR) / ( 1 + dR + dR^2/4 + dZ^2/4 ) ;
		// k_squared ~= (1 + dR)*( 1 - (dR + dR^2/4 + dZ^2/4) + (dR + dR^2/4 + dZ^2/4)^2))
		//   ~= ( 1 + dR ) * ( 1 - dR + .75 * dR^2 / 4 - dZ^2/4 );
		//   ~ = 1 - .75*dR^2 - .25*dZ^2;
		// Thus, m = 1 - dR*dR*.25 - dZ*dZ*.25;
		// K(m) ~= -(1/2) ln( 1 - m ) + ln(4) + O(|m-1|);
		// Set k = k^2 = 1 outside of the argument of K(). Keep linear terms, as doubles can detect those, but the squares are 
		// beyond double resolution.
		double K_approx;
		if ( dR == 0 )
		{
			K_approx = -::log( ::fabs( dZ ) ) - .5*::log( .25 ) + ::log( 4 );
		}
		else if ( dZ == 0 )
		{
			K_approx = -::log( ::fabs( dR ) ) - .5*::log( .25 ) + ::log( 4 );
		}
		else
		{
			if ( ::fabs( dR ) >= ::fabs( dZ ) )
			{
				double t = dZ/dR;
				K_approx = -::log( ::fabs( dR ) ) - .5*::log( 1 + t*t ) - .5*::log( .25 ) + ::log( 4 );
			}
			else
			{
				double t = dR/dZ;
				K_approx = -::log( ::fabs( dZ ) ) - .5*::log( 1 + t*t ) - .5*::log( .25 ) + ::log( 4 );
			}
		}

		answer *= R*( 1. + dR  ) * ( K_approx - 2.0 );
	}

	return answer;
}

const double xk[] = {0, 0.1943570033249354, 0.3772097381640342, 0.5391467053879678, 0.6742714922484358, 0.7806074389832003, 0.8595690586898966, 0.9148792632645746, 0.9513679640727469, 0.9739668681956774, 0.9870405605073769, 0.9940555066314021, 0.9975148564572244, 0.9990651964557858, 0.9996882640283532, 0.9999093846951440, 0.9999774771924616, 0.9999953160412205, 0.9999992047371147, 0.9999998927816124, 0.9999999888756649, 0.9999999991427051, 0.999999999952856, 0.999999999998232, 0.999999999999957};

const double wk[] = {0.1963495408493621,0.1904104648293382,0.173701844905907,0.1491828782311446,0.1207470724265376,0.09217973104519348,0.06638478442850675,0.04505767730866796,0.02875279931434859,0.01717776346664597,0.009548217946354038,0.004896875686700097,0.00229289587374098,0.0009678251282580301,0.0003628147184876642,0.0001187433505354336,0.00003327506421908961,7.81031990509301e-6,1.49796267039634e-6,2.282915074213832e-7,2.67890056961788e-8,2.335910283592051e-9,1.453895726781973e-10,6.172317347078991e-12,1.697723034317386e-13};

double BoundaryPsi( mfem::FiniteElementSpace *q_space, mfem::Vector & zero_soln, mfem::Vector const& r )
{

	mfem::Mesh *mesh = q_space->GetMesh();
	double Answer = 0;

	GridFunction q_fn;
	q_fn.MakeRef( q_space, zero_soln, 0 );

	for (int i = 0; i < mesh->GetNBE(); i++)
	{
		mfem::FaceElementTransformations *tr = mesh->GetBdrFaceTransformations(i);
		if (tr != NULL)
		{
			mfem::Vector favect;
			mfem::FiniteElement const  &bdr_cell = *q_space->GetFE(tr->Elem1No);
			double w, f_val;

			if ( tr->Elem2No >= 0 )
			{
				// Interior Face. Do Nothing.
				continue;
			}

			int dim = tr->Face->GetSpaceDim();
			int ndof_cell = bdr_cell.GetDof();

			mfem::Vector shape_f( ndof_cell );
			favect.SetSize( ndof_cell * dim );
			favect = 0.0;


			// Boundary face
			// check if it's the singular one.
			
			mfem::Vector a( 2 );
			mfem::Vector b( 2 );
			mfem::IntegrationPoint test_point;
			test_point.x = 0; test_point.y = 0; test_point.z = 0; test_point.weight = 0;
			tr->Face->Transform( test_point, a );
			test_point.x = 1; test_point.y = 0; test_point.z = 0; test_point.weight = 0;
			tr->Face->Transform( test_point, b );

			mfem::Vector n( 2 );
			n = b;
			n -= a;

			mfem::Vector x( 2 );
			x = r;
			x -= a;

			double s;
			if ( n( 0 ) == 0 )
				s = x( 1 ) / n( 1 );
			else
				s = x( 0 ) / n( 0 );

			if ( n.Norml2() == 0 )
				throw std::logic_error( "This is not on!!" );

			double SinhTanhEps = 1e-6;
			if ( ( x( 0 )*n( 1 ) == x( 1 )*n( 0 ) ) && ( s >= -1e-3 && s <= 1.001 ) )
			{
				if ( s >= SinhTanhEps && s <= 1.0 - SinhTanhEps )
				{
					double LeftInt = 0.0,RightInt = 0.0;
					int N_sinhtanh = 16;
					// 'left' integral is on [0,s]
					// 'right' integral is on [s,1]
					for ( int i = -N_sinhtanh; i <= N_sinhtanh; i++ )
					{
						mfem::IntegrationPoint ipL,ipR;
						ipL.y = 0; ipL.z = 0;
						ipR.y = 0; ipR.z = 0;
						if ( i < 0 )
						{
							ipL.x = s*( 1.0 - xk[ abs( i ) ] )/2.0; 
							ipR.x = ( 1.0 + s - ( 1.0 - s )*xk[ abs( i ) ] )/2.0; 
						}
						else
						{
							ipL.x = s*( 1.0 + xk[ abs( i ) ] )/2.0; 
							ipR.x = ( 1.0 + s + ( 1.0 - s )*xk[ abs( i ) ] )/2.0; 
						}
						// eip_L is inside the boundary cell, but is on the edge of it that corresponds
						// to the Boundary Face that we are integrating over
						mfem::IntegrationPoint eip_L;
						mfem::IntegrationPoint eip_R;
						tr->Loc1.Transform(ipL, eip_L);
						tr->Loc1.Transform(ipR, eip_R);
						
						w = wk[ abs( i ) ] / 2.0;
						mfem::Vector normal( dim );
						mfem::Vector r_star( dim );
						mfem::Vector q_value( 2 );
						// For ipL
						tr->Face->SetIntPoint(&ipL);
						mfem::CalcOrtho( tr->Face->Jacobian(), normal );
						tr->Face->Transform( ipL, r_star );

						f_val = GreensFunction( r, r_star );
						q_fn.GetVectorValue( tr->Elem1No, eip_L, q_value );
						LeftInt += ( w * s ) * f_val * ( q_value * normal );

						// For ipR
						tr->Face->SetIntPoint(&ipR);
						mfem::CalcOrtho( tr->Face->Jacobian(), normal );
						tr->Face->Transform( ipR, r_star );

						f_val = GreensFunction( r, r_star );
						q_fn.GetVectorValue( tr->Elem1No, eip_R, q_value );
						RightInt += ( w * ( 1.0 - s ) ) * f_val * ( q_value * normal );
					}
					Answer += LeftInt + RightInt;
				}
				else 
				{
					double tmp_ans = 0;
					int N_sinhtanh = 16;
					for ( int i = -N_sinhtanh; i <= N_sinhtanh; i++ )
					{
						mfem::IntegrationPoint ip;
						ip.y = 0; ip.z = 0;
						if ( i < 0 )
							ip.x = ( 1.0 - xk[ abs( i ) ] )/2.0; 
						else
							ip.x = ( 1.0 + xk[ abs( i ) ] )/2.0;
						// eip_L is inside the boundary cell, but is on the edge of it that corresponds
						// to the Boundary Face that we are integrating over
						mfem::IntegrationPoint eip_L;
						tr->Loc1.Transform(ip, eip_L);
						tr->Face->SetIntPoint(&ip);
						mfem::Vector normal( dim );
						// Normal contains the factor of det(J) that arises in the change of 
						// variables in the integration
						mfem::CalcOrtho( tr->Face->Jacobian(), normal );
						mfem::Vector r_star( dim );
						tr->Face->Transform( ip, r_star );

						f_val = GreensFunction( r, r_star );

						w = wk[ abs( i ) ] / 2.0;

						mfem::Vector q_value( 2 );
						q_fn.GetVectorValue( tr->Elem1No, eip_L, q_value );

						tmp_ans += w * f_val * ( q_value * normal );

					}

					Answer += tmp_ans;
				}
			}
			else
			{
				// Not singular, just use MFEM's gaussian quadrature.
				const mfem::IntegrationRule *ir = nullptr;
				int order = 2 * bdr_cell.GetOrder() + 2;
				if (bdr_cell.GetMapType() == FiniteElement::VALUE)
				{
					order += tr->Face->OrderW();
				}
				ir = &IntRules.Get(tr->FaceGeom, order);

				double tmp_ans = 0.0;

				for (int p = 0; p < ir->GetNPoints(); p++)
				{
					const mfem::IntegrationPoint &ip = ir->IntPoint(p);


					// eip_L is inside the boundary cell, but is on the edge of it that corresponds
					// to the Boundary Face that we are integrating over
					mfem::IntegrationPoint eip_L;
					tr->Loc1.Transform(ip, eip_L);

					tr->Face->SetIntPoint(&ip);
					mfem::Vector normal( dim );
					// Normal contains the factor of det(J) that arises in the change of 
					// variables in the integration
					mfem::CalcOrtho( tr->Face->Jacobian(), normal );
					mfem::Vector r_star( dim );
					tr->Face->Transform( ip, r_star );

					f_val = GreensFunction( r, r_star );

					w = ip.weight;

					mfem::Vector q_value( 2 );
					q_fn.GetVectorValue( tr->Elem1No, eip_L, q_value );

					tmp_ans += w * f_val * ( q_value * normal );

				}

				Answer += tmp_ans;
			}
		}
	}

	return Answer;
}

double GreensFunctionPsi( mfem::Mesh * mesh, mfem::Vector r, std::function<double( const mfem::Vector& )> const& j_coil )
{
	double Answer = 0;
	for (int i = 0; i < mesh->GetNE(); i++)
	{
		mfem::ElementTransformation *T = mesh->GetElementTransformation( i );
		mfem::IntegrationRule const *ir = &IntRules.Get(T->GetGeometryType(), 16);
		double CellAnswer = 0;
		for (int p = 0; p < ir->GetNPoints(); p++)
		{
			mfem::Vector pt( 3 );
			const mfem::IntegrationPoint &ip = ir->IntPoint(p);
			T->Transform( ip, pt );
			if ( ::fabs( j_coil( pt ) ) < 1e-9 )
				continue;
			double tmp = ip.weight  * T->Weight() * GreensFunction( r, pt ) * j_coil( pt );

			CellAnswer += tmp;
		}

		Answer += CellAnswer;
	}
	return Answer;
}

double PsiFromZeroBC( mfem::FiniteElementSpace *q_space, mfem::FiniteElementSpace *u_space, mfem::Vector & zero_soln, mfem::Vector const& r )
{

	mfem::Mesh *mesh = q_space->GetMesh();
	double Answer = 0;

	GridFunction q_fn,u_fn;
	q_fn.MakeRef( q_space, zero_soln, 0 );
	u_fn.MakeRef( u_space, zero_soln, q_space->GetVSize() );

	Answer  = BoundaryPsi( q_space, zero_soln, r );

	double psiVal = 0.0;

	{
		mfem::DenseMatrix ptMat( 2, 1 );
		ptMat( 0, 0 ) = r( 0 );
		ptMat( 1, 0 ) = r( 1 );
		mfem::Array<int> elem_ids( 1 );
		mfem::Array< mfem::IntegrationPoint > ips( 1 );
		mesh->FindPoints( ptMat, elem_ids, ips, false );

		psiVal = u_fn.GetValue( elem_ids[ 0 ], ips[ 0 ] );
	}

	return Answer + psiVal;
}

double GreensFunctionBoundaryCoefficient::Eval( mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip )
{
	if ( T.GetGeometryType() != mfem::Geometry::Type::SEGMENT )
		throw new std::logic_error( "This only works with Segements!" );
	mfem::Vector pt( 3 );
	T.Transform( ip, pt );
	return 2*BoundaryPsi( Q_space, psi_hat, pt );
}


}
