#include <functional>
#include "StdFnCoeffs.hpp"

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

	if ( dR == 0 && dZ == 0 )
		throw std::logic_error( "Logarithmic divergence detected. Caution is advised." );

	if ( ( ::fabs( dR ) > 1e-7 ) || ( ::fabs( dZ ) > 1e-7 ) )
	{
		answer *= ::sqrt( ( R + R_star )*( R + R_star ) + ( Z - Z_star )*( Z - Z_star ) ) * ( ( 1.0 - 0.5*k_squared )*std::comp_ellint_1( k ) - std::comp_ellint_2( k ) );
		if ( !std::isfinite( answer ) )
			throw std::logic_error( "Wat" );
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

using namespace mfem;

const unsigned int N = 20;
const double xk[] = {0, 0.1943570033249354, 0.3772097381640342, 0.5391467053879678, 0.6742714922484358, 0.7806074389832003, 0.8595690586898966, 0.9148792632645746, 0.9513679640727469, 0.9739668681956774, 0.9870405605073769, 0.9940555066314021, 0.9975148564572244, 0.9990651964557858, 0.9996882640283532, 0.9999093846951440, 0.9999774771924616, 0.9999953160412205, 0.9999992047371147, 0.9999998927816124, 0.9999999888756649, 0.9999999991427051, 0.999999999952856, 0.999999999998232, 0.999999999999957};

const double wk[] = {0.1963495408493621,0.1904104648293382,0.173701844905907,0.1491828782311446,0.1207470724265376,0.09217973104519348,0.06638478442850675,0.04505767730866796,0.02875279931434859,0.01717776346664597,0.009548217946354038,0.004896875686700097,0.00229289587374098,0.0009678251282580301,0.0003628147184876642,0.0001187433505354336,0.00003327506421908961,7.81031990509301e-6,1.49796267039634e-6,2.282915074213832e-7,2.67890056961788e-8,2.335910283592051e-9,1.453895726781973e-10,6.172317347078991e-12,1.697723034317386e-13};


double SinhTanhQuad( double a, double b, std::function<double( double )> F )
{
	double L = ( b - a ) /2.;
	double quad_ans = 0;
	quad_ans = F( L*( 1. + xk[ 0 ] ) + a )*wk[ 0 ];
	for ( unsigned int i=1; i < N; i++ )
	{
		quad_ans += wk[ i ]*( F( L*( 1. + xk[ i ] ) + a ) + F( L*( 1. - xk[ i ] ) + a ) );
	}
	return quad_ans * L;
};


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


			

			b -= a;

			mfem::Vector x( 2 );
			x = r;
			x -= a;

			double s;
			if ( b( 0 ) == 0 )
				s = x( 1 ) / b( 1 );
			else
				s = x( 0 ) / b( 0 );

			if ( b.Norml2() == 0 )
				throw std::logic_error( "This is not on!!" );

			if ( ( x( 0 )*b( 1 ) == x( 1 )*b( 0 ) ) && ( s >= -1e-3 && s <= 1.001 ) )
			{

				
				
				// Takes values in [0,1]
				auto IntFunc = [&]( double y ) {

					IntegrationPoint ip_x;
					ip_x.x = y; ip_x.y = 0; ip_x.z = 0;
					ip_x.weight = 1;
					tr->Face->SetIntPoint(&ip_x);
				
					mfem::Vector normal( 2 );
					mfem::CalcOrtho( tr->Face->Jacobian(), normal );
				
					mfem::IntegrationPoint eip_L;
					tr->Loc1.Transform(ip_x, eip_L);
					
					mfem::Vector r_star( 2 );

					tr->Face->Transform( ip_x, r_star );
					double GF_val = GreensFunction( r, r_star );
					
					mfem::Vector q_value( 2 );
					q_fn.GetVectorValue( tr->Elem1No, eip_L, q_value );

					return GF_val * ( q_value * normal );
				};
					

				if ( s > 0.0 && s < 1.0 )
				{
					double ExpInt;
					double LInt = SinhTanhQuad( 0.0, s, IntFunc );
					double RInt = SinhTanhQuad( s, 1.0, IntFunc );
					ExpInt = LInt + RInt;
					Answer += ExpInt;
				}
				else if ( s <= 0 || s >= 1 )
				{
					double ExpInt = SinhTanhQuad( 0.0, 1, IntFunc );
					Answer += ExpInt;
				}
				else {
					throw std::logic_error( "Kapow" );
				}
			}
			else 
			{
				const mfem::IntegrationRule *ir = nullptr;
				int order = 2 * bdr_cell.GetOrder() + 10;
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

class GreensFunctionBoundaryCoefficient : public mfem::Coefficient
{
	protected:
		mfem::Mesh const *mesh;
		mfem::FiniteElementSpace *Q_space;
		mfem::Vector & psi_hat;
	public:
		GreensFunctionBoundaryCoefficient( mfem::Mesh const* mesh_r, mfem::FiniteElementSpace *q_space, mfem::Vector & zero_soln )
			: mesh( mesh_r ), Q_space( q_space ), psi_hat( zero_soln )
		{
		}

		virtual double Eval( mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip )
		{
			if ( T.GetGeometryType() != mfem::Geometry::Type::SEGMENT )
				throw new std::logic_error( "This only works with Segements!" );
			mfem::Vector pt( 3 );
			T.Transform( ip, pt );
			return BoundaryPsi( Q_space, psi_hat, pt );
		}
		
};

