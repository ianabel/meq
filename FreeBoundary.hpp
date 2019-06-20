

double GreensFunction( mfem::Vector const& r, mfem::Vector const& r_star )
{
	double R = r[ 0 ];
	double R_star = r_star[ 0 ];
	double Z = r[ 1 ];
	double Z_star = r_star[ 1 ];

	double k_squared = 4 * R * R_star / ( ( R + R_star )*( R + R_star ) + ( Z - Z_star )*( Z - Z_star ) );

	double answer = 1.0 / M_PI;

	answer *= ::sqrt( ( R*R_star + Z*Z_star ) / ( k_squared ) );
	
	double k = ::sqrt( k_squared );
	answer *= ( ( 1.0 - 0.5*k_squared )*std::comp_ellint_1( k ) - std::comp_ellint_2( k ) );
	return answer;
}

double BoundaryPsi( FiniteElementSpace *q_space, mfem::Vector & zero_soln, mfem::Vector const& r )
{

	LinerForm *lf = new LinearForm( q_space );

	auto bdGF = std::bind( GreensFunction, std::placeholders::_1, r );
	StdFunctionCoefficient GreensFunctionCoefficient( bdGF );
	lf.AddBoundaryIntegrator( new VectorBoundaryFluxLFIntegrator( GreensFunctionCoefficient ) );

	GridFunction gradPsi;
	gradPsi.MakeRef( q_space, zero_soln, 0 );
	return lf( gradPsi );
	
}

class GreensFunctionBoundaryCoefficient : public mfem::Coefficient
{
	protected:
		mfem::Mesh const *mesh;
		mfem::FiniteElementSpace *Q_space;
		mfem::Vector & psi_hat;
	public:
		GreensFunctionBoundaryCoefficient( mfem::Mesh const* mesh_r, FiniteElementSpace *q_space, mfem::Vector & zero_soln )
			: mesh( mesh_r ), Q_space( q_space ), psi_hat( zero_soln )
		{
		}

		virtual override double Eval( mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip )
		{
			if ( T.GetGeometryType() != mfem::Geometry::Type::SEGMENT )
				throw new std::logic_error( "This only works with Segements!" );
			mfem::Vector pt( 3 );
			T.Transform( ip, pt );
			return BoundaryPsi( Q_space, psi_hat, pt );
		}
		
}
