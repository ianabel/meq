
// Analytic manufactured solution
// from <citation>

class TSVSoln1 {

	private:
		double r0,kr,kz;
	public:
		TSVSoln1( double r, double k_r, double k_z )
			: r0( r ), kr( k_r ), kz( k_z )
		{
		}
		~TSVSoln1() {};

		double Psi( const mfem::Vector& pt ) const
		{
			double r(pt(0));
			double z(pt(1));
			return ::sin( kr*( r + r0 ) )*::cos( kz*z );
		}

		void GradPsi( const mfem::Vector &pt, mfem::Vector &q ) const
		{
			double r(pt(0));
			double z(pt(1));
			
			q.SetSize( 2 );
			q( 0 ) =  kr*::cos( kr*( r + r0 ) )*::cos( kz*z );
			q( 1 ) = -kz*::sin( kr*( r + r0 ) )*::sin( kz*z );
		}

		double operator()( mfem::Vector const& pt, double psi ) const
		{
			double r = pt( 0 );
			double z = pt( 1 );
			double sc = ::sin( kr*( r + r0 ) )*::cos( kz*z );
			return ( ( kr*kr + kz*kz )*psi + ( kr/r )*::cos( kr*( r + r0 ) )*::cos( kz*z ) 
				+ r * ( sc*sc - psi*psi + ::exp( -sc ) - ::exp( -psi ) ) )/r;
		}
};


