
// Analytic Equilibrium with dissimilar source functions 
// from <citation>

class McCarthyEquilibrium {

	private:
		double T;
		std::vector<double> c;
	public:
		McCarthyEquilibrium( double T_in, std::initializer_list<double> l )
			: T( T_in ), c( l )
		{
			if ( c.size() != 18 )
				throw std::invalid_argument( "Needs precisely 18 free parameters" );
		}
		~McCarthyEquilibrium() {};

		double Psi( const mfem::Vector& pt ) const
		{
			double r(pt(0));
			double z(pt(1));
			double p = ::sqrt( T );
			double q = p/2.0;
			double nu = ::sqrt( .75 )*p;

			double s = ::sqrt( r*r + z*z );

			double psi_ret;

			psi_ret = c[ 0 ] + c[ 1 ]*r*r + r*std::cyl_bessel_j( 1, p*r )*( c[ 2 ] + c[ 3 ]*z ) + c[ 4 ]*::cos( p*z ) + c[ 5 ]*::sin( p*z )
				+ r*r*( c[ 6 ]*::cos( p*z ) + c[ 7 ]*::sin( p*z ) ) + c[ 8 ]*::cos( p*s ) + c[ 9 ]*::sin( p*s ) 
				+ r*std::cyl_bessel_j( 1, nu*r )*( c[ 10 ]*::cos( q*z ) + c[ 11 ]*::sin( q*z ) ) + r*std::cyl_bessel_j( 1, q*r )*( c[ 12 ]*::cos( nu*z ) + c[ 13 ]*::sin( nu*z ) )
				+ r*std::cyl_neumann( 1, nu*r )*( c[ 14 ]*::cos( q*z ) + c[ 15 ]*::sin( q*z ) ) + r*std::cyl_neumann( 1, q*r )*( c[ 16 ]*::cos( nu*z ) + c[ 17 ]*::sin( nu*z ) );

			return psi_ret;
		}

		void GradPsi( const mfem::Vector &pt, mfem::Vector &q_vec ) const
		{
			double r(pt(0));
			double z(pt(1));
			double p = ::sqrt( T );
			double q = p/2.0;
			double nu = ::sqrt( .75 )*p;

			double s = ::sqrt( r*r + z*z );
			q_vec.SetSize( 2 );
			q_vec( 0 ) = 0;
			q_vec( 1 ) = 0;

			// d( r J_1(p*r) ) / dr 
			double drJpdr = std::cyl_bessel_j( 1, p*r ) + ( p*r/2.0 )*( std::cyl_bessel_j( 0, p*r ) - std::cyl_bessel_j( 2, p*r ) );
			double drJqdr = std::cyl_bessel_j( 1, q*r ) + ( q*r/2.0 )*( std::cyl_bessel_j( 0, q*r ) - std::cyl_bessel_j( 2, q*r ) );
			double drJnudr = std::cyl_bessel_j( 1, nu*r ) + ( nu*r/2.0 )*( std::cyl_bessel_j( 0, nu*r ) - std::cyl_bessel_j( 2, nu*r ) );

			double drYqdr = std::cyl_neumann( 1, q*r ) + ( q*r/2.0 )*( std::cyl_neumann( 0, q*r ) - std::cyl_neumann( 2, q*r ) );
			double drYnudr = std::cyl_neumann( 1, nu*r ) + ( nu*r/2.0 )*( std::cyl_neumann( 0, nu*r ) - std::cyl_neumann( 2, nu*r ) );

			// d psi / d r
			q_vec( 0 ) = 2*c[ 1 ]*r + drJpdr * ( c[ 2 ] + c[ 3 ]*z ) + 2*r*( c[ 6 ]*::cos( p*z ) + c[ 7 ]*::sin( p*z ) ) - ( c[ 8 ] * p * r / s )*::sin( p*s ) + ( c[ 9 ] * p * r / s ) * ::sin( p*s ) 
				+ drJnudr*( c[ 10 ]*::cos( q*z ) + c[ 11 ]*::sin( q*z ) ) + drJqdr*( c[ 12 ]*::cos( nu*z ) + c[ 13 ]*::sin( nu*z ) )
				+ drYnudr*( c[ 14 ]*::cos( q*z ) + c[ 15 ]*::sin( q*z ) ) + drYqdr*( c[ 16 ]*::cos( nu*z ) + c[ 17 ]*::sin( nu*z ) );

			q_vec( 1 ) = r*std::cyl_bessel_j( 1, p*r )*c[ 3 ] - p*c[ 4 ]*::sin( p*z ) + p*c[ 5 ]*::cos( p*z )
				+ r*r*( -p * c[ 6 ]*::sin( p*z ) + p * c[ 7 ]*::cos( p*z ) ) - ( c[ 8 ] * p * z / s )*::sin( p*s ) + ( c[ 9 ] * p * z / s )*::cos( p*s ) 
				+ r*std::cyl_bessel_j( 1, nu*r )*( -q * c[ 10 ]*::sin( q*z ) + q * c[ 11 ]*::cos( q*z ) ) + r*std::cyl_bessel_j( 1, q*r )*( -nu*c[ 12 ]*::sin( nu*z ) + nu * c[ 13 ]*::cos( nu*z ) )
				+ r*std::cyl_neumann( 1, nu*r )*( -q * c[ 14 ]*::sin( q*z ) + q * c[ 15 ]*::cos( q*z ) ) + r*std::cyl_neumann( 1, q*r )*(  -nu * c[ 16 ]*::sin( nu*z ) + nu * c[ 17 ]*::cos( nu*z ) );

		}

		double Pprime( double psi ) const
		{
			return -c[ 0 ]*T;
		}

		double FFprime( double psi ) const
		{
			return T*psi - c[ 1 ];
		}

		double operator()( mfem::Vector const& pt, double psi ) const
		{
			double R = pt( 0 );
			return T*( psi - c[ 0 ] - c[ 1 ]*R*R )/R;
		}
};


