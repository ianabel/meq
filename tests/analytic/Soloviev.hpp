
class SolovievEquilibrium {

	private:
		static double psi_P( double A, double C, double R, double z )
		{
			double R4 = R*R*R*R;

			return C*R4/8.0 + A*( ( R*R/2.0 )*::log( R ) );
		}
		static double psi_1( double R, double z ) {
			return 1;
		};
		static double psi_2( double R, double z ) {
			return R*R;
		};
		static double psi_3( double R, double z ) {
			return z*z - R*R*::log( R );
		};
		static double psi_4( double R, double z ) {
			return R*R*( R*R - 4 * z*z );
		};

		static void grad_psi_P( double A, double C, double R, double z, mfem::Vector& grad )
		{
			grad( 0 ) = C * R * R * R/2.0 + A * R * ( ::log( R ) + 1/2.0 );
			grad( 1 ) = 0.0;
		}

		static void grad_psi_1( double R, double z, mfem::Vector& grad ) {
			grad( 0 ) = 0.0;
			grad( 1 ) = 0.0;
		};
		static void grad_psi_2( double R, double z, mfem::Vector& grad ) {
			grad( 0 ) = 2.0*R;
			grad( 1 ) = 0.0;
		};
		static void grad_psi_3( double R, double z, mfem::Vector& grad ) {
			grad( 0 ) = - R*( 2.0*::log( R ) + 1.0 );
			grad( 1 ) = 2*z;
		};
		static void grad_psi_4( double R, double z, mfem::Vector& grad ) {
			grad( 0 ) = 4.0 * R *( R*R - 2.0*z*z );
			grad( 1 ) = - 8.0 * R * R * z;
		};


		double A,C;
		double c1,c2,c3,c4;
	public:
		SolovievEquilibrium( double A_in, double C_in, double c1_in, double c2_in, double c3_in, double c4_in ) 
		{
			A = A_in;
			C = C_in;
			c1 = c1_in;
			c2 = c2_in;
			c3 = c3_in;
			c4 = c4_in;
		}
		~SolovievEquilibrium() {};

		double Psi( const mfem::Vector& pt ) const
		{
			double R(pt(0));
			double Z(pt(1));
			return psi_P( A, C, R, Z ) + c1 * psi_1( R, Z ) + c2 * psi_2( R, Z ) + c3 * psi_3( R, Z ) + c4 * psi_4( R, Z );
		}

		void GradPsi( const mfem::Vector &pt, mfem::Vector &q ) const
		{
			double R(pt(0));
			double Z(pt(1));
			mfem::Vector tmp( 2 );
			q.SetSize( 2 );
			q( 0 ) = 0;
			q( 1 ) = 0;

			grad_psi_P( A, C, R, Z, tmp );
			q += tmp;

			grad_psi_1( R, Z, tmp );
			tmp *= c1;
			q += tmp;

			grad_psi_2( R, Z, tmp );
			tmp *= c2;
			q += tmp;

			grad_psi_3( R, Z, tmp );
			tmp *= c3;
			q += tmp;

			grad_psi_4( R, Z, tmp );
			tmp *= c4;
			q += tmp;
		}

		double Pprime( double psi ) const
		{
			return C;
		}

		double FFprime( double psi ) const
		{
			return A;
		}

		double operator()( const mfem::Vector &pt, double psi ) const
		{
			double R = pt( 0 );
			return -( R*R*Pprime( psi ) + FFprime( psi ) )/R;
		};
};

