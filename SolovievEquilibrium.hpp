
class SolovievEq {

	private:
		static double psi_P( double A, double C, double R, double z )
		{
			double R4 = R*R*R*R;

			return ( C+A )*( R4/8.0 ) + A*( ( R*R/2.0 )*::log( R ) - R4/8.0 );
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
			return R*R*( R*R - 4 * Z*Z );
		};



		double A,C;
		double c1,c2,c3,c4;
	public:
		SolovievEq( double A_in, double C_in, double c1_in, double c2_in, double c3_in, double c4_in ) 
		{
			A = A_in;
			C = C_in;
			c1 = c1_in;
			c2 = c2_in;
			c3 = c3_in;
			c4 = c4_in;
		}
		~SolovievEq();

		double psi( const mfem::Vector& pt )
		{
			double R(pt(0));
			double Z(pt(1));
			return psi_P( A, C, R, Z ) + c1 * psi_1( R, Z ) + c2 * psi_2( R, Z ) + c3 * psi_3( R, Z ) + c4 * psi_4( R, Z );
		}

		void q( const mfem::Vector &pt, mfem::Vector &q )
