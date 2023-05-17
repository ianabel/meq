#include <cmath>
double PlasmaModel::Density( PlasmaModel::Index s, double R, double z )
{
	double PsiVal = Psi( R , z );
	double T =  Temperature( s, R, z );
	double Te = Temperature( ElectronIdx, R, z ); 
	double Zs = Content[ s ].Z;
	double m = Content[ s ].mass;
	double omega = Omega( R, z );
	return N( s, PsiVal ) * std::exp( -Zs * Te * phi0( R, z )/ T + m*omega*omega*R*R / ( 2.0 * T )  );
}	

double PlasmaModel::Temperature( PlasmaModel::Index s, double R, double z )
{
	return Content[ s ].Temperature( Psi( R, z ) );
}

double PlasmaModel::Omega( double R, double z )
{
	return Omega_psi( Psi( R, z ) );
}

// Returns e phi0 / T_e
double phi0( double R, double z )
{
	double PsiVal = Psi( R, z );
	double m = Content[ IonIdx ].mass;
	double Z = Content[ IonIdx ].Z;
	double omega = Omega( R, z );
	double Ti = Temperature( IonIdx, R, z ), Te = Temperature( ElectronIdx, R, z );
	return ( 1.0 / ( 1.0 + Z*Te/Ti ) ) * ( m * omega * omega * R * R )/( 2.0 * Ti );
}

double N( Index s, double psi )
{
	double T =  Content[ s ].Temperature( s, psi );
	double Te = Content[ s ].Temperature( ElectronIdx, psi ); 
	double Zs = Content[ s ].Z;
	double m = Content[ s ].mass;
	double omega = Omega_psi( psi );
	double R = RMid_psi( psi );
	return Content[ s ].MidplaneDensity( psi ) * std::exp( Zs * Te * phi0( R, z )/T - m*omega*omega*R*R/( 2.0 * T );

}
