
// Test harness for 
// integrating f(x) * K(x) over [0,1], where K is the compelte
// elliptic integral of the first kind
//


unsigned int N_GL = 5;
constexpr double xk[] = {
5.6522282050800971359e-3,
7.3430371742652273406e-2,
2.8495740446255815371e-1,
6.1948226408477838141e-1,
9.1575808300469833378e-1
};
constexpr double wk[] = {
2.1046945791854629119e-2,
1.3070554074444669759e-1,
2.8970230167131415684e-1,
3.5022037012039871029e-1,
2.0832484167198580616e-1
};
#include <iomanip>
#include <boost/math/special_functions/ellint_1.hpp>

int main( int argc, char **argv )
{
	// Uses the extended Gauss-log quadrature
	// to compute Integrals of K(x)
	double ans = 0.0;
	auto F = []( double x ) {
		return boost::math::ellint_1( ::sqrt( 1 - x ) );
	};
	for ( unsigned int i=0; i<N_GL; i++ )
	{
		ans += wk[ i ] * F( xk[ i ] );
	}
	std::cerr << std::setprecision( 16 ) << ans << std::endl;
}
