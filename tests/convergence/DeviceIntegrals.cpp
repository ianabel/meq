/*
 * THE POINT ROUTES UNDER A DEVICE (MEASUREMENTS.md M-DA-5).
 *
 * addSourceColumns() -- and so assembleSourceDerivative() and every q-driven
 * pin row built on it -- integrateOverDomain(), deliveredPlasmaCurrent() and
 * sourceTableReach() take the potential at every quadrature point from a
 * kernel and scatter on the device, evaluating their integrand on the host.
 * Each against its host loop on one solved NormalisedMHDSource equilibrium,
 * route on and route off: the same numbers to the round-off of the two
 * geometries (the tables' radius and weight against the host's Transform()
 * and Weight(), the same to the bit in double, M-322), and the device route
 * taken. Needs a Device (MEQ_TEST_DEVICE) and a build with device kernels,
 * and skips loudly where either is missing.
 */
#define BOOST_TEST_MODULE DeviceIntegrals
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "convergence/ConvergenceHarness.hpp"
#include "unit/DeviceKernels.hpp"

namespace
{
	struct DeviceFromEnvironment
	{
		DeviceFromEnvironment()
		{
			char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
			if ( requested != nullptr && *requested != '\0' )
				static mfem::Device device( requested );
		}
	};

	double const eps = std::numeric_limits<mfem::real_t>::epsilon();

	/// The q-driven shape: p' constant, gg' linear in Psi, on the standard box.
	struct Solved
	{
		mfem::Mesh mesh;
		meq::NormalisedMHDSource source;
		mfem::ConstantCoefficient zero;
		mfem::FunctionCoefficient guess;
		meq::GradShafranovSolver solver;

		explicit Solved( int cells )
			: mesh( meq::tests::makeMesh( meq::tests::standardBox(), cells ) ),
			  source( std::make_shared<meq::ConstantProfile const>( 0.45 ),
			          std::make_shared<meq::SplineProfile const>( std::vector<meq::Knot>{
			              { 0.0, 0.55*2.20, 0.55*0.55 }, { 1.0, 0.55*( 2.20 + 0.55 ), 0.55*0.55 } } ),
			          0.30, 1.0 ),
			  zero( 0.0 ),
			  guess( []( mfem::Vector const &x )
			  {
				  meq::tests::Rectangle const box = meq::tests::standardBox();
				  return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				         *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			  } ),
			  solver( mesh, 2 )
		{
			solver.setAssemblyMode( meq::GradShafranovSolver::AssemblyMode::Batched );
			solver.setLocalFactorMode( meq::GradShafranovSolver::LocalFactorMode::Batched );
			solver.setSource( source, 0.30 );
			solver.setBoundaryData( zero );
			solver.setInitialGuess( guess );
			solver.setNewtonControl( 1.0e-12, 1.0e-14, 60 );
			solver.solve();
		}
	};

	bool deviceConfigured( char const *name )
	{
		if ( mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return true;
		std::printf( "\n  %s: SKIPPED -- no Device is configured; MEQ_TEST_DEVICE=debug or cuda runs it\n",
		             name );
		return false;
	}

	double relative( double a, double b )
	{
		return std::abs( a - b )/std::max( std::abs( b ), std::numeric_limits<double>::min() );
	}
}

BOOST_TEST_GLOBAL_FIXTURE( DeviceFromEnvironment );

/*
 * THE COLUMNS: assembleSourceDerivative() of a function of psi alone -- the
 * shape of the q-driven pin rows' integrands -- and of one reading R and z.
 */
BOOST_AUTO_TEST_CASE( theSourceColumnsAreTheHostLoops, *meq::tests::requiresDeviceKernels() )
{
	if ( !deviceConfigured( "theSourceColumnsAreTheHostLoops" ) )
		return;
	Solved s( 8 );
	mfem::Vector const state( s.solver.state() );
	auto column = [ & ]( bool device, std::function<double( double, double, double )> const &f )
	{
		s.solver.setDevicePointRoutes( device );
		long const before = s.solver.devicePointRoutesTaken();
		mfem::Vector out;
		s.solver.assembleSourceDerivative( state, f, out );
		BOOST_TEST( ( s.solver.devicePointRoutesTaken() > before ) == device,
		            "the " << ( device ? "device" : "host" ) << " route was not the one taken" );
		out.HostRead();
		return out;
	};
	std::vector<std::function<double( double, double, double )>> const integrands {
		[]( double, double, double psi ) { return psi > 0.0 ? std::sqrt( psi ) : 0.0; },
		[]( double radius, double z, double psi ) { return radius*z + psi*psi; } };
	for ( auto const &g : integrands )
	{
		mfem::Vector const host = column( false, g );
		mfem::Vector const device = column( true, g );
		double worst = 0.0, scale = 0.0;
		for ( int i = 0; i < host.Size(); ++i )
		{
			worst = std::max( worst, std::abs( double( host( i ) - device( i ) ) ) );
			scale = std::max( scale, std::abs( double( host( i ) ) ) );
		}
		std::printf( "\n  a source column, device against host: %.2e of its largest entry %.3e\n",
		             worst/scale, scale );
		BOOST_TEST( scale > 0.0 );
		BOOST_TEST( worst <= 64.0*eps*scale, "the device route's column differs by " << worst/scale
		            << " of its largest entry" );
	}
	s.solver.setDevicePointRoutes( true );
}

/*
 * THE INTEGRALS: the delivered flux's quadrature, the delivered current and
 * the table reach.
 */
BOOST_AUTO_TEST_CASE( theDeliveredIntegralsAreTheHostLoops, *meq::tests::requiresDeviceKernels() )
{
	if ( !deviceConfigured( "theDeliveredIntegralsAreTheHostLoops" ) )
		return;
	Solved s( 8 );
	double const axis = s.solver.psiAxis();
	double const boundary = s.solver.psiBoundary();
	auto flux = [ axis, boundary ]( double radius, double, double psi )
	{
		double const x = ( psi - boundary )/( axis - boundary );
		if ( !( radius > 0.0 ) || !( x > 0.0 ) )
			return 0.0;
		double const value = 4.84 + 2.42*x + 0.3025*x*x;
		return std::sqrt( value )/radius;
	};
	struct Read
	{
		double flux, current, low, high;
		long taken;
	};
	auto read = [ & ]( bool device )
	{
		s.solver.setDevicePointRoutes( device );
		long const before = s.solver.devicePointRoutesTaken();
		Read out{};
		out.flux = s.solver.integrateOverDomain( flux );
		out.current = s.solver.deliveredPlasmaCurrent();
		meq::TableReach const reach = s.solver.sourceTableReach();
		out.low = reach.below;
		out.high = reach.above;
		out.taken = s.solver.devicePointRoutesTaken() - before;
		return out;
	};
	Read const host = read( false );
	Read const device = read( true );
	s.solver.setDevicePointRoutes( true );
	std::printf( "\n  device against host: Phi %.12g against %.12g (%.2e), mu0 I_p %.12g against %.12g "
	             "(%.2e); %ld device routes taken, %ld on the host route\n",
	             device.flux, host.flux, relative( device.flux, host.flux ), device.current, host.current,
	             relative( device.current, host.current ), device.taken, host.taken );
	BOOST_TEST( host.taken == 0L );
	BOOST_TEST( device.taken == 3L, "of the three routes, " << device.taken << " ran on the device" );
	BOOST_TEST( relative( device.flux, host.flux ) <= 64.0*eps );
	BOOST_TEST( relative( device.current, host.current ) <= 64.0*eps );
	BOOST_TEST( std::abs( device.low - host.low ) <= 64.0*eps );
	BOOST_TEST( std::abs( device.high - host.high ) <= 64.0*eps );
}

/*
 * THE INTEGRANDS AS DATA: the delivered flux's g and the three functions of
 * it the q-driven pin rows differentiate into, carried as a polynomial or as
 * a spline's knots (hermiteKnots()), integrated and scattered in the kernel.
 * Each against its host loop -- the same arithmetic on the host -- and the
 * column against assembleSourceDerivative() of the same function as a
 * lambda, which is how the pin rows write it today.
 */
BOOST_AUTO_TEST_CASE( theIntegrandsAsDataAreTheHostLoops, *meq::tests::requiresDeviceKernels() )
{
	if ( !deviceConfigured( "theIntegrandsAsDataAreTheHostLoops" ) )
		return;
	Solved s( 8 );
	mfem::Vector const state( s.solver.state() );
	std::vector<double> const g2 { 4.84, 2.42, 0.3025 };
	// G as a spline: a quadratic is reproduced by its Hermite knots exactly.
	class Quadratic : public meq::Profile
	{
		public:
			explicit Quadratic( std::vector<double> c ) : c( std::move( c ) ) {}
			double operator()( double x ) const override { return c[ 0 ] + x*( c[ 1 ] + x*c[ 2 ] ); }
			double prime( double x ) const override { return c[ 1 ] + 2.0*x*c[ 2 ]; }
			double doublePrime( double ) const override { return 2.0*c[ 2 ]; }
		private:
			std::vector<double> c;
	};
	Quadratic const exact( g2 );
	meq::SplineProfile const square( meq::hermiteKnots( exact, { 0.0, 0.25, 0.5, 0.75, 1.0 } ) );
	meq::SplineProfile const hat( std::vector<meq::Knot>{ { 0.0, 0.0, 0.0 }, { 0.5, 1.0, 0.0 },
	                                                     { 1.0, 0.0, 0.0 } } );

	// The delivered flux: polynomial and spline, device against host.
	meq::PsiIntegrand polynomial;
	polynomial.polynomial = g2;
	meq::PsiIntegrand spline;
	spline.square = &square;
	auto flux = [ & ]( meq::PsiIntegrand const &f, bool device )
	{
		s.solver.setDevicePointRoutes( device );
		long const before = s.solver.devicePointRoutesTaken();
		double const value = s.solver.integratePsiIntegrand( f );
		BOOST_TEST( ( s.solver.devicePointRoutesTaken() > before ) == device );
		return value;
	};
	double const polyHost = flux( polynomial, false ), polyDevice = flux( polynomial, true );
	double const splineHost = flux( spline, false ), splineDevice = flux( spline, true );
	std::printf( "\n  Phi: polynomial %.12g device against %.12g host (%.2e); spline %.12g against %.12g "
	             "(%.2e); the two G's %.2e apart\n", polyDevice, polyHost, relative( polyDevice, polyHost ),
	             splineDevice, splineHost, relative( splineDevice, splineHost ),
	             relative( splineHost, polyHost ) );
	BOOST_TEST( relative( polyDevice, polyHost ) <= 64.0*eps );
	BOOST_TEST( relative( splineDevice, splineHost ) <= 64.0*eps );
	BOOST_TEST( relative( splineHost, polyHost ) <= 64.0*eps, "the Hermite spline of G is not G" );

	// The four kinds as columns, device against host and against the lambda.
	double const axis = s.solver.psiAxis(), boundary = s.solver.psiBoundary(), span = axis - boundary;
	for ( meq::PsiIntegrandKind kind : { meq::PsiIntegrandKind::Root, meq::PsiIntegrandKind::SlopeOverRoot,
	                                     meq::PsiIntegrandKind::PsiSlopeOverRoot,
	                                     meq::PsiIntegrandKind::RatioOverRoot } )
	{
		meq::PsiIntegrand f;
		f.kind = kind;
		f.square = &square;
		f.numerator = &hat;
		auto column = [ & ]( bool device )
		{
			s.solver.setDevicePointRoutes( device );
			mfem::Vector out;
			s.solver.assemblePsiIntegrandColumn( state, f, out );
			out.HostRead();
			return out;
		};
		mfem::Vector const host = column( false );
		mfem::Vector const device = column( true );
		s.solver.setDevicePointRoutes( false );
		mfem::Vector lambda;
		s.solver.assembleSourceDerivative( state,
			[ & ]( double, double, double psi )
			{
				double const x = ( psi - boundary )/span;
				if ( !( x > 0.0 ) )
					return 0.0;
				double const value = exact( x );
				if ( !( value > 0.0 ) )
					return 0.0;
				double const g = std::sqrt( value );
				switch ( kind )
				{
					case meq::PsiIntegrandKind::Root: return g;
					case meq::PsiIntegrandKind::SlopeOverRoot: return exact.prime( x )/( 2.0*g*span );
					case meq::PsiIntegrandKind::PsiSlopeOverRoot: return x*exact.prime( x )/( 2.0*g*span );
					case meq::PsiIntegrandKind::RatioOverRoot: return hat( x )/( 2.0*g );
				}
				return 0.0;
			}, lambda );
		lambda.HostRead();
		double worst = 0.0, against = 0.0, scale = 0.0;
		for ( int i = 0; i < host.Size(); ++i )
		{
			worst = std::max( worst, std::abs( double( host( i ) - device( i ) ) ) );
			against = std::max( against, std::abs( double( host( i ) - lambda( i ) ) ) );
			scale = std::max( scale, std::abs( double( host( i ) ) ) );
		}
		std::printf( "  kind %d: device against host %.2e, host against the lambda %.2e, of %.3e\n",
		             static_cast<int>( kind ), worst/scale, against/scale, scale );
		BOOST_TEST( scale > 0.0 );
		BOOST_TEST( worst <= 64.0*eps*scale );
		BOOST_TEST( against <= 64.0*eps*scale, "the column as data is not the column as a lambda" );
	}
	s.solver.setDevicePointRoutes( true );
}
