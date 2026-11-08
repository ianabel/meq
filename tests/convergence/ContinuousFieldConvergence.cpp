#define BOOST_TEST_MODULE ContinuousFieldConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "mfem.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "meq/ContinuousField.hpp"
#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"
#include "ConvergenceHarness.hpp"

/*
 * psi~ AND q~, THE CONTINUOUS PAIR THE FLUX-SURFACE EXTRACTION ROOTS AND
 * WEIGHTS BY (ContinuousField.hpp; MANTA-EQUILIBRIUM-PLAN.md section 4.5;
 * M-364). What the extraction's smoothness and accuracy rest on, asserted
 * on the field itself:
 *
 *   1. the moment interpolant reproduces polynomials of its degree, on
 *      triangles and on quadrilaterals -- the unisolvence, the moments'
 *      orientation and the dof bookkeeping, which a wrong edge orientation or
 *      a miscounted cell basis would each break;
 *   2. psi~ and q~ are continuous across every interior face, to round-off,
 *      and psi~'s trace has psihat's P_k moments on every interior edge;
 *   3. the threaded build is the serial one bit for bit;
 *   4. psi~ converges at psi*'s k+2 and q~ at q*'s k+1, against Solov'ev's
 *      closed form, and is within a small factor of psi*'s own error.
 */
namespace
{
	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	/// Machine epsilon at the build's precision: the round-off bounds below
	/// are multiples of it, so that a single build is asserted at its own
	/// floor and not at double's.
	double const eps = std::numeric_limits<mfem::real_t>::epsilon();

	/// The larger, and NaN if either is: std::max( a, NaN ) is a, which
	/// would let a singular moment system pass as an accurate one.
	double worse( double a, double b )
	{
		return ( std::isnan( a ) || std::isnan( b ) ) ? std::nan( "" ) : std::max( a, b );
	}

	double rate( double coarse, double fine )
	{
		return std::log( coarse/fine )/std::log( 2.0 );
	}

	/// A polynomial of total degree @a p in ( R, z ), every monomial present.
	double polynomial( int p, double radius, double z )
	{
		double value = 0.0;
		for ( int total = 0; total <= p; ++total )
			for ( int a = 0; a <= total; ++a )
				value += ( 1.0 + 0.1*a + 0.01*total )*std::pow( radius, a )*std::pow( z, total - a );
		return value;
	}

	double physicalPolynomial( int p, mfem::ElementTransformation &transform, mfem::IntegrationPoint const &ip )
	{
		mfem::Vector x( 2 );
		transform.Transform( ip, x );
		return polynomial( p, x( 0 ), x( 1 ) );
	}

	/// One Solov'ev solve, kept alive and post-processed. Member order is
	/// load bearing; FluxSurfaceConvergence.cpp's SolvedEquilibrium.
	class SolvedEquilibrium
	{
		public:
			SolvedEquilibrium( Equilibrium const &eqIn, Rectangle const &box, int order, int n )
				: eq( eqIn ),
				  mesh( meq::tests::makeMesh( box, n ) ),
				  sourceCoeff( [ this ]( mfem::Vector const &x ) { return eq.f( x( 0 ), x( 1 ), 0.0 ); } ),
				  psiCoeff( [ this ]( mfem::Vector const &x ) { return eq.psi( x( 0 ), x( 1 ) ); } ),
				  solver( mesh, order )
			{
				solver.setSource( sourceCoeff );
				solver.setBoundaryData( psiCoeff );
				solver.solve();
				solver.postProcess();
			}

			SolvedEquilibrium( SolvedEquilibrium const & ) = delete;
			SolvedEquilibrium &operator=( SolvedEquilibrium const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			mfem::FunctionCoefficient &exactPsi()
			{
				return psiCoeff;
			}

			Equilibrium const &equilibrium() const
			{
				return eq;
			}

		private:
			Equilibrium eq;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			meq::GradShafranovSolver solver;
	};
}

/*
 * 1. THE INTERPOLANT REPRODUCES ITS DEGREE. Fed the traces and cell values of
 * a polynomial of degree p, the moment interpolant of degree p is that
 * polynomial: on triangles ( P_(p-3) cell moments ) and on quadrilaterals
 * ( Q_(p-2) ), at p = 3, 4, 5. Mutated -- a cell monomial in x where it
 * should be in y, or psi~'s interior edges from one side's psi* in place of
 * psihat -- the case and its neighbours fail (M-366).
 */
BOOST_AUTO_TEST_CASE( theInterpolantReproducesItsDegree )
{
	Rectangle const box = meq::tests::standardBox();
	for ( mfem::Element::Type type : { mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL } )
		for ( int p = 3; p <= 5; ++p )
		{
			mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 3, 4, type, false, box.width(), box.height() );
			double const minRadius = box.minRadius;
			double const zMin = box.zMin;
			mesh.Transform( [ minRadius, zMin ]( mfem::Vector const &in, mfem::Vector &out )
			{
				// A shear, so that no edge is aligned with an axis.
				out( 0 ) = in( 0 ) + minRadius + 0.1*in( 1 );
				out( 1 ) = in( 1 ) + zMin;
			} );
			mfem::H1_FECollection collection( p, 2 );
			mfem::FiniteElementSpace space( &mesh, &collection );
			mfem::GridFunction field( &space );
			meq::momentInterpolant( field,
				[ & ]( int e, mfem::IntegrationPoint const &ip )
				{
					thread_local mfem::IsoparametricTransformation transform;
					mesh.GetElementTransformation( e, &transform );
					return physicalPolynomial( p, transform, ip );
				},
				[ & ]( int, mfem::FaceElementTransformations &transform, mfem::IntegrationPoint const &ip )
				{
					return physicalPolynomial( p, transform, ip );
				} );

			double worst = 0.0, scale = 0.0;
			mfem::IsoparametricTransformation transform;
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				mesh.GetElementTransformation( e, &transform );
				mfem::IntegrationRule const &rule = mfem::IntRules.Get( mesh.GetElementBaseGeometry( e ), 2*p );
				for ( int q = 0; q < rule.GetNPoints(); ++q )
				{
					double const exact = physicalPolynomial( p, transform, rule.IntPoint( q ) );
					worst = worse( worst, std::abs( field.GetValue( e, rule.IntPoint( q ) ) - exact ) );
					scale = std::max( scale, std::abs( exact ) );
				}
			}
			// Relative to the polynomial's size: the cell systems are
			// monomial moments, conditioned like a degree-p Vandermonde.
			std::printf( "  %-15s p = %d: worst | interpolant - polynomial | = %.2e of its size\n",
			             type == mfem::Element::TRIANGLE ? "triangles," : "quadrilaterals,", p, worst/scale );
			BOOST_TEST( worst < 65536.0*eps*scale, "the moment interpolant of degree " << p << " does "
			            "not reproduce a polynomial of that degree: " << worst/scale << " of its size" );
		}
}

/*
 * 2. CONTINUOUS, AND CONSISTENT WITH psihat. psi~ and q~ read from the two
 * elements either side of every interior face agree to round-off; psi~'s
 * trace has psihat's P_k moments on every interior edge.
 */
BOOST_AUTO_TEST_CASE( theFieldsAreContinuousAndTakePsihatsMoments )
{
	int const order = 2;
	SolvedEquilibrium solved( Equilibrium::iterExample2(), meq::tests::standardBox(), order, 8 );
	meq::GradShafranovSolver const &solver = solved.theSolver();
	mfem::GridFunction const &psi = solver.continuousPotential();
	mfem::GridFunction const &flux = solver.continuousFlux();
	mfem::GridFunction const &trace = solver.trace();
	mfem::FiniteElementSpace const &traceSpace = *trace.FESpace();
	mfem::Mesh &mesh = *psi.FESpace()->GetMesh();

	double psiScale = 0.0, fluxScale = 0.0;
	for ( int i = 0; i < psi.Size(); ++i )
		psiScale = std::max( psiScale, static_cast<double>( std::abs( psi( i ) ) ) );
	for ( int i = 0; i < flux.Size(); ++i )
		fluxScale = std::max( fluxScale, static_cast<double>( std::abs( flux( i ) ) ) );

	double psiJump = 0.0, fluxJump = 0.0, moment = 0.0;
	mfem::FaceElementTransformations faceTransform;
	mfem::IsoparametricTransformation one, two;
	mfem::Vector a( 2 ), b( 2 ), local, shape;
	mfem::Array<int> dofs;
	for ( int f = 0; f < mesh.GetNumFaces(); ++f )
	{
		mesh.GetFaceElementTransformations( f, faceTransform, one, two );
		if ( faceTransform.Elem2No < 0 )
			continue;
		mfem::FiniteElement const *traceFe = traceSpace.GetFaceElement( f );
		traceSpace.GetFaceVDofs( f, dofs );
		trace.GetSubVector( dofs, local );
		shape.SetSize( traceFe->GetDof() );
		std::vector<double> moments( static_cast<std::size_t>( order + 1 ), 0.0 );
		mfem::IntegrationRule const &rule = mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*order + 6 );
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			faceTransform.SetAllIntPoints( &ip );
			double const left = psi.GetValue( faceTransform.Elem1No, faceTransform.GetElement1IntPoint() );
			double const right = psi.GetValue( faceTransform.Elem2No, faceTransform.GetElement2IntPoint() );
			psiJump = worse( psiJump, std::abs( left - right ) );
			flux.GetVectorValue( faceTransform.Elem1No, faceTransform.GetElement1IntPoint(), a );
			flux.GetVectorValue( faceTransform.Elem2No, faceTransform.GetElement2IntPoint(), b );
			fluxJump = worse( fluxJump, std::hypot( a( 0 ) - b( 0 ), a( 1 ) - b( 1 ) ) );
			traceFe->CalcShape( ip, shape );
			double const hat = shape*local;
			double const x = 2.0*ip.x - 1.0;
			double const legendre[ 3 ] = { 1.0, x, 1.5*x*x - 0.5 };
			for ( int m = 0; m <= order; ++m )
				moments[ static_cast<std::size_t>( m ) ] += ip.weight*( left - hat )*legendre[ m ];
		}
		for ( double value : moments )
			moment = worse( moment, std::abs( value ) );
	}
	std::printf( "\n  across interior faces: | psi~ jump | %.2e, | q~ jump | %.2e of their scales; "
	             "psi~'s trace less psihat, worst P_k moment, %.2e of psi~'s scale\n",
	             psiJump/psiScale, fluxJump/fluxScale, moment/psiScale );
	BOOST_TEST( psiJump < 512.0*eps*psiScale, "psi~ is not continuous across a face: " << psiJump );
	BOOST_TEST( fluxJump < 512.0*eps*fluxScale, "q~ is not continuous across a face: " << fluxJump );
	BOOST_TEST( moment < 512.0*eps*psiScale, "psi~'s trace does not have psihat's P_k moments: " << moment );
}

/*
 * 3. THE THREADED BUILD IS THE SERIAL ONE, BIT FOR BIT: every dof is written
 * by one entity, and the vertex means are summed in element order.
 */
BOOST_AUTO_TEST_CASE( theThreadedBuildIsTheSerialOne )
{
#ifdef _OPENMP
	SolvedEquilibrium solved( Equilibrium::iterExample2(), meq::tests::standardBox(), 3, 8 );
	meq::GradShafranovSolver const &solver = solved.theSolver();
	int const threads = omp_get_max_threads();
	omp_set_num_threads( 1 );
	meq::ContinuousFields const serial = meq::buildContinuousFields(
		solver.postProcessedPotential(), solver.postProcessedFlux(), solver.trace() );
	omp_set_num_threads( std::max( threads, 4 ) );
	meq::ContinuousFields const threaded = meq::buildContinuousFields(
		solver.postProcessedPotential(), solver.postProcessedFlux(), solver.trace() );
	omp_set_num_threads( threads );
	std::size_t differ = 0;
	for ( int i = 0; i < serial.potential->Size(); ++i )
		differ += ( *serial.potential )( i ) != ( *threaded.potential )( i );
	for ( int i = 0; i < serial.flux->Size(); ++i )
		differ += ( *serial.flux )( i ) != ( *threaded.flux )( i );
	std::printf( "\n  one thread against %d: %zu of %d dofs differ\n", std::max( threads, 4 ), differ,
	             serial.potential->Size() + serial.flux->Size() );
	BOOST_TEST( differ == 0u, "the threaded build differs from the serial one at " << differ << " dofs" );
#else
	std::printf( "\n  built without OpenMP: one thread is all there is, and nothing to compare\n" );
#endif
}

/*
 * 4. psi~ KEEPS psi*'S ORDER AND q~ KEEPS q*'S. Against Solov'ev's closed
 * form, k = 1, 2, 3 on n = 8, 16: psi~ at k+2 and q~ at k+1, the same
 * two-tier slack as the rest of the convergence suite; and psi~'s error within
 * 3x of psi*'s on every mesh -- the interpolant moves the field by the face
 * jumps, which are psi*'s own error.
 */
BOOST_AUTO_TEST_CASE( theContinuousFieldsKeepTheirSourcesOrder )
{
	// THE RATES ARE ASSERTED WHERE psi* ITSELF STILL CONVERGES. In single,
	// psi* reaches the round-off floor at k = 2 on n = 16 ( 1.7e-7 ) and its
	// error GROWS from n = 8 to 16 at k = 3, so psi~ cannot be asked for k+2
	// there; it is still asked to stay within 3x of psi* at every order.
	auto rated = []( int order )
	{
#ifdef MFEM_USE_SINGLE
		return order == 1;
#else
		return order >= 1;
#endif
	};
	Rectangle const box = meq::tests::standardBox();
	Equilibrium const eq = Equilibrium::iterExample2();
	mfem::VectorFunctionCoefficient exactFlux( 2, [ &eq ]( mfem::Vector const &x, mfem::Vector &q )
	{
		double qR = 0.0, qZ = 0.0;
		eq.flux( x( 0 ), x( 1 ), qR, qZ );
		q( 0 ) = qR;
		q( 1 ) = qZ;
	} );
	std::printf( "\n  %3s %4s %12s %12s %8s %12s %12s %8s\n", "k", "n", "psi* L2", "psi~ L2", "rate",
	             "q* L2", "q~ L2", "rate" );
	for ( int order : { 1, 2, 3 } )
	{
		double previousPsi = 0.0, previousFlux = 0.0;
		for ( int n : { 8, 16 } )
		{
			SolvedEquilibrium solved( eq, box, order, n );
			meq::GradShafranovSolver const &solver = solved.theSolver();
			double const starPsi = solver.postProcessedPotential().ComputeL2Error( solved.exactPsi() );
			double const tildePsi = solver.continuousPotential().ComputeL2Error( solved.exactPsi() );
			double const starFlux = solver.postProcessedFlux().ComputeL2Error( exactFlux );
			double const tildeFlux = solver.continuousFlux().ComputeL2Error( exactFlux );
			if ( previousPsi > 0.0 && rated( order ) )
			{
				double const psiRate = rate( previousPsi, tildePsi );
				double const fluxRate = rate( previousFlux, tildeFlux );
				std::printf( "  %3d %4d %12.3e %12.3e %8.3f %12.3e %12.3e %8.3f\n", order, n, starPsi,
				             tildePsi, psiRate, starFlux, tildeFlux, fluxRate );
				BOOST_TEST( psiRate > order + 2 - 0.3, "psi~ converges at " << psiRate
				            << " at k = " << order << ", not psi*'s k+2" );
				BOOST_TEST( fluxRate > order + 1 - 0.3, "q~ converges at " << fluxRate
				            << " at k = " << order << ", not q*'s k+1" );
			}
			else
				std::printf( "  %3d %4d %12.3e %12.3e %8s %12.3e %12.3e %8s\n", order, n, starPsi,
				             tildePsi, previousPsi > 0.0 ? "(floor)" : "-", starFlux, tildeFlux,
				             previousPsi > 0.0 ? "(floor)" : "-" );
			BOOST_TEST( tildePsi < 3.0*starPsi, "psi~'s error " << tildePsi << " is more than 3x psi*'s "
			            << starPsi << " at k = " << order << ", n = " << n );
			previousPsi = tildePsi;
			previousFlux = tildeFlux;
		}
	}
}
