/*
 * DEVICE AUDIT TIMINGS -- not a ctest; every number is a timing.
 *
 * Times, under a configured mfem::Device, each of the accuracy pieces that
 * landed between c4ee278 and 1e8a3ad, called the way the q-driven loop and
 * the write phase call them, on a normalised MHD source of the q-driven
 * shape (constant p', gg' linear in Psi) on the standard box:
 *
 *   continuous fields   meq::buildContinuousFields( psi*, q*, psihat )
 *   labelledAxis        psi~ at the located axis
 *   extraction          extractFluxSurfaces on psi~, q~ (24 x 128)
 *   integrateOverDomain the delivered Phi_edge's quadrature
 *   delivered I_p       deliveredPlasmaCurrent()
 *   table reach         sourceTableReach()
 *   residual (host)     evaluateResidual(): a plain DarcyNPCOperator::Mult,
 *                       the operator linearise() and applyRetainedJacobian()
 *                       build
 *   addSourceColumns    one pass of the source loop
 *   linearise           the retained factorisation, with psi_ax's border
 *   back-solve          solveLinearised( passes = 0 )
 *   applyRetainedJacobian
 *   assembleSourceDerivative
 *   tangent passes=p    linearisedTangent( ..., p )
 *   inverse norm        linearisedInverseNorm( 20 )
 *
 * Usage: DeviceAudit [--device cpu|debug|cuda] [--sizes 16,32,64] [--order 2]
 *                    [--repeats 3]
 *
 * Under a Device it takes the driver's device defaults (batched assembly and
 * local factorisation, cuDSS under cuda); AUDIT_SERIAL=1 keeps the host's
 * keys and AUDIT_PARDISO=1 keeps PARDISO, the two controls M-383 ran.
 * MEASUREMENTS.md M-384 is what it measured.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "mfem.hpp"
#include "meq/ContinuousField.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/FluxExtraction.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/Source.hpp"

#if __has_include( <nvtx3/nvToolsExt.h> ) && defined( MFEM_USE_CUDA )
#include <nvtx3/nvToolsExt.h>
#define AUDIT_NVTX 1
#endif

namespace meq::tests
{
	struct Rectangle
	{
		double minRadius, maxRadius, zMin, zMax;
		double width() const { return maxRadius - minRadius; }
		double height() const { return zMax - zMin; }
	};
	inline Rectangle standardBox() { return Rectangle{ 0.6, 1.4, -0.6, 0.6 }; }
	inline mfem::Mesh makeMesh( Rectangle const &box, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE,
		                                               false, box.width(), box.height() );
		double const minRadius = box.minRadius, zMin = box.zMin;
		mesh.Transform( [ minRadius, zMin ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + minRadius;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}
}

namespace
{
	using Solver = meq::GradShafranovSolver;

	double now()
	{
		return std::chrono::duration<double>(
			std::chrono::steady_clock::now().time_since_epoch() ).count();
	}

	struct Range
	{
		explicit Range( char const *name )
		{
#ifdef AUDIT_NVTX
			nvtxRangePushA( name );
#else
			static_cast<void>( name );
#endif
		}
		~Range()
		{
#ifdef AUDIT_NVTX
			nvtxRangePop();
#endif
		}
	};

	std::vector<int> parseList( char const *text )
	{
		std::vector<int> out;
		std::stringstream in( text );
		std::string item;
		while ( std::getline( in, item, ',' ) )
			out.push_back( std::atoi( item.c_str() ) );
		return out;
	}

	/// Best and median of @a repeats calls, in seconds.
	struct Timing
	{
		double best = 0.0, median = 0.0;
	};

	Timing timeIt( char const *name, int repeats, std::function<void()> const &f )
	{
		std::vector<double> t;
		for ( int r = 0; r < repeats; ++r )
		{
			Range const range( name );
			double const t0 = now();
			f();
			t.push_back( now() - t0 );
		}
		std::sort( t.begin(), t.end() );
		return { t.front(), t[ t.size()/2 ] };
	}

	void row( char const *name, Timing const &t, double reference )
	{
		std::printf( "  %-28s %10.3f ms  (median %10.3f ms)  %7.2f x residual\n",
		             name, 1e3*t.best, 1e3*t.median,
		             reference > 0.0 ? t.best/reference : 0.0 );
		std::fflush( stdout );
	}
}

int main( int argc, char **argv )
{
	std::string device = "cpu";
	std::vector<int> sizes { 16, 32, 64 };
	int order = 2;
	int repeats = 3;
	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--device" && i + 1 < argc )
			device = argv[ ++i ];
		else if ( arg == "--sizes" && i + 1 < argc )
			sizes = parseList( argv[ ++i ] );
		else if ( arg == "--order" && i + 1 < argc )
			order = std::atoi( argv[ ++i ] );
		else if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
	}

	std::unique_ptr<mfem::Device> deviceHandle;
	if ( device != "cpu" )
	{
		deviceHandle = std::make_unique<mfem::Device>( device.c_str() );
		deviceHandle->Print();
	}
	bool const onDevice = device != "cpu";

	double const g0 = 2.20, g1 = 0.55, mu0 = 1.0;
	std::vector<double> const g2 { g0*g0, 2.0*g0*g1, g1*g1 };

	for ( int n : sizes )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );

		std::vector<meq::Knot> knots;
		for ( double x : { 0.0, 1.0 } )
			knots.push_back( meq::Knot{ x, g1*( g0 + g1*x ), g1*g1 } );
		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		auto ggPrime = std::make_shared<meq::SplineProfile const>( knots );
		meq::NormalisedMHDSource source( pPrime, ggPrime, 0.30, mu0 );

		Solver solver( mesh, order );
		if ( onDevice && !std::getenv( "AUDIT_SERIAL" ) )
		{
			solver.setAssemblyMode( Solver::AssemblyMode::Batched );
			solver.setLocalFactorMode( Solver::LocalFactorMode::Batched );
			if ( device == "cuda" && !std::getenv( "AUDIT_PARDISO" )
			     && Solver::traceSolverAvailable( Solver::TraceSolver::cuDSS ) )
				solver.setTraceSolver( Solver::TraceSolver::cuDSS );
		}
		mfem::FunctionCoefficient guess(
			[ box ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			} );
		mfem::ConstantCoefficient zero( 0.0 );
		solver.setSource( source, 0.30 );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 60 );

		std::printf( "\n=== n = %d (%d triangles), k = %d, device %s ===\n",
		             n, mesh.GetNE(), order, device.c_str() );

		double const t0 = now();
		{
			Range const range( "solve" );
			solver.solve();
		}
		double const solveTime = now() - t0;
		Solver::StepProfile const &profile = solver.stepProfile();
		double const residualPerCall = profile.residualCalls > 0
		                               ? profile.residualSeconds/profile.residualCalls : 0.0;
		std::printf( "  solve %.3f s, %d Newton steps, psi_ax %.7e; device integrand %s;"
		             " Newton residual %.3f ms/call (%ld calls), gradient %.3f ms/call,"
		             " trace solve %.3f ms/call\n",
		             solveTime, solver.newtonIterations(), solver.psiAxis(),
		             solver.deviceIntegrandTaken() ? "yes" : "NO",
		             1e3*residualPerCall, profile.residualCalls,
		             profile.gradientCalls ? 1e3*profile.gradientSeconds/profile.gradientCalls : 0.0,
		             profile.traceSolveCalls ? 1e3*profile.traceSolveSeconds/profile.traceSolveCalls : 0.0 );

		Timing const post = timeIt( "postProcess", 1, [ & ] { solver.postProcess(); } );
		row( "postProcess (Reconstruct)", post, residualPerCall );

		// ---- the continuous extraction ----
		Timing const continuous = timeIt( "continuousFields", repeats, [ & ]
		{
			meq::ContinuousFields fields = meq::buildContinuousFields(
				solver.postProcessedPotential(), solver.postProcessedFlux(), solver.trace() );
			static_cast<void>( fields.potential->HostRead() );
		} );
		row( "buildContinuousFields", continuous, residualPerCall );

		// The interpolant's two halves (M-385): its quadratures, once per set
		// of spaces, and its application, per extraction, on either route.
		std::unique_ptr<meq::ContinuousInterpolant> interpolant;
		Timing const quadratures = timeIt( "interpolantBuild", repeats, [ & ]
		{
			interpolant = std::make_unique<meq::ContinuousInterpolant>( *solver.postProcessedPotential().FESpace(),
				*solver.postProcessedFlux().FESpace(), *solver.trace().FESpace() );
		} );
		row( "  interpolant quadratures", quadratures, residualPerCall );
		Timing const applyHost = timeIt( "interpolantHost", repeats, [ & ]
		{
			static_cast<void>( interpolant->apply( solver.postProcessedPotential(), solver.postProcessedFlux(),
			                                       solver.trace(), meq::ContinuousInterpolant::Route::Host ) );
		} );
		row( "  interpolant apply, host", applyHost, residualPerCall );
		if ( onDevice )
		{
			Timing const applyDevice = timeIt( "interpolantDevice", repeats, [ & ]
			{
				static_cast<void>( interpolant->apply( solver.postProcessedPotential(), solver.postProcessedFlux(),
				                                       solver.trace(), meq::ContinuousInterpolant::Route::Device ) );
			} );
			row( "  interpolant apply, device", applyDevice, residualPerCall );
			std::printf( "      the device route ran: %s\n", interpolant->lastRanOnDevice() ? "yes" : "NO" );
		}

		meq::ContourTracer const tracer( solver );   // builds and caches psi~, q~ once
		meq::CriticalPointFinder const finder( solver );
		meq::CriticalPoint located;
		Timing const findAxis = timeIt( "findAxis", repeats, [ & ] { located = finder.findAxis(); } );
		row( "findAxis (pre-existing)", findAxis, residualPerCall );
		meq::CriticalPoint labelled;
		Timing const label = timeIt( "labelledAxis", repeats, [ & ]
		{
			labelled = meq::labelledAxis( solver, tracer, located );
		} );
		row( "labelledAxis", label, residualPerCall );

		meq::FluxFamilyOptions options;
		std::size_t surfaces = 0;
		Timing const extract = timeIt( "extractFluxSurfaces", repeats, [ & ]
		{
			meq::FluxSurfaceFamily const family =
				meq::extractFluxSurfaces( tracer, labelled, solver.psiBoundary(), options );
			surfaces = family.surfaces.size();
		} );
		row( "extractFluxSurfaces 24x128", extract, residualPerCall );

		// ---- the delivered scale pins ----
		double const axis = solver.psiAxis(), boundary = solver.psiBoundary();
		double const span = axis - boundary;
		double flux = 0.0;
		Timing const integrate = timeIt( "integrateOverDomain", repeats, [ & ]
		{
			flux = solver.integrateOverDomain(
				[ & ]( double radius, double, double psi )
				{
					if ( !( radius > 0.0 ) || span == 0.0 )
						return 0.0;
					double const x = ( psi - boundary )/span;
					if ( !( x > 0.0 ) )
						return 0.0;
					double const value = g2[ 0 ] + x*( g2[ 1 ] + x*g2[ 2 ] );
					return value > 0.0 ? std::sqrt( value )/radius : 0.0;
				} );
		} );
		row( "integrateOverDomain (Phi)", integrate, residualPerCall );

		long const takenBefore = solver.iterateRoutesTakenOnDevice();
		long const declinedBefore = solver.iterateRoutesDeclinedOnDevice();
		double current = 0.0;
		Timing const delivered = timeIt( "deliveredPlasmaCurrent", repeats, [ & ]
		{
			current = solver.deliveredPlasmaCurrent();
		} );
		row( "deliveredPlasmaCurrent", delivered, residualPerCall );
		std::printf( "      device route: %ld taken, %ld declined over %d calls\n",
		             solver.iterateRoutesTakenOnDevice() - takenBefore,
		             solver.iterateRoutesDeclinedOnDevice() - declinedBefore, repeats );

		Timing const reach = timeIt( "sourceTableReach", repeats, [ & ]
		{
			static_cast<void>( solver.sourceTableReach() );
		} );
		row( "sourceTableReach", reach, residualPerCall );

		// ---- the linearised solve ----
		mfem::Vector const state( solver.state() );
		mfem::Vector residual;
		Timing const hostResidual = timeIt( "evaluateResidual", repeats, [ & ]
		{
			solver.evaluateResidual( state, residual );
			static_cast<void>( residual.HostRead() );
		} );
		row( "evaluateResidual (plain NPC)", hostResidual, residualPerCall );

		mfem::Vector column( state.Size() );
		Timing const sourceLoop = timeIt( "addSourceColumns", repeats, [ & ]
		{
			column.HostWrite();
			for ( int i = 0; i < column.Size(); ++i )
				column( i ) = 0.0;
			std::vector<mfem::Vector *> out { &column };
			solver.addSourceColumns( state, 1,
				[ &source ]( double radius, double z, double psi, double *g )
				{
					g[ 0 ] = source.f( radius, z, psi );
					return true;
				}, out );
		} );
		row( "addSourceColumns (one pass)", sourceLoop, residualPerCall );

		Timing const lin = timeIt( "linearise", repeats, [ & ] { solver.linearise( state ); } );
		row( "linearise", lin, residualPerCall );

		mfem::Vector dRdp, dx, applied;
		std::vector<double> ds;
		Timing const derivative = timeIt( "assembleSourceDerivative", repeats, [ & ]
		{
			solver.assembleSourceDerivative( state,
				[ &source ]( double radius, double z, double psi )
				{
					return source.dFdPsi( radius, z, psi )*0.01 + 1.0/radius;
				}, dRdp );
		} );
		row( "assembleSourceDerivative", derivative, residualPerCall );

		Timing const backSolve = timeIt( "solveLinearised0", repeats, [ & ]
		{
			solver.solveLinearised( dRdp, { 0.0 }, dx, ds, 0 );
			static_cast<void>( dx.HostRead() );
		} );
		row( "back-solve (passes = 0)", backSolve, residualPerCall );

		Timing const apply = timeIt( "applyRetainedJacobian", repeats, [ & ]
		{
			solver.applyRetainedJacobian( dx, applied );
			static_cast<void>( applied.HostRead() );
		} );
		row( "applyRetainedJacobian", apply, residualPerCall );

		for ( int passes : { 1, 2, 4 } )
		{
			std::string const name = "tangent passes=" + std::to_string( passes );
			Timing const tangent = timeIt( name.c_str(), repeats, [ & ]
			{
				solver.linearisedTangent( dRdp, {}, nullptr, dx, ds, passes );
				static_cast<void>( dx.HostRead() );
			} );
			row( name.c_str(), tangent, residualPerCall );
			dx.HostRead();
			std::printf( "      dpsi_ax %.9e, max|dx| %.9e\n", ds.empty() ? 0.0 : ds[ 0 ], dx.Normlinf() );
		}

		Timing const inverse = timeIt( "linearisedInverseNorm", 1, [ & ]
		{
			std::printf( "      1/sigma_min estimate %.9e\n", solver.linearisedInverseNorm( 20 ) );
		} );
		row( "linearisedInverseNorm(20)", inverse, residualPerCall );

		std::printf( "  checks: %zu surfaces, Phi %.9e, mu0 I_p %.9e, psi_ax %.9e labelled %.9e\n",
		             surfaces, flux, current, axis, labelled.psi );
	}
	return 0;
}
