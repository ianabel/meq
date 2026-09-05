/*
 * NPC thread scaling: what AssemblyMode::Threaded is worth on a NONLINEAR
 * solve, and whether PARDISO's MKL threads are reachable at the same time.
 *
 * NOT A CTEST, per the standing rule for tests/performance/: every number
 * below is a timing, and a threaded timing on this machine is a measurement
 * about the machine. It exits non-zero only for the CORRECTNESS properties that
 * make the timings mean anything -- threaded assembly reproducing serial
 * assembly bit for bit on a nonlinear source, and the trace solvers agreeing.
 *
 * WHY THIS EXISTS BESIDE TraceSolverScaling, WHICH ALREADY SWEEPS BOTH THREAD
 * AXES. That harness times a LINEAR solve: assembly, reduction and one trace
 * factorisation. It was written when AssemblyMode::Threaded meant ComputeH()'s
 * element loop and nothing else, so a linear problem exercised the whole of it.
 *
 * It does not any more. MFEM now threads MultNL() as well -- the residual and
 * the Jacobian assembly, and therefore NPCResidual() and NPCGradient(), which
 * is every NPC step -- and measures a whole NPC step at 1.9-2.1x on eight
 * threads. None of that appears on a linear path, because a linear problem
 * never calls MultNL at all. So the option's value has moved to a loop the
 * existing harness structurally cannot reach, and MEQ's recorded reason for
 * defaulting to Serial was measured against the version that lacked it.
 *
 * THE THREE QUESTIONS, in the order the answers matter:
 *
 *   1. Is Threaded now worth taking by default? CLAUDE.md records the
 *      measurement that said no -- HighBetaConvergence went 21.5 s to 39 s
 *      under an automatic gate, because MFEM forks a team and buffers PER CALL
 *      and a caller that assembles hundreds of times inside a bordered Newton
 *      pays it every time. That argument is about ComputeH(). MultNL() is
 *      called far more often than assembly is, so the same caller now has
 *      something to gain from the same flag. --highbeta is that case.
 *
 *   2. Does threaded assembly defuse MKL_NUM_THREADS? Measured in isolation it
 *      does, completely: MKL suppresses its own threading inside an active
 *      OpenMP parallel region, so the element-local dgetrs and dgemm cost the
 *      same at MKL=8 as at MKL=1 -- while the SAME kernels in a serial element
 *      loop cost 5.9x more at k=3. If that holds in situ, then
 *      MKL_NUM_THREADS=1 stops being a requirement and becomes a requirement
 *      of the SERIAL mode only.
 *
 *   3. Which makes PARDISO's threads reachable. CLAUDE.md's *What to do* item 0
 *      says MEQ cannot spend them because MKL_NUM_THREADS is process-wide and
 *      the setting PARDISO wants is the setting that ruins ComputeH(). If (2)
 *      holds, the trace solve runs on the master thread OUTSIDE any parallel
 *      region and gets all of them, while the element loop is nested and gets
 *      none -- and no mkl_set_num_threads_local() plumbing is needed at all.
 *
 * ONE THREAD-COUNT WARNING IS BUILT IN, because it is a configuration a user
 * can reach by accident and it is catastrophic rather than merely slow. Threaded
 * assembly with OMP_NUM_THREADS=1 and MKL_NUM_THREADS>1 measured 12.4 s against
 * 0.069 s serial on the isolated kernels -- a team of one thread with MKL
 * threading live. CLAUDE.md records Threaded-at-one-thread as 0.86x; that is
 * true only at MKL=1.
 *
 * Usage:  NpcThreadScaling [--orders k,k] [--sizes n,n] [--repeats N]
 *                          [--highbeta] [--pedestal]
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#ifdef MFEM_USE_OPENMP
	#include <omp.h>
#endif

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/ManufacturedNonlinear.hpp"
#include "analytic/PressurePedestal.hpp"

namespace
{
	using Solver = meq::GradShafranovSolver;
	using AM = Solver::AssemblyMode;
	using TS = Solver::TraceSolver;

	double now()
	{
		using namespace std::chrono;
		return duration<double>( steady_clock::now().time_since_epoch() ).count();
	}

	struct Box { double rMin, rMax, zMin, zMax; };

	Box box() { return { 0.6, 1.4, -0.6, 0.6 }; }

	mfem::Mesh makeMesh( Box const &b, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false,
			b.rMax - b.rMin, b.zMax - b.zMin );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			double *c = mesh.GetVertex( v );
			c[ 0 ] += b.rMin;
			c[ 1 ] += b.zMin;
		}
		return mesh;
	}

	/// A meq::Source over any fixture with f() and dFdPsi(). The same adapter
	/// tests/convergence/ConvergenceHarness.hpp carries; duplicated rather than
	/// shared because that header pulls in Boost.Test, which this binary is
	/// deliberately not linked against.
	template<typename Equilibrium>
	class EquilibriumSource : public meq::Source
	{
		public:
			explicit EquilibriumSource( Equilibrium const &eqIn ) : eq( eqIn ) {}

			double f( double r, double z, double psi ) const override
			{
				return eq.f( r, z, psi );
			}

			double dFdPsi( double r, double z, double psi ) const override
			{
				return eq.dFdPsi( r, z, psi );
			}

		private:
			Equilibrium const &eq;
	};

	/// What one configuration cost, and what it produced.
	struct Run
	{
		double prepareTime = 0.0;
		double solveTime = 0.0;
		int newtonIterations = 0;
		bool converged = false;
		std::vector<double> psi;
		std::vector<double> flux;
	};

	std::vector<double> copyOf( mfem::GridFunction const &g )
	{
		std::vector<double> out( static_cast<size_t>( g.Size() ) );
		for ( int i = 0; i < g.Size(); ++i )
			out[ static_cast<size_t>( i ) ] = g( i );
		return out;
	}

	double worstDifference( std::vector<double> const &a,
	                        std::vector<double> const &b )
	{
		if ( a.size() != b.size() )
			return std::numeric_limits<double>::infinity();
		double worst = 0.0;
		for ( size_t i = 0; i < a.size(); ++i )
			worst = std::max( worst, std::fabs( a[ i ] - b[ i ] ) );
		return worst;
	}

	/// One solve, timed. `prepare()` is timed separately from the rest of
	/// `solve()` because they are different loops: prepare() is ComputeH()'s and
	/// the remainder is MultNL()'s, and the whole point of this harness is that
	/// the second one is new.
	/// `datum` is the fixture's own boundary condition, and it is NOT a free
	/// choice in either case.
	///
	/// For a fixture with an exact solution it must be that solution's trace, or
	/// the solve is of a different problem -- a valid one, but not the one whose
	/// Newton history is being timed. For the GS-2 pedestal there is no exact
	/// solution, and the datum is section 4.2's sign-changing ramp: it has to
	/// change sign so that the pedestal layer lies INSIDE the mesh rather than
	/// pressed against the boundary, and it has to be non-zero because that
	/// source vanishes at psi = 0 -- with homogeneous data psi == 0 SOLVES the
	/// problem and Newton stops on it in no iterations at all.
	/// tests/convergence/PedestalConvergence.cpp gives the full account.
	template<typename Fixture>
	Run runOnce( Box const &b, int order, int n, Fixture const &fixture,
	             mfem::Coefficient &datum, bool guessFromDatum, AM mode, TS trace )
	{
		Run out;
		mfem::Mesh mesh = makeMesh( b, n );
		Solver solver( mesh, order );
		solver.setAssemblyMode( mode );
		solver.setTraceSolver( trace );

		EquilibriumSource<Fixture> const source( fixture );
		solver.setSource( source );
		solver.setBoundaryData( datum );

		// Only the trivial-branch cases need a guess, and giving one where it is
		// not needed would change the iteration count this harness is comparing.
		if ( guessFromDatum )
			solver.setInitialGuess( datum );

		double const t0 = now();
		solver.prepare();
		double const t1 = now();
		out.prepareTime = t1 - t0;

		try
		{
			solver.solve();
			out.converged = true;
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}
		out.solveTime = now() - t1;
		out.newtonIterations = solver.newtonIterations();

		if ( out.converged )
		{
			out.psi = copyOf( solver.potential() );
			out.flux = copyOf( solver.flux() );
		}
		return out;
	}

	std::vector<int> parseList( char const *text )
	{
		std::vector<int> out;
		std::string item;
		for ( char const *c = text; ; ++c )
		{
			if ( *c == ',' || *c == '\0' )
			{
				if ( !item.empty() )
					out.push_back( std::atoi( item.c_str() ) );
				item.clear();
				if ( *c == '\0' )
					break;
			}
			else
				item.push_back( *c );
		}
		return out;
	}

	char const *traceName( TS t )
	{
		switch ( t )
		{
			case TS::UMFPack: return "UMFPack";
			case TS::Pardiso: return "PARDISO";
			case TS::cuDSS:   return "cuDSS";
		}
		return "?";
	}

	int envInt( char const *name, int fallback )
	{
		char const *v = std::getenv( name );
		return v ? std::atoi( v ) : fallback;
	}
}

int main( int argc, char **argv )
{
	std::vector<int> orders = { 2, 3 };
	std::vector<int> sizes = { 32, 48 };
	int repeats = 3;
	bool wantPedestal = true;
	bool wantExample5 = true;
	// Which trace solvers to include. Restricting matters here in a way it does
	// not in TraceSolverScaling: UMFPACK and PARDISO respond to MKL_NUM_THREADS
	// in OPPOSITE directions, so a row that sweeps both at MKL > 1 is dominated
	// by UMFPACK's collapse and says nothing about PARDISO.
	std::string wantTrace = "all";

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--orders" && i + 1 < argc )
			orders = parseList( argv[ ++i ] );
		else if ( arg == "--sizes" && i + 1 < argc )
			sizes = parseList( argv[ ++i ] );
		else if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
		else if ( arg == "--pedestal" )
			wantExample5 = false;
		else if ( arg == "--example5" )
			wantPedestal = false;
		else if ( arg == "--trace" && i + 1 < argc )
			wantTrace = argv[ ++i ];
	}

	std::printf( "\n=== MEQ NPC thread scaling ===\n" );

	int const ompThreads =
#ifdef MFEM_USE_OPENMP
		omp_get_max_threads();
#else
		1;
#endif
	int const mklRequested = envInt( "MKL_NUM_THREADS", 0 );

	std::printf( "  OMP_NUM_THREADS=%-5s MKL_NUM_THREADS=%-5s "
	             "omp_get_max_threads()=%d, best of %d\n",
	             std::getenv( "OMP_NUM_THREADS" ) ? std::getenv( "OMP_NUM_THREADS" ) : "unset",
	             std::getenv( "MKL_NUM_THREADS" ) ? std::getenv( "MKL_NUM_THREADS" ) : "unset",
	             ompThreads, repeats );
	std::printf( "  build:" );
#ifdef MFEM_USE_OPENMP
	std::printf( " OPENMP" );
#endif
#ifdef MFEM_THREAD_SAFE
	std::printf( " THREAD_SAFE" );
#endif
#ifdef MFEM_USE_LAPACK
	std::printf( " LAPACK" );
#endif
#ifdef MFEM_USE_MKL_PARDISO
	std::printf( " MKL_PARDISO" );
#endif
#ifdef MFEM_USE_SUITESPARSE
	std::printf( " SUITESPARSE" );
#endif
	std::printf( "\n" );

	bool const canThread = Solver::assemblyModeAvailable( AM::Threaded );
	if ( !canThread )
		std::printf( "\n  *** this MFEM cannot thread the element loop "
		             "(needs MFEM_USE_OPENMP and MFEM_THREAD_SAFE); the "
		             "threaded columns are omitted\n" );

	// The one configuration warning worth building in. See the file comment.
	if ( canThread && ompThreads == 1 && mklRequested > 1 )
		std::printf( "\n  *** WARNING: OMP_NUM_THREADS=1 with MKL_NUM_THREADS=%d.\n"
		             "      Measured on the isolated element-local kernels, a\n"
		             "      one-thread OpenMP team with MKL threading live cost\n"
		             "      12.4 s against 0.069 s serial at k=3 -- 178x. Either\n"
		             "      raise OMP_NUM_THREADS or set MKL_NUM_THREADS=1.\n",
		             mklRequested );

	std::vector<TS> traceSolvers;
	for ( TS t : { TS::UMFPack, TS::Pardiso } )
	{
		if ( !Solver::traceSolverAvailable( t ) )
			continue;
		if ( wantTrace == "umfpack" && t != TS::UMFPack )
			continue;
		if ( wantTrace == "pardiso" && t != TS::Pardiso )
			continue;
		traceSolvers.push_back( t );
	}
	if ( traceSolvers.empty() )
	{
		std::printf( "\n  no trace solver matches --trace %s in this build\n",
		             wantTrace.c_str() );
		return 1;
	}

	int failures = 0;
	double worstThreadedPsi = 0.0;
	double worstThreadedFlux = 0.0;
	int newtonMismatches = 0;
	double worstTraceAgreement = 0.0;

	meq::analytic::ManufacturedNonlinear const example5
		= meq::analytic::ManufacturedNonlinear::example5();
	meq::analytic::PressurePedestal const pedestal
		= meq::analytic::PressurePedestal::pedestal();

	std::printf( "\n  a whole nonlinear solve, serial assembly against threaded\n" );
	std::printf( "    %-10s %2s %5s %8s %9s %9s %8s %9s %9s %8s %5s\n",
	             "case", "k", "n", "solver",
	             "prep_ser", "prep_thr", "prep_x",
	             "solve_ser", "solve_thr", "solve_x", "its" );
	// Flushed, because a row at MKL > 1 with UMFPACK can take minutes and an
	// unflushed header makes that look like a hang rather than a measurement.
	std::fflush( stdout );

	for ( int pass = 0; pass < 2; ++pass )
	{
		char const *name = ( pass == 0 ) ? "example5" : "pedestal";
		if ( pass == 0 && !wantExample5 ) continue;
		if ( pass == 1 && !wantPedestal ) continue;

		// example5's datum is the trace of its own exact solution; the
		// pedestal's is section 4.2's sign-changing ramp, which is also its
		// initial guess. See runOnce().
		double const zMax = box().zMax;
		mfem::FunctionCoefficient exactDatum( [ &example5 ]( mfem::Vector const &x )
		{
			return example5.psi( x( 0 ), x( 1 ) );
		} );
		mfem::FunctionCoefficient rampDatum( [ zMax ]( mfem::Vector const &x )
		{
			return 0.3*x( 1 )/zMax;
		} );
		mfem::Coefficient &datum = ( pass == 0 )
			? static_cast<mfem::Coefficient &>( exactDatum )
			: static_cast<mfem::Coefficient &>( rampDatum );
		bool const guessFromDatum = ( pass != 0 );

		for ( int order : orders )
		{
			for ( int n : sizes )
			{
				for ( TS trace : traceSolvers )
				{
					Run best_s, best_t;
					best_s.prepareTime = best_s.solveTime = 1e30;
					best_t.prepareTime = best_t.solveTime = 1e30;

					for ( int r = 0; r < repeats; ++r )
					{
						Run s = ( pass == 0 )
							? runOnce( box(), order, n, example5, datum, guessFromDatum, AM::Serial, trace )
							: runOnce( box(), order, n, pedestal, datum, guessFromDatum, AM::Serial, trace );
						if ( s.prepareTime < best_s.prepareTime ) best_s.prepareTime = s.prepareTime;
						if ( s.solveTime < best_s.solveTime ) best_s.solveTime = s.solveTime;
						best_s.newtonIterations = s.newtonIterations;
						best_s.converged = s.converged;
						if ( best_s.psi.empty() ) { best_s.psi = s.psi; best_s.flux = s.flux; }
					}

					if ( canThread )
					{
						for ( int r = 0; r < repeats; ++r )
						{
							Run t = ( pass == 0 )
								? runOnce( box(), order, n, example5, datum, guessFromDatum, AM::Threaded, trace )
								: runOnce( box(), order, n, pedestal, datum, guessFromDatum, AM::Threaded, trace );
							if ( t.prepareTime < best_t.prepareTime ) best_t.prepareTime = t.prepareTime;
							if ( t.solveTime < best_t.solveTime ) best_t.solveTime = t.solveTime;
							best_t.newtonIterations = t.newtonIterations;
							best_t.converged = t.converged;
							if ( best_t.psi.empty() ) { best_t.psi = t.psi; best_t.flux = t.flux; }
						}
					}

					std::printf( "    %-10s %2d %5d %8s %9.4f %9.4f %8.2f %9.4f %9.4f %8.2f %2d/%-2d",
					             name, order, n, traceName( trace ),
					             best_s.prepareTime,
					             canThread ? best_t.prepareTime : 0.0,
					             canThread && best_t.prepareTime > 0.0
					               ? best_s.prepareTime/best_t.prepareTime : 0.0,
					             best_s.solveTime,
					             canThread ? best_t.solveTime : 0.0,
					             canThread && best_t.solveTime > 0.0
					               ? best_s.solveTime/best_t.solveTime : 0.0,
					             best_s.newtonIterations,
					             canThread ? best_t.newtonIterations : 0 );

					if ( !best_s.converged || ( canThread && !best_t.converged ) )
						std::printf( "  (did not converge)" );
					std::printf( "\n" );
					std::fflush( stdout );

					if ( canThread && best_s.converged && best_t.converged )
					{
						worstThreadedPsi = std::max( worstThreadedPsi,
							worstDifference( best_s.psi, best_t.psi ) );
						worstThreadedFlux = std::max( worstThreadedFlux,
							worstDifference( best_s.flux, best_t.flux ) );
						if ( best_s.newtonIterations != best_t.newtonIterations )
							++newtonMismatches;
					}
				}

				// The trace solvers, against the first available one, on the
				// serial path so that only one variable moves.
				if ( traceSolvers.size() > 1 )
				{
					Run reference = ( pass == 0 )
						? runOnce( box(), order, n, example5, datum, guessFromDatum, AM::Serial, traceSolvers[ 0 ] )
						: runOnce( box(), order, n, pedestal, datum, guessFromDatum, AM::Serial, traceSolvers[ 0 ] );
					for ( size_t i = 1; i < traceSolvers.size(); ++i )
					{
						Run other = ( pass == 0 )
							? runOnce( box(), order, n, example5, datum, guessFromDatum, AM::Serial, traceSolvers[ i ] )
							: runOnce( box(), order, n, pedestal, datum, guessFromDatum, AM::Serial, traceSolvers[ i ] );
						if ( reference.converged && other.converged )
						{
							double scale = 0.0;
							for ( double v : reference.psi )
								scale = std::max( scale, std::fabs( v ) );
							if ( scale > 0.0 )
								worstTraceAgreement = std::max( worstTraceAgreement,
									worstDifference( reference.psi, other.psi )/scale );
						}
					}
				}
			}
		}
	}

	std::printf( "\n  correctness, which is what makes the timings mean anything\n" );
	if ( canThread )
	{
		std::printf( "    threaded vs serial, psi     : worst difference %.3e%s\n",
		             worstThreadedPsi, worstThreadedPsi == 0.0 ? "  (exact)" : "  *** NOT EXACT ***" );
		std::printf( "    threaded vs serial, flux    : worst difference %.3e%s\n",
		             worstThreadedFlux, worstThreadedFlux == 0.0 ? "  (exact)" : "  *** NOT EXACT ***" );
		std::printf( "    threaded vs serial, Newton  : %d case(s) took a different "
		             "number of iterations%s\n", newtonMismatches,
		             newtonMismatches == 0 ? "  (exact)" : "  *** DIFFERS ***" );
	}
	else
		std::printf( "    threaded assembly           : unavailable in this build\n" );

	if ( traceSolvers.size() > 1 )
		std::printf( "    trace solvers agree to      : %.3e relative\n",
		             worstTraceAgreement );

	// The exit code. Only correctness, never a timing.
	if ( canThread && worstThreadedPsi != 0.0 )
	{
		std::printf( "\n  *** threaded assembly is NOT bit for bit in psi (%.3e).\n"
		             "      MFEM documents exactness on both its threaded loops, so\n"
		             "      this is a real change rather than a tolerance to widen.\n"
		             "      On a NONLINEAR source the loop is MultNL's, and the\n"
		             "      caller's own integrators sit on it: check whether\n"
		             "      anything MEQ installs holds per-point scratch as a\n"
		             "      member. meq::SourceIntegrator's `shape` is guarded on\n"
		             "      MFEM_THREAD_SAFE for exactly that reason.\n",
		             worstThreadedPsi );
		++failures;
	}
	if ( canThread && worstThreadedFlux != 0.0 )
	{
		std::printf( "\n  *** threaded assembly is NOT bit for bit in the flux (%.3e).\n",
		             worstThreadedFlux );
		++failures;
	}
	if ( canThread && newtonMismatches > 0 )
	{
		std::printf( "\n  *** %d case(s) where threaded assembly changed the NEWTON\n"
		             "      ITERATION COUNT while reaching the same answer. That is\n"
		             "      the signature of a corrupted JACOBIAN rather than a\n"
		             "      corrupted residual -- Newton reaches the same root by a\n"
		             "      longer path, so no error norm and no convergence rate can\n"
		             "      see it. AssembleElementGrad is where to look.\n",
		             newtonMismatches );
		++failures;
	}
	if ( traceSolvers.size() > 1 && worstTraceAgreement > 1.0e-10 )
	{
		std::printf( "\n  *** the trace solvers disagree by %.3e, which is past the\n"
		             "      1e-10 that makes setTraceSolver() a performance choice.\n",
		             worstTraceAgreement );
		++failures;
	}

	std::printf( "\n" );
	return failures == 0 ? 0 : 1;
}
