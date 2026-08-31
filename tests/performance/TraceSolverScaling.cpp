/*
 * TraceSolverScaling -- the performance harness meq did not have.
 *
 * This is a BENCHMARK, not a test: it is built but not registered with ctest,
 * because everything it reports is a timing and CLAUDE.md's standing rule is
 * that a threaded timing on this machine is a measurement about the machine.
 * Asserting on one would make the suite fail for reasons that have nothing to
 * do with the code. What it *does* assert is the two correctness properties
 * that make the timings meaningful at all -- that threaded assembly reproduces
 * serial assembly bit for bit, and that PARDISO reproduces UMFPACK -- and it
 * says so loudly if either breaks, because a faster wrong answer is not a
 * result.
 *
 * Three things are measured, and they are three different questions:
 *
 *   1. ASSEMBLY, serial against threaded. MFEM gained
 *      DarcyHybridization::SetAssemblyMode(); the element-local half of
 *      ComputeH() threads and the scatter into the trace matrix cannot. This
 *      measures what that is worth on meq's own problem rather than on the
 *      library's benchmark, and the serial scatter is the Amdahl ceiling.
 *
 *   2. THE TRACE SOLVE, UMFPACK against MKL PARDISO. Both are direct. They
 *      reach the same answer or one of them is being handed a different matrix.
 *
 *   3. THREADS. The scan is driven from outside, by running this binary once
 *      per thread count, because MKL fixes its threading at first use and an
 *      in-process sweep would measure the first setting several times.
 *
 * Run it through tests/performance/scan.sh, which sets the environment each
 * thread count needs and collects the table.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "analytic/Soloviev.hpp"

#ifdef MFEM_USE_OPENMP
#include <omp.h>
#endif

#ifdef MFEM_USE_CUDA
#include <cuda_runtime.h>
#endif

namespace
{
	using Clock = std::chrono::steady_clock;

	double seconds( Clock::time_point start )
	{
		return std::chrono::duration<double>( Clock::now() - start ).count();
	}

	/*
	 * Stop the clock only when the GPU has actually finished.
	 *
	 * WITHOUT THIS EVERY DEVICE NUMBER IN THIS FILE IS FICTION, and the fiction
	 * is convincing rather than obviously broken: cuDSS queues its work on a
	 * stream and returns, so a host-side timer around SetOperator() measures the
	 * enqueue and not the factorisation. Measured, the difference is not
	 * marginal -- a warm cuDSS setup timed without a sync reads 2.0e-04 s for a
	 * factorisation of 148,224 unknowns, which is a millisecond-scale answer to
	 * a second-scale question, and the cost merely reappears somewhere later and
	 * unattributed. The CPU solvers are synchronous, so this is a no-op for them
	 * and is called unconditionally rather than only on the device path.
	 */
	void syncDevice()
	{
#ifdef MFEM_USE_CUDA
		cudaDeviceSynchronize();
#endif
	}

	/// The benchmark box, and a triangulated mesh on it. Duplicated from
	/// tests/convergence/ConvergenceHarness.hpp rather than included, because
	/// that header pulls in Boost.Test and this binary is not a test. Ten lines
	/// of mesh construction is a smaller cost than linking a test framework into
	/// a benchmark, and it keeps the two independent -- the harness is free to
	/// grow assertions without changing what this measures.
	struct Box
	{
		double rMin, rMax, zMin, zMax;
	};

	Box box()
	{
		return Box{ 0.6, 1.4, -0.6, 0.6 };
	}

	mfem::Mesh makeMesh( Box const &b, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE,
		                                               false, b.rMax - b.rMin,
		                                               b.zMax - b.zMin );
		double const rMin = b.rMin;
		double const zMin = b.zMin;
		mesh.Transform( [ rMin, zMin ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	/// Best of `repeats`, which is the right statistic for a timing whose noise
	/// is one-sided: a run can be slowed by the machine and cannot be sped up by
	/// it, so the minimum is the closest thing to the cost of the work.
	template <typename F>
	double bestOf( int repeats, F &&f )
	{
		double best = 1.0e300;
		for ( int i = 0; i < repeats; ++i )
		{
			auto start = Clock::now();
			f();
			syncDevice();
			best = std::min( best, seconds( start ) );
		}
		return best;
	}

	/// Every entry of two assembled trace matrices, compared exactly. Returns
	/// the largest absolute difference, or -1 if the two do not even have the
	/// same structure -- which is a different and worse failure.
	double compareExactly( mfem::SparseMatrix const &a, mfem::SparseMatrix const &b )
	{
		if ( a.Height() != b.Height() || a.NumNonZeroElems() != b.NumNonZeroElems() )
			return -1.0;

		int const nnz = a.NumNonZeroElems();
		double const *av = a.GetData();
		double const *bv = b.GetData();
		int const *aj = a.GetJ();
		int const *bj = b.GetJ();

		double worst = 0.0;
		for ( int i = 0; i < nnz; ++i )
		{
			if ( aj[ i ] != bj[ i ] )
				return -1.0;
			worst = std::max( worst, std::fabs( av[ i ] - bv[ i ] ) );
		}
		return worst;
	}

	struct Built
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		mfem::SparseMatrix *matrix = nullptr;
	};

	/// Mesh, spaces and coefficients, but NOT assembled. Kept out of the timed
	/// region deliberately: building a Mesh and three FiniteElementSpaces is
	/// serial work that no assembly mode touches, and leaving it in would dilute
	/// the very ratio this benchmark exists to report.
	Built prepareUnassembled( int order, int n, mfem::Coefficient &source,
	                          mfem::Coefficient &psi )
	{
		Built out;
		out.mesh = std::make_unique<mfem::Mesh>( makeMesh( box(), n ) );
		out.solver = std::make_unique<meq::GradShafranovSolver>( *out.mesh, order );
		out.solver->setSource( source );
		out.solver->setBoundaryData( psi );
		return out;
	}

	/// Assemble and reduce, in one mode. setAssemblyMode() clears the solver's
	/// `built` flag whatever mode is asked for, which is what makes prepare()
	/// re-runnable and so what makes a best-of-N timing possible at all.
	void assembleInto( Built &b, meq::GradShafranovSolver::AssemblyMode mode )
	{
		b.solver->setAssemblyMode( mode );
		b.solver->prepare();
		b.matrix = dynamic_cast<mfem::SparseMatrix *>( &b.solver->reducedOperator() );
	}
}

int main( int argc, char **argv )
{
	int repeats = 3;
	std::string device = "cpu";
	std::vector<int> orders = { 2, 3 };
	std::vector<int> sizes = { 16, 32, 48, 64 };
	bool assemblyTiming = true;
	bool reuse = false;

	// Comma-separated lists, so a scan can restrict itself. That is not a
	// convenience: under threaded MKL the element-local dense factorisation at
	// k = 3 degrades by a factor of forty, which makes a full sweep cost hours
	// and buries the solver columns it was run for. Restricting to k = 2 is how
	// the direct-solver scan gets measured at all.
	auto parseList = []( char const *text )
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
	};

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
		else if ( arg == "--device" && i + 1 < argc )
			device = argv[ ++i ];
		else if ( arg == "--orders" && i + 1 < argc )
			orders = parseList( argv[ ++i ] );
		else if ( arg == "--sizes" && i + 1 < argc )
			sizes = parseList( argv[ ++i ] );
		else if ( arg == "--no-assembly-timing" )
			assemblyTiming = false;
		else if ( arg == "--reuse" )
			reuse = true;
	}

	/*
	 * The device, if one was asked for. CONSTRUCTED FIRST AND LEFT ALIVE for the
	 * whole run, which mfem::Device requires -- it is global state that decides
	 * where every Vector afterwards allocates.
	 *
	 * It is a flag rather than a default because it is not a free choice. cuDSS
	 * reads its matrix through SparseMatrix::ReadI/ReadJ/ReadData and its vectors
	 * through Read()/Write(), which return DEVICE pointers only when a device
	 * backend is configured -- so cuDSS cannot be measured without it. But
	 * configuring one changes where the rest of meq allocates too, and meq's HDG
	 * assembly has no device kernels, so the CPU columns taken under "cuda" would
	 * be CPU columns paying a memory manager they do not use. Two runs, two
	 * configurations, and the comparison between UMFPACK and cuDSS is made
	 * WITHIN one of them rather than across both.
	 */
	std::unique_ptr<mfem::Device> deviceHandle;
	if ( device != "cpu" )
	{
#ifndef MFEM_USE_CUDA
		// Exit 0, not 1. Asking for a device in a build that has none is a
		// SKIP: this binary doubles as the `cuDSSTraceSolver` ctest, which must
		// pass on a machine with no GPU rather than fail for having none.
		std::printf( "\n  device \"%s\" requested, but this MFEM has no CUDA -- "
		             "skipped\n\n", device.c_str() );
		return 0;
#endif
#ifndef MFEM_USE_CUDSS
		std::printf( "\n  note: CUDA is on but MFEM_USE_CUDSS is off, so the cuDSS "
		             "columns are absent and this run checks only the CPU solvers\n" );
#endif
		deviceHandle = std::make_unique<mfem::Device>( device.c_str() );
		std::printf( "\n" );
		deviceHandle->Print();
	}

	// What the environment actually gave us, printed rather than assumed: the
	// whole point of the scan is that these differ between runs, and a table
	// with the requested thread count in it rather than the delivered one is
	// how a scan comes to describe a machine it never ran on.
	char const *ompEnv = std::getenv( "OMP_NUM_THREADS" );
	char const *mklEnv = std::getenv( "MKL_NUM_THREADS" );
	char const *layer  = std::getenv( "MKL_THREADING_LAYER" );

	int ompMax = 1;
#ifdef MFEM_USE_OPENMP
	ompMax = omp_get_max_threads();
#endif

	std::printf( "\n=== meq trace-solver scaling ===\n" );
	std::printf( "  OMP_NUM_THREADS=%-6s MKL_NUM_THREADS=%-6s MKL_THREADING_LAYER=%s\n",
	             ompEnv ? ompEnv : "(unset)", mklEnv ? mklEnv : "(unset)",
	             layer ? layer : "(unset)" );
	std::printf( "  omp_get_max_threads() = %d, best of %d, device \"%s\", setup is %s\n",
	             ompMax, repeats, device.c_str(),
	             reuse ? "WARM (symbolic/reordering reused)" : "cold" );
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
#ifdef MFEM_USE_CUDA
	std::printf( " CUDA" );
#endif
#ifdef MFEM_USE_CUDSS
	std::printf( " CUDSS" );
#endif
#ifdef MFEM_USE_SUITESPARSE
	std::printf( " SUITESPARSE" );
#endif
	std::printf( "\n" );

	meq::analytic::SolovievEquilibrium const eq =
		meq::analytic::SolovievEquilibrium::nstx();
	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psi( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	using AM = meq::GradShafranovSolver::AssemblyMode;

	bool threadedAvailable = true;
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	threadedAvailable = false;
#endif

	std::printf( "\n  assembly + reduction, and the trace solve\n" );
	std::printf( "    %2s %5s %8s %9s %10s %10s %7s %11s %11s %11s %11s"
#ifdef MFEM_USE_CUDSS
	             " %11s %11s"
#endif
	             "\n",
	             "k", "n", "trace", "nnz/row", "asm ser", "asm thr", "speedup",
	             "UMF setup", "UMF solve", "PAR setup", "PAR solve"
#ifdef MFEM_USE_CUDSS
	             , "cuDSS setup", "cuDSS solve"
#endif
	             );

	int structureFailures = 0;
	int agreementFailures = 0;
	double worstAssemblyDiff = 0.0;
	double worstSolverDiff = 0.0;
	double worstCudssDiff = -1.0;

	for ( int order : orders )
	{
		for ( int n : sizes )
		{
			Built serial = prepareUnassembled( order, n, source, psi );
			assembleInto( serial, AM::Serial );
			if ( serial.matrix == nullptr )
			{
				std::printf( "    k=%d n=%d: reduced operator is not a SparseMatrix\n",
				             order, n );
				continue;
			}

			int const size = serial.matrix->Height();
			double const nnzPerRow =
				static_cast<double>( serial.matrix->NumNonZeroElems() )/size;

			double serialAsm = assemblyTiming
				? bestOf( repeats, [ & ]() { assembleInto( serial, AM::Serial ); } )
				: -1.0;

			double threadedAsm = -1.0;
			double assemblyDiff = 0.0;
			if ( threadedAvailable )
			{
				// A SEPARATE solver for the threaded matrix, so both are live at
				// once and can be compared entry by entry. Comparing against a copy
				// taken earlier from the same solver would not do: prepare()
				// reassembles in place, so the "serial" matrix would be whatever the
				// last timed repeat left behind.
				Built threaded = prepareUnassembled( order, n, source, psi );
				assembleInto( threaded, AM::Threaded );
				assemblyDiff = compareExactly( *serial.matrix, *threaded.matrix );
				if ( assemblyDiff < 0.0 )
					++structureFailures;
				else
					worstAssemblyDiff = std::max( worstAssemblyDiff, assemblyDiff );

				if ( assemblyTiming )
					threadedAsm = bestOf( repeats, [ & ]()
					{
						assembleInto( threaded, AM::Threaded );
					} );
			}

			// Back to serial, so the matrix the solvers below are handed is the
			// one this row's UMFPACK/PARDISO columns claim to be about.
			assembleInto( serial, AM::Serial );

			mfem::Vector const &rhs = serial.solver->reducedRhs();

			// ---- UMFPACK, configured exactly as GradShafranov.cpp does ----
			mfem::Vector umfX( size );
			umfX = 0.0;
			double umfSetup = 0.0, umfSolve = 0.0;
			{
				mfem::UMFPackSolver umf;
				umf.Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
				// With --reuse the reported setup is the WARM one: the symbolic
				// analysis is kept and only the numeric factorisation is redone.
				// That is the steady-state cost of a Newton step, which re-forms
				// the same sparsity every iterate, and is the number to compare
				// across solvers for meq's actual workload. One cold call first,
				// outside the timing, so the analysis really is already done.
				if ( reuse )
				{
					umf.SetReuseSymbolic( true );
					umf.SetOperator( *serial.matrix );
				}
				umfSetup = bestOf( repeats, [ & ]() { umf.SetOperator( *serial.matrix ); } );
				umfSolve = bestOf( repeats, [ & ]() { umf.Mult( rhs, umfX ); } );
			}

			double parSetup = -1.0, parSolve = -1.0, agreement = -1.0;
#ifdef MFEM_USE_MKL_PARDISO
			{
				mfem::Vector parX( size );
				parX = 0.0;
				mfem::PardisoSolver par;
				// Structurally symmetric, not symmetric in value: on the extension
				// path HDGExtensionIntegrator puts an outer product into the flux
				// block and the values lose symmetry, while the sparsity keeps it.
				par.SetMatrixType( mfem::PardisoSolver::REAL_STRUCTURE_SYMMETRIC );
				par.SetPrintLevel( 0 );
				if ( reuse )
				{
					par.SetReuseSymbolic( true );
					par.SetOperator( *serial.matrix );
				}
				parSetup = bestOf( repeats, [ & ]() { par.SetOperator( *serial.matrix ); } );
				parSolve = bestOf( repeats, [ & ]() { par.Mult( rhs, parX ); } );

				mfem::Vector diff( parX );
				diff -= umfX;
				agreement = diff.Norml2()/std::max( 1.0e-300, umfX.Norml2() );
				worstSolverDiff = std::max( worstSolverDiff, agreement );
				if ( agreement >= 1.0e-10 )
					++agreementFailures;
			}
#endif

			double cudssSetup = -1.0, cudssSolve = -1.0;
#ifdef MFEM_USE_CUDSS
			{
				// NONSYMMETRIC and FULL, deliberately, and not because the matrix
				// always is. On a fitted mesh the trace matrix is symmetric to
				// 2e-16 and negative definite; on the extension path -- meq's
				// headline configuration -- HDGExtensionIntegrator puts an outer
				// product into the flux block and the relative asymmetry is 5.4e-1.
				// The sparsity is symmetric on both. So the general setting is the
				// only correct one, and it is what PARDISO's
				// REAL_STRUCTURE_SYMMETRIC amounts to as well.
				mfem::Vector cudssX( size );
				cudssX = 0.0;
				mfem::CuDSSSolver cudss;
				cudss.SetMatrixSymType( mfem::CuDSSSolver::NONSYMMETRIC );
				cudss.SetMatrixViewType( mfem::CuDSSSolver::FULL );
				// cuDSS spells it SetReorderingReuse, and requires it BEFORE the
				// first SetOperator -- it verifies against Ac == nullptr.
				if ( reuse )
					cudss.SetReorderingReuse( true );
				if ( reuse )
					cudss.SetOperator( *serial.matrix );
				cudssSetup = bestOf( repeats, [ & ]()
				{
					cudss.SetOperator( *serial.matrix );
				} );
				cudssSolve = bestOf( repeats, [ & ]() { cudss.Mult( rhs, cudssX ); } );

				mfem::Vector diff( cudssX );
				diff -= umfX;
				double const rel = diff.Norml2()/std::max( 1.0e-300, umfX.Norml2() );
				worstCudssDiff = std::max( worstCudssDiff, rel );
				if ( rel >= 1.0e-10 )
					++agreementFailures;
			}
#endif

			std::printf( "    %2d %5d %8d %9.1f %10.4f %10.4f %7s %11.4f %11.4f %11.4f %11.4f"
#ifdef MFEM_USE_CUDSS
			             " %11.4f %11.4f"
#endif
			             "\n",
			             order, n, size, nnzPerRow, serialAsm, threadedAsm,
			             threadedAsm > 0.0
			                 ? ( std::to_string( serialAsm/threadedAsm ).substr( 0, 5 ) + "x" ).c_str()
			                 : "-",
			             umfSetup, umfSolve, parSetup, parSolve
#ifdef MFEM_USE_CUDSS
			             , cudssSetup, cudssSolve
#endif
			             );
			std::fflush( stdout );
		}
	}

	std::printf( "\n  correctness, which is what makes the timings mean anything\n" );
	if ( threadedAvailable )
		std::printf( "    threaded assembly vs serial : worst entry difference %.3e%s\n",
		             worstAssemblyDiff,
		             worstAssemblyDiff == 0.0 ? "   (bit for bit, as documented)" : "" );
	else
		std::printf( "    threaded assembly           : unavailable in this build\n" );
#ifdef MFEM_USE_MKL_PARDISO
	std::printf( "    PARDISO vs UMFPACK          : worst relative difference %.3e\n",
	             worstSolverDiff );
#endif
#ifdef MFEM_USE_CUDSS
	if ( worstCudssDiff >= 0.0 )
		std::printf( "    cuDSS vs UMFPACK            : worst relative difference %.3e\n",
		             worstCudssDiff );
#endif

	if ( structureFailures > 0 )
	{
		std::printf( "\n  *** %d case(s) where threaded assembly changed the SPARSITY, "
		             "not merely the values. That is not a rounding question.\n",
		             structureFailures );
		return 1;
	}
	if ( worstAssemblyDiff != 0.0 )
	{
		std::printf( "\n  *** threaded assembly is NOT bit for bit (%.3e). MFEM "
		             "documents that it is, so this is a real regression.\n",
		             worstAssemblyDiff );
		return 1;
	}
	if ( agreementFailures > 0 )
	{
		std::printf( "\n  *** %d case(s) where PARDISO and UMFPACK disagree. Check the "
		             "matrix type before the timings.\n", agreementFailures );
		return 1;
	}

	std::printf( "\n" );
	return 0;
}
