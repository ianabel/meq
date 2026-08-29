#define BOOST_TEST_MODULE TraceSolverComparison
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * THE GLOBAL TRACE SOLVE: UMFPACK AGAINST PARDISO.
 *
 * meq solves the hybridized trace system with UMFPackSolver and
 * UMFPACK_ORDERING_METIS. CLAUDE.md records an earlier attempt to compare it
 * against MKL PARDISO which had to be abandoned: Debian's intel-mkl 2020.4.304
 * returned error -3 above n ~ 3000, demonstrated on a plain 5-point Laplacian
 * with no MFEM linked at all, so the verdict was on the packaging and not on the
 * solver. That file left the question open with a list of what a proper look
 * would need. This is that look, against oneAPI MKL 2026.1.
 *
 * WHAT IS ASSERTED AND WHAT IS ONLY PRINTED, because the difference matters
 * here more than usual. The two solvers must agree -- that is a statement about
 * the code and it is asserted. The timings are NOT asserted, for the reason
 * CLAUDE.md gives twice over: PARDISO is threaded, threaded reduction order
 * depends on the thread count, and a measurement taken on this machine is not a
 * measurement about the code. They are printed so a reader can see them and
 * repeat them somewhere that matters.
 *
 * BOTH PATHS SHIP. PARDISO needs oneMKL, whose licence is not everybody's to
 * accept, and MFEM_USE_MKL_PARDISO is off in most builds. So this file is
 * written to compile and pass either way: without it, the PARDISO rows are
 * skipped and the UMFPACK ones still run.
 */

namespace
{
	using meq::tests::EquilibriumSource;

	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}

	double seconds( std::chrono::steady_clock::time_point start )
	{
		return std::chrono::duration<double>( std::chrono::steady_clock::now()
		                                      - start ).count();
	}

	/// The assembled trace system for one mesh and degree, on the LINEAR path so
	/// that the operator really is a matrix rather than a residual.
	struct TraceSystem
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		mfem::SparseMatrix *matrix = nullptr;
	};

	TraceSystem assemble( int order, int n, meq::analytic::SolovievEquilibrium const &eq,
	                      mfem::Coefficient &sourceCoeff, mfem::Coefficient &psiCoeff )
	{
		(void)eq;
		TraceSystem system;
		system.mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( box(), n ) );
		system.solver = std::make_unique<meq::GradShafranovSolver>( *system.mesh, order );
		system.solver->setSource( sourceCoeff );
		system.solver->setBoundaryData( psiCoeff );
		system.solver->prepare();
		system.matrix = dynamic_cast<mfem::SparseMatrix *>(
			&system.solver->reducedOperator() );
		return system;
	}
}

BOOST_AUTO_TEST_CASE( umfpackAndPardisoAgreeOnTheTraceSystem )
{
	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

#ifdef MFEM_USE_MKL_PARDISO
	std::printf( "\n  the hybridized trace system: UMFPACK against MKL PARDISO\n" );
#else
	std::printf( "\n  the hybridized trace system: UMFPACK only "
	             "(MFEM_USE_MKL_PARDISO is off)\n" );
#endif
	std::printf( "    %2s %6s %8s %9s %11s %11s %11s %11s\n",
	             "k", "n", "trace", "nnz/row", "UMF setup", "UMF solve",
	             "PAR setup", "PAR solve" );

	// Four sizes, the largest past where Debian's MKL used to fail with error -3.
	for ( int order : { 2, 3 } )
	{
		for ( int n : { 16, 32, 48 } )
		{
			TraceSystem system = assemble( order, n, eq, sourceCoeff, psiCoeff );
			BOOST_TEST_REQUIRE( system.matrix != nullptr,
			                    "the reduced operator is not a SparseMatrix, so "
			                    "neither direct solver can be handed it" );

			int const size = system.matrix->Height();
			double const nnzPerRow =
				static_cast<double>( system.matrix->NumNonZeroElems() )/size;

			mfem::Vector const &rhs = system.solver->reducedRhs();

			// ---- UMFPACK, configured exactly as GradShafranov.cpp does ----
			mfem::Vector umfpackX( size );
			umfpackX = 0.0;
			double umfpackSetup = 0.0, umfpackSolve = 0.0;
			{
				mfem::UMFPackSolver umfpack;
				umfpack.Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
				auto start = std::chrono::steady_clock::now();
				umfpack.SetOperator( *system.matrix );
				umfpackSetup = seconds( start );

				start = std::chrono::steady_clock::now();
				umfpack.Mult( rhs, umfpackX );
				umfpackSolve = seconds( start );
			}

			double pardisoSetup = -1.0, pardisoSolve = -1.0;
			double agreement = -1.0;
#ifdef MFEM_USE_MKL_PARDISO
			{
				mfem::Vector pardisoX( size );
				pardisoX = 0.0;
				mfem::PardisoSolver pardiso;
				// STRUCTURE symmetric, not symmetric: CLAUDE.md measures the trace
				// matrix as symmetric to 2e-16 on a FITTED mesh and asymmetric at
				// 5.4e-1 on the extension path, where HDGExtensionIntegrator puts
				// an outer product into the flux block. The sparsity is symmetric
				// on both. Type 1 is therefore the one that is right for meq's
				// headline configuration as well as this one.
				pardiso.SetMatrixType( mfem::PardisoSolver::REAL_STRUCTURE_SYMMETRIC );
				pardiso.SetPrintLevel( 0 );

				auto start = std::chrono::steady_clock::now();
				pardiso.SetOperator( *system.matrix );
				pardisoSetup = seconds( start );

				start = std::chrono::steady_clock::now();
				pardiso.Mult( rhs, pardisoX );
				pardisoSolve = seconds( start );

				mfem::Vector difference( pardisoX );
				difference -= umfpackX;
				agreement = difference.Norml2()
				            /std::max( 1.0e-300, umfpackX.Norml2() );
			}
#endif

			std::printf( "    %2d %6d %8d %9.1f %11.4f %11.4f %11.4f %11.4f",
			             order, n, size, nnzPerRow, umfpackSetup, umfpackSolve,
			             pardisoSetup, pardisoSolve );
			if ( agreement >= 0.0 )
				std::printf( "   agree to %.2e", agreement );
			std::printf( "\n" );
			std::fflush( stdout );

#ifdef MFEM_USE_MKL_PARDISO
			// THE ONE ASSERTION. Two direct solvers on the same matrix must reach
			// the same answer; anything else is one of them being handed a
			// different matrix, or a matrix type that misdescribes it.
			BOOST_TEST( agreement < 1.0e-10,
			            "k = " << order << ", n = " << n << ": PARDISO and UMFPACK "
			            "differ by " << agreement << " relative on a trace system of "
			            << size << " unknowns. Check the matrix type first -- meq's "
			            "trace matrix is structurally symmetric and, on the "
			            "extension path, not symmetric in its values" );
#endif
		}
	}
}
