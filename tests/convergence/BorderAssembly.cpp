#define BOOST_TEST_MODULE BorderAssembly
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"

// omp_get_max_threads() and omp_set_num_threads(). Guarded on the same pair
// src/meq/Threading.hpp is: with either absent meq's own element loops are
// serial by construction and there is nothing to compare.
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
#include <omp.h>
#endif

#include "DivertedMachine.hpp"

using namespace meqtest;

/*
 * ===========================================================================
 * THE BORDER'S ELEMENT LOOPS: COMPRESSED EXACTLY, AND THREADED EXACTLY
 * ===========================================================================
 *
 * THREADING-PLAN.md items A and B, on the one fixture in this tree that
 * carries every border row at once -- psi_ax, psi_bnd, XP-3's two, the current
 * scale and the exterior modes. examples/diverted-tokamak.toml is freegs4e's
 * A_testtokamak_classic, and DivertedMachine.hpp assembles it exactly as
 * apps/meq.cpp does.
 *
 * WHAT THE TWO ITEMS HAVE IN COMMON IS THE ACCEPTANCE AND NOT THE MECHANISM.
 * Neither is allowed to move a digit:
 *
 *   * item A compresses eleven of the fifteen border rows onto the indices
 *     they are actually non-zero on, and dropping a term that is EXACTLY zero
 *     leaves every partial sum bit-unchanged. So the assertion is
 *     `0.000e+00` against the dense contraction, not a tolerance.
 *   * item B threads five plasma element loops, and three of them write into
 *     an L2 space where no two elements share a dof. So the assertion is
 *     `0.000e+00` between one thread and eight, not a tolerance.
 *
 * AND THE TWO SCALAR ACCUMULATORS ARE THE EXCEPTION THAT IS WORTH STATING.
 * assemblePlasmaCurrent() and assembleCurrentNormalisationCorner() carried one
 * running total across every quadrature point of every element, and a flat
 * left fold is not reproducible from per-element partials by any grouping --
 * `( a + b ) + ( c + d )` is not `( ( a + b ) + c ) + d`. THREADING-PLAN.md
 * item B says per-element partials are "bit-identical to the serial loop" and
 * that is false; what IS true, and what is asserted below, is that they give
 * the same bits at every thread count. The re-association against the
 * PREVIOUS serial code is a one-off of order 1e-16 and is recorded rather than
 * hidden.
 */

namespace
{
	/// A deterministic filling of a state vector, so that two runs contract
	/// against the same numbers without either depending on a solve.
	///
	/// NOT std::rand and not a single constant. A constant column makes a dot
	/// product insensitive to which entries it visited -- the very thing item A
	/// changes -- so the values have to vary with the index; and a library
	/// generator would make the case depend on somebody else's sequence.
	void fillDeterministic( mfem::Vector &v, double seed )
	{
		for ( int i = 0; i < v.Size(); ++i )
			v( i ) = std::sin( seed + 0.7*static_cast<double>( i ) )
			         *( 1.0 + 0.001*static_cast<double>( i % 17 ) );
	}

	/// The dense contraction the compressed one replaces, written out here so
	/// that the comparison is against an independent loop rather than against
	/// another call into the thing being tested.
	double denseDot( mfem::Vector const &row, mfem::Vector const &v )
	{
		double total = 0.0;
		int const n = std::min( row.Size(), v.Size() );
		for ( int i = 0; i < n; ++i )
			total += row( i )*v( i );
		return total;
	}

	/// The three block sizes, in the order `solution` holds them: flux,
	/// potential, trace. GradShafranovSolver::blockOffsets is private, and this
	/// layout is what its own documentation states.
	int solutionSize( meq::GradShafranovSolver &solver )
	{
		return solver.fluxSpace().GetVSize()
		       + solver.potentialSpace().GetVSize()
		       + solver.traceSpace().GetVSize();
	}
}

/*
 * ITEM A. THE COMPRESSED ROWS CONTRACT EXACTLY, AND THE SUPPORT IS CHECKED
 * RATHER THAN TRUSTED.
 *
 * Two claims, and they are different claims. The first is arithmetic: summing
 * the terms on the support in ascending index order reproduces the dense
 * loop's every partial sum, because the terms it drops were exactly `0.0`.
 * The second is bookkeeping: the support really does list every index the
 * scatter wrote. The first is a theorem and the second is a property of the
 * code, so the second is where a defect would live -- which is why
 * exteriorTransmissionRows() checks it before returning and why this case
 * checks it again from outside.
 *
 * **AND THE COMPRESSION HAS TO BE A REAL ONE.** A support list containing
 * every index would pass both assertions and buy nothing, so the ratio is
 * asserted too: the plan sizes it at 984 of about 109,200, a factor of 111.
 */
BOOST_AUTO_TEST_CASE( theCompressedBorderRowsContractExactly )
{
	Machine m = buildMachine();
	completeMachine( m );
	m.solver->prepare();

	std::vector<int> support;
	std::vector<mfem::Vector> const rows =
		m.solver->exteriorTransmissionRows( *m.exterior, &support );

	BOOST_TEST_REQUIRE( static_cast<int>( rows.size() ) == m.exterior->modeCount(),
	                    "the sweep did not return one row per exterior mode" );
	BOOST_TEST_REQUIRE( !support.empty(),
	                    "the sweep recorded an EMPTY support, so every "
	                    "contraction below would be zero and would agree "
	                    "trivially" );

	int const n = solutionSize( *m.solver );

	std::printf( "\n  THE TRANSMISSION ROWS, COMPRESSED\n" );
	std::printf( "    %zu modes, support %zu of %d, a factor of %.1f\n",
	             rows.size(), support.size(), n,
	             static_cast<double>( n )/static_cast<double>( support.size() ) );
	std::fflush( stdout );

	// The support is a real compression rather than a relabelling of the whole
	// vector. 111 on the plan's arithmetic; 10 is the bar, so a mesh whose
	// Gamma_h is a larger share of the domain does not fail this.
	BOOST_TEST( static_cast<double>( n )/static_cast<double>( support.size() )
	            > 10.0,
	            "the recorded support is " << support.size() << " of " << n
	            << ", which is not a compression worth having" );

	// ASCENDING AND UNIQUE, which is what the exactness argument rests on: the
	// surviving terms have to arrive in the same order the dense loop visited
	// them in.
	for ( std::size_t k = 1; k < support.size(); ++k )
		BOOST_TEST_REQUIRE( support[ k ] > support[ k - 1 ],
		                    "the support is not strictly ascending at " << k );

	// AND IT COVERS EVERY NON-ZERO. This is the bookkeeping half, asked of the
	// rows from outside the function that built them.
	for ( std::size_t mode = 0; mode < rows.size(); ++mode )
		BOOST_TEST_REQUIRE( meq::firstNonzeroOutside( rows[ mode ], support ) < 0,
		                    "mode " << mode
		                    << " is non-zero at an index the recorded support "
		                    "does not list" );

	meq::CompressedRows compressed;
	compressed.compress( rows, support );
	BOOST_TEST_REQUIRE( compressed.rowCount()
	                    == static_cast<int>( rows.size() ) );

	// FIVE COLUMNS, each varying with the index, plus the two degenerate ones
	// -- all zeros, and all ones -- because a compression bug that dropped the
	// whole support would agree with the dense loop on a zero column and a
	// support that listed extra indices would agree on a constant one.
	std::vector<mfem::Vector> columns;
	for ( int c = 0; c < 5; ++c )
	{
		columns.emplace_back( n );
		fillDeterministic( columns.back(), 0.31*static_cast<double>( c + 1 ) );
	}
	columns.emplace_back( n );
	columns.back() = 0.0;
	columns.emplace_back( n );
	columns.back() = 1.0;

	double worst = 0.0;
	std::vector<double> gathered;
	for ( std::size_t mode = 0; mode < rows.size(); ++mode )
	{
		for ( std::size_t c = 0; c < columns.size(); ++c )
		{
			double const dense = denseDot( rows[ mode ], columns[ c ] );
			double const sparse =
				compressed.dot( static_cast<int>( mode ), columns[ c ],
				                gathered );
			worst = std::max( worst, std::abs( dense - sparse ) );

			// ZERO, NOT A TOLERANCE. Exactness is the claim.
			BOOST_TEST( dense - sparse == 0.0,
				"mode " << mode << ", column " << c << ": the dense contraction "
				"gives " << dense << " and the compressed one " << sparse
				<< ", a difference of " << ( dense - sparse )
				<< ". Dropping terms that are EXACTLY zero cannot move a sum, so "
				"a difference here means the support and the row disagree about "
				"which entries are non-zero." );

			// And the pre-gathered route, which is what the elimination uses
			// once per column rather than once per ( row, column ) pair.
			compressed.gather( columns[ c ], gathered );
			double const preGathered =
				compressed.dot( static_cast<int>( mode ), gathered );
			BOOST_TEST( preGathered - sparse == 0.0,
				"mode " << mode << ", column " << c
				<< ": gathering the column once and gathering it per row give "
				"different answers" );
		}
	}

	std::printf( "    worst |dense - compressed| over %zu modes x %zu columns: "
	             "%.3e\n", rows.size(), columns.size(), worst );
	std::fflush( stdout );

	/*
	 * AND THE CURRENT ROW, WHICH IS THE OTHER LONG ONE AND HAS ITS OWN
	 * SUPPORT. It lives on the potential dofs of the elements XP-1's fill
	 * reached, which moves with the plasma, so it is rebuilt every Newton step
	 * where the transmission rows are built once per solve.
	 */
	mfem::Vector state( n );
	fillDeterministic( state, 1.7 );

	mfem::Vector currentRow( n );
	std::vector<int> currentSupport;
	m.solver->assembleCurrentRow( state, currentRow, &currentSupport );

	BOOST_TEST_REQUIRE( !currentSupport.empty(),
	                    "the current row recorded an EMPTY support, so the "
	                    "plasma fill reached nothing and the comparison below "
	                    "is vacuous" );
	BOOST_TEST_REQUIRE( meq::firstNonzeroOutside( currentRow, currentSupport ) < 0,
	                    "the current row is non-zero at an index its recorded "
	                    "support does not list" );

	std::printf( "    current row: support %zu of %d, a factor of %.1f\n",
	             currentSupport.size(), n,
	             static_cast<double>( n )
	             /static_cast<double>( currentSupport.size() ) );
	std::fflush( stdout );

	meq::CompressedRows compressedCurrent;
	compressedCurrent.compress( currentRow, currentSupport );
	for ( std::size_t c = 0; c < columns.size(); ++c )
	{
		double const dense = denseDot( currentRow, columns[ c ] );
		double const sparse = compressedCurrent.dot( 0, columns[ c ], gathered );
		BOOST_TEST( dense - sparse == 0.0,
			"current row, column " << c << ": " << dense << " against " << sparse );
	}
}

/*
 * ITEM B. THE FIVE PLASMA ELEMENT LOOPS GIVE THE SAME NUMBERS AT ONE THREAD
 * AND AT EIGHT.
 *
 * WHY THIS IS NOT A WHOLE-SOLVE PIN. A race is not a deterministic wrong
 * answer: one that does not fire on this mesh at this thread count leaves a
 * solve looking perfect, and the same binary can pass a hundred times and fail
 * on the hundred and first. What CAN be asserted deterministically is the
 * thing the design claims -- that the five assemblers are functions of their
 * arguments and of nothing else -- and the way to ask it is to call them
 * directly. That is why they are public; see the header.
 *
 * THE HAZARD THIS GUARDS IS NOT SUBTLE ONCE IT FIRES. `dofs`, `shape` and
 * `point` were declared ABOVE the element loop. Shared across threads they are
 * exactly the meq::SourceIntegrator defect CLAUDE.md records: CalcShape writes
 * into a buffer another thread is reading, and Vector::SetSize is a
 * REALLOCATION, so the failure mode is memory corruption rather than a stale
 * read. Checked to discriminate: hoisting one of them back out of the body
 * makes this case fail by a wide margin on every one of the five.
 *
 * A FOURTH PIECE OF SCRATCH IS GUARDED HERE THAT NOTHING HAD NAMED.
 * FiniteElementSpace::GetElementDofs( int, Array<int> & ) -- the two-argument
 * overload all five loops used -- writes into FiniteElementSpace::DoFTrans, a
 * MUTABLE MEMBER of the space, on every call. It is the same species as
 * Mesh::GetElementTransformation's shared IsoparametricTransformation, and
 * MFEM's own header says the caller-provided overload is the way to avoid it.
 * For an L2 space the only write is a null pointer, so the race is benign
 * today -- by accident of the space rather than by contract.
 */
BOOST_AUTO_TEST_CASE( theBorderAssemblersAreThreadSafe )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	BOOST_TEST_MESSAGE( "MFEM is built without MFEM_USE_OPENMP and "
	                    "MFEM_THREAD_SAFE together, so meq's element loops are "
	                    "serial by construction and there is nothing to "
	                    "compare" );
#else
	Machine m = buildMachine();
	completeMachine( m );
	m.solver->prepare();

	/*
	 * A REAL PLASMA SUPPORT, because every one of the five skips an element
	 * the fill did not reach -- so without one the loops would run over
	 * everything or over nothing and the uneven trip count the schedule has to
	 * cope with would not be exercised.
	 */
	m.solver->setPlasmaSupportFrozen( false );
	m.solver->refreshPlasmaComponent( *m.carried );
	m.solver->setPlasmaSupportFrozen( true );

	int const n = solutionSize( *m.solver );
	int const fluxSize = m.solver->fluxSpace().GetVSize();
	int const potentialSize = m.solver->potentialSpace().GetVSize();

	// The transferred freegs4e guess in the potential block, which is the only
	// block any of the five reads, and a deterministic filling elsewhere so
	// that a loop reading the wrong block would not read zeros.
	mfem::Vector state( n );
	fillDeterministic( state, 2.9 );
	for ( int i = 0; i < potentialSize; ++i )
		state( fluxSize + i ) = ( *m.carried )( i );

	struct Result
	{
		double current = 0.0;
		double againstAxis = 0.0;
		double againstBoundary = 0.0;
		mfem::Vector column;
		mfem::Vector row;
		mfem::Vector normalisationAxis;
		mfem::Vector normalisationBoundary;
		bool axisAssembled = false;
		bool boundaryAssembled = false;
	};

	auto runAll = [ & ]( int threads, Result &out )
	{
		omp_set_num_threads( threads );
		out.current = m.solver->assemblePlasmaCurrent( state );
		m.solver->assembleCurrentColumn( state, out.column );
		m.solver->assembleCurrentRow( state, out.row );
		m.solver->assembleCurrentNormalisationCorner( state, out.againstAxis,
		                                              out.againstBoundary );
		out.axisAssembled =
			m.solver->assembleNormalisationColumn( state, true,
			                                       out.normalisationAxis );
		out.boundaryAssembled =
			m.solver->assembleNormalisationColumn( state, false,
			                                       out.normalisationBoundary );
	};

	int const ambient = omp_get_max_threads();
	// AT LEAST TWO, whatever the environment says. A case that silently ran
	// both arms at one thread would pass on a machine with OMP_NUM_THREADS=1
	// and assert nothing -- and this project has a recorded instance of
	// exactly that shape, a precondition quietly failing and taking the
	// comparison with it.
	int const many = std::max( ambient, 2 );

	Result serial;
	Result threaded;
	runAll( 1, serial );
	runAll( many, threaded );
	omp_set_num_threads( ambient );

	std::printf( "\n  THE FIVE PLASMA ELEMENT LOOPS AT 1 THREAD AND AT %d\n",
	             many );
	std::printf( "    plasma elements %d of %d, n = %d\n",
	             m.solver->plasmaComponentElements(),
	             m.solver->plasmaCandidateElements(), n );
	std::printf( "    int F/R            %.17e\n", serial.current );
	std::printf( "                       %.17e\n", threaded.current );
	std::printf( "    dI/d(psi_ax)       %.17e\n", serial.againstAxis );
	std::printf( "    dI/d(psi_bnd)      %.17e\n", serial.againstBoundary );
	std::printf( "    analytic columns   axis %s, boundary %s\n",
	             serial.axisAssembled ? "yes" : "no",
	             serial.boundaryAssembled ? "yes" : "no" );
	std::fflush( stdout );

	// The fixture has to reach the loops at all. A plasma of zero elements
	// would make every comparison below a comparison of zeros.
	BOOST_TEST_REQUIRE( m.solver->plasmaComponentElements() > 0,
	                    "the plasma fill reached no elements, so none of the "
	                    "five loops has a body to run" );
	BOOST_TEST_REQUIRE( serial.current != 0.0,
	                    "int F/R is exactly zero, so assemblePlasmaCurrent() "
	                    "accumulated nothing and the comparison is vacuous" );

	BOOST_TEST( serial.current - threaded.current == 0.0,
		"assemblePlasmaCurrent() gives " << serial.current << " at one thread "
		"and " << threaded.current << " at " << many << ", a difference of "
		<< ( serial.current - threaded.current )
		<< ". The per-element partial sums are indexed by ELEMENT and summed in "
		"element order precisely so that this cannot happen; a reduction( + : ) "
		"would leave the association to the runtime." );

	BOOST_TEST( serial.againstAxis - threaded.againstAxis == 0.0,
		"assembleCurrentNormalisationCorner()'s axis entry moved by "
		<< ( serial.againstAxis - threaded.againstAxis ) );
	BOOST_TEST( serial.againstBoundary - threaded.againstBoundary == 0.0,
		"assembleCurrentNormalisationCorner()'s boundary entry moved by "
		<< ( serial.againstBoundary - threaded.againstBoundary ) );

	BOOST_TEST( serial.axisAssembled == threaded.axisAssembled );
	BOOST_TEST( serial.boundaryAssembled == threaded.boundaryAssembled );

	struct Comparison
	{
		char const *name;
		mfem::Vector const *serial;
		mfem::Vector const *threaded;
	};
	std::vector<Comparison> const vectors = {
		{ "assembleCurrentColumn", &serial.column, &threaded.column },
		{ "assembleCurrentRow", &serial.row, &threaded.row },
		{ "assembleNormalisationColumn( axis )", &serial.normalisationAxis,
		  &threaded.normalisationAxis },
		{ "assembleNormalisationColumn( boundary )",
		  &serial.normalisationBoundary, &threaded.normalisationBoundary },
	};

	for ( Comparison const &one : vectors )
	{
		BOOST_TEST_REQUIRE( one.serial->Size() == one.threaded->Size(),
		                    one.name << " returned vectors of different sizes" );

		// AND IT HAS TO HAVE ASSEMBLED SOMETHING. A loop that wrote nothing
		// agrees with itself perfectly.
		double magnitude = 0.0;
		int worstIndex = -1;
		double worst = 0.0;
		for ( int i = 0; i < one.serial->Size(); ++i )
		{
			magnitude = std::max( magnitude, std::abs( ( *one.serial )( i ) ) );
			double const difference =
				( *one.serial )( i ) - ( *one.threaded )( i );
			if ( std::abs( difference ) > std::abs( worst ) )
			{
				worst = difference;
				worstIndex = i;
			}
		}

		std::printf( "    %-40s max |entry| %.3e, worst difference %.3e\n",
		             one.name, magnitude, worst );
		std::fflush( stdout );

		BOOST_TEST_REQUIRE( magnitude > 0.0,
			one.name << " returned an identically zero vector, so agreeing at "
			"every entry says nothing" );

		// ENTRY FOR ENTRY, AT ZERO. The potential space is L2, so no two
		// elements share a dof and each dof's accumulation order is untouched
		// by the threading -- these four are bit-identical to the serial loop
		// as well as to themselves at any thread count.
		BOOST_TEST( worst == 0.0,
			one.name << " differs between one thread and " << many
			<< " by " << worst << " at entry " << worstIndex
			<< " of " << one.serial->Size()
			<< ". The scatter is disjoint by construction -- an L2 potential "
			"dof belongs to exactly one element -- so a difference here is "
			"shared scratch in the loop body and not the arithmetic." );
	}
#endif
}
