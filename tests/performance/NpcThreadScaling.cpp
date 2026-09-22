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
 *                          [--pedestal] [--example5]   the two defaults
 *                          [--soloviev]   the LINEAR-source arm
 *                          [--highbeta]   the BORDERED arm, psi_ax an unknown
 *                          [--trace umfpack|pardiso|cudss|all]
 *                          [--device cpu|cuda|debug]
 *
 * THE CASE FLAGS COMPOSE. With none given the two defaults run; the first one
 * given clears them and each one after it adds. `--soloviev` and `--highbeta`
 * are opt-in because each changes what the table means -- one is linear, the
 * other bordered -- and neither should be averaged with the defaults.
 *
 * `--device debug` IS THE INSTRUMENT, NOT A SLOWER cuda. mfem::Device( "debug" )
 * has device memory semantics with host arithmetic and mprotect's the host page,
 * so a raw host read of a device-valid buffer is a NAMED FAULT WITH A BACKTRACE
 * rather than a wrong number. Every link of the alias chain behind M-79 was
 * found with it and none of them was visible under "cuda", which protects
 * nothing and hands back a stale host copy instead. Reach for it first.
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
#include "analytic/HighBetaPoloidal.hpp"
#include "analytic/Soloviev.hpp"

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

	struct Box { double minRadius, maxRadius, zMin, zMax; };

	Box box() { return { 0.6, 1.4, -0.6, 0.6 }; }

	mfem::Mesh makeMesh( Box const &b, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false,
			b.maxRadius - b.minRadius, b.zMax - b.zMin );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			double *c = mesh.GetVertex( v );
			c[ 0 ] += b.minRadius;
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

			double f( double radius, double z, double psi ) const override
			{
				return eq.f( radius, z, psi );
			}

			double dFdPsi( double radius, double z, double psi ) const override
			{
				return eq.dFdPsi( radius, z, psi );
			}

		private:
			Equilibrium const &eq;
	};

	/*
	 * THE SAME ADAPTER FOR A SOURCE WHOSE psi_ax IS AN UNKNOWN, which the
	 * high-beta arm needs and the other three do not. setSource( source, guess )
	 * closes the pair by a bordered Newton, so this must be a
	 * meq::NormalisedSource rather than a meq::Source.
	 *
	 * It holds the equilibrium BY VALUE where EquilibriumSource holds a
	 * reference: setNormalisation() mutates it once per residual evaluation, so
	 * a shared one would have two solves writing to the same psi_ax.
	 */
	template<typename Equilibrium>
	class NormalisedEquilibriumSource : public meq::NormalisedSource
	{
		public:
			explicit NormalisedEquilibriumSource( Equilibrium const &eqIn ) : eq( eqIn ) {}

			double f( double radius, double z, double psi ) const override
			{
				return eq.f( radius, z, psi );
			}

			double dFdPsi( double radius, double z, double psi ) const override
			{
				return eq.dFdPsi( radius, z, psi );
			}

			void setNormalisation( double psiAxis, double psiBoundary ) override
			{
				// The analytic fixtures are written for psi_bnd = 0. Refused
				// rather than ignored, which is this harness's own rule.
				if ( psiBoundary != 0.0 )
					throw std::invalid_argument(
						"NpcThreadScaling: the high-beta fixture cannot represent "
						"a non-zero boundary flux" );
				eq.setPsiAxis( psiAxis );
			}
			using meq::NormalisedSource::setNormalisation;

			double normalisation() const override { return eq.psiAxis(); }
			double boundaryNormalisation() const override { return 0.0; }

		private:
			Equilibrium eq;
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

		/*
		 * THE TWO FIELDS THAT MAKE A DEVICE RUN READABLE, and neither was here
		 * before upstream pointed out why they had to be.
		 *
		 * `finiteFailures` is mfem::Vector::CheckFinite() taken BEFORE any norm.
		 * Norml2() guards its reduction with fabs(v) > 0, so an ALL-NaN vector
		 * reports a norm of ZERO -- which would make a threaded-vs-serial
		 * comparison of two NaN fields read as perfect agreement, and would make
		 * a non-zero one meaningless in the other direction. It has to be asked
		 * separately and first.
		 *
		 * `exactError` is the L2 error against a closed form, and it exists
		 * because an ITERATION COUNT CANNOT DISCRIMINATE. Under a Device this
		 * harness reported 0/0 Newton iterations on a case needing four; zero
		 * iterations is also what a correct linear solve reports. Only a
		 * comparison against an answer known before the code runs separates
		 * them -- which is how upstream built their own reproduction, and the
		 * reason theirs could see what a two-arm comparison could not.
		 * Negative where the fixture has no closed form.
		 */
		int finiteFailures = -1;
		double exactError = -1.0;

		/// max |psi_h|, which separates "returned garbage" from "returned the
		/// INITIAL ITERATE" -- and the second is what a residual evaluation
		/// that silently produces nothing looks like from outside, since
		/// Newton then stops at iteration zero on an apparently converged
		/// residual and hands back what it was given.
		double psiPeak = -1.0;
	};

	std::vector<double> copyOf( mfem::GridFunction const &g )
	{
		// operator()( int ) IS A RAW HOST READ. Under mfem::Device( "debug" ) the
		// host page of a device-valid buffer is mprotect'd, so this faults rather
		// than returning a stale number -- which is the whole point of that
		// device and is how the alias chain behind M-79 was found. HostRead()
		// first. It is a no-op with no Device configured.
		g.HostRead();
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
	/*
	 * The L2 error against a closed form, where the fixture has one.
	 *
	 * AN OVERLOAD PAIR RATHER THAN A RUNTIME PREDICATE, and the difference is
	 * not style: a `if ( hasClosedForm( f ) )` inside one template still
	 * INSTANTIATES the psi() call for every fixture, and PressurePedestal has
	 * no psi() to instantiate. Overload resolution picks the concrete one for
	 * Solov'ev and the template for everything else, so only the reachable body
	 * is ever compiled.
	 */
	template<typename Fixture>
	double exactErrorOf( Solver &, Fixture const & )
	{
		return -1.0;
	}

	double exactErrorOf( Solver &solver,
	                     meq::analytic::SolovievEquilibrium const &fixture )
	{
		mfem::FunctionCoefficient exact( [ &fixture ]( mfem::Vector const &x )
		{
			return fixture.psi( x( 0 ), x( 1 ) );
		} );
		return solver.potential().ComputeL2Error( exact );
	}

	/// example5 has a closed form too, and the FAILING arm needs it as much as
	/// the passing one: "0/0 Newton iterations" says the solve stopped, and only
	/// this says what it stopped on.
	double exactErrorOf( Solver &solver,
	                     meq::analytic::ManufacturedNonlinear const &fixture )
	{
		mfem::FunctionCoefficient exact( [ &fixture ]( mfem::Vector const &x )
		{
			return fixture.psi( x( 0 ), x( 1 ) );
		} );
		return solver.potential().ComputeL2Error( exact );
	}

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
			// FINITENESS FIRST, before anything takes a norm. See Run.
			out.finiteFailures = solver.potential().CheckFinite()
			                   + solver.flux().CheckFinite();
			out.exactError = exactErrorOf( solver, fixture );
			out.psiPeak = solver.potential().Normlinf();
			out.psi = copyOf( solver.potential() );
			out.flux = copyOf( solver.flux() );
		}
		return out;
	}

	/*
	 * THE BORDERED ARM -- pass 3, --highbeta.
	 *
	 * psi_ax IS AN UNKNOWN HERE AND IS NOT ON THE OTHER THREE, which is what
	 * earns it a pass of its own rather than another row. The border spends
	 * extra residual evaluations and extra backsolves against the SAME
	 * factorisation, so the balance between element-local work and the trace
	 * solve is a different one -- and it is the configuration MEQ actually runs
	 * on a physical problem, where example5 and the pedestal hold psi_ax fixed.
	 *
	 * OPT-IN rather than default, for the reason --soloviev is: it changes the
	 * shape of the table, and a thread-scaling figure averaged over a bordered
	 * case and two unbordered ones describes neither.
	 *
	 * IT HAS NO CLOSED FORM, so `exactError` stays negative and the row prints
	 * max|psi_h| alone. That is weaker than the other three arms and is the
	 * reason not to read this one on its own under a device: an iteration count
	 * cannot discriminate a wrong answer from a right one. See Run.
	 */
	Run runBordered( Box const &b, int order, int n, int nu, double amplitude,
	                 AM mode, TS trace )
	{
		Run out;

		double const w = b.maxRadius - b.minRadius;
		double const h = b.zMax - b.zMin;
		double const eigenvalue = M_PI*M_PI*( 1.0/( w*w ) + 1.0/( h*h ) );
		double const estimate = std::sqrt( nu*amplitude/eigenvalue );

		mfem::Mesh mesh = makeMesh( b, n );
		Solver solver( mesh, order );
		solver.setAssemblyMode( mode );
		solver.setTraceSolver( trace );

		NormalisedEquilibriumSource<meq::analytic::HighBetaPoloidal> source(
			meq::analytic::HighBetaPoloidal::peaked( nu, amplitude, estimate ) );
		mfem::ConstantCoefficient zero( 0.0 );

		// The guess is PART OF THE PROBLEM STATEMENT here, not an optimisation:
		// at a fixed normalisation this equation has a small positive solution
		// and a large one, and Newton from the Dirichlet datum walks onto the
		// small branch. HighBetaConvergence records the same bump for the reason.
		double const minRadius = b.minRadius;
		double const zMin = b.zMin;
		mfem::FunctionCoefficient guess(
			[ estimate, minRadius, zMin, w, h ]( mfem::Vector const &x )
			{
				return estimate*std::sin( M_PI*( x( 0 ) - minRadius )/w )
				       *std::sin( M_PI*( x( 1 ) - zMin )/h );
			} );

		solver.setSource( source, estimate );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-10, 1.0e-14, 30 );

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
			// FINITENESS FIRST, before anything takes a norm. See Run.
			out.finiteFailures = solver.potential().CheckFinite()
			                   + solver.flux().CheckFinite();
			out.psiPeak = solver.potential().Normlinf();
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
	// OPT-IN, because it is a LINEAR problem and this harness is about the
	// nonlinear path: including it by default would put a one-step solve in
	// every thread-scaling table, where it measures assembly and nothing else.
	bool wantSoloviev = false;
	// OPT-IN for a different reason -- it is the BORDERED case, psi_ax an
	// unknown. See runBordered().
	bool wantHighBeta = false;
	/*
	 * THE CASE FLAGS COMPOSE, AND THEY USED NOT TO.
	 *
	 * Each was written as "turn the other one off", which reads correctly for
	 * one flag and silently selects NOTHING for two: `--example5 --pedestal`
	 * cleared both and the binary printed a header, a correctness block and no
	 * rows. Nothing said so, because running no cases is not an error.
	 *
	 * So the first case flag clears the defaults and every flag after it ADDS.
	 * One flag behaves exactly as before; two now mean what they say.
	 */
	bool selectionMade = false;
	// Which trace solvers to include. Restricting matters here in a way it does
	// not in TraceSolverScaling: UMFPACK and PARDISO respond to MKL_NUM_THREADS
	// in OPPOSITE directions, so a row that sweeps both at MKL > 1 is dominated
	// by UMFPACK's collapse and says nothing about PARDISO.
	std::string wantTrace = "all";
	/*
	 * The device, if one is asked for. A FLAG AND NOT A DEFAULT, for the reason
	 * TraceSolverScaling.cpp gives at the same place: mfem::Device is global
	 * state deciding where every Vector afterwards allocates, cuDSS cannot be
	 * measured without one -- it reads its matrix through
	 * SparseMatrix::ReadI/ReadJ/ReadData and its vectors through Read()/Write(),
	 * which hand back host pointers otherwise and abort inside CUDA -- and
	 * configuring one changes where the rest of meq allocates too.
	 *
	 * WHAT IS NEW HERE, AND WHY IT IS WORTH A SECOND HARNESS. TraceSolverScaling
	 * times cuDSS on an EXTRACTED trace matrix: assembly, reduction and one
	 * factorisation, with no MultNL() anywhere. That is group 4 of MFEM's
	 * doc/HDG-DEVICE-OFFLOAD.md measured on its own, which is the thing that
	 * plan says not to do -- and it cannot show the cost of not doing groups 2
	 * and 3, because a linear solve calls the element-local nonlinear path zero
	 * times. A whole NPC solve calls it once per residual and once per Jacobian,
	 * so this is where "the trace solve alone is worse than nothing" would
	 * become a number rather than a prediction.
	 *
	 * IT DOES NOT GET THAT FAR, AND WHAT IT FINDS INSTEAD IS WORTH MORE THAN THE
	 * TIMING WOULD HAVE BEEN. A whole NPC solve with an mfem::Device configured
	 * for CUDA does not compute the right answer, and cuDSS is a BYSTANDER
	 * rather than the cause -- measured on all three trace solvers at
	 * --device cuda, k=2, n=32:
	 *
	 *   OMP_NUM_THREADS=8   aborts, and from several threads at once, in
	 *                       MemoryManager::CheckHostMemoryType_ reached through
	 *                       Vector::AddElementVector inside
	 *                       DarcyHybridization::MultNL's own OpenMP region --
	 *                       "host pointer is not registered". MFEM_USE_EXCEPTIONS
	 *                       makes it a throw, which cannot leave a parallel
	 *                       region, so it lands as `terminate called recursively`.
	 *   OMP_NUM_THREADS=1   NO abort, and this is the dangerous one: example5
	 *                       reports 0/0 Newton iterations on a case that needs
	 *                       four, and the flux disagrees between assembly modes
	 *                       by 4.006e-02 -- the SAME value to four figures on
	 *                       UMFPack, PARDISO and cuDSS, which is what says there
	 *                       is one common fault rather than three.
	 *
	 * THE 0/0 IS REPAIRED AND THE ABORT IS NOT, AND THEY WERE NEVER ONE
	 * FAULT. M-79 found a device-memory alias chain at four sites, two of them
	 * MEQ's own, and with those fixed the OMP_NUM_THREADS=1 symptom above is
	 * gone: the same configuration now runs a real Newton to a residual floor
	 * of 2.14e-15 against the host's 7.03e-17, rather than reporting 0/0 and
	 * an identically zero potential.
	 *
	 * WHAT REMAINS, RE-MEASURED RATHER THAN REMEMBERED, IS TWO THINGS AND
	 * NEITHER IS UNIVERSAL. At OMP_NUM_THREADS=8 a device still aborts on the
	 * PLAIN path -- examples/mhd-rectangle.toml, the message above, from
	 * several threads at once -- while the BORDERED path survives it and
	 * reaches the host's answer, 12 iterations and psi_ax 9.484400e-02 on
	 * examples/limited-tokamak.toml. And at one thread the device degrades
	 * Newton's RATE by a case-dependent amount: 5 steps to 9 on the plain
	 * case, none at all on that bordered one. MEASUREMENTS.md M-129 is the
	 * table.
	 *
	 * ONE ROW OF THAT TABLE HAS SINCE BEEN WALKED AND WAS NOT THE DEVICE.
	 * examples/machine-f-diiid failed its bordered Newton outright and
	 * finished through the driver's ladder; that was two unsynced reads of
	 * MEQ's own -- prepare()'s alias seeding and GetEssentialTrueDofs() -- and
	 * with them fixed the case solves on a device in 2 steps to the host's
	 * every digit. MEASUREMENTS.md M-130. Read the remaining plain-path
	 * degradation as a suspect rather than a datum; nobody has walked it.
	 *
	 * So the refusal apps/meq.cpp carries rests on the trade argument, which
	 * is unchanged: the integrators are 46-53% of an NPC step and have no
	 * kernels. "The solve does not survive a Device at all" was the second
	 * reason it used to rest on and it is too strong -- the bordered path,
	 * which is what MEQ runs, survives one.
	 *
	 * Re-run this the day the integrators get device kernels. The flag and the
	 * cuDSS column are here so that it costs one command rather than an
	 * afternoon.
	 */
	std::string device = "cpu";

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--orders" && i + 1 < argc )
			orders = parseList( argv[ ++i ] );
		else if ( arg == "--sizes" && i + 1 < argc )
			sizes = parseList( argv[ ++i ] );
		else if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
		else if ( arg == "--pedestal" || arg == "--example5"
		          || arg == "--soloviev" || arg == "--highbeta" )
		{
			if ( !selectionMade )
			{
				selectionMade = true;
				wantExample5 = wantPedestal = false;
				wantSoloviev = wantHighBeta = false;
			}
			if ( arg == "--pedestal" )      wantPedestal = true;
			else if ( arg == "--example5" ) wantExample5 = true;
			else if ( arg == "--soloviev" ) wantSoloviev = true;
			else                            wantHighBeta = true;
		}
		else if ( arg == "--trace" && i + 1 < argc )
			wantTrace = argv[ ++i ];
		else if ( arg == "--device" && i + 1 < argc )
			device = argv[ ++i ];
	}

	/*
	 * CONSTRUCTED FIRST AND LEFT ALIVE for the whole run, which mfem::Device
	 * requires. Before the meshes, before the sources, before any Vector.
	 */
	std::unique_ptr<mfem::Device> deviceHandle;
	if ( device != "cpu" )
	{
#ifndef MFEM_USE_CUDA
		// Exit 0, not 1: asking for a device in a build that has none is a skip,
		// not a failure. This binary is not a ctest, but the rule is the same
		// one TraceSolverScaling follows and there is no reason to differ.
		std::printf( "\n  device \"%s\" requested, but this MFEM has no CUDA -- "
		             "skipped\n\n", device.c_str() );
		return 0;
#else
		deviceHandle = std::make_unique<mfem::Device>( device.c_str() );
		std::printf( "\n" );
		deviceHandle->Print();
#endif
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
	for ( TS t : { TS::UMFPack, TS::Pardiso, TS::cuDSS } )
	{
		if ( !Solver::traceSolverAvailable( t ) )
			continue;
		if ( wantTrace == "umfpack" && t != TS::UMFPack )
			continue;
		if ( wantTrace == "pardiso" && t != TS::Pardiso )
			continue;
		if ( wantTrace == "cudss" && t != TS::cuDSS )
			continue;
		// cuDSS IS OPT-IN EVEN WHEN THE BUILD HAS IT, and it is the only one of
		// the three that is. Without an mfem::Device it does not fall back --
		// it reads host pointers through the device-aware accessors and aborts
		// inside CUDA with a message naming cudaMemcpyDeviceToDevice and
		// nothing about the solver -- so an "all" sweep on a CPU run must not
		// pick it up. --trace cudss --device cuda is the way to reach it.
		if ( t == TS::cuDSS && ( wantTrace != "cudss" || device == "cpu" ) )
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
	/*
	 * THE LINEAR-SOURCE ARM -- pass 2, --soloviev.
	 *
	 * The experiment HDG-DEVICE-AND-LEVEL-2-FROM-HDGDEV.md asks for: the same
	 * solve with the source "replaced by a linear one", to separate a genuinely
	 * non-linear integrand from the other three differences between MEQ's
	 * failing device run and upstream's passing reproduction.
	 *
	 * SolovievEquilibrium IS THE RIGHT SUBSTITUTION AND A TRULY LINEAR SOURCE
	 * WOULD NOT BE, which is the whole subtlety of the request. `dFdPsi` is
	 * identically zero, so the problem is linear and Newton takes one exact
	 * step -- but it still arrives through setSource( Source const & ), so
	 * `nonlinearSource` is non-null, `usesNonlinearForms()` is true, and the
	 * FORM ROUTING IS UNCHANGED: meq::SourceIntegrator on the potential-mass
	 * non-linear form's domain, both HDGDiffusionIntegrators on its faces, M_p
	 * null, and therefore the same c_bfi_p branch of EnableHybridization().
	 * `AssembleElementVector` still runs per element per residual.
	 *
	 * Handing MEQ a psi-independent source instead would take the OTHER branch
	 * of buildForms(), construct M_p, and change which arm of
	 * EnableHybridization() fires -- so a pass would implicate the routing as
	 * readily as the integrand and the experiment would decide nothing. This is
	 * the same construction upstream used for their own reproduction, where an
	 * inert domain integrator sits on the non-linear form so that "the routing
	 * is yours".
	 *
	 * And it has a CLOSED FORM, which the other two fixtures do not. That is
	 * what makes the arm readable at all: under a Device this harness reported
	 * 0/0 Newton iterations, and zero iterations is also what a correct linear
	 * solve reports.
	 */
	meq::analytic::SolovievEquilibrium const soloviev
		= meq::analytic::SolovievEquilibrium::nstx();

	std::printf( "\n  a whole nonlinear solve, serial assembly against threaded\n" );
	std::printf( "    %-10s %2s %5s %8s %9s %9s %8s %9s %9s %8s %5s\n",
	             "case", "k", "n", "solver",
	             "prep_ser", "prep_thr", "prep_x",
	             "solve_ser", "solve_thr", "solve_x", "its" );
	// Flushed, because a row at MKL > 1 with UMFPACK can take minutes and an
	// unflushed header makes that look like a hang rather than a measurement.
	std::fflush( stdout );

	for ( int pass = 0; pass < 4; ++pass )
	{
		char const *name = ( pass == 0 ) ? "example5"
		                 : ( pass == 1 ) ? "pedestal"
		                 : ( pass == 2 ) ? "soloviev" : "highbeta";
		if ( pass == 0 && !wantExample5 ) continue;
		if ( pass == 1 && !wantPedestal ) continue;
		if ( pass == 2 && !wantSoloviev ) continue;
		if ( pass == 3 && !wantHighBeta ) continue;

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
		// Solov'ev's datum is the trace of ITS exact solution, as example5's is
		// of its own. The source is F and not F/R; the solver applies the 1/R.
		mfem::FunctionCoefficient solovievDatum(
			[ &soloviev ]( mfem::Vector const &x )
			{
				return soloviev.psi( x( 0 ), x( 1 ) );
			} );
		// Pass 3 binds solovievDatum and never reads it: runBordered() carries
		// its own homogeneous datum and its own guess, both being part of that
		// problem's statement rather than of this harness's.
		mfem::Coefficient &datum = ( pass == 0 )
			? static_cast<mfem::Coefficient &>( exactDatum )
			: ( pass == 1 ) ? static_cast<mfem::Coefficient &>( rampDatum )
			: static_cast<mfem::Coefficient &>( solovievDatum );
		// Only the pedestal needs one: its source vanishes at psi = 0, so
		// homogeneous data would make psi == 0 solve the problem.
		bool const guessFromDatum = ( pass == 1 );

		for ( int order : orders )
		{
			for ( int n : sizes )
			{
				/*
				 * ONE DISPATCH FOR EVERY CALL SITE, AND THE DUPLICATION IT
				 * REPLACES CARRIED A BUG.
				 *
				 * The four sites below -- serial, threaded, and the two the
				 * trace-solver comparison runs -- each spelled the pass out as
				 * a ternary chain. The two in the comparison stopped at
				 * `pass == 0 ? example5 : pedestal`, so on the --soloviev arm
				 * it ran the PEDESTAL fixture carrying Solov'ev's datum: a
				 * combination that is neither case, and it still printed an
				 * agreement figure. Adding a fourth pass to four chains would
				 * have been a fifth chance to do it again.
				 */
				auto runCase = [ & ]( AM mode, TS trace ) -> Run
				{
					switch ( pass )
					{
						case 0:  return runOnce( box(), order, n, example5, datum, guessFromDatum, mode, trace );
						case 1:  return runOnce( box(), order, n, pedestal, datum, guessFromDatum, mode, trace );
						case 2:  return runOnce( box(), order, n, soloviev, datum, guessFromDatum, mode, trace );
						// nu = 2, amplitude 1: HighBetaConvergence's own
						// converging point, so a non-convergence here is the
						// harness rather than the physics.
						default: return runBordered( box(), order, n, 2, 1.0, mode, trace );
					}
				};

				for ( TS trace : traceSolvers )
				{
					Run best_s, best_t;
					best_s.prepareTime = best_s.solveTime = 1e30;
					best_t.prepareTime = best_t.solveTime = 1e30;

					for ( int repeat = 0; repeat < repeats; ++repeat )
					{
						Run s = runCase( AM::Serial, trace );
						if ( s.prepareTime < best_s.prepareTime ) best_s.prepareTime = s.prepareTime;
						if ( s.solveTime < best_s.solveTime ) best_s.solveTime = s.solveTime;
						best_s.newtonIterations = s.newtonIterations;
						best_s.converged = s.converged;
						best_s.finiteFailures = s.finiteFailures;
						best_s.exactError = s.exactError;
						best_s.psiPeak = s.psiPeak;
						if ( best_s.psi.empty() ) { best_s.psi = s.psi; best_s.flux = s.flux; }
					}

					if ( canThread )
					{
						for ( int repeat = 0; repeat < repeats; ++repeat )
						{
							Run t = runCase( AM::Threaded, trace );
							if ( t.prepareTime < best_t.prepareTime ) best_t.prepareTime = t.prepareTime;
							if ( t.solveTime < best_t.solveTime ) best_t.solveTime = t.solveTime;
							best_t.newtonIterations = t.newtonIterations;
							best_t.converged = t.converged;
							best_t.finiteFailures = t.finiteFailures;
							best_t.exactError = t.exactError;
							best_t.psiPeak = t.psiPeak;
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

					// NON-FINITE FIRST. A NaN field reports an L2 error and a
					// worst-difference of whatever the guarded reductions make
					// of it, so this has to be said before either number is
					// read rather than inferred from them afterwards.
					if ( best_s.finiteFailures > 0
					     || ( canThread && best_t.finiteFailures > 0 ) )
						std::printf( "  *** NON-FINITE ENTRIES: %d serial, %d threaded",
						             best_s.finiteFailures,
						             canThread ? best_t.finiteFailures : 0 );
					else if ( best_s.exactError >= 0.0 )
						std::printf( "  L2 vs exact: %.3e ser", best_s.exactError );
					if ( best_s.psiPeak >= 0.0 )
						std::printf( "  max|psi_h| %.6e", best_s.psiPeak );
					if ( canThread && best_s.exactError >= 0.0
					     && best_s.finiteFailures == 0 )
						std::printf( ", %.3e thr", best_t.exactError );
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
					Run reference = runCase( AM::Serial, traceSolvers[ 0 ] );
					for ( size_t i = 1; i < traceSolvers.size(); ++i )
					{
						Run other = runCase( AM::Serial, traceSolvers[ i ] );
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

	/*
	 * The exit code. Only correctness, never a timing.
	 *
	 * AND EXACTNESS IS ONLY CLAIMED AT MKL_NUM_THREADS=1, WHICH THE DIAGNOSTIC
	 * HAS TO SAY OR IT SENDS THE READER HUNTING A RACE THAT IS NOT THERE.
	 * The two assembly modes hand MKL DIFFERENT THREAD COUNTS at MKL > 1: the
	 * serial element loop is outside any parallel region and its element-local
	 * dgetrs and dgemm get all of them, while the threaded loop is an active
	 * OpenMP region and MKL suppresses its own threading inside one, so the same
	 * kernels run sequentially. A blocked BLAS-3 sums in a different order from
	 * an unblocked loop, so the two modes are then computing the same thing by
	 * different associations -- arithmetic reassociation inside MKL, not a race
	 * in MFEM and not shared scratch in MEQ. Measured here: 1.4e-15 in psi and
	 * 1.8e-13 in the flux at MKL=8, against 0.0e+00 in both at MKL=1.
	 *
	 * GradShafranov.cpp's setAssemblyMode() documentation records the same
	 * mechanism and the same magnitudes, which is why the message points there
	 * rather than repeating the argument.
	 */
	bool const mklMayReassociate = ( mklRequested > 1 );
	char const *reassociationNote = mklMayReassociate
		? "      BUT MKL_NUM_THREADS > 1, AND EXACTNESS IS ONLY CLAIMED AT 1.\n"
		  "      The serial element loop gets every MKL thread and the threaded\n"
		  "      one gets none -- MKL suppresses itself inside an active OpenMP\n"
		  "      region -- so a blocked BLAS-3 reassociates against an unblocked\n"
		  "      loop. That is arithmetic, not a race. Re-run at\n"
		  "      MKL_NUM_THREADS=1 before reading anything into this.\n"
		: "      MFEM documents exactness on both its threaded loops, so\n"
		  "      this is a real change rather than a tolerance to widen.\n"
		  "      On a NONLINEAR source the loop is MultNL's, and the\n"
		  "      caller's own integrators sit on it: check whether\n"
		  "      anything MEQ installs holds per-point scratch as a\n"
		  "      member. meq::SourceIntegrator's `shape` is guarded on\n"
		  "      MFEM_THREAD_SAFE for exactly that reason.\n";

	if ( canThread && worstThreadedPsi != 0.0 )
	{
		std::printf( "\n  *** threaded assembly is NOT bit for bit in psi (%.3e).\n%s",
		             worstThreadedPsi, reassociationNote );
		// Not a failure at MKL > 1: the two modes are entitled to differ there,
		// and counting it would make the binary's exit code a statement about
		// the environment rather than about the code.
		if ( !mklMayReassociate )
			++failures;
	}
	if ( canThread && worstThreadedFlux != 0.0 )
	{
		std::printf( "\n  *** threaded assembly is NOT bit for bit in the flux (%.3e).\n%s",
		             worstThreadedFlux, mklMayReassociate ? reassociationNote : "" );
		if ( !mklMayReassociate )
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
