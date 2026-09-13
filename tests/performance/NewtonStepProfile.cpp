/*
 * WHERE A NEWTON STEP'S TIME GOES: the four-leg split of
 * HDG-NEWTON-STEP-PROFILE-FROM-HDGDEV.md, plus that request's level-3 shape
 * parameters.
 *
 * NOT A CTEST, per the standing rule for tests/performance/: every number below
 * is a timing. It exits non-zero only if a case fails to converge, because a
 * leg share taken from a diverged solve is a share of the wrong thing.
 *
 * WHAT IT IS FOR. Upstream has two candidate pieces of work in
 * DarcyHybridization, about 4-5x apart in cost, and which one is worth building
 * turns on where a MEQ Newton step actually spends its time:
 *
 *   (a) cache the state-independent half of the local condensation -- A^-1 B^T,
 *       B A^-1 B^T, A^-1 C^T, B A^-1 C^T, C A^-1 B^T, all constant across the
 *       whole Newton loop under LocalOpType::PotNL and rebuilt every step. By
 *       flop count ~70-75% of ComputeElementH. One to two sessions.
 *   (b) CCSZ interpolatory HDG. Five stages, five to eight sessions, and its
 *       performance payoff is bounded by the residual integrators plus
 *       ConstructGrad.
 *
 * Their own profile and gffp's disagree by enough that the answer flips: on
 * 128x128 quads at order 2 they read the trace solve at 54-59% of a whole solve,
 * where gffp reads 152 ms a step against a 1.4 ms solve. **If the trace solve is
 * over half, neither (a) nor (b) is the next thing to build**, and that answer
 * costs an afternoon instead of two weeks. This harness exists to give it.
 *
 * THE LEGS ARE SPLIT BY CALL SITE AND THEY DO NOT OVERLAP. See
 * GradShafranovSolver::StepProfile for where each timer sits and why. The
 * short form:
 *
 *   residual   DarcyNPCOperator::Mult
 *   gradient   DarcyNPCOperator::GetGradient -- "assemble and factor the
 *              Jacobian at x" in MFEM's own words, so ComputeH() is IN HERE
 *   factor     the trace package's SetOperator: the numeric factorisation
 *   backsolve  the trace package's Mult
 *   component  refreshPlasmaComponent(), the XP-1 flood fill, zero on these
 *              fixtures and reported so that its absence is on the record
 *   other      the REMAINDER, total minus the above. Source refresh, the
 *              critical-point search, the boundary datum, and MFEM's own vector
 *              arithmetic inside NewtonSolver
 *
 * `other` being a remainder rather than a measurement is deliberate and it is
 * the honest form: anything the five timed legs miss lands there instead of
 * silently inflating one of them.
 *
 * THREE THINGS THE REQUEST NAMES AS SPOILERS, AND WHAT IS DONE ABOUT EACH.
 *
 *   "Separate allocation from arithmetic." The timing wrappers allocate nothing
 *   and copy no vector, so no leg carries a malloc that is not the wrapped
 *   call's own. What each leg DOES carry is written above. One thing is worth
 *   naming: prepare() is OUTSIDE solve() and therefore outside `total`
 *   entirely, so the form assembly and the essential-BC elimination are in
 *   neither the legs nor the remainder -- they are reported separately as
 *   `prep`, since they are paid once per mesh rather than once per step.
 *
 *   "Pin MKL_NUM_THREADS." newton-step-profile.sh does, and this binary prints
 *   what it was given rather than what it wanted.
 *
 *   "Medians of five, not one run." --repeats defaults to 5 and the statistic
 *   is the MEDIAN, not the best: a share is not a share if it comes from one
 *   run, and for a SHARE the median is the right centre where for a minimum
 *   wall time the best-of is. Both are printed, so the spread is visible.
 *
 * Usage:  NewtonStepProfile [--repeats N] [--serial] [--threaded]
 *                           [--case pedestal|example5|all]
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

namespace
{
	using Solver = meq::GradShafranovSolver;
	using AM = Solver::AssemblyMode;
	using Profile = Solver::StepProfile;

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
	/// NpcThreadScaling carries, duplicated for the same reason: the shared one
	/// pulls in Boost.Test and these binaries are deliberately not linked
	/// against it.
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

	/*
	 * The BORDERED path's source: psi_ax is an unknown, so the fixture is held
	 * BY VALUE and setNormalisation() mutates it. The same adapter
	 * tests/convergence/ConvergenceHarness.hpp carries, duplicated for the same
	 * reason as EquilibriumSource above -- that header includes Boost.Test.
	 */
	template<typename Equilibrium>
	class NormalisedEquilibriumSource : public meq::NormalisedSource
	{
		public:
			explicit NormalisedEquilibriumSource( Equilibrium const &eqIn ) : eq( eqIn ) {}

			double f( double r, double z, double psi ) const override
			{
				return eq.f( r, z, psi );
			}

			double dFdPsi( double r, double z, double psi ) const override
			{
				return eq.dFdPsi( r, z, psi );
			}

			void setNormalisation( double psiAxis, double psiBoundary ) override
			{
				// The analytic fixtures are written for psi_bnd = 0. Refused
				// rather than ignored, which is the harness's own rule.
				if ( psiBoundary != 0.0 )
					throw std::invalid_argument(
						"NewtonStepProfile: the high-beta fixture cannot represent "
						"a non-zero boundary flux" );
				eq.setPsiAxis( psiAxis );
			}
			using meq::NormalisedSource::setNormalisation;

			double normalisation() const override { return eq.psiAxis(); }
			double boundaryNormalisation() const override { return 0.0; }

		private:
			Equilibrium eq;
	};

	/// Everything one run reports: the legs, prep, and the level-3 shape.
	struct Run
	{
		char const *name = "";
		int n = 0;
		bool bordered = false;

		Profile profile;
		double prepSeconds = 0.0;
		int newtonIterations = 0;
		long numericFactorisations = 0;
		long symbolicFactorisations = 0;
		bool converged = false;

		// Level 3, which is shape rather than timing and so is taken once.
		int elements = 0;
		int order = 0;
		int fluxDofsPerElement = 0;
		int potentialDofsPerElement = 0;
		int traceDofsPerFace = 0;
		int facesPerElement = 0;
		int traceDofs = 0;
		int essentialTraceDofs = 0;
	};

	/// Level 3, taken once per run because it is shape rather than timing.
	void captureShape( Run &out, Solver &solver, mfem::Mesh &mesh, int order )
	{
		out.elements = mesh.GetNE();
		out.order = order;
		out.fluxDofsPerElement = solver.fluxSpace().GetFE( 0 )->GetDof()
		                       * solver.fluxSpace().GetVDim();
		out.potentialDofsPerElement = solver.potentialSpace().GetFE( 0 )->GetDof();
		out.traceDofsPerFace = solver.traceSpace().GetFaceElement( 0 )->GetDof();
		// In two dimensions an element's faces are its edges, which is what
		// DarcyHybridization's nf counts.
		out.facesPerElement = mesh.GetElement( 0 )->GetNEdges();
		out.traceDofs = solver.traceSpace().GetVSize();
		out.essentialTraceDofs = solver.essentialTraceDofs().Size();
	}

	template<typename Fixture>
	Run runOnce( Box const &b, int order, int n, Fixture const &fixture,
	             mfem::Coefficient &datum, bool guessFromDatum, AM mode )
	{
		Run out;
		mfem::Mesh mesh = makeMesh( b, n );
		Solver solver( mesh, order );
		solver.setAssemblyMode( mode );

		EquilibriumSource<Fixture> const source( fixture );
		solver.setSource( source );
		solver.setBoundaryData( datum );
		if ( guessFromDatum )
			solver.setInitialGuess( datum );

		// prepare() IS OUTSIDE THE LEGS AND IS TIMED SEPARATELY. It is the form
		// assembly and the essential-BC elimination, paid once per mesh, so
		// folding it into a per-step share would inflate whichever leg it landed
		// in -- and it would land in the remainder, which is where an unexplained
		// cost is hardest to notice.
		double const t0 = now();
		solver.prepare();
		out.prepSeconds = now() - t0;

		try
		{
			solver.solve();
			out.converged = true;
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}

		out.profile = solver.stepProfile();
		out.newtonIterations = solver.newtonIterations();
		out.numericFactorisations = solver.numericFactorisations();
		out.symbolicFactorisations = solver.symbolicFactorisations();

		captureShape( out, solver, mesh, order );
		return out;
	}

	/*
	 * THE BORDERED PATH -- psi_ax as an unknown, which is what MEQ's free-boundary
	 * and high-beta work actually runs and the one upstream asks about by name
	 * ("whether psi_ax's refresh forces any re-assembly we have not spotted").
	 *
	 * It is a separate function rather than a flag on runOnce because it is a
	 * different solve: setSource( NormalisedSource &, double ) closes the system
	 * by a bordered Newton with a hand-rolled loop, so the legs are timed at
	 * different call sites inside the library -- and the border spends EXTRA
	 * residual evaluations and EXTRA backsolves against the same factorisation,
	 * which is precisely what the call counts below are for.
	 */
	Run runBordered( Box const &b, int order, int n, int nu, double amplitude,
	                 AM mode )
	{
		Run out;
		out.bordered = true;

		double const w = b.rMax - b.rMin;
		double const h = b.zMax - b.zMin;
		double const eigenvalue = M_PI*M_PI*( 1.0/( w*w ) + 1.0/( h*h ) );
		double const estimate = std::sqrt( nu*amplitude/eigenvalue );

		mfem::Mesh mesh = makeMesh( b, n );
		Solver solver( mesh, order );
		solver.setAssemblyMode( mode );

		NormalisedEquilibriumSource<meq::analytic::HighBetaPoloidal> source(
			meq::analytic::HighBetaPoloidal::peaked( nu, amplitude, estimate ) );
		mfem::ConstantCoefficient zero( 0.0 );

		// The guess is PART OF THE PROBLEM STATEMENT here, not an optimisation:
		// at a fixed normalisation this equation has a small positive solution
		// and a large one, and Newton from the Dirichlet datum walks onto the
		// small branch. HighBetaConvergence records the same bump for the reason.
		double const rMin = b.rMin;
		double const zMin = b.zMin;
		mfem::FunctionCoefficient guess(
			[ estimate, rMin, zMin, w, h ]( mfem::Vector const &x )
			{
				return estimate*std::sin( M_PI*( x( 0 ) - rMin )/w )
				       *std::sin( M_PI*( x( 1 ) - zMin )/h );
			} );

		solver.setSource( source, estimate );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-10, 1.0e-14, 30 );

		double const t0 = now();
		solver.prepare();
		out.prepSeconds = now() - t0;

		try
		{
			solver.solve();
			out.converged = true;
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}

		out.profile = solver.stepProfile();
		out.newtonIterations = solver.newtonIterations();
		out.numericFactorisations = solver.numericFactorisations();
		out.symbolicFactorisations = solver.symbolicFactorisations();
		captureShape( out, solver, mesh, order );
		return out;
	}

	double median( std::vector<double> v )
	{
		if ( v.empty() )
			return 0.0;
		std::sort( v.begin(), v.end() );
		std::size_t const m = v.size()/2;
		return ( v.size() % 2 == 1 )
		       ? v[ m ] : 0.5*( v[ m - 1 ] + v[ m ] );
	}

	double smallest( std::vector<double> const &v )
	{
		return v.empty() ? 0.0 : *std::min_element( v.begin(), v.end() );
	}

	char const *modeName( AM m ) { return ( m == AM::Threaded ) ? "threaded" : "serial"; }
}

int main( int argc, char **argv )
{
	int repeats = 5;
	bool wantSerial = true;
	bool wantThreaded = true;
	std::string wantCase = "all";

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
		else if ( arg == "--serial" )
			wantThreaded = false;
		else if ( arg == "--threaded" )
			wantSerial = false;
		else if ( arg == "--case" && i + 1 < argc )
			wantCase = argv[ ++i ];
	}

	int const ompThreads =
#ifdef MFEM_USE_OPENMP
		omp_get_max_threads();
#else
		1;
#endif

	std::printf( "\n=== MEQ Newton-step leg profile ===\n" );
	std::printf( "  OMP_NUM_THREADS=%-5s MKL_NUM_THREADS=%-5s "
	             "omp_get_max_threads()=%d, MEDIAN of %d\n",
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
#ifdef MFEM_USE_MKL_PARDISO
	std::printf( " MKL_PARDISO" );
#endif
#ifdef MFEM_USE_SUITESPARSE
	std::printf( " SUITESPARSE" );
#endif
	std::printf( "\n" );

	{
		// Which package the trace leg is, printed rather than assumed: the
		// default moved to PARDISO where the build has it, and the two do not
		// respond to MKL_NUM_THREADS in the same direction.
		mfem::Mesh probe = makeMesh( box(), 2 );
		Solver const probeSolver( probe, 1 );
		char const *which =
			( probeSolver.traceSolver() == Solver::TraceSolver::Pardiso ) ? "PARDISO"
			: ( probeSolver.traceSolver() == Solver::TraceSolver::cuDSS ) ? "cuDSS"
			: "UMFPack";
		std::printf( "  trace solver: %s (the library default), symbolic reuse ON\n",
		             which );
	}

	meq::analytic::ManufacturedNonlinear const example5
		= meq::analytic::ManufacturedNonlinear::example5();
	meq::analytic::PressurePedestal const pedestal
		= meq::analytic::PressurePedestal::pedestal();

	int failures = 0;
	std::vector<Run> shapes;

	std::printf( "\n  LEVEL 1 -- the legs, as a share of one whole solve()\n" );
	std::printf( "    %-9s %-9s %2s %4s %6s %8s %8s %8s %8s %8s %8s\n",
	             "case", "assembly", "k", "n", "its",
	             "total_s", "residual", "gradient", "factor", "backsolv", "other" );
	std::fflush( stdout );

	enum class Kind { Example5, Pedestal, HighBeta };
	struct Case
	{
		char const *name;
		int order;
		int n;
		Kind kind;
	};
	/*
	 * The cheap one and the heavy one the request asks for, PLUS a bordered one
	 * it does not ask for and should have.
	 *
	 * example5 at k = 2, n = 24 is exactly what was named. The pedestal at
	 * k = 2, n = 32 is GS-2 section 4.2's internal transport barrier, the
	 * stiffest source in MEQ's suite that still converges undamped, and it
	 * stands in for "PedestalConvergence" -- that binary sweeps four orders and
	 * three meshes and its total would be a mixture rather than a profile.
	 *
	 * highbeta at k = 2, n = 16 is psi_ax AS AN UNKNOWN. It is here because
	 * without it every row would be the plain path, and MEQ's headline
	 * configuration is not the plain path: the free-boundary and high-beta work
	 * closes the system by a bordered Newton, which spends extra residual
	 * evaluations and extra backsolves against the SAME factorisation. A profile
	 * that omitted it would answer for the fixtures rather than for MEQ.
	 */
	std::vector<Case> const cases = {
		{ "example5", 2, 24, Kind::Example5 },
		{ "pedestal", 2, 32, Kind::Pedestal },
		{ "highbeta", 2, 16, Kind::HighBeta },
	};

	for ( Case const &c : cases )
	{
		if ( wantCase != "all" && wantCase != c.name )
			continue;

		double const zMax = box().zMax;
		mfem::FunctionCoefficient exactDatum( [ &example5 ]( mfem::Vector const &x )
		{
			return example5.psi( x( 0 ), x( 1 ) );
		} );
		mfem::FunctionCoefficient rampDatum( [ zMax ]( mfem::Vector const &x )
		{
			return 0.3*x( 1 )/zMax;
		} );
		mfem::Coefficient &datum = ( c.kind == Kind::Pedestal )
			? static_cast<mfem::Coefficient &>( rampDatum )
			: static_cast<mfem::Coefficient &>( exactDatum );

		for ( AM mode : { AM::Serial, AM::Threaded } )
		{
			if ( mode == AM::Serial && !wantSerial ) continue;
			if ( mode == AM::Threaded && !wantThreaded ) continue;
			if ( mode == AM::Threaded && !Solver::assemblyModeAvailable( AM::Threaded ) )
				continue;

			std::vector<double> total, residual, gradient, factor, backsolve,
			                    component, computeH, other, prep;
			Run last;

			for ( int r = 0; r < repeats; ++r )
			{
				Run run;
				switch ( c.kind )
				{
					case Kind::Pedestal:
						run = runOnce( box(), c.order, c.n, pedestal, datum, true, mode );
						break;
					case Kind::Example5:
						run = runOnce( box(), c.order, c.n, example5, datum, false, mode );
						break;
					case Kind::HighBeta:
						// nu = 2, amplitude 1: HighBetaConvergence's own converging
						// point, so a non-convergence here would be the harness
						// rather than the physics.
						run = runBordered( box(), c.order, c.n, 2, 1.0, mode );
						break;
				}
				run.name = c.name;
				run.n = c.n;
				if ( !run.converged )
				{
					std::printf( "    %-9s %-9s %2d %4d   DID NOT CONVERGE\n",
					             c.name, modeName( mode ), c.order, c.n );
					++failures;
					break;
				}
				total.push_back( run.profile.totalSeconds );
				residual.push_back( run.profile.residualSeconds );
				gradient.push_back( run.profile.gradientSeconds );
				factor.push_back( run.profile.traceFactorSeconds );
				backsolve.push_back( run.profile.traceSolveSeconds );
				component.push_back( run.profile.componentSeconds );
				computeH.push_back( run.profile.computeHSeconds );
				other.push_back( run.profile.otherSeconds() );
				prep.push_back( run.prepSeconds );
				last = run;
			}
			if ( total.empty() )
				continue;

			double const t = median( total );
			auto share = [ t ]( double v ) { return ( t > 0.0 ) ? 100.0*v/t : 0.0; };

			std::printf( "    %-9s %-9s %2d %4d %6d %8.4f %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
			             c.name, modeName( mode ), c.order, c.n,
			             last.newtonIterations, t,
			             share( median( residual ) ), share( median( gradient ) ),
			             share( median( factor ) ), share( median( backsolve ) ),
			             share( median( other ) ) );
			std::printf( "    %-9s %-9s %2s %4s %6s %8s %8.4f %8.4f %8.4f %8.4f %8.4f   seconds\n",
			             "", "", "", "", "", "",
			             median( residual ), median( gradient ), median( factor ),
			             median( backsolve ), median( other ) );

			// Per step, which is what the request wants to fall out, plus the
			// call counts that say what "per step" means for each leg.
			int const its = std::max( 1, last.newtonIterations );
			std::printf( "    %-9s %-9s per step (%d): residual %.4f x%ld, "
			             "gradient %.4f x%ld, factor %.4f x%ld, backsolve %.4f x%ld\n",
			             "", "", its,
			             median( residual )/its, last.profile.residualCalls,
			             median( gradient )/its, last.profile.gradientCalls,
			             median( factor )/its, last.profile.traceFactorCalls,
			             median( backsolve )/its, last.profile.traceSolveCalls );
			std::printf( "    %-9s %-9s prepare() %.4f (outside total), "
			             "component fill %.4f x%ld, spread total %.4f-%.4f\n",
			             "", "",
			             median( prep ), median( component ),
			             last.profile.componentCalls,
			             smallest( total ), *std::max_element( total.begin(), total.end() ) );

			/*
			 * LEVEL 2 -- ComputeH() INSIDE the gradient leg, from upstream's
			 * static accumulator. Printed as a share of `gradient` AND of the
			 * whole solve, because those answer different questions: the first
			 * says how much of assembling the Jacobian is the element-local
			 * condensation, the second bounds what caching its state-independent
			 * half could ever be worth. It is NOT added to the leg shares above;
			 * it is already inside one of them.
			 */
			double const g = median( gradient );
			double const ch = median( computeH );
			std::printf( "    %-9s %-9s LEVEL 2: ComputeH %.4f x%ld -- %.1f%% of "
			             "the gradient leg, %.1f%% of the whole solve\n\n",
			             "", "", ch, last.profile.computeHCalls,
			             ( g > 0.0 ) ? 100.0*ch/g : 0.0, share( ch ) );
			std::fflush( stdout );

			if ( mode == AM::Threaded || !wantThreaded )
				shapes.push_back( last );
		}
	}

	/*
	 * LEVEL 3 -- shape, not timing. The request says they can compute the
	 * cacheable fraction themselves from the dimensions but cannot guess them.
	 */
	std::printf( "  LEVEL 3 -- the dimensions, so the cacheable fraction can be computed\n" );
	std::printf( "    %-9s %2s %4s %8s %6s %6s %6s %4s %9s %9s %7s\n",
	             "case", "k", "n", "elements", "na", "nd", "nc", "nf",
	             "trace_dof", "ess_trace", "numfact" );
	for ( Run const &shape : shapes )
		std::printf( "    %-9s %2d %4d %8d %6d %6d %6d %4d %9d %9d %4ld/%d\n",
		             shape.name, shape.order, shape.n, shape.elements,
		             shape.fluxDofsPerElement, shape.potentialDofsPerElement,
		             shape.traceDofsPerFace, shape.facesPerElement,
		             shape.traceDofs, shape.essentialTraceDofs,
		             shape.numericFactorisations, shape.newtonIterations );

	std::printf( "\n    na = flux dofs per element (vector), nd = potential dofs per element,\n"
	             "    nc = trace dofs per face, nf = faces per element. numfact is\n"
	             "    GetNumNumericFactorizations() against the Newton step count, which is\n"
	             "    the ratio the symbolic reuse is for -- one analysis, one factorisation\n"
	             "    per step.\n" );

	/*
	 * The two level-3 questions that are about the CONFIGURATION rather than the
	 * dimensions, answered here because the answer is a property of these
	 * fixtures and a reader should not have to infer it from their absence.
	 */
	std::printf( "\n    THE EXTENSION PATH IS OFF ON ALL THREE ROWS, and that does NOT\n"
	             "    limit what this table says about (a). These are fitted meshes with\n"
	             "    the Dirichlet datum on the mesh boundary -- no [boundary.shape], no\n"
	             "    exterior coupling. The tempting inference is that MEQ's CURVED and\n"
	             "    FREE-BOUNDARY cases are the other side of a line, because\n"
	             "    HDGExtensionIntegrator makes A asymmetric there. IT IS WRONG, and\n"
	             "    asymmetry read as state dependence is the error. Bnl_data is written\n"
	             "    at ONE site, ConstructGrad(), and only from grad_Aup -- the (0,1)\n"
	             "    block of a BLOCK NON-LINEAR integrator's element gradient, i.e. a\n"
	             "    flux law depending on the potential. HDGExtensionIntegrator is a\n"
	             "    BilinearFormIntegrator on the flux mass BILINEAR form, so it lands\n"
	             "    in A and is assembled once per mesh. A^-1 B^T is constant whether or\n"
	             "    not A is symmetric, so (a)'s precondition survives the extension\n"
	             "    path and these shares carry over to the problems MEQ exists for.\n"
	             "\n    psi_ax IS an unknown on the highbeta row and is not on the other\n"
	             "    two, which is what makes that row worth reading separately. Its\n"
	             "    border costs EXTRA RESIDUAL EVALUATIONS and EXTRA BACKSOLVES against\n"
	             "    the same factorisation -- visible in the per-step call counts above,\n"
	             "    where residual and backsolve calls exceed the iteration count while\n"
	             "    factor calls do not. It forces no re-assembly: setNormalisation()\n"
	             "    changes what the source integrator EVALUATES, and the source sits on\n"
	             "    the non-linear form, so the next residual picks it up without\n"
	             "    darcy->Assemble() being called again. Nothing in the loop calls\n"
	             "    DarcyForm::Update() or Reset().\n" );

	if ( failures > 0 )
		std::printf( "\n  *** %d case(s) did not converge; their legs are not reported.\n",
		             failures );

	return ( failures > 0 ) ? 1 : 0;
}
