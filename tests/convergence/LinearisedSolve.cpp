/*
 * J^{-1} AT A RETAINED LINEARISATION, WITH THE BORDER.
 * GradShafranovSolver::linearise() factorises the bordered Newton
 * matrix of a converged state -- the field Jacobian, psi_ax's column, its row
 * and corner -- and keeps it; solveLinearised() and linearisedTangent() reuse
 * it, one back-solve a direction; linearisedInverseNorm() estimates 1/sigma_min
 * from it. MEASUREMENTS.md M-375 and M-376.
 *
 *   * the bordered tangent -- dx and d psi_ax -- against a Richardson
 *     difference of two solves converged to round-off;
 *   * every border in it: psi_bnd pinned at a limiter point, and at the
 *     two-lobe mode's X-point with ( R_X, z_X ), each border unknown against
 *     the same difference (M-401);
 *   * the retained factorisation reused across directions, to the bit of a
 *     fresh one, and refused once anything has evaluated a residual since;
 *   * a NORMALISED closed form through the bordered path: a manufactured
 *     source F = ( lambda^2/span ) g0 ( 1 + c ( Psi - phi ) ) + e ( psi^2 -
 *     lambda^2 phi^2 ) g0, g0 = -Delta* phi, whose psi = lambda phi has its
 *     axis where phi peaks, so psi_ax = lambda and the exact tangent in
 *     lambda is phi with d psi_ax/d lambda = 1 -- and F reads the
 *     normalisation, so the border is load-bearing;
 *   * sigma_min: the estimate against the distance to the least eigenvalue
 *     of the linearised operator, M-374's near-fold case.
 */
#define BOOST_TEST_MODULE LinearisedSolve
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "convergence/ConvergenceHarness.hpp"
#include "unit/DeviceKernels.hpp"

/*
 * MEQ_TEST_DEVICE, WHEN SET, CONFIGURES AN mfem::Device BEFORE ANY CASE RUNS
 * (NewtonConvergence.cpp's pattern; a function-local static, since a
 * namespace-scope mfem::Device may be constructed before MFEM's own
 * singleton). Unset, there is none, which is how the plain ctest runs. Set,
 * it is the device arm tests/CMakeLists.txt registers: the retained linearisation
 * read its back-solves and refinement residuals raw while they were
 * device-valid, and the bordered tangent came out wrong under --device cuda
 * (MEASUREMENTS.md M-383).
 */
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
}

BOOST_TEST_GLOBAL_FIXTURE( DeviceFromEnvironment );

namespace
{
	using meq::tests::standardBox;

	/// phi = sin( a ( R - R_0 ) ) sin( b ( z - z_0 ) ) on the standard box:
	/// zero on Gamma, one at its centre, and g0 = -Delta* phi.
	struct Mode
	{
		double radiusLow = 0.6, zLow = -0.6, width = 0.8, height = 1.2;
		/// THE TWO-LOBE MODE, where mu or nu is non-zero: sin( b( z - z_0 ) )
		/// becomes S( t ) = sin t + mu sin 3t + nu sin 2t, two maxima in z with
		/// a saddle between them -- an X-point for the border to carry. Zero
		/// leaves the single mode, expression for expression.
		double mu = 0.0, nu = 0.0;
		double a() const { return M_PI/width; }
		double b() const { return M_PI/height; }
		bool lobed() const { return mu != 0.0 || nu != 0.0; }
		double lobes( double z ) const
		{
			double const t = b()*( z - zLow );
			return std::sin( t ) + mu*std::sin( 3.0*t ) + nu*std::sin( 2.0*t );
		}
		double lobesSecond( double z ) const
		{
			double const t = b()*( z - zLow );
			return -b()*b()*( std::sin( t ) + 9.0*mu*std::sin( 3.0*t ) + 4.0*nu*std::sin( 2.0*t ) );
		}
		double phi( double radius, double z ) const
		{
			if ( lobed() )
				return std::sin( a()*( radius - radiusLow ) )*lobes( z );
			return std::sin( a()*( radius - radiusLow ) )*std::sin( b()*( z - zLow ) );
		}
		void gradPhi( double radius, double z, double &dr, double &dz ) const
		{
			dr = a()*std::cos( a()*( radius - radiusLow ) )*std::sin( b()*( z - zLow ) );
			dz = b()*std::sin( a()*( radius - radiusLow ) )*std::cos( b()*( z - zLow ) );
		}
		/// -Delta* phi, Delta* = d_RR - ( 1/R ) d_R + d_zz.
		double g0( double radius, double z ) const
		{
			if ( lobed() )
			{
				double const sr = std::sin( a()*( radius - radiusLow ) );
				double const cr = std::cos( a()*( radius - radiusLow ) );
				return a()*a()*sr*lobes( z ) + ( a()/radius )*cr*lobes( z ) - sr*lobesSecond( z );
			}
			double const cr = std::cos( a()*( radius - radiusLow ) );
			double const sz = std::sin( b()*( z - zLow ) );
			return ( a()*a() + b()*b() )*phi( radius, z ) + ( a()/radius )*cr*sz;
		}
	};

	/*
	 * THE MANUFACTURED NORMALISED SOURCE. F = ( lambda^2/span ) g0 ( 1 + c ( Psi
	 * - phi ) ) + e ( psi^2 - lambda^2 phi^2 ) g0, Psi = ( psi - psi_bnd )/span:
	 * at psi = lambda phi, psi_ax = lambda, psi_bnd = 0 it is lambda g0 =
	 * -Delta*( lambda phi ), so that is a solution with its axis at phi's peak;
	 * and F reads the span, so psi_ax is a genuine unknown whose column and row
	 * the tangent needs.
	 *
	 * THE e TERM IS THERE TO BREAK A SYMMETRY. Without it F is homogeneous of
	 * degree one in ( psi, psi_ax, lambda ), so the DISCRETE solution is
	 * exactly lambda phi_h and the discrete tangent exactly psi_h/lambda: the
	 * tangent's rates would be the solution's by construction. With it they
	 * are measured.
	 */
	class ManufacturedNormalised : public meq::NormalisedSource
	{
		public:
			ManufacturedNormalised( double lambdaIn, double cIn, double eIn, double psiAxis )
				: lambda( lambdaIn ), c( cIn ), e( eIn ), axis( psiAxis ) {}

			double f( double radius, double z, double psi ) const override
			{
				if ( !insidePlasma( radius, z, psi ) )
					return 0.0;
				double const span = axis - boundary;
				double const x = ( psi - boundary )/span;
				double const phi = mode.phi( radius, z );
				return lambda*lambda/span*mode.g0( radius, z )*( 1.0 + c*( x - phi ) )
				       + e*( psi*psi - lambda*lambda*phi*phi )*mode.g0( radius, z );
			}
			double dFdPsi( double radius, double z, double psi ) const override
			{
				if ( !insidePlasma( radius, z, psi ) )
					return 0.0;
				double const span = axis - boundary;
				return lambda*lambda/( span*span )*mode.g0( radius, z )*c + 2.0*e*psi*mode.g0( radius, z );
			}
			bool normalisationDerivatives( double radius, double z, double psi, double &dAxis,
			                               double &dBoundary ) const override
			{
				double const span = axis - boundary;
				double const x = ( psi - boundary )/span;
				double const k = lambda*lambda*mode.g0( radius, z );
				double const bracket = 1.0 + c*( x - mode.phi( radius, z ) );
				// span = a - b; dx/da = -x/span, dx/db = ( x - 1 )/span.
				dAxis = -k*bracket/( span*span ) + k/span*c*( -x/span );
				dBoundary = k*bracket/( span*span ) + k/span*c*( ( x - 1.0 )/span );
				return true;
			}
			/// dF/d lambda at fixed psi and labels.
			double dFdLambda( double radius, double z, double psi ) const
			{
				double const span = axis - boundary;
				double const x = ( psi - boundary )/span;
				double const phi = mode.phi( radius, z );
				return 2.0*lambda/span*mode.g0( radius, z )*( 1.0 + c*( x - phi ) )
				       - 2.0*e*lambda*phi*phi*mode.g0( radius, z );
			}
			void setNormalisation( double psiAxis, double psiBoundary ) override
			{
				axis = psiAxis;
				boundary = psiBoundary;
			}
			using meq::NormalisedSource::setNormalisation;
			double normalisation() const override { return axis; }
			double boundaryNormalisation() const override { return boundary; }

			Mode mode;

		private:
			double lambda, c, e, axis, boundary = 0.0;
	};

	struct Solved
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<ManufacturedNormalised> source;
		std::unique_ptr<mfem::ConstantCoefficient> zero;
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		double h = 0.0;
	};

	/// Which psi_bnd the solve closes: the box's zero (psi_bnd is not an
	/// unknown), a limiter point inside the box, or the two-lobe mode's
	/// X-point.
	enum class Edge { Fixed, Limiter, XPoint };

	/// The two-lobe mode's saddle, ( R, z ): the sine's peak in R and
	/// S'( t ) = 0 between the lobes in z, by Newton from t = pi/2.
	void saddleOf( Mode const &mode, double &radius, double &z )
	{
		double t = 0.5*M_PI;
		for ( int i = 0; i < 50; ++i )
		{
			double const d1 = std::cos( t ) + 3.0*mode.mu*std::cos( 3.0*t ) + 2.0*mode.nu*std::cos( 2.0*t );
			double const d2 = -std::sin( t ) - 9.0*mode.mu*std::sin( 3.0*t ) - 4.0*mode.nu*std::sin( 2.0*t );
			t -= d1/d2;
		}
		radius = mode.radiusLow + 0.5*mode.width;
		z = mode.zLow + t/mode.b();
	}

	Mode modeFor( Edge edge )
	{
		Mode mode;
		if ( edge == Edge::XPoint )
		{
			mode.mu = 0.4;
			mode.nu = 0.15;
		}
		return mode;
	}

	std::unique_ptr<Solved> solveManufactured( double lambda, int order, int cells, Edge edge = Edge::Fixed )
	{
		auto s = std::make_unique<Solved>();
		s->mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( standardBox(), cells ) );
		s->h = standardBox().width()/cells;
		Mode const mode = modeFor( edge );
		double const peak = edge == Edge::XPoint ? 1.1 : 1.0;
		s->source = std::make_unique<ManufacturedNormalised>( lambda, 0.3, 0.25, 0.9*peak*lambda );
		s->source->mode = mode;
		s->zero = std::make_unique<mfem::ConstantCoefficient>( 0.0 );
		s->guess = std::make_unique<mfem::FunctionCoefficient>(
			[ mode, lambda ]( mfem::Vector const &x ) { return 0.9*lambda*mode.phi( x( 0 ), x( 1 ) ); } );
		s->solver = std::make_unique<meq::GradShafranovSolver>( *s->mesh, order );
		if ( edge == Edge::Limiter )
			s->solver->setBoundaryFluxPoint( 1.3, 0.0 );
		if ( edge == Edge::XPoint )
		{
			double radius = 0.0, z = 0.0;
			saddleOf( mode, radius, z );
			s->solver->setXPointBoundary( radius, z );
			s->solver->setBoundaryFluxInitial( 0.9*lambda*mode.phi( radius, z ) );
		}
		s->solver->setSource( *s->source, 0.9*peak*lambda );
		s->solver->setBoundaryData( *s->zero );
		s->solver->setInitialGuess( *s->guess );
		s->solver->setNewtonControl( 1.0e-14, 1.0e-15, 40 );
		if ( edge == Edge::XPoint )
		{
			/*
			 * XP-2'S START: psi_bnd pinned at the saddle, then the border from
			 * there. From the guess itself the bordered Newton's first dense
			 * corner is singular -- every entry finite, rows 0 and 3 nearly
			 * parallel -- on every lobe shape tried (TODO, *The X-point
			 * border from a cold start*); the tangent is the converged
			 * state's either way.
			 */
			double radius = 0.0, z = 0.0;
			saddleOf( mode, radius, z );
			meq::GradShafranovSolver pinned( *s->mesh, order );
			pinned.setBoundaryFluxPoint( radius, z );
			pinned.setSource( *s->source, 0.9*peak*lambda );
			pinned.setBoundaryData( *s->zero );
			pinned.setInitialGuess( *s->guess );
			pinned.setNewtonControl( 1.0e-14, 1.0e-15, 40 );
			meq::tests::namedRoutesFromEnvironment( pinned );
			pinned.solve();
			s->solver->setInitialState( pinned.state(), {} );
			s->solver->setBoundaryFluxInitial( pinned.psiBoundary() );
			s->source->setNormalisation( pinned.psiAxis(), pinned.psiBoundary() );
		}
		meq::tests::namedRoutesFromEnvironment( *s->solver );
		s->solver->solve();
		return s;
	}

	/// The bordered tangent in lambda at a solved state.
	void tangent( Solved &s, mfem::Vector &dx, double &dAxis,
	              int passes = meq::GradShafranovSolver::linearisedRefinementPasses )
	{
		s.solver->linearise( s.solver->state() );
		mfem::Vector dRdp;
		ManufacturedNormalised const *source = s.source.get();
		s.solver->assembleSourceDerivative( s.solver->state(),
			[ source ]( double radius, double z, double psi ) { return source->dFdLambda( radius, z, psi ); },
			dRdp );
		std::vector<double> ds;
		s.solver->linearisedTangent( dRdp, {}, nullptr, dx, ds, passes );
		dAxis = ds.front();
		// Every digit, so that two arms -- the host and MEQ_TEST_DEVICE --
		// can be compared to the bit by their output (M-383).
		std::printf( "  tangent bits: d psi_ax %.17g, |dx| %.17g\n", dAxis, dx.Norml2() );
		// Read raw by the caller: asked for on the host, as any raw reader
		// must, and after the norm, whose device read protects the host page
		// under --device debug.
		dx.HostRead();
	}

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}
}

/*
 * THE BORDERED TANGENT IS THE DERIVATIVE OF THE DISCRETE SOLUTION, psi_ax
 * included: against ( 4 D( h/2 ) - D( h ) )/3 of two solves at lambda +- h,
 * converged to round-off on the same mesh.
 */
BOOST_AUTO_TEST_CASE( theBorderedTangentIsADifferenceOfTwoConvergedSolves )
{
	double const lambda = 0.2;
	int const order = 2, cells = 8;
	auto base = solveManufactured( lambda, order, cells );
	BOOST_REQUIRE( base->solver->linearisedBorders() == 0 );
	mfem::Vector dx;
	double dAxis = 0.0;
	tangent( *base, dx, dAxis );
	BOOST_TEST( base->solver->linearisedBorders() == 1 );

	auto at = [ & ]( double by, double &axisOut )
	{
		auto s = solveManufactured( lambda + by, order, cells );
		axisOut = s->solver->psiAxis();
		return mfem::Vector( s->solver->state() );
	};
	double const h = 1.0e-3*lambda;
	auto central = [ & ]( double step, double &axisDifference )
	{
		double ap = 0.0, am = 0.0;
		mfem::Vector d = at( step, ap );
		d -= at( -step, am );
		d /= 2.0*step;
		axisDifference = ( ap - am )/( 2.0*step );
		return d;
	};
	double axisHalf = 0.0, axisFull = 0.0;
	mfem::Vector richardson = central( h/2.0, axisHalf );
	richardson *= 4.0;
	richardson -= central( h, axisFull );
	richardson /= 3.0;
	double const axisRichardson = ( 4.0*axisHalf - axisFull )/3.0;

	mfem::Vector error( richardson );
	error -= dx;
	double const relative = error.Norml2()/dx.Norml2();
	double const axisRelative = std::abs( axisRichardson - dAxis )/std::abs( dAxis );
	std::printf( "  the bordered tangent against two converged solves: state %.3e relative, "
	             "d psi_ax %.10f against %.10f (%.3e)\n", relative, dAxis, axisRichardson, axisRelative );
	/*
	 * A PRECONDITION OF THE FORMAT, NOT A TOLERANCE. The reference is a
	 * difference of two solves, and in single each solve stops at its own
	 * floor, so the reference itself is good to about 2e-3 over this step
	 * (M-378) -- it cannot judge a tangent to 1e-8. Single checks the
	 * tangent against the closed form instead, at the solution's rates.
	 */
#ifdef MFEM_USE_SINGLE
	BOOST_TEST_MESSAGE( "  single: a difference of two single solves cannot resolve the tangent to "
	                    "1e-8; not asserted (measured " << relative << ")" );
#else
	BOOST_TEST( relative < 1.0e-8, "the bordered tangent misses the difference by " << relative );
	BOOST_TEST( axisRelative < 1.0e-8, "d psi_ax misses the difference by " << axisRelative );
#endif
}

/*
 * EVERY BORDER IS IN THE TANGENT: psi_bnd pinned at a limiter point, and at
 * the two-lobe mode's X-point with ( R_X, z_X ) two more unknowns -- the
 * bordered Newton's blocks, re-assembled at the state by linearise() -- each
 * against the Richardson difference of two solves converged to round-off,
 * the state and every border unknown.
 */
namespace
{
	void bordersAgainstSolves( Edge edge, char const *name )
	{
		double const lambda = 0.2;
		// Nine cells, not eight: on eight the saddle at R = 1 sits on a mesh
		// line, where q_h jumps across the face and the X-point's Newton
		// alternates between the two elements' roots and does not converge.
		int const order = 2, cells = 9;
		auto base = solveManufactured( lambda, order, cells, edge );
		base->solver->linearise( base->solver->state() );
		meq::GradShafranovSolver::LinearisedBorders const layout = base->solver->linearisedLayout();
		BOOST_TEST_REQUIRE( layout.axis == 0 );
		BOOST_TEST_REQUIRE( layout.boundary == 1 );
		BOOST_TEST_REQUIRE( ( layout.xPoint == 2 ) == ( edge == Edge::XPoint ) );
		BOOST_TEST( layout.count == ( edge == Edge::XPoint ? 4 : 2 ) );
		mfem::Vector dRdp, dx;
		ManufacturedNormalised const *source = base->source.get();
		base->solver->assembleSourceDerivative( base->solver->state(),
			[ source ]( double radius, double z, double psi ) { return source->dFdLambda( radius, z, psi ); },
			dRdp );
		std::vector<double> ds;
		base->solver->linearisedTangent( dRdp, {}, nullptr, dx, ds );

		/// The state and the border unknowns, as the solve leaves them.
		auto at = [ & ]( double by, std::vector<double> &borders )
		{
			auto s = solveManufactured( lambda + by, order, cells, edge );
			borders = { s->solver->psiAxis(), s->solver->psiBoundary() };
			if ( edge == Edge::XPoint )
			{
				borders.push_back( s->solver->xPointR() );
				borders.push_back( s->solver->xPointZ() );
			}
			return mfem::Vector( s->solver->state() );
		};
		auto central = [ & ]( double step, std::vector<double> &borders )
		{
			std::vector<double> plus, minus;
			mfem::Vector d = at( step, plus );
			d -= at( -step, minus );
			d /= 2.0*step;
			borders.resize( plus.size() );
			for ( std::size_t i = 0; i < plus.size(); ++i )
				borders[ i ] = ( plus[ i ] - minus[ i ] )/( 2.0*step );
			return d;
		};
		double const h = 1.0e-3*lambda;
		std::vector<double> half, full;
		mfem::Vector richardson = central( h/2.0, half );
		richardson *= 4.0;
		richardson -= central( h, full );
		richardson /= 3.0;
		mfem::Vector error( richardson );
		error -= dx;
		double const relative = error.Norml2()/dx.Norml2();
		std::printf( "\n  %s: the bordered tangent against two converged solves, state %.3e relative\n",
		             name, relative );
		char const *labels[] = { "d psi_ax", "d psi_bnd", "d R_X", "d z_X" };
		double worst = 0.0;
		for ( std::size_t i = 0; i < half.size(); ++i )
		{
			double const reference = ( 4.0*half[ i ] - full[ i ] )/3.0;
			// Against the border's own size or the state's: d z_X may be
			// zero by the mode's symmetry in R.
			double const apart = std::abs( reference - ds[ i ] )
			                     /std::max( std::abs( reference ), static_cast<double>( dx.Normlinf() ) );
			worst = std::max( worst, apart );
			std::printf( "    %-10s %+.10e against %+.10e (%.3e)\n", labels[ i ], ds[ i ], reference, apart );
		}
		std::fflush( stdout );
#ifdef MFEM_USE_SINGLE
		BOOST_TEST_MESSAGE( "  single: a difference of two single solves cannot resolve the tangent to "
		                    "1e-8; not asserted (measured " << relative << ", borders " << worst << ")" );
#else
		BOOST_TEST( relative < 1.0e-8, name << ": the bordered tangent misses the difference by " << relative );
		BOOST_TEST( worst < 1.0e-8, name << ": a border unknown misses the difference by " << worst );
#endif
	}
}

BOOST_AUTO_TEST_CASE( theLimiterBorderIsInTheTangent ) { bordersAgainstSolves( Edge::Limiter, "limiter point" ); }
BOOST_AUTO_TEST_CASE( theXPointBorderIsInTheTangent ) { bordersAgainstSolves( Edge::XPoint, "X-point" ); }

/*
 * setInitialState()'S ITERATE IS THE ONE THE SOLVE STARTS FROM, on any Device:
 * a converged state handed to a fresh solver of the same problem starts it
 * near its answer. Under --device debug it was not -- the copy into the
 * state went to the device while the trace datum, written through an alias
 * whose flags still said host, put the guess back -- and every warm start, the
 * q-driven loop's and a restart's, began from the guess (M-402).
 */
BOOST_AUTO_TEST_CASE( theInitialStateIsWhereTheSolveStarts )
{
	auto base = solveManufactured( 0.2, 2, 8 );
	double const scale = base->solver->newtonResiduals().front();
	meq::GradShafranovSolver again( *base->mesh, 2 );
	again.setSource( *base->source, base->solver->psiAxis() );
	again.setBoundaryData( *base->zero );
	again.setInitialGuess( *base->guess );
	again.setNewtonControl( 1.0e-14, 1.0e-15, 40 );
	meq::tests::namedRoutesFromEnvironment( again );
	again.setInitialState( base->solver->state(), {} );
	again.solve();
	double const first = again.newtonResiduals().front();
	std::printf( "  its first step: field %.3e, axis constraint %.3e\n",
	             again.borderSteps().empty() ? -1.0 : again.borderSteps().front().fieldNorm,
	             again.borderSteps().empty() ? -1.0 : again.borderSteps().front().axis );
	std::printf( "\n  setInitialState(): the first residual %.3e against the cold start's %.3e, %d iterations "
	             "against %d\n", first, scale, again.newtonIterations(), base->solver->newtonIterations() );
	/*
	 * NOT AT ZERO ON THE HOST, AND NOTHING IS DROPPED: the field residual of
	 * the first step is 1.5e-16 on every arm, and psi_ax's label the base's.
	 * The 1.9e-4 is the axis constraint alone. On eight cells the axis sits
	 * at a mesh vertex (R = 1, z = 0), the fresh solver's first location is
	 * a cold sweep, and on the host it roots q_h in another element than the
	 * converged solve tracked -- psi_h jumps by O( h^( k+1 ) ) across the
	 * faces there -- so it starts 1.9e-4 off and converges back in 4; the
	 * card's device route happens to root the same element, 1.5e-16 in none.
	 * On 7 and 9 cells, the axis off every face, the host starts at 1.8e-15
	 * and 6.2e-16 in 0 iterations; on 10, a vertex again, 8.0e-5 in 4
	 * (M-405). From the guess the first residual is the cold start's own.
	 */
	BOOST_TEST( first < 1.0e-2*scale, "a converged state handed to setInitialState() starts at "
	            << first << ", the cold start's " << scale << ": the solve began from the guess" );
	BOOST_TEST( again.newtonIterations() < base->solver->newtonIterations() );
}

/*
 * THE FACTORISATION IS RETAINED, AND RETAINED TO THE BIT: a direction solved
 * again after another answers as it did, and as a fresh linearisation of the
 * same state answers. And it is refused once a residual has been evaluated
 * on the solver, since that overwrote the blocks it lives in.
 */
BOOST_AUTO_TEST_CASE( theRetainedFactorisationAnswersAsAFreshOneDoes )
{
	auto s = solveManufactured( 0.2, 2, 8 );
	s->solver->linearise( s->solver->state() );
	int const size = s->solver->state().Size();
	mfem::Vector b1( size ), b2( size ), x1, x2, x1Again, x1Fresh;
	for ( int i = 0; i < size; ++i )
	{
		b1( i ) = std::sin( 0.31*i );
		b2( i ) = std::cos( 0.17*i + 0.4 );
	}
	std::vector<double> d1, d2, d1Again, d1Fresh;
	s->solver->solveLinearised( b1, { 0.3 }, x1, d1 );
	s->solver->solveLinearised( b2, { -0.1 }, x2, d2 );
	s->solver->solveLinearised( b1, { 0.3 }, x1Again, d1Again );

	s->solver->linearise( s->solver->state() );
	s->solver->solveLinearised( b1, { 0.3 }, x1Fresh, d1Fresh );

	std::printf( "  retained bits: ds %.17g, |x| %.17g\n", d1.front(), x1.Norml2() );
	for ( mfem::Vector const *x : { &x1, &x1Again, &x1Fresh } )
		x->HostRead();   // read raw below, after the norm's device read
	bool again = d1Again == d1, fresh = d1Fresh == d1;
	for ( int i = 0; i < size; ++i )
	{
		again = again && x1Again( i ) == x1( i );
		fresh = fresh && x1Fresh( i ) == x1( i );
	}
	BOOST_TEST( again, "the same direction through the retained factorisation answered differently" );
	BOOST_TEST( fresh, "a fresh linearisation of the same state answered differently from the retained one" );

	mfem::Vector residual;
	s->solver->evaluateResidual( s->solver->state(), residual );
	BOOST_TEST( !s->solver->isLinearised() );
	BOOST_CHECK_THROW( s->solver->solveLinearised( b1, { 0.3 }, x1, d1 ), std::logic_error );
}

/*
 * THE CLOSED FORM, THROUGH THE BORDER. delta psi_h and delta q_h against phi
 * and grad phi/R, d psi_ax against 1, k = 1...3 on 4 to 32 cells -- k = 3
 * from 8, since its solution is pre-asymptotic on 4 (2.9 in L2 from 4 to 8,
 * the tangent with it), and its finest rung about 1e-6 in L-infinity.
 * Asserted as the solution's rates are: k + 1 - 0.2 in L2 on every pair and
 * in L-infinity from the second on. M-375.
 */
namespace
{
	void borderedRates( int order, std::vector<int> const &meshes )
	{
		double const lambda = 0.2;
		Mode const mode;
		mfem::FunctionCoefficient phi( [ mode ]( mfem::Vector const &x ) { return mode.phi( x( 0 ), x( 1 ) ); } );
		mfem::VectorFunctionCoefficient q( 2, [ mode ]( mfem::Vector const &x, mfem::Vector &v )
		{
			double dr = 0.0, dz = 0.0;
			mode.gradPhi( x( 0 ), x( 1 ), dr, dz );
			v( 0 ) = dr/x( 0 );
			v( 1 ) = dz/x( 0 );
		} );
		std::vector<mfem::IntegrationRule const *> irs( mfem::Geometry::NumGeom );
		for ( int i = 0; i < mfem::Geometry::NumGeom; ++i )
			irs[ i ] = &mfem::IntRules.Get( i, 2*order + 6 );

		struct Row
		{
			double h, psiL2, psiMax, dL2, dMax, qL2, qMax, axis;
			/// psi_h's own flux errors, against lambda grad phi/R.
			double solutionQL2 = 0.0, solutionQMax = 0.0;
		};
		std::vector<Row> rows;
		for ( int cells : meshes )
		{
			auto s = solveManufactured( lambda, order, cells );
			mfem::Vector dx;
			double dAxis = 0.0;
			tangent( *s, dx, dAxis );
			int const fluxSize = s->solver->fluxSpace().GetVSize();
			mfem::GridFunction dpsi( &s->solver->potentialSpace() ), dq( &s->solver->fluxSpace() );
			for ( int i = 0; i < dpsi.Size(); ++i )
				dpsi( i ) = dx( fluxSize + i );
			for ( int i = 0; i < fluxSize; ++i )
				dq( i ) = -dx( i );
			mfem::FunctionCoefficient psiExact( [ mode, lambda ]( mfem::Vector const &x )
			                                    { return lambda*mode.phi( x( 0 ), x( 1 ) ); } );
			Row row;
			row.h = s->h;
			row.psiL2 = s->solver->potential().ComputeL2Error( psiExact, irs.data() );
			row.psiMax = s->solver->potential().ComputeMaxError( psiExact, irs.data() );
			row.dL2 = dpsi.ComputeL2Error( phi, irs.data() );
			row.dMax = dpsi.ComputeMaxError( phi, irs.data() );
			row.qL2 = dq.ComputeL2Error( q, irs.data() );
			row.qMax = dq.ComputeMaxError( q, irs.data() );
			row.axis = std::abs( dAxis - 1.0 );
			mfem::GridFunction qh( &s->solver->fluxSpace() );
			for ( int i = 0; i < fluxSize; ++i )
				qh( i ) = -s->solver->state()( i );
			mfem::VectorFunctionCoefficient qSolution( 2, [ mode, lambda ]( mfem::Vector const &x, mfem::Vector &v )
			{
				double dr = 0.0, dz = 0.0;
				mode.gradPhi( x( 0 ), x( 1 ), dr, dz );
				v( 0 ) = lambda*dr/x( 0 );
				v( 1 ) = lambda*dz/x( 0 );
			} );
			row.solutionQL2 = qh.ComputeL2Error( qSolution, irs.data() );
			row.solutionQMax = qh.ComputeMaxError( qSolution, irs.data() );
			rows.push_back( row );
		}

		std::printf( "\n  k = %d: psi_h = lambda phi and its bordered tangent in lambda, errors (rates)\n", order );
		std::printf( "    %6s | %-20s %-20s | %-20s %-20s | %-20s %-20s | %-20s\n", "h", "psi L2", "psi max",
		             "dpsi L2", "dpsi max", "dq L2", "dq max", "|d psi_ax - 1|" );
		for ( std::size_t i = 0; i < rows.size(); ++i )
		{
			auto cell = [ & ]( double Row::*field )
			{
				double const e = rows[ i ].*field;
				if ( i == 0 )
					std::printf( " %9.3e (    )   ", e );
				else
					std::printf( " %9.3e (%4.2f)   ", e, rate( rows[ i - 1 ].*field, e, rows[ i - 1 ].h/rows[ i ].h ) );
			};
			std::printf( "    %6.4f |", rows[ i ].h );
			cell( &Row::psiL2 );
			cell( &Row::psiMax );
			std::printf( "|" );
			cell( &Row::dL2 );
			cell( &Row::dMax );
			std::printf( "|" );
			cell( &Row::qL2 );
			cell( &Row::qMax );
			std::printf( "|" );
			cell( &Row::axis );
			std::printf( "\n" );
		}
		std::fflush( stdout );

		double const wanted = order + 1.0 - meq::tests::rateSlack;
		for ( std::size_t i = 1; i < rows.size(); ++i )
		{
			double const ratio = rows[ i - 1 ].h/rows[ i ].h;
			/*
			 * Held to the solution's rate where the solution has one: a pair
			 * on which psi_h's own error in that field and norm misses the
			 * rate is past the format's floor for the solution, and the
			 * tangent, refined against the same residual, can be no better.
			 * It holds on every pair in double (M-378).
			 */
			auto check = [ & ]( char const *name, double Row::*tangentField, double Row::*solutionField )
			{
				double const solutionRate = rate( rows[ i - 1 ].*solutionField, rows[ i ].*solutionField, ratio );
				if ( solutionRate < wanted )
				{
					std::printf( "    k = %d, h = %.4f: %s not asserted -- the solution's own rate there is %.2f\n",
					             order, rows[ i ].h, name, solutionRate );
					return;
				}
				double const r = rate( rows[ i - 1 ].*tangentField, rows[ i ].*tangentField, ratio );
				BOOST_TEST( r >= wanted, "k = " << order << ": " << name << " at " << r );
			};
			check( "d psi L2", &Row::dL2, &Row::psiL2 );
			check( "d q L2", &Row::qL2, &Row::solutionQL2 );
			if ( i >= 2 )
			{
				check( "d psi max", &Row::dMax, &Row::psiMax );
				check( "d q max", &Row::qMax, &Row::solutionQMax );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theBorderedTangentConvergesAtDegreeOne ) { borderedRates( 1, { 4, 8, 16, 32 } ); }
BOOST_AUTO_TEST_CASE( theBorderedTangentConvergesAtDegreeTwo ) { borderedRates( 2, { 4, 8, 16, 32 } ); }
BOOST_AUTO_TEST_CASE( theBorderedTangentConvergesAtDegreeThree ) { borderedRates( 3, { 8, 16, 32 } ); }

/*
 * sigma_min FROM THE FACTORISATION THE TANGENT ALREADY HAS: the estimate of
 * 1/sigma_min against the relative distance delta to the least eigenvalue of
 * the linearised operator, on M-374's near-fold case (F = mu psi + g,
 * no border; k = 2 on 8 cells, whose least eigenvalue M-374 fitted at
 * 23.074977827). Close in it must grow as 1/delta. M-376.
 */
namespace
{
	class AffineSource : public meq::Source
	{
		public:
			explicit AffineSource( double muIn ) : mu( muIn ) {}
			double f( double radius, double z, double psi ) const override
			{
				return mu*psi + mode.g0( radius, z ) - mu*mode.phi( radius, z );
			}
			double dFdPsi( double, double, double ) const override
			{
				return mu;
			}
		private:
			double mu;
			Mode mode;
	};
}

BOOST_AUTO_TEST_CASE( theInverseNormEstimateGrowsAsTheOperatorNearsItsFold )
{
	double const eigen = 23.074977827233;
	std::printf( "\n  1/sigma_min FROM THE RETAINED FACTORISATION, k = 2 on 8 cells\n" );
	std::printf( "    %10s %14s %14s %12s\n", "delta", "estimate", "x delta", "last change" );
	std::vector<double> scaled;
	for ( double delta : { 1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5 } )
	{
		AffineSource source( eigen*( 1.0 - delta ) );
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), 8 );
		mfem::ConstantCoefficient zero( 0.0 );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( zero );
		solver.setNewtonControl( 1.0e-14, 1.0e-15, 40 );
		meq::tests::namedRoutesFromEnvironment( solver );
		solver.solve();
		solver.linearise( solver.state() );
		double change = 0.0;
		double const estimate = solver.linearisedInverseNorm( 20, &change );
		std::printf( "    %10.1e %14.4e %14.4e %12.2e   bits %.17g\n", delta, estimate, estimate*delta, change,
		             estimate );
		if ( delta <= 1.0e-2 )
			scaled.push_back( estimate*delta );
	}
	std::fflush( stdout );
	double const low = *std::min_element( scaled.begin(), scaled.end() );
	double const high = *std::max_element( scaled.begin(), scaled.end() );
	BOOST_TEST( high/low < 1.5, "the estimate times delta varied by " << high/low
	            << " over delta = 1e-2 to 1e-5; near the fold it should be a constant" );
}

/*
 * WITHOUT A BORDER THE RETAINED PATH IS solveTangent()'s, so TangentConvergence's
 * closed-form rates (M-373), measured through solveTangent(), are this
 * path's too: on the affine case away from its fold, the same tangent from
 * both, with a non-zero Dirichlet datum derivative to exercise the lifting.
 */
BOOST_AUTO_TEST_CASE( withoutABorderTheRetainedTangentIsSolveTangents )
{
	AffineSource source( 10.0 );
	Mode const mode;
	mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), 8 );
	mfem::ConstantCoefficient zero( 0.0 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( zero );
	solver.setNewtonControl( 1.0e-14, 1.0e-15, 40 );
	meq::tests::namedRoutesFromEnvironment( solver );
	solver.solve();
	mfem::Vector dRdp;
	solver.assembleSourceDerivative( solver.state(),
		[ mode ]( double radius, double z, double ) { return mode.g0( radius, z ); }, dRdp );
	mfem::FunctionCoefficient datum( []( mfem::Vector const &x ) { return 0.1*x( 0 )*x( 1 ); } );
	mfem::Vector fresh, retainedTangent;
	solver.solveTangent( solver.state(), dRdp, &datum, fresh );
	solver.linearise( solver.state() );
	BOOST_TEST( solver.linearisedBorders() == 0 );
	std::vector<double> ds;
	solver.linearisedTangent( dRdp, {}, &datum, retainedTangent, ds );
	mfem::Vector difference( fresh );
	difference -= retainedTangent;
	double const relative = difference.Normlinf()/fresh.Normlinf();
	std::printf( "\n  no border: the retained tangent against solveTangent(), %.3e relative; bits |dx| %.17g\n",
	             relative, retainedTangent.Norml2() );
	BOOST_TEST( relative < 1.0e-13, "the retained path's tangent differs from solveTangent()'s by " << relative );
}

/*
 * ITERATIVE REFINEMENT, BY PASS: the bordered tangent's error against the
 * closed form after 0...4 passes of solveLinearised()'s refinement -- each a
 * residual of the bordered system with J applied, and one back-solve -- at
 * k = 2 and 3 on 8 to 32 cells, beside the solution's own error. MEASUREMENT
 * ONLY (disabled; run by name, in each precision). M-378.
 */
BOOST_AUTO_TEST_CASE( theTangentErrorByRefinementPass, * boost::unit_test::disabled() )
{
	double const lambda = 0.2;
	Mode const mode;
	mfem::FunctionCoefficient phi( [ mode ]( mfem::Vector const &x ) { return mode.phi( x( 0 ), x( 1 ) ); } );
	mfem::FunctionCoefficient psiExact( [ mode, lambda ]( mfem::Vector const &x )
	                                    { return lambda*mode.phi( x( 0 ), x( 1 ) ); } );
	std::printf( "\n  THE BORDERED TANGENT BY REFINEMENT PASS: d psi L2 against phi, |d psi_ax - 1|\n" );
	for ( int order : { 2, 3 } )
		for ( int cells : { 8, 16, 32 } )
		{
			auto s = solveManufactured( lambda, order, cells );
			std::vector<mfem::IntegrationRule const *> irs( mfem::Geometry::NumGeom );
			for ( int i = 0; i < mfem::Geometry::NumGeom; ++i )
				irs[ i ] = &mfem::IntRules.Get( i, 2*order + 6 );
			double const psiError = s->solver->potential().ComputeL2Error( psiExact, irs.data() );
			std::printf( "    k = %d, %2d cells: psi_h %.3e (relative %.3e) | d psi:", order, cells, psiError,
			             psiError/lambda );
			for ( int passes = 0; passes <= 4; ++passes )
			{
				mfem::Vector dx;
				double dAxis = 0.0;
				tangent( *s, dx, dAxis, passes );
				int const fluxSize = s->solver->fluxSpace().GetVSize();
				mfem::GridFunction dpsi( &s->solver->potentialSpace() );
				for ( int i = 0; i < dpsi.Size(); ++i )
					dpsi( i ) = dx( fluxSize + i );
				std::printf( " %d: %.3e (%.1e)", passes, dpsi.ComputeL2Error( phi, irs.data() ),
				             std::abs( dAxis - 1.0 ) );
			}
			std::printf( "\n" );
			std::fflush( stdout );
		}
}

/*
 * THE RETAINED LINEARISATION'S DEVICE ROUTE IS ITS HOST ROUTE (M-409): a
 * NormalisedMHDSource, whose POD form lets the retained operator take the
 * device integrand and the source columns the plasma loops' kernels, solved
 * on 8 cells at k = 2. The same tangent -- along a gg' direction, the q
 * loop's derivative -- with the route on and off, and that direction's dR/du
 * both ways. Under `debug` in double the kernels are the host loops on the
 * host's arithmetic and the columns agree to the bit; in single, and on the
 * card, to the format's round-off; the tangents to the round-off of two
 * operators that sum their elements in different orders.
 * Needs a Device (MEQ_TEST_DEVICE) and device kernels, and says so where
 * either is missing.
 */
BOOST_AUTO_TEST_CASE( theDeviceRouteIsTheHostRoute, *meq::tests::requiresDeviceKernels() )
{
	if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
	{
		std::printf( "\n  theDeviceRouteIsTheHostRoute: SKIPPED -- no Device is configured; "
		             "MEQ_TEST_DEVICE=debug or cuda runs it\n" );
		return;
	}
	meq::tests::Rectangle const box = standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, 8 );
	double const g0 = 2.20, g1 = 0.55;
	std::vector<meq::Knot> knots;
	for ( double x : { 0.0, 1.0 } )
		knots.push_back( meq::Knot{ x, g1*( g0 + g1*x ), g1*g1 } );
	meq::NormalisedMHDSource source( std::make_shared<meq::ConstantProfile const>( 0.45 ),
	                                 std::make_shared<meq::SplineProfile const>( knots ), 0.30, 1.0 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setAssemblyMode( meq::GradShafranovSolver::AssemblyMode::Batched );
	solver.setLocalFactorMode( meq::GradShafranovSolver::LocalFactorMode::Batched );
	mfem::FunctionCoefficient guess( [ box ]( mfem::Vector const &x )
	{
		return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
		       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
	} );
	mfem::ConstantCoefficient zero( 0.0 );
	solver.setSource( source, 0.30 );
	solver.setBoundaryData( zero );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-12, 1.0e-14, 60 );
	solver.solve();
	mfem::Vector const state( solver.state() );
	BOOST_TEST_REQUIRE( solver.deviceIntegrandTaken(), "the solve did not take the device integrand, so "
	                    "there is no device route for the linearisation to follow" );

	// A cardinal-like direction in gg': a hat at Psi = 0.5, as a Hermite spline.
	meq::SplineProfile const direction( std::vector<meq::Knot>{ { 0.0, 0.0, 0.0 }, { 0.5, 1.0, 0.0 },
	                                                           { 1.0, 0.0, 0.0 } } );
	auto tangentBy = [ & ]( bool device, mfem::Vector &column, mfem::Vector &dx, double &dAxis )
	{
		solver.setLinearisedDeviceRoute( device );
		long const before = solver.linearisedDeviceTaken();
		solver.linearise( state );
		solver.assembleGGPrimeDerivative( state, direction, column );
		std::vector<double> ds;
		solver.linearisedTangent( column, {}, nullptr, dx, ds );
		dAxis = ds.empty() ? 0.0 : ds.front();
		static_cast<void>( column.HostRead() );
		static_cast<void>( dx.HostRead() );
		return solver.linearisedDeviceTaken() - before;
	};
	mfem::Vector hostColumn, deviceColumn, hostDx, deviceDx;
	double hostAxis = 0.0, deviceAxis = 0.0;
	long const hostTaken = tangentBy( false, hostColumn, hostDx, hostAxis );
	long const deviceTaken = tangentBy( true, deviceColumn, deviceDx, deviceAxis );
	solver.setLinearisedDeviceRoute( true );

	double columnWorst = 0.0, columnScale = 0.0, dxWorst = 0.0, dxScale = 0.0;
	std::size_t columnDiffer = 0;
	for ( int i = 0; i < hostColumn.Size(); ++i )
	{
		columnDiffer += hostColumn( i ) != deviceColumn( i );
		columnWorst = std::max( columnWorst, std::abs( double( hostColumn( i ) - deviceColumn( i ) ) ) );
		columnScale = std::max( columnScale, std::abs( double( hostColumn( i ) ) ) );
	}
	for ( int i = 0; i < hostDx.Size(); ++i )
	{
		dxWorst = std::max( dxWorst, std::abs( double( hostDx( i ) - deviceDx( i ) ) ) );
		dxScale = std::max( dxScale, std::abs( double( hostDx( i ) ) ) );
	}
	bool const card = mfem::Device::Allows( mfem::Backend::CUDA_MASK | mfem::Backend::HIP_MASK );
	double const eps = std::numeric_limits<mfem::real_t>::epsilon();
	std::printf( "\n  the retained route %s: %ld device steps against %ld; dR/du differs at %zu entries, "
	             "by %.2e of %.3e; the tangent by %.2e of %.3e, d psi_ax %.10g against %.10g; "
	             "%d refinement passes kept\n",
	             card ? "on the card" : "under debug", deviceTaken, hostTaken, columnDiffer, columnWorst,
	             columnScale, dxWorst, dxScale, deviceAxis, hostAxis, solver.linearisedPassesTaken() );
	BOOST_TEST( hostTaken == 0L, "the host route took the device " << hostTaken << " times" );
	BOOST_TEST( deviceTaken > 0L, "the device route never ran: linearise() did not take the device "
	            "integrand's operator, or the plasma loops refused" );
	/*
	 * dR/du: the kernel against the host loop. In double under `debug` they
	 * are one arithmetic and agree to the bit; in single the host evaluates
	 * ggPrimeDerivative() in double and the kernel in real_t, on radii the
	 * tables round to real_t (M-322), so the two agree to single's round-off
	 * on either arm, which is what this asserts there.
	 */
#ifdef MFEM_USE_SINGLE
	BOOST_TEST( columnWorst <= 64.0*eps*columnScale, "dR/du on the device route differs from the "
	            "host's by " << columnWorst/columnScale << " of its largest entry" );
#else
	if ( card )
		BOOST_TEST( columnWorst <= 64.0*eps*columnScale, "dR/du on the card differs from the host's by "
		            << columnWorst/columnScale << " of its largest entry" );
	else
		BOOST_TEST( columnDiffer == 0u, "under debug the source kernel is the host loop on the host's "
		            "arithmetic and must agree to the bit; " << columnDiffer << " entries differ" );
#endif
	// The tangent: two operators, the device integrand's and MFEM's, whose
	// element sums run in different orders, through one factorisation each.
	BOOST_TEST( dxWorst <= std::max( 1.0e-11, 1.0e3*eps )*dxScale, "the device route's tangent differs "
	            "from the host's by " << dxWorst/dxScale << " of its largest entry" );
	BOOST_TEST( std::abs( deviceAxis - hostAxis ) <= std::max( 1.0e-11, 1.0e3*eps )*std::abs( hostAxis ),
	            "d psi_ax is " << deviceAxis << " on the device route against " << hostAxis );
}
