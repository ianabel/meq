/*
 * J^{-1} AT A RETAINED LINEARISATION, WITH THE BORDER.
 * GradShafranovSolver::linearise() factorises the bordered Newton
 * matrix of a converged state -- the field Jacobian, psi_ax's column, its row
 * and corner -- and keeps it; solveLinearised() and linearisedTangent() reuse
 * it, one back-solve a direction; linearisedInverseNorm() estimates 1/sigma_min
 * from it. MEASUREMENTS.md M-LIN-3 and M-LIN-4.
 *
 *   * the bordered tangent -- dx and d psi_ax -- against a Richardson
 *     difference of two solves converged to round-off;
 *   * the retained factorisation reused across directions, to the bit of a
 *     fresh one, and refused once anything has evaluated a residual since;
 *   * a NORMALISED closed form through the bordered path: a manufactured
 *     source F = ( lambda^2/span ) g0 ( 1 + c ( Psi - phi ) ) + e ( psi^2 -
 *     lambda^2 phi^2 ) g0, g0 = -Delta* phi, whose psi = lambda phi has its
 *     axis where phi peaks, so psi_ax = lambda and the exact tangent in
 *     lambda is phi with d psi_ax/d lambda = 1 -- and F reads the
 *     normalisation, so the border is load-bearing;
 *   * sigma_min: the estimate against the distance to the least eigenvalue
 *     of the linearised operator, M-LIN-2's near-fold case.
 */
#define BOOST_TEST_MODULE LinearisedSolve
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "convergence/ConvergenceHarness.hpp"

namespace
{
	using meq::tests::standardBox;

	/// phi = sin( a ( R - R_0 ) ) sin( b ( z - z_0 ) ) on the standard box:
	/// zero on Gamma, one at its centre, and g0 = -Delta* phi.
	struct Mode
	{
		double radiusLow = 0.6, zLow = -0.6, width = 0.8, height = 1.2;
		double a() const { return M_PI/width; }
		double b() const { return M_PI/height; }
		double phi( double radius, double z ) const
		{
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

	std::unique_ptr<Solved> solveManufactured( double lambda, int order, int cells )
	{
		auto s = std::make_unique<Solved>();
		s->mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( standardBox(), cells ) );
		s->h = standardBox().width()/cells;
		s->source = std::make_unique<ManufacturedNormalised>( lambda, 0.3, 0.25, 0.9*lambda );
		s->zero = std::make_unique<mfem::ConstantCoefficient>( 0.0 );
		Mode const mode;
		s->guess = std::make_unique<mfem::FunctionCoefficient>(
			[ mode, lambda ]( mfem::Vector const &x ) { return 0.9*lambda*mode.phi( x( 0 ), x( 1 ) ); } );
		s->solver = std::make_unique<meq::GradShafranovSolver>( *s->mesh, order );
		s->solver->setSource( *s->source, 0.9*lambda );
		s->solver->setBoundaryData( *s->zero );
		s->solver->setInitialGuess( *s->guess );
		s->solver->setNewtonControl( 1.0e-14, 1.0e-15, 40 );
		s->solver->solve();
		return s;
	}

	/// The bordered tangent in lambda at a solved state.
	void tangent( Solved &s, mfem::Vector &dx, double &dAxis )
	{
		s.solver->linearise( s.solver->state() );
		mfem::Vector dRdp;
		ManufacturedNormalised const *source = s.source.get();
		s.solver->assembleSourceDerivative( s.solver->state(),
			[ source ]( double radius, double z, double psi ) { return source->dFdLambda( radius, z, psi ); },
			dRdp );
		std::vector<double> ds;
		s.solver->linearisedTangent( dRdp, {}, nullptr, dx, ds );
		dAxis = ds.front();
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
	// The harness's step for a difference: eps^(1/3) in single, the
	// historical 1e-3 in double.
	double const h = meq::tests::jacobianDifferenceStep( 1.0e-3 )*lambda;
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
	 * IN SINGLE A DIFFERENCE OF TWO SOLVES RESOLVES THE DERIVATIVE TO ABOUT
	 * 1e-3, NOT 1e-8: each solve is good to a few eps and the step is
	 * eps^(1/3). Measured 4.6e-04 in the state (M-LIN-5), so the bound there
	 * is 1e-2 -- which the border's column dropped, 0.62, still exceeds
	 * sixty-fold. Double keeps 1e-8.
	 */
#ifdef MFEM_USE_SINGLE
	double const bound = 1.0e-2;
#else
	double const bound = 1.0e-8;
#endif
	BOOST_TEST( relative < bound, "the bordered tangent misses the difference by " << relative );
	BOOST_TEST( axisRelative < bound, "d psi_ax misses the difference by " << axisRelative );
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
 * in L-infinity from the second on. M-LIN-3.
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
		// The scales of the tangent's round-off floor: max | phi | = 1 and
		// max | grad phi/R | = a/R_0 at the inner wall.
		double const psiScale = 1.0, qScale = mode.a()/mode.radiusLow;
		for ( std::size_t i = 1; i < rows.size(); ++i )
		{
			double const ratio = rows[ i - 1 ].h/rows[ i ].h;
			auto check = [ & ]( char const *name, double Row::*field, double scale )
			{
				// Under the floor the floor is the bound; in double it is
				// zero and every rate is asserted (M-LIN-5).
				if ( rows[ i ].*field < meq::tests::tangentRoundOffFloor( scale, rows[ i ].h ) )
				{
					BOOST_TEST_MESSAGE( "  k = " << order << ", h = " << rows[ i ].h << ": " << name
					                    << " is under the tangent's round-off floor, rate not asserted" );
					return;
				}
				double const r = rate( rows[ i - 1 ].*field, rows[ i ].*field, ratio );
				BOOST_TEST( r >= wanted, "k = " << order << ": " << name << " at " << r );
			};
			check( "d psi L2", &Row::dL2, psiScale );
			check( "d q L2", &Row::qL2, qScale );
			if ( i >= 2 )
			{
				check( "d psi max", &Row::dMax, psiScale );
				check( "d q max", &Row::qMax, qScale );
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
 * the linearised operator, on M-LIN-2's near-fold case (F = mu psi + g,
 * no border; k = 2 on 8 cells, whose least eigenvalue M-LIN-2 fitted at
 * 23.074977827). Close in it must grow as 1/delta. M-LIN-4.
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
		solver.solve();
		solver.linearise( solver.state() );
		double change = 0.0;
		double const estimate = solver.linearisedInverseNorm( 20, &change );
		std::printf( "    %10.1e %14.4e %14.4e %12.2e\n", delta, estimate, estimate*delta, change );
		// In single the fold is resolved to about delta = 1e-4: closer, the
		// factorisation's own round-off is the size of sigma_min, and the
		// estimate falls short (19 at 1e-5, M-LIN-5). The range is the
		// format's precondition; double holds it to 1e-5.
#ifdef MFEM_USE_SINGLE
		bool const resolved = delta >= 1.0e-4;
#else
		bool const resolved = true;
#endif
		if ( delta <= 1.0e-2 && resolved )
			scaled.push_back( estimate*delta );
	}
	std::fflush( stdout );
	double const low = *std::min_element( scaled.begin(), scaled.end() );
	double const high = *std::max_element( scaled.begin(), scaled.end() );
	BOOST_TEST( high/low < 1.5, "the estimate times delta varied by " << high/low
	            << " over the resolved deltas from 1e-2; near the fold it should be a constant" );
}

/*
 * WITHOUT A BORDER THE RETAINED PATH IS solveTangent()'s, so TangentConvergence's
 * closed-form rates (M-LIN-1), measured through solveTangent(), are this
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
	std::printf( "\n  no border: the retained tangent against solveTangent(), %.3e relative\n", relative );
	BOOST_TEST( relative < 1.0e-13, "the retained path's tangent differs from solveTangent()'s by " << relative );
}
