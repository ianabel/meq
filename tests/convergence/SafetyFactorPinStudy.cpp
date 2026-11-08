#define BOOST_TEST_MODULE SafetyFactorPinStudy
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "mfem.hpp"

#include "meq/FluxExtraction.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorSolve.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * WHICH q CONDITION GIVES WAY TO EACH SCALE PIN -- A STUDY, RUN ONCE, BUILT
 * AND NOT REGISTERED WITH ctest. MEASUREMENTS.md M-355 and M-356.
 *
 * A q-driven equilibrium needs one scalar beside p( Psi_N ) and q( Psi_N ) to
 * select a member of the family they leave (M-353). Each pin adds a
 * condition, or removes an unknown, so to keep the outer system square one q
 * condition gives way: a direction u of the fit's residual space, in the
 * orthonormal coordinates y = R c of the fit (R from the QR of its design
 * matrix on the surfaces' labels). Candidates:
 *
 *   y_k     the fit's k-th orthonormal row: the misfit left is the k-th
 *           orthogonal polynomial over the surfaces;
 *   axis    the coefficient of Psi^d: the misfit left is proportional to
 *           Psi^d, zero to order d at the edge and largest on the axis --
 *           "q at the axis limit".
 *
 * For each pin and each u, the smallest singular value of the square
 * Jacobian at SafetyFactorSolver's reference, columns in the same orthonormal
 * coordinates. Then I_p as [source] PlasmaCurrent poses it -- the profile
 * scale a border unknown of the inner solve -- around the free loop.
 */
namespace
{
	using meq::GradShafranovSolver;

	double const g0 = 2.20;
	double const g1 = 0.55;
	double const pPrimeValue = 0.45;
	double const rowWeight = 10.0;

	meq::FluxFamilyOptions familyOptions()
	{
		meq::FluxFamilyOptions options;
		options.surfaces = 12;
		options.angles = 96;
		return options;
	}

	enum class Inner
	{
		Free,      ///< psi_ax a border unknown
		Held,      ///< the normalisation fixed at ( held, 0 ), plain Newton
		Current    ///< psi_ax and the profile scale border unknowns, mu0 I_p given
	};

	struct Solved
	{
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<GradShafranovSolver> solver;
		double psiAxis = 0.0;
		double current = 0.0;
		double scale = 1.0;
		bool converged = true;
	};

	Solved solve( mfem::Mesh &mesh, int order,
	              std::shared_ptr<meq::Profile const> const &ggPrime,
	              Inner inner, double parameter, Solved const *warm )
	{
		Solved out;
		meq::tests::Rectangle const box = meq::tests::standardBox();
		auto pPrime = std::make_shared<meq::ConstantProfile const>( pPrimeValue );
		out.source = std::make_unique<meq::NormalisedMHDSource>(
			pPrime, ggPrime, inner == Inner::Held ? parameter : 1.0, 1.0 );
		out.solver = std::make_unique<GradShafranovSolver>( mesh, order );

		mfem::FunctionCoefficient guess(
			[ box ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			} );
		mfem::ConstantCoefficient zero( 0.0 );

		if ( inner == Inner::Held )
		{
			out.source->setNormalisation( parameter, 0.0 );
			out.solver->setSource( static_cast<meq::Source const &>( *out.source ) );
		}
		else
			out.solver->setSource( *out.source, warm ? warm->psiAxis : 0.30 );
		if ( inner == Inner::Current )
		{
			out.solver->setPlasmaCurrent( parameter );
			if ( warm != nullptr )
				out.solver->setPlasmaCurrentScaleInitial( warm->scale );
		}
		out.solver->setBoundaryData( zero );
		if ( warm != nullptr )
			out.solver->setInitialGuess( warm->solver->potential() );
		else
			out.solver->setInitialGuess( guess );
		out.solver->setNewtonControl( 1.0e-12, 1.0e-14, 60 );
		try
		{
			out.solver->solve();
		}
		catch ( std::exception const & )
		{
			out.converged = false;
			return out;
		}
		out.solver->postProcess();
		out.psiAxis = inner == Inner::Held ? parameter : out.solver->psiAxis();
		out.current = out.solver->deliveredPlasmaCurrent();
		out.scale = out.solver->plasmaCurrentScale();
		return out;
	}

	std::shared_ptr<meq::Profile const> ggPrimeOf( std::vector<double> const &c )
	{
		return std::make_shared<meq::SplineProfile const>(
			meq::ggPrimeKnotsFromCoefficients( c, 17, meq::EdgeExtension::Clamp ) );
	}

	struct Step
	{
		bool ok = false;
		Eigen::VectorXd fitted;
		double located = 0.0;
		double current = 0.0;
		double flux = 0.0;
		double scale = 1.0;
	};

	/// The delivered toroidal flux, int_plasma g/R dA with g^2 = sum c_j Psi^j,
	/// as the driver's ToroidalFlux pin reads it (M-362).
	double deliveredFlux( GradShafranovSolver const &solver, Eigen::VectorXd const &c )
	{
		double const boundary = solver.psiBoundary();
		double const span = solver.psiAxis() - boundary;
		return solver.integrateOverDomain( [ & ]( double radius, double, double psi )
		{
			double const normalised = ( psi - boundary )/span;
			if ( !( radius > 0.0 ) || !( normalised > 0.0 ) )
				return 0.0;
			double value = 0.0, power = 1.0;
			for ( Eigen::Index j = 0; j < c.size(); ++j )
			{
				value += c( j )*power;
				power *= normalised;
			}
			return value > 0.0 ? std::sqrt( value )/radius : 0.0;
		} );
	}

	struct Fixture
	{
		int const order = 2;
		mfem::Mesh mesh;
		std::vector<double> label, targetQ;
		double referenceAxis = 0.0;
		double referenceCurrent = 0.0;
		double referenceFlux = 0.0;
		std::unique_ptr<Solved> base;
		int solves = 0;

		explicit Fixture( int cells = 24 )
			: mesh( meq::tests::makeMesh( meq::tests::standardBox(), cells ) )
		{
			std::vector<meq::Knot> knots;
			for ( double x : { 0.0, 1.0 } )
				knots.push_back( meq::Knot{ x, g1*( g0 + g1*x ), g1*g1 } );
			base = std::make_unique<Solved>( solve(
				mesh, order, std::make_shared<meq::SplineProfile const>( knots ),
				Inner::Free, 0.0, nullptr ) );
			BOOST_TEST_REQUIRE( base->converged );
			referenceAxis = base->psiAxis;
			referenceCurrent = base->current;
			{
				Eigen::VectorXd reference( 3 );
				reference << g0*g0, 2.0*g0*g1, g1*g1;
				referenceFlux = deliveredFlux( *base->solver, reference );
			}
			meq::FluxFamilyOptions options = familyOptions();
			double const axis = referenceAxis;
			options.toroidalField = [ axis ]( double psi ) { return g0 + g1*psi/axis; };
			meq::FluxSurfaceFamily const family =
				meq::extractFluxSurfaces( *base->solver, options );
			BOOST_TEST_REQUIRE( family.safetyFactorAvailable );
			for ( meq::FluxSurface const &s : family.surfaces )
			{
				label.push_back( s.normalisedFlux );
				targetQ.push_back( s.safetyFactor );
			}
		}

		double target( double psiN ) const
		{
			std::size_t best = 0;
			for ( std::size_t i = 1; i < label.size(); ++i )
				if ( std::abs( label[ i ] - psiN ) < std::abs( label[ best ] - psiN ) )
					best = i;
			return targetQ[ best ];
		}

		Eigen::MatrixXd orthonormal( int degree ) const
		{
			Eigen::MatrixXd v( static_cast<Eigen::Index>( label.size() ), degree + 1 );
			for ( std::size_t i = 0; i < label.size(); ++i )
			{
				double power = 1.0;
				for ( int k = 0; k <= degree; ++k )
				{
					v( static_cast<Eigen::Index>( i ), k ) = power;
					power *= 1.0 - label[ i ];
				}
			}
			Eigen::HouseholderQR<Eigen::MatrixXd> qr( v );
			Eigen::MatrixXd r = qr.matrixQR().topRows( degree + 1 )
			                      .triangularView<Eigen::Upper>();
			return r/std::sqrt( static_cast<double>( label.size() ) );
		}

		Step evaluate( Eigen::VectorXd const &c, Inner inner, double parameter,
		               bool keep = false )
		{
			Step out;
			std::vector<double> coefficients( c.data(), c.data() + c.size() );
			++solves;
			auto attempt = std::make_unique<Solved>( solve(
				mesh, order, ggPrimeOf( coefficients ), inner, parameter, base.get() ) );
			if ( !attempt->converged )
				return out;
			meq::FluxSurfaceFamily family;
			try
			{
				family = meq::extractFluxSurfaces( *attempt->solver, familyOptions() );
			}
			catch ( std::exception const & )
			{
				return out;
			}
			meq::ToroidalField const inverted = meq::invertSafetyFactor(
				family, [ this ]( double psiN ) { return target( psiN ); } );
			std::vector<double> const fitted = meq::fitToroidalFieldSquared(
				inverted, static_cast<unsigned int>( c.size() - 1 ) );
			out.fitted = Eigen::Map<Eigen::VectorXd const>(
				fitted.data(), static_cast<Eigen::Index>( fitted.size() ) );
			out.located = family.psiAxis - family.psiBoundary;
			out.current = attempt->current;
			out.flux = deliveredFlux( *attempt->solver, c );
			out.scale = attempt->scale;
			out.ok = true;
			if ( keep )
			{
				attempt->psiAxis = attempt->solver->psiAxis();
				base = std::move( attempt );
			}
			return out;
		}
	};

	Eigen::VectorXd answer( int degree )
	{
		Eigen::VectorXd c = Eigen::VectorXd::Zero( degree + 1 );
		c( 0 ) = g0*g0;
		c( 1 ) = 2.0*g0*g1;
		c( 2 ) = g1*g1;
		return c;
	}

	/// d( G - c )/dc in c, and the scalars' derivatives, centrally.
	struct Jacobian
	{
		Eigen::MatrixXd fit;
		Eigen::RowVectorXd located, current, flux;
	};

	Jacobian differentiate( Fixture &fixture, Eigen::VectorXd const &c, Inner inner,
	                        double parameter, double step )
	{
		int const n = static_cast<int>( c.size() );
		Jacobian out;
		out.fit.resize( n, n );
		out.located.resize( n );
		out.current.resize( n );
		out.flux.resize( n );
		for ( int j = 0; j < n; ++j )
		{
			Eigen::VectorXd plus = c, minus = c;
			plus( j ) += step;
			minus( j ) -= step;
			Step const up = fixture.evaluate( plus, inner, parameter );
			Step const down = fixture.evaluate( minus, inner, parameter );
			BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
			out.fit.col( j ) = ( up.fitted - down.fitted )/( 2.0*step );
			out.fit( j, j ) -= 1.0;
			out.located( j ) = ( up.located - down.located )/( 2.0*step );
			out.current( j ) = ( up.current - down.current )/( 2.0*step );
			out.flux( j ) = ( up.flux - down.flux )/( 2.0*step );
		}
		return out;
	}

	/// An orthonormal basis of the complement of @a u, as rows.
	Eigen::MatrixXd complementOf( Eigen::VectorXd u )
	{
		u.normalize();
		Eigen::Index const n = u.size();
		Eigen::MatrixXd projector = Eigen::MatrixXd::Identity( n, n ) - u*u.transpose();
		Eigen::JacobiSVD<Eigen::MatrixXd> svd( projector, Eigen::ComputeFullU );
		return svd.matrixU().leftCols( n - 1 ).transpose();
	}

	struct Sigma
	{
		double smallest, largest;
	};

	Sigma sigmaOf( Eigen::MatrixXd const &m )
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd( m );
		return { svd.singularValues()( svd.singularValues().size() - 1 ),
		         svd.singularValues()( 0 ) };
	}
}

BOOST_AUTO_TEST_CASE( whichQConditionGivesWayToEachPin )
{
	Fixture fixture;
	std::printf( "\n  M-355: smallest singular value of each pinned square system "
	             "at the reference (box n = 24, k = 2), weight %.0f on the pin row\n",
	             rowWeight );
	std::printf( "    psi_ax* %.9e, mu0 I_p* %.9e, Phi* %.9e\n", fixture.referenceAxis,
	             fixture.referenceCurrent, fixture.referenceFlux );
	for ( int degree : { 2, 4 } )
	{
		Eigen::VectorXd const c = answer( degree );
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		Eigen::MatrixXd const rInverse = r.inverse();
		Jacobian const free = differentiate( fixture, c, Inner::Free, 0.0, 1.0e-3 );
		Jacobian const held = differentiate( fixture, c, Inner::Held,
		                                     fixture.referenceAxis, 1.0e-3 );
		Sigma const plain = sigmaOf( r*free.fit*rInverse );
		std::printf( "    degree %d: free (no pin) sigma_min %.3e, sigma_max %.3e\n", degree,
		             plain.smallest, plain.largest );
		std::printf( "      %-8s %22s %22s %22s %22s %22s\n", "gives way", "R0B0 (c0 held)",
		             "Phi_edge, psi_ax free", "Phi_edge, psi_ax held", "I_p, psi_ax free",
		             "Phi delivered" );

		Eigen::MatrixXd const columnsFree = r.bottomRightCorner( degree, degree ).inverse();
		for ( int k = 0; k <= degree + 1; ++k )
		{
			Eigen::VectorXd u = Eigen::VectorXd::Zero( degree + 1 );
			if ( k <= degree )
				u( k ) = 1.0;
			else
				u = r.col( degree );
			Eigen::MatrixXd const keep = complementOf( u );

			Eigen::MatrixXd edge = keep*( r*free.fit ).rightCols( degree )*columnsFree;
			Eigen::MatrixXd flux( degree + 1, degree + 1 ), fluxHeld( degree + 1, degree + 1 ),
				current( degree + 1, degree + 1 ), delivered( degree + 1, degree + 1 );
			flux.topRows( degree ) = keep*r*free.fit;
			flux.row( degree ) = rowWeight*free.located/fixture.referenceAxis;
			fluxHeld.topRows( degree ) = keep*r*held.fit;
			fluxHeld.row( degree ) = rowWeight*held.located/fixture.referenceAxis;
			current.topRows( degree ) = keep*r*free.fit;
			current.row( degree ) = rowWeight*free.current/fixture.referenceCurrent;
			delivered.topRows( degree ) = keep*r*free.fit;
			delivered.row( degree ) = rowWeight*free.flux/fixture.referenceFlux;
			Sigma const a = sigmaOf( edge );
			Sigma const b = sigmaOf( flux*rInverse );
			Sigma const h = sigmaOf( fluxHeld*rInverse );
			Sigma const i = sigmaOf( current*rInverse );
			Sigma const f = sigmaOf( delivered*rInverse );
			char name[ 16 ];
			if ( k <= degree )
				std::snprintf( name, sizeof name, "y_%d", k );
			else
				std::snprintf( name, sizeof name, "axis" );
			std::printf( "      %-8s %10.3e/%10.3e %10.3e/%10.3e %10.3e/%10.3e %10.3e/%10.3e %10.3e/%10.3e\n",
			             name, a.smallest, a.largest, b.smallest, b.largest, h.smallest,
			             h.largest, i.smallest, i.largest, f.smallest, f.largest );
		}
		std::printf( "      d ln psi_ax/dc = " );
		for ( Eigen::Index j = 0; j <= degree; ++j )
			std::printf( " %10.3e", free.located( j )/fixture.referenceAxis );
		std::printf( "\n      d ln Phi/dc    = " );
		for ( Eigen::Index j = 0; j <= degree; ++j )
			std::printf( " %10.3e", free.flux( j )/fixture.referenceFlux );
		std::printf( "\n      d ln I_p/dc    = " );
		for ( Eigen::Index j = 0; j <= degree; ++j )
			std::printf( " %10.3e", free.current( j )/fixture.referenceCurrent );
		std::printf( "\n" );
	}
	std::printf( "    %d solves\n", fixture.solves );
}

/*
 * [source] PlasmaCurrent AS IT IS POSED TODAY, beside the free loop: the
 * profile scale alpha a border unknown of every inner solve. alpha multiplies
 * p' and gg' together, so a converged loop at alpha != 1 has a pressure alpha p
 * and a g^2 whose varying part is alpha times the one it fitted.
 */
BOOST_AUTO_TEST_CASE( thePlasmaCurrentBorderBesideTheFreeLoop )
{
	Fixture fixture;
	std::printf( "\n  M-356: [source] PlasmaCurrent's border inside the free q loop\n" );
	for ( double factor : { 1.0, 1.1 } )
	{
		double const target = factor*fixture.referenceCurrent;
		std::vector<double> start = { 2.0*2.0, 0.0, 0.0 };
		int evaluations = 0;
		auto map = [ & ]( std::vector<double> const &c ) -> std::vector<double>
		{
			++evaluations;
			Eigen::VectorXd point = Eigen::Map<Eigen::VectorXd const>(
				c.data(), static_cast<Eigen::Index>( c.size() ) );
			for ( int i = 0; i <= 64; ++i )
			{
				double const psi = i/64.0;
				double value = c[ 0 ] + c[ 1 ]*psi + c[ 2 ]*psi*psi;
				if ( !( value > 0.0 ) )
					return { c[ 0 ] + 1.0, c[ 1 ] + 1.0, c[ 2 ] + 1.0 };
			}
			Step const at = fixture.evaluate( point, Inner::Current, target, true );
			if ( !at.ok )
				return { c[ 0 ] + 1.0, c[ 1 ] + 1.0, c[ 2 ] + 1.0 };
			return { at.fitted( 0 ), at.fitted( 1 ), at.fitted( 2 ) };
		};
		meq::OuterNewtonOptions options;
		options.maxIterations = 15;
		meq::OuterNewtonResult const result = meq::solveForToroidalField( start, map, options );
		std::vector<double> const &c = result.coefficients;
		Step const at = fixture.evaluate(
			Eigen::Map<Eigen::VectorXd const>( c.data(), 3 ), Inner::Current, target, true );
		double worst = 0.0;
		for ( double psi : fixture.label )
		{
			double const s = 1.0 - psi;
			double const varying = c[ 1 ]*s + c[ 2 ]*s*s;
			worst = std::max( worst, std::abs( std::sqrt( ( c[ 0 ] + at.scale*varying )
			                                              /( c[ 0 ] + varying ) ) - 1.0 ) );
		}
		std::printf( "    I_p = %.2f I_p*: %s\n", factor, result.status.c_str() );
		std::printf( "      alpha %.9f, psi_ax %.9e (reference %.9e), c = %.6e %.6e %.6e, "
		             "conditioning %.3e\n", at.scale, fixture.base->psiAxis,
		             fixture.referenceAxis, c[ 0 ], c[ 1 ], c[ 2 ], result.jacobianConditioning );
		std::printf( "      the equilibrium's q against the target, with g^2's varying part "
		             "alpha times the fitted one: worst %.3e\n", worst );
	}
	std::printf( "    %d solves\n", fixture.solves );
}
