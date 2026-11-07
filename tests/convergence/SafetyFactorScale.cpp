#define BOOST_TEST_MODULE SafetyFactorScale
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
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * D1'S TORUS CHECK -- MANTA-EQUILIBRIUM-PLAN.md section 5, D1.
 *
 * In a cylinder, fixed-boundary equilibria with the same p( Psi_N ) and
 * q( Psi_N ) form a one-parameter family: the toroidal field's scale is free
 * (check_cylinder_scale.py). In a torus the same is VMEC's statement that
 * p( s ), iota( s ) and the boundary need PHIEDGE beside them -- and with q
 * and Psi_N fixed, Phi_edge = 2 pi ( psi_ax - psi_bnd ) int_0^1 q dPsi_N, so
 * the free scalar is psi_ax.
 *
 * MEQ's q-driven solve leaves psi_ax free and fits ALL of g^2, its constant
 * included. This file measures what that does, on SafetyFactorSolver's
 * fixture (g = 2.20 + 0.55 Psi, so g^2 is exactly quadratic):
 *
 *   1. the outer Jacobian G' - I at the answer, psi_ax free, by degree, in an
 *      orthonormal basis of the fit's own space -- whether it is near-singular
 *      along the continuum, and whether that falls with the degree;
 *   2. psi_ax HELD at S = ( 1 + delta ) psi_ax*, the least-squares fixed point
 *      of [ G_S( c ) - c ; w ( psi_ax,located - S )/S ] -- whether a member
 *      of the family exists off the reference, to what residual, and how well
 *      the held problem determines c;
 *   3. the same with a source scale alpha as a square system's last unknown,
 *      which is the plan's rejected alternative.
 */
namespace
{
	using meq::GradShafranovSolver;

	double const g0 = 2.20;
	double const g1 = 0.55;
	double const pPrimeValue = 0.45;

	meq::FluxFamilyOptions familyOptions()
	{
		meq::FluxFamilyOptions options;
		options.surfaces = 12;
		options.angles = 96;
		return options;
	}

	struct Solved
	{
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<GradShafranovSolver> solver;
		double psiAxis = 0.0;
		bool converged = true;
	};

	/// psi_ax FREE (held < 0): the bordered Newton, as SafetyFactorSolver.
	/// HELD (held > 0): the normalisation fixed at ( held, 0 ) and the source
	/// solved as a plain one, scaled by alpha.
	Solved solve( mfem::Mesh &mesh, int order,
	              std::shared_ptr<meq::Profile const> const &ggPrime,
	              double held, double alpha, Solved const *warm )
	{
		Solved out;
		meq::tests::Rectangle const box = meq::tests::standardBox();
		auto pPrime = std::make_shared<meq::ConstantProfile const>( pPrimeValue );
		out.source = std::make_unique<meq::NormalisedMHDSource>(
			pPrime, ggPrime, held > 0.0 ? held : 1.0, 1.0 );
		out.source->setCurrentScale( alpha );
		out.solver = std::make_unique<GradShafranovSolver>( mesh, order );

		mfem::FunctionCoefficient guess(
			[ box ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			} );
		mfem::ConstantCoefficient zero( 0.0 );

		if ( held > 0.0 )
		{
			out.source->setNormalisation( held, 0.0 );
			out.solver->setSource( static_cast<meq::Source const &>( *out.source ) );
		}
		else
			out.solver->setSource( *out.source, warm ? warm->psiAxis : 0.30 );
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
		out.psiAxis = held > 0.0 ? held : out.solver->psiAxis();
		return out;
	}

	/// The polynomial g^2( Psi ) = sum c_k Psi^k as gg' knots.
	std::shared_ptr<meq::Profile const> ggPrimeOf( std::vector<double> const &c )
	{
		return std::make_shared<meq::SplineProfile const>(
			meq::ggPrimeKnotsFromCoefficients( c, 17, meq::EdgeExtension::Clamp ) );
	}

	/// One map evaluation: the fitted g^2 and the located axis's flux.
	struct Step
	{
		bool ok = false;
		Eigen::VectorXd fitted;
		double located = 0.0;
	};

	struct Fixture
	{
		int const n;
		int const order = 2;
		mfem::Mesh mesh;
		std::vector<double> label, targetQ;
		double referenceAxis = 0.0;
		std::unique_ptr<Solved> base;
		int solves = 0;

		explicit Fixture( int cells = 24 )
			: n( cells ), mesh( meq::tests::makeMesh( meq::tests::standardBox(), cells ) )
		{
			// THE REFERENCE AND ITS q, exactly as SafetyFactorSolver builds them.
			std::vector<meq::Knot> knots;
			for ( double x : { 0.0, 1.0 } )
				knots.push_back( meq::Knot{ x, g1*( g0 + g1*x ), g1*g1 } );
			base = std::make_unique<Solved>( solve(
				mesh, order, std::make_shared<meq::SplineProfile const>( knots ),
				-1.0, 1.0, nullptr ) );
			BOOST_TEST_REQUIRE( base->converged );
			referenceAxis = base->psiAxis;
			meq::FluxFamilyOptions options = familyOptions();
			double const axis = referenceAxis;
			options.toroidalField = [ axis ]( double psi )
			{
				return g0 + g1*psi/axis;
			};
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

		/// R of the QR of the fit's design matrix on the labels, in the
		/// source's Psi: y = R c is an orthonormal coordinate on P_d in the
		/// discrete L2 norm over the surfaces, so singular values in y are
		/// the map's and not the monomial basis's.
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
			r /= std::sqrt( static_cast<double>( label.size() ) );
			return r;
		}

		Step evaluate( Eigen::VectorXd const &c, double held, double alpha,
		               bool keep = false )
		{
			Step out;
			std::vector<double> coefficients( c.data(), c.data() + c.size() );
			for ( int i = 0; i <= 64; ++i )
			{
				double const psi = i/64.0;
				double value = 0.0, power = 1.0;
				for ( double ck : coefficients )
				{
					value += ck*power;
					power *= psi;
				}
				if ( !( value > 0.0 ) )
					return out;
			}
			++solves;
			auto attempt = std::make_unique<Solved>( solve(
				mesh, order, ggPrimeOf( coefficients ), held, alpha, base.get() ) );
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
			out.located = family.psiAxis;
			out.ok = true;
			if ( keep )
			{
				attempt->psiAxis = family.psiAxis;
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

	void printVector( char const *name, Eigen::VectorXd const &v )
	{
		std::printf( "      %-18s", name );
		for ( Eigen::Index i = 0; i < v.size(); ++i )
			std::printf( " %11.4e", v( i ) );
		std::printf( "\n" );
	}
}

namespace
{
	struct Held
	{
		bool solved = false;
		Eigen::VectorXd c;
		Eigen::MatrixXd jacobian;
		Eigen::VectorXd residual;
	};

	/// The least-squares fixed point with psi_ax held at @a held: Gauss-Newton
	/// on [ R( G_S( c ) - c ) ; w ( psi_ax,located - S )/S ], R the fit's
	/// orthonormalising factor.
	Held solveHeld( Fixture &fixture, int degree, double held, double weight,
	                Eigen::VectorXd c )
	{
		Held out;
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		auto assemble = [ & ]( Step const &at, Eigen::VectorXd const &point )
		{
			Eigen::VectorXd value( degree + 2 );
			value.head( degree + 1 ) = r*( at.fitted - point );
			value( degree + 1 ) = weight*( at.located - held )/held;
			return value;
		};
		out.jacobian.resize( degree + 2, degree + 1 );
		for ( int iteration = 0; iteration < 8; ++iteration )
		{
			Step const at = fixture.evaluate( c, held, 1.0, true );
			if ( !at.ok )
				return out;
			out.residual = assemble( at, c );
			for ( int j = 0; j <= degree; ++j )
			{
				double const step = 1.0e-3;
				Eigen::VectorXd plus = c, minus = c;
				plus( j ) += step;
				minus( j ) -= step;
				Step const up = fixture.evaluate( plus, held, 1.0 );
				Step const down = fixture.evaluate( minus, held, 1.0 );
				if ( !up.ok || !down.ok )
					return out;
				out.jacobian.col( j ) = ( assemble( up, plus ) - assemble( down, minus ) )
				                        /( 2.0*step );
			}
			Eigen::VectorXd const dc =
				out.jacobian.colPivHouseholderQr().solve( -out.residual );
			c += dc;
			if ( dc.norm() < 1.0e-9*c.norm() )
				break;
		}
		Step const at = fixture.evaluate( c, held, 1.0, true );
		if ( !at.ok )
			return out;
		out.residual = assemble( at, c );
		out.c = c;
		out.solved = true;
		return out;
	}
}

namespace
{
	/// G' - I at @a c, psi_ax free, in the fit's orthonormal coordinates; and
	/// d psi_ax / d c, column by column.
	struct FreeJacobian
	{
		Eigen::MatrixXd inY;
		Eigen::VectorXd located;
	};

	FreeJacobian freeJacobian( Fixture &fixture, Eigen::VectorXd const &c, double step )
	{
		int const degree = static_cast<int>( c.size() ) - 1;
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		Eigen::MatrixXd jacobian( degree + 1, degree + 1 );
		FreeJacobian out;
		out.located.resize( degree + 1 );
		for ( int j = 0; j <= degree; ++j )
		{
			Eigen::VectorXd plus = c, minus = c;
			plus( j ) += step;
			minus( j ) -= step;
			Step const up = fixture.evaluate( plus, -1.0, 1.0 );
			Step const down = fixture.evaluate( minus, -1.0, 1.0 );
			BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
			jacobian.col( j ) = ( up.fitted - down.fitted )/( 2.0*step );
			jacobian( j, j ) -= 1.0;
			out.located( j ) = ( up.located - down.located )/( 2.0*step );
		}
		out.inY = r*jacobian*r.inverse();
		return out;
	}
}

/*
 * 1. psi_ax FREE, AS MEQ SOLVES q-DRIVEN TODAY: G' - I at the answer has ONE
 *    near-null direction, and it is the scale -- psi_ax and g moving together,
 *    d psi_ax/psi_ax = ( 1/2 ) d c_0/c_0 -- smaller with every degree.
 */
BOOST_AUTO_TEST_CASE( theFreeOuterJacobiansNearNullDirectionIsTheScale )
{
	Fixture fixture;
	std::printf( "\n  D1 TORUS CHECK 1: G' - I at the answer, psi_ax free (M-SESSION-2)\n" );
	std::printf( "    %6s %11s %11s %11s %13s\n", "degree", "sigma_min", "sigma_max",
	             "ratio", "dlnpsi/dlnc0" );
	for ( int degree : { 2, 4 } )
	{
		Eigen::VectorXd const c = answer( degree );
		FreeJacobian const at = freeJacobian( fixture, c, 1.0e-3 );
		Eigen::JacobiSVD<Eigen::MatrixXd> svd( at.inY, Eigen::ComputeFullV );
		Eigen::VectorXd const sigma = svd.singularValues();
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		Eigen::VectorXd const nullC = r.inverse()*svd.matrixV().col( degree );
		double const scale = ( at.located.dot( nullC )/fixture.referenceAxis )
		                     /( nullC( 0 )/c( 0 ) );
		std::printf( "    %6d %11.3e %11.3e %11.3e %13.4f\n", degree, sigma( degree ),
		             sigma( 0 ), sigma( degree )/sigma( 0 ), scale );
		BOOST_TEST( sigma( degree )/sigma( 0 ) < 1.0e-3,
			"at degree " << degree << " the free q-driven Jacobian's smallest singular "
			"value is " << sigma( degree )/sigma( 0 ) << " of its largest: the scale "
			"continuum of fixed p( Psi_N ), q( Psi_N ) equilibria is gone, and D1's "
			"reason for holding psi_ax with it" );
		BOOST_TEST( std::abs( scale - 0.5 ) < 0.1,
			"the near-null direction moves psi_ax by " << scale << " of c_0 relatively, "
			"where the scale continuum, g proportional to psi_ax, moves it by 1/2" );
	}
	std::printf( "    %d solves\n", fixture.solves );
}

/*
 * 2. psi_ax HELD 10% off the reference: a member of the family exists, to a
 *    residual that falls with the degree of g^2 (8x from 2 to 4, 29x to 5,
 *    M-SESSION-2) -- the representation pins
 *    the scale, not the torus -- and holding psi_ax raises the smallest
 *    singular value by two orders. alpha, a source scale, as the square
 *    system's last unknown is printed: the plan rejects it as near-singular.
 */
BOOST_AUTO_TEST_CASE( holdingPsiAxisSelectsAMemberOfTheFamily )
{
	Fixture fixture;
	double const weight = 10.0;
	double const delta = 0.10;
	double const held = ( 1.0 + delta )*fixture.referenceAxis;
	std::printf( "\n  D1 TORUS CHECK 2: psi_ax held at 1.10 psi_ax*, least squares "
	             "(M-SESSION-2)\n" );
	std::printf( "    %6s %11s %11s %11s %11s %11s %11s\n", "degree", "|fit res|",
	             "axis res", "s_min", "s_max", "alpha ratio", "free s_min" );
	std::vector<double> residuals;
	for ( int degree : { 2, 4 } )
	{
		Held const at = solveHeld( fixture, degree, held, weight,
		                           answer( degree )*( 1.0 + delta )*( 1.0 + delta ) );
		BOOST_TEST_REQUIRE( at.solved );
		Eigen::JacobiSVD<Eigen::MatrixXd> tall( at.jacobian );

		Eigen::MatrixXd withAlpha( degree + 2, degree + 2 );
		withAlpha.leftCols( degree + 1 ) = at.jacobian;
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		{
			double const step = 1.0e-4;
			Step const up = fixture.evaluate( at.c, held, 1.0 + step );
			Step const down = fixture.evaluate( at.c, held, 1.0 - step );
			BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
			withAlpha.col( degree + 1 ).head( degree + 1 ) =
				r*( up.fitted - down.fitted )/( 2.0*step );
			withAlpha( degree + 1, degree + 1 ) =
				weight*( up.located - down.located )/( 2.0*step*held );
		}
		Eigen::JacobiSVD<Eigen::MatrixXd> alpha( withAlpha );

		FreeJacobian const free = freeJacobian( fixture, answer( degree ), 1.0e-3 );
		Eigen::JacobiSVD<Eigen::MatrixXd> freeSvd( free.inY );

		double const residual = at.residual.head( degree + 1 ).norm();
		residuals.push_back( residual );
		std::printf( "    %6d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n", degree, residual,
		             at.residual( degree + 1 )/weight, tall.singularValues()( degree ),
		             tall.singularValues()( 0 ),
		             alpha.singularValues()( degree + 1 )/alpha.singularValues()( 0 ),
		             freeSvd.singularValues()( degree ) );
		BOOST_TEST( tall.singularValues()( degree ) > 30.0*freeSvd.singularValues()( degree ),
			"at degree " << degree << ", holding psi_ax leaves the smallest singular value "
			<< tall.singularValues()( degree ) << " against " << freeSvd.singularValues()( degree )
			<< " free: holding it no longer determines the scale" );
	}
	BOOST_TEST( residuals[ 1 ] < 0.25*residuals[ 0 ],
		"10% along the continuum the held fixed point's residual is " << residuals[ 0 ]
		<< " at degree 2 and " << residuals[ 1 ] << " at degree 4: if it does not fall with "
		"the representation, the torus pins the scale and D1's continuum is not there" );
	std::printf( "    %d solves\n", fixture.solves );
}

/*
 * 3. AND THE MESH DOES NOT PIN IT: the same residual, 5% along, at n = 12 and
 *    24 (and 48 in M-SESSION-2), each against its own mesh's reference.
 */
BOOST_AUTO_TEST_CASE( theMeshDoesNotPinTheScale )
{
	std::printf( "\n  D1 TORUS CHECK 3: the residual 5%% along the continuum, by resolution "
	             "(M-SESSION-2)\n" );
	std::vector<double> residuals;
	for ( int cells : { 12, 24 } )
	{
		Fixture fixture( cells );
		double const held = 0.95*fixture.referenceAxis;
		Held const at = solveHeld( fixture, 2, held, 10.0, answer( 2 )*0.95*0.95 );
		BOOST_TEST_REQUIRE( at.solved );
		residuals.push_back( at.residual.head( 3 ).norm() );
		std::printf( "    n = %2d: |fit res| %.3e, axis res %.3e\n", cells, residuals.back(),
		             at.residual( 3 )/10.0 );
	}
	BOOST_TEST( std::abs( residuals[ 1 ] - residuals[ 0 ] ) < 0.1*residuals[ 1 ],
		"the held residual moves from " << residuals[ 0 ] << " to " << residuals[ 1 ]
		<< " with the mesh: what pins the scale is then the discretisation" );
}
