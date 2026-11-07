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
		int const n = 24;
		int const order = 2;
		mfem::Mesh mesh;
		std::vector<double> label, targetQ;
		double referenceAxis = 0.0;
		std::unique_ptr<Solved> base;
		int solves = 0;

		Fixture()
			: mesh( meq::tests::makeMesh( meq::tests::standardBox(), 24 ) )
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

BOOST_AUTO_TEST_CASE( theFreeOuterJacobianAlongTheContinuum )
{
	Fixture fixture;
	std::printf( "\n  D1 TORUS CHECK 1: G' - I at the answer, psi_ax free\n" );
	std::printf( "    reference psi_ax %.10e\n", fixture.referenceAxis );
	for ( int degree : { 2, 3, 4, 5 } )
	{
		Eigen::VectorXd const c = answer( degree );
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		for ( double step : { 1.0e-3, 1.0e-2 } )
		{
			Step const at = fixture.evaluate( c, -1.0, 1.0 );
			BOOST_TEST_REQUIRE( at.ok );
			Eigen::MatrixXd jac( degree + 1, degree + 1 );
			Eigen::VectorXd located( degree + 1 );
			for ( int j = 0; j <= degree; ++j )
			{
				Eigen::VectorXd plus = c, minus = c;
				plus( j ) += step;
				minus( j ) -= step;
				Step const up = fixture.evaluate( plus, -1.0, 1.0 );
				Step const down = fixture.evaluate( minus, -1.0, 1.0 );
				BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
				jac.col( j ) = ( up.fitted - down.fitted )/( 2.0*step );
				jac( j, j ) -= 1.0;
				located( j ) = ( up.located - down.located )/( 2.0*step );
			}
			Eigen::MatrixXd const inY = r*jac*r.inverse();
			Eigen::JacobiSVD<Eigen::MatrixXd> svd( inY, Eigen::ComputeFullV );
			Eigen::VectorXd const sigma = svd.singularValues();
			std::printf( "    degree %d, step %.0e: |G(c*) - c*|_y %.3e\n", degree,
			             step, ( r*( at.fitted - c ) ).norm() );
			printVector( "sigma (y)", sigma );
			std::printf( "      sigma_min/sigma_max %.4e\n",
			             sigma( degree )/sigma( 0 ) );
			Eigen::VectorXd const nullY = svd.matrixV().col( degree );
			Eigen::VectorXd const nullC = r.inverse()*nullY;
			printVector( "near-null (c)", nullC/nullC( 0 ) );
			std::printf( "      d psi_ax along it (per unit c_0) %.4e, psi_ax %.4e\n",
			             located.dot( nullC )/nullC( 0 ), fixture.referenceAxis );
		}
	}
	std::printf( "    %d solves\n", fixture.solves );
}

BOOST_AUTO_TEST_CASE( holdingPsiAxisAlongTheContinuum )
{
	Fixture fixture;
	double const reference = fixture.referenceAxis;
	double const weight = 10.0;
	std::printf( "\n  D1 TORUS CHECK 2: psi_ax HELD at ( 1 + delta ) psi_ax*, least "
	             "squares\n" );
	std::printf( "    reference psi_ax %.10e; axis row weighted %.1f\n", reference,
	             weight );
	for ( int degree : { 2, 3, 4, 5 } )
	{
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		std::printf( "    degree %d\n", degree );
		std::printf( "      %7s %11s %11s %11s %11s %11s %11s %11s\n", "delta",
		             "|fit res|", "axis res", "s_min/s_max", "sq s_min/max",
		             "alpha s_min", "g(1)", "g(0)" );
		for ( double delta : { -0.10, -0.05, 0.0, 0.05, 0.10 } )
		{
			double const held = ( 1.0 + delta )*reference;
			// From the scale the cylinder suggests, g^2 ~ S^2.
			Eigen::VectorXd c = answer( degree )*( 1.0 + delta )*( 1.0 + delta );
			// Warm-start the held solves from the reference.
			Eigen::MatrixXd jac( degree + 2, degree + 1 );
			Eigen::VectorXd residual( degree + 2 );
			auto assemble = [ & ]( Step const &at, Eigen::VectorXd const &point )
			{
				Eigen::VectorXd out( degree + 2 );
				out.head( degree + 1 ) = r*( at.fitted - point );
				out( degree + 1 ) = weight*( at.located - held )/held;
				return out;
			};
			bool solved = true;
			for ( int iteration = 0; iteration < 8; ++iteration )
			{
				Step const at = fixture.evaluate( c, held, 1.0, true );
				if ( !at.ok )
				{
					solved = false;
					break;
				}
				residual = assemble( at, c );
				for ( int j = 0; j <= degree; ++j )
				{
					double const step = 1.0e-3;
					Eigen::VectorXd plus = c, minus = c;
					plus( j ) += step;
					minus( j ) -= step;
					Step const up = fixture.evaluate( plus, held, 1.0 );
					Step const down = fixture.evaluate( minus, held, 1.0 );
					if ( !up.ok || !down.ok )
					{
						solved = false;
						break;
					}
					jac.col( j ) = ( assemble( up, plus ) - assemble( down, minus ) )
					               /( 2.0*step );
				}
				if ( !solved )
					break;
				Eigen::VectorXd const dc =
					jac.colPivHouseholderQr().solve( -residual );
				c += dc;
				if ( dc.norm() < 1.0e-9*c.norm() )
					break;
			}
			if ( !solved )
			{
				std::printf( "      %7.3f  did not solve\n", delta );
				continue;
			}
			Step const at = fixture.evaluate( c, held, 1.0, true );
			residual = assemble( at, c );
			Eigen::JacobiSVD<Eigen::MatrixXd> tall( jac );
			Eigen::JacobiSVD<Eigen::MatrixXd> square(
				Eigen::MatrixXd( jac.topRows( degree + 1 ) ) );

			// THE SQUARE SYSTEM WITH A SOURCE SCALE alpha AS ITS LAST UNKNOWN.
			Eigen::MatrixXd withAlpha( degree + 2, degree + 2 );
			withAlpha.leftCols( degree + 1 ) = jac;
			{
				double const step = 1.0e-4;
				Step const up = fixture.evaluate( c, held, 1.0 + step );
				Step const down = fixture.evaluate( c, held, 1.0 - step );
				BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
				withAlpha.col( degree + 1 ) =
					( assemble( up, c ) - assemble( down, c ) )/( 2.0*step );
			}
			Eigen::JacobiSVD<Eigen::MatrixXd> alpha( withAlpha );
			double gEdge = std::sqrt( c( 0 ) );
			double gAxis = std::sqrt( c.sum() );
			std::printf( "      %7.3f %11.3e %11.3e %11.3e %11.3e %11.3e %11.6f %11.6f\n",
			             delta, residual.head( degree + 1 ).norm(),
			             residual( degree + 1 )/weight,
			             tall.singularValues()( degree )/tall.singularValues()( 0 ),
			             square.singularValues()( degree )
			                 /square.singularValues()( 0 ),
			             alpha.singularValues()( degree + 1 )
			                 /alpha.singularValues()( 0 ),
			             gEdge, gAxis );
		}
	}
	std::printf( "    %d solves\n", fixture.solves );
}
