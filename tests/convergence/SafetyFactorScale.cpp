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
 * WHAT FIXES THE SCALE OF A q-DRIVEN EQUILIBRIUM -- MEASUREMENTS.md M-353,
 * and MANTA-EQUILIBRIUM-PLAN.md's D1 (on manta/design), which rests on it.
 *
 * In a cylinder, fixed-boundary equilibria with the same p( Psi_N ) and
 * q( Psi_N ) form a one-parameter family: the toroidal field's scale is free.
 * In a torus the same is VMEC's statement that p( s ), iota( s ) and the
 * boundary need PHIEDGE beside them -- and with q and Psi_N fixed,
 * Phi_edge = 2 pi ( psi_ax - psi_bnd ) int_0^1 q dPsi_N, so the free scalar
 * is psi_ax.
 *
 * MEQ's q-driven solve leaves psi_ax free and fits ALL of g^2, its constant
 * included. This file measures what that does, on SafetyFactorSolver's
 * fixture (g = 2.20 + 0.55 Psi, so g^2 is exactly quadratic):
 *
 *   1. the outer Jacobian G' - I at the answer, psi_ax free, in an
 *      orthonormal basis of the fit's own space -- whether it is near-singular
 *      along the continuum;
 *   2. psi_ax HELD at S = ( 1 + delta ) psi_ax*, the least-squares fixed point
 *      of [ G_S( c ) - c ; w ( psi_ax,located - S )/S ] -- whether a member
 *      of the family exists off the reference, to what residual, and how well
 *      the held problem determines c; a source scale alpha as a square
 *      system's last unknown, which D1 rejects; q's own scale s as that
 *      unknown, q = s q_target, which is a continuum again; and the axis row
 *      in place of one fit row, square;
 *   3. whether the mesh pins the scale.
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

		/// meq::toroidalMapFloorOrthonormal() here: the noise on @ref fitted
		/// in the fit's orthonormal coordinates.
		double floor = 0.0;
	};

	/*
	 * THE DIFFERENCES ARE TAKEN IN THE FIT'S ORTHONORMAL COORDINATES, y = R c,
	 * and the step is the production loop's: the least asked for, raised to
	 * 100 map floors measured in those coordinates (meq::solveForToroidalField
	 * with meq::PinnedToroidalFieldMap::differenceBasis()). There the map's
	 * noise is uniform; in the monomial coefficients the top one's is 1/R_dd
	 * times the rest -- 3.5e-3 at degree 4 in single, where a fixed 1e-3 in c
	 * lifted the free sigma_min this file exists to show is near-null from
	 * 3.6e-4 to 1.5e-3 (MEASUREMENTS.md M-358, M-359).
	 */
	double differenceStep( double least, Step const &at )
	{
		return std::max( least, 100.0*at.floor );
	}

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
			out.floor = meq::toroidalMapFloorOrthonormal( inverted,
				static_cast<unsigned int>( c.size() - 1 ), meq::newtonTolerance( 1.0e-12 ) );
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
			double const step = differenceStep( 1.0e-3, at );
			Eigen::MatrixXd const toC = r.inverse();
			Eigen::MatrixXd alongY( degree + 2, degree + 1 );
			for ( int k = 0; k <= degree; ++k )
			{
				Eigen::VectorXd const plus = c + step*toC.col( k );
				Eigen::VectorXd const minus = c - step*toC.col( k );
				Step const up = fixture.evaluate( plus, held, 1.0 );
				Step const down = fixture.evaluate( minus, held, 1.0 );
				if ( !up.ok || !down.ok )
					return out;
				alongY.col( k ) = ( assemble( up, plus ) - assemble( down, minus ) )
				                  /( 2.0*step );
			}
			// Columns in c, as the tall system's singular values are quoted.
			out.jacobian = alongY*r;
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

	FreeJacobian freeJacobian( Fixture &fixture, Eigen::VectorXd const &c, double least )
	{
		int const degree = static_cast<int>( c.size() ) - 1;
		Eigen::MatrixXd const r = fixture.orthonormal( degree );
		Eigen::MatrixXd const toC = r.inverse();
		Step const here = fixture.evaluate( c, -1.0, 1.0 );
		BOOST_TEST_REQUIRE( here.ok );
		double const step = differenceStep( least, here );
		FreeJacobian out;
		out.inY.resize( degree + 1, degree + 1 );
		Eigen::VectorXd locatedY( degree + 1 );
		for ( int k = 0; k <= degree; ++k )
		{
			Step const up = fixture.evaluate( c + step*toC.col( k ), -1.0, 1.0 );
			Step const down = fixture.evaluate( c - step*toC.col( k ), -1.0, 1.0 );
			BOOST_TEST_REQUIRE( ( up.ok && down.ok ) );
			out.inY.col( k ) = r*( up.fitted - down.fitted )/( 2.0*step );
			out.inY( k, k ) -= 1.0;
			locatedY( k ) = ( up.located - down.located )/( 2.0*step );
		}
		// d psi_ax/dc = R^T d psi_ax/dy.
		out.located = r.transpose()*locatedY;
		return out;
	}
}

/*
 * 1 AND 2 ON ONE FIXTURE, n = 12, AT DEGREES 2 AND 4 -- the suite's version;
 * MEASUREMENTS.md M-353 has the longer study at n = 24 and degrees 2 to 5.
 *
 * psi_ax FREE, as MEQ solves q-driven today: G' - I at the answer has ONE
 * near-null direction, and it is the scale -- psi_ax and g moving together,
 * d psi_ax/psi_ax = ( 1/2 ) d c_0/c_0.
 *
 * psi_ax HELD 10% off the reference: a member of the family exists, to a
 * residual that falls with the degree of g^2 -- the representation pins the
 * scale, not the torus -- and holding psi_ax raises the smallest singular
 * value by orders, as does the axis row in place of the fit's highest row,
 * which is square. Printed and not asserted: alpha, a source scale, and q's
 * own scale s as the square system's last unknown -- each a continuum again.
 * s's column costs no solve: the fitted g^2 is s^2 times s = 1's, and at
 * s = 1 the other columns are the least-squares Jacobian's.
 */
BOOST_AUTO_TEST_CASE( theScaleIsFreeAndHoldingPsiAxisSelectsAMember )
{
	Fixture fixture( 12 );
	double const weight = 10.0;
	double const delta = 0.10;
	double const held = ( 1.0 + delta )*fixture.referenceAxis;
	std::printf( "\n  WHAT FIXES THE SCALE, n = 12 (M-353)\n" );
	std::vector<double> residuals;
	for ( int degree : { 2, 4 } )
	{
		Eigen::MatrixXd const r = fixture.orthonormal( degree );

		// ---- psi_ax free ----------------------------------------------------
		Eigen::VectorXd const answerC = answer( degree );
		FreeJacobian const free = freeJacobian( fixture, answerC, 1.0e-3 );
		Eigen::JacobiSVD<Eigen::MatrixXd> freeSvd( free.inY, Eigen::ComputeFullV );
		Eigen::VectorXd const sigma = freeSvd.singularValues();
		Eigen::VectorXd const nullC = r.inverse()*freeSvd.matrixV().col( degree );
		double const scale = ( free.located.dot( nullC )/fixture.referenceAxis )
		                     /( nullC( 0 )/answerC( 0 ) );
		std::printf( "    degree %d, psi_ax free: sigma_min %.3e, sigma_max %.3e, ratio %.3e, "
		             "d ln psi_ax/d ln c_0 along it %.4f\n", degree, sigma( degree ), sigma( 0 ),
		             sigma( degree )/sigma( 0 ), scale );
		BOOST_TEST( sigma( degree )/sigma( 0 ) < 1.0e-3,
			"at degree " << degree << " the free q-driven Jacobian's smallest singular "
			"value is " << sigma( degree )/sigma( 0 ) << " of its largest: the scale "
			"continuum of fixed p( Psi_N ), q( Psi_N ) equilibria is gone" );
		BOOST_TEST( std::abs( scale - 0.5 ) < 0.1,
			"the near-null direction moves psi_ax by " << scale << " of c_0 relatively, "
			"where the scale continuum, g proportional to psi_ax, moves it by 1/2" );

		// ---- psi_ax held ----------------------------------------------------
		Held const at = solveHeld( fixture, degree, held, weight,
		                           answerC*( 1.0 + delta )*( 1.0 + delta ) );
		BOOST_TEST_REQUIRE( at.solved );
		Eigen::JacobiSVD<Eigen::MatrixXd> tall( at.jacobian );
		double const residual = at.residual.head( degree + 1 ).norm();
		residuals.push_back( residual );

		Eigen::MatrixXd withAlpha( degree + 2, degree + 2 );
		withAlpha.leftCols( degree + 1 ) = at.jacobian;
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

		Eigen::MatrixXd withScale( degree + 2, degree + 2 );
		withScale.leftCols( degree + 1 ) = at.jacobian;
		{
			Step const here = fixture.evaluate( at.c, held, 1.0 );
			BOOST_TEST_REQUIRE( here.ok );
			withScale.col( degree + 1 ).head( degree + 1 ) = r*( 2.0*here.fitted );
			withScale( degree + 1, degree + 1 ) = 0.0;
		}
		Eigen::JacobiSVD<Eigen::MatrixXd> sSvd( withScale );

		double bestDropped = 0.0;
		int bestRow = -1;
		for ( int k = 0; k <= degree; ++k )
		{
			Eigen::MatrixXd dropped( degree + 1, degree + 1 );
			int row = 0;
			for ( int i = 0; i < degree + 2; ++i )
				if ( i != k )
					dropped.row( row++ ) = at.jacobian.row( i );
			Eigen::JacobiSVD<Eigen::MatrixXd> droppedSvd( dropped );
			if ( droppedSvd.singularValues()( degree ) > bestDropped )
			{
				bestDropped = droppedSvd.singularValues()( degree );
				bestRow = k;
			}
		}

		std::printf( "    degree %d, psi_ax held: residual %.3e, axis %.3e; least squares "
		             "s_min %.3e s_max %.3e; axis row for y_%d s_min %.3e; alpha ratio %.3e; "
		             "q's scale s_min %.3e ratio %.3e\n", degree, residual,
		             at.residual( degree + 1 )/weight, tall.singularValues()( degree ),
		             tall.singularValues()( 0 ), bestRow, bestDropped,
		             alpha.singularValues()( degree + 1 )/alpha.singularValues()( 0 ),
		             sSvd.singularValues()( degree + 1 ),
		             sSvd.singularValues()( degree + 1 )/sSvd.singularValues()( 0 ) );
		BOOST_TEST( tall.singularValues()( degree ) > 30.0*sigma( degree ),
			"at degree " << degree << ", holding psi_ax leaves the smallest singular value "
			<< tall.singularValues()( degree ) << " against " << sigma( degree )
			<< " free: holding it no longer determines the scale" );
		BOOST_TEST( bestDropped > 30.0*sigma( degree ),
			"at degree " << degree << ", the axis row in place of fit row y_" << bestRow
			<< " leaves the square system's smallest singular value " << bestDropped
			<< " against " << sigma( degree ) << " free" );
	}
	BOOST_TEST( residuals[ 1 ] < 0.25*residuals[ 0 ],
		"10% along the continuum the held fixed point's residual is " << residuals[ 0 ]
		<< " at degree 2 and " << residuals[ 1 ] << " at degree 4: if it does not fall with "
		"the representation, the torus pins the scale and the continuum is not there" );
	std::printf( "    %d solves\n", fixture.solves );
}

/*
 * 3. AND THE MESH DOES NOT PIN IT: the same residual, 5% along, at n = 12 and
 *    18 (and 24 and 48 in M-353), each against its own mesh's reference.
 */
BOOST_AUTO_TEST_CASE( theMeshDoesNotPinTheScale )
{
	std::printf( "\n  WHAT FIXES THE SCALE, 3: the residual 5%% along the continuum, by resolution "
	             "(M-353)\n" );
	std::vector<double> residuals;
	for ( int cells : { 12, 18 } )
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
