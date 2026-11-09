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

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorTangent.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * WHAT FIXES THE SCALE OF A q-DRIVEN EQUILIBRIUM -- MEASUREMENTS.md M-353,
 * M-QR-6.
 *
 * Fixed-boundary equilibria with the same p( Psi_N ) and q( Psi_N ) form a
 * one-parameter family: the toroidal field's scale is free, VMEC's PHIEDGE.
 * With q matched at every knot, that family is the near-null direction of the
 * q rows' Jacobian -- g^2 and psi_ax moving together, d ln psi_ax = ( 1/2 )
 * d ln g^2 -- and each scale pin is a row that closes it. This file measures
 * both on SafetyFactorSolver's fixture at the closed form g = 2.20 + 0.55 Psi,
 * the Jacobian the exact tangent's (SafetyFactorEquilibrium):
 *
 *   1. the q rows alone: sigma_min/sigma_max, and what the near-null vector
 *      does to psi_ax;
 *   2. with each pin's row, R0B0, ToroidalFlux and PlasmaCurrent, unweighted
 *      and relative: sigma_min raised by orders under every pin -- the pins'
 *      conditioning, which the delivered flux and current rows must not lose
 *      against R0B0's;
 *
 * Measured at n = 8 and 16 (10 and 19 knots): the near-null vector moves
 * psi_ax by 0.487 and 0.484 of g^2, and the pins lift sigma_min from 1.4e-5
 * and 9.5e-3 to 3.8e-2 and 3.4e-2 (R0B0), 1.5e-1 and 1.1e-1 (ToroidalFlux),
 * 1.1e-2 and 1.2e-2 (PlasmaCurrent): the current is the weakest pin, about a
 * third of R0B0, the flux the strongest (M-QR-6). The near-null singular value
 * itself does not fall monotonically with the mesh (1.4e-5, 6.2e-3, 9.5e-3 at
 * n = 8, 12, 16), so nothing is asserted of it but that it is the scale.
 */
namespace
{
	double const g0 = 2.20;
	double const g1 = 0.55;

	double referenceGSquared( double psi )
	{
		return ( g0 + g1*psi )*( g0 + g1*psi );
	}

	struct Measured
	{
		double ratio = 0.0;          // sigma_min/sigma_max of the q rows
		double sigmaMin = 0.0;
		double sigmaNext = 0.0;      // the next one up
		std::vector<double> along;   // | pin . v |/| pin |, v the near-null vector
		double scale = 0.0;          // d ln psi_ax/d ln u along the near-null vector
		std::vector<double> pinned;  // sigma_min with each pin row
	};

	Measured measure( int cells, int order )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, cells );
		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		auto flat = std::make_shared<meq::ConstantProfile const>( 0.0 );
		meq::NormalisedMHDSource source( pPrime, flat, 0.30, 1.0 );
		meq::GradShafranovSolver solver( mesh, order );
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

		std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( mesh ), order );
		meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
		Measured out;
		Eigen::MatrixXd q;
		Eigen::VectorXd nullVector;
		Eigen::Index const m = static_cast<Eigen::Index>( knots );
		for ( meq::ScalePin kind : { meq::ScalePin::EdgeField, meq::ScalePin::ToroidalFlux,
		                             meq::ScalePin::PlasmaCurrent } )
		{
			meq::SafetyFactorEquilibrium::Options options;
			options.family.surfaces = cut.surfaces;
			options.family.innerCut = cut.innerCut;
			options.family.outerCut = cut.outerCut;
			options.pin = kind;
			meq::SafetyFactorEquilibrium equilibrium( solver, source, options );
			meq::ToroidalFieldSpline spline( equilibrium.labels() );
			std::vector<double> values;
			for ( double label : spline.labels() )
				values.push_back( referenceGSquared( 1.0 - label ) );
			spline.setValues( values );
			meq::SafetyFactorEvaluation evaluation;
			BOOST_TEST_REQUIRE( equilibrium.evaluate( spline, evaluation ) );
			meq::SafetyFactorDerivative derivative;
			BOOST_TEST_REQUIRE( equilibrium.differentiate( spline, derivative ) );

			// The q rows, relative: d( q_i/q_i )/d ln u_j.
			q.resize( m, m );
			for ( Eigen::Index i = 0; i < m; ++i )
				for ( Eigen::Index j = 0; j < m; ++j )
				{
					std::size_t const ii = static_cast<std::size_t>( i ), jj = static_cast<std::size_t>( j );
					double d = derivative.qOverG[ ii*knots + jj ]/evaluation.qOverG[ ii ];
					if ( i == j )
						d += 0.5/values[ ii ];
					q( i, j ) = d*values[ jj ];
				}
			Eigen::RowVectorXd pin( m );
			for ( Eigen::Index j = 0; j < m; ++j )
			{
				std::size_t const jj = static_cast<std::size_t>( j );
				pin( j ) = kind == meq::ScalePin::EdgeField
				           ? spline.cardinal( jj, 0.0 )*values[ jj ]/spline( 0.0 )
				           : derivative.measured[ jj ]*values[ jj ]/evaluation.measured;
			}
			if ( kind == meq::ScalePin::EdgeField )
			{
				Eigen::JacobiSVD<Eigen::MatrixXd> svd( q, Eigen::ComputeFullV );
				Eigen::VectorXd const sigma = svd.singularValues();
				out.sigmaMin = sigma( m - 1 );
				out.sigmaNext = sigma( m - 2 );
				nullVector = svd.matrixV().col( m - 1 );
				out.ratio = sigma( m - 1 )/sigma( 0 );
				Eigen::VectorXd const v = svd.matrixV().col( m - 1 );   // in ln u
				double dAxis = 0.0;
				for ( Eigen::Index j = 0; j < m; ++j )
					dAxis += equilibrium.axisDerivative()[ static_cast<std::size_t>( j ) ]
					         *v( j )*values[ static_cast<std::size_t>( j ) ];
				out.scale = ( dAxis/solver.psiAxis() )/v.mean();
			}
			out.along.push_back( std::abs( pin.dot( nullVector ) )/pin.norm() );
			Eigen::MatrixXd closed( m + 1, m );
			closed.topRows( m ) = q;
			closed.row( m ) = pin;
			Eigen::JacobiSVD<Eigen::MatrixXd> svd( closed );
			out.pinned.push_back( svd.singularValues()( m - 1 ) );
		}
		std::printf( "    n = %2d, %zu knots: q rows sigma_min %.3e, next %.3e, ratio %.3e, d ln psi_ax/d ln g^2 "
		             "along it %.4f; sigma_min under R0B0 %.3e, ToroidalFlux %.3e, PlasmaCurrent %.3e; "
		             "|pin.v|/|pin| %.3f %.3f %.3f\n",
		             cells, knots, out.sigmaMin, out.sigmaNext, out.ratio, out.scale, out.pinned[ 0 ], out.pinned[ 1 ],
		             out.pinned[ 2 ], out.along[ 0 ], out.along[ 1 ], out.along[ 2 ] );
		return out;
	}
}

BOOST_AUTO_TEST_CASE( theQRowsLeaveTheScaleFreeAndEachPinFixesIt )
{
	std::printf( "\n  WHAT FIXES THE SCALE (M-353, M-QR-6)\n" );
	char const *names[] = { "R0B0", "ToroidalFlux", "PlasmaCurrent" };
	for ( int cells : { 8, 16 } )
	{
		Measured const at = measure( cells, 2 );
		BOOST_TEST( at.ratio < 1.0e-3,
			"n = " << cells << ": the q rows' smallest singular value is " << at.ratio
			<< " of their largest: the scale continuum of fixed p( Psi_N ), q( Psi_N ) "
			"equilibria is gone, which the torus does not do" );
		BOOST_TEST( std::abs( at.scale - 0.5 ) < 0.1,
			"n = " << cells << ": the near-null direction moves psi_ax by " << at.scale
			<< " of g^2 relatively, where the scale continuum, g proportional to psi_ax, "
			"moves it by 1/2: it is not the scale" );
		for ( std::size_t k = 0; k < 3; ++k )
		{
			BOOST_TEST( at.pinned[ k ] > at.sigmaMin,
				"n = " << cells << ": the " << names[ k ] << " row leaves sigma_min at "
				<< at.pinned[ k ] << ", not above the q rows' " << at.sigmaMin
				<< ": it does not reach the scale" );
			BOOST_TEST( at.pinned[ k ] > 0.1*at.pinned[ 0 ],
				"n = " << cells << ": the " << names[ k ] << " row closes the scale at sigma_min "
				<< at.pinned[ k ] << ", ten times worse than R0B0's " << at.pinned[ 0 ]
				<< ": its conditioning is bad, and the loop under it inherits that" );
		}
	}
}
