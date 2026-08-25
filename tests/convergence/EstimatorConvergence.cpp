#define BOOST_TEST_MODULE EstimatorConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/Estimator.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"

/*
 * The stage-6 acceptance test: the residual error estimator of
 * refs/HDG-GradShafranov-Adaptive.pdf eq (20), measured rather than trusted.
 *
 * "The adaptive loop converges" is weak evidence and this project has twice been
 * bitten by things that converged beautifully to the wrong answer, so the
 * estimator is pinned here against a closed form instead. Three claims, all
 * measured on the same fitted Solov'ev rectangle that
 * tests/convergence/SolovievConvergence.cpp uses:
 *
 *   1. EFFECTIVITY. eta / || q - q_h || is bounded above and below and tends to
 *      a constant. That is what an estimator being efficient and reliable means,
 *      and it is the one property a rate table alone cannot show.
 *   2. EVERY COMPONENT SEPARATELY. eta_1 ... eta_5 each converge at k+1, which
 *      is what GS-2 Table 1 reports ("nearly optimal" for the estimator and all
 *      its separate components). A single total can hide one term being wrong;
 *      five rates cannot -- and one term WAS wrong, see below.
 *   3. THE PSI-STAR CLAIM. eta_2 built on the raw psi_h loses exactly one order,
 *      at every k, which converts the paper's stated reason for post-processing
 *      into a measurement.
 *
 * The denominator in 1 is || q - q_h ||, the L2 error in the flux. That is the
 * quantity a residual estimator for a mixed method controls -- Cockburn and
 * Zhang's analysis, which eq (20) is built from, is an equivalence with the flux
 * error -- and it is also the physically interesting one here, since the magnetic
 * field is built from q. Measured, it is the right choice: the ratio settles to
 * 10.5 at k = 1, 15.3 at k = 2 and 21.9 at k = 3, moving by less than one per
 * cent over the last three meshes at every order. Against || psi - psi_h || it
 * is nearly as flat, and against || psi - psi*_h || it doubles with every
 * refinement, which is what a mismatch of rates looks like: psi* converges at
 * k+2 and the estimator does not.
 *
 * ===================================================================
 * ETA_5 AS PRINTED IN EQ (20) IS WRONG, AND THE WAY IT IS WRONG IS
 * PROVABLE RATHER THAN A MATTER OF TASTE.
 * ===================================================================
 *
 * eq (20)'s fifth term is
 *
 *     eta_5^2 = sum_e h_e^-1 || psihat_h - psi*_h ||_e^2.
 *
 * psihat_h lives in M_h = P_k( e ). psi*_h lives in P_(k+1)( K ), so its trace on
 * e is a polynomial of degree k+1, and NO element of M_h can represent it. Write
 * P_M for the L2( e ) projection onto M_h. Then, by orthogonality,
 *
 *     || psihat_h - psi*_h ||_e^2 = || psihat_h - P_M psi*_h ||_e^2
 *                                 + || ( I - P_M ) psi*_h ||_e^2.
 *
 * The second piece has nothing to do with the error. Put the EXACT solution in
 * for psi*_h and P_M psi for psihat_h -- the best either space can do -- and it
 * is still there, at || ( I - P_M ) psi ||_e = O( h^(k+3/2) ), which the h_e^-1
 * and the O( h^-2 ) edges turn into O( h^k ). So the term as printed does not
 * vanish on the exact solution, and it has an O( h^k ) floor that the rest of the
 * estimator, which converges at k+1, eventually sits underneath.
 *
 * MEASURED, and this is the whole of the case:
 *
 *              eta_5 as printed      eta_5 inside M_h     the floor,
 *                                                       (I-P_M) psi exact
 *     k = 1     1.00 1.00 1.00        1.97 1.99 2.00     1.00 1.00 1.00
 *     k = 2     2.00 2.00 2.00        2.93 2.97 2.98     2.00 2.00 2.00
 *     k = 3     3.00 3.00 3.00        3.98 3.99 4.00     3.00 3.00 3.00
 *
 * The printed form tracks the floor exactly -- at k = 1, h = 0.025, the floor is
 * 2.058e-3 and the printed eta_5 is 2.058e-3, agreeing to four figures -- so it
 * is not merely polluted by the floor, it IS the floor. And it drags the total
 * down with it: eta then converges at 1.73, 1.43, 1.17 at k = 1, heading for 1
 * instead of 2, and the effectivity index grows without bound.
 *
 * Taking the difference inside M_h fixes both, and reproduces the k+1 that GS-2
 * Table 1 reports for eta_5 -- 1.81, 1.5, 1.81, 1.94 at k = 1 and 2.82, 3.93,
 * 3.48, 3.71 at k = 3. So the paper's own numbers agree with the projected form
 * and not with its printed formula. It is also the reading its own prose asks
 * for: eta_5 "estimates the rate of convergence of the hybrid variable and post
 * processed solution -- restricted to the element boundaries -- as approximations
 * to the local trace", and comparing a P_k( e ) function with a degree-(k+1)
 * trace as approximations to the same thing means comparing them in P_k( e ).
 *
 * meq therefore computes the projected form by default and keeps the printed one
 * behind ResidualEstimator::TraceComparison::Literal, so that this stays a
 * measurement. This is the third erratum found in the same pair of papers, after
 * the sign of the Solov'ev source and the Solov'ev coefficients themselves.
 *
 * WHAT THE PRE-MODERNISATION ESTIMATOR DID, for the record, since it is the thing
 * being replaced. It used the raw psi_h everywhere psi*_h belongs (an order in
 * eta_2, measured below); it dropped the quadrature weight from every edge
 * integral, so the edge terms were unweighted point sums; it dropped the Jacobian
 * weight from eta_3 entirely; and it used the unnormalised CalcOrtho normal in
 * eta_3, multiplying that jump by the edge length. None of that is inherited
 * here.
 *
 * WHY THE FITTED RECTANGLE. Because the three claims above are claims about
 * RATES against a closed form, and on the extension path the sequence of
 * computational domains is not a smooth function of h -- stage 5 records the
 * error not even being monotone in the mesh count at k = 3. The estimator is
 * exercised on the curved boundary in tests/convergence/AdaptiveRefinement.cpp,
 * where eta and all five components do converge at k+1 once eta_5 is told to
 * leave Gamma_h alone; the reason it has to be is
 * meq::ResidualEstimator::setTransferredBoundary().
 */

namespace
{

	double const rMin = 0.6;
	double const rMax = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	using Term = meq::ResidualEstimator::Term;
	using Potential = meq::ResidualEstimator::Potential;
	using TraceComparison = meq::ResidualEstimator::TraceComparison;

	int const termCount = meq::ResidualEstimator::termCount;

	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq
			= meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	/// The same rectangle SolovievConvergence.cpp measures on, so that the errors
	/// in the effectivity denominator are the ones already regression-tested
	/// there. Triangles, r well away from zero.
	mfem::Mesh makeMesh( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE, false,
		                                               rMax - rMin, zMax - zMin );
		mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	struct Measurement
	{
		double h;
		int elements;
		double eta;
		double component[ 5 ];
		double etaLiteral;
		double componentLiteral[ 5 ];
		double etaTwoRaw;
		double floorTerm;
		double errorPsi;
		double errorFlux;
		double errorStar;
	};

	double faceLength( mfem::Mesh &mesh, int face )
	{
		mfem::Array<int> vertices;
		mesh.GetFaceVertices( face, vertices );
		return mfem::Distance( mesh.GetVertex( vertices[ 0 ] ),
		                       mesh.GetVertex( vertices[ 1 ] ),
		                       mesh.SpaceDimension() );
	}

	/// sqrt( sum_e h_e^-1 || ( I - P_M ) psi ||_e^2 ) for the EXACT psi, summed
	/// over both sides of every interior edge the way eta_5 is.
	///
	/// This is eta_5 as eq (20) prints it, evaluated at the exact solution with
	/// psihat_h at its best possible value P_M psi: everything an error indicator
	/// is supposed to measure has been removed, and what is left is the term's
	/// floor. Nothing in the solver enters it -- it is a statement about the two
	/// finite element spaces and the exact solution, which is why it can be
	/// asserted as a property of the formula rather than of a computation.
	double traceProjectionFloor( mfem::Mesh &mesh, int order )
	{
		mfem::DG_Interface_FECollection traceColl( order, mesh.Dimension() );
		mfem::FiniteElementSpace traceFes( &mesh, &traceColl );

		meq::analytic::SolovievEquilibrium const &eq = equilibrium();

		double total = 0.0;
		mfem::Vector shape;
		mfem::Vector point;
		mfem::DenseMatrix mass;
		mfem::Vector load;
		mfem::Vector coefficients;

		for ( int f = 0; f < mesh.GetNumFaces(); ++f )
		{
			mfem::FaceElementTransformations *ftr = mesh.GetFaceElementTransformations( f );
			if ( !ftr )
				continue;

			mfem::FiniteElement const *fe = traceFes.GetFaceElement( f );
			int const dof = fe->GetDof();
			shape.SetSize( dof );
			mass.SetSize( dof );
			mass = 0.0;
			load.SetSize( dof );
			load = 0.0;

			mfem::IntegrationRule const &ir =
				mfem::IntRules.Get( ftr->GetGeometryType(), 2*( order + 1 ) + 6 );

			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				ftr->SetAllIntPoints( &ip );
				ftr->Transform( ip, point );
				double const weight = ip.weight*ftr->Weight();
				double const exact = eq.psi( point( 0 ), point( 1 ) );

				fe->CalcShape( ip, shape );
				mfem::AddMult_a_VVt( weight, shape, mass );
				load.Add( weight*exact, shape );
			}

			mfem::DenseMatrixInverse inverse( mass );
			coefficients.SetSize( dof );
			inverse.Mult( load, coefficients );

			// A second pass, integrating ( psi - P_M psi )^2 directly. The
			// algebraic shortcut || psi ||^2 - load . ( M^-1 load ) is the same
			// number and is not usable: at k = 3 on the finest mesh the remainder
			// is 1e-16 of || psi ||^2, so the subtraction loses every significant
			// figure and the floor comes out 2.6 times too large with a rate of
			// 1.6 instead of 3. Measured, and it is exactly the kind of thing a
			// single-mesh check would not have caught.
			double remainder = 0.0;
			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				ftr->SetAllIntPoints( &ip );
				ftr->Transform( ip, point );
				double const weight = ip.weight*ftr->Weight();
				double const exact = eq.psi( point( 0 ), point( 1 ) );

				fe->CalcShape( ip, shape );
				double const projected = shape*coefficients;
				remainder += weight*( exact - projected )*( exact - projected );
			}

			double const he = faceLength( mesh, f );
			int const sides = ( ftr->Elem2No >= 0 ) ? 2 : 1;
			total += sides*remainder/he;
		}
		return std::sqrt( total );
	}

	Measurement measure( int order, int n )
	{
		meq::analytic::SolovievEquilibrium const &eq = equilibrium();
		meq::SolovievSource const source( eq.getA() );

		mfem::Mesh mesh = makeMesh( n );

		mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.f( x( 0 ), x( 1 ), 0.0 );
		} );
		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		solver.solve();
		solver.postProcess();

		// The meq::Source overload, so that eta_1 evaluates F at psi* rather than
		// at a frozen coefficient. For Solov'ev the two agree to round-off, which
		// is exactly why this benchmark can be used to check the semi-linear path
		// of the estimator against the linear one.
		meq::ResidualEstimator estimator( solver, source );

		Measurement point;
		point.h = ( rMax - rMin )/static_cast<double>( n );
		point.elements = mesh.GetNE();

		point.eta = estimator.GetTotalError();
		for ( int t = 0; t < termCount; ++t )
			point.component[ t ] = estimator.component( static_cast<Term>( t ) );

		estimator.setTraceComparison( TraceComparison::Literal );
		point.etaLiteral = estimator.GetTotalError();
		for ( int t = 0; t < termCount; ++t )
			point.componentLiteral[ t ] = estimator.component( static_cast<Term>( t ) );
		estimator.setTraceComparison( TraceComparison::Projected );

		estimator.setPotential( Potential::Raw );
		point.etaTwoRaw = estimator.component( Term::Constitutive );
		estimator.setPotential( Potential::PostProcessed );

		point.floorTerm = traceProjectionFloor( mesh, order );

		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.errorStar = solver.postProcessedPotentialError( psiCoeff );
		return point;
	}

	double rate( double coarseError, double fineError, double refinementRatio )
	{
		return std::log( coarseError/fineError )/std::log( refinementRatio );
	}

	/// The same four dyadic meshes as SolovievConvergence.cpp, so three measured
	/// rates per quantity per order.
	std::vector<int> const meshSizes = { 4, 8, 16, 32 };

	/// k+1, less the slack allowed for a two-mesh rate estimate. Larger than
	/// SolovievConvergence.cpp's 0.15 because eta_3 and eta_4 approach their rate
	/// from below on this sequence -- the tightest measured is eta_4 at k = 2,
	/// h = 0.1, which comes out at 2.83.
	double const rateSlack = 0.25;

	std::vector<Measurement> study( int order )
	{
		std::vector<Measurement> points;
		points.reserve( meshSizes.size() );
		for ( int n : meshSizes )
			points.push_back( measure( order, n ) );
		return points;
	}

	void printTable( int order, std::vector<Measurement> const &points )
	{
		std::printf( "\n  eq (20) on the fitted Solov'ev rectangle, k = %d\n", order );
		std::printf( "  %8s %7s %11s %5s", "h", "elem", "eta", "rate" );
		for ( int t = 0; t < termCount; ++t )
			std::printf( " %11s %5s", meq::ResidualEstimator::name( static_cast<Term>( t ) ),
			             "rate" );
		std::printf( "\n" );

		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			double const ratio = ( i == 0 ) ? 0.0 : points[ i - 1 ].h/p.h;
			std::printf( "  %8.5f %7d %11.4e", p.h, p.elements, p.eta );
			if ( i == 0 )
				std::printf( " %5s", "-" );
			else
				std::printf( " %5.2f", rate( points[ i - 1 ].eta, p.eta, ratio ) );
			for ( int t = 0; t < termCount; ++t )
			{
				std::printf( " %11.4e", p.component[ t ] );
				if ( i == 0 )
					std::printf( " %5s", "-" );
				else
					std::printf( " %5.2f", rate( points[ i - 1 ].component[ t ],
					                             p.component[ t ], ratio ) );
			}
			std::printf( "\n" );
		}

		std::printf( "    effectivity, eta / L2 error:\n" );
		std::printf( "    %8s %10s %10s %10s\n", "h", "/L2(q)", "/L2(psi)", "/L2(psi*)" );
		for ( Measurement const &p : points )
			std::printf( "    %8.5f %10.3f %10.3f %10.3e\n",
			             p.h, p.eta/p.errorFlux, p.eta/p.errorPsi, p.eta/p.errorStar );
		std::fflush( stdout );
	}

}

/// Every component of eq (20) at k+1, and the total with them. GS-2 Table 1
/// reports "nearly optimal" rates for the estimator and all its separate
/// components; this is that claim, one rate at a time, because a single total
/// can hide one wrong term and five rates cannot.
BOOST_AUTO_TEST_CASE( everyComponentConvergesAtKPlusOne )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order );
		printTable( order, points );

		double const expected = order + 1.0 - rateSlack;

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;

			double const total = rate( points[ i - 1 ].eta, points[ i ].eta, ratio );
			BOOST_TEST( total >= expected,
			            "k = " << order << ", h = " << points[ i ].h << ": eta converged at "
			            << total << ", wanted " << expected );

			for ( int t = 0; t < termCount; ++t )
			{
				double const measured = rate( points[ i - 1 ].component[ t ],
				                              points[ i ].component[ t ], ratio );
				BOOST_TEST( measured >= expected,
				            "k = " << order << ", h = " << points[ i ].h << ": "
				            << meq::ResidualEstimator::name( static_cast<Term>( t ) )
				            << " converged at " << measured << ", wanted " << expected );
			}
		}

		// A rate is blind to an estimator that is uniformly too large or too
		// small, so the absolute value is pinned too. The ceilings sit at about
		// three times the measured finest-mesh eta: 6.391e-04, 3.820e-06,
		// 1.927e-08 at k = 1, 2, 3.
		double const ceiling[ 3 ] = { 2.0e-3, 1.2e-5, 6.0e-8 };
		BOOST_TEST( points.back().eta < ceiling[ order - 1 ],
		            "k = " << order << ": eta is " << points.back().eta
		            << " on the finest mesh, above the ceiling " << ceiling[ order - 1 ] );
	}
}

/// The effectivity index. An estimator is useful when eta / || error || is
/// bounded above (reliable) and below (efficient) and does not drift, which is a
/// stronger statement than "eta converges at the right rate" -- two quantities
/// can share a rate and still differ by a factor that grows with k or with the
/// mesh. The denominator is the L2 error in the flux; see the file comment.
BOOST_AUTO_TEST_CASE( theEffectivityIndexIsBoundedAndSettles )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order );

		std::printf( "\n  effectivity, k = %d\n", order );
		std::vector<double> index;
		for ( Measurement const &p : points )
		{
			index.push_back( p.eta/p.errorFlux );
			std::printf( "    h = %8.5f   eta = %11.4e   L2(q) = %11.4e   ratio = %8.4f\n",
			             p.h, p.eta, p.errorFlux, index.back() );
		}
		std::fflush( stdout );

		for ( double value : index )
		{
			BOOST_TEST( value > 1.0,
			            "k = " << order << ": the effectivity index is " << value
			            << ", so eta is smaller than the error it is supposed to bound" );
			BOOST_TEST( value < 60.0,
			            "k = " << order << ": the effectivity index is " << value
			            << ", so eta overestimates the flux error by more than the "
			            "measured 10-22" );
		}

		// And it settles. Measured, the ratio between consecutive indices is
		// within 1.5 per cent of one over the last three meshes at every order --
		// 10.72 -> 10.61 -> 10.55 -> 10.52 at k = 1. A drifting index is what a
		// term of the wrong order looks like, and it is exactly what the printed
		// eta_5 produces; see theLiteralTraceTermIsItsOwnFloor below.
		for ( std::size_t i = 2; i < index.size(); ++i )
		{
			double const drift = index[ i ]/index[ i - 1 ];
			BOOST_TEST( ( drift > 0.95 && drift < 1.05 ),
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the effectivity index moved by a factor " << drift
			            << " on this refinement, so it is not settling" );
		}
	}
}

/// eta_2 on the raw psi_h loses exactly one order, at every k. This is the
/// paper's stated reason for post-processing at all, and CLAUDE.md's reason for
/// bringing the dropped stage-3 work back for stage 6, turned into a
/// measurement. It is cheap, and it is the difference between reproducing the
/// published estimator and reproducing meq's own degraded copy of it.
BOOST_AUTO_TEST_CASE( theRawPotentialCostsEtaTwoExactlyOneOrder )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order );

		std::printf( "\n  eta_2 on psi* against psi_h, k = %d\n", order );
		std::printf( "    %8s %12s %6s %12s %6s\n",
		             "h", "eta_2(psi*)", "rate", "eta_2(psi_h)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			std::printf( "    %8.5f %12.4e", p.h, p.component[ 1 ] );
			if ( i == 0 )
			{
				std::printf( " %6s %12.4e %6s\n", "-", p.etaTwoRaw, "-" );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( " %6.3f %12.4e %6.3f\n",
				             rate( points[ i - 1 ].component[ 1 ], p.component[ 1 ], ratio ),
				             p.etaTwoRaw,
				             rate( points[ i - 1 ].etaTwoRaw, p.etaTwoRaw, ratio ) );
			}
		}
		std::fflush( stdout );

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const withStar = rate( points[ i - 1 ].component[ 1 ],
			                              points[ i ].component[ 1 ], ratio );
			double const withRaw = rate( points[ i - 1 ].etaTwoRaw,
			                             points[ i ].etaTwoRaw, ratio );

			BOOST_TEST( withStar >= order + 1.0 - rateSlack,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": eta_2 on psi* converged at " << withStar );
			BOOST_TEST( withRaw < order + 0.5,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": eta_2 on the raw psi_h converged at " << withRaw
			            << ", which is better than the reduced order the paper predicts "
			            "-- so either psi* is not doing what it is here for, or the two "
			            "branches are not measuring the same term" );
			BOOST_TEST( withRaw > order - 0.5,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": eta_2 on the raw psi_h converged at " << withRaw
			            << ", worse than the one order the paper predicts losing" );
		}

		// And it is not a small effect: at k = 1 the raw form is 124 times larger
		// on the finest mesh, at k = 2 it is 178 times, at k = 3 it is 407 times.
		BOOST_TEST( points.back().etaTwoRaw > 50.0*points.back().component[ 1 ],
		            "k = " << order << ": eta_2 is " << points.back().etaTwoRaw
		            << " on psi_h against " << points.back().component[ 1 ]
		            << " on psi*, which is not the difference an order should make "
		            "over four refinements" );
	}
}

/// eta_5 as eq (20) prints it IS its own floor, and the floor carries no
/// information about the error.
///
/// The floor is computed with no solver in it at all: the exact psi, the best
/// psihat_h that M_h admits, and eq (20)'s own weights. If the printed formula
/// were a residual it would be zero. It is not, it converges at k rather than
/// k+1, and it agrees with the measured printed eta_5 to three or four
/// significant figures -- so the printed term is not a residual with a floor
/// under it, it is the floor.
BOOST_AUTO_TEST_CASE( theLiteralTraceTermIsItsOwnFloor )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order );

		std::printf( "\n  eta_5, printed against projected, k = %d\n", order );
		std::printf( "    %8s %12s %6s %12s %6s %12s %6s %9s\n",
		             "h", "projected", "rate", "printed", "rate", "floor", "rate",
		             "print/flr" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			std::printf( "    %8.5f %12.4e", p.h, p.component[ 4 ] );
			if ( i == 0 )
			{
				std::printf( " %6s %12.4e %6s %12.4e %6s %9.5f\n",
				             "-", p.componentLiteral[ 4 ], "-", p.floorTerm, "-",
				             p.componentLiteral[ 4 ]/p.floorTerm );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( " %6.3f %12.4e %6.3f %12.4e %6.3f %9.5f\n",
				             rate( points[ i - 1 ].component[ 4 ], p.component[ 4 ], ratio ),
				             p.componentLiteral[ 4 ],
				             rate( points[ i - 1 ].componentLiteral[ 4 ],
				                   p.componentLiteral[ 4 ], ratio ),
				             p.floorTerm,
				             rate( points[ i - 1 ].floorTerm, p.floorTerm, ratio ),
				             p.componentLiteral[ 4 ]/p.floorTerm );
			}
		}
		std::fflush( stdout );

		for ( Measurement const &p : points )
		{
			// The floor is not zero, which is the whole argument: eq (20)'s fifth
			// term does not vanish on the exact solution.
			BOOST_TEST( p.floorTerm > 0.0,
			            "k = " << order << ", h = " << p.h
			            << ": the projection floor came out zero, so either the trace "
			            "space represents a degree-(k+1) trace exactly or the "
			            "measurement is wrong" );

			// And it is the printed term, to within a few per cent.
			double const share = p.componentLiteral[ 4 ]/p.floorTerm;
			BOOST_TEST( ( share > 0.9 && share < 1.1 ),
			            "k = " << order << ", h = " << p.h << ": the printed eta_5 is "
			            << p.componentLiteral[ 4 ] << " and its solution-independent floor "
			            << p.floorTerm << ", a ratio of " << share
			            << " -- if that has moved away from one, the printed term has "
			            "stopped being dominated by the floor and this argument needs "
			            "revisiting" );
		}

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const floorRate = rate( points[ i - 1 ].floorTerm,
			                               points[ i ].floorTerm, ratio );
			double const literalRate = rate( points[ i - 1 ].componentLiteral[ 4 ],
			                                 points[ i ].componentLiteral[ 4 ], ratio );
			double const projectedRate = rate( points[ i - 1 ].component[ 4 ],
			                                   points[ i ].component[ 4 ], ratio );

			BOOST_TEST( std::abs( floorRate - order ) < 0.25,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the projection floor converged at " << floorRate
			            << ", not the k it should" );
			BOOST_TEST( std::abs( literalRate - order ) < 0.25,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the printed eta_5 converged at " << literalRate
			            << ", not the k the floor forces on it" );
			BOOST_TEST( projectedRate >= order + 1.0 - rateSlack,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the projected eta_5 converged at " << projectedRate
			            << ", wanted " << order + 1.0 - rateSlack );
		}

		// The consequence for the estimator as a whole: with the printed eta_5
		// the total loses an order, which is what makes this worth fixing rather
		// than noting.
		double const span = points.front().h/points.back().h;
		double const printedTotal = rate( points.front().etaLiteral,
		                                  points.back().etaLiteral, span );
		double const projectedTotal = rate( points.front().eta, points.back().eta, span );
		std::printf( "    across the sequence: eta printed %.3f, eta projected %.3f\n",
		             printedTotal, projectedTotal );
		std::fflush( stdout );

		// Measured across the sequence: 1.99 against 1.44 at k = 1, 2.99 against
		// 2.32 at k = 2, 3.97 against 3.55 at k = 3. The gap narrows with k on a
		// fixed mesh sequence rather than widening, and that is not a weakening of
		// the argument: the printed eta_5 is O( h^k ) while every other term is
		// O( h^(k+1) ), so it always wins eventually, but at k = 3 eta_1 starts
		// out fifty times larger and four dyadic meshes are not enough for the
		// crossover. Chasing it would need finer meshes than round-off allows.
		BOOST_TEST( projectedTotal - printedTotal > 0.3,
		            "k = " << order << ": eta converges at " << printedTotal
		            << " with the printed eta_5 and " << projectedTotal
		            << " with the projected one, so the printed form is no longer "
		            "costing the estimator an order and this test has stopped saying "
		            "anything" );
	}
}

/// The post-processing itself: DarcyForm::Reconstruct() gives psi* at k+2, so
/// nothing had to be written by hand for stage 6. Measured here rather than in
/// SolovievConvergence.cpp because it is the estimator that needs it, and because
/// a psi* that had quietly stopped superconverging would show up first as eta_2
/// losing an order -- which is the same symptom as using psi_h, and worth being
/// able to tell apart.
BOOST_AUTO_TEST_CASE( reconstructGivesThePostProcessedPotentialAtKPlusTwo )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order );

		std::printf( "\n  psi* against psi_h, k = %d\n", order );
		std::printf( "    %8s %12s %6s %12s %6s\n",
		             "h", "L2(psi)", "rate", "L2(psi*)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			std::printf( "    %8.5f %12.4e", p.h, p.errorPsi );
			if ( i == 0 )
			{
				std::printf( " %6s %12.4e %6s\n", "-", p.errorStar, "-" );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( " %6.3f %12.4e %6.3f\n",
				             rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorStar,
				             rate( points[ i - 1 ].errorStar, p.errorStar, ratio ) );
			}
		}
		std::fflush( stdout );

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const measured = rate( points[ i - 1 ].errorStar,
			                              points[ i ].errorStar, ratio );
			BOOST_TEST( measured >= order + 2.0 - rateSlack,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi* converged at " << measured << ", wanted "
			            << order + 2.0 - rateSlack
			            << " -- if this has dropped to k+1 then Reconstruct() is no "
			            "longer post-processing and eta_2 cannot be right either" );
		}

		// Measured: 3.495e-07, 9.302e-10, 1.555e-12 on the finest mesh at
		// k = 1, 2, 3, against 3.595e-05, 1.913e-07 and 5.493e-10 for psi_h.
		BOOST_TEST( points.back().errorStar < 0.1*points.back().errorPsi,
		            "k = " << order << ": psi* is only " << points.back().errorStar
		            << " against " << points.back().errorPsi << " for psi_h, which is "
		            "not the gain an extra order should give over four refinements" );
	}
}

/// The two marking strategies, on a made-up indicator where the answer can be
/// worked out by hand. Doerfler is what the convergence proof assumes; maximum is
/// what GS-2's own experiments used, at gamma = 0.3.
BOOST_AUTO_TEST_CASE( theMarkingStrategiesMarkWhatTheyClaimTo )
{
	// Squares 16, 9, 4, 1, summing to 30.
	mfem::Vector indicator( 4 );
	indicator( 0 ) = 1.0;
	indicator( 1 ) = 4.0;
	indicator( 2 ) = 2.0;
	indicator( 3 ) = 3.0;

	mfem::Array<int> marked;

	// Doerfler at gamma = 0.5 wants 15 of 30: 16 from the largest alone does it.
	meq::markDoerfler( indicator, 0.5, marked );
	BOOST_TEST( marked.Size() == 1, "Doerfler at 0.5 marked " << marked.Size()
	            << " elements, not the one whose square is already past half the total" );
	BOOST_TEST( marked[ 0 ] == 1, "Doerfler marked element " << marked[ 0 ]
	            << ", not the largest" );

	// At 0.9 it wants 27: 16 + 9 = 25 is not enough, 16 + 9 + 4 = 29 is.
	meq::markDoerfler( indicator, 0.9, marked );
	BOOST_TEST( marked.Size() == 3, "Doerfler at 0.9 marked " << marked.Size()
	            << " elements, not the three needed to reach 27 of 30" );

	// At 1.0 it must take everything, and a minimal set is still minimal: the
	// loop must not stop one short of the total.
	meq::markDoerfler( indicator, 1.0, marked );
	BOOST_TEST( marked.Size() == 4, "Doerfler at 1.0 marked " << marked.Size()
	            << " of 4 elements" );

	// Maximum at 0.5 takes eta_K >= 2, which is three of the four.
	meq::markMaximum( indicator, 0.5, marked );
	BOOST_TEST( marked.Size() == 3, "maximum marking at 0.5 marked " << marked.Size()
	            << " elements, not the three at or above half of 4" );

	// GS-2's own choice. At 0.3 the threshold is 1.2, so the element at 1 is out.
	meq::markMaximum( indicator, 0.3, marked );
	BOOST_TEST( marked.Size() == 3, "maximum marking at 0.3 marked " << marked.Size()
	            << " elements" );

	// And at 1.0 only the maximum itself.
	meq::markMaximum( indicator, 1.0, marked );
	BOOST_TEST( marked.Size() == 1, "maximum marking at 1.0 marked " << marked.Size()
	            << " elements, not just the largest" );

	// A zero indicator marks nothing rather than everything. GS-2 section 3.3
	// names exactly this as the condition under which the computational domains
	// stop exhausting Omega, so it is a state to report rather than to paper over.
	mfem::Vector zero( 4 );
	zero = 0.0;
	meq::markDoerfler( zero, 0.5, marked );
	BOOST_TEST( marked.Size() == 0, "Doerfler marked " << marked.Size()
	            << " elements on a zero indicator" );
	meq::markMaximum( zero, 0.3, marked );
	BOOST_TEST( marked.Size() == 0, "maximum marking marked " << marked.Size()
	            << " elements on a zero indicator" );
}
