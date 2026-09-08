#define BOOST_TEST_MODULE LimiterCurve
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * THE LIMITER AS A CURVE, RESTRICTED TO THE POLYGON THE MESH IS FITTED TO.
 *
 * FB-6's remaining capability gap, and FREE-BOUNDARY-PLAN.md section 7.20 is
 * why it is one. `psi_bnd` is the flux at the limiter CONTACT, and until now
 * MEQ had to be TOLD where that contact is: setBoundaryFluxPoint() takes a
 * point. The physical problem is `psi_bnd = max psi` over the limiter CURVE,
 * with the contact an output of the solve rather than an input to it -- which
 * is what every production free-boundary code does, and what takes another
 * code's grid artefact out of MEQ's comparison against it. Section 7.20
 * measures the artefact: `freegs4e`'s 513^2 run puts its contact 0.6 m round
 * the same limiter circle from where its own 129^2 run puts it, because a
 * contact read off a grid is a maximum over CELLS.
 *
 * THE RESTRICTION IS THE POLYGON AND IT COSTS ALMOST NOTHING, WHICH IS THE
 * POINT. A limiter is PRESCRIBED INPUT -- it does not move with the solution --
 * so the mesh can be fitted to it, exactly as section 7.9's conductors are.
 * Fitted, the limiter IS a union of mesh faces, and then:
 *
 *   - there is no cut element and no cut quadrature anywhere in it;
 *   - there is no separate limiter geometry to keep consistent with the mesh;
 *   - the contact lies on a face, so "which element" has a two-valued answer
 *     with a KNOWN difference -- the face jump, O( h^{k+1} ) -- rather than an
 *     ambiguous one.
 *
 * What it gives up is that the polygon INSCRIBES whatever smooth curve it was
 * fitted to. On a gmsh mesh built by tools/mesh/halfdisc.py --limiter that is
 * all it gives up: measured, the polygon's vertices sit on the true circle to
 * 3.3e-16 and its length falls short only by the chords.
 *
 * THIS FILE PAINTS ITS POLYGON ONTO A CARTESIAN MESH RATHER THAN MESHING ONE,
 * AND THAT IS DELIBERATE. Any set of elements defines a fitted polygon -- its
 * own boundary -- so the constraint can be exercised with no external mesher
 * in the loop, which keeps this a test of the CONSTRAINT rather than of gmsh.
 * The staircase that painting gives is a worse limiter than a fragmented
 * circle and a perfectly good polygon, and every property asserted below is a
 * property of the polygon rather than of the curve behind it.
 */
namespace
{
	using meq::GradShafranovSolver;

	/*
	 * THE LIMITER IS A RECTANGLE ON MESH LINES, AND THAT IS THE WHOLE OF THE
	 * RESTRICTION.
	 *
	 * standardBox() is [ 0.6, 1.4 ] x [ -0.6, 0.6 ] with n cells a side, so
	 * these four edges land on mesh lines at n = 12, 24 and 48 exactly -- the
	 * offsets are 3/12, 9/12, 3/12, 9/12 of the way across and every one of
	 * them is an integer number of cells at all three. So the polygon is the
	 * SAME CURVE at every mesh, its perimeter is exactly 2.0 at every mesh, and
	 * there is no geometric error between the limiter and the mesh at all.
	 *
	 * A CIRCLE PAINTED THE SAME WAY WOULD NOT BE, AND MEASURING IT IS WHAT
	 * SETTLED THIS. Painting a disc of radius 0.30 by centroid gives a closed
	 * polygon at every mesh -- and its perimeter reads 2.828, 2.368, 2.611
	 * against the circle's 1.885 at n = 12, 24, 48, converging to nothing. That
	 * is the staircase: the polygon converges to the circle as a SET, at
	 * O( h ) in Hausdorff distance, while its LENGTH converges to the wrong
	 * number, and `max psi` over it would then converge at O( h ) rather than
	 * at the field's own order. Fitting a mesh to a curve is what
	 * tools/mesh/halfdisc.py --limiter does instead, and it is a stronger thing
	 * than painting: fragmenting the circle into the geometry puts the polygon's
	 * VERTICES on the true circle -- measured, to 3.3e-16 -- so the polygon
	 * inscribes it at O( h^2 ). Painting is fitting the polygon to the mesh;
	 * fragmenting is fitting the mesh to the curve, and only the second
	 * converges at the rate a limiter deserves.
	 *
	 * Choosing a polygon that is exactly representable takes that whole question
	 * out of this file, which is about the CONSTRAINT.
	 */
	double const limiterRMin = 0.8;
	double const limiterRMax = 1.2;
	double const limiterZMin = -0.3;
	double const limiterZMax = 0.3;

	double limiterPerimeter()
	{
		return 2.0*( ( limiterRMax - limiterRMin )
		             + ( limiterZMax - limiterZMin ) );
	}

	int const limiterAttribute = 20;

	/// Paint the region the limiter encloses. Every element whose centroid is
	/// inside the rectangle takes `limiterAttribute`; the polygon is then that
	/// region's own boundary, and because the rectangle is on mesh lines it is
	/// the rectangle exactly rather than a staircase approximating it.
	int paintLimiter( mfem::Mesh &mesh )
	{
		int painted = 0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::Vector centre( 2 );
			mesh.GetElementCenter( e, centre );
			bool const inside = centre( 0 ) > limiterRMin
			                    && centre( 0 ) < limiterRMax
			                    && centre( 1 ) > limiterZMin
			                    && centre( 1 ) < limiterZMax;
			mesh.SetAttribute( e, inside ? limiterAttribute : 1 );
			if ( inside )
				++painted;
		}
		mesh.SetAttributes();
		return painted;
	}

	/// The faces of the painted polygon, as ( face, element-on-the-inside )
	/// pairs -- the same rule the solver uses, re-derived here rather than
	/// shared, so that this file checks the solver rather than agreeing with
	/// it by construction.
	struct PolygonFace
	{
		int face;
		int element;
		bool second;
	};

	std::vector<PolygonFace> polygonFaces( mfem::Mesh &mesh )
	{
		std::vector<PolygonFace> out;
		for ( int f = 0; f < mesh.GetNumFaces(); ++f )
		{
			int first = -1, second = -1;
			mesh.GetFaceElements( f, &first, &second );
			if ( first < 0 || second < 0 )
				continue;
			bool const a = mesh.GetAttribute( first ) == limiterAttribute;
			bool const b = mesh.GetAttribute( second ) == limiterAttribute;
			if ( a == b )
				continue;
			out.push_back( { f, a ? first : second, !a } );
		}
		return out;
	}

	/// `max psi_h` over the polygon by DENSE SAMPLING, and where it is.
	///
	/// The independent route: the solver closes on the contact by a bracketed
	/// golden section per face, this walks a uniform grid over every face. Two
	/// different searches over the same set, which is what makes the comparison
	/// below a statement about the answer rather than about the algorithm.
	double sampledMaximum( mfem::Mesh &mesh, mfem::GridFunction const &potential,
	                       int samplesPerFace, double &r, double &z )
	{
		double best = -std::numeric_limits<double>::infinity();
		for ( PolygonFace const &entry : polygonFaces( mesh ) )
		{
			mfem::FaceElementTransformations faceScratch;
			mfem::IsoparametricTransformation e1, e2;
			mesh.GetFaceElementTransformations( entry.face, faceScratch, e1, e2 );
			if ( faceScratch.GetGeometryType() == mfem::Geometry::INVALID )
				continue;

			for ( int i = 0; i < samplesPerFace; ++i )
			{
				double const t = static_cast<double>( i )/( samplesPerFace - 1 );
				mfem::IntegrationPoint ip;
				ip.Set1w( t, 1.0 );
				faceScratch.SetAllIntPoints( &ip );

				mfem::IntegrationPoint const &inside =
					entry.second ? faceScratch.GetElement2IntPoint()
					             : faceScratch.GetElement1IntPoint();
				double const value = potential.GetValue( entry.element, inside );
				if ( value > best )
				{
					best = value;
					mfem::Vector point( 2 );
					faceScratch.Transform( ip, point );
					r = point( 0 );
					z = point( 1 );
				}
			}
		}
		return best;
	}

	/// The guess of HighBetaConvergence's FB-3 case, which is where this
	/// fixture's profiles and amplitude come from as well.
	mfem::FunctionCoefficient bump( double height )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		double const rMin = box.rMin;
		double const zMin = box.zMin;
		double const width = box.width();
		double const depth = box.height();
		return mfem::FunctionCoefficient(
			[ height, rMin, zMin, width, depth ]( mfem::Vector const &x )
			{
				return height*std::sin( M_PI*( x( 0 ) - rMin )/width )
				       *std::sin( M_PI*( x( 1 ) - zMin )/depth );
			} );
	}
}

/*
 * THE POLYGON IS RECOVERED FROM THE ELEMENT ATTRIBUTES, AND IT IS CLOSED.
 *
 * The cheapest possible statement, and it guards the one thing that would
 * otherwise fail silently: an empty polygon. `max` over an empty set is minus
 * infinity and its border row is all zeroes, which does NOT diverge -- it
 * converges, to a psi_bnd pinned by nothing at all.
 *
 * Closedness is asserted by counting: a closed polygon of straight segments
 * has every vertex in exactly two of them. A rule that picked up a stray face
 * -- the failure halfdisc.py's own docstring records, where conductor outlines
 * joined Gamma -- would leave a vertex with one or three.
 */
BOOST_AUTO_TEST_CASE( theLimiterPolygonIsClosedAndFittedToTheMesh )
{
	meq::tests::Rectangle const box = meq::tests::standardBox();

	std::printf( "\n  THE PAINTED LIMITER POLYGON\n" );
	std::printf( "    %5s %9s %9s %9s %14s %14s\n",
	             "n", "enclosed", "faces", "odd ends", "length", "exact" );

	for ( int n : { 12, 24, 48 } )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );
		int const painted = paintLimiter( mesh );
		std::vector<PolygonFace> const faces = polygonFaces( mesh );

		BOOST_TEST_REQUIRE( painted > 0,
			"nothing was painted at n = " << n );
		BOOST_TEST_REQUIRE( !faces.empty(),
			"the polygon is empty at n = " << n << ", which is the failure that "
			"does not diverge: max over an empty set is -infinity and its "
			"border row is zero" );

		// Every vertex in exactly two segments.
		std::map<int, int> uses;
		double length = 0.0;
		for ( PolygonFace const &entry : faces )
		{
			mfem::Array<int> v;
			mesh.GetFaceVertices( entry.face, v );
			BOOST_TEST_REQUIRE( v.Size() == 2,
				"a face of a 2-D mesh has " << v.Size() << " vertices" );
			for ( int i = 0; i < 2; ++i )
				++uses[ v[ i ] ];
			double const *a = mesh.GetVertex( v[ 0 ] );
			double const *b = mesh.GetVertex( v[ 1 ] );
			length += std::hypot( b[ 0 ] - a[ 0 ], b[ 1 ] - a[ 1 ] );
		}

		int odd = 0;
		for ( auto const &use : uses )
			if ( use.second != 2 )
				++odd;

		std::printf( "    %5d %9d %9zu %9d %14.6f %14.6f\n", n, painted,
		             faces.size(), odd, length, limiterPerimeter() );

		// EXACTLY the prescribed rectangle, at every mesh, because its edges
		// are mesh lines. A staircase would read something larger and would not
		// settle -- see the note beside limiterRMin for what a painted circle
		// does instead.
		BOOST_TEST( std::fabs( length - limiterPerimeter() ) < 1.0e-12,
			"the polygon's length at n = " << n << " is " << length
			<< " against the prescribed rectangle's " << limiterPerimeter()
			<< ". The rectangle is on mesh lines, so the two must agree to "
			"round-off; a gap means the paint is not landing on the cells the "
			"edges were chosen for." );

		BOOST_TEST( odd == 0,
			"the polygon at n = " << n << " has " << odd << " vertices that are "
			"not in exactly two segments, so it is not a closed curve. A stray "
			"face is halfdisc.py's own recorded failure, one dimension down." );
	}
}

/*
 * THE CONTACT IS FOUND WHERE THE FIELD IS LARGEST ON THE POLYGON.
 *
 * The property the constraint is NAMED for, checked by a second search over
 * the same set: the solver closes on the contact by a bracketed golden section
 * per face, this walks a uniform grid over every face and takes the best. Two
 * different searches, so agreement is a statement about the answer rather than
 * about the algorithm -- and a dense sample cannot be beaten by a maximiser
 * that is looking on the wrong side of the face, in the wrong element, or at
 * the wrong field, which are the three ways this could be wrong and still
 * converge.
 *
 * THE COMPARISON IS ONE-SIDED AND THAT IS THE SHARP DIRECTION. A sampler on a
 * finite grid can only UNDER-state a maximum, so the located value must be
 * greater than or equal to every sample, to round-off; it may legitimately
 * exceed the best sample, by however much a degree-k polynomial rises between
 * two sample points. Asserting equality would be asserting that the sampler is
 * as good as the search.
 */
BOOST_AUTO_TEST_CASE( theContactIsFoundWhereTheFieldIsLargestOnThePolygon )
{
	int const order = 2;
	meq::tests::Rectangle const box = meq::tests::standardBox();

	std::printf( "\n  THE LOCATED CONTACT AGAINST A DENSE SAMPLE OF THE POLYGON\n" );
	std::printf( "    %5s %7s %14s %14s %14s %12s %12s\n",
	             "n", "newton", "psi_ax", "psi_bnd", "sampled max",
	             "contact r", "contact z" );

	for ( int n : { 12, 24, 48 } )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );
		paintLimiter( mesh );

		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		auto ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
		meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );

		mfem::ConstantCoefficient zero( 0.0 );
		mfem::FunctionCoefficient guess = bump( 0.30 );

		GradShafranovSolver solver( mesh, order );
		solver.setLimiterSurface( limiterAttribute );
		solver.setSource( source, 0.30 );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 40 );
		solver.solve();

		BOOST_TEST_REQUIRE( solver.limiterContactWasLocated(),
			"no contact was located at n = " << n );

		double sampledR = 0.0, sampledZ = 0.0;
		double const sampled = sampledMaximum( mesh, solver.potential(), 401,
		                                       sampledR, sampledZ );

		std::printf( "    %5d %7d %14.6e %14.6e %14.6e %12.6f %12.6f\n",
		             n, solver.newtonIterations(), solver.psiAxis(),
		             solver.psiBoundary(), sampled, solver.limiterContactR(),
		             solver.limiterContactZ() );

		BOOST_TEST( solver.psiBoundary() >= sampled - 1.0e-11,
			"psi_bnd at n = " << n << " is " << solver.psiBoundary()
			<< " but a dense sample of the polygon found " << sampled
			<< " at ( " << sampledR << ", " << sampledZ << " ). The constraint "
			"is a MAXIMUM over the polygon, so nothing on it may beat the "
			"located contact -- a value that is beaten is the search reading "
			"the wrong side of a face, the wrong element, or the wrong field." );

		// AND IT MUST BE ON THE POLYGON, which is the other half of "max over
		// the polygon" and is not implied by the value being large.
		double const dR = std::min( std::fabs( solver.limiterContactR() - limiterRMin ),
		                            std::fabs( solver.limiterContactR() - limiterRMax ) );
		double const dZ = std::min( std::fabs( solver.limiterContactZ() - limiterZMin ),
		                            std::fabs( solver.limiterContactZ() - limiterZMax ) );
		BOOST_TEST( std::min( dR, dZ ) < 1.0e-12,
			"the contact at n = " << n << " is ( " << solver.limiterContactR()
			<< ", " << solver.limiterContactZ() << " ), which is not on the "
			"rectangle's boundary" );
	}
}

namespace
{
	/// One solve of the fixture, with whatever limiter constraint the caller
	/// installs. Returns psi_bnd; reports psi_ax, the contact and the Newton
	/// history through the solver it was handed.
	struct Solved
	{
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		double contactR = 0.0;
		double contactZ = 0.0;
		int newton = 0;
		std::vector<double> residuals;
	};

	/// `located` chooses the constraint: the polygon when true, the prescribed
	/// point ( r, z ) when false. EVERYTHING ELSE IS HELD, which is what makes
	/// the difference between two calls a statement about the constraint.
	Solved solveWith( int n, int order, bool located, double r, double z )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );
		paintLimiter( mesh );

		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		auto ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
		meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );

		mfem::ConstantCoefficient zero( 0.0 );
		mfem::FunctionCoefficient guess = bump( 0.30 );

		GradShafranovSolver solver( mesh, order );
		if ( located )
			solver.setLimiterSurface( limiterAttribute );
		else
			solver.setBoundaryFluxPoint( r, z );
		solver.setSource( source, 0.30 );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 40 );
		solver.solve();

		Solved out;
		out.psiAxis = solver.psiAxis();
		out.psiBoundary = solver.psiBoundary();
		out.contactR = solver.limiterContactR();
		out.contactZ = solver.limiterContactZ();
		out.newton = solver.newtonIterations();
		out.residuals = solver.newtonResiduals();
		return out;
	}
}

/*
 * THE FOUND CONTACT AGAINST THE EXACT ONE, WITH A WRONG PRESCRIBED POINT AS THE
 * CONTROL.
 *
 * THIS FIXTURE KNOWS ITS OWN ANSWER, WHICH IS WHY IT WAS CHOSEN. The domain,
 * the boundary condition and the source are all symmetric in `z`, so the exact
 * `psi` is, and the maximum over the rectangle's outboard edge is at `z = 0`
 * EXACTLY. So the contact is a known point and the located one can be measured
 * against it rather than against another run.
 *
 * THE DISCRETE CONTACT IS NOT AT z = 0 AND THAT IS THE MESH, NOT AN ERROR.
 * MakeCartesian2D splits every cell along one diagonal, which is not symmetric
 * in `z`, so `psi_h` is not either. The located contact converges to the
 * midplane rather than sitting on it, and the rate is what says the search is
 * finding a property of the FIELD rather than of the mesh.
 *
 * THE CONTROL IS A PRESCRIBED POINT THAT IS WRONG, and it is wrong the way a
 * real one is: FREE-BOUNDARY-PLAN.md section 7.20 records `freegs4e`'s own
 * contact moving 0.6 m round its limiter circle between two grids, because a
 * contact read off a grid is a maximum over CELLS. A contact prescribed
 * somewhere else on the SAME polygon is that failure in miniature, and it is
 * the column that says the located constraint is doing work: without it every
 * assertion below is satisfied by a solver that quietly ignores the polygon.
 */
BOOST_AUTO_TEST_CASE( theFoundContactAgreesWithTheExactOneAndAWrongOneDoesNot )
{
	int const order = 2;

	// Known by symmetry: the outboard edge, on the midplane.
	double const exactR = limiterRMax;
	double const exactZ = 0.0;
	// Half way up the same edge -- still ON the polygon, so this is not a
	// question about whether the point is legal.
	double const wrongZ = 0.15;

	std::printf( "\n  THE FOUND CONTACT, THE EXACT ONE, AND A WRONG ONE\n" );
	std::printf( "    %5s %14s %14s %14s %12s %12s\n",
	             "n", "found", "exact point", "wrong point", "|found-ex|",
	             "|wrong-ex|" );

	std::vector<double> gap, control, contactOffset;
	for ( int n : { 12, 24, 48 } )
	{
		Solved const found = solveWith( n, order, true, 0.0, 0.0 );
		Solved const exact = solveWith( n, order, false, exactR, exactZ );
		Solved const wrong = solveWith( n, order, false, exactR, wrongZ );

		gap.push_back( std::fabs( found.psiBoundary - exact.psiBoundary ) );
		control.push_back( std::fabs( wrong.psiBoundary - exact.psiBoundary ) );
		contactOffset.push_back( std::fabs( found.contactZ - exactZ ) );

		std::printf( "    %5d %14.6e %14.6e %14.6e %12.3e %12.3e\n",
		             n, found.psiBoundary, exact.psiBoundary,
		             wrong.psiBoundary, gap.back(), control.back() );
	}

	std::printf( "\n    the located contact's own distance from the midplane,"
	             " which is exact by symmetry:\n" );
	for ( std::size_t i = 0; i < contactOffset.size(); ++i )
	{
		std::printf( "      %12.3e", contactOffset[ i ] );
		if ( i )
			std::printf( "   rate %6.3f",
			             std::log2( contactOffset[ i - 1 ]/contactOffset[ i ] ) );
		std::printf( "\n" );
	}

	// THE CONTACT IS AN OUTPUT AND IT CONVERGES. A search that returned a mesh
	// artefact -- an element centre, a dof, a face midpoint -- would sit at a
	// fixed fraction of h from the midplane and read a rate of 1, or would not
	// move at all.
	double const contactRate =
		std::log2( contactOffset.front()/contactOffset.back() )/2.0;
	std::printf( "\n    contact position converges at %.3f\n", contactRate );
	BOOST_TEST( contactRate > 1.8,
		"the located contact approaches the midplane at only " << contactRate
		<< ", where the field it is a stationary point of converges at "
		<< order + 1 << ". A rate near 1 is a search returning a mesh entity "
		"rather than a property of psi_h." );

	// AND FINDING IT BEATS BEING TOLD IT WRONG, BY ORDERS. The gap against the
	// exactly-prescribed point is the envelope theorem paying: the contact is a
	// STATIONARY point, so an O( h^2 ) error in WHERE it is costs O( h^4 ) in
	// the VALUE there.
	for ( std::size_t i = 0; i < gap.size(); ++i )
		BOOST_TEST( gap[ i ] < 0.05*control[ i ],
			"at mesh " << i << " the located contact is " << gap[ i ]
			<< " from the exactly-prescribed one while a contact prescribed "
			"0.15 m away on the same polygon is " << control[ i ]
			<< ". The located constraint has to be much closer to the right "
			"answer than a plausible wrong point, or it is not finding "
			"anything." );

	// THE CONTROL MUST NOT CONVERGE, which is what makes it a control. A wrong
	// contact is wrong by dist x |grad psi| however fine the mesh.
	BOOST_TEST( control.back() > 0.5*control.front(),
		"the wrongly-prescribed contact's error fell from " << control.front()
		<< " to " << control.back() << " under refinement. It must not: it is a "
		"different point of the same curve, so its error is a property of the "
		"SOLUTION and not of the mesh. If it converges, this column is not a "
		"control and the comparison above is empty." );
}

/*
 * NEWTON'S OWN ORDER, WHICH IS THE ONLY THING THAT CAN SEE THE BORDER ROW.
 *
 * The row this constraint contributes is the contact element's potential shape
 * functions AT the contact, undifferenced, with a corner of exactly 1 -- and
 * that is only correct because of the envelope theorem. `x*` moves with the
 * solution, so the chain rule carries a term `grad psi . d( x* )/dlambda`, and
 * the claim is that it vanishes: on the interior of an edge because the
 * TANGENTIAL derivative is zero at a maximum of the restriction while
 * `d( x* )/dlambda` is tangential, and at a vertex because `x*` does not move
 * at all.
 *
 * NOTHING IN AN ERROR NORM CAN CHECK THAT. This file's own table above would
 * be unchanged by a row missing the position term: Newton converges to the same
 * discrete solution whatever Jacobian carried it there, which is the finding
 * CLAUDE.md records under *A wrong Jacobian is invisible to a convergence
 * table* -- perturbing dF/dpsi by 5% leaves every error and every rate
 * unchanged to six figures and drops the observed order to exactly 1.000. So
 * the order is the assertion, and a missing envelope term would show up here
 * and nowhere else.
 *
 * BOUNDED BOTH SIDES, per this tree's standing caution: 1 is a broken Jacobian
 * and anything near 4 is the iterate walking into its basin rather than a rate.
 * The order is taken on the TAIL -- the last triple above the round-off floor
 * -- and not on the best triple anywhere in the history.
 */
BOOST_AUTO_TEST_CASE( theBorderRowKeepsNewtonsQuadraticOrder )
{
	int const order = 2;

	std::printf( "\n  NEWTON ON THE LOCATED-CONTACT BORDER\n" );

	for ( int n : { 24, 48 } )
	{
		Solved const found = solveWith( n, order, true, 0.0, 0.0 );

		std::printf( "    n = %d, %d iterations\n", n, found.newton );
		for ( std::size_t i = 0; i < found.residuals.size(); ++i )
		{
			std::printf( "      %2zu  %14.6e", i, found.residuals[ i ] );
			if ( i >= 2 )
				std::printf( "   order %6.3f",
				             meq::tests::newtonOrder( found.residuals[ i - 2 ],
				                                      found.residuals[ i - 1 ],
				                                      found.residuals[ i ] ) );
			std::printf( "\n" );
		}

		// The tail: the last triple whose middle and last are above the floor.
		double const floorLevel = 1.0e-13*found.residuals.front();
		double tail = 0.0;
		for ( std::size_t i = 2; i < found.residuals.size(); ++i )
			if ( found.residuals[ i ] > floorLevel )
				tail = meq::tests::newtonOrder( found.residuals[ i - 2 ],
				                                found.residuals[ i - 1 ],
				                                found.residuals[ i ] );

		std::printf( "      tail order %.3f\n", tail );

		BOOST_TEST( tail > 1.5,
			"the bordered Newton's tail order at n = " << n << " is " << tail
			<< ". A rate near 1 on a Newton method is a Jacobian statement, and "
			"the Jacobian entry this constraint owns is the border row -- so "
			"the first thing to suspect is the envelope argument that lets the "
			"row be the shape functions with no sensitivity of the contact in "
			"it. See LimiterConstraint::LocatedContact." );
		BOOST_TEST( tail < 3.5,
			"the tail order at n = " << n << " reads " << tail
			<< ", which is not a rate -- it is a short or non-monotone history "
			"being read as one. Print the whole history before believing it." );
	}
}

/*
 * WHAT IT REFUSES, AND THE CONTROL THAT STOPS THE REFUSALS BEING VACUOUS.
 *
 * An EMPTY polygon is the failure worth catching, because it is the one that
 * does not announce itself: `max` over an empty set is minus infinity and its
 * border row is all zeroes, so the bordered solve does not diverge -- it
 * converges, to a `psi_bnd` pinned by nothing at all. Asking for an attribute
 * the mesh does not carry is how a caller reaches that, and the commonest way
 * is a mesh built without `--limiter`.
 *
 * AND A PRESCRIBED CONTACT BESIDE A FOUND ONE IS REFUSED RATHER THAN RANKED.
 * Both constraints converge, to equilibria differing by the O( h ) the point
 * version costs -- section 7.20 measures that -- so a precedence rule would
 * decide which equilibrium is reported on the strength of call order.
 *
 * THE CONTROL IS THE LAST BLOCK: a setter that refused everything would satisfy
 * every assertion above it.
 */
BOOST_AUTO_TEST_CASE( anEmptyPolygonAndADoublyPinnedContactAreBothRefused )
{
	meq::tests::Rectangle const box = meq::tests::standardBox();

	{
		// A mesh with no limiter region at all -- every element attribute 1,
		// which is what a mesh built without --limiter looks like.
		mfem::Mesh mesh = meq::tests::makeMesh( box, 12 );
		GradShafranovSolver solver( mesh, 2 );
		BOOST_CHECK_THROW( solver.setLimiterSurface( limiterAttribute ),
		                   std::runtime_error );
	}

	{
		mfem::Mesh mesh = meq::tests::makeMesh( box, 12 );
		paintLimiter( mesh );
		GradShafranovSolver solver( mesh, 2 );
		solver.setBoundaryFluxPoint( 1.2, 0.0 );
		BOOST_CHECK_THROW( solver.setLimiterSurface( limiterAttribute ),
		                   std::logic_error );
	}

	{
		mfem::Mesh mesh = meq::tests::makeMesh( box, 12 );
		paintLimiter( mesh );
		GradShafranovSolver solver( mesh, 2 );
		BOOST_CHECK_THROW( solver.setLimiterSurface( 0 ),
		                   std::invalid_argument );
	}

	{
		// THE CONTROL. The same call on the same mesh with the region present
		// must be accepted, and must report the constraint it installed.
		mfem::Mesh mesh = meq::tests::makeMesh( box, 12 );
		paintLimiter( mesh );
		GradShafranovSolver solver( mesh, 2 );
		BOOST_CHECK_NO_THROW( solver.setLimiterSurface( limiterAttribute ) );
		BOOST_TEST( ( solver.limiterConstraint()
		              == GradShafranovSolver::LimiterConstraint::LocatedContact ),
			"setLimiterSurface() did not install LocatedContact, so the three "
			"refusals above are compatible with a setter that refuses "
			"everything" );
	}
}
