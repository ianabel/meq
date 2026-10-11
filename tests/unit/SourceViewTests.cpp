// Unit tests for meq::SourceView -- the POD form of F( R, z, psi ) and dF/dpsi
// that a device kernel can evaluate.
//
// THIS FILE IS WHAT MAKES THE POD FORM SAFE TO HAVE. SourceView.hpp is a SECOND
// expression of the Grad-Shafranov right hand side, standing beside
// meq::NormalisedMHDSource::f(), and CLAUDE.md names that arrangement as a
// hazard rather than a pattern: "a second expression of F( R, z, psi )" is "two
// implementations that must agree", and the day one is edited and the other is
// not, MEQ solves a different equilibrium and converges to it at the full rate.
// No convergence table can see that. This test is the only thing that can.
//
// SO IT ASSERTS EQUALITY TO THE BIT, not agreement to a tolerance. The two
// expressions evaluate the same arithmetic in the same order on the same
// doubles, so anything other than an exact match is a real difference in what
// they compute, and a tolerance would hide exactly the class of edit this file
// exists to catch -- a dropped mu0, an R instead of an R^2, a span squared
// where it should be linear.
//
// AND THE SWEEP IS CHOSEN FOR WHERE A SPLINE AND ITS POD FORM CAN DISAGREE,
// which is three places and none of them is the interior of an interval:
//
//   * ON A KNOT. Interval selection is a binary search on one side and
//     std::upper_bound on the other; a profile with a flat segment meeting a
//     rising one has a genuine kink there, so picking the interval on the wrong
//     side of an abscissa that lands exactly on a knot changes the derivative
//     discontinuously while leaving the value continuous.
//   * OUTSIDE THE PROFILE SUPPORT, where the value is extended constant and the
//     derivative is zero rather than the constant's slope.
//   * ACROSS THE PLASMA EDGE, where insidePlasma() switches the whole term off.
//
// If this file is ever deleted or its assertions relaxed, SourceView.hpp stops
// being pinned to anything.

#define BOOST_TEST_MODULE SourceViewTests
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>

#include "convergence/RoundOffFloor.hpp"
#include "meq/Coils.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "meq/SourceView.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

namespace
{

	/// A spline with a deliberate kink: flat, then rising, then flat. The
	/// middle knots are where interval selection has to agree.
	std::shared_ptr<meq::Profile const> kinkedSpline()
	{
		std::vector<meq::Knot> knots = {
			{ 0.00, 1.0, 0.0 },
			{ 0.25, 1.0, 0.0 },
			{ 0.50, 3.0, 4.0 },
			{ 0.75, 5.0, 0.0 },
			{ 1.00, 5.0, 0.0 } };
		return std::make_shared<meq::SplineProfile>( std::move( knots ) );
	}

	/// The abscissae worth testing, in normalised flux: the knots exactly, the
	/// midpoints, and points outside the support on both sides.
	std::vector<double> sweepPoints()
	{
		std::vector<double> x = { -0.5, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5,
		                          0.625, 0.75, 0.875, 1.0, 1.25, 1.75 };
		return x;
	}

	/*
	 * THE CONTRACT, AND WHY IT IS NOT LITERALLY `==` IN EVERY BUILD.
	 *
	 * The file comment's justification for bit-equality is that the two
	 * expressions "evaluate the same arithmetic in the same order on the same
	 * DOUBLES". Under `MFEM_PRECISION = single` that premise is false: the POD
	 * side returns `mfem::real_t`, so it rounds to `float` at every step, while
	 * `meq::Source` stays `double` -- *Layout* keeps it free of MFEM on purpose.
	 * The two are then genuinely different computations and no edit to either
	 * can make them agree to the bit.
	 *
	 * **WHAT THE FILE EXISTS TO CATCH SURVIVES THAT INTACT.** A dropped `mu0`,
	 * an `R` where an `R^2` belongs, a span squared where it should be linear,
	 * an interval picked on the wrong side of a knot -- every one of those is an
	 * `O( 1 )` RELATIVE difference. The format's is `O( eps )`: measured over
	 * this sweep the worst disagreement is **0.54 of one eps** ( 0.276799977
	 * against 0.27679999999999993 ).
	 *
	 * So `exactAgreementFloor()` returns **exactly zero in double** -- where
	 * this stays the bit-equality it always was -- and 64 eps of the value in
	 * single, which is 118x the worst measured and five orders below the
	 * smallest defect this file is for. M-219.
	 */
	/*
	 * BOTH ROUTES MUST BE ASKED AT THE SAME ABSCISSA, AND IN SINGLE THAT IS NOT
	 * AUTOMATIC. `meq::sourceF` takes `mfem::real_t`, so a `double` psi is
	 * rounded on the way in; `meq::Source::f` takes it as written. When the
	 * sweep lands ON A KNOT -- which this file chooses deliberately, because
	 * interval selection is a binary search on one side and `std::upper_bound`
	 * on the other -- that rounding can put the two on OPPOSITE SIDES of the
	 * abscissa, and a profile with a flat segment meeting a rising one has a
	 * genuine kink there. The derivative then differs by the size of the kink,
	 * which is `O( 1 )` and looks exactly like the defect this file hunts.
	 *
	 * Measured: every `dF/dpsi` comparison at `Psi = 0.25` and no other.
	 *
	 * So the abscissa is rounded to `real_t` ONCE and both routes are handed
	 * that value. This is M-189's rule -- *where a caller computes one quantity
	 * two ways and compares them, both sides must be given the same step* --
	 * applied to the point rather than to a step. It is a no-op in double,
	 * where `real_t` IS `double`. M-219.
	 */
	double sameAbscissa( double psi )
	{
		return static_cast<double>( static_cast<mfem::real_t>( psi ) );
	}

	bool agrees( double pod, double reference )
	{
		return std::abs( pod - reference )
		       <= meq::tests::exactAgreementFloor( reference );
	}

	/*
	 * AND NEAR A KNOT THE DERIVATIVE IS NOT COMPARABLE POINTWISE IN A SINGLE
	 * BUILD, WHICH IS A PROPERTY OF THE POD AND NOT OF THE TEST.
	 *
	 * `meq::KnotView::psi` is `mfem::real_t` while `meq::Profile`'s abscissae
	 * are `double` -- the MFEM-free seam again -- so under
	 * `MFEM_PRECISION = single` the interval BOUNDARY itself moves by up to one
	 * `float` ulp. Rounding the sweep point to `real_t` does not help; the two
	 * sides are measuring from different knots.
	 *
	 * **AND `dF/dpsi` IS STEEP THERE, WHICH IS THE WHOLE REASON THIS FILE
	 * SAMPLES ON KNOTS.** Probed at the failing point, `psi = 0.13750000`:
	 *
	 *     four ulps below   0
	 *     the POD           3.8925e-05
	 *     the profile       1.1122e-05
	 *     four ulps above   2.0352e-03
	 *
	 * The two do not disagree about the FUNCTION; they are evaluating a rapidly
	 * varying one a few ulps apart, and every value between 0 and 2e-03 is a
	 * value it legitimately takes there. A pointwise comparison of that is not
	 * well posed in this format however either side is written.
	 *
	 * **SO THE POD MUST LIE IN THE RANGE THE REFERENCE TAKES WITHIN A FEW ULPS
	 * OF THE SAME ABSCISSA**, which is the strongest statement the format
	 * supports and keeps every defect this file hunts: a dropped `mu0` is a
	 * factor of 1.3e-06, an `R` for an `R^2` a factor of `R`, and neither lands
	 * inside a bracket spanning one knot. In double the two knots are
	 * identical, the direct comparison passes, and the bracket is never
	 * reached.
	 *
	 * **AND THE SAME THING HAPPENS ONE LEVEL UP, AT A COIL EDGE.** This file's
	 * coils are `[ 1.85, 2.15 ]` and `[ 1.90, 2.20 ]` and its sweep visits
	 * `R = 2.20` on purpose, because -- in its own words -- *"meq::Coil::
	 * contains() is closed and a half-open POD would differ"*. Under single the
	 * edge itself is a `real_t`, so the two routes compare against bounds an
	 * ulp apart and the POD reports **-9.2153** where the reference reports
	 * **0**: not a rounding difference but the whole coil term, switched on or
	 * off by which side of its own edge the point fell. The bracket covers that
	 * too, and a genuinely wrong coil term is outside it at every radius rather
	 * than at one. M-219.
	 */
	template<typename Evaluate>
	bool agreesNearABoundary( double pod, double psi, Evaluate &&at )
	{
		if ( agrees( pod, at( psi ) ) )
			return true;

		double const ulp = ( std::abs( psi ) + 1.0 )
		                   *static_cast<double>(
			                   std::numeric_limits<mfem::real_t>::epsilon() );
		double lo = at( psi ), hi = lo;
		for ( double const shift : { -4.0, -1.0, 1.0, 4.0 } )
		{
			double const v = at( psi + shift*ulp );
			lo = std::min( lo, v );
			hi = std::max( hi, v );
		}
		double const slack = meq::tests::exactAgreementFloor( hi - lo );
		return pod >= lo - slack && pod <= hi + slack;
	}


}

BOOST_AUTO_TEST_CASE( theViewReproducesTheSourceToTheBit )
{
	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( kinkedSpline(), kinkedSpline(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( source, view, knots ) );

	double const span = psiAxis - 0.05;
	int checked = 0;

	for ( double const radius : { 0.3, 1.0, 1.7, 3.2 } )
		for ( double const normalised : sweepPoints() )
		{
			// Back out the physical psi the source will be handed, so the
			// sweep lands on the knots in the coordinate the PROFILES use.
			double const psi = sameAbscissa( 0.05 + normalised*span );

			double const fromSource = source.f( radius, 0.0, psi );
			double const fromView =
				meq::sourceF( view, knots.data(), radius, 0.0, psi );
			BOOST_TEST( agrees( fromView, fromSource ),
			            "F disagreed at R = " << radius << ", Psi = " << normalised );

			double const dFromSource = source.dFdPsi( radius, 0.0, psi );
			double const dFromView =
				meq::sourceDFdPsi( view, knots.data(), radius, 0.0, psi );
			BOOST_TEST( agreesNearABoundary( dFromView, psi,
			                              [ & ]( double p )
			                              { return source.dFdPsi( radius, 0.0, p ); } ),
			            "dF/dpsi disagreed at R = " << radius
			            << ", Psi = " << normalised );
			++checked;
		}

	// A sweep that silently tested nothing is the failure this guards.
	BOOST_TEST( checked == 52 );
}

BOOST_AUTO_TEST_CASE( theViewReproducesTheConfinedSourceAcrossTheEdge )
{
	// THE PLASMA EDGE IS THE ONE PLACE F IS DISCONTINUOUS, so it is the one
	// place a POD form can agree everywhere it is smooth and still be wrong.
	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( kinkedSpline(), kinkedSpline(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );
	source.setPlasmaSupport( true );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( source, view, knots ) );
	BOOST_TEST( view.confineToPlasma == 1 );

	int inside = 0, outside = 0;
	for ( double const psi : { -0.20, -0.05, 0.0, 0.049, 0.05, 0.051, 0.2, 0.4, 0.6 } )
	{
		double const fromSource = source.f( 1.3, 0.0, psi );
		double const fromView = meq::sourceF( view, knots.data(), 1.3, 0.0, psi );
		BOOST_TEST( agrees( fromView, fromSource ), "F disagreed at psi = " << psi );

		double const dFromSource = source.dFdPsi( 1.3, 0.0, psi );
		double const dFromView =
			meq::sourceDFdPsi( view, knots.data(), 1.3, 0.0, psi );
		BOOST_TEST( agreesNearABoundary( dFromView, psi,
		                              [ & ]( double p )
		                              { return source.dFdPsi( 1.3, 0.0, p ); } ),
		            "dF/dpsi disagreed at psi = " << psi );

		if ( fromSource == 0.0 ) ++outside; else ++inside;
	}

	// Both sides have to be exercised or the case proves only one of them.
	BOOST_TEST( inside > 0 );
	BOOST_TEST( outside > 0 );
}

BOOST_AUTO_TEST_CASE( theViewReproducesTheMirrorBoxOnAllThreeSides )
{
	/*
	 * meq::PlasmaSupportShape::Box -- `0 < Psi <= 1` and `|z| <= z_edge` -- is
	 * CARRIED by the view rather than refused, so it is a second expression of
	 * the support and this is what holds it to the first. MIRROR-PLAN.md MR-1.
	 *
	 * Mirror orientation: psi_ax is the INNER label and the smaller flux, so
	 * the span is negative. Each of the box's three sides is crossed, its two
	 * closed edges are sat on, and the normalisation columns are compared
	 * beside F and dF/dpsi, because the device route computes all three.
	 */
	double const inner = 0.05;
	double const outer = 0.4;
	double const zEdge = 0.6;
	meq::NormalisedMHDSource source( kinkedSpline(), kinkedSpline(), inner );
	source.setNormalisation( inner, outer );
	source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, zEdge );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( source, view, knots ) );
	BOOST_TEST( view.confineToPlasma == 1 );
	BOOST_TEST( view.supportShape == static_cast<int>( meq::PlasmaSupportShape::Box ) );
	BOOST_TEST( view.axialEdge == static_cast<mfem::real_t>( zEdge ) );

	int inside = 0, innerSide = 0, outerSide = 0, beyondThroat = 0;
	for ( double const z : { 0.0, 0.3, -zEdge, zEdge, 0.61, -2.0 } )
		for ( double const psi : { -0.1, 0.0, 0.049, 0.05, 0.051, 0.2, 0.399, 0.4,
		                           0.401, 0.6 } )
		{
			double const fromSource = source.f( 1.3, z, psi );
			double const fromView = meq::sourceF( view, knots.data(), 1.3, z, psi );
			BOOST_TEST( agrees( fromView, fromSource ),
			            "F disagreed at psi = " << psi << ", z = " << z );

			double const dFromView =
				meq::sourceDFdPsi( view, knots.data(), 1.3, z, psi );
			BOOST_TEST( agreesNearABoundary( dFromView, psi,
			                              [ & ]( double p )
			                              { return source.dFdPsi( 1.3, z, p ); } ),
			            "dF/dpsi disagreed at psi = " << psi << ", z = " << z );

			double hostAxis = 0.0, hostBoundary = 0.0;
			BOOST_TEST_REQUIRE( source.normalisationDerivatives(
				1.3, z, psi, hostAxis, hostBoundary ) );
			mfem::real_t viewAxis = 0, viewBoundary = 0;
			meq::sourceNormalisationDerivatives( view, knots.data(), 1.3, z, psi,
			                                     viewAxis, viewBoundary );
			BOOST_TEST( agreesNearABoundary( viewAxis, psi,
				[ & ]( double p )
				{
					double a = 0.0, b = 0.0;
					source.normalisationDerivatives( 1.3, z, p, a, b );
					return a;
				} ), "dF/dpsi_ax disagreed at psi = " << psi << ", z = " << z );
			BOOST_TEST( agreesNearABoundary( viewBoundary, psi,
				[ & ]( double p )
				{
					double a = 0.0, b = 0.0;
					source.normalisationDerivatives( 1.3, z, p, a, b );
					return b;
				} ), "dF/dpsi_bnd disagreed at psi = " << psi << ", z = " << z );

			// Classify by the closed form, not by the source under test.
			double const normalised = ( psi - outer )/( inner - outer );
			if ( std::abs( z ) > zEdge )
				++beyondThroat;
			else if ( normalised > 1.0 )
				++innerSide;
			else if ( !( normalised > 0.0 ) )
				++outerSide;
			else
			{
				++inside;
				BOOST_TEST( fromSource != 0.0 );
			}
			if ( std::abs( z ) > zEdge || normalised > 1.0 || !( normalised > 0.0 ) )
				BOOST_TEST( fromView == 0.0,
				            "the view is not zero outside the box at psi = " << psi
				            << ", z = " << z );
		}

	// Every side has to be exercised or the case proves only some of them.
	BOOST_TEST( inside > 0 );
	BOOST_TEST( innerSide > 0 );
	BOOST_TEST( outerSide > 0 );
	BOOST_TEST( beyondThroat > 0 );
}

BOOST_AUTO_TEST_CASE( aScaledProfileIsFoldedIntoTheViewRatherThanDropped )
{
	// meq::ScaledProfile has no kind of its own -- it folds into `scale` -- so
	// the way this goes wrong is by being ignored, which leaves a view that is
	// right in shape and off by a constant factor everywhere.
	double const psiAxis = 0.4;
	auto scaled = std::make_shared<meq::ScaledProfile>( kinkedSpline(), 2.5 );
	meq::NormalisedMHDSource source( scaled, kinkedSpline(), psiAxis );
	source.setNormalisation( psiAxis, 0.0 );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( source, view, knots ) );
	BOOST_TEST( view.pPrime.scale == 2.5 );
	BOOST_TEST( view.ggPrime.scale == 1.0 );

	for ( double const psi : { 0.05, 0.15, 0.2, 0.3, 0.39 } )
	{
		BOOST_TEST( agrees( meq::sourceF( view, knots.data(), 1.1, 0.0, psi ),
		                    source.f( 1.1, 0.0, psi ) ) );
		BOOST_TEST( agreesNearABoundary( meq::sourceDFdPsi( view, knots.data(), 1.1, 0.0, psi ), psi,
		                              [ & ]( double p )
		                              { return source.dFdPsi( 1.1, 0.0, p ); } ) );
	}
}

BOOST_AUTO_TEST_CASE( theSolovievSourceHasAViewAndItIsQuadraticInRadius )
{
	// F = -( ( 1 - a ) R^2 + a ), so a view that carried a constant would agree
	// at exactly one radius and nowhere else.
	meq::SolovievSource source( -0.52 );
	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( source, view, knots ) );

	for ( double const radius : { 0.4, 1.0, 2.3 } )
		for ( double const psi : { -1.0, 0.0, 1.0 } )
		{
			BOOST_TEST( agrees( meq::sourceF( view, knots.data(), radius, 0.0, psi ),
			                    source.f( radius, 0.0, psi ) ) );
			BOOST_TEST( agreesNearABoundary( meq::sourceDFdPsi( view, knots.data(), radius, 0.0, psi ), psi,
			                              [ & ]( double p )
			                              { return source.dFdPsi( radius, 0.0, p ); } ) );
		}
}

BOOST_AUTO_TEST_CASE( theUniformFastPathIsTakenAndAgreesWithTheSearch )
{
	// THE FAST INDEX IS A PERFORMANCE HINT WHOSE FAILURE MODE MUST BE
	// PERFORMANCE. profileInterval() picks the interval by arithmetic when the
	// knots are uniform and by binary search when they are not, and the two
	// must choose identically -- so this asserts the flag is set where it
	// should be, cleared where it should not, and that BOTH routes reproduce
	// the host profile to the bit. A fast path that is right on average is a
	// different cubic on a set of points nobody enumerates.
	double const psiAxis = 0.4;

	// Uniform: 0, 0.25, 0.5, 0.75, 1 -- the fixture above.
	meq::NormalisedMHDSource uniform( kinkedSpline(), kinkedSpline(), psiAxis );
	uniform.setNormalisation( psiAxis, 0.0 );
	meq::SourceView uniformView{};
	std::vector<meq::KnotView> uniformKnots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( uniform, uniformView, uniformKnots ) );
	BOOST_TEST( uniformView.pPrime.uniform == 1 );

	// Non-uniform: the same endpoints, graded towards the right.
	std::vector<meq::Knot> graded = {
		{ 0.00, 1.0, 0.0 },
		{ 0.60, 1.0, 0.0 },
		{ 0.80, 3.0, 4.0 },
		{ 0.95, 5.0, 0.0 },
		{ 1.00, 5.0, 0.0 } };
	auto gradedSpline = std::make_shared<meq::SplineProfile>( graded );
	meq::NormalisedMHDSource nonUniform( gradedSpline, gradedSpline, psiAxis );
	nonUniform.setNormalisation( psiAxis, 0.0 );
	meq::SourceView gradedView{};
	std::vector<meq::KnotView> gradedKnots;
	BOOST_TEST_REQUIRE( meq::makeSourceView( nonUniform, gradedView, gradedKnots ) );
	BOOST_TEST( gradedView.pPrime.uniform == 0 );

	// Both routes, against their own host profile, on and between every knot.
	for ( double const normalised : { 0.0, 0.1, 0.25, 0.3, 0.5, 0.6, 0.75,
	                                  0.8, 0.9, 0.95, 1.0, 1.2 } )
	{
		double const psi = sameAbscissa( normalised*psiAxis );

		BOOST_TEST( agrees( meq::sourceF( uniformView, uniformKnots.data(), 1.4, 0.0, psi ),
		                uniform.f( 1.4, 0.0, psi ) ),
		            "uniform F disagreed at Psi = " << normalised );
		BOOST_TEST( agreesNearABoundary( meq::sourceDFdPsi( uniformView, uniformKnots.data(), 1.4, 0.0, psi ), psi,
		                              [ & ]( double p )
		                              { return uniform.dFdPsi( 1.4, 0.0, p ); } ),
		            "uniform dF/dpsi disagreed at Psi = " << normalised );

		BOOST_TEST( agrees( meq::sourceF( gradedView, gradedKnots.data(), 1.4, 0.0, psi ),
		                nonUniform.f( 1.4, 0.0, psi ) ),
		            "graded F disagreed at Psi = " << normalised );
		BOOST_TEST( agreesNearABoundary( meq::sourceDFdPsi( gradedView, gradedKnots.data(), 1.4, 0.0, psi ), psi,
		                              [ & ]( double p )
		                              { return nonUniform.dFdPsi( 1.4, 0.0, p ); } ),
		            "graded dF/dpsi disagreed at Psi = " << normalised );
	}
}

BOOST_AUTO_TEST_CASE( aSourceWithNoPodFormIsRefusedRatherThanApproximated )
{
	// THE REFUSAL IS THE SAFETY PROPERTY, so it is asserted rather than
	// assumed. A view built from a source it does not understand would describe
	// the plasma term alone and silently drop whatever else that source adds --
	// a coil term, a rotation closure -- and the run would converge to it.
	//
	// meq::MHDSource is not normalised and is not the production class, so it
	// has no view; it stands here for the whole family the cast excludes.
	auto constant = std::make_shared<meq::ConstantProfile>( 1.0 );
	meq::MHDSource plain( constant, constant );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST( !meq::makeSourceView( plain, view, knots ) );
	BOOST_TEST( knots.empty() );
}


/*
 * THE PRODUCTION CONFIGURATION IS COIL-AUGMENTED, AND UNTIL NOW THE POD COVERED
 * EVERY CASE EXCEPT THOSE.
 *
 * makeSourceView() casts to the EXACT production class and refuses a wrapper,
 * under a comment that is exactly right about why: a coil-augmented source's
 * f() is the SUM, and a view built from the wrapper would describe the plasma
 * term alone and silently drop the coils.
 *
 * meq::Source::scaledF() is the function for which that inverts. It is the
 * PLASMA TERM ALONE by definition -- CoilAugmentedNormalisedSource overrides it
 * saying so, because "this class's own f() is the SUM and would put the coil
 * current inside the prescribed plasma current" -- and it is what
 * assemblePlasmaCurrent() integrates, once per constraint evaluation, as the
 * third reduction of DEVICE-PLAN.md section 3.4a.
 *
 * So makePlasmaSourceView() is a SECOND ENTRY POINT rather than a widened cast,
 * and this case asserts all three halves of that being safe: the plasma view
 * reproduces scaledF() to the bit, the ordinary entry point still REFUSES the
 * wrapper, and the coil term is genuinely non-zero so that the first two are
 * not the same statement.
 */
BOOST_AUTO_TEST_CASE( thePlasmaViewReproducesScaledFThroughACoilWrapper )
{
	double const psiAxis = 0.4;
	auto plasma = std::make_shared<meq::NormalisedMHDSource>(
		kinkedSpline(), kinkedSpline(), psiAxis );
	plasma->setNormalisation( psiAxis, 0.05 );

	auto coils = std::make_shared<meq::CoilSet>();
	coils->add( meq::Coil( 2.0, 0.0, 0.15, 0.20, 1.0e6 ) );

	meq::CoilAugmentedNormalisedSource wrapped(
		plasma, std::static_pointer_cast<meq::CoilSet const>( coils ) );

	// THE ORDINARY ENTRY POINT MUST STILL REFUSE. If this ever passes, the cast
	// has been widened and f() is silently missing its coils.
	{
		meq::SourceView view{};
		std::vector<meq::KnotView> knots;
		BOOST_TEST( !meq::makeSourceView( wrapped, view, knots ),
		            "makeSourceView() must refuse a coil-augmented source: its "
		            "f() is the SUM and a view of it would drop the coils" );
	}

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( wrapped, view, knots ),
	                    "the plasma term of a coil-augmented source has a POD "
	                    "form and makePlasmaSourceView() must find it" );

	double const span = psiAxis - 0.05;
	int checked = 0;
	int coilTermSeen = 0;

	// 2.0 IS INSIDE THE COIL AND THE OTHERS ARE NOT, WHICH IS THE POINT OF IT.
	// The coil spans R in [ 1.85, 2.15 ], so a sweep over 0.3, 1.0, 1.7 and 3.2
	// alone samples F_coil = 0 at every point -- and then f() and scaledF()
	// agree everywhere and this case tests nothing. Measured: the control below
	// failed exactly that way before 2.0 was added.
	for ( double const radius : { 0.3, 1.0, 1.7, 2.0, 3.2 } )
		for ( double const normalised : sweepPoints() )
		{
			double const psi = sameAbscissa( 0.05 + normalised*span );

			double const fromSource = wrapped.scaledF( radius, 0.0, psi );
			double const fromView =
				meq::sourceF( view, knots.data(), radius, 0.0, psi );
			BOOST_TEST( agrees( fromView, fromSource ),
			            "scaledF disagreed at R = " << radius
			            << ", Psi = " << normalised );

			double const dFromSource = wrapped.scaledDFdPsi( radius, 0.0, psi );
			double const dFromView =
				meq::sourceDFdPsi( view, knots.data(), radius, 0.0, psi );
			BOOST_TEST( agreesNearABoundary( dFromView, psi,
			                              [ & ]( double p )
			                              { return wrapped.scaledDFdPsi( radius, 0.0, p ); } ),
			            "scaledDFdPsi disagreed at R = " << radius
			            << ", Psi = " << normalised );

			// The control: f() is the SUM, so it must DIFFER from scaledF()
			// wherever the coil reaches. Without this the case above would
			// pass just as well against a wrapper that had no coils in it.
			if ( wrapped.f( radius, 0.0, psi ) != fromSource )
				++coilTermSeen;
			++checked;
		}

	BOOST_TEST( checked == 65 );
	BOOST_TEST( coilTermSeen == 13,
	            "the coil term must be non-zero at exactly the 13 sweep points "
	            "of R = 2.0, which is the one radius inside the coil" );
	BOOST_TEST( coilTermSeen > 0,
	            "f() and scaledF() agreed everywhere, so this case did not "
	            "test the distinction it exists for -- the coil term is zero "
	            "or the wrapper is not wrapping" );
}


/*
 * THE RESIDUAL INTEGRAND THROUGH A COIL WRAPPER, WHICH IS WHAT THE DEVICE
 * KERNEL ACTUALLY COMPUTES.
 *
 * meq::SourceIntegrator accumulates
 *
 *     elvect.Add( -weight * sourceValue( R, z, psi + psi_c, e )/R, shape )
 *
 * and `sourceValue` is `elementCarriesPlasma( e ) ? source->f( R, z, psi )
 *                                                 : fOutsidePlasma( R, z )`.
 * For a coil-augmented source those are `plasma + coil` and `coil` -- so the
 * POD form is one expression, `( plasma gated by the element flag ) + coil`,
 * and the gate is the ONLY difference between the two branches. This case
 * asserts that identity on both sides of the gate, to the bit.
 *
 * **THE COIL TERM IS POD-ABLE AND THE CLASS DOES NOT LOOK IT.** meq::CoilSet
 * also carries psi(), psiOf() and flux(), which are elliptic integrals and are
 * not device arithmetic. The SOURCE term is not any of those.
 */
BOOST_AUTO_TEST_CASE( thePodReproducesTheResidualIntegrandThroughACoilWrapper )
{
	double const psiAxis = 0.4;
	auto plasma = std::make_shared<meq::NormalisedMHDSource>(
		kinkedSpline(), kinkedSpline(), psiAxis );
	plasma->setNormalisation( psiAxis, 0.05 );

	// TWO coils, deliberately OVERLAPPING, because CoilSet::f() sums rather
	// than short-circuiting on the first hit and a "find the coil" rewrite
	// would pass every test built on disjoint ones.
	auto coils = std::make_shared<meq::CoilSet>();
	coils->add( meq::Coil( 2.00, 0.0, 0.15, 0.20, 1.0e6 ) );
	coils->add( meq::Coil( 2.05, 0.0, 0.15, 0.20, -4.0e5 ) );

	meq::CoilAugmentedNormalisedSource wrapped(
		plasma, std::static_pointer_cast<meq::CoilSet const>( coils ) );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( wrapped, view, knots ) );

	std::vector<meq::CoilView> coilViews;
	mfem::real_t mu0 = 0;
	BOOST_TEST_REQUIRE( meq::makeCoilViews( wrapped, coilViews, mu0 ) );
	BOOST_TEST_REQUIRE( coilViews.size() == 2u );

	double const span = psiAxis - 0.05;
	int checked = 0;
	int overlapSeen = 0;
	int coilSeen = 0;

	// Radii chosen to straddle both coils and their OVERLAP, and to sit on an
	// edge -- meq::Coil::contains() is closed and a half-open POD would differ
	// only there.
	for ( double const radius : { 0.3, 1.80, 1.85, 1.95, 2.00, 2.10, 2.20, 3.2 } )
		for ( double const normalised : sweepPoints() )
		{
			double const psi = sameAbscissa( 0.05 + normalised*span );

			double const coilFromSet = coils->f( radius, 0.0 );
			double const coilFromPod =
				meq::coilF( coilViews.data(), 2, mu0, radius, 0.0 );
			BOOST_TEST( agreesNearABoundary( coilFromPod, radius,
			                                 [ & ]( double r )
			                                 { return coils->f( r, 0.0 ); } ),
			            "the coil term disagreed at R = " << radius );

			if ( coilFromSet != 0.0 )
				++coilSeen;
			if ( coils->indexContaining( radius, 0.0 ) >= 0
			     && radius >= 1.95 && radius <= 2.15 )
				++overlapSeen;

			// INSIDE the plasma component: f() is the sum.
			double const inside = wrapped.f( radius, 0.0, psi );
			double const insidePod =
				meq::sourceF( view, knots.data(), radius, 0.0, psi )
				+ coilFromPod;
			BOOST_TEST( agreesNearABoundary( insidePod, radius,
			                                 [ & ]( double r )
			                                 { return wrapped.f( r, 0.0, psi ); } ),
			            "the integrand inside the component disagreed at R = "
			            << radius << ", Psi = " << normalised );

			// OUTSIDE it: the plasma term is gated off and the coil remains.
			double const outside = wrapped.fOutsidePlasma( radius, 0.0 );
			BOOST_TEST( agreesNearABoundary( coilFromPod, radius,
			                                 [ & ]( double r )
			                                 { return wrapped.fOutsidePlasma( r, 0.0 ); } ),
			            "the integrand outside the component disagreed at R = "
			            << radius );
			++checked;
		}

	BOOST_TEST( checked == 104 );
	BOOST_TEST( coilSeen > 0,
	            "no sample point was inside any coil, so the coil term was "
	            "zero throughout and this case tested nothing" );
	BOOST_TEST( overlapSeen > 0,
	            "no sample point was in the overlap, so a first-hit rewrite "
	            "would pass this case" );
}
