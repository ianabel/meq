#define BOOST_TEST_MODULE AxisAgreement
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/GradShafranov.hpp"

#include "analytic/HighBetaPoloidal.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * IS psi_ax A MAGNETIC AXIS, OR MERELY THE LARGEST NUMBER IN THE POTENTIAL
 * VECTOR? THE DEFINITION DOES NOT SAY, AND ON 2026-09-06 THAT PRODUCED A
 * COMPLETELY WRONG EQUILIBRIUM WITH EVERY DIAGNOSTIC GREEN.
 *
 * GradShafranovSolver::psiAxis() is the largest NODAL value of psi_h, and that
 * is deliberate: one nodal value is one entry of the discrete unknown, so the
 * bordered Newton's row is exactly -e_j. The constraint it closes,
 * G( lambda, s ) = s - max psi_h, is satisfied at machine zero by a spurious
 * nodal spike exactly as it is by an axis, so NOTHING the solver reports can
 * tell them apart -- and on a free-boundary machine case at k = 2 nothing did:
 * 17 Newton steps, psi_ax - max psi_h reading 0.000e+00, the prescribed plasma
 * current delivered to seven figures, and psi_ax = 2.734289e+00 against a peak
 * of 8.64e-02.
 *
 * meq::CriticalPointFinder::checkAxis() is the guard, and what it compares is
 * the NORMALISED FLUX at a genuine O-point -- a zero of q_h, which is a solved
 * field carrying the potential's own order rather than a derivative of one.
 * Psi at the magnetic axis is 1 by definition when psi_ax is the axis flux, so
 * the comparison is against a known number in the units the profiles actually
 * consume.
 *
 * WHAT THIS FILE ASSERTS IS THE GUARD, AND THE FIXTURES ARE CHEAP ON PURPOSE.
 * The failure that motivated it is a free-boundary solve with a PRESCRIBED
 * PLASMA CURRENT, an exterior coupling and a limiter point, and the first of
 * those does not exist here to be driven -- so the original case is not
 * reproducible in this tree at any price, and inventing a fixture that happened
 * to pass would be worse than saying so.
 *
 * What the guard actually needs to meet is a psi_h whose largest nodal value is
 * not its axis, and a SPIKED DOF is exactly that -- the shape of the observed
 * defect term for term: the potential carries a huge isolated nodal value while
 * q_h, a separately solved field, still has its zero where the plasma is. Two
 * rungs, in the ladder this tree uses everywhere:
 *
 *   * a PROJECTED pair ( psi, q ), no solver at all, where the axis is known in
 *     closed form and the guard's answer can be checked against arithmetic;
 *   * a REAL bordered-Newton solve of the high-beta source, which exercises
 *     CriticalPointFinder( solver ) -- the ctor the driver calls -- against the
 *     solver's own psiAxis().
 *
 * Both rungs are run healthy first, which is the half that says the guard does
 * not simply refuse everything.
 *
 * AND SINCE 2026-09-07 THE SOLVED RUNG IS PINNED TO THE OLD DEFINITION OF
 * psi_ax, WHICH IS WHAT KEEPS IT A TEST OF THE GUARD.
 *
 * GradShafranovSolver::AxisConstraint now defaults to LocatedAxis: psi_ax is
 * constrained AT a zero of q_h. Under that default Psi at the located axis is 1
 * BY CONSTRUCTION, so a guard case run on it would be checking a solve against
 * the formula it used -- and would pass while testing nothing. So
 * theSolversOwnAxisFluxIsCheckedAgainstAZeroOfTheFlux asks for
 * AxisConstraint::NodalMaximum explicitly, which is the configuration in which
 * psi_ax CAN be wrong and in which the spiked-dof half means anything.
 *
 * theLocatedAxisConstraintPutsPsiAxisOnTheAxis is the separate case for the new
 * default, and it asserts the things the construction does NOT give away: that
 * the constraint really was applied at a located axis rather than falling back
 * to the nodal maximum, and that the two definitions actually differ on this
 * field -- without which everything else about it would be vacuous.
 *
 * WHICH IS NOT TO SAY THE GUARD IS EMPTY UNDER THE NEW DEFAULT. It re-locates
 * the axis INDEPENDENTLY, by a full sweep taking the largest Psi, where the
 * solver's constraint follows one root from a warm start -- and the two land on
 * different elements' versions of it, an O( h^{k+1} ) apart. Measured here, the
 * solver constrains at ( 1.0918, 0.0003 ) and the guard reports ( 1.0913,
 * -0.0003 ), so Psi comes back at 0.99949 rather than at exactly 1.
 */

namespace
{
	using meq::analytic::HighBetaPoloidal;
	using meq::tests::NormalisedEquilibriumSource;
	using meq::tests::standardBox;

	/// The centre of the paraboloid, and its half-widths. Well inside
	/// standardBox() so that the axis is an interior maximum and the corners are
	/// comfortably negative.
	double const centreR = 1.0;
	double const centreZ = 0.0;
	double const halfR = 0.55;
	double const halfZ = 0.55;
	double const peak = 0.25;

	/// psi = peak ( 1 - ( ( r - R0 )/a )^2 - ( z/b )^2 ), a quadratic, so P_k
	/// represents it exactly at k >= 2 and the only error in the fixture is
	/// round-off. Its maximum is `peak` at ( R0, 0 ) and nowhere else.
	double paraboloid( mfem::Vector const &x )
	{
		double const u = ( x( 0 ) - centreR )/halfR;
		double const v = ( x( 1 ) - centreZ )/halfZ;
		return peak*( 1.0 - u*u - v*v );
	}

	/// q = ( 1/r ) grad_bar( psi ), in MEQ's sign convention -- what
	/// GradShafranovSolver::flux() returns, NOT the raw block, which holds -q
	/// and would silently turn every maximum into a minimum.
	void paraboloidFlux( mfem::Vector const &x, mfem::Vector &value )
	{
		value.SetSize( 2 );
		value( 0 ) = -2.0*peak*( x( 0 ) - centreR )/( halfR*halfR*x( 0 ) );
		value( 1 ) = -2.0*peak*( x( 1 ) - centreZ )/( halfZ*halfZ*x( 0 ) );
	}

	/// A field with no interior extremum at all: psi rising monotonically in r,
	/// which is the wall-hugging annulus branch a free-boundary solve can settle
	/// on. q never vanishes, so there is no axis to find.
	double monotone( mfem::Vector const &x )
	{
		return 0.1*x( 0 );
	}

	void monotoneFlux( mfem::Vector const &x, mfem::Vector &value )
	{
		value.SetSize( 2 );
		value( 0 ) = 0.1/x( 0 );
		value( 1 ) = 0.0;
	}

	/// The dof carrying the largest value of @a field, and the value there.
	/// Independent of CriticalPointFinder's own search on purpose: a test that
	/// located the spike with the routine under test would be checking a solve
	/// against the formula it used.
	int largestNodalDof( mfem::GridFunction const &field, double &value )
	{
		int best = -1;
		value = -std::numeric_limits<double>::infinity();
		for ( int i = 0; i < field.Size(); ++i )
			if ( field( i ) > value )
			{
				value = field( i );
				best = i;
			}
		return best;
	}

	/// A dof of @a field at least @a away metres from ( r, z ), so that a spike
	/// planted there is unambiguously somewhere else. Returns -1 if there is
	/// none, which no mesh in this file produces.
	int dofAwayFrom( mfem::GridFunction const &field, double r, double z,
	                 double away )
	{
		mfem::FiniteElementSpace const *space = field.FESpace();
		mfem::Mesh *mesh = space->GetMesh();
		mfem::Array<int> dofs;

		for ( int e = 0; e < mesh->GetNE(); ++e )
		{
			mfem::Vector centre;
			mesh->GetElementCenter( e, centre );
			double const dr = centre( 0 ) - r;
			double const dz = centre( 1 ) - z;
			if ( std::sqrt( dr*dr + dz*dz ) < away )
				continue;

			space->GetElementDofs( e, dofs );
			if ( dofs.Size() > 0 )
				return dofs[ 0 ] >= 0 ? dofs[ 0 ] : -1 - dofs[ 0 ];
		}
		return -1;
	}

	void report( char const *what, meq::AxisAgreement const &found )
	{
		std::printf( "    %-22s psi_ax = %12.6e at ( %7.4f, %7.4f )\n",
		             what, found.psiAxis, found.nodeR, found.nodeZ );
		if ( !found.located )
		{
			std::printf( "    %-22s NO O-point of that sense anywhere on the mesh "
			             "( %d saddles )\n", "", found.saddles );
			return;
		}
		std::printf( "    %-22s O-point psi = %12.6e at ( %7.4f, %7.4f ), "
		             "|q| = %8.2e\n",
		             "", found.axis.psi, found.axis.r, found.axis.z,
		             found.axis.fluxResidual );
		std::printf( "    %-22s Psi there = %10.4e   separation = %8.2e m = "
		             "%6.2f element diameters   %d extrema, %d saddles   %s\n",
		             "", found.normalisedFlux, found.separation,
		             found.separationInElements, found.extrema, found.saddles,
		             found.agrees ? "AGREES" : "DISAGREES" );
		std::fflush( stdout );
	}
}

/*
 * THE PROJECTED RUNG. A quadratic psi and its exact flux, so the magnetic axis
 * is ( 1.0, 0.0 ) by arithmetic and the guard's answer can be checked against
 * that rather than against another routine.
 *
 * Psi at the located axis reads slightly ABOVE 1 on the healthy column and that
 * is the correct behaviour rather than slack: psi_ax is the largest NODAL value
 * and the peak of the polynomial over a closed element is at least that, so the
 * healthy reading approaches 1 from above. The guard is one sided for exactly
 * this reason -- see meq::AxisAgreement.
 */
BOOST_AUTO_TEST_CASE( aSpikedNodalValueIsNotAMagneticAxis )
{
	int const order = 2;
	int const n = 12;

	mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
	mfem::L2_FECollection potentialCollection( order, mesh.Dimension(),
	                                           mfem::BasisType::GaussLobatto );
	mfem::L2_FECollection fluxCollection( order, mesh.Dimension(),
	                                      mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace potentialSpace( &mesh, &potentialCollection );
	mfem::FiniteElementSpace fluxSpace( &mesh, &fluxCollection, 2 );

	mfem::GridFunction potential( &potentialSpace );
	mfem::GridFunction flux( &fluxSpace );
	mfem::FunctionCoefficient exactPotential( paraboloid );
	mfem::VectorFunctionCoefficient exactFlux( 2, paraboloidFlux );
	potential.ProjectCoefficient( exactPotential );
	flux.ProjectCoefficient( exactFlux );

	meq::CriticalPointFinder finder( flux, potential );

	std::printf( "\n  a projected paraboloid, k = %d, n = %d: the axis is "
	             "( %.4f, %.4f ), psi there %.6e\n",
	             order, n, centreR, centreZ, peak );

	double nodal = 0.0;
	largestNodalDof( potential, nodal );
	meq::AxisAgreement const healthy = finder.checkAxis( nodal );
	report( "healthy", healthy );

	BOOST_TEST( healthy.located,
		"checkAxis() found no interior maximum of a paraboloid whose maximum is "
		"at ( " << centreR << ", " << centreZ << " ). The flux was projected from "
		"a closed form, so this is the root finder and not the field." );
	BOOST_TEST( std::abs( healthy.axis.r - centreR ) < 1.0e-8,
		"the located axis is at r = " << healthy.axis.r << " against an exact "
		<< centreR );
	BOOST_TEST( std::abs( healthy.axis.z - centreZ ) < 1.0e-8,
		"the located axis is at z = " << healthy.axis.z << " against an exact "
		<< centreZ );
	BOOST_TEST( healthy.agrees,
		"the guard fires on a field whose largest nodal value IS its axis: Psi "
		"read " << healthy.normalisedFlux << ". A guard that refuses a healthy "
		"run is worse than none." );
	BOOST_TEST( healthy.normalisedFlux >= 1.0,
		"Psi at the axis read " << healthy.normalisedFlux << ", below 1. psi_ax "
		"is the largest NODAL value and the peak of a polynomial over a closed "
		"element is at least its largest nodal value, so a healthy field reaches "
		"1 from above." );

	/*
	 * THE DEFECT, PLANTED: one dof of one element far from the axis raised to a
	 * multiple of the true peak, with the flux left alone. That is the observed
	 * failure's shape exactly -- psi_h carries a huge isolated nodal value while
	 * q_h still has its zero where the plasma is.
	 *
	 * The two multipliers are the two measured cases: 29x is the free-boundary
	 * machine run ( psi_ax = 2.734289e+00 against a peak of 8.64e-02, Psi 0.032 )
	 * and 3.2x is the half-disc that latched onto the corner where Gamma meets
	 * the axis ( Psi 0.31 ). A guard that caught only the first would be tuned to
	 * one number.
	 */
	for ( double multiplier : { 29.0, 3.2 } )
	{
		mfem::GridFunction spiked( potential );
		int const victim = dofAwayFrom( spiked, centreR, centreZ, 0.3 );
		BOOST_REQUIRE( victim >= 0 );
		spiked( victim ) = multiplier*peak;

		double spikedNodal = 0.0;
		int const found = largestNodalDof( spiked, spikedNodal );
		BOOST_REQUIRE_EQUAL( found, victim );

		meq::CriticalPointFinder spikedFinder( flux, spiked );
		meq::AxisAgreement const bad = spikedFinder.checkAxis( spikedNodal );
		std::printf( "    spiked %.1fx:\n", multiplier );
		report( "", bad );

		BOOST_TEST( bad.located,
			"the spike removed the axis from the search, which it must not: the "
			"flux was not touched." );
		BOOST_TEST( !bad.agrees,
			"a psi_ax " << multiplier << " times the true peak, planted on a "
			"single dof away from the axis, was ACCEPTED at Psi = "
			<< bad.normalisedFlux << ". This is the whole of what the guard "
			"exists to catch." );
		BOOST_TEST( std::abs( bad.normalisedFlux - 1.0/multiplier ) < 1.0e-6,
			"Psi read " << bad.normalisedFlux << " against the 1/" << multiplier
			<< " = " << 1.0/multiplier << " arithmetic demands. The guard is "
			"reporting something other than psi( O-point )/psi_ax." );
		BOOST_TEST( bad.separationInElements > 2.0,
			"the spike was planted at least 0.3 m from the axis and the "
			"separation reads only " << bad.separationInElements << " element "
			"diameters, so the reported node is not where the spike is." );
	}
}

/*
 * A CONDUCTOR'S O-POINT IS AN EXTREMUM OF psi AND WINS EVERY SCORE THIS GUARD
 * KNOWS, AND ON A MACHINE WHOSE COILS ARE MESHED INSIDE Omega IT IS THERE TO BE
 * FOUND.
 *
 * Every machine MEQ can pose has them inside: no semicircle centred on the axis
 * both encloses a plasma and excludes the coils near it, so the conductors are
 * in the domain and enter through the SOURCE. Any one of them carrying current
 * of the plasma's own sign contributes a local maximum of psi, and one carrying
 * a large enough fraction of I_p contributes the LARGEST one. checkAxis()
 * competes candidates on normalised flux, so it hands that back -- correctly,
 * and uselessly, because a plasma has no magnetic axis inside a conductor.
 *
 * THE FIXTURE IS TWO GAUSSIAN BUMPS WITH AN EXACT FLUX, so nothing here depends
 * on a solve: the plasma's peaks at 0.25 at ( 1.00, 0.00 ) and the conductor's
 * at 0.40 at ( 1.60, 0.30 ). That ordering is the observed one on the diverted
 * machine of FREE-BOUNDARY-PLAN.md section 10, where P1L carries +1.37e+05 A
 * against an I_p of 2.0e+05 A and its O-point reads psi = 1.21e-01 against the
 * reference equilibrium's axis at 8.3e-02.
 *
 * IT ASSERTS BOTH HALVES AND THE FIRST IS WHAT KEEPS THE SECOND HONEST:
 * unfiltered, the guard must return the CONDUCTOR. If it did not, the exclusion
 * would be untested whatever the filtered half reported.
 */
namespace
{
	double const plasmaR = 1.00;
	double const plasmaZ = 0.00;
	double const plasmaPeak = 0.25;
	double const bumpWidth = 0.22;

	double const conductorR = 1.60;
	double const conductorZ = 0.30;
	double const conductorPeak = 0.40;
	double const conductorHalf = 0.06;

	double twoBumps( mfem::Vector const &x )
	{
		double const a = ( x( 0 ) - plasmaR )/bumpWidth;
		double const b = ( x( 1 ) - plasmaZ )/bumpWidth;
		double const c = ( x( 0 ) - conductorR )/bumpWidth;
		double const d = ( x( 1 ) - conductorZ )/bumpWidth;
		return plasmaPeak*std::exp( -( a*a + b*b ) )
		       + conductorPeak*std::exp( -( c*c + d*d ) );
	}

	/// q = ( 1/r ) grad_bar( psi ), the same convention paraboloidFlux() uses.
	void twoBumpsFlux( mfem::Vector const &x, mfem::Vector &value )
	{
		value.SetSize( 2 );
		double const w = bumpWidth*bumpWidth;
		double const a = ( x( 0 ) - plasmaR )/bumpWidth;
		double const b = ( x( 1 ) - plasmaZ )/bumpWidth;
		double const c = ( x( 0 ) - conductorR )/bumpWidth;
		double const d = ( x( 1 ) - conductorZ )/bumpWidth;
		double const one = plasmaPeak*std::exp( -( a*a + b*b ) );
		double const two = conductorPeak*std::exp( -( c*c + d*d ) );
		value( 0 ) = ( -2.0*one*( x( 0 ) - plasmaR )/w
		               - 2.0*two*( x( 0 ) - conductorR )/w )/x( 0 );
		value( 1 ) = ( -2.0*one*( x( 1 ) - plasmaZ )/w
		               - 2.0*two*( x( 1 ) - conductorZ )/w )/x( 0 );
	}
}

BOOST_AUTO_TEST_CASE( aConductorsOwnOPointIsNotAMagneticAxis )
{
	int const order = 3;
	int const n = 24;
	meq::tests::Rectangle const box{ 0.60, 2.00, -0.60, 0.90 };

	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	mfem::L2_FECollection potentialCollection( order, mesh.Dimension(),
	                                           mfem::BasisType::GaussLobatto );
	mfem::L2_FECollection fluxCollection( order, mesh.Dimension(),
	                                      mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace potentialSpace( &mesh, &potentialCollection );
	mfem::FiniteElementSpace fluxSpace( &mesh, &fluxCollection, 2 );

	mfem::GridFunction potential( &potentialSpace );
	mfem::GridFunction flux( &fluxSpace );
	mfem::FunctionCoefficient exactPotential( twoBumps );
	mfem::VectorFunctionCoefficient exactFlux( 2, twoBumpsFlux );
	potential.ProjectCoefficient( exactPotential );
	flux.ProjectCoefficient( exactFlux );

	// The conductor as a machine would carry it: the predicate the driver
	// installs is a lambda over CoilSet::indexContaining(), so this uses the
	// same call rather than an open-coded rectangle.
	meq::CoilSet coils;
	coils.add( meq::Coil( conductorR, conductorZ, conductorHalf,
	                      conductorHalf, 1.0 ) );

	double nodal = 0.0;
	largestNodalDof( potential, nodal );

	std::printf( "\n  two bumps, k = %d, n = %d: the plasma at ( %.2f, %.2f ) "
	             "peaking %.2f, the conductor at ( %.2f, %.2f ) peaking %.2f\n",
	             order, n, plasmaR, plasmaZ, plasmaPeak,
	             conductorR, conductorZ, conductorPeak );

	meq::CriticalPointFinder open( flux, potential );
	meq::AxisAgreement const unfiltered = open.checkAxis( nodal );
	report( "unfiltered", unfiltered );

	BOOST_TEST_REQUIRE( unfiltered.located,
		"the unfiltered guard found no interior extremum at all on a field that "
		"has two by construction, so neither half of this case means anything." );
	BOOST_TEST_REQUIRE( coils.indexContaining( unfiltered.axis.r,
	                                           unfiltered.axis.z ) >= 0,
		"the unfiltered guard returned ( " << unfiltered.axis.r << ", "
		<< unfiltered.axis.z << " ), which is NOT inside the conductor. The "
		"filtered half asserts that the exclusion MOVES the answer, and that can "
		"only mean something if the answer starts on the conductor -- so this is "
		"the fixture's premise rather than a claim about the code." );

	// AND NOW THE SAME FIELD WITH THE CONDUCTOR EXCLUDED.
	meq::CriticalPointFinder filtered( flux, potential );
	filtered.setExcluded(
		[ &coils ]( double r, double z )
		{ return coils.indexContaining( r, z ) >= 0; } );
	meq::AxisAgreement const guarded = filtered.checkAxis( nodal );
	report( "conductor excluded", guarded );

	BOOST_TEST_REQUIRE( guarded.located,
		"the exclusion took the conductor's O-point and left the guard with "
		"nothing, on a field carrying a second extremum 0.6 m away. "
		"setExcluded() filters WHICH candidate is competed and must not empty "
		"the competition." );
	BOOST_TEST( std::abs( guarded.axis.r - plasmaR ) < 2.0e-2,
		"with the conductor excluded the guard returned r = " << guarded.axis.r
		<< " against the plasma bump's " << plasmaR );
	BOOST_TEST( std::abs( guarded.axis.z - plasmaZ ) < 2.0e-2,
		"with the conductor excluded the guard returned z = " << guarded.axis.z
		<< " against the plasma bump's " << plasmaZ );
	BOOST_TEST( std::abs( guarded.axis.psi - plasmaPeak ) < 1.0e-3,
		"the excluded guard's axis carries psi = " << guarded.axis.psi
		<< " against the plasma bump's peak of " << plasmaPeak << ". A value "
		"near " << conductorPeak << " would mean it reached the conductor by "
		"another route." );

	// AND THE EXCLUDED POINT IS NOT COUNTED EITHER, which is what setExcluded()
	// documents: the count is of CANDIDATES, and a point that cannot be an axis
	// is not one. Without it the diagnostic would report a competition it did
	// not hold.
	BOOST_TEST( guarded.extrema == unfiltered.extrema - 1,
		"the unfiltered sweep counted " << unfiltered.extrema
		<< " extrema and the filtered one " << guarded.extrema
		<< ", where exactly one candidate was excluded." );
}

/*
 * NO INTERIOR EXTREMUM IS A RESULT, NOT AN ABSENCE OF ONE. A monotone psi is the
 * wall-hugging annulus branch: every constraint a bordered Newton imposes can be
 * satisfied by it -- they constrain the current and the normalisations, and none
 * says the plasma is a core -- so it is a converged answer with no magnetic axis
 * in it, and psi_ax is then the edge of nothing.
 */
BOOST_AUTO_TEST_CASE( aFieldWithNoInteriorExtremumHasNoAxisToAgreeWith )
{
	int const order = 2;
	int const n = 8;

	mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
	mfem::L2_FECollection collection( order, mesh.Dimension(),
	                                  mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace potentialSpace( &mesh, &collection );
	mfem::FiniteElementSpace fluxSpace( &mesh, &collection, 2 );

	mfem::GridFunction potential( &potentialSpace );
	mfem::GridFunction flux( &fluxSpace );
	mfem::FunctionCoefficient exactPotential( monotone );
	mfem::VectorFunctionCoefficient exactFlux( 2, monotoneFlux );
	potential.ProjectCoefficient( exactPotential );
	flux.ProjectCoefficient( exactFlux );

	double nodal = 0.0;
	largestNodalDof( potential, nodal );

	meq::CriticalPointFinder finder( flux, potential );
	meq::AxisAgreement const found = finder.checkAxis( nodal );

	std::printf( "\n  a monotone psi, the annulus branch:\n" );
	report( "", found );

	BOOST_TEST( !found.located,
		"an O-point was reported on a field that rises monotonically in r, where "
		"q never vanishes. sweep() found " << found.extrema << " extrema." );
	BOOST_TEST( !found.agrees,
		"the guard accepted a psi_ax on a field with no magnetic axis at all." );
}

/*
 * THE SOLVED RUNG, THROUGH THE CTOR THE DRIVER CALLS. The high-beta source with
 * psi_ax an unknown of the bordered Newton -- so psiAxis() is a genuine solved
 * quantity here rather than an argument, and CriticalPointFinder( solver ) takes
 * flux() and potential() itself, which is where handing it the raw block instead
 * would turn the maximum into a minimum.
 *
 * AND IT IS DRIVEN AT AxisConstraint::NodalMaximum DELIBERATELY, WHICH IS NOT
 * THE SOLVER'S DEFAULT SINCE 2026-09-07.
 *
 * The default is now AxisConstraint::LocatedAxis: psi_ax is constrained AT a
 * zero of q_h. Under it, Psi at the located axis is 1 by construction, so
 * checkAxis().agrees is very nearly true whatever the field does -- the guard
 * would be checking a solve against the formula it used, which is this tree's
 * own recorded trap. A case that passed for that reason while looking like a
 * test of the guard would be worse than no case.
 *
 * So this one pins the OLD definition, where psi_ax CAN be something that is not
 * an axis, and it is the configuration the guard exists for: the spiked-dof half
 * below depends on psi_ax being the largest nodal value, since that is exactly
 * what the defect delivers. theLocatedAxisConstraintPutsPsiAxisOnTheAxis is the
 * separate case for the new default, and it asserts a different thing.
 */
BOOST_AUTO_TEST_CASE( theSolversOwnAxisFluxIsCheckedAgainstAZeroOfTheFlux )
{
	int const order = 2;
	int const n = 8;
	int const nu = 2;
	double const amplitude = 1.0;

	// The dimensional estimate sqrt( nu A / lambda_1 ) on standardBox(), which is
	// what HighBetaConvergence seeds this source with. It is a starting point for
	// the border and for the bump, not an answer.
	meq::tests::Rectangle const box = standardBox();
	double const width = box.rMax - box.rMin;
	double const height = box.zMax - box.zMin;
	double const lambda = M_PI*M_PI*( 1.0/( width*width ) + 1.0/( height*height ) );
	double const estimate = std::sqrt( nu*amplitude/lambda );

	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	HighBetaPoloidal equilibrium =
		HighBetaPoloidal::peaked( nu, amplitude, estimate );
	NormalisedEquilibriumSource<HighBetaPoloidal> source( equilibrium );

	mfem::ConstantCoefficient zero( 0.0 );
	// A separable sine bump of about the right height. The trivial branch and a
	// small-amplitude second solution are both in reach from the Dirichlet datum
	// on this source, so the guess is part of the problem statement.
	double const rMin = box.rMin;
	double const zMin = box.zMin;
	mfem::FunctionCoefficient guess(
		[ estimate, rMin, zMin, width, height ]( mfem::Vector const &x )
		{
			return estimate*std::sin( M_PI*( x( 0 ) - rMin )/width )
			       *std::sin( M_PI*( x( 1 ) - zMin )/height );
		} );

	meq::GradShafranovSolver solver( mesh, order );
	// THE OLD DEFINITION, ON PURPOSE. See the header: under the default the
	// guard's verdict is nearly self-referential, and the spike below needs
	// psi_ax to BE the largest nodal value.
	solver.setAxisConstraint(
		meq::GradShafranovSolver::AxisConstraint::NodalMaximum );
	solver.setSource( source, estimate );
	solver.setBoundaryData( zero );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 40 );
	solver.solve();

	// The bordered path drives || ( R, gamma G ) ||; the constraint itself is
	// reported separately because the two are in different units. Both are
	// required before the guard's answer means anything.
	BOOST_REQUIRE( !solver.newtonResiduals().empty() );
	BOOST_REQUIRE( solver.newtonResiduals().back() < 1.0e-8 );
	BOOST_REQUIRE( std::abs( solver.normalisationResidual() ) < 1.0e-10 );

	std::printf( "\n  the high-beta bordered Newton, k = %d, n = %d, "
	             "nu = %d, A = %.1f\n", order, n, nu, amplitude );

	meq::CriticalPointFinder finder( solver );
	meq::AxisAgreement const healthy =
		finder.checkAxis( solver.psiAxis(), solver.psiBoundary() );
	report( "solved", healthy );

	BOOST_TEST( healthy.located,
		"no interior maximum of q_h was found on a converged high-beta solve, "
		"whose psi is a single hump vanishing on the box. If this fails on the "
		"CriticalPointFinder( solver ) ctor and not on the two-field one, the "
		"suspect is flux() against the raw block: the raw one holds -q, and in "
		"even dimension that turns every Maximum into a Minimum silently." );
	BOOST_TEST( healthy.agrees,
		"the guard refuses a converged bordered-Newton solve at Psi = "
		<< healthy.normalisedFlux << ", where the constraint psi_ax - max psi_h "
		"reads " << solver.normalisationResidual() << "." );
	BOOST_TEST( std::abs( healthy.nodalExtreme - solver.psiAxis() )
	            <= 1.0e-10*std::abs( solver.psiAxis() ),
		"checkAxis() recomputed the largest nodal value as "
		<< healthy.nodalExtreme << " where the solver reports psi_ax = "
		<< solver.psiAxis() << ". Same field, same rule, so a disagreement is "
		"itself a finding." );
	BOOST_TEST( healthy.separationInElements < 3.0,
		"the largest nodal value sits " << healthy.separationInElements
		<< " element diameters from the axis on a HEALTHY solve. It should be a "
		"node of the axis element or of a near neighbour." );

	/*
	 * The same solve with the potential spiked and psi_ax set to the spike, which
	 * is what the defect delivers: the constraint is still satisfied to machine
	 * zero, because psi_ax IS the largest nodal value. Every number the solver
	 * prints is unchanged.
	 */
	mfem::GridFunction &potential = solver.potential();
	int const victim = dofAwayFrom( potential, healthy.axis.r, healthy.axis.z,
	                                0.3 );
	BOOST_REQUIRE( victim >= 0 );
	double const spike = 29.0*solver.psiAxis();
	potential( victim ) = spike;

	meq::CriticalPointFinder spikedFinder( solver );
	meq::AxisAgreement const bad = spikedFinder.checkAxis( spike,
	                                                      solver.psiBoundary() );
	std::printf( "    spiked 29x:\n" );
	report( "", bad );

	BOOST_TEST( !bad.agrees,
		"a psi_ax 29 times the axis flux, planted on one dof of a converged "
		"solve, was ACCEPTED at Psi = " << bad.normalisedFlux << "." );
	BOOST_TEST( bad.located,
		"the spike removed the axis from the search. The roots are a property of "
		"q_h alone, which was not touched." );
	BOOST_TEST( std::abs( bad.axis.psi - healthy.axis.psi )
	            <= 1.0e-12*std::abs( healthy.axis.psi ),
		"the located O-point's psi moved from " << healthy.axis.psi << " to "
		<< bad.axis.psi << " when the POTENTIAL was spiked, which it can only do "
		"if the spiked element is the axis element." );
}

/*
 * THE NEW DEFAULT: psi_ax IS THE FLUX AT THE LOCATED AXIS, SO Psi THERE IS 1.
 *
 * AxisConstraint::LocatedAxis constrains psi_ax at a zero of q_h rather than at
 * the largest nodal value of psi_h -- FREE-BOUNDARY-PLAN.md section 11.5,
 * option 3. What that buys is that psi_ax MEANS the flux at a magnetic axis
 * rather than merely being the largest number in the potential vector.
 *
 * WHAT THIS CASE CAN AND CANNOT ASSERT, SAID PLAINLY BECAUSE THE DISTINCTION IS
 * THE WHOLE POINT OF SPLITTING IT FROM THE ONE ABOVE. Psi = 1 at the located
 * axis is very nearly TRUE BY CONSTRUCTION here, so it is NOT evidence that the
 * guard works -- theSolversOwnAxisFluxIsCheckedAgainstAZeroOfTheFlux is where
 * that is measured, at the old definition, where psi_ax can be wrong. What this
 * case is for is the two things the construction does NOT give for free:
 *
 *   * that the constraint was actually applied at a located axis rather than
 *     falling back to the nodal maximum, which axisWasLocated() reports and
 *     which a silent fallback would hide;
 *   * that the two definitions DISAGREE on this field, which is what says the
 *     choice is not cosmetic. If they agreed to round-off the constraint would
 *     be doing nothing and Psi = 1 would prove nothing.
 *
 * The gap is small here because nu = 2 is a mild source on a resolved mesh;
 * HighBetaConvergence's header measures the same gap at eight per cent on a
 * stiff one, and records that the two constraints then converge to genuinely
 * different equilibria.
 */
BOOST_AUTO_TEST_CASE( theLocatedAxisConstraintPutsPsiAxisOnTheAxis )
{
	int const order = 2;
	int const n = 8;
	int const nu = 2;
	double const amplitude = 1.0;

	meq::tests::Rectangle const box = standardBox();
	double const width = box.rMax - box.rMin;
	double const height = box.zMax - box.zMin;
	double const lambda = M_PI*M_PI*( 1.0/( width*width ) + 1.0/( height*height ) );
	double const estimate = std::sqrt( nu*amplitude/lambda );

	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	HighBetaPoloidal equilibrium =
		HighBetaPoloidal::peaked( nu, amplitude, estimate );
	NormalisedEquilibriumSource<HighBetaPoloidal> source( equilibrium );

	mfem::ConstantCoefficient zero( 0.0 );
	double const rMin = box.rMin;
	double const zMin = box.zMin;
	mfem::FunctionCoefficient guess(
		[ estimate, rMin, zMin, width, height ]( mfem::Vector const &x )
		{
			return estimate*std::sin( M_PI*( x( 0 ) - rMin )/width )
			       *std::sin( M_PI*( x( 1 ) - zMin )/height );
		} );

	meq::GradShafranovSolver solver( mesh, order );
	// Explicit rather than relied upon: this case is ABOUT the choice, so a
	// change of default must not silently turn it into a copy of the one above.
	solver.setAxisConstraint(
		meq::GradShafranovSolver::AxisConstraint::LocatedAxis );
	solver.setSource( source, estimate );
	solver.setBoundaryData( zero );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 40 );
	solver.solve();

	BOOST_REQUIRE( !solver.newtonResiduals().empty() );
	BOOST_REQUIRE( solver.newtonResiduals().back() < 1.0e-8 );

	double nodal = 0.0;
	largestNodalDof( solver.potential(), nodal );

	meq::CriticalPointFinder finder( solver );
	meq::AxisAgreement const found =
		finder.checkAxis( solver.psiAxis(), solver.psiBoundary() );

	std::printf( "\n  the located-axis constraint, k = %d, n = %d, nu = %d, "
	             "A = %.1f\n", order, n, nu, amplitude );
	std::printf( "    psi_ax %13.6e at ( %.4f, %.4f ), largest nodal value "
	             "%13.6e, gap %.2e relative\n",
	             solver.psiAxis(), solver.axisR(), solver.axisZ(), nodal,
	             std::abs( solver.psiAxis() - nodal )
	                 /std::max( std::abs( solver.psiAxis() ), 1.0e-300 ) );
	report( "solved", found );

	// IT WAS CONSTRAINED WHERE IT SAYS. Without this every assertion below would
	// pass just as well on a silent fallback to the nodal maximum.
	BOOST_TEST( solver.axisWasLocated(),
		"psi_ax fell back to the largest nodal value: no O-point of q_h was "
		"reachable on a single-hump source over a box, where there is one." );

	// The border still closes, at the new point. 1e-9 rather than round-off
	// because the constraint point comes from a root search.
	BOOST_TEST( std::abs( solver.normalisationResidual() )
	            < 1.0e-9*std::abs( solver.psiAxis() ),
		"the normalisation constraint reads " << solver.normalisationResidual()
		<< " against psi_ax = " << solver.psiAxis() );

	BOOST_TEST( found.located,
		"no O-point of q_h was found by checkAxis() on a field the solver just "
		"constrained AT one." );

	// Psi = 1 at the axis. Nearly true by construction -- see the header -- so
	// this is a consistency check between the constraint and the guard, not
	// evidence that the guard works.
	BOOST_TEST( found.agrees,
		"the guard refuses a solve whose psi_ax was constrained AT a zero of "
		"q_h, reading Psi = " << found.normalisedFlux << " where 1 is expected. "
		"The two disagree about WHICH O-point, or about the field." );

	// AND THE CHOICE IS NOT COSMETIC. If the located axis and the nodal maximum
	// agreed to round-off, the constraint would be doing nothing and everything
	// above would be vacuous. The gap is O( h^2 ) and small on this mild source;
	// what is ruled out is that it is zero.
	double const gap = std::abs( solver.psiAxis() - nodal )
	                   /std::max( std::abs( solver.psiAxis() ), 1.0e-300 );
	BOOST_TEST( gap > 1.0e-6,
		"psi_ax and the largest nodal value agree to " << gap << " relative, so "
		"this case cannot distinguish the two constraints and proves nothing "
		"about either." );
	BOOST_TEST( gap < 1.0e-2,
		"psi_ax and the largest nodal value differ by " << gap << " relative on "
		"a MILD source at nu = 2, where the gap should be the O( h^2 ) distance "
		"between two definitions. A gap this large means the located search is "
		"following something that is not the peak." );
}
