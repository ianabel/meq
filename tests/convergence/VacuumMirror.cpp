/*
 * MR-3: THE VACUUM MIRROR, THROUGH MEQ. MIRROR-PLAN.md, stage MR-3.
 *
 * No plasma. geq's two-pack mirror machine, a half-disc made by
 * tools/mesh/halfdisc.py, and the exterior coupling live -- the rung between
 * the analytic fixtures and a full machine solve, and the one that validates
 * the geometry, the mesh, the coil model and the exterior TOGETHER and without
 * the plasma. Three acceptances, one per case:
 *
 *   (a) psi against meq::coilPsi() -- the exact vacuum field inside Gamma -- at
 *       the rate the discretisation promises, k+1, over an h-ladder at three
 *       orders;
 *   (b) the mode count N, swept on MEQ's own solve until the error stops
 *       falling, which is the number MR-0 bounded and this stage measures;
 *   (c) psi against geq's own vacuum field, with the disagreement ATTRIBUTED by
 *       measuring geq's field against meq::coilPsi() directly.
 *
 * The machine, the mesh, D_h and both vacuum routes are MirrorMachine.hpp's,
 * shared with RotatingMirror; its comment says why the coils are inside Gamma
 * and why Gamma is at 3.0.
 *
 * THIS FILE TAKES THE MESHED ROUTE FOR (a) AND (b), because a subtracted
 * vacuum has psi_p = 0 identically and nothing to converge --
 * COIL-SUBTRACTION-PLAN.md section 10.1 says a rate there would be wrong in
 * both directions -- while a meshed one makes the coil field the solve's own
 * answer, which is the claim. Subtraction appears in (c), where it is the
 * instrument rather than the claim.
 */

#define BOOST_TEST_MODULE VacuumMirror
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "RoundOffFloor.hpp"

#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/Sampler.hpp"

#include "convergence/ConvergenceHarness.hpp"
#include "convergence/GeqReference.hpp"
#include "convergence/MirrorMachine.hpp"

namespace
{
	using namespace meqtest::mirror;

	/// meshFile()'s stem for this test.
	char const *const meshStem = "VacuumMirror";

	std::unique_ptr<Domain> buildDomain( int refinements )
	{
		return meqtest::mirror::buildDomain( refinements, meshStem );
	}

	/**
	 * The relative tail of psi_coil's own trace on Gamma past N retained
	 * modes, in the operator's weight dGamma/R = dmu/( 1 - mu^2 ) -- what M-169
	 * computed for the plasma, computed for the copper. It is the best any
	 * N-mode coupling can do ON Gamma, and entry N is the tail after degrees
	 * 2 .. N + 1.
	 *
	 * **ITS OWN 400-POINT RULE IN mu, AND NOT ExteriorDtN::coefficients()**,
	 * so that the reference this case reads N against shares no quadrature
	 * with the class under test. A fixed 30-point rule in coefficients() read
	 * this tail as a flat 0.82 out to N = 48 at 96 modes -- M-260, the
	 * quadrature and not the trace -- which is why that rule now scales with
	 * the mode count; this one resolves the highest mode with room to spare.
	 * The tail is summed DOWN from the last mode rather than subtracted from
	 * the norm, which would floor it at the square root of round-off.
	 */
	std::vector<double> traceTail( meq::CoilSet const &coils, int maxModes )
	{
		meq::ExteriorDtN const wide( 0.0, rhoGamma, maxModes );
		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 799 );

		std::vector<double> a( static_cast<std::size_t>( maxModes ), 0.0 );
		double norm2 = 0.0;
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			double const mu = 2.0*ip.x - 1.0;
			double const w = 2.0*ip.weight/( ( 1.0 - mu )*( 1.0 + mu ) );
			double const radius = rhoGamma*std::sqrt( ( 1.0 - mu )*( 1.0 + mu ) );
			double const z = rhoGamma*mu;
			double const psi = coils.psi( radius, z );
			norm2 += w*psi*psi;
			for ( int i = 0; i < maxModes; ++i )
				a[ static_cast<std::size_t>( i ) ] +=
					w*psi*wide.basis( wide.firstMode() + i, radius, z );
		}

		std::vector<double> tail( static_cast<std::size_t>( maxModes ) + 1, 0.0 );
		double rest = 0.0;
		for ( int i = maxModes; i >= 0; --i )
		{
			tail[ static_cast<std::size_t>( i ) ] = std::sqrt( rest/norm2 );
			if ( i > 0 )
			{
				int const n = wide.firstMode() + i - 1;
				double const coefficient = a[ static_cast<std::size_t>( i - 1 ) ]
				                           /wide.mass( n );
				rest += coefficient*coefficient*wide.mass( n );
			}
		}
		return tail;
	}

	/// geq's vacuum field on its own grid, R z psi per row.
	struct GeqPoint
	{
		double radius;
		double z;
		double psi;
	};

	std::vector<GeqPoint> readGeqPoints()
	{
		std::string const path = answerDirectory() + "/points.dat";
		std::ifstream in = meqtest::geq::openOutput( answerDirectory(), "points.dat" );
		std::vector<GeqPoint> points;
		std::string line;
		while ( std::getline( in, line ) )
		{
			if ( line.empty() || line[ 0 ] == '#' )
				continue;
			std::istringstream row( line );
			GeqPoint p{};
			BOOST_TEST_REQUIRE( static_cast<bool>( row >> p.radius >> p.z >> p.psi ),
			                    "malformed row in " << path << ": " << line );
			points.push_back( p );
		}
		BOOST_TEST_REQUIRE( !points.empty(), path << " has no points" );
		return points;
	}

	double rms( std::vector<double> const &v )
	{
		double sum = 0.0;
		for ( double x : v )
			sum += x*x;
		return std::sqrt( sum/static_cast<double>( v.size() ) );
	}
}


/*
 * (a) MEQ'S psi AGAINST meq::coilPsi(), AT k+1.
 *
 * With no plasma and every conductor inside Gamma, the exact answer is the
 * packs' own field: Delta* psi_coil = -mu0 R J inside the copper and 0 outside,
 * psi_coil = 0 on the axis, and its exterior a decaying Gegenbauer series that
 * the coupling represents to the mode truncation. So psi_h converges on
 * coilPsi() and the order is the discretisation's.
 *
 * TWO ROUTES TO ONE FIELD, AND THAT IS WHAT MAKES IT A TEST. The solve sees the
 * packs only as a current density, CoilSet::f() = mu0 R I/area, sampled at
 * quadrature points of the elements meshed to them; the reference integrates
 * the filament Green's function over each rectangle. A factor, a sign or a
 * misplaced conductor in either is a disagreement that does not fall with h:
 * each of M-260's four mutations -- one pack moved by its own half-width, one
 * pack's current reversed, 0.1% on CoilSet::f(), 1% on the DtN's block entry --
 * turns every order red here, the smallest of them flattening the k = 3 rate
 * to 0.03.
 *
 * WHERE IT IS ASSERTED, AND WHY TWO NORMS. Over the MACHINE REGION -- geq's own
 * box, R <= 0.63 and |z| <= 1.5, where MR-4's plasma lives and 0.15 m clear of
 * the copper -- with the two-tier rule CLAUDE_HDGGS.md gives for the extension
 * path: 0.30 of slack on every pair and 0.15 over the sequence. Over the WHOLE
 * of D_h on the finest pair only: the background is 0.4 m out to the arc, so
 * level 0 resolves the vacuum near Gamma with three elements across a metre,
 * and k = 1's first pair there reads 1.27 before reaching 2.14. Both norms reach
 * k+1 at k = 3 as well, which says the r^2 log r at the packs' corners -- which
 * caps a cut conductor's rate -- is not limiting a conductor meshed to.
 */
BOOST_AUTO_TEST_CASE( theVacuumMirrorConvergesOnTheCoilFieldAtKPlusOne )
{
	std::shared_ptr<meq::CoilSet> const reference = machineCoils();

	// THE GEOMETRY DECISION, ASKED OF THE TYPES THAT OWN IT. Every pack wholly
	// inside Gamma, measured to its FARTHEST point -- the precondition the
	// subtracted route checks and the meshed route needs just as much, since
	// a source beyond Gamma is a current the exterior series cannot hold. And
	// at M-169's rho_Gamma = 1.7 the copper straddles Gamma, so NEITHER route
	// admits it: that radius was chosen for the plasma alone.
	meq::ConductorField copper( machine().mu0 );
	meq::ExteriorCoilSet outside( machine().mu0 );
	for ( double sign : { -1.0, 1.0 } )
	{
		copper.add( pack( sign ) );
		outside.add( pack( sign ) );
	}
	double const containment = copper.containment( 0.0, rhoGamma );
	BOOST_TEST_REQUIRE( containment > 0.0,
		"a pack reaches beyond Gamma = " << rhoGamma << " by " << -containment
		<< " m, so its field is not the exterior series' to represent" );
	BOOST_TEST( copper.containment( 0.0, 1.7 ) < 0.0 );
	BOOST_TEST( outside.clearance( 0.0, 1.7 ) < 0.0 );

	std::printf( "\n  MR-3 (a): geq's mirror, no plasma, packs MESHED inside "
	             "Gamma = %.2f ( %.4f m clear ), N = %d\n",
	             rhoGamma, containment, defaultModes );
	std::printf( "    pack R [%.6f, %.6f], z +-[%.3f, %.3f], I = %.6e A each\n\n",
	             machine().innerRadius, machine().outerRadius,
	             machine().centreZ - 0.5*machine().length,
	             machine().centreZ + 0.5*machine().length, machine().packCurrent );
	std::printf( "    %2s %5s %7s %8s %13s %7s %13s %7s\n", "k", "level", "NE",
	             "h", "L2 D_h", "rate", "L2 machine", "rate" );

	int const levels = 3;
	int const orders = 3;
	std::vector<std::vector<double>> global( orders );
	std::vector<std::vector<double>> region( orders );
	double normGlobal = 0.0;
	double normRegion = 0.0;
	double controlRatio = 0.0;

	for ( int level = 0; level < levels; ++level )
	{
		std::unique_ptr<Domain> const d = buildDomain( level );
		for ( int order = 1; order <= orders; ++order )
		{
			std::unique_ptr<VacuumSolve> const s =
				solveMeshed( *d, order, defaultModes, machineCoils() );
			ReferenceTable const table( *d->sub, order, *reference );
			Errors const e = table.errors( s->solver->potential() );
			std::size_t const k = static_cast<std::size_t>( order - 1 );
			global[ k ].push_back( e.global );
			region[ k ].push_back( e.region );
			normGlobal = e.globalNorm;
			normRegion = e.regionNorm;

			double gRate = std::numeric_limits<double>::quiet_NaN();
			double rRate = std::numeric_limits<double>::quiet_NaN();
			if ( level > 0 )
			{
				gRate = meq::tests::rate( global[ k ][ level - 1 ], e.global, 2.0 );
				rRate = meq::tests::rate( region[ k ][ level - 1 ], e.region, 2.0 );
			}
			std::printf( "    %2d %5d %7d %8.4f %13.4e %7.3f %13.4e %7.3f\n",
			             order, level, d->sub->GetNE(), d->h, e.global, gRate,
			             e.region, rRate );
			std::fflush( stdout );

			// AFFINE, SO ONE STEP REACHES THE FLOOR. The coils enter as data
			// and the interior source is empty, so dF/dpsi = 0 and the
			// residual is affine in ( x, a ). Asserted on the residual history
			// rather than the count, which is M-219's lesson: whether a
			// second step is taken to notice is a property of the target.
			BOOST_TEST( meq::tests::oneExactStepReachedTheFloor(
			                s->solver->newtonResiduals() ),
				"k = " << order << ", level " << level << ": "
				<< meq::tests::stepReductionNotice( s->solver->newtonResiduals() ) );

			// THE INSTRUMENT, ONCE: the table's own L2 against the solver's.
			// It integrates the same rule with the same reference, so a
			// disagreement is the table -- a transposed point, a dropped
			// weight -- and every rate below would inherit it.
			if ( level == 0 && order == 1 )
			{
				mfem::FunctionCoefficient exact( [ &reference ]( mfem::Vector const &x )
				{
					return reference->psi( x( 0 ), x( 1 ) );
				} );
				double const direct = s->solver->potentialError( exact );
				// The table sums in double and potentialError() in real_t, so
				// in single they part at a few eps of the error. M-272.
				BOOST_TEST_REQUIRE( std::abs( direct - e.global )
				                    <= meq::tests::scalarAgreementFloor( 1.0e-10*direct,
				                                                         direct ),
					"the error table reads " << e.global << " where "
					"potentialError() reads " << direct );
			}

			// THE CONTROL, ONCE: the same solve with the coupling REMOVED, so
			// Gamma_h carries psi = 0. A coupling that did nothing would leave
			// the packs' field truncated at the boundary and every rate above
			// would still be a rate -- towards a different function.
			if ( level == 0 && order == 2 )
			{
				std::unique_ptr<VacuumSolve> const bare =
					solveMeshed( *d, order, defaultModes, machineCoils(), false );
				Errors const c = table.errors( bare->solver->potential() );
				controlRatio = c.region/e.region;
				std::printf( "       control, coupling removed: L2 D_h %.4e, "
				             "machine %.4e ( %.0fx this row )\n",
				             c.global, c.region, controlRatio );
				std::fflush( stdout );
			}
		}
	}

	std::printf( "\n    |psi_coil|: %.4e over D_h, %.4e over the machine region\n",
	             normGlobal, normRegion );

	for ( int order = 1; order <= orders; ++order )
	{
		std::size_t const k = static_cast<std::size_t>( order - 1 );
		double const wanted = order + 1.0;
		double const overall = meq::tests::rate( region[ k ].front(),
		                                         region[ k ].back(),
		                                         std::pow( 2.0, levels - 1 ) );
		double const finestGlobal = meq::tests::rate(
			global[ k ][ levels - 2 ], global[ k ][ levels - 1 ], 2.0 );
		std::printf( "    k = %d: machine region %.3f over the sequence, "
		             "D_h %.3f on the finest pair, k+1 = %.0f\n",
		             order, overall, finestGlobal, wanted );

		if ( !meq::tests::rateOfOrderIsMeasurable( wanted ) )
		{
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice(
				order, "psi against the packs' own field", wanted ) );
			continue;
		}
		for ( int level = 1; level < levels; ++level )
		{
			double const pair = meq::tests::rate( region[ k ][ level - 1 ],
			                                      region[ k ][ level ], 2.0 );
			BOOST_TEST( pair > wanted - 0.30,
				"k = " << order << ": psi converges on the coil field at " << pair
				<< " over the machine region between levels " << level - 1
				<< " and " << level << ", against k+1 = " << wanted << ". The "
				"exact answer is coilPsi() and the region is clear of the "
				"copper, so a short rate is the source term, the coupling or the "
				"transfer -- not a singularity." );
		}
		BOOST_TEST( overall > wanted - 0.15,
			"k = " << order << ": the machine-region rate over the whole ladder "
			"is " << overall << " against k+1 = " << wanted );
		BOOST_TEST( finestGlobal > wanted - 0.30,
			"k = " << order << ": over the whole of D_h the finest pair reads "
			<< finestGlobal << " against k+1 = " << wanted );
	}

	BOOST_TEST( controlRatio > 50.0,
		"removing the exterior coupling moved the machine-region error by only "
		<< controlRatio << "x, so the coupling is not what is carrying the "
		"packs' field across Gamma" );
}

/*
 * (b) THE MODE COUNT, SWEPT ON MEQ'S OWN SOLVE -- AND IT IS THE COPPER'S.
 *
 * MR-0 bounded N from the PLASMA's current: M-169 needed 32 modes at
 * rho_Gamma = 1.7 for 2.65e-05 at M = 2. With the packs inside Gamma the
 * plasma is not what sets it. The copper reaches rho = 2.176, so the packs'
 * trace on Gamma = 3.0 carries mode n at about ( 2.176/3.0 )^n = 0.725^n, and
 * that is the spectrum the coupling must hold. The trace tail column is that
 * spectrum, computed here from coilPsi() alone with no solve.
 *
 * WHAT IS ASSERTED. That the error over D_h falls by more than two decades and
 * then STOPS -- the last two rows agree to 2% -- because past that point the
 * truncation is below the discretisation error and N buys nothing; that
 * defaultModes is already on the plateau, so (a) and (c) are not measuring the
 * truncation; and that the machine region feels the truncation too, at N = 4,
 * but plateaus far earlier than D_h does, the modes past N being harmonic
 * fields that decay inward from Gamma as ( rho/rho_Gamma )^n.
 */
BOOST_AUTO_TEST_CASE( theModeCountIsSetByTheCopperAndNotByThePlasma )
{
	int const order = 3;
	int const level = 1;
	std::vector<int> const sweep = { 4, 8, 12, 16, 20, 24, 32, 48 };

	std::shared_ptr<meq::CoilSet> const reference = machineCoils();
	std::vector<double> const tail = traceTail( *reference, 96 );

	std::unique_ptr<Domain> const d = buildDomain( level );
	ReferenceTable const table( *d->sub, order, *reference );

	std::printf( "\n  MR-3 (b): the mode count, k = %d, level %d ( %d elements ), "
	             "Gamma = %.2f\n", order, level, d->sub->GetNE(), rhoGamma );
	std::printf( "    copper to rho = %.4f, so mode n falls as about %.3f^n on Gamma\n\n",
	             std::hypot( machine().outerRadius,
	                         machine().centreZ + 0.5*machine().length ),
	             std::hypot( machine().outerRadius,
	                         machine().centreZ + 0.5*machine().length )/rhoGamma );
	std::printf( "    %3s %13s %13s %13s %13s %13s\n", "N", "trace tail",
	             "L2 D_h", "relative", "L2 machine", "relative" );

	std::vector<double> global;
	std::vector<double> globalNorm;
	std::vector<double> region;
	for ( int modes : sweep )
	{
		std::unique_ptr<VacuumSolve> const s =
			solveMeshed( *d, order, modes, machineCoils() );
		Errors const e = table.errors( s->solver->potential() );
		global.push_back( e.global );
		globalNorm.push_back( e.globalNorm );
		region.push_back( e.region );
		std::printf( "    %3d %13.4e %13.4e %13.4e %13.4e %13.4e\n", modes,
		             tail[ static_cast<std::size_t>( modes ) ], e.global,
		             e.global/e.globalNorm, e.region, e.region/e.regionNorm );
		std::fflush( stdout );
	}

	std::size_t const last = sweep.size() - 1;
	std::size_t plateau = last;
	for ( std::size_t i = 0; i <= last; ++i )
		if ( global[ i ] <= 1.05*global[ last ] )
		{
			plateau = i;
			break;
		}
	std::printf( "\n    the error over D_h is within 5%% of N = %d's from N = %d "
	             "on; the trace tail there is %.3e\n\n", sweep[ last ],
	             sweep[ plateau ], tail[ static_cast<std::size_t>( sweep[ plateau ] ) ] );

	auto at = [ &sweep ]( std::vector<double> const &v, int modes )
	{
		std::size_t const i = static_cast<std::size_t>(
			std::find( sweep.begin(), sweep.end(), modes ) - sweep.begin() );
		BOOST_TEST_REQUIRE( i < v.size() );
		return v[ i ];
	};

	// IT FALLS: the modes are what carries the field across Gamma.
	BOOST_TEST( global.front() > 100.0*at( global, defaultModes ),
		"the error over D_h falls only from " << global.front() << " at N = "
		<< sweep.front() << " to " << at( global, defaultModes ) << " at N = "
		<< defaultModes << " -- the coupling is not resolving the packs' trace" );
	BOOST_TEST( region.front() > 100.0*at( region, defaultModes ),
		"the machine region does not feel the truncation at N = "
		<< sweep.front() << ": " << region.front() << " against "
		<< at( region, defaultModes ) );

	// AND IT STOPS: past the plateau N buys nothing, so the error left is the
	// discretisation's. A rising last row would be the coupling degrading with
	// N; a still-falling one would mean (a) is reading the truncation.
	//
	// ONLY WHERE THE PLATEAU IS THE DISCRETISATION'S. The transfer differences
	// its path endpoints at extensionFdStep(), so it carries eps/fd_step of the
	// field: 2.4e-05 relative in single, 6.6e-05 on this || psi ||, where the
	// single plateau sits at 6.8e-05 and moves by up to 20% from one N to the
	// next. That is round-off with nothing left to resolve, and the 2% band
	// cannot be asked of it. In double the floor is eleven orders down and
	// this always runs. M-275.
	bool const plateauIsResolved = global[ last ]
		> 4.0*globalNorm[ last ]
		  *static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() )
		  /meq::extensionFdStep();
	if ( !plateauIsResolved )
	{
		std::printf( "    *** SKIPPED \"it stops\": the plateau, %.3e, is under 4 eps/fd_step "
		             "x || psi || = %.3e, the extension's own floor in this format. M-275\n",
		             global[ last ], 4.0*globalNorm[ last ]
		             *static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() )
		             /meq::extensionFdStep() );
	}
	else
	{
		BOOST_TEST( std::abs( global[ last ]/at( global, defaultModes ) - 1.0 ) < 0.02,
			"between N = " << defaultModes << " and N = " << sweep[ last ]
			<< " the error over D_h moved from " << at( global, defaultModes )
			<< " to " << global[ last ] << ", so it has not stopped falling" );
	}
	BOOST_TEST( sweep[ plateau ] <= defaultModes,
		"the plateau starts at N = " << sweep[ plateau ] << ", beyond the "
		<< defaultModes << " every other case here runs at" );

	// AND THE MACHINE REGION PLATEAUS FIRST, which is the statement MR-4 can
	// spend: the plasma sees the truncation through ( rho/rho_Gamma )^n.
	BOOST_TEST( region[ 1 ] < 1.2*region[ last ],
		"the machine region is still " << region[ 1 ]/region[ last ]
		<< "x its plateau at N = " << sweep[ 1 ] );
}

/*
 * (c) AGAINST geq'S OWN VACUUM FIELD, WITH THE DISAGREEMENT ATTRIBUTED.
 *
 * tools/geq-benchmark/reference/vacuum/points.dat is geq's field on every
 * other node of its 65 x 65 grid. With no plasma that field is a closed-form
 * sum over geq's twenty point filaments -- export_vacuum.py asserts so to
 * 5e-16 -- so MEQ is being compared with a coil MODEL, and four numbers say
 * how much of the difference is whose:
 *
 *   1. geq against meq::filamentPsi() over meq::filamentStack( pack, 1, 10 ):
 *      the same twenty filaments in MEQ's own kernel. Round-off, which settles
 *      every convention at once -- amperes per filament, turns, the sign of
 *      each pack, weber per radian against weber.
 *   2. geq against coilPsi(): THE MODEL FLOOR. Ten filaments at the radial
 *      centre are the midpoint rule for the rectangle, so this is geq's
 *      quadrature of its own winding, measured with no solve anywhere -- the
 *      reference measured against itself, M-167's lesson.
 *   3. MEQ, packs meshed, against geq: equal to 2, because
 *   4. MEQ, packs meshed, against coilPsi() is hundreds of times smaller.
 *
 * AND THEN THE MODEL IS REMOVED. Handing MEQ geq's own twenty filaments,
 * SUBTRACTED ( CS-3; a filament has no other route ), the floor is gone and
 * what is left is the mode truncation of a trace with a logarithm 0.88 m from
 * Gamma. So the residual is attributed in both directions: MEQ reproduces geq
 * when given geq's machine, and differs from it by exactly geq's model when
 * given the winding.
 *
 * Points outside D_h are DROPPED, and so is any inside a conductor; geq's box
 * is wholly inside both Gamma and clear of the copper, so none is, and that is
 * asserted rather than assumed.
 */
BOOST_AUTO_TEST_CASE( theVacuumMirrorAgreesWithGeqToGeqsOwnCoilModel )
{
	int const order = 3;
	int const level = 1;

	std::shared_ptr<meq::CoilSet> const reference = machineCoils();
	std::vector<GeqPoint> const points = readGeqPoints();
	int const count = static_cast<int>( points.size() );

	meq::ConductorField const geqModel = geqFilaments();
	BOOST_TEST_REQUIRE( geqModel.containment( 0.0, rhoGamma ) > 0.0 );

	std::unique_ptr<Domain> const d = buildDomain( level );
	std::unique_ptr<VacuumSolve> const meshed =
		solveMeshed( *d, order, defaultModes, machineCoils() );
	std::unique_ptr<VacuumSolve> const subtracted =
		solveSubtracted( *d, order, defaultModes, geqModel );

	mfem::Vector positions( 2*count );
	for ( int i = 0; i < count; ++i )
	{
		positions( i ) = points[ static_cast<std::size_t>( i ) ].radius;
		positions( count + i ) = points[ static_cast<std::size_t>( i ) ].z;
	}
	meq::PointLocator locator( *d->sub );
	locator.locate( positions );
	mfem::Vector meshedValues;
	mfem::Vector remainder;
	locator.interpolate( meshed->solver->potential(), meshedValues );
	locator.interpolate( subtracted->solver->potential(), remainder );

	std::vector<double> self;
	std::vector<double> model;
	std::vector<double> meqGeq;
	std::vector<double> meqExact;
	std::vector<double> matched;
	std::vector<double> scale;
	int dropped = 0;
	for ( int i = 0; i < count; ++i )
	{
		GeqPoint const &p = points[ static_cast<std::size_t>( i ) ];
		if ( locator.codes()[ i ] == 2 || reference->indexContaining( p.radius, p.z ) >= 0 )
		{
			++dropped;
			continue;
		}
		double const filaments = geqModel.psi( p.radius, p.z );
		double const winding = reference->psi( p.radius, p.z );
		self.push_back( filaments - p.psi );
		model.push_back( winding - p.psi );
		meqGeq.push_back( meshedValues( i ) - p.psi );
		meqExact.push_back( meshedValues( i ) - winding );
		matched.push_back( remainder( i ) + filaments - p.psi );
		scale.push_back( p.psi );
	}
	BOOST_TEST_REQUIRE( dropped == 0,
		dropped << " of geq's " << count << " points are outside D_h or inside a "
		"pack, where its box should be inside Gamma and clear of the copper" );

	double const norm = rms( scale );
	double const conventions = rms( self )/norm;
	double const floor = rms( model )/norm;
	double const meqAgainstGeq = rms( meqGeq )/norm;
	double const meqAgainstExact = rms( meqExact )/norm;
	double const geqMachine = rms( matched )/norm;

	std::printf( "\n  MR-3 (c): against geq's vacuum field, %d points, k = %d, "
	             "level %d, N = %d -- rms relative to rms psi = %.4e\n\n",
	             count - dropped, order, level, defaultModes, norm );
	std::printf( "    1. geq's filaments in MEQ's kernel  vs geq   %11.3e\n",
	             conventions );
	std::printf( "    2. the winding, coilPsi()          vs geq   %11.3e  "
	             "<- geq's coil model\n", floor );
	std::printf( "    3. MEQ, packs meshed               vs geq   %11.3e\n",
	             meqAgainstGeq );
	std::printf( "    4. MEQ, packs meshed               vs 2.    %11.3e\n",
	             meqAgainstExact );
	std::printf( "    5. MEQ, geq's filaments subtracted vs geq   %11.3e\n\n",
	             geqMachine );
	std::fflush( stdout );

	// What this build can resolve, relative: the double build never reaches it.
	double const resolvable = meq::tests::roundOffFloor( norm )/norm;

	// 1. THE CONVENTIONS. Two implementations of one Green's function, scipy's
	// elliptic integrals against Carlson's, so round-off and not a tolerance.
	BOOST_TEST( conventions < 1.0e-12,
		"geq's twenty filaments evaluated by MEQ's own kernel differ from geq "
		"by " << conventions << " relative -- a convention is wrong: current per "
		"filament, turns, a pack's sign, or Wb against Wb/rad" );

	// 2 AND 3. THE DISAGREEMENT IS geq'S COIL MODEL.
	BOOST_TEST( floor > 1.0e-5,
		"geq's filaments and the winding agree to " << floor << ", so there is "
		"no model floor to attribute and this case has stopped measuring one" );
	BOOST_TEST( std::abs( meqAgainstGeq/floor - 1.0 ) < 0.05,
		"MEQ differs from geq by " << meqAgainstGeq << " where geq's own coil "
		"model differs from the winding by " << floor << " -- the disagreement "
		"is not the model, so some of it is MEQ's" );

	// 4. AND IT IS NOT MEQ'S.
	BOOST_TEST( meqAgainstExact < std::max( 0.02*floor, resolvable ),
		"MEQ's own error at geq's points is " << meqAgainstExact << ", not small "
		"against the model floor " << floor );

	// 5. AND IT GOES AWAY WHEN MEQ IS GIVEN geq'S MACHINE.
	BOOST_TEST( geqMachine < std::max( 1.0e-3*floor, resolvable ),
		"handed geq's own filaments, MEQ still differs from geq by "
		<< geqMachine << " against a model floor of " << floor );
}
