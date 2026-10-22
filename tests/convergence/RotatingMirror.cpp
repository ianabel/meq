/*
 * MR-4: THE FREE-BOUNDARY ROTATING MIRROR, THROUGH MEQ, AGAINST geq.
 * MIRROR-PLAN.md MR-4b, MR-4c and MR-4d.
 *
 * Everything at once: meq::NormalisedRotatingSource at two species, both labels
 * pinned at prescribed points on the midplane, the box support, geq's twenty
 * filaments subtracted, the exterior coupling live, and the bordered Newton
 * closing psi_ax, psi_bnd and the exterior modes together. VacuumMirror is the
 * same machine without the plasma, and MirrorMachine.hpp is what the two share.
 *
 * THE PROBLEM IS geq'S, READ FROM THE COMMITTED REFERENCE.
 * tools/geq-benchmark/reference/rotating/ -- export_rotating.py, MR-4a -- holds
 * what MEQ is GIVEN: meta.txt ( R_ref, the two label points, the box end, the
 * species ) and the tables, `Psi f df/dPsi` in MEQ's Psi and SI with n_s the
 * density on R = R_ref, read by GeqReference.hpp's one reader. Every profile is
 *
 *     f = f_edge + ( f_peak - f_edge ) C psi_n^3 ( 1 - psi_n )^4,
 *
 * so in MEQ's Psi = 1 - psi_n its slope vanishes as Psi^3 at the outer label
 * and as ( 1 - Psi )^2 at the inner: with g g' = 0 quasineutrality has no
 * psi-dependence at either label, phi_0' = 0 there, and F = mu0 R^2 dp/dpsi
 * vanishes at BOTH edges of the box -- at power j = 3 outside and j = 2 inside.
 * That is the box's precondition, and j >= 2 rather than 1 is what keeps the
 * edge from capping the order at k >= 2 -- CLAUDE_FB.md, *The plasma edge caps
 * the order*. The shape is ASYMMETRIC, so a table read in psi_n rather than Psi
 * gives different numbers: the abscissa's sense is something this case catches.
 * SplineProfile holds each table at its end value outside [ 0, 1 ].
 *
 * geq's ANSWER IS NOT TRACKED. solution.txt and psi_*.dat are export_rotating.py's
 * output and are gitignored; MR-4d reads them at run time and fails, naming the
 * command, when they are absent. No number of geq's is written in this file.
 *
 * THE SEQUENCE IS THE ONE THE SETTERS IMPOSE, MIRROR-PLAN.md MR-4b:
 *
 *   1. the vacuum on its own solver -- one solver holds one source -- by the
 *      subtracted route, so its state is the remainder the plasma solve's
 *      state is, and the labels' starting values are the TOTAL flux at the
 *      points;
 *   2. the plasma source, g g' = 0, box to geq's |z| = z_edge;
 *   3. both points, and setBoundaryFluxInitial( psi_vac( outer point ) ): its
 *      default 0 makes the first span psi_ax - 0 > 0 and puts the box around
 *      the axis -- measured, a cost of two Newton steps here rather than a
 *      failure, and asserted against below;
 *   4. setSource( source, psi_vac( inner point ) ), the datum, the extension,
 *      the DtN, and the vacuum GridFunction as the guess, which seeds q as well.
 *
 * THE FILAMENTS ARE SUBTRACTED AND NOT WRAPPED. MR-4b's step 2 names
 * meq::CoilAugmentedNormalisedSource, but that wrapper carries MESHED coils, and
 * a filament has no route but subtraction; with every conductor in psi_c the
 * wrapper would hold an empty set and add nothing. apps/meq.cpp takes the same
 * branch: under a subtracted conductor model it releases the meshed set and
 * hands the solver the bare plasma source. And it puts geq's own coil model in,
 * which M-260 measured MEQ reproducing to 3.9e-10 with no plasma -- so the coil
 * model leaves the comparison MR-4d makes.
 *
 * THE AXIAL CUTOFF IS MESHED TO. The box support's source JUMPS at |z| = z_edge
 * -- by 8% of its peak here, which checkAxialCutoff() reports -- and a jump
 * across a line that cuts elements caps the order like a source stepping at a
 * plasma edge, j = 0. The line is prescribed, so it is aligned to, as the
 * coils are: MirrorMachine.hpp's PlasmaMesh makes z = +-z_edge element faces
 * at every level. theAlignedCutoffKeepsTheOrder measures the difference
 * against the same mesh with the fragmentation dropped, which cuts it.
 *
 * WHERE psi IS READ, AND THE TRAP IN IT. The total flux at a point is the
 * remainder's polynomial AT the point plus psi_c evaluated EXACTLY there, which
 * is what the label rows do. totalPotential() carries the NODAL INTERPOLANT of
 * psi_c instead, off by its own interpolation error -- about 7e-06 here, two
 * hundred times geq's resolution floor.
 */

#define BOOST_TEST_MODULE RotatingMirror
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "mfem.hpp"

#include "meq/Config.hpp"
#include "meq/ConductorField.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"

#include "convergence/ConvergenceHarness.hpp"
#include "convergence/GeqReference.hpp"
#include "convergence/MirrorMachine.hpp"
#include "convergence/RoundOffFloor.hpp"

#include "SourceTree.hpp"

namespace
{
	using namespace meqtest::mirror;

	/// meshFile()'s stems, one per ctest since the two run in parallel; the
	/// mesh that cuts the cutoff adds "Cut".
	char const *const solveStem = "RotatingMirror";
	char const *const ladderStem = "RotatingMirrorLadder";
	/// The midplane suite's ladder, its own ctest and so its own files.
	char const *const midplaneStem = "RotatingMirrorMidplane";
	/// The stem sampled() meshes under: the resolution suite's, unless the
	/// midplane suite -- which runs in its own process -- has claimed it.
	char const *ladderStemInUse = ladderStem;

	/// The committed problem -- meta.txt, filaments.dat and the tables, in the
	/// source tree -- or MEQ_GEQ_ROTATING_REFERENCE for a regenerated one, a
	/// directory holding a whole export: the override GeqSourceComparison and
	/// MirrorMachine.hpp take.
	std::string rotatingDirectory()
	{
		char const *override = std::getenv( "MEQ_GEQ_ROTATING_REFERENCE" );
		return override ? std::string( override )
		                : meqtest::geq::inputDirectory( "rotating" );
	}

	/// geq's answer -- the psi files and solution.txt, untracked, in the runs
	/// directory -- or the same override as rotatingDirectory().
	std::string rotatingAnswerDirectory()
	{
		char const *override = std::getenv( "MEQ_GEQ_ROTATING_REFERENCE" );
		return override ? std::string( override )
		                : meqtest::geq::answerDirectory( "rotating" );
	}

	meqtest::geq::GeqCase readProblem()
	{
		meqtest::geq::GeqCase const p = meqtest::geq::readCase( rotatingDirectory() );
		BOOST_TEST_REQUIRE( p.species.size() == 2u,
			p.directory << "meta.txt names " << p.species.size() << " species; "
			"the rotating mirror is electrons and deuterons" );
		BOOST_TEST_REQUIRE( p.referenceRadius > 0.0 );
		BOOST_TEST_REQUIRE( p.mu0 == machine().mu0,
			"the rotating case's mu0 is not the vacuum case's: two machines" );
		BOOST_TEST_REQUIRE( std::isfinite( p.axisPointR ) );
		BOOST_TEST_REQUIRE( std::isfinite( p.boundaryPointR ) );
		BOOST_TEST_REQUIRE( std::isfinite( p.axialEdge ) );
		BOOST_TEST_REQUIRE( p.gridCommon > 1 );
		return p;
	}

	meqtest::geq::GeqCase const &problem()
	{
		static meqtest::geq::GeqCase const p = readProblem();
		return p;
	}

	/// What the source is built from: the species, omega, the gauge they are
	/// tabulated in, and R_ref -- zero under the midplane gauge, which has none.
	struct MirrorProfiles
	{
		std::vector<meq::Species> species;
		std::shared_ptr<meq::Profile const> omega;
		double referenceRadius = 0.0;
		meq::RotatingSource::Gauge gauge = meq::RotatingSource::Gauge::FixedRadius;
	};

	/// The reference's tables in MEQ's gauge: n_s on R = R_ref, never the raw
	/// Nraw_* that are geq's own.
	MirrorProfiles mirrorProfiles()
	{
		MirrorProfiles out;
		out.species = meqtest::geq::speciesFrom( problem(), false );
		out.omega = meqtest::geq::omegaFrom( problem() );
		out.referenceRadius = problem().referenceRadius;
		return out;
	}

	/// The least-squares slope of log v against log d: the power at which
	/// |F|/R vanishes at a label.
	double fittedPower( std::vector<double> const &d, std::vector<double> const &v )
	{
		double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
		double const n = static_cast<double>( d.size() );
		for ( std::size_t i = 0; i < d.size(); ++i )
		{
			double const x = std::log( d[ i ] );
			double const y = std::log( v[ i ] );
			sx += x;
			sy += y;
			sxx += x*x;
			sxy += x*y;
		}
		return ( n*sxy - sx*sy )/( n*sxx - sx*sx );
	}

	/// The same source with dF/dpsi scaled -- a Jacobian wrong by a known
	/// factor and a residual left exact -- for the control on what the Newton
	/// history can see.
	class SkewedJacobian : public meq::NormalisedRotatingSource
	{
		public:
			SkewedJacobian( double factorIn, MirrorProfiles const &p, double psiAxis,
			                double mu0 )
				: meq::NormalisedRotatingSource( p.species, p.omega,
				      std::make_shared<meq::ConstantProfile const>( 0.0 ),
				      p.referenceRadius, psiAxis, mu0 ),
				  factor( factorIn )
			{
			}

			double dFdPsi( double radius, double z, double psi ) const override
			{
				return factor*meq::NormalisedRotatingSource::dFdPsi( radius, z, psi );
			}

		private:
			double factor;
	};

	/// @a field at the points; NaN where one is not found.
	std::vector<double> valuesAt( mfem::GridFunction const &field,
	                              std::vector<double> const &radius,
	                              std::vector<double> const &z )
	{
		int const count = static_cast<int>( radius.size() );
		mfem::Vector positions( 2*count );
		for ( int i = 0; i < count; ++i )
		{
			positions( i ) = radius[ static_cast<std::size_t>( i ) ];
			positions( count + i ) = z[ static_cast<std::size_t>( i ) ];
		}
		meq::PointLocator locator( *field.FESpace()->GetMesh() );
		locator.locate( positions );
		mfem::Vector values;
		locator.interpolate( field, values );
		std::vector<double> out( static_cast<std::size_t>( count ),
		                         std::numeric_limits<double>::quiet_NaN() );
		for ( int i = 0; i < count; ++i )
			if ( locator.codes()[ i ] != 2 )
				out[ static_cast<std::size_t>( i ) ] = values( i );
		return out;
	}

	/// A whole solved state to start a run from in place of the vacuum:
	/// setInitialState()'s iterate and coefficients, the labels, and the
	/// crossings the bordered gauge's knots start at -- null for the start's
	/// own, found on its field.
	struct Start
	{
		mfem::Vector state;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		std::vector<double> coefficients;
		std::shared_ptr<meq::Profile const> crossing;
	};

	/// The filaments' own midplane crossing, defined with the vacuum-gauge
	/// suite.
	std::shared_ptr<meq::SplineProfile const> vacuumCrossing( unsigned int intervals );

	/// The midplane suite's three species, defined with it.
	MirrorProfiles threeSpecies( double peak, meq::RotatingSource::Gauge gauge, double floor );

	/// How one run is driven. The defaults are MR-4b's sequence as written.
	struct Options
	{
		int order = 2;
		int level = 0;
		int modes = defaultModes;
		bool aligned = true;           ///< the cutoff on element faces
		std::string stem = solveStem;  ///< meshFile()'s
		bool seedBoundary = true;      ///< setBoundaryFluxInitial() from the vacuum
		double relativeTolerance = -1.0;  ///< < 0 keeps the solver's default
		double absoluteTolerance = -1.0;
		int maxIterations = 40;
		double jacobianFactor = 1.0;   ///< != 1 builds SkewedJacobian
		/// A half-disc to cut D_h from in place of meshFile()'s -- the
		/// driver's own, for the driver suite. `aligned` and `stem` are then
		/// unread, since the file is whatever it is.
		std::string meshPath;
		/// No vacuum solve and no initial guess: the remainder starts at zero,
		/// so the labels start at the filaments' own flux at their points --
		/// the start a [[coils]] file with no [initialguess] gives the driver.
		bool coldStart = false;
		/// psi_ax's starting value under coldStart.
		double axisGuess = std::numeric_limits<double>::quiet_NaN();
		/// The deuterons' density built from the electrons' by quasineutrality,
		/// meq::neutralisingDensity(), as SourceFactory builds the species a
		/// file marks Neutralising -- rather than read from the reference's own
		/// table, which agrees with it to 4e-16.
		bool neutraliseDeuterium = false;
		/// Profiles to solve with in place of the reference's own -- another
		/// gauge's tables of the same plasma, or another plasma -- on the same
		/// machine, labels and box. Null takes mirrorProfiles().
		MirrorProfiles const *profiles = nullptr;
		/// A fixed-curve source's curve, installed before the solve so that
		/// the solver's own -- the vacuum's, label-checked -- is never asked
		/// for: what reproduces a plasma posed on one curve and gauged on
		/// another. Null leaves the solver to install the vacuum's.
		std::shared_ptr<meq::Profile const> gaugeCurve;
		/// The border columns differenced rather than analytic, the control
		/// that sees what the analytic columns leave out.
		bool differencedColumns = false;
		/// setMidplaneGaugeKnots(), the midplane gauge's knots; zero keeps
		/// the solver's own count.
		int gaugeKnots = 0;
		/// A whole state to start from, with no vacuum solve.
		Start const *start = nullptr;
		/// The profiles posed on a HELD table, @a heldCrossing or the
		/// filaments' own crossing -- Gauge::FixedCurve on it, so F is
		/// pointwise and there are no knots: what the midplane gauge's
		/// crossings would give if they never moved. Not a gauge for the
		/// midplane densities, and the control on the bordered one.
		bool holdCrossing = false;
		std::shared_ptr<meq::Profile const> heldCrossing;
		/// setMidplaneCrossingInitial() with no state: the bordered gauge's
		/// knots start on this table, and the iterate where it would anyway.
		std::shared_ptr<meq::Profile const> crossingStart;
		/// The plasma solve's HDG stabilisation: the constructor's tau, read
		/// as tau0 in tau = tau0/max( R, R_f ) when diffusionScaledTau. THE
		/// MIRROR'S, meq::MirrorStabilisation -- what the driver gives a box
		/// file -- unless a case says otherwise; the vacuum solve is
		/// mirrorSolver()'s, at the same numbers, either way.
		double tau = meq::MirrorStabilisation::tau0;
		bool diffusionScaledTau = true;
		double tauFloorRadius = meq::MirrorStabilisation::floorRadius;
	};

	/// Everything one plasma solve aliases, kept alive with it.
	struct MirrorRun
	{
		std::unique_ptr<Domain> domain;
		meq::ConductorField filaments{ geqFilaments() };
		std::unique_ptr<VacuumSolve> vacuum;
		std::unique_ptr<meq::NormalisedRotatingSource> source;
		std::unique_ptr<meq::ExteriorDtN> dtn;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<meq::GradShafranovSolver> solver;

		double vacuumAxis = 0.0;      ///< total vacuum psi at the inner point
		double vacuumBoundary = 0.0;  ///< and at the outer
		bool converged = false;
		std::string failure;
		double vacuumSeconds = 0.0;
		double plasmaSeconds = 0.0;
	};

	/// The total flux, psi_p + psi_c, at the points: @a remainder's polynomial
	/// there and psi_c exactly there. See the file comment for why not
	/// totalPotential().
	std::vector<double> totalAt( MirrorRun const &run,
	                             mfem::GridFunction const &remainder,
	                             std::vector<double> const &radius,
	                             std::vector<double> const &z )
	{
		std::vector<double> values = valuesAt( remainder, radius, z );
		for ( std::size_t i = 0; i < values.size(); ++i )
			values[ i ] += run.filaments.psi( radius[ i ], z[ i ] );
		return values;
	}

	/// totalAt() for the flux: q_p + q_c at the points, R component first,
	/// each the remainder's polynomial there plus ConductorField::flux()
	/// exactly there. MEQ's sign convention, which flux() and the filaments'
	/// flux() share.
	std::pair<std::vector<double>, std::vector<double>>
	totalFluxAt( MirrorRun const &run, mfem::GridFunction const &remainder,
	             std::vector<double> const &radius, std::vector<double> const &z )
	{
		mfem::FiniteElementSpace const &space = *remainder.FESpace();
		BOOST_TEST_REQUIRE( space.GetVDim() == 2 );
		BOOST_TEST_REQUIRE( space.GetOrdering() == mfem::Ordering::byNODES );
		mfem::FiniteElementSpace scalar( space.GetMesh(), space.FEColl() );
		std::vector<double> component[ 2 ];
		for ( int c = 0; c < 2; ++c )
		{
			mfem::GridFunction part( &scalar );
			for ( int i = 0; i < part.Size(); ++i )
				part( i ) = remainder( c*part.Size() + i );
			component[ c ] = valuesAt( part, radius, z );
		}
		for ( std::size_t i = 0; i < radius.size(); ++i )
		{
			double qR = 0.0, qZ = 0.0;
			run.filaments.flux( radius[ i ], z[ i ], qR, qZ );
			component[ 0 ][ i ] += qR;
			component[ 1 ][ i ] += qZ;
		}
		return { component[ 0 ], component[ 1 ] };
	}

	std::unique_ptr<MirrorRun> runMirror( Options const &o )
	{
		meqtest::geq::GeqCase const &p = problem();
		auto run = std::make_unique<MirrorRun>();
		if ( !o.meshPath.empty() )
			run->domain = buildDomainFromFile( o.meshPath, o.level );
		else
		{
			PlasmaMesh plasma;
			plasma.axialEdge = p.axialEdge;
			plasma.aligned = o.aligned;
			run->domain = buildDomain( o.level, o.aligned ? o.stem : o.stem + "Cut",
			                           plasma );
		}

		if ( o.start )
		{
			run->vacuumAxis = o.start->psiAxis;
			run->vacuumBoundary = o.start->psiBoundary;
		}
		else if ( o.coldStart )
		{
			run->vacuumAxis = run->filaments.psi( p.axisPointR, p.axisPointZ );
			run->vacuumBoundary = run->filaments.psi( p.boundaryPointR, p.boundaryPointZ );
		}
		else
		{
			auto const t0 = std::chrono::steady_clock::now();
			run->vacuum = solveSubtracted( *run->domain, o.order, o.modes, run->filaments );
			auto const t1 = std::chrono::steady_clock::now();
			run->vacuumSeconds = std::chrono::duration<double>( t1 - t0 ).count();

			std::vector<double> const labels = totalAt( *run,
				run->vacuum->solver->potential(), { p.axisPointR, p.boundaryPointR },
				{ p.axisPointZ, p.boundaryPointZ } );
			run->vacuumAxis = labels[ 0 ];
			run->vacuumBoundary = labels[ 1 ];
		}
		BOOST_TEST_REQUIRE( std::isfinite( run->vacuumAxis ) );
		BOOST_TEST_REQUIRE( std::isfinite( run->vacuumBoundary ) );
		double const axisStart = o.coldStart ? o.axisGuess : run->vacuumAxis;
		BOOST_TEST_REQUIRE( std::isfinite( axisStart ) );

		MirrorProfiles profiles = o.profiles ? *o.profiles : mirrorProfiles();
		if ( o.neutraliseDeuterium )
		{
			BOOST_TEST_REQUIRE( problem().species[ 1 ].name == "deuterium" );
			profiles.species[ 1 ].density = meq::neutralisingDensity( profiles.species, 1 );
		}
		if ( o.jacobianFactor != 1.0 )
			run->source = std::make_unique<SkewedJacobian>(
				o.jacobianFactor, profiles, axisStart, machine().mu0 );
		else
			run->source = std::make_unique<meq::NormalisedRotatingSource>(
				profiles.species, profiles.omega,
				std::make_shared<meq::ConstantProfile const>( 0.0 ),
				profiles.referenceRadius, axisStart, machine().mu0,
				meq::RotatingSource::Closure::Automatic,
				o.holdCrossing ? meq::RotatingSource::Gauge::FixedCurve : profiles.gauge );
		run->source->setPlasmaSupportShape( meq::PlasmaSupportShape::Box, p.axialEdge );
		if ( o.gaugeCurve )
			run->source->unnormalised().setMidplaneCrossing( o.gaugeCurve );
		else if ( o.holdCrossing )
			run->source->unnormalised().setMidplaneCrossing(
				o.heldCrossing ? o.heldCrossing : vacuumCrossing( 1024 ) );

		run->dtn = std::make_unique<meq::ExteriorDtN>( 0.0, rhoGamma, o.modes );
		run->solver = std::make_unique<meq::GradShafranovSolver>( *run->domain->sub,
		                                                          o.order, o.tau );
		meq::GradShafranovSolver &solver = *run->solver;
		if ( o.diffusionScaledTau )
			solver.setTauScaling( meq::GradShafranovSolver::TauScaling::Diffusion,
			                      o.tauFloorRadius );
		solver.setAxisFluxPoint( p.axisPointR, p.axisPointZ );
		solver.setBoundaryFluxPoint( p.boundaryPointR, p.boundaryPointZ );
		if ( o.seedBoundary )
			solver.setBoundaryFluxInitial( run->vacuumBoundary );
		solver.setSource( *run->source, axisStart );
		solver.setBoundaryData( run->zero );
		solver.setExtension( *run->domain->path, run->domain->gammaHMarker );
		solver.setConductorField( run->filaments );
		solver.setExteriorCoupling( *run->dtn );
		if ( run->vacuum )
			solver.setInitialGuess( run->vacuum->solver->potential() );
		if ( o.gaugeKnots > 0 )
			solver.setMidplaneGaugeKnots( o.gaugeKnots );
		if ( o.start )
		{
			solver.setInitialState( o.start->state, o.start->coefficients );
			solver.setMidplaneCrossingInitial( o.start->crossing );
		}
		else if ( o.crossingStart )
			solver.setMidplaneCrossingInitial( o.crossingStart );
		if ( o.relativeTolerance > 0.0 )
			solver.setNewtonControl( o.relativeTolerance, o.absoluteTolerance,
			                         o.maxIterations );
		if ( o.differencedColumns )
			solver.setBorderColumn( meq::GradShafranovSolver::BorderColumn::Differenced );

		auto const t2 = std::chrono::steady_clock::now();
		try
		{
			solver.solve();
			run->converged = true;
		}
		catch ( std::exception const &e )
		{
			run->failure = e.what();
		}
		auto const t3 = std::chrono::steady_clock::now();
		run->plasmaSeconds = std::chrono::duration<double>( t3 - t2 ).count();
		return run;
	}

	void printHistory( MirrorRun const &run )
	{
		meq::GradShafranovSolver const &solver = *run.solver;
		std::vector<double> const &history = solver.newtonResiduals();
		std::vector<meq::GradShafranovSolver::BorderStep> const &steps =
			solver.borderSteps();
		std::printf( "    %4s %14s %9s %11s %11s %11s %11s %8s %5s\n", "it", "merit",
		             "order", "|R|", "axis", "boundary", "exterior", "damping",
		             "threw" );
		for ( std::size_t i = 0; i < history.size(); ++i )
		{
			double order = std::numeric_limits<double>::quiet_NaN();
			if ( i >= 2 )
				order = meq::tests::newtonOrder( history[ i - 2 ], history[ i - 1 ],
				                                 history[ i ] );
			std::printf( "    %4zu %14.6e %9.3f", i, history[ i ], order );
			if ( i < steps.size() )
			{
				meq::GradShafranovSolver::BorderStep const &s = steps[ i ];
				std::printf( " %11.3e %11.3e %11.3e %11.3e %8.4f %5d", s.fieldNorm,
				             s.axis, s.boundary, s.exterior, s.damping, s.trialsThrew );
			}
			std::printf( "\n" );
		}
		std::fflush( stdout );
	}

	/// What one converged run measures.
	struct Measured
	{
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		double span = 0.0;
		double fieldAtAxisPoint = 0.0;
		double fieldAtBoundaryPoint = 0.0;
		double peak = 0.0;            ///< max |F|/R over the nodes
		/// |F|/R on the midplane just inside each label, at Psi = 1 - delta and
		/// Psi = delta, for each delta in deltas(); and AT the inner label.
		std::vector<double> innerSide;
		std::vector<double> outerSide;
		double atInnerLabel = 0.0;
		double innerPower = 0.0;      ///< fittedPower() over powerDeltas()
		double outerPower = 0.0;
		/// max |F|/R over nearDeltas(), both labels, against the peak.
		double nearLabels = 0.0;
		std::vector<double> amplitudes;  ///< the exterior's energy amplitudes
		/// The PLASMA's own: the energy amplitudes of a - a_vacuum. Under the
		/// split the coupling's coefficients are the TOTAL trace's, and the
		/// copper's is most of it.
		std::vector<double> plasmaAmplitudes;
		double truncation = 0.0;
		meq::GradShafranovSolver::AxialCutoffCheck cutoff;
		meq::GradShafranovSolver::AxisSourceCheck axis;
		meq::GradShafranovSolver::PlasmaSupportNodes support;
	};

	/**
	 * WHERE F IS READ AGAINST ITS EDGE POWER, AND WHY NOT AT THE LABEL ITSELF.
	 * The tables are cubic Hermite on 2049 knots, a spacing of 4.9e-4 in Psi,
	 * and inside the FIRST interval at the outer label the interpolant of
	 * f ~ Psi^4 is the cubic through ( 0, 0 ) and ( h, h^4 ) with zero slope at
	 * 0 -- whose slope is LINEAR in Psi, at a coefficient of order h^2. So F
	 * there falls as d^1 at about 1e-11 of its peak, and a fit through it
	 * reads 2.1 where the profile is 3. Measured: the outer label's |F|/R reads
	 * 6.5e-12 and 2.7e-11 of the peak at d = 1e-5 and 1e-4 -- the
	 * interpolant -- then 4.6e-08, 3.7e-07, 2.9e-06 at 2e-3, 4e-3, 8e-3, eight
	 * times a doubling, d^3 to 0.02. The inner label's f ~ ( 1 - Psi )^3 is
	 * reproduced by the cubic exactly and reads d^2 at every d.
	 *
	 * So the power is fitted from four to sixteen knot spacings in, where F is
	 * the profile's, and the knot-scale remainder is held to the plan's
	 * amplitude bound instead: below 1e-06 of the peak within 1e-4 of each
	 * label, which is MR-4a's own acceptance of geq's state.
	 */
	std::vector<double> const &powerDeltas()
	{
		static std::vector<double> const values = { 2.0e-3, 4.0e-3, 8.0e-3 };
		return values;
	}

	std::vector<double> const &nearDeltas()
	{
		static std::vector<double> const values = { 1.0e-5, 1.0e-4 };
		return values;
	}

	/// Every d printed: nearDeltas() then powerDeltas().
	std::vector<double> deltas()
	{
		std::vector<double> all = nearDeltas();
		all.insert( all.end(), powerDeltas().begin(), powerDeltas().end() );
		return all;
	}

	Measured measure( MirrorRun &run )
	{
		meqtest::geq::GeqCase const &p = problem();
		meq::GradShafranovSolver &solver = *run.solver;
		Measured m;
		m.psiAxis = solver.psiAxis();
		m.psiBoundary = solver.psiBoundary();
		m.span = m.psiAxis - m.psiBoundary;

		// The labels against the field AT THE POINTS, on the row's own terms.
		std::vector<double> const atPoints = totalAt( run, solver.potential(),
			{ p.axisPointR, p.boundaryPointR }, { p.axisPointZ, p.boundaryPointZ } );
		m.fieldAtAxisPoint = atPoints[ 0 ];
		m.fieldAtBoundaryPoint = atPoints[ 1 ];

		m.cutoff = solver.checkAxialCutoff();
		m.axis = solver.checkAxisSource();
		m.support = solver.plasmaSupportNodes();
		m.peak = m.cutoff.peak;

		// F at the source's own converged normalisation, which the solver set
		// before its last residual and is set again here so that nothing
		// depends on which trial that was.
		meq::NormalisedRotatingSource &source = *run.source;
		source.setNormalisation( m.psiAxis, m.psiBoundary );
		auto density = [ & ]( double radius, double normalised )
		{
			double const psi = m.psiBoundary + normalised*m.span;
			return std::abs( source.f( radius, 0.0, psi ) )/radius;
		};
		for ( double delta : deltas() )
		{
			m.innerSide.push_back( density( p.axisPointR, 1.0 - delta ) );
			m.outerSide.push_back( density( p.boundaryPointR, delta ) );
		}
		m.atInnerLabel = density( p.axisPointR, 1.0 );
		std::vector<double> inner, outer;
		for ( double delta : powerDeltas() )
		{
			inner.push_back( density( p.axisPointR, 1.0 - delta ) );
			outer.push_back( density( p.boundaryPointR, delta ) );
		}
		m.innerPower = fittedPower( powerDeltas(), inner );
		m.outerPower = fittedPower( powerDeltas(), outer );
		for ( double delta : nearDeltas() )
			m.nearLabels = std::max( { m.nearLabels,
				density( p.axisPointR, 1.0 - delta )/m.peak,
				density( p.boundaryPointR, delta )/m.peak } );

		m.amplitudes = run.dtn->modeAmplitudes( solver.exteriorCoefficients() );
		// The plasma's own needs the vacuum's on the same mesh and order, which
		// a cold start never solved for; it then stays empty.
		if ( run.vacuum )
		{
			std::vector<double> own = solver.exteriorCoefficients();
			std::vector<double> const &vacuum = run.vacuum->solver->exteriorCoefficients();
			for ( std::size_t i = 0; i < own.size(); ++i )
				own[ i ] -= vacuum[ i ];
			m.plasmaAmplitudes = run.dtn->modeAmplitudes( own );
		}
		m.truncation = run.dtn->truncationRatio( solver.exteriorCoefficients() );
		return m;
	}

	void printMeasured( Measured const &m, MirrorRun const &run )
	{
		meqtest::geq::GeqCase const &p = problem();
		std::printf( "    psi_ax %.12e, field at ( %.2f, 0 ) %.12e, difference "
		             "%.3e\n", m.psiAxis, p.axisPointR, m.fieldAtAxisPoint,
		             m.psiAxis - m.fieldAtAxisPoint );
		std::printf( "    psi_bnd %.12e, field at ( %.2f, 0 ) %.12e, difference "
		             "%.3e\n", m.psiBoundary, p.boundaryPointR, m.fieldAtBoundaryPoint,
		             m.psiBoundary - m.fieldAtBoundaryPoint );
		std::printf( "    span %.6e; the plasma moved psi_ax by %.6e and psi_bnd by "
		             "%.6e from the vacuum\n", m.span, m.psiAxis - run.vacuumAxis,
		             m.psiBoundary - run.vacuumBoundary );
		std::printf( "    peak |F|/R %.6e at ( %.4f, %.4f )\n", m.peak, m.cutoff.peakR,
		             m.cutoff.peakZ );
		std::printf( "    |F|/R / peak on the midplane:   delta" );
		for ( double delta : deltas() )
			std::printf( " %10.0e", delta );
		std::printf( "\n      inner, Psi = 1 - delta            " );
		for ( double v : m.innerSide )
			std::printf( " %10.3e", v/m.peak );
		std::printf( "   at Psi = 1: %.3e\n      outer, Psi = delta                " ,
		             m.atInnerLabel/m.peak );
		for ( double v : m.outerSide )
			std::printf( " %10.3e", v/m.peak );
		std::printf( "\n      fitted power over d = 2e-3 .. 8e-3: inner %.4f, outer "
		             "%.4f; within 1e-4 of either label at most %.3e of the peak\n",
		             m.innerPower, m.outerPower, m.nearLabels );
		std::printf( "    axial cutoff |z| = %.2f: %.6e at ( %.4f, %.4f ), %.3e of "
		             "the peak, %d samples\n", m.cutoff.axialEdge, m.cutoff.atCutoff,
		             m.cutoff.cutoffR, m.cutoff.cutoffZ, m.cutoff.relative,
		             m.cutoff.cutoffSamples );
		std::printf( "    axis: Psi on R = 0 %.6f, inside the plasma %d, source "
		             "vanishes there %d; support holds %d of %d nodes\n",
		             m.axis.normalisedFluxOnAxis, m.axis.axisInsidePlasma,
		             m.axis.sourceVanishesOnAxis, m.support.inside, m.support.total );
		auto tail = [ & ]( char const *label, std::vector<double> const &amplitudes )
		{
			double norm2 = 0.0;
			for ( double a : amplitudes )
				norm2 += a*a;
			std::printf( "    %s, N = %zu: norm %.4e; energy tail past N':", label,
			             amplitudes.size(), std::sqrt( norm2 ) );
			for ( std::size_t cut : { 4u, 8u, 16u, 24u, 28u, 30u } )
			{
				double rest = 0.0;
				for ( std::size_t i = cut; i < amplitudes.size(); ++i )
					rest += amplitudes[ i ]*amplitudes[ i ];
				std::printf( " %zu:%.2e", cut, std::sqrt( rest/norm2 ) );
			}
			std::printf( "\n" );
		};
		tail( "exterior, total trace ", m.amplitudes );
		if ( !m.plasmaAmplitudes.empty() )
			tail( "exterior, plasma's own", m.plasmaAmplitudes );
		std::printf( "    truncation ratio, total trace %.3e\n", m.truncation );
		std::fflush( stdout );
	}

	/**
	 * IS THE TAIL QUADRATIC? The LAST triple's order is above 1.5, or its last
	 * residual is at round-off, 2e-15 of the first -- the triple's first
	 * member never the guess, whose step measures how far the guess was.
	 *
	 * WHY NOT bestNewtonOrder(), WHICH IS WHAT THE PLAN NAMES. From the vacuum
	 * the first step takes the residual down five decades -- the plasma is a
	 * perturbation of a field the guess already holds -- and the third lands on
	 * round-off, so no triple lies above the floor in the asymptotic regime:
	 * bestNewtonOrder() reads 0.66 to 1.17 on these runs, which says nothing
	 * about the Jacobian.
	 *
	 * WHY ROUND-OFF COUNTS AS A PASS, which is PrescribedAxisPoint's argument
	 * about where a tail lands. The stopping target is 1e-12 of the cold
	 * residual, so a LINEAR tail stops just under it, near 1e-13 of the first
	 * residual; reaching 2e-15 from above the target takes a contraction
	 * better than 1e-3 in one step, which is the quadratic step.
	 * aWrongJacobianIsVisibleInTheHistory keeps the control standing: dF/dpsi
	 * scaled by 0.9, the residual exact, must fail this.
	 */
	bool quadraticTail( std::vector<double> const &history, std::string &why )
	{
		if ( history.size() < 4 )
		{
			why = "fewer than four residuals, so there is no triple past the "
			      "guess's own step";
			return false;
		}
		std::size_t const n = history.size() - 1;
		double const order = meq::tests::newtonOrder( history[ n - 2 ],
		                                              history[ n - 1 ], history[ n ] );
		// Ending within one eps of the first residual is ending at the
		// format's floor, where the last triple straddles it: in single a
		// four-residual history does, 3.1e-08 at three species and 5e19
		// (M-294). In double eps is under the 2e-15 this always
		// asked, so nothing moves there.
		double const relative = history[ n ]/history.front();
		double const floor = std::max( 2.0e-15, static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() ) );
		if ( order > 1.5 || relative <= floor )
			return true;
		char buffer[ 256 ];
		std::snprintf( buffer, sizeof buffer, "the last triple's order is %.3f and "
		               "it ends at %.3e of the first residual, above round-off; ",
		               order, relative );
		why = buffer;
		return false;
	}

	/// Whether this build leaves @a history a tail for quadraticTail() to
	/// read. In single a solve reaches the format's floor in three residuals,
	/// so there is no triple past the guess's own step: that is the format and
	/// not the Jacobian, and it is said loudly instead.
	/// aWrongJacobianIsVisibleInTheHistory runs in both builds, because a
	/// wrong Jacobian's history is long. M-272.
	bool tailIsReadable( std::vector<double> const &history )
	{
		if ( history.size() < 4 && meq::tests::newtonOrderIsUnmeasurable( 0.0 ) )
		{
			std::printf( "%s\n", meq::tests::newtonOrderNotice().c_str() );
			return false;
		}
		return true;
	}

	// ------------------------------------------------------------------
	// MR-4c and MR-4d: psi at geq's nodes
	// ------------------------------------------------------------------

	/// geq's common nodes, the 65 x 65 every resolution it was solved at
	/// shares, from meta.txt's grid line: z fastest, as the psi files list them.
	struct Nodes
	{
		std::vector<double> radius;
		std::vector<double> z;
	};

	Nodes const &geqNodes()
	{
		static Nodes const nodes = []
		{
			meqtest::geq::GeqCase const &p = problem();
			int const n = p.gridCommon;
			Nodes out;
			for ( int i = 0; i < n; ++i )
				for ( int j = 0; j < n; ++j )
				{
					out.radius.push_back( p.gridRMin
						+ ( p.gridRMax - p.gridRMin )*i/( n - 1.0 ) );
					out.z.push_back( p.gridZMin
						+ ( p.gridZMax - p.gridZMin )*j/( n - 1.0 ) );
				}
			for ( std::size_t i = 0; i < out.radius.size(); ++i )
				BOOST_TEST_REQUIRE( inMachineRegion( out.radius[ i ], out.z[ i ] ) );
			return out;
		}();
		return nodes;
	}

	/// One converged run, reduced to what MR-4c and MR-4d read from it.
	struct Sampled
	{
		int order = 0;
		int level = 0;
		bool aligned = true;
		int elements = 0;
		double h = 0.0;
		double seconds = 0.0;
		int newtonSteps = 0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		std::vector<double> psi;      ///< psi_h + psi_c at geqNodes(), NaN where not found
		std::vector<double> psiStar;  ///< psi* + psi_c
		std::vector<double> plasma;   ///< psi* alone: the remainder, the plasma's own flux
		std::vector<double> fluxR;    ///< q_h + q_c, the R component, at geqNodes()
		std::vector<double> fluxZ;    ///< and the z component
		/// psi* + psi_c at the two label points, which psiAxis() and
		/// psiBoundary() read by psi_h.
		double starAxis = 0.0;
		double starBoundary = 0.0;
		/// The merit history, and how often the midplane gauge tabulated
		/// its crossings off an iterate -- zero under the fixed radius.
		std::vector<double> history;
		int midplaneRefreshes = 0;
	};

	/// A converged @a run of @a o, reduced to a Sampled. The run is left
	/// post-processed and otherwise as it was.
	Sampled reduce( std::unique_ptr<MirrorRun> const &run, Options const &o )
	{
		if ( !run->converged )
			printHistory( *run );
		BOOST_TEST_REQUIRE( run->converged, "k = " << o.order << ", level " << o.level
			<< ( o.aligned ? ", aligned" : ", cut" ) << ": " << run->failure );

		Sampled s;
		s.order = o.order;
		s.level = o.level;
		s.aligned = o.aligned;
		s.elements = run->domain->sub->GetNE();
		s.h = run->domain->h;
		s.newtonSteps = static_cast<int>( run->solver->newtonResiduals().size() ) - 1;
		s.history = run->solver->newtonResiduals();
		s.midplaneRefreshes = run->solver->midplaneGaugeRefreshes();
		s.psiAxis = run->solver->psiAxis();
		s.psiBoundary = run->solver->psiBoundary();
		Nodes const &nodes = geqNodes();
		s.psi = totalAt( *run, run->solver->potential(), nodes.radius, nodes.z );
		std::tie( s.fluxR, s.fluxZ ) = totalFluxAt( *run, run->solver->flux(), nodes.radius,
		                                            nodes.z );
		run->solver->postProcess();
		mfem::GridFunction const &star = run->solver->postProcessedPotential();
		s.psiStar = totalAt( *run, star, nodes.radius, nodes.z );
		s.plasma = valuesAt( star, nodes.radius, nodes.z );
		meqtest::geq::GeqCase const &p = problem();
		std::vector<double> const labels = totalAt( *run, star,
			{ p.axisPointR, p.boundaryPointR }, { p.axisPointZ, p.boundaryPointZ } );
		s.starAxis = labels[ 0 ];
		s.starBoundary = labels[ 1 ];
		s.seconds = run->vacuumSeconds + run->plasmaSeconds;
		return s;
	}

	/// One converged run of @a o, reduced to a Sampled.
	Sampled sample( Options const &o )
	{
		return reduce( runMirror( o ), o );
	}

	/// sampled()'s runs, by ( level, k, aligned ): a case that has had to
	/// keep a run of its own -- the midplane transfer reads the best one's
	/// field -- files its reduction here rather than solving it twice.
	std::map<std::tuple<int, int, bool>, Sampled> &sampledCache()
	{
		static std::map<std::tuple<int, int, bool>, Sampled> cache;
		return cache;
	}

	/// Each ( level, k, aligned ) solved once per process, since MR-4d's best
	/// run is MR-4c's finest.
	Sampled const &sampled( int level, int order, bool aligned )
	{
		std::map<std::tuple<int, int, bool>, Sampled> &cache = sampledCache();
		auto const key = std::make_tuple( level, order, aligned );
		auto const found = cache.find( key );
		if ( found != cache.end() )
			return found->second;

		Options o;
		o.level = level;
		o.order = order;
		o.aligned = aligned;
		o.stem = ladderStemInUse;
		return cache.emplace( key, sample( o ) ).first->second;
	}

	/// max and rms of a - b over the nodes where both are finite.
	struct Difference
	{
		double max = 0.0;
		double rms = 0.0;
		int count = 0;
		std::size_t worst = 0;  ///< the node where max is
	};

	Difference difference( std::vector<double> const &a, std::vector<double> const &b )
	{
		BOOST_TEST_REQUIRE( a.size() == b.size() );
		Difference d;
		double sum = 0.0;
		for ( std::size_t i = 0; i < a.size(); ++i )
		{
			if ( !std::isfinite( a[ i ] ) || !std::isfinite( b[ i ] ) )
				continue;
			double const e = std::abs( a[ i ] - b[ i ] );
			if ( e > d.max )
			{
				d.max = e;
				d.worst = i;
			}
			sum += e*e;
			++d.count;
		}
		d.rms = d.count > 0 ? std::sqrt( sum/d.count ) : 0.0;
		return d;
	}

	/// difference() for the flux: max and rms of | q_a - q_b |, the
	/// Euclidean length of the difference, over the nodes where both are
	/// finite.
	Difference nodalFluxDifference( Sampled const &a, Sampled const &b )
	{
		BOOST_TEST_REQUIRE( a.fluxR.size() == b.fluxR.size() );
		std::vector<double> length( a.fluxR.size() );
		for ( std::size_t i = 0; i < length.size(); ++i )
			length[ i ] = std::hypot( a.fluxR[ i ] - b.fluxR[ i ], a.fluxZ[ i ] - b.fluxZ[ i ] );
		return difference( length, std::vector<double>( length.size(), 0.0 ) );
	}

	int finiteCount( std::vector<double> const &v )
	{
		return static_cast<int>( std::count_if( v.begin(), v.end(),
			[]( double x ) { return std::isfinite( x ); } ) );
	}

	/// One of geq's untracked psi files, R z psi on geqNodes() in order.
	std::vector<double> readGeqPsi( std::string const &name )
	{
		std::string const path = rotatingAnswerDirectory() + "/" + name;
		std::ifstream in = meqtest::geq::openOutput( rotatingAnswerDirectory(), name );
		Nodes const &nodes = geqNodes();
		std::vector<double> psi;
		std::string line;
		while ( std::getline( in, line ) )
		{
			if ( line.empty() || line[ 0 ] == '#' )
				continue;
			std::istringstream row( line );
			double radius = 0.0, z = 0.0, value = 0.0;
			BOOST_TEST_REQUIRE( static_cast<bool>( row >> radius >> z >> value ),
			                    "malformed row in " << path << ": " << line );
			std::size_t const i = psi.size();
			BOOST_TEST_REQUIRE( i < nodes.radius.size(),
				path << " has more rows than meta.txt's grid has nodes" );
			BOOST_TEST_REQUIRE( std::abs( radius - nodes.radius[ i ] ) <= 1.0e-12,
				path << " row " << i << " is at R = " << radius << " where meta.txt's "
				"grid puts node " << i << " at " << nodes.radius[ i ] );
			BOOST_TEST_REQUIRE( std::abs( z - nodes.z[ i ] ) <= 1.0e-12,
				path << " row " << i << " is at z = " << z << " where meta.txt's "
				"grid puts node " << i << " at " << nodes.z[ i ] );
			psi.push_back( value );
		}
		BOOST_TEST_REQUIRE( psi.size() == nodes.radius.size(),
			path << " has " << psi.size() << " rows for " << nodes.radius.size()
			<< " nodes" );
		return psi;
	}

	/// solution.txt's labels, `key psi_axis psi_bndry`.
	struct Labels
	{
		double axis = std::numeric_limits<double>::quiet_NaN();
		double boundary = std::numeric_limits<double>::quiet_NaN();
	};

	Labels readGeqLabels( std::string const &key )
	{
		std::string const path = rotatingAnswerDirectory() + "/solution.txt";
		std::ifstream in = meqtest::geq::openOutput( rotatingAnswerDirectory(), "solution.txt" );
		std::string line;
		while ( std::getline( in, line ) )
		{
			std::istringstream row( line );
			std::string name;
			row >> name;
			if ( name != key )
				continue;
			Labels out;
			BOOST_TEST_REQUIRE( static_cast<bool>( row >> out.axis >> out.boundary ),
			                    "malformed " << key << " in " << path );
			return out;
		}
		BOOST_TEST_REQUIRE( false, path << " has no " << key << "; regenerate it "
			"with " << meqtest::geq::regenerationCommand( rotatingAnswerDirectory() ) );
		return Labels();
	}
}

/*
 * SEVEN SUITES, ONE ctest EACH, so that the solve's acceptance does not wait on
 * the resolution study: `solve` is MR-4b and its controls, `resolution` is
 * MR-4c and MR-4d, which share their runs through sampled(), and `driver` is
 * MR-5 -- the shipped example through the shipped binary, against the library.
 * `adaptive` judges the estimator on that example and `device` runs it on the
 * card; both need the binary, as `driver` does. `midplane` and `vacuum_gauge`
 * are the two crossing gauges, MIRROR-PLAN.md section 5.
 */
BOOST_AUTO_TEST_SUITE( solve )

/*
 * THE FILAMENTS ARE THE REFERENCE'S OWN. geqFilaments() builds geq's twenty
 * from the VACUUM case's pack, and the rotating case writes the twenty it was
 * solved with; a filament moved or a current changed between the two exports
 * would be a disagreement MR-4d could only read as the plasma's.
 */
BOOST_AUTO_TEST_CASE( theSubtractedFilamentsAreTheRotatingCasesOwn )
{
	std::string const path = rotatingDirectory() + "/filaments.dat";
	std::ifstream in( path );
	BOOST_TEST_REQUIRE( in.good(), "cannot open " << path << ", which is tracked" );
	meq::ConductorField written( machine().mu0 );
	int count = 0;
	std::string line;
	while ( std::getline( in, line ) )
	{
		if ( line.empty() || line[ 0 ] == '#' )
			continue;
		std::istringstream row( line );
		double radius = 0.0, z = 0.0, current = 0.0;
		BOOST_TEST_REQUIRE( static_cast<bool>( row >> radius >> z >> current ),
		                    "malformed row in " << path << ": " << line );
		written.add( meq::CurrentFilament( radius, z, current ) );
		++count;
	}
	BOOST_TEST_REQUIRE( count == 2*machine().filamentsPerPack );

	meq::ConductorField const built = geqFilaments();
	Nodes const &nodes = geqNodes();
	double worst = 0.0;
	for ( std::size_t i = 0; i < nodes.radius.size(); ++i )
	{
		double const a = written.psi( nodes.radius[ i ], nodes.z[ i ] );
		double const b = built.psi( nodes.radius[ i ], nodes.z[ i ] );
		worst = std::max( worst, std::abs( a - b )/std::abs( a ) );
	}
	std::printf( "\n  filaments.dat against geqFilaments(), at geq's %zu nodes: worst "
	             "relative %.3e\n", nodes.radius.size(), worst );
	BOOST_TEST( worst < 1.0e-12,
		"the subtracted filaments differ from the rotating case's own by " << worst
		<< " relative, so the coil model is back in MR-4d's comparison" );
}

/*
 * MR-4b: THE SOLVE, AT k = 2 AND 3 ON THE FIRST TWO MESHES, AND WHAT IT IS
 * ASKED.
 *
 *   CONVERGED FROM THE VACUUM, with a negative span -- MIRROR-PLAN.md's first
 *   unknown.
 *
 *   QUADRATICALLY, by quadraticTail(); see it for why not the best triple.
 *
 *   THE SEED WAS THE VACUUM'S: both label rows read zero at iteration 0. With
 *   setBoundaryFluxInitial() left out the boundary row starts large and the
 *   first span is positive, putting the box around the axis; the bordered
 *   Newton recovers -- the boundary row is linear in the state, so the first
 *   step pins psi_bnd -- in more steps. So the plan's MANDATORY is a cost here
 *   and not a failure, and this is the assertion that sees it dropped.
 *
 *   F VANISHES AT EACH LABEL AT THE PROFILES' OWN POWER: |F|/R on the midplane
 *   at Psi = 1 - d and Psi = d falls as d^2 at the inner label and d^3 at the
 *   outer. A fixed threshold at one d would be a statement about the
 *   profile's amplitude; the power is the statement about the edge.
 *
 *   THE LABELS ARE THE FIELD AT THE POINTS, read on the row's own terms -- the
 *   remainder by its shape functions there and psi_c exactly there -- to
 *   round-off of the span.
 *
 *   AND THE DIAGNOSTICS: span < 0, the symmetry axis outside the box
 *   ( Psi( R = 0 ) > 1, MR-4's correction to checkAxisSource() ), the support
 *   non-empty, the axial cutoff and the exterior's trace tail printed. The
 *   tail is printed twice, because under the split the coupling's
 *   coefficients are the TOTAL trace's and the copper's is nearly all of it:
 *   the plasma's own is a - a_vacuum, on the same mesh and order.
 */
BOOST_AUTO_TEST_CASE( theRotatingMirrorConvergesFromTheVacuum )
{
	for ( int level : { 0, 1 } )
	{
		for ( int order : { 2, 3 } )
		{
			Options o;
			o.order = order;
			o.level = level;
			std::unique_ptr<MirrorRun> const run = runMirror( o );

			std::printf( "\n  MR-4b: k = %d, level %d ( %d elements ), N = %d -- vacuum "
			             "%.1f s, plasma %.1f s\n", order, level,
			             run->domain->sub->GetNE(), o.modes, run->vacuumSeconds,
			             run->plasmaSeconds );
			std::printf( "    vacuum labels: psi( %.2f, 0 ) = %.12e, psi( %.2f, 0 ) = "
			             "%.12e\n", problem().axisPointR, run->vacuumAxis,
			             problem().boundaryPointR, run->vacuumBoundary );
			printHistory( *run );
			if ( !run->converged )
				std::printf( "    NOT CONVERGED: %s\n", run->failure.c_str() );
			std::fflush( stdout );
			BOOST_TEST_REQUIRE( run->converged,
				"k = " << order << ", level " << level << ": the bordered Newton did "
				"not converge from the vacuum with a negative span: " << run->failure );

			Measured const m = measure( *run );
			printMeasured( m, *run );
			std::vector<double> const &history = run->solver->newtonResiduals();
			std::printf( "    best-triple order %.3f, the plan's measure, which a "
			             "history this short cannot carry\n",
			             meq::tests::bestNewtonOrder( history ) );
			std::fflush( stdout );

			// In single the solve reaches the format's floor in three residuals,
			// so there is no triple past the guess's own step to read a tail
			// from. That is the format and not the Jacobian, and it is skipped
			// loudly; aWrongJacobianIsVisibleInTheHistory still runs, because a
			// wrong Jacobian's history is long. M-272.
			if ( history.size() < 4 && meq::tests::newtonOrderIsUnmeasurable( 0.0 ) )
				std::printf( "%s\n", meq::tests::newtonOrderNotice().c_str() );
			else
			{
				std::string why;
				BOOST_TEST( quadraticTail( history, why ),
					"k = " << order << ", level " << level << ": " << why << "a tail "
					"that does not land where a quadratic model puts it is a Jacobian "
					"that does not match the residual" );
			}

			std::vector<meq::GradShafranovSolver::BorderStep> const &steps =
				run->solver->borderSteps();
			BOOST_TEST_REQUIRE( !steps.empty() );
			BOOST_TEST( std::abs( steps.front().boundary ) <= 1.0e-12*steps.front().fieldNorm,
				"the psi_bnd row reads " << steps.front().boundary << " at iteration 0. "
				"Seeded from the vacuum at the outer point it is zero; "
				"setBoundaryFluxInitial() was not given the vacuum's value, and the "
				"first span then puts the box around the axis" );
			BOOST_TEST( std::abs( steps.front().axis ) <= 1.0e-12*steps.front().fieldNorm,
				"the psi_ax row reads " << steps.front().axis << " at iteration 0, "
				"so setSource() was not given the vacuum's value at the inner point" );

			BOOST_TEST( m.span < 0.0,
				"a mirror's span is negative -- the inner label is the smaller "
				"flux -- and this one is " << m.span );
			// Two routes to one number: the border unknown, and the field read
			// at its point. Round-off of the format on the span. M-272.
			double const labelTolerance = meq::tests::scalarAgreementFloor(
				1.0e-12*std::abs( m.span ), m.span );
			BOOST_TEST( std::abs( m.psiAxis - m.fieldAtAxisPoint ) <= labelTolerance,
				"psi_ax " << m.psiAxis << " is not the field at its point, "
				<< m.fieldAtAxisPoint );
			BOOST_TEST( std::abs( m.psiBoundary - m.fieldAtBoundaryPoint ) <= labelTolerance,
				"psi_bnd " << m.psiBoundary << " is not the field at its point, "
				<< m.fieldAtBoundaryPoint );

			BOOST_TEST( std::abs( m.innerPower - 2.0 ) < 0.05,
				"|F|/R vanishes at the inner label as d^" << m.innerPower << " where "
				"the profiles' ( 1 - Psi )^2 slope puts it at d^2. Anything near d^0 "
				"is a source that STEPS at Psi = 1, which the box's precondition "
				"forbids" );
			BOOST_TEST( std::abs( m.outerPower - 3.0 ) < 0.05,
				"|F|/R vanishes at the outer label as d^" << m.outerPower
				<< " where the profiles' Psi^3 slope puts it at d^3" );
			BOOST_TEST( m.nearLabels < 1.0e-6,
				"|F|/R reaches " << m.nearLabels << " of its peak within 1e-4 of a "
				"label, where the plan's bound -- MR-4a's on geq's own state -- is "
				"1e-06" );

			BOOST_TEST( m.axis.reachesAxis );
			BOOST_TEST( m.axis.normalisedFluxOnAxis > 1.0 );
			BOOST_TEST( !m.axis.axisInsidePlasma,
				"checkAxisSource() puts the symmetry axis inside a box it is "
				"outside of, at Psi = " << m.axis.normalisedFluxOnAxis );
			BOOST_TEST( m.support.inside > 0 );
			BOOST_TEST( m.cutoff.applies );
			BOOST_TEST( m.cutoff.cutoffSamples > 0 );
		}
	}
}

/*
 * THE COLD-REFERENCE TARGET IS NOT LOOSE AGAINST THE PLASMA'S OWN FLUX.
 *
 * MIRROR-PLAN.md's third unknown: the stopping target is scaled from the cold
 * residual, which the 5.9 MA packs would dominate, while the plasma moves
 * psi_bnd by about 1e-03 of the label span's 6e-02. Under the split the coils
 * are not in the residual at all, and the quadratic tail lands at round-off
 * anyway -- so the same solve with the target taken three orders tighter must
 * move nothing the plasma is measured by.
 */
BOOST_AUTO_TEST_CASE( theStoppingTargetIsNotLooseAgainstThePlasmasOwnFlux )
{
#ifdef MFEM_USE_SINGLE
	// The default target is meq::newtonTolerance(), which in single IS the
	// format's floor, so there is no tighter target to compare it with: a
	// solve asked for 1e-15 cannot converge. A double-build acceptance. M-272.
	std::printf( "\n  *** SKIPPED: in single the default stopping target is the "
	             "format's own floor, so no tighter target exists to compare it "
	             "against. M-272\n" );
	return;
#endif
	Options base;
	std::unique_ptr<MirrorRun> const reference = runMirror( base );
	Options tight = base;
	tight.relativeTolerance = 1.0e-15;
	tight.absoluteTolerance = 1.0e-18;
	tight.maxIterations = 12;
	std::unique_ptr<MirrorRun> const tighter = runMirror( tight );
	BOOST_TEST_REQUIRE( reference->converged, reference->failure );
	BOOST_TEST_REQUIRE( tighter->converged, tighter->failure );

	mfem::GridFunction difference( tighter->solver->potential() );
	difference -= reference->solver->potential();
	double const plasmaFlux = std::abs( reference->solver->psiBoundary()
	                                    - reference->vacuumBoundary );
	double const moved = std::max( { static_cast<double>( difference.Normlinf() ),
		std::abs( tighter->solver->psiAxis() - reference->solver->psiAxis() ),
		std::abs( tighter->solver->psiBoundary() - reference->solver->psiBoundary() ) } );

	std::printf( "\n  THE STOPPING TARGET, k = %d, level %d: %zu steps at the default, "
	             "%zu at 1e-15\n", base.order, base.level,
	             reference->solver->newtonResiduals().size() - 1,
	             tighter->solver->newtonResiduals().size() - 1 );
	std::printf( "    the tighter solve moves psi by %.3e at most, against the "
	             "plasma's own %.3e at psi_bnd: %.3e of it\n", moved, plasmaFlux,
	             moved/plasmaFlux );
	std::fflush( stdout );

	BOOST_TEST( moved < 1.0e-8*plasmaFlux,
		"three orders on the stopping target moved psi by " << moved
		<< ", which is " << moved/plasmaFlux << " of the plasma's own flux" );
}

/*
 * THE CONTROL ON quadraticTail(): A JACOBIAN WRONG BY TEN PER CENT.
 *
 * dF/dpsi scaled by 0.9 and the residual left exact, so the solve converges
 * to the SAME equilibrium -- linearly, at about 1e-01 a step. If this ever
 * passed quadraticTail(), the main case's assertion could not see a wrong
 * Jacobian and its green would mean nothing.
 *
 * NOT ONE PER CENT, because that is a coin toss rather than a control. A
 * contraction of about 1e-03 a step stops anywhere from the 1e-12 target down
 * to about 1e-15 of the first residual, which straddles quadraticTail()'s
 * round-off branch at 2e-15 -- measured at 4.4e-15 on this problem, 2.2 times
 * clear. At 1e-01 a step the last residual cannot fall below 1e-13.
 */
BOOST_AUTO_TEST_CASE( aWrongJacobianIsVisibleInTheHistory )
{
	Options exact;
	std::unique_ptr<MirrorRun> const reference = runMirror( exact );
	Options skewed = exact;
	skewed.jacobianFactor = 0.9;
	std::unique_ptr<MirrorRun> const wrong = runMirror( skewed );

	std::printf( "\n  CONTROL: dF/dpsi x 0.9, k = %d, level %d\n", skewed.order,
	             skewed.level );
	printHistory( *wrong );
	BOOST_TEST_REQUIRE( reference->converged, reference->failure );
	BOOST_TEST_REQUIRE( wrong->converged,
		"the control did not converge, so its failing the tail test would be "
		"a failure to solve rather than a linear tail: " << wrong->failure );

	mfem::GridFunction difference( wrong->solver->potential() );
	difference -= reference->solver->potential();
	// Two solves stopped at the same target by different routes: in single
	// each is good to the format's floor, which a whole solve reaches at a few
	// thousand eps (exactAgreementFloor()'s argument). M-272.
	double const scale = std::abs( reference->solver->psiBoundary() );
	BOOST_TEST( difference.Normlinf()
	            < std::max( 1.0e-9*scale, meq::tests::exactAgreementFloor( scale ) ),
		"a Jacobian error moved the ANSWER by " << difference.Normlinf()
		<< "; it should move only the route" );

	std::string why;
	bool const passed = quadraticTail( wrong->solver->newtonResiduals(), why );
	std::printf( "    quadraticTail(): %s\n", passed ? "PASSED" : why.c_str() );
	std::fflush( stdout );
	BOOST_TEST( !passed,
		"a ten per cent Jacobian error passes quadraticTail(), so that test "
		"cannot see one" );
}


BOOST_AUTO_TEST_SUITE_END()

namespace
{
	/// MR-4c's reference and MR-4d's best run: the finest level at the highest
	/// order, read by psi*.
	int const bestLevel = 2;
	int const bestOrder = 4;

	Sampled const &reference()
	{
		return sampled( bestLevel, bestOrder, true );
	}

	/// The profiles' edge order at the INNER label, the smaller of the two,
	/// so the one that caps: p' ~ d^2 there and d^3 at the outer.
	double const edgeOrder = 2.0;

	/// M-09's caps for MEQ at the plain Gauss rule: psi_h at
	/// min( k+1, j+1.5 ), psi* at min( k+2, j+2 ).
	double potentialRate( int order )
	{
		return std::min( order + 1.0, edgeOrder + 1.5 );
	}

	double postProcessedRate( int order )
	{
		return std::min( order + 2.0, edgeOrder + 2.0 );
	}

	/**
	 * THE RUNGS EACH ORDER IS ASSERTED ON, AND WHY k <= 2 STARTS A LEVEL
	 * LATER. Level 0 puts 0.1 m elements across a plasma annulus 0.25 m wide
	 * at the midplane and 0.16 m at the throat, two or three elements, and
	 * at k <= 2 that is not yet the asymptotic regime: MEASURED, psi_h reads
	 * 1.29 then 1.81 at k = 1 and 2.38 then 2.83 at k = 2 over levels 0..2,
	 * and psi* at k = 1 RISES from level 0 to 1. So k = 1 and 2 are taken
	 * over levels 1..3 -- the precondition relaxed, the slack not -- and
	 * k = 3, which is asymptotic from level 0, over 0..2.
	 */
	std::vector<int> rungs( int order )
	{
		if ( order <= 2 )
			return { 1, 2, 3 };
		return { 0, 1, 2 };
	}

	/**
	 * The reference's own error, bounded from the one rung below it at its
	 * order: | psi*( level 1 ) - psi*( level 2 ) | / ( 2^p - 1 ), p the
	 * DESIGN rate psi* is capped at by the inner edge -- not the observed one,
	 * which is higher, so this over-states it.
	 */
	Difference const &referenceError()
	{
		static Difference const e = []
		{
			Difference d = difference( sampled( bestLevel - 1, bestOrder,
			                                    true ).psiStar,
			                           reference().psiStar );
			double const factor = std::pow( 2.0, postProcessedRate( bestOrder ) )
			                      - 1.0;
			d.max /= factor;
			d.rms /= factor;
			return d;
		}();
		return e;
	}

	/// geq's four grid corners, by index into geqNodes() -- z fastest.
	std::vector<std::size_t> gridCorners()
	{
		std::size_t const n = static_cast<std::size_t>( problem().gridCommon );
		return { 0, n - 1, ( n - 1 )*n, n*n - 1 };
	}

	/// @a v with the corners set to NaN, so difference() skips them.
	std::vector<double> withoutCorners( std::vector<double> v )
	{
		for ( std::size_t i : gridCorners() )
			v[ i ] = std::numeric_limits<double>::quiet_NaN();
		return v;
	}

	void printDifference( char const *label, Difference const &d )
	{
		std::printf( "    %-34s max %10.3e at ( %.3f, %+.3f )   rms %10.3e   %d nodes\n",
		             label, d.max, geqNodes().radius[ d.worst ], geqNodes().z[ d.worst ],
		             d.rms, d.count );
	}
}

BOOST_AUTO_TEST_SUITE( resolution )

/*
 * MR-4c: MEQ AGAINST ITSELF, AT geq'S NODES.
 *
 * No closed form exists, so the error at each rung is taken against
 * reference(): level 2 at k = 4, read by psi*. psi -- remainder at the point
 * plus psi_c exactly there -- on geq's 65 x 65 common nodes, every one of them
 * in the machine region and in D_h, which is the set MR-4d compares on.
 *
 * THE RATE THE METHOD PROMISES IS THE EDGE'S, NOT k+1. F vanishes as d^2 at
 * the inner label and d^3 at the outer, and M-09 caps MEQ at
 * min( k+1, j+1.5 ) for psi_h and min( k+2, j+2 ) for psi* at the shipped
 * rule, with j = 2 the smaller: psi_h 2, 3, 3.5 and psi* 3, 4, 4 at k = 1, 2,
 * 3. The axial cutoff is the third edge, and it is aligned to -- see
 * theAlignedCutoffKeepsTheOrder for what that is worth.
 *
 * ASSERTED ON THE RMS OVER THE NODES, by the two-tier rule CLAUDE_HDGGS.md
 * gives for the extension path: every pair within 0.30 of the design and the
 * whole ladder within 0.15. The maximum is printed with where it is, and not
 * asserted: a pointwise value of a discontinuous psi_h depends on where the
 * node sits in its element, and its pairs scatter.
 *
 * psi* IS ASSERTED AT EVERY ORDER AND psi_h ONLY WHERE THE LADDER REACHES
 * ITS ASYMPTOTE, which is k = 3. psi* is what the driver reports and what
 * MR-4d compares. psi_h's error is local -- psi_h minus its OWN psi* reads
 * the same digits as psi_h minus the reference -- and it sits in the plasma's
 * interior, 0.1 < Psi < 0.9, not at an edge, where it is still pre-asymptotic
 * at k <= 2 on every level this can afford. MEASURED: at k = 1, 1.81, 1.81
 * and then 1.94 between levels 3 and 4 ( 310,000 elements, two and a half
 * minutes ); at k = 2, 2.83 and 2.86. Those rows are printed against their
 * design and not asserted -- the precondition, the asymptotic regime, is
 * what is missing, and level 4 is what it costs.
 *
 * q_h IS PRINTED BESIDE THEM: | q_h - q_ref |, q = grad psi/R with the
 * filaments' own q_c added exactly, against the reference run's q_h -- what a
 * choice of the stabilisation trades psi_h against. ITS FINEST RUNGS READ THE
 * REFERENCE: the reference's q_h differs from level 1's at k = 4 by about
 * 3e-06 rms, and k = 2 at level 3 and k = 3 at level 2 both stop near
 * 1.1e-06 rms under constant tau = 1 and the mirror's tau = 3/max( R, 0.1 m )
 * alike. grad_bar( psi* )/R reads within 1.2e-07 rms of the reference's q_h,
 * so it is no better a reference; a finer one is ( MEASUREMENTS.md
 * M-292 ). So q_h is asserted where psi_h is -- k = 3, the asymptotic
 * rung -- and on the one pair the reference does not reach, levels 0 and 1,
 * at the edge's rate for psi_h with the pair's 0.30 of slack: 3.56 against
 * 3.5, MEASURED ( M-294 ). It does NOT see the axis: with the
 * plasma solve's floor at 1e-6 m, the unbounded tau0/R's stand-in, no error
 * column here moves by more than 4%, since geq's nodes stay at
 * R >= 0.05 -- VacuumMirror (a)'s q_h is the check that catches that. The
 * ladders at k <= 2 and the k = 3 rung at level 2 are printed and not
 * asserted.
 */
BOOST_AUTO_TEST_CASE( psiConvergesAtTheRateItsEdgesAllow )
{
	Sampled const &ref = reference();
	Nodes const &nodes = geqNodes();
	std::printf( "\n  MR-4c: against level %d, k = %d, psi* ( %d elements, %.1f s ): "
	             "%d of %zu nodes in D_h\n", ref.level, ref.order, ref.elements,
	             ref.seconds, finiteCount( ref.psiStar ), nodes.radius.size() );
	BOOST_TEST_REQUIRE( finiteCount( ref.psiStar ) == static_cast<int>( nodes.radius.size() ),
		"a node of geq's grid is not in D_h, which covers the machine region" );
	Difference const &own = referenceError();
	std::printf( "    the reference's own error, bounded from level %d at rate %.0f: "
	             "max %.3e, rms %.3e\n", bestLevel - 1,
	             postProcessedRate( bestOrder ), own.max, own.rms );
	Difference const ownFlux = nodalFluxDifference( sampled( bestLevel - 1, bestOrder, true ), ref );
	std::printf( "    the reference's q_h, level %d against %d at k = %d: max %.3e, rms %.3e\n",
	             bestLevel - 1, bestLevel, bestOrder, ownFlux.max, ownFlux.rms );

	std::printf( "\n    %2s %5s %7s %6s | %10s %10s %6s | %10s %10s %6s | %10s %10s %6s | %10s %10s\n",
	             "k", "level", "NE", "s", "psi_h max", "rms", "rate", "psi* max", "rms",
	             "rate", "q_h max", "rms", "rate", "psi*( ax )", "psi*( bnd )" );
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<double> potential, post, flux;
		std::vector<int> const levels = rungs( order );
		for ( int level : levels )
		{
			Sampled const &s = sampled( level, order, true );
			Difference const h = difference( s.psi, ref.psiStar );
			Difference const star = difference( s.psiStar, ref.psiStar );
			Difference const q = nodalFluxDifference( s, ref );
			potential.push_back( h.rms );
			post.push_back( star.rms );
			flux.push_back( q.rms );
			double hRate = std::numeric_limits<double>::quiet_NaN();
			double starRate = std::numeric_limits<double>::quiet_NaN();
			double qRate = std::numeric_limits<double>::quiet_NaN();
			std::size_t const n = potential.size();
			if ( n > 1 )
			{
				hRate = meq::tests::rate( potential[ n - 2 ], potential[ n - 1 ], 2.0 );
				starRate = meq::tests::rate( post[ n - 2 ], post[ n - 1 ], 2.0 );
				qRate = meq::tests::rate( flux[ n - 2 ], flux[ n - 1 ], 2.0 );
			}
			std::printf( "    %2d %5d %7d %6.1f | %10.3e %10.3e %6.3f | %10.3e %10.3e "
			             "%6.3f | %10.3e %10.3e %6.3f | %10.2e %10.2e\n", order, level, s.elements,
			             s.seconds, h.max, h.rms, hRate, star.max, star.rms, starRate, q.max, q.rms,
			             qRate, s.starAxis - ref.starAxis, s.starBoundary - ref.starBoundary );
			std::fflush( stdout );
		}
		std::printf( "       %-5s %.3f over the ladder, printed, not asserted\n", "q_h",
		             meq::tests::rate( flux.front(), flux.back(),
		                               std::pow( 2.0, levels.back() - levels.front() ) ) );
		if ( order >= 3 && meq::tests::rateOfOrderIsMeasurable( potentialRate( order ) ) )
		{
			double const firstPair = meq::tests::rate( flux[ 0 ], flux[ 1 ], 2.0 );
			std::printf( "       %-5s %.3f between levels %d and %d against %.1f, the pair the "
			             "reference does not reach\n", "q_h", firstPair, levels[ 0 ], levels[ 1 ],
			             potentialRate( order ) );
			BOOST_TEST( firstPair > potentialRate( order ) - 0.30,
				"k = " << order << ": q_h converges at " << firstPair << " between levels "
				<< levels[ 0 ] << " and " << levels[ 1 ] << " against the " << potentialRate( order )
				<< " the edges allow psi_h. A change to the stabilisation, the support or the "
				"coupling is the first suspect; the axis is not, since geq's nodes stay at "
				"R >= 0.05 and VacuumMirror (a) is the axis's check" );
		}

		double const span = std::pow( 2.0, levels.back() - levels.front() );
		struct Quantity
		{
			char const *name;
			std::vector<double> const *errors;
			double wanted;
			bool asserted;
		};
		for ( Quantity const &q : {
			Quantity{ "psi_h", &potential, potentialRate( order ), order >= 3 },
			Quantity{ "psi*", &post, postProcessedRate( order ), true } } )
		{
			std::vector<double> const &e = *q.errors;
			double const overall = meq::tests::rate( e.front(), e.back(), span );
			std::printf( "       %-5s %.3f over the ladder against %.1f%s\n", q.name,
			             overall, q.wanted, q.asserted ? "" : ", pre-asymptotic: "
			             "printed, not asserted" );
			if ( !q.asserted )
				continue;
			if ( !meq::tests::rateOfOrderIsMeasurable( q.wanted ) )
			{
				BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice(
					order, q.name, q.wanted ) );
				continue;
			}
			for ( std::size_t i = 1; i < e.size(); ++i )
			{
				double const pair = meq::tests::rate( e[ i - 1 ], e[ i ], 2.0 );
				BOOST_TEST( pair > q.wanted - 0.30,
					"k = " << order << ": " << q.name << " converges at " << pair
					<< " between levels " << levels[ i - 1 ] << " and " << levels[ i ]
					<< " against the " << q.wanted << " the edge orders allow. A "
					"short rate is the source, the support, the coupling or the "
					"transfer -- or the cutoff no longer on element faces" );
			}
			BOOST_TEST( overall > q.wanted - 0.15,
				"k = " << order << ": " << q.name << " converges at " << overall
				<< " over levels " << levels.front() << ".." << levels.back()
				<< " against " << q.wanted );
		}
	}
	std::fflush( stdout );
}

/*
 * THE AXIAL CUTOFF, MESHED TO AND NOT.
 *
 * The box support's source jumps at |z| = z_edge by several per cent of its
 * peak ( checkAxialCutoff(), MR-4b ), so where that line cuts elements the
 * exact psi carries |d|^2 across it inside elements -- a j = 0 edge. The same
 * ladder on the same mesh with the fragmentation dropped, against the same
 * reference, at k = 3, where the inner edge's cap on psi* is 4 and a j = 0
 * band would cap it near 2.
 *
 * MEASURED, AND IT IS THE POST-PROCESSED FIELD THAT SEES IT. At k = 3 psi*'s
 * rms reads 5.5 and 6.4 aligned against 4.2 and 3.6 cut, a factor 14 apart at
 * level 2; psi_h is 3.4 and 3.5 against 3.2 and 3.7, below the inner edge's
 * own cap either way. At k = 2 the two meshes read alike, psi* 4.4 and 4.3
 * against 4.1 and 4.6 -- at that order the edge caps first. The cutoff
 * crosses about 5% of the peak current, which is why the band costs a
 * constant here before it costs an order at these levels.
 */
BOOST_AUTO_TEST_CASE( theAlignedCutoffKeepsTheOrder )
{
	int const order = 3;
	Sampled const &ref = reference();
	std::printf( "\n  THE CUTOFF, k = %d: psi* against level %d, k = %d\n", order,
	             ref.level, ref.order );
	std::vector<double> aligned, cut;
	for ( int level : rungs( order ) )
	{
		Sampled const &a = sampled( level, order, true );
		Sampled const &c = sampled( level, order, false );
		Difference const da = difference( a.psiStar, ref.psiStar );
		Difference const dc = difference( c.psiStar, ref.psiStar );
		aligned.push_back( da.rms );
		cut.push_back( dc.rms );
		std::printf( "    level %d: aligned ( %5d elements ) max %.3e rms %.3e; cut "
		             "( %5d ) max %.3e rms %.3e\n", level, a.elements, da.max, da.rms,
		             c.elements, dc.max, dc.rms );
	}
	std::printf( "    rms rates: aligned %.3f %.3f, cut %.3f %.3f\n",
	             meq::tests::rate( aligned[ 0 ], aligned[ 1 ], 2.0 ),
	             meq::tests::rate( aligned[ 1 ], aligned[ 2 ], 2.0 ),
	             meq::tests::rate( cut[ 0 ], cut[ 1 ], 2.0 ),
	             meq::tests::rate( cut[ 1 ], cut[ 2 ], 2.0 ) );
	std::fflush( stdout );
	// In single both columns end on the format's floor, so the factor between
	// them is round-off: a double-build acceptance there, skipped loudly. The
	// floor rises with tau. At constant tau = 1 the aligned column ended at
	// 5.5e-08, under 16 eps x psi_bnd = 1.25e-07 ( M-272 ); under the mirror's
	// stabilisation every psi* rung from level 1 reads 2.4e-07 to 2.7e-07 rms,
	// flat, 35 eps x psi_bnd, so the guard is 128 eps x psi_bnd = 1.0e-06,
	// 3.6 above it, where double's aligned column is 5e-09
	// ( M-294 ).
	if ( !meq::tests::aboveTheFormatsFloor( aligned.back(), ref.psiBoundary, 128.0 ) )
	{
		std::printf( "    *** SKIPPED: the aligned column ends at %.3e, under this "
		             "build's floor on psi*, so the factor to the cut one is not "
		             "measurable. M-272\n", aligned.back() );
		return;
	}
	BOOST_TEST( aligned.back() < 0.2*cut.back(),
		"aligning the cutoff to element faces leaves psi* at " << aligned.back()
		<< " against " << cut.back() << " cut, so the alignment is not what "
		"MR-4c's k = 3 rate rests on -- or the line no longer cuts the other mesh" );
}

/*
 * MR-4d: MEQ AGAINST geq, AS A MATRIX.
 *
 * MEQ's best answer -- reference(), psi* -- against geq's direct free-boundary
 * solve at 129, 257 and 513 on geq's common nodes, and the labels against
 * solution.txt's. Only a difference that FALLS as geq refines rules out a
 * conversion error ( M-158's argument ): a wrong gauge, Psi's sense, eV for J
 * or a stray filament is a difference that stays. So:
 *
 *   it falls at every step, in both norms;
 *   at about geq's own rate -- each step's rate within 0.30 of geq's
 *   observed log2( |g129 - g257| / |g257 - g513| );
 *   and at 513 it is within geq's own 257-against-513 difference, which is
 *   what a geq converging at order 1.6 against an exact answer puts it at
 *   ( 1/( 2^1.6 - 1 ) = 0.49 of it ).
 *
 * MEQ's own error, MR-4c's bound on the reference, is printed beside it and
 * asserted below geq's 257-against-513, so the residual is geq's.
 *
 * THE CONTROL: shipped von Hagenow at 257, biased by its fixed eps = 1e-2 m.
 * The same MEQ answer must read that bias at its own size -- the comparison
 * resolves an error of about 1e-05 -- and far above the direct difference.
 *
 * geq'S FOUR GRID CORNERS ARE LEFT OUT, AND FOR A MEASURED REASON. There geq's
 * psi carries exactly HALF the plasma's own flux -- the remainder's, psi*
 * without psi_c -- at every resolution and under both boundary routes, where
 * the neighbouring nodes carry all of it to 3e-06: MEASURED, 0.500000 at
 * ( 0.05, +-1.5 ) and 0.500002 at ( 0.63, +-1.5 ), the latter 2.4e-05, more
 * than the Hagenow bias. A resolution-independent defect of the reference's
 * at four nodes, so it is asserted as found on the three direct solves -- if
 * geq's corners change, this says so and they rejoin the comparison -- and
 * the matrix is the other 4221.
 */
BOOST_AUTO_TEST_CASE( meqAgreesWithGeqAsGeqRefines )
{
	// geq's answer first, so that a clone without it fails in a second and
	// names the command, rather than after the solve.
	std::vector<double> const g129 = readGeqPsi( "psi_direct_129.dat" );
	std::vector<double> const g257 = readGeqPsi( "psi_direct_257.dat" );
	std::vector<double> const g513 = readGeqPsi( "psi_direct_513.dat" );
	std::vector<double> const hagenow = readGeqPsi( "psi_hagenow_257.dat" );
	for ( char const *key : { "labels_direct_129", "labels_direct_257",
	                          "labels_direct_513", "labels_hagenow_257" } )
		readGeqLabels( key );
	Sampled const &best = reference();

	std::printf( "\n  MR-4d: MEQ at level %d, k = %d, psi* ( %d elements ) against geq; "
	             "%d of %zu nodes in D_h\n", best.level, best.order, best.elements,
	             finiteCount( best.psiStar ), best.psiStar.size() );
	BOOST_TEST_REQUIRE( finiteCount( best.psiStar ) == static_cast<int>( best.psiStar.size() ) );

	/*
	 * IN SINGLE, geq'S OWN RESOLUTION IS UNDER MEQ'S FLOOR. geq's 257-against-
	 * 513 difference is 3.5e-08 at most, and MEQ's psi* in single is good to
	 * about 3e-07 there, so MEQ - geq stops falling at MEQ's floor and the
	 * corners' shares are MEQ's error over a small plasma flux. Every claim
	 * below that needs geq's resolution is then a double-build acceptance and
	 * is skipped loudly. What a conversion error would produce -- 1e-03 and up
	 * -- is still asserted, as is the Hagenow control's rms. M-272.
	 */
	double const fluxScale = std::abs( best.psiBoundary );
	bool const resolvable = meq::tests::aboveTheFormatsFloor(
		difference( withoutCorners( g257 ), g513 ).max, fluxScale );
	if ( !resolvable )
		std::printf( "    *** SKIPPED: geq's own 257 - 513 difference is under this "
		             "build's floor on psi*, so the matrix, its rates and the corners "
		             "are not measurable. The conversion bound and the Hagenow control "
		             "still run. M-272\n" );

	// The corners, asserted as found, on the direct solves: the Hagenow
	// control's carry its bias besides.
	for ( std::size_t i : resolvable ? gridCorners() : std::vector<std::size_t>() )
		for ( auto const *g : { &g129, &g257, &g513 } )
		{
			double const share = 1.0 - ( best.psiStar[ i ] - ( *g )[ i ] )/best.plasma[ i ];
			BOOST_TEST( std::abs( share - 0.5 ) < 1.0e-3,
				"at geq's corner ( " << geqNodes().radius[ i ] << ", " << geqNodes().z[ i ]
				<< " ) geq carries " << share << " of the plasma's own flux, where "
				"it carried one half. Its corners have changed: if they now carry all "
				"of it, put them back in the comparison" );
		}
	std::printf( "    geq's grid corners carry %.6f of the plasma's own flux "
	             "( %.3e at ( %.2f, %.2f ) ) and are left out\n",
	             1.0 - ( best.psiStar[ gridCorners().back() ] - g513[ gridCorners().back() ] )
	                   /best.plasma[ gridCorners().back() ],
	             best.plasma[ gridCorners().back() ], geqNodes().radius[ gridCorners().back() ],
	             geqNodes().z[ gridCorners().back() ] );

	std::vector<double> const meq = withoutCorners( best.psiStar );
	Difference const d129 = difference( meq, g129 );
	Difference const d257 = difference( meq, g257 );
	Difference const d513 = difference( meq, g513 );
	Difference const dControl = difference( meq, hagenow );
	Difference const geqCoarse = difference( withoutCorners( g129 ), g257 );
	Difference const geqFine = difference( withoutCorners( g257 ), g513 );
	Difference const geqBias = difference( withoutCorners( hagenow ), g257 );
	Difference const &own = referenceError();
	Difference const potential513 = difference( withoutCorners( best.psi ), g513 );

	std::printf( "\n" );
	printDifference( "MEQ - geq direct 129", d129 );
	printDifference( "MEQ - geq direct 257", d257 );
	printDifference( "MEQ - geq direct 513", d513 );
	printDifference( "geq 129 - geq 257", geqCoarse );
	printDifference( "geq 257 - geq 513", geqFine );
	printDifference( "MEQ's own error, bounded ( MR-4c )", own );
	printDifference( "MEQ psi_h - geq 513, for scale", potential513 );
	printDifference( "CONTROL: MEQ - geq Hagenow 257", dControl );
	printDifference( "geq Hagenow 257 - geq direct 257", geqBias );

	double const geqMax = std::log2( geqCoarse.max/geqFine.max );
	double const geqRms = std::log2( geqCoarse.rms/geqFine.rms );
	std::printf( "\n    rates, max: MEQ - geq %.3f then %.3f; geq's own %.3f\n",
	             std::log2( d129.max/d257.max ), std::log2( d257.max/d513.max ), geqMax );
	std::printf( "    rates, rms: MEQ - geq %.3f then %.3f; geq's own %.3f\n",
	             std::log2( d129.rms/d257.rms ), std::log2( d257.rms/d513.rms ), geqRms );
	std::printf( "    at 513: %.3f ( max ) and %.3f ( rms ) of geq's own 257 - 513\n",
	             d513.max/geqFine.max, d513.rms/geqFine.rms );

	std::printf( "\n    labels %22s %22s %12s %12s %12s %12s\n", "psi_axis", "psi_bndry",
	             "psi_ax - ", "psi_bnd - ", "psi*( ax ) -", "psi*( bnd ) -" );
	for ( char const *key : { "labels_direct_129", "labels_direct_257",
	                          "labels_direct_513", "labels_hagenow_257" } )
	{
		Labels const l = readGeqLabels( key );
		std::printf( "    %-19s %.15e %.15e %12.3e %12.3e %12.3e %12.3e\n", key, l.axis,
		             l.boundary, best.psiAxis - l.axis, best.psiBoundary - l.boundary,
		             best.starAxis - l.axis, best.starBoundary - l.boundary );
	}
	std::fflush( stdout );

	if ( !resolvable )
	{
		BOOST_TEST( d513.max < meq::tests::exactAgreementFloor( fluxScale ),
			"MEQ - geq 513 reads " << d513.max << ", which is a conversion and not "
			"a resolution: the labels, R_ref and the gauge, Psi's sense, J against eV" );
		BOOST_TEST( std::abs( dControl.rms - geqBias.rms ) < 0.05*geqBias.rms,
			"MEQ reads the Hagenow reference " << dControl.rms << " rms away where geq's "
			"own direct-to-Hagenow difference is " << geqBias.rms );
		return;
	}

	for ( auto const &norm : { std::make_tuple( "max", d129.max, d257.max, d513.max,
	                                            geqFine.max, geqMax, own.max ),
	                           std::make_tuple( "rms", d129.rms, d257.rms, d513.rms,
	                                            geqFine.rms, geqRms, own.rms ) } )
	{
		char const *name = std::get<0>( norm );
		double const a = std::get<1>( norm ), b = std::get<2>( norm ), c = std::get<3>( norm );
		double const floor = std::get<4>( norm ), geqRate = std::get<5>( norm );
		BOOST_TEST( ( a > b && b > c ),
			name << ": MEQ - geq reads " << a << ", " << b << ", " << c << " at 129, "
			"257, 513 and does not fall. A difference that stays is a conversion -- "
			"the labels, R_ref and the gauge, Psi's sense, J against eV, the cutoff, "
			"the filaments -- and not geq's resolution" );
		for ( auto const &step : { std::make_pair( a, b ), std::make_pair( b, c ) } )
		{
			double const rate = std::log2( step.first/step.second );
			BOOST_TEST( rate > geqRate - 0.30,
				name << ": MEQ - geq falls at " << rate << " where geq converges at "
				<< geqRate );
		}
		BOOST_TEST( c < floor,
			name << ": at 513 MEQ - geq is " << c << ", above geq's own 257 - 513 "
			"difference " << floor );
		BOOST_TEST( std::get<6>( norm ) < floor,
			name << ": MEQ's own error, bounded at " << std::get<6>( norm ) << ", is "
			"not below geq's 257 - 513 difference " << floor << ", so the residual "
			"cannot be put down to geq" );
	}

	// The labels: psi* at the points against geq's, by psiRZ at the same
	// points. geq's own converge to 1e-11, so these are within MEQ's error.
	Labels const l513 = readGeqLabels( "labels_direct_513" );
	BOOST_TEST( std::abs( best.starAxis - l513.axis ) < geqFine.max,
		"psi* at the inner point differs from geq's psi_axis by "
		<< best.starAxis - l513.axis );
	BOOST_TEST( std::abs( best.starBoundary - l513.boundary ) < geqFine.max,
		"psi* at the outer point differs from geq's psi_bndry by "
		<< best.starBoundary - l513.boundary );

	// The control: the bias read at its own size, far above the direct
	// difference.
	BOOST_TEST( std::abs( dControl.rms - geqBias.rms ) < 0.05*geqBias.rms,
		"MEQ reads the Hagenow reference " << dControl.rms << " rms away where geq's "
		"own direct-to-Hagenow difference is " << geqBias.rms );
	BOOST_TEST( dControl.max > 50.0*d257.max,
		"the shipped-Hagenow control reads " << dControl.max << " against the direct "
		"257's " << d257.max << ", so the comparison does not resolve a 1e-05 bias" );
}

BOOST_AUTO_TEST_SUITE_END()

#ifdef MEQ_RUN_PATH
/*
 * MR-5: THE SHIPPED MIRROR, THROUGH THE SHIPPED BINARY, AGAINST THE LIBRARY.
 *
 * examples/mirror-rotating.toml is what a user runs, and what rots about an
 * example is the FILE -- a key renamed, a table's path, a default moved -- so
 * this runs that file through meq-run, exactly as a user would, with only its
 * output redirected (MachineFixedBoundary's text substitution), and holds the
 * answer against the library solving geq's problem.
 *
 * THE LIBRARY SIDE IS BUILT FROM THE REFERENCE, NOT FROM THE TOML: the
 * species, tables, R_ref, labels, edge, filaments, Gamma and N are
 * runMirror()'s, read from tools/geq-benchmark/reference/. Only the mesh --
 * the driver's own, cut the driver's way -- and psi_ax's starting value are
 * taken from the file, the second being an initial value rather than the
 * problem. So agreement dof for dof says the file poses geq's problem, and
 * the start, mesh and iteration being the same, the bound is round-off.
 *
 * WHAT IT CANNOT SEE is the alignment of the cutoff, since both sides read
 * one mesh -- so the generator's command is asserted on its own.
 */
namespace
{
	std::string const mirrorExample = meq::test::sourcePath( "examples/mirror-rotating.toml" );
	/// The example GENERATES its mesh, and a generated [mesh] File is the run's
	/// product, so meq-run writes it relative to the working directory -- the
	/// build tree -- where the library side reads it back.
	char const *const mirrorMesh = "mirror-rotating.msh";
	/// Relative to the working directory, which is the build tree.
	char const *const driverScratch = "mirror-driver-scratch";

	std::string readText( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( std::istreambuf_iterator<char>( stream ),
		                    std::istreambuf_iterator<char>() );
	}

	int runLogged( std::string const &command, std::string const &log )
	{
		int const raw = std::system( ( command + " > " + log + " 2>&1" ).c_str() );
		return WIFEXITED( raw ) ? WEXITSTATUS( raw ) : -1;
	}

	/// The shipped example with [output] redirected into the scratch directory
	/// under @a prefix, and @a solverLines added under [solver]. TEXT
	/// SUBSTITUTION ON THE SHIPPED FILE, because a copy would drift and a
	/// drifted copy passing is worse than no test; shippedConfig() makes the
	/// tables it reads absolute, since the variant is written away from them.
	///
	/// @param lines   further lines replaced whole: the first line beginning
	///                with each key, which must exist, becomes its text.
	/// @param tables  text inserted before [output], for a table the example
	///                does not carry.
	std::string configureMirror( std::string const &prefix,
	                             std::string const &solverLines,
	                             std::vector<std::pair<std::string, std::string>> const &lines = {},
	                             std::string const &tables = "" )
	{
		std::string text = meq::test::shippedConfig( "examples/mirror-rotating.toml" );
		BOOST_TEST_REQUIRE( !text.empty(), mirrorExample << " is missing or empty" );

		auto replaceLine = [ &text ]( std::string const &key, std::string const &line )
		{
			std::size_t at = text.find( "\n" + key );
			BOOST_TEST_REQUIRE( at != std::string::npos,
			                    "the shipped example has no " << key << " line" );
			++at;
			std::size_t const end = text.find( '\n', at );
			text.replace( at, end - at, line );
		};
		replaceLine( "Directory = ", "Directory = \"" + std::string( driverScratch ) + "\"" );
		replaceLine( "Prefix = ", "Prefix = \"" + prefix + "\"" );
		for ( auto const &line : lines )
			replaceLine( line.first, line.second );
		if ( !solverLines.empty() )
		{
			std::size_t const at = text.find( "\n[solver]\n" );
			BOOST_TEST_REQUIRE( at != std::string::npos );
			text.insert( at + 10, solverLines );
		}
		if ( !tables.empty() )
		{
			std::size_t const at = text.find( "\n[output]\n" );
			BOOST_TEST_REQUIRE( at != std::string::npos );
			text.insert( at + 1, tables );
		}

		std::filesystem::create_directories( driverScratch );
		std::string const path = std::string( driverScratch ) + "/" + prefix + ".toml";
		std::ofstream out( path );
		out << text;
		return path;
	}

	struct DriverRun
	{
		int status = -1;
		std::string log;
		std::string stem;     ///< the output's path without its extension
		double seconds = 0.0;
	};

	/// meq-run on the shipped file, once per process: two cases read it.
	DriverRun const &shippedDriverRun()
	{
		static DriverRun const run = []
		{
			DriverRun out;
			std::string const prefix = "mirror-driver";
			std::string const config = configureMirror( prefix, "" );
			out.stem = std::string( driverScratch ) + "/" + prefix;
			auto const t0 = std::chrono::steady_clock::now();
			out.status = runLogged( "MKL_NUM_THREADS=1 " + std::string( MEQ_RUN_PATH )
			                        + " " + config + " --gf", out.stem + ".log" );
			out.seconds = std::chrono::duration<double>(
				std::chrono::steady_clock::now() - t0 ).count();
			out.log = readText( out.stem + ".log" );
			return out;
		}();
		return run;
	}

	/// The number after @a marker in @a text, or NaN.
	double numberAfter( std::string const &text, std::string const &marker )
	{
		std::size_t const at = text.find( marker );
		if ( at == std::string::npos )
			return std::numeric_limits<double>::quiet_NaN();
		return std::strtod( text.c_str() + at + marker.size(), nullptr );
	}

	bool says( std::string const &text, std::string const &phrase )
	{
		return text.find( phrase ) != std::string::npos;
	}

	/// meq-run on a variant of the shipped file, @a arguments after it, with
	/// @a environment before it.
	DriverRun runVariant( std::string const &prefix,
	                      std::vector<std::pair<std::string, std::string>> const &lines,
	                      std::string const &tables = "",
	                      std::string const &arguments = "",
	                      std::string const &environment = "" )
	{
		DriverRun out;
		std::string const config = configureMirror( prefix, "", lines, tables );
		out.stem = std::string( driverScratch ) + "/" + prefix;
		out.status = runLogged( environment + " MKL_NUM_THREADS=1 "
		                        + std::string( MEQ_RUN_PATH ) + " " + config + " --gf"
		                        + ( arguments.empty() ? "" : " " + arguments ),
		                        out.stem + ".log" );
		out.log = readText( out.stem + ".log" );
		return out;
	}

	/// A run's exact files read back: Omega_h, psi_h and the HDG flux q_h --
	/// `_grad_psi.gf` is the solved q = grad psi / R, not B (docs/output.rst).
	/// Under [conductors] both are the REMAINDER's, and so is what the
	/// estimator measures.
	struct StoredAnswer
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<mfem::GridFunction> psi;
		std::unique_ptr<mfem::GridFunction> flux;
		/// The physical flux psi_c + psi_h, where the run wrote one -- under
		/// [conductors] -- and null otherwise.
		std::unique_ptr<mfem::GridFunction> total;
	};

	StoredAnswer readAnswer( std::string const &stem )
	{
		StoredAnswer out;
		out.mesh = std::make_unique<mfem::Mesh>( ( stem + ".mesh" ).c_str(), 1, 1 );
		std::ifstream psi( stem + "_psi.gf" );
		BOOST_TEST_REQUIRE( psi.good(), "no " << stem << "_psi.gf" );
		out.psi = std::make_unique<mfem::GridFunction>( out.mesh.get(), psi );
		std::ifstream flux( stem + "_grad_psi.gf" );
		BOOST_TEST_REQUIRE( flux.good(), "no " << stem << "_grad_psi.gf" );
		out.flux = std::make_unique<mfem::GridFunction>( out.mesh.get(), flux );
		BOOST_TEST_REQUIRE( out.flux->FESpace()->GetVDim() == 2 );
		std::ifstream total( stem + "_psi_total.gf" );
		if ( total.good() )
			out.total = std::make_unique<mfem::GridFunction>( out.mesh.get(), total );
		return out;
	}

	/// || q_a - q_b || in L2 over Omega_h of @a a, at the rule of order 2k + 4
	/// GradShafranovSolver::fluxError() takes, and || q_b || on the same points.
	///
	/// THE TWO MESHES ARE DIFFERENT POLYGONS INSIDE ONE Gamma. Each is cut from
	/// its own background, so a sliver of Omega_h( a ) next to Gamma_h can lie
	/// outside Omega_h( b ); its points are dropped, and their measure returned
	/// so that a caller can say how much was.
	struct FluxDifference
	{
		double error = 0.0;
		double norm = 0.0;
		double area = 0.0;
		double uncovered = 0.0;
	};

	FluxDifference fluxDifference( StoredAnswer const &a, StoredAnswer const &b, int order )
	{
		mfem::Mesh &meshA = *a.mesh;
		std::vector<double> radius, z, weight, valueR, valueZ;
		mfem::Vector value( 2 );
		for ( int e = 0; e < meshA.GetNE(); ++e )
		{
			mfem::ElementTransformation *transformation = meshA.GetElementTransformation( e );
			mfem::IntegrationRule const &rule =
				mfem::IntRules.Get( meshA.GetElementBaseGeometry( e ), 2*order + 4 );
			for ( int i = 0; i < rule.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &point = rule.IntPoint( i );
				transformation->SetIntPoint( &point );
				mfem::Vector x;
				transformation->Transform( point, x );
				a.flux->GetVectorValue( *transformation, point, value );
				radius.push_back( x( 0 ) );
				z.push_back( x( 1 ) );
				weight.push_back( point.weight*transformation->Weight() );
				valueR.push_back( value( 0 ) );
				valueZ.push_back( value( 1 ) );
			}
		}

		// q_b one component at a time, since PointLocator interpolates a
		// scalar field: the component's coefficients copied into the scalar
		// space of the same collection.
		mfem::FiniteElementSpace const &vector = *b.flux->FESpace();
		mfem::FiniteElementSpace scalar( b.mesh.get(), vector.FEColl(), 1 );
		mfem::GridFunction componentR( &scalar ), componentZ( &scalar );
		for ( int i = 0; i < scalar.GetNDofs(); ++i )
		{
			componentR( i ) = ( *b.flux )( vector.DofToVDof( i, 0 ) );
			componentZ( i ) = ( *b.flux )( vector.DofToVDof( i, 1 ) );
		}

		int const count = static_cast<int>( radius.size() );
		mfem::Vector positions( 2*count );
		for ( int i = 0; i < count; ++i )
		{
			positions( i ) = radius[ static_cast<std::size_t>( i ) ];
			positions( count + i ) = z[ static_cast<std::size_t>( i ) ];
		}
		meq::PointLocator locator( *b.mesh );
		locator.locate( positions );
		mfem::Vector atR, atZ;
		locator.interpolate( componentR, atR );
		locator.interpolate( componentZ, atZ );

		FluxDifference out;
		for ( int i = 0; i < count; ++i )
		{
			std::size_t const p = static_cast<std::size_t>( i );
			out.area += weight[ p ];
			if ( locator.codes()[ i ] == 2 )
			{
				out.uncovered += weight[ p ];
				continue;
			}
			double const dR = valueR[ p ] - atR( i );
			double const dZ = valueZ[ p ] - atZ( i );
			out.error += weight[ p ]*( dR*dR + dZ*dZ );
			out.norm += weight[ p ]*( atR( i )*atR( i ) + atZ( i )*atZ( i ) );
		}
		out.error = std::sqrt( out.error );
		out.norm = std::sqrt( out.norm );
		return out;
	}

	/// One row of the driver's adaptive table.
	struct CycleRow
	{
		int elements = 0;
		int trace = 0;
		int marked = 0;
		double eta = 0.0;
	};

	/// The adaptive table the driver prints, `cycle elem trace marked wide eta
	/// it`, one row per cycle.
	std::vector<CycleRow> cycleTable( std::string const &log )
	{
		std::vector<CycleRow> rows;
		std::size_t const at = log.find( "cycle     elem" );
		if ( at == std::string::npos )
			return rows;
		std::istringstream stream( log.substr( at ) );
		std::string line;
		std::getline( stream, line );
		while ( std::getline( stream, line ) )
		{
			std::istringstream fields( line );
			int cycle = 0, wide = 0;
			CycleRow row;
			if ( !( fields >> cycle >> row.elements >> row.trace >> row.marked >> wide
			        >> row.eta ) )
				break;
			rows.push_back( row );
		}
		return rows;
	}
}

BOOST_AUTO_TEST_SUITE( driver )

/*
 * THE SHIPPED MIRROR SOLVES THROUGH THE DRIVER, AND IT IS geq'S PROBLEM.
 *
 *   THE MESH IS ALIGNED TO THE CUTOFF: `meq --mesh-command` carries
 *   --align-plasma and the box -1.5 .. 1.5, and no --coil, since the
 *   filaments are subtracted.
 *
 *   THE BINARY EXITS 0, and says what a mirror run says -- psi_ax a label at a
 *   prescribed point, psi_bnd seeded at the outer one, the axial cutoff -- and
 *   none of what a tokamak's axis check says of a mirror, nor a retry.
 *
 *   THE DRIVER'S psi_h IS THE LIBRARY'S, dof for dof, to round-off of its
 *   size, and in the same number of Newton steps -- the second is what sees
 *   the psi_bnd seed dropped, MR-4b's 3 against 5.
 *
 *   THE LABELS AND DIAGNOSTICS ARE MR-4b's: a negative span, each label the
 *   field at its point, the symmetry axis outside the box, a non-empty
 *   support and a cutoff that lies in the mesh. And the .nc says so.
 *
 *   AND psi_inner STARTS AT THE FILAMENTS' OWN FLUX AT THE ELECTRODE, read
 *   off the cold state there by the driver -- the file has no PsiAxis to give
 *   it -- so the inner label row reads zero at the first step.
 */
BOOST_AUTO_TEST_CASE( theShippedMirrorSolvesThroughTheDriver )
{
	std::filesystem::create_directories( driverScratch );
	meqtest::geq::GeqCase const &p = problem();
	meq::Configuration const file( mirrorExample );

	// ---- the generator's command -----------------------------------------
	{
		std::string const log = std::string( driverScratch ) + "/mesh-command.txt";
		int const status = runLogged( std::string( MEQ_DRIVER_PATH ) + " --mesh-command "
		                              + mirrorExample, log );
		std::string const command = readText( log );
		std::printf( "\n  MR-5: meq --mesh-command %s\n    %s", mirrorExample.c_str(),
		             command.c_str() );
		BOOST_TEST_REQUIRE( status == 0, command );
		BOOST_TEST( says( command, "--align-plasma" ),
			"the cutoff is not aligned to: the mesh would cut elements across "
			"|z| = AxialEdge, where the source jumps, and cost the order" );
		BOOST_TEST( says( command, "--plasma 0 -1.5 0.75 3 " ) );
		BOOST_TEST( says( command, "--rho 3.4 " ) );
		BOOST_TEST( !says( command, "--coil" ),
			"a subtracted conductor was emitted to the mesher" );
	}

	// ---- the file poses the reference's problem --------------------------
	meq::SourceConfig const &source = file.getSource();
	double const filamentsAtAxis = geqFilaments().psi( p.axisPointR, p.axisPointZ );
	BOOST_TEST( source.axisPoint.given );
	BOOST_TEST( source.axisPoint.radius == p.axisPointR );
	BOOST_TEST( source.axisPoint.z == p.axisPointZ );
	BOOST_TEST( file.getBoundary().limiter.radius == p.boundaryPointR );
	BOOST_TEST( file.getBoundary().limiter.z == p.boundaryPointZ );
	BOOST_TEST( source.axialEdge == p.axialEdge );
	BOOST_TEST( source.getRotating().referenceRadius == p.referenceRadius );

	// ---- the driver --------------------------------------------------------
	DriverRun const &run = shippedDriverRun();
	BOOST_TEST_REQUIRE( run.status == 0,
		"meq-run " << mirrorExample << " exited " << run.status << ":\n" << run.log );
	// THE INNER LABEL IS REPORTED ONCE, AS psi_inner AT THE ELECTRODE. A mirror
	// has no magnetic axis, so the tokamak's prescribed-point line -- psi_ax,
	// and the axis check that does not apply -- is not printed beside it.
	for ( char const *phrase : { "converged in ", "constraint psi_inner - psi_h(",
	                             "psi_inner starts at", "psi_bnd starts at", "axial cutoff",
	                             "the plasma support is the BOX" } )
		BOOST_TEST( says( run.log, phrase ), "the driver did not say '" << phrase << "'" );
	for ( char const *phrase : { "no interior extremum", "wall-hugging",
	                             "NOT a magnetic axis", "NO interior extremum",
	                             "did not converge", "Retrying",
	                             "support is EMPTY", "at the prescribed point",
	                             "     psi_ax = " } )
		BOOST_TEST( !says( run.log, phrase ),
			"the driver said '" << phrase << "' of a mirror:\n" << run.log );
	int const driverSteps = static_cast<int>( numberAfter( run.log, "converged in " ) );
	double const driverAxis = numberAfter( run.log, "     psi_inner = " );
	double const innerStart = numberAfter( run.log, "MEQ: psi_inner starts at " );
	BOOST_TEST( std::abs( innerStart - filamentsAtAxis ) <= 1.0e-6*std::abs( filamentsAtAxis ),
		"psi_inner started at " << innerStart << " where the filaments' flux at the "
		"electrode is " << filamentsAtAxis << ": the driver did not read the cold "
		"state at [boundary.electrode]" );

	// ---- the library, on geq's problem and the driver's mesh ---------------
	Options o;
	o.meshPath = mirrorMesh;
	o.level = 0;
	o.order = 2;
	o.modes = defaultModes;
	o.coldStart = true;
	o.neutraliseDeuterium = true;
	o.axisGuess = filamentsAtAxis;
	std::unique_ptr<MirrorRun> const library = runMirror( o );
	printHistory( *library );
	BOOST_TEST_REQUIRE( library->converged, library->failure );
	Measured const m = measure( *library );
	printMeasured( m, *library );
	int const librarySteps = static_cast<int>( library->solver->newtonResiduals().size() ) - 1;

	// ---- dof for dof -------------------------------------------------------
	mfem::Mesh storedMesh( ( run.stem + ".mesh" ).c_str(), 1, 1 );
	std::ifstream stream( run.stem + "_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good(), "the driver wrote no " << run.stem << "_psi.gf" );
	mfem::GridFunction const stored( &storedMesh, stream );
	mfem::GridFunction const &solved = library->solver->potential();
	BOOST_TEST_REQUIRE( storedMesh.GetNE() == library->domain->sub->GetNE() );
	BOOST_TEST_REQUIRE( stored.Size() == solved.Size() );
	double worst = 0.0;
	double largest = 0.0;
	for ( int i = 0; i < stored.Size(); ++i )
	{
		worst = std::max( worst, static_cast<double>( std::abs( stored( i ) - solved( i ) ) ) );
		largest = std::max( largest, static_cast<double>( std::abs( solved( i ) ) ) );
	}
	std::printf( "\n  MR-5: driver against library on %d elements: psi_h to %.3e over %d "
	             "dofs ( max |psi_h| %.3e, %.3e of it ); Newton %d and %d steps; driver "
	             "wall %.2f s, library solve %.2f s\n", storedMesh.GetNE(), worst,
	             stored.Size(), largest, worst/largest, driverSteps, librarySteps,
	             run.seconds, library->plasmaSeconds );
	std::printf( "    psi_inner printed %.6e, library %.12e\n", driverAxis, m.psiAxis );
	std::fflush( stdout );

	// IN SINGLE THE TWO ROUTES MEET AT THE EXTENSION'S OWN FLOOR. The transfer
	// differences its path endpoints at extensionFdStep(), so it carries
	// eps/fd_step of the field it moves -- 2.4e-05 relative in single -- and
	// the driver and the library reach their endpoints through two
	// expressions of one level set. Measured 1.578e-06, peaking at Gamma, against
	// eps/fd_step x |psi_bnd| = 1.556e-06: a ratio of 1.01. Applied in single
	// only, since in double that floor would loosen the bit-identity this
	// asserts there. M-272.
	double driverBound = meq::tests::scalarAgreementFloor( 1.0e-12, largest );
#ifdef MFEM_USE_SINGLE
	driverBound = std::max( driverBound, 4.0*std::abs( m.psiBoundary )
	                        *static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() )
	                        /meq::extensionFdStep() );
#endif
	BOOST_TEST( worst < driverBound,
		"the driver's psi_h differs from the library's on geq's problem by " << worst
		<< ", so the shipped file does not pose that problem" );
	BOOST_TEST( driverSteps == librarySteps,
		"the driver took " << driverSteps << " Newton steps and the library "
		<< librarySteps << "; a mirror's psi_bnd seeded anywhere but at the outer "
		"label costs steps, so look first at that" );
	BOOST_TEST( std::abs( driverAxis - m.psiAxis ) <= 1.0e-6*std::abs( m.psiAxis ),
		"the driver printed psi_inner = " << driverAxis << " against " << m.psiAxis );

	// ---- MR-4b's labels and diagnostics, on the library run ----------------
	BOOST_TEST( m.span < 0.0 );
	double const labelTolerance = meq::tests::scalarAgreementFloor(
		1.0e-12*std::abs( m.span ), m.span );
	BOOST_TEST( std::abs( m.psiAxis - m.fieldAtAxisPoint ) <= labelTolerance );
	BOOST_TEST( std::abs( m.psiBoundary - m.fieldAtBoundaryPoint ) <= labelTolerance );
	BOOST_TEST( !m.axis.axisInsidePlasma );
	BOOST_TEST( m.cutoff.applies );
	BOOST_TEST( m.cutoff.cutoffSamples > 0 );
	BOOST_TEST( m.support.inside > 0 );
	std::vector<meq::GradShafranovSolver::BorderStep> const &steps =
		library->solver->borderSteps();
	BOOST_TEST_REQUIRE( !steps.empty() );
	BOOST_TEST( std::abs( steps.front().axis ) <= 1.0e-12*steps.front().fieldNorm,
		"the psi_ax row reads " << steps.front().axis << " at iteration 0" );
	BOOST_TEST( std::abs( steps.front().boundary ) <= 1.0e-12*steps.front().fieldNorm,
		"the psi_bnd row reads " << steps.front().boundary << " at iteration 0" );

	// ---- and the interchange file says what the run was --------------------
	{
		std::string const header = std::string( driverScratch ) + "/header.txt";
		BOOST_TEST_REQUIRE( runLogged( "ncdump -h " + run.stem + ".nc", header ) == 0,
		                    readText( header ) );
		std::string const nc = readText( header );
		for ( char const *attribute : { ":axis_constraint = \"prescribed point\"",
		                                ":axis_point_r = 0.1 ;",
		                                ":plasma_support = \"box",
		                                ":plasma_support_axial_edge = 1.5 ;",
		                                ":axial_cutoff_relative = " } )
			BOOST_TEST( says( nc, attribute ), "the .nc lacks " << attribute );
		BOOST_TEST( !says( nc, ":axis_normalised_flux" ),
			"the .nc reports a magnetic axis a mirror does not have" );
	}
}

/*
 * A MIRROR WITH PicardSweeps IS REFUSED BEFORE ANYTHING IS SOLVED. Each sweep
 * solves the unbordered problem at frozen labels under Picard-then-Newton, and
 * on this file every one of them diverges -- a non-finite first residual --
 * while the bordered Newton needs no initialiser; so Config refuses the pair,
 * and this asserts the binary carries that through as exit 1 naming the key.
 */
BOOST_AUTO_TEST_CASE( aMirrorWithPicardSweepsIsRefusedAtParse )
{
	std::string const config = configureMirror( "mirror-driver-picard", "PicardSweeps = 2\n" );
	std::string const log = std::string( driverScratch ) + "/mirror-driver-picard.log";
	int const status = runLogged( std::string( MEQ_DRIVER_PATH ) + " " + config, log );
	std::string const text = readText( log );
	BOOST_TEST( status == 1, "exit " << status << ":\n" << text );
	BOOST_TEST( says( text, "solver.PicardSweeps" ), text );
}

/*
 * THE SHIPPED MIRROR IN THE MIDPLANE GAUGE, THROUGH THE DRIVER: its fixed-radius
 * tables read as midplane densities, which is another plasma and one the box's
 * precondition does not suit. ReferenceRadius out, Gauge = "midplane" in.
 *
 *   NEWTON IS QUADRATIC, in the fixed radius's steps or one more: the crossings
 *   are border unknowns. [M-284](MEASUREMENTS.md#m-284).
 *
 *   AND THE DRIVER WARNS ON THE ANSWER. geq's densities are flat at the outer
 *   label at fixed R_ref, so at fixed R_0 they carry -n m omega^2 R_0 R_0'/T
 *   there and F does not vanish at the label: the midplane gauge's
 *   precondition can only be read on the solved crossings, and the run says
 *   so -- exit 0, since a jump at a label costs the order, not the answer.
 */
BOOST_AUTO_TEST_CASE( theShippedMirrorInTheMidplaneGaugeIsQuadraticAndWarnsOnItsLabel )
{
	DriverRun const midplane = runVariant( "mirror-driver-midplane",
		{ { "ReferenceRadius = ", "Gauge = \"midplane\"" } } );
	BOOST_TEST_REQUIRE( midplane.status == 0, midplane.log );
	DriverRun const &fixed = shippedDriverRun();
	BOOST_TEST_REQUIRE( fixed.status == 0, fixed.log );
	double const steps = numberAfter( midplane.log, "converged in " );
	double const fixedSteps = numberAfter( fixed.log, "converged in " );
	std::printf( "\n  THE SHIPPED MIRROR IN THE MIDPLANE GAUGE: %g Newton steps, the fixed radius %g\n",
	             steps, fixedSteps );
	BOOST_TEST( steps <= fixedSteps + 1,
		"the midplane gauge took " << steps << " steps against the fixed radius's " << fixedSteps
		<< ", so its Newton is not quadratic:\n" << midplane.log );
	BOOST_TEST( says( midplane.log, "on the solved midplane crossings F at the outer label" ),
		"the driver did not warn that F stands at a label on the answer:\n" << midplane.log );
	BOOST_TEST( !says( fixed.log, "on the solved midplane crossings" ), fixed.log );
}

/*
 * A MIRROR WHOSE FIRST ATTEMPT FAILS IS NOT RETRIED: beside the box the
 * bordered Picard-then-Newton rung is skipped and the failure exits 2 as it
 * is. MEASUREMENTS.md M-277.
 *
 * WHY NOT: started from psi_inner = 0.2 rather than at the electrode's
 * flux, the rung converges in eight to a SECOND DISCRETE ROOT, well
 * converged, with one element of the box at -7.9e-02 where the shipped
 * answer's largest |psi_h| is 6.6e-03 ( M-277, at constant tau = 1 ). A
 * file cannot pose that start -- [source] PsiAxis is refused beside the box
 * and psi_inner starts at the state's flux at [boundary.electrode] -- but an
 * [initialguess] far from the answer can, and the rung would then report the
 * root it reaches. So the failure here is made the plainest way, a Newton
 * cap one below the steps the shipped file converges in, and exit 2 is the
 * assertion that matters: putting the rung back turns it into an exit 0.
 *
 * THE CAP IS READ FROM THE SHIPPED RUN, NOT WRITTEN DOWN. That count is three
 * in double and two in single, whose Newton stops at a looser residual, so a
 * fixed cap of 2 converged in single and the case failed for the wrong reason.
 *
 * Of the rung's two triggers this is the failure. The annulus cannot fire on
 * a box run, since [boundary.xpoint] is refused beside it, and
 * theShippedMirrorSolvesThroughTheDriver asserts that on a converged one.
 */
BOOST_AUTO_TEST_CASE( aMirrorWhoseFirstAttemptFailsIsNotRetried )
{
	DriverRun const &shipped = shippedDriverRun();
	BOOST_TEST_REQUIRE( shipped.status == 0, shipped.log );
	int const shippedSteps
		= static_cast<int>( numberAfter( shipped.log, "converged in " ) );
	BOOST_TEST_REQUIRE( shippedSteps >= 2,
		"the shipped mirror converges in " << shippedSteps << " Newton steps, so "
		"no smaller cap can make its first attempt fail:\n" << shipped.log );
	int const cap = shippedSteps - 1;
	std::printf( "\n  A MIRROR WHOSE FIRST ATTEMPT FAILS: the shipped file converges in "
	             "%d, capped at %d\n", shippedSteps, cap );

	std::string const prefix = "mirror-driver-no-retry";
	std::string const config = configureMirror( prefix, "",
	                                            { { "NewtonMaxIterations = ",
	                                                "NewtonMaxIterations = "
	                                                + std::to_string( cap ) } } );
	std::string const stem = std::string( driverScratch ) + "/" + prefix;
	std::remove( ( stem + ".nc" ).c_str() );
	int const status = runLogged( "MKL_NUM_THREADS=1 " + std::string( MEQ_RUN_PATH )
	                              + " " + config + " --gf", stem + ".log" );
	std::string const text = readText( stem + ".log" );
	BOOST_TEST( status == 2,
		"meq-run exited " << status << " on a mirror whose bordered Newton is "
		"capped at " << cap << " steps. An exit 0 is the retry rung run beside the box, which "
		"from a poor start reports a root that is not the shipped equilibrium:\n"
		<< text );
	BOOST_TEST( says( text, "the bordered Newton did not converge" ),
		"the first attempt's failure is not reported:\n" << text );
	BOOST_TEST( says( text, "retry is not run beside" ),
		"the driver did not say why the rung was skipped:\n" << text );
	for ( char const *phrase : { "Retrying", "bordered-picard-then-newton",
	                             "did not converge either", "NO interior extremum" } )
		BOOST_TEST( !says( text, phrase ),
			"the driver said '" << phrase << "' of a mirror:\n" << text );
	BOOST_TEST( !std::filesystem::exists( stem + ".nc" ),
		"a failed mirror wrote " << stem << ".nc" );
}
/*
 * THE TRUE MIDPLANE IN THE .nc, AND THAT IT POSES THE PLASMA AGAIN.
 * [M-282](MEASUREMENTS.md#m-282).
 *
 * examples/mirror-three-species.toml through meq-run three ways, on one
 * generated mesh:
 *
 *   1. AS SHIPPED, in the midplane gauge: the .nc's midplane_n_* -- the
 *      densities where each SOLVED surface crosses z = 0 -- are the DensityFile
 *      tables back EXACTLY AT THE GAUGE'S KNOTS, where the solver's crossings
 *      are the solved ones; between them the potential is pinned on the knots'
 *      interpolant, within a micron of the crossing, and the tables come back
 *      to what that gap moves them by -- the knots' own statement of the
 *      gauge, [M-284](MEASUREMENTS.md#m-284);
 *   2. THE SAME TABLES ON THE VACUUM'S CURVE, Gauge = "vacuum-midplane": a
 *      different plasma, whose .nc relabels it onto its own solved midplane;
 *   3. THOSE RELABELLED COLUMNS WRITTEN AS DensityFile TABLES and solved in the
 *      midplane gauge: the same plasma as 2, to what the knots' interpolant
 *      moves a density by -- which is what makes the columns a statement of
 *      the plasma rather than a picture of it.
 *
 * Read through ncdump, FluxGridConvergence's reason: the check does not share
 * code with the writer it checks.
 */
namespace
{
	char const *const threeSpeciesExample = "examples/mirror-three-species.toml";

	/// The values of one variable out of an `ncdump -v` dump.
	std::vector<double> dumpedValues( std::string const &text, std::string const &name )
	{
		std::size_t const data = text.find( "\ndata:" );
		BOOST_TEST_REQUIRE( data != std::string::npos, "ncdump printed no data section" );
		std::size_t const at = text.find( " " + name + " =", data );
		BOOST_TEST_REQUIRE( at != std::string::npos, "ncdump printed no values for " << name );
		std::size_t const begin = text.find( '=', at ) + 1;
		std::size_t const end = text.find( ';', begin );
		BOOST_TEST_REQUIRE( end != std::string::npos );
		std::vector<double> values;
		std::string const body = text.substr( begin, end - begin );
		std::size_t cursor = 0;
		while ( cursor < body.size() )
		{
			while ( cursor < body.size() && ( body[ cursor ] == ','
			        || std::isspace( static_cast<unsigned char>( body[ cursor ] ) ) ) )
				++cursor;
			if ( cursor >= body.size() )
				break;
			std::size_t used = 0;
			values.push_back( std::stod( body.substr( cursor ), &used ) );
			cursor += used;
		}
		return values;
	}

	/// The three-species example with [output] and its generated mesh in the
	/// scratch directory, each of @a lines' first line beginning with its key
	/// replaced whole, and each DensityFile whose path names @a files' first
	/// pointed at its second.
	std::string configureThreeSpecies( std::string const &prefix,
	                                   std::vector<std::pair<std::string, std::string>> lines,
	                                   std::vector<std::pair<std::string, std::string>> const &files = {} )
	{
		std::string text = meq::test::shippedConfig( threeSpeciesExample );
		BOOST_TEST_REQUIRE( !text.empty() );
		lines.emplace_back( "Directory = ", "Directory = \"" + std::string( driverScratch ) + "\"" );
		lines.emplace_back( "Prefix = ", "Prefix = \"" + prefix + "\"" );
		lines.emplace_back( "File = \"mirror-three-species.msh\"",
		                    "File = \"" + std::string( driverScratch ) + "/three-species.msh\"" );
		for ( auto const &line : lines )
		{
			std::size_t at = text.find( "\n" + line.first );
			BOOST_TEST_REQUIRE( at != std::string::npos, threeSpeciesExample << " has no " << line.first << " line" );
			++at;
			text.replace( at, text.find( '\n', at ) - at, line.second );
		}
		for ( auto const &file : files )
		{
			std::size_t const name = text.find( file.first );
			BOOST_TEST_REQUIRE( name != std::string::npos, threeSpeciesExample << " reads no " << file.first );
			std::size_t const at = text.rfind( '\n', name ) + 1;
			text.replace( at, text.find( '\n', at ) - at, "DensityFile = \"" + file.second + "\"" );
		}
		std::filesystem::create_directories( driverScratch );
		std::string const path = std::string( driverScratch ) + "/" + prefix + ".toml";
		std::ofstream out( path );
		out << text;
		return path;
	}

	/// One run of a three-species variant, and its .nc's midplane columns.
	struct MidplaneColumns
	{
		DriverRun run;
		int newtonIterations = -1;
		std::vector<double> psi, radius, radiusSlope;
		std::map<std::string, std::vector<double>> value, slope;
	};

	/// A .nc column and its slope column as a Hermite table on the .nc's knots.
	meq::SplineProfile columnTable( std::vector<double> const &psi, std::vector<double> const &value,
	                                std::vector<double> const &slope )
	{
		std::vector<meq::Knot> knots;
		for ( std::size_t i = 0; i < psi.size(); ++i )
			knots.push_back( meq::Knot{ psi[ i ], value[ i ], slope[ i ] } );
		return meq::SplineProfile( std::move( knots ) );
	}

	MidplaneColumns runThreeSpecies( std::string const &prefix,
	                                 std::vector<std::pair<std::string, std::string>> const &lines,
	                                 std::vector<std::pair<std::string, std::string>> const &files = {} )
	{
		MidplaneColumns out;
		std::string const config = configureThreeSpecies( prefix, lines, files );
		out.run.stem = std::string( driverScratch ) + "/" + prefix;
		out.run.status = runLogged( "MKL_NUM_THREADS=1 " + std::string( MEQ_RUN_PATH ) + " " + config + " --gf",
		                            out.run.stem + ".log" );
		out.run.log = readText( out.run.stem + ".log" );
		if ( out.run.status != 0 )
			return out;
		std::string variables = "midplane_Psi,midplane_R,midplane_dR_dPsi";
		for ( char const *s : { "electron", "deuterium", "carbon" } )
			variables += std::string( ",midplane_n_" ) + s + ",midplane_dn_" + s + "_dPsi";
		std::string const dump = out.run.stem + ".dump";
		BOOST_TEST_REQUIRE( runLogged( "ncdump -v " + variables + " " + out.run.stem + ".nc", dump ) == 0,
		                    readText( dump ) );
		std::string const text = readText( dump );
		out.newtonIterations = static_cast<int>( numberAfter( text, ":newton_iterations = " ) );
		out.psi = dumpedValues( text, "midplane_Psi" );
		out.radius = dumpedValues( text, "midplane_R" );
		out.radiusSlope = dumpedValues( text, "midplane_dR_dPsi" );
		for ( char const *s : { "electron", "deuterium", "carbon" } )
		{
			out.value[ s ] = dumpedValues( text, std::string( "midplane_n_" ) + s );
			out.slope[ s ] = dumpedValues( text, std::string( "midplane_dn_" ) + s + "_dPsi" );
			BOOST_TEST_REQUIRE( out.value[ s ].size() == out.psi.size() );
			BOOST_TEST_REQUIRE( out.slope[ s ].size() == out.psi.size() );
		}
		return out;
	}

	/// One species' columns as a DensityFile, `Psi f df/dPsi`.
	std::string writeDensityTable( MidplaneColumns const &c, std::string const &species, std::string const &name )
	{
		std::string const path = std::filesystem::absolute(
			std::filesystem::path( driverScratch ) / name ).string();
		std::ofstream out( path );
		out << "# midplane_n_" << species << " and its slope, relabelled from " << c.run.stem << ".nc\n";
		char line[ 128 ];
		for ( std::size_t i = 0; i < c.psi.size(); ++i )
		{
			std::snprintf( line, sizeof line, "%.17e %.17e %.17e\n", c.psi[ i ],
			               c.value.at( species )[ i ], c.slope.at( species )[ i ] );
			out << line;
		}
		return path;
	}
}

BOOST_AUTO_TEST_CASE( theNcCarriesTheTrueMidplaneAndItPosesThePlasmaAgain )
{
	// 1. As shipped: the midplane gauge's own tables, back.
	MidplaneColumns const shipped = runThreeSpecies( "three-midplane", {} );
	BOOST_TEST_REQUIRE( shipped.run.status == 0, shipped.run.log );
	std::printf( "\n  THE TRUE MIDPLANE IN THE .nc: %zu crossing knots, R_0 %.6f .. %.6f m; midplane gauge "
	             "%d Newton steps\n", shipped.psi.size(), shipped.radius.front(), shipped.radius.back(),
	             shipped.newtonIterations );
	BOOST_TEST( std::abs( shipped.radius.front() - 0.35 ) < 1.0e-4 );
	BOOST_TEST( std::abs( shipped.radius.back() - 0.1 ) < 1.0e-4 );
	BOOST_TEST( shipped.newtonIterations <= 4,
		"the midplane gauge took " << shipped.newtonIterations << " steps, so its Newton is not quadratic" );

	// THE GAP THE KNOTS LEAVE: the solved crossing, a Hermite table on the
	// .nc's own knots and slopes, against the interpolant through its values
	// at the solver's knots -- which is the table the potential was pinned on.
	meq::SplineProfile const crossing = columnTable( shipped.psi, shipped.radius, shipped.radiusSlope );
	auto const basis = std::make_shared<meq::CardinalSpline const>( 0.0, 1.0,
		static_cast<std::size_t>( meq::GradShafranovSolver::defaultMidplaneGaugeKnots ) );
	std::vector<double> atKnots;
	for ( std::size_t k = 0; k < basis->size(); ++k )
		atKnots.push_back( crossing( basis->knot( k ) ) );
	std::shared_ptr<meq::Profile const> const interpolant = meq::crossingThroughKnots( basis, atKnots );
	double knotGap = 0.0;
	for ( std::size_t i = 0; i < shipped.psi.size(); ++i )
		knotGap = std::max( knotGap, std::abs( ( *interpolant )( shipped.psi[ i ] ) - shipped.radius[ i ] ) );
	std::printf( "    the crossings against their %zu knots' interpolant: %.3e m at most\n", basis->size(), knotGap );

	BOOST_TEST( knotGap < 1.0e-5, "the crossings sit " << knotGap << " m from their knots' interpolant" );

	/*
	 * BETWEEN THE KNOTS THE TABLES COME BACK MOVED BY THE GAP, EXACTLY: a
	 * density read at the solved crossing R_0 is its table times exp( A_s ),
	 * A_s the exponent at R_0 of the source pinned on the interpolant -- which
	 * the test's closed-form statement of the same plasma, threeSpecies(),
	 * evaluates without the .nc. So the columns are the tables, transferred
	 * from the interpolant to the crossing, and nothing else.
	 */
	MirrorProfiles const closed = threeSpecies( 5.0e18, meq::RotatingSource::Gauge::Midplane, 1.0e-6 );
	meq::RotatingSource pinned( closed.species, closed.omega,
		std::make_shared<meq::ConstantProfile const>( 0.0 ), 0.0, machine().mu0,
		meq::RotatingSource::Closure::Automatic, meq::RotatingSource::Gauge::Midplane );
	pinned.setMidplaneCrossing( interpolant );
	std::map<std::string, double> worstOf;
	for ( std::size_t index : { std::size_t( 0 ), std::size_t( 2 ) } )
	{
		std::string const s = index == 0 ? "electron" : "carbon";
		meq::SplineProfile const table = meq::SplineProfile::fromFile(
			meq::test::sourcePath( "examples/mirror-three-species-n-" + s + ".dat" ) );
		double worstValue = 0.0, worstAgainstGap = 0.0;
		for ( std::size_t i = 0; i < shipped.psi.size(); ++i )
		{
			double const x = shipped.psi[ i ];
			double const stated = shipped.value.at( s )[ i ]/table( x ) - 1.0;
			double const fromGap = std::expm1( pinned.densityExponent( index, shipped.radius[ i ], x ) );
			worstValue = std::max( worstValue, std::abs( stated ) );
			worstAgainstGap = std::max( worstAgainstGap, std::abs( stated - fromGap ) );
		}
		worstOf[ s ] = worstValue;
		std::printf( "    %-9s midplane_n against its DensityFile: %.3e relative at most; against the table "
		             "moved by the gap, %.3e\n", s.c_str(), worstValue, worstAgainstGap );
		// A RATIO NEAR ONE IS KNOWN TO eps, so in single the two statements of
		// it meet at the format's floor: 4.2e-08 for the electrons and 1.9e-07
		// for carbon, against 16 eps = 1.9e-06. Skipped loudly there. M-272.
		if ( meq::tests::aboveTheFormatsFloor( worstAgainstGap, 1.0 ) )
			BOOST_TEST( worstAgainstGap < 1.0e-3*worstValue + 1.0e-9,
				s << ": the .nc's midplane density is not its table moved from the knots' interpolant to the "
				"solved crossing -- off by " << worstAgainstGap << " where the gap moves it by " << worstValue );
		else
			std::printf( "    *** SKIPPED: %s's two statements of the gap meet under this build's floor on a "
			             "ratio near one. M-272\n", s.c_str() );
	}

	// 2. The same tables on the vacuum's curve: another plasma, and quadratic.
	MidplaneColumns const vacuum = runThreeSpecies( "three-vacuum",
		{ { "Gauge = ", "Gauge = \"vacuum-midplane\"" } } );
	BOOST_TEST_REQUIRE( vacuum.run.status == 0, vacuum.run.log );
	std::printf( "    vacuum-midplane gauge: %d Newton steps\n", vacuum.newtonIterations );
	BOOST_TEST( vacuum.newtonIterations <= 4 );
	double moved = 0.0;
	for ( std::size_t i = 0; i < vacuum.psi.size(); ++i )
		moved = std::max( moved, std::abs( std::log( vacuum.value.at( "carbon" )[ i ]
			/meq::SplineProfile::fromFile( meq::test::sourcePath(
				"examples/mirror-three-species-n-carbon.dat" ) )( vacuum.psi[ i ] ) ) ) );
	std::printf( "    its true-midplane carbon against the table it was given: | ln ratio | up to %.3e\n", moved );
	BOOST_TEST( moved > 1.0e-3, "the vacuum curve's plasma relabels onto its own midplane as the "
		"table it was given, so the two gauges are not being told apart" );

	// 3. Its relabelled columns as tables, in the midplane gauge: the same plasma.
	std::string const electron = writeDensityTable( vacuum, "electron", "three-relabelled-n-electron.dat" );
	std::string const carbon = writeDensityTable( vacuum, "carbon", "three-relabelled-n-carbon.dat" );
	MidplaneColumns const again = runThreeSpecies( "three-relabelled", {},
		{ { "mirror-three-species-n-electron.dat", electron }, { "mirror-three-species-n-carbon.dat", carbon } } );
	BOOST_TEST_REQUIRE( again.run.status == 0, again.run.log );
	// psi*, the field the run reports, and psi_h beside it, dof for dof on the
	// one mesh; against the shipped run, which is another plasma.
	auto largestGap = []( std::string const &first, std::string const &second, double &scale )
	{
		StoredAnswer const a = readAnswer( first );
		StoredAnswer const b = readAnswer( second );
		std::ifstream starA( first + "_psistar.gf" ), starB( second + "_psistar.gf" );
		BOOST_TEST_REQUIRE( ( starA.good() && starB.good() ), "no _psistar.gf beside " << first );
		mfem::GridFunction const a2( a.mesh.get(), starA ), b2( b.mesh.get(), starB );
		BOOST_TEST_REQUIRE( a2.Size() == b2.Size(), "the two runs are not on one mesh" );
		double worst = 0.0;
		scale = 0.0;
		for ( int i = 0; i < a2.Size(); ++i )
		{
			worst = std::max( worst, double( std::abs( a2( i ) - b2( i ) ) ) );
			scale = std::max( scale, double( std::abs( a2( i ) ) ) );
		}
		return worst;
	};
	double scale = 0.0, shippedScale = 0.0;
	double const worst = largestGap( vacuum.run.stem, again.run.stem, scale );
	double const other = largestGap( vacuum.run.stem, shipped.run.stem, shippedScale );
	double worstColumn = 0.0;
	for ( std::size_t i = 0; i < again.psi.size(); ++i )
		worstColumn = std::max( worstColumn, std::abs( again.value.at( "carbon" )[ i ] - vacuum.value.at( "carbon" )[ i ] )
			/vacuum.value.at( "carbon" )[ i ] );
	std::printf( "    posed again in the midplane gauge, %d Newton steps: psi* differs by %.3e of its largest "
	             "%.3e, where the shipped tables' plasma differs by %.3e; the relabelled carbon by %.3e "
	             "relative\n", again.newtonIterations, worst/scale, scale, other/scale, worstColumn );
	// The relabelled tables are stated at the vacuum run's crossings, and the
	// midplane gauge pins on its knots' interpolant: the densities the two
	// plasmas carry differ by what part 1's gap moves a density by, and the
	// field by no more of its own scale than the electrons' share of that.
	// In single the two plasmas meet at the format's floor on psi*, 16 eps of
	// psi_bnd = 6.6e-02, which is 1.3e-07 against the 4.5e-08 they differ by:
	// the knots' share and the control's ratio are then round-off, and both
	// are skipped loudly. M-272.
	double const boundary = numberAfter( vacuum.run.log, "     psi_bnd = " );
	BOOST_TEST_REQUIRE( std::isfinite( boundary ), "the driver printed no psi_bnd" );
	if ( !meq::tests::aboveTheFormatsFloor( worst, boundary ) )
	{
		std::printf( "    *** SKIPPED: the relabelled plasma lands %.3e from the one it was read off, "
		             "under this build's floor on psi* of psi_bnd = %.3e. M-272\n", worst, boundary );
		return;
	}
	BOOST_TEST( worst <= worstOf.at( "electron" )*scale,
		"the relabelled densities posed in the midplane gauge solve to a different psi*, " << worst/scale
		<< " of its largest, beyond the " << worstOf.at( "electron" ) << " the knots' interpolant moves the "
		"electrons by -- so the .nc's midplane columns do not state the plasma that was solved" );
	BOOST_TEST( other > 1.0e4*worst,
		"the two gauges' plasmas on the same tables solve within " << other << " of each other, so this "
		"case cannot tell one plasma from another" );
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( adaptive )

/*
 * THE ESTIMATOR, JUDGED ON THE MIRROR: eta AGAINST THE TRUE ERROR, CYCLE BY
 * CYCLE, THROUGH THE DRIVER.
 *
 * DriverAcceptance's everyAdaptiveFrameSaysWhatEachElementIsMadeOf runs this
 * file adaptively for its frames and asks nothing of eta. Everything that
 * judges the estimator -- EstimatorConvergence, AdaptiveRefinement,
 * FreeBoundaryCoupling -- does it against a closed form, and none of them has
 * the split, the box or the labels. So here: the shipped file with three
 * Doerfler cycles at 0.5, DriverAcceptance's block, and each cycle's answer
 * against MR-4c's reference -- MEQ against itself, RefinementLevels = 2 at
 * k = 4 through the same binary -- in the quantity eta estimates, the flux:
 *
 *     effectivity = eta / || q - q_h ||,
 *
 * which AdaptiveRefinement holds in ( 5, 40 ) at k = 2 on graded meshes, 15.3
 * on uniform ones. The same band here, at every cycle and on the uniform L1
 * mesh too, whose eta is a one-cycle adaptive run's; it is not widened for a
 * mirror. And eta and the true error both fall at every cycle.
 *
 * EACH CYCLE'S ANSWER IS A RUN STOPPED THERE. A run writes the last cycle's
 * exact files and only VTK frames before it, so MaxIterations = 1, 2 and 3
 * are three runs; the loop is deterministic, and the element counts being the
 * long run's table is asserted, so the fields measured are its cycles.
 *
 * BOTH SIDES ARE THE REMAINDER: [conductors] subtracts the filaments, so the
 * solved q_h is psi_p's, and so is the residual eta is built on.
 *
 * WHAT IT FOUND, TWICE, AND BOTH WERE THE SPLIT. MEASUREMENTS.md M-279.
 *   eta_5's datum was the exterior's modal sum, the PHYSICAL flux on Gamma,
 *   where the load the solve imposes is that less psi_c: eta_5 read 2.68,
 *   3.58, 4.49 over the cycles, eta ROSE, and the marking chased it.
 *   apps/meq.cpp shifts the datum.
 *   And eta_1 read F at psi_p rather than at psi_p + psi_c, so its residual
 *   was O( 1 ) and eta fell at FIRST order: 1.48e-02, 7.46e-03, 3.69e-03
 *   under uniform refinement against a flux error falling at 3, an
 *   effectivity of 17, 68, 308. ResidualEstimator reads the total flux.
 * With the first repaired alone the cycles read 16.9, 32.7, 41.2 -- out of
 * the band only at the last -- which is why the uniform L1 index is asserted
 * beside them: it reads 68 there.
 *
 * THE REFERENCE RESOLVES WHAT IS MEASURED AGAINST IT. Its own flux error is
 * bounded from L1 at k = 4 as referenceError() bounds psi*'s, and must sit
 * below a tenth of every error asserted on. The uniform L2 row is printed
 * and not asserted: its error is within a factor three of that bound.
 *
 * THE COST AGAINST UNIFORM REFINEMENT IS PRINTED AND NOT ASSERTED, since
 * nothing in this tree asserts an adaptive loop against uniform refinement:
 * the last cycle's error beside the uniform L0 -> L1 line, log-log in trace
 * dofs, at the same count.
 */
BOOST_AUTO_TEST_CASE( theEstimatorTracksTheErrorOnTheAdaptiveMirror )
{
	int const order = 2;   // the shipped degree
	int const cycles = 3;
	// AdaptiveRefinement's band at k = 2.
	double const lowest = 5.0;
	double const highest = 40.0;
	// Its own mesh, so that this ctest and RotatingMirrorDriver never generate
	// one file at once.
	std::string const mesh = "File = \"mirror-adaptive.msh\"";
	auto adaptivity = []( int n )
	{
		return "[adaptivity]\nEnabled = true\nMaxIterations = " + std::to_string( n )
		       + "\nTheta = 0.5\n\n";
	};
	auto uniform = [ & ]( int level, int degree )
	{
		return std::vector<std::pair<std::string, std::string>>{
			{ "File = ", mesh + "\nRefinementLevels = " + std::to_string( level ) },
			{ "PolynomialDegree = ", "PolynomialDegree = " + std::to_string( degree ) } };
	};

	// ---- the runs ------------------------------------------------------------
	std::vector<DriverRun> stopped;
	for ( int n = 1; n <= cycles; ++n )
	{
		stopped.push_back( runVariant( "mirror-adaptive-" + std::to_string( n ),
		                               { { "File = ", mesh } }, adaptivity( n ) ) );
		BOOST_TEST_REQUIRE( stopped.back().status == 0,
			"the adaptive mirror stopped at " << n << " cycles exited "
			<< stopped.back().status << ":\n" << stopped.back().log );
	}
	DriverRun const reference = runVariant( "mirror-adaptive-reference", uniform( 2, 4 ) );
	DriverRun const below = runVariant( "mirror-adaptive-reference-L1", uniform( 1, 4 ) );
	// One cycle each, which estimates and does not refine: the uniform meshes' eta.
	std::vector<DriverRun> uniformRuns;
	for ( int level = 1; level <= 2; ++level )
		uniformRuns.push_back( runVariant( "mirror-adaptive-uniform-L"
		                                   + std::to_string( level ), uniform( level, order ),
		                                   adaptivity( 1 ) ) );
	for ( DriverRun const *run : std::vector<DriverRun const *>{
		&reference, &below, &uniformRuns[ 0 ], &uniformRuns[ 1 ] } )
		BOOST_TEST_REQUIRE( run->status == 0, run->stem << " exited " << run->status
		                    << ":\n" << run->log );

	// ---- the reference, and its own error -----------------------------------
	StoredAnswer const fine = readAnswer( reference.stem );
	StoredAnswer const coarse = readAnswer( below.stem );
	FluxDifference const step = fluxDifference( coarse, fine, 4 );
	// The flux converges with psi_h, capped by the inner edge at
	// potentialRate( 4 ) = 3.5 (MR-4c), so this over-states the reference's
	// error, as referenceError() does for psi*.
	double const referenceBound = step.error/( std::pow( 2.0, potentialRate( 4 ) ) - 1.0 );

	// ---- each answer against it ---------------------------------------------
	struct Row
	{
		std::string name;
		int elements = 0;
		int trace = 0;
		int marked = 0;
		double eta = 0.0;
		FluxDifference truth;
	};
	// @a cycle of @a run's table, against the reference.
	auto measured = [ & ]( std::string const &name, DriverRun const &run,
	                       std::size_t cycle )
	{
		std::vector<CycleRow> const table = cycleTable( run.log );
		BOOST_TEST_REQUIRE( table.size() == cycle + 1, "no last cycle " << cycle
		                    << " in " << run.stem << "'s table:\n" << run.log );
		StoredAnswer const answer = readAnswer( run.stem );
		Row row;
		row.name = name;
		row.elements = answer.mesh->GetNE();
		row.trace = table[ cycle ].trace;
		row.eta = table[ cycle ].eta;
		row.truth = fluxDifference( answer, fine, order );
		BOOST_TEST_REQUIRE( row.elements == table[ cycle ].elements,
			run.stem << " wrote " << row.elements << " elements where its last "
			"cycle had " << table[ cycle ].elements );
		return row;
	};
	std::vector<CycleRow> const longRun = cycleTable( stopped.back().log );
	BOOST_TEST_REQUIRE( longRun.size() == static_cast<std::size_t>( cycles ),
	                    "no " << cycles << "-row table in:\n" << stopped.back().log );
	std::vector<Row> rows;
	for ( std::size_t c = 0; c < longRun.size(); ++c )
	{
		rows.push_back( measured( std::to_string( c ), stopped[ c ], c ) );
		BOOST_TEST_REQUIRE( rows.back().elements == longRun[ c ].elements,
			"the run stopped at " << c + 1 << " cycles ends on " << rows.back().elements
			<< " elements where the long run's cycle " << c << " had "
			<< longRun[ c ].elements << ", so its answer is not that cycle's" );
		// What the long run marked, since a stopped run marks nothing last.
		rows.back().marked = longRun[ c ].marked;
	}
	std::vector<Row> uniformRows;
	for ( std::size_t i = 0; i < uniformRuns.size(); ++i )
		uniformRows.push_back( measured( "L" + std::to_string( i + 1 ), uniformRuns[ i ], 0 ) );

	std::printf( "\n  THE ESTIMATOR ON THE ADAPTIVE MIRROR, k = %d, Doerfler 0.5, against "
	             "L2 k = 4 ( %d elements )\n", order, fine.mesh->GetNE() );
	std::printf( "    the reference's own flux error, bounded from L1 at rate %.1f: %.3e, "
	             "%.3e of || q ||\n", potentialRate( 4 ), referenceBound,
	             referenceBound/step.norm );
	std::printf( "    %6s %7s %8s %7s %11s %11s %11s %9s %10s\n", "cycle", "elem", "trace",
	             "marked", "eta", "L2(q-q_h)", "relative", "eta/err", "uncovered" );
	auto print = []( Row const &r )
	{
		std::printf( "    %6s %7d %8d %7d %11.4e %11.4e %11.4e %9.3f %10.3e\n",
		             r.name.c_str(), r.elements, r.trace, r.marked, r.eta, r.truth.error,
		             r.truth.error/r.truth.norm, r.eta/r.truth.error,
		             r.truth.uncovered/r.truth.area );
	};
	for ( Row const &r : rows )
		print( r );
	std::printf( "    uniform refinement, one cycle each:\n" );
	for ( Row const &r : uniformRows )
		print( r );

	// The cost: the uniform L0 -> L1 line at the last cycle's trace count.
	Row const &first = rows.front();
	Row const &last = rows.back();
	double const slope = std::log( uniformRows[ 0 ].truth.error/first.truth.error )
	                     /std::log( static_cast<double>( uniformRows[ 0 ].trace )/first.trace );
	double const uniformThere = first.truth.error
	                            *std::pow( static_cast<double>( last.trace )/first.trace, slope );
	double const uniformCount = first.trace
	                            *std::pow( last.truth.error/first.truth.error, 1.0/slope );
	std::printf( "    cost: cycle %d reaches %.3e on %d trace dofs; uniform refinement, "
	             "error ~ dofs^%.2f from L0 to L1, reads %.3e there and needs %.0f dofs "
	             "for %.3e\n", cycles - 1, last.truth.error, last.trace, slope,
	             uniformThere, uniformCount, last.truth.error );
	std::fflush( stdout );

	// ---- what holds ------------------------------------------------------------
	std::vector<Row> asserted = rows;
	asserted.push_back( uniformRows[ 0 ] );
	for ( Row const &r : asserted )
	{
		BOOST_TEST( referenceBound < 0.1*r.truth.error,
			r.name << ": the reference's own error, bounded at " << referenceBound
			<< ", is not below a tenth of the " << r.truth.error << " measured "
			"against it, so the effectivity index is not this answer's" );
		BOOST_TEST( r.truth.uncovered < 1.0e-3*r.truth.area,
			r.name << ": " << r.truth.uncovered << " of Omega_h's " << r.truth.area
			<< " lies outside the reference's, so the true error is not over the "
			"domain eta is" );
		double const index = r.eta/r.truth.error;
		BOOST_TEST( ( index > lowest && index < highest ),
			r.name << ": the effectivity index is " << index << " against "
			"AdaptiveRefinement's ( " << lowest << ", " << highest << " ) at k = 2. "
			"Above it with eta FALLING AT FIRST ORDER under uniform refinement, F "
			"is being read at the remainder and not at psi_p + psi_c; with eta_5 "
			"carrying it and RISING over the cycles, eta_5's datum is the physical "
			"modal sum and not the remainder's. MEASUREMENTS.md M-279" );
	}
	for ( std::size_t c = 1; c < rows.size(); ++c )
	{
		BOOST_TEST( rows[ c ].eta < rows[ c - 1 ].eta,
			"eta went from " << rows[ c - 1 ].eta << " to " << rows[ c ].eta
			<< " at cycle " << c );
		BOOST_TEST( rows[ c ].truth.error < rows[ c - 1 ].truth.error,
			"the flux error went from " << rows[ c - 1 ].truth.error << " to "
			<< rows[ c ].truth.error << " at cycle " << c << " -- eta came down "
			"and the true error did not, which is what a misleading indicator "
			"looks like" );
	}
	// And the loop is adaptive rather than uniform.
	for ( std::size_t c = 0; c + 1 < rows.size(); ++c )
		BOOST_TEST( ( rows[ c ].marked > 0 && rows[ c ].marked < rows[ c ].elements ),
			"cycle " << c << " marked " << rows[ c ].marked << " of "
			<< rows[ c ].elements << " elements" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( device )

/*
 * THE SHIPPED MIRROR ON THE CARD IS THE HOST'S, dof for dof.
 *
 * MEASUREMENTS.md M-276 raced this file under --device cuda through
 * tools/geq-benchmark/race_geq.py and found device double reproducing the
 * host to 2.0e-15 in psi* at geq's nodes, at L2, k = 4. Nothing asserted it,
 * and DeviceParity's cases are tokamaks: none has a box support, two labels
 * at points or a subtracted filament, which are the paths this file is the
 * only route to. So the shipped file, through meq-run, twice --
 *
 *   --device cpu and --device cuda, OMP_NUM_THREADS = 1 on BOTH, which is
 *   DeviceParity's reason: the driver refuses threaded assembly above one
 *   thread under a Device, and pinning the host too leaves the device the
 *   only difference;
 *
 * and the exact psi_h written by each compared over every dof, with the
 * Newton count beside it (M-187: a device can be right on more work).
 *
 * THE BOUND IS THE ONE THE DRIVER-AGAINST-LIBRARY COMPARISON ABOVE USES,
 * 1e-12 of the largest value, on psi_h and on q_h. The two arms are not
 * bit-identical and are not run to be -- the card takes the driver's device
 * defaults, cuDSS and the batched local factorisation -- and the card is not
 * even bit-reproducible run to run, so a printed-digit comparison is the
 * wrong instrument; what separates the two is round-off in two
 * factorisations of one system. M-276 read 2.0e-15 in psi* against a flux of
 * order 7e-02, 3e-14 of it; three takes here read 4.8e-14 to 1.2e-13 of
 * max |psi_h| and 5.7e-14 to 1.0e-13 of max |q_h| (M-280). So the bound
 * sits a factor eight above the worst take and twelve orders below the
 * defects the device path has had: M-206 was a factor of two. In single it takes the
 * driver comparison's floor too, unmeasured on a card here.
 *
 * SKIPPED LOUDLY, NOT FAILED, without CUDA in the build or a card that MFEM
 * configures -- read off MFEM's own `Device configuration:` line, which is
 * DeviceParity's test for the same question.
 */
BOOST_AUTO_TEST_CASE( theShippedMirrorOnTheCardIsTheHosts )
{
#ifndef MFEM_USE_CUDA
	BOOST_TEST_MESSAGE( "  *** SKIPPED: this build has no CUDA" );
	std::printf( "\n  THE MIRROR ON THE CARD: SKIPPED -- this build has no CUDA\n" );
	return;
#else
	// Its own mesh, as the adaptive suite has, so the three ctests on this
	// file never generate one at once.
	std::vector<std::pair<std::string, std::string>> const lines =
		{ { "File = ", "File = \"mirror-device.msh\"" } };
	DriverRun const host = runVariant( "mirror-device-cpu", lines, "", "--device cpu",
	                                   "OMP_NUM_THREADS=1" );
	BOOST_TEST_REQUIRE( host.status == 0,
		"the mirror under --device cpu exited " << host.status << ":\n" << host.log );
	DriverRun const card = runVariant( "mirror-device-cuda", lines, "", "--device cuda",
	                                   "OMP_NUM_THREADS=1" );
	if ( !says( card.log, "Device configuration: cuda" ) )
	{
		std::printf( "\n  THE MIRROR ON THE CARD: SKIPPED -- built with CUDA, no device "
		             "configured (exit %d)\n", card.status );
		BOOST_TEST_MESSAGE( "  *** SKIPPED: built with CUDA, and MFEM configured no "
		                    "cuda device on this machine" );
		return;
	}
	BOOST_TEST_REQUIRE( card.status == 0,
		"the mirror exits " << card.status << " under --device cuda and 0 on the "
		"host:\n" << card.log );

	StoredAnswer const onHost = readAnswer( host.stem );
	StoredAnswer const onCard = readAnswer( card.stem );
	BOOST_TEST_REQUIRE( onHost.mesh->GetNE() == onCard.mesh->GetNE() );
	BOOST_TEST_REQUIRE( onHost.psi->Size() == onCard.psi->Size() );
	double worst = 0.0;
	double largest = 0.0;
	for ( int i = 0; i < onHost.psi->Size(); ++i )
	{
		worst = std::max( worst, static_cast<double>(
			std::abs( ( *onCard.psi )( i ) - ( *onHost.psi )( i ) ) ) );
		largest = std::max( largest, static_cast<double>( std::abs( ( *onHost.psi )( i ) ) ) );
	}
	double fluxWorst = 0.0;
	double fluxLargest = 0.0;
	for ( int i = 0; i < onHost.flux->Size(); ++i )
	{
		fluxWorst = std::max( fluxWorst, static_cast<double>(
			std::abs( ( *onCard.flux )( i ) - ( *onHost.flux )( i ) ) ) );
		fluxLargest = std::max( fluxLargest,
		                        static_cast<double>( std::abs( ( *onHost.flux )( i ) ) ) );
	}
	int const hostSteps = static_cast<int>( numberAfter( host.log, "converged in " ) );
	int const cardSteps = static_cast<int>( numberAfter( card.log, "converged in " ) );
	double const hostAxis = numberAfter( host.log, "     psi_inner = " );
	double const cardAxis = numberAfter( card.log, "     psi_inner = " );

	std::printf( "\n  THE MIRROR ON THE CARD against --device cpu, %d elements, one host "
	             "thread each\n", onHost.mesh->GetNE() );
	std::printf( "    psi_h: %.3e over %d dofs, %.3e of max |psi_h| %.3e\n", worst,
	             onHost.psi->Size(), worst/largest, largest );
	std::printf( "    q_h:   %.3e over %d dofs, %.3e of max |q_h| %.3e\n", fluxWorst,
	             onHost.flux->Size(), fluxWorst/fluxLargest, fluxLargest );
	std::printf( "    Newton %d on the host, %d on the card; psi_inner printed %.6e and %.6e\n",
	             hostSteps, cardSteps, hostAxis, cardAxis );
	std::fflush( stdout );

	double bound = meq::tests::scalarAgreementFloor( 1.0e-12*largest, largest );
	double fluxBound = meq::tests::scalarAgreementFloor( 1.0e-12*fluxLargest, fluxLargest );
#ifdef MFEM_USE_SINGLE
	// The driver comparison's own floor in single, the extension's
	// eps/fd_step of the field it moves (theShippedMirrorSolvesThroughTheDriver).
	double const extension = 4.0*static_cast<double>(
		std::numeric_limits<mfem::real_t>::epsilon() )/meq::extensionFdStep();
	bound = std::max( bound, extension*largest );
	fluxBound = std::max( fluxBound, extension*fluxLargest );
	// AND psi_h IS A REMAINDER, SO ITS ROUND-OFF IS THE PHYSICAL FLUX'S. The
	// conductors are out of the mesh, psi_h = psi - psi_c is about 1e-3 at its
	// largest, and the terms it is made from are the size of psi, 1.73. On a
	// card-against-host pair the physical flux differs by 1.2e-07, which is
	// 0.6 eps of it, and psi_h by up to 1.1e-07, which against psi_h's own size
	// is 550 to 900 eps and failed one run in three against the floor above.
	if ( onHost.total )
	{
		double totalLargest = 0.0;
		for ( int i = 0; i < onHost.total->Size(); ++i )
			totalLargest = std::max( totalLargest,
			                         static_cast<double>( std::abs( ( *onHost.total )( i ) ) ) );
		bound = std::max( bound, meq::tests::scalarAgreementFloor( 0.0, totalLargest ) );
	}
#endif
	BOOST_TEST( worst < bound,
		"psi_h under --device cuda differs from the host's by " << worst << " of a "
		"largest " << largest << ", above round-off in two factorisations. Walk it "
		"with MEQ_STAGE_CHECKSUM=1 on both arms and diff: the first stage that "
		"disagrees names it (CLAUDE_DEVICE.md)" );
	BOOST_TEST( fluxWorst < fluxBound,
		"q_h under --device cuda differs from the host's by " << fluxWorst
		<< " of a largest " << fluxLargest );
	BOOST_TEST( cardSteps == hostSteps,
		cardSteps << " Newton steps on the card against " << hostSteps << " on the "
		"host -- the answer may be right and the work is not (MEASUREMENTS.md M-187)" );
#endif
}

BOOST_AUTO_TEST_SUITE_END()
#endif

/*
 * THE MIDPLANE GAUGE ON THE MIRROR. MIRROR-PLAN.md section 5.
 *
 * meq::RotatingSource::Gauge::Midplane pins phi_0 where each flux surface
 * crosses z = 0, R_0( psi ), which the solver re-tabulates once per Newton
 * iteration and freezes within it. Three things are asked of it here, one
 * ctest -- RotatingMirrorMidplane -- for the three:
 *
 *   1. THE SAME EQUILIBRIUM IN BOTH GAUGES. geq's problem, given to MEQ in the
 *      fixed-radius gauge as MR-4 gives it, is transferred to the midplane
 *      gauge on MEQ's own best answer and solved again: the two must be one
 *      equilibrium to within MEQ's discretisation error, and the difference
 *      must fall as the mesh refines. MEQ-only, so it needs nothing of geq's.
 *   2. AGAINST geq, IN THE MIDPLANE GAUGE. export_rotating.py's midplane
 *      tables, the same plasma transferred on geq's own 513^2 crossings, and
 *      MR-4d's matrix against geq's direct solves -- to M-263's resolution.
 *   3. THE CASE THE FIXED RADIUS CANNOT POSE: electrons, deuterium and carbon
 *      at M = 6, whose fixed-radius densities span hundreds of decades. It
 *      solves, its n_s0 are the physical midplane densities, and it
 *      self-converges as the mesh refines.
 */
namespace
{
	/// The two densities of the reference's two-species plasma, tabulated on
	/// the curve @a crossing in @a gauge, at the knots of the fixed-radius
	/// tables -- which the crossing must share. Two species share one
	/// exponent, the closed form's
	///
	///     A = C( Psi )( R^2 - R_ref^2 )/2,
	///     C = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ),
	///
	/// so the density where a surface meets the curve is n_s0 e^A at R = R_c,
	/// and its slope is the chain rule through C and R_c:
	///
	///     d/dPsi ( n e^A ) = e^A [ n' + n( C' ( R_c^2 - R_ref^2 )/2 + C R_c R_c' ) ].
	MirrorProfiles transferOnto( MirrorProfiles const &fixed, meq::SplineProfile const &crossing,
	                             meq::RotatingSource::Gauge gauge )
	{
		BOOST_TEST_REQUIRE( fixed.species.size() == 2u );
		BOOST_TEST_REQUIRE( ( fixed.gauge == meq::RotatingSource::Gauge::FixedRadius ) );
		auto const *electron = dynamic_cast<meq::SplineProfile const *>(
			fixed.species[ 0 ].density.get() );
		BOOST_TEST_REQUIRE( electron != nullptr );
		std::vector<meq::Knot> const &knots = electron->knots();
		BOOST_TEST_REQUIRE( crossing.knots().size() == knots.size() );

		meq::Species const &one = fixed.species[ 0 ];
		meq::Species const &two = fixed.species[ 1 ];
		meq::Profile const &omega = *fixed.omega;
		double const k = one.charge*two.mass - two.charge*one.mass;
		double const rRef2 = fixed.referenceRadius*fixed.referenceRadius;

		std::vector<meq::Knot> table;
		for ( std::size_t i = 0; i < knots.size(); ++i )
		{
			double const x = knots[ i ].psi;
			double const w = omega( x )*omega( x );
			double const wPrime = 2.0*omega( x )*omega.prime( x );
			double const d = one.charge*( *two.temperature )( x ) - two.charge*( *one.temperature )( x );
			double const dPrime = one.charge*two.temperature->prime( x ) - two.charge*one.temperature->prime( x );
			double const c = k*w/d;
			double const cPrime = k*( wPrime/d - w*dPrime/( d*d ) );

			meq::Knot const &r = crossing.knots()[ i ];
			BOOST_TEST_REQUIRE( r.psi == x );
			double const half = ( r.value*r.value - rRef2 )/2.0;
			double const factor = std::exp( c*half );
			double const n = ( *one.density )( x );
			double const nPrime = one.density->prime( x );
			table.push_back( meq::Knot{ x, n*factor,
				factor*( nPrime + n*( cPrime*half + c*r.value*r.derivative ) ) } );
		}

		MirrorProfiles out = fixed;
		out.species[ 0 ].density = std::make_shared<meq::SplineProfile const>( std::move( table ) );
		out.species[ 1 ].density = meq::neutralisingDensity( out.species, 1 );
		out.referenceRadius = 0.0;
		out.gauge = gauge;
		return out;
	}

	/// @a run's own midplane crossings at its answer, tabulateMidplaneCrossing()
	/// on the solver's midplane field -- the same reading the solver gives a
	/// midplane-gauged source at every iterate -- on @a intervals.
	std::shared_ptr<meq::SplineProfile const> answerCrossing( MirrorRun &run, unsigned int intervals )
	{
		meq::GradShafranovSolver &solver = *run.solver;
		double const axis = solver.psiAxis();
		double const boundary = solver.psiBoundary();
		double rMin = 0.0, rMax = 0.0;
		meq::MidplaneField const physical = solver.midplaneField( rMin, rMax );
		meq::MidplaneField const normalised = [ & ]( double radius, double &x, double &dXdR )
		{
			double psi = 0.0, slope = 0.0;
			if ( !physical( radius, psi, slope ) )
				return false;
			x = ( psi - boundary )/( axis - boundary );
			dXdR = slope/( axis - boundary );
			return true;
		};
		return meq::tabulateMidplaneCrossing( normalised, rMin, rMax, 0.0, 1.0, intervals );
	}

	/// The knot count less one of the reference's fixed-radius tables.
	unsigned int tableIntervals( MirrorProfiles const &fixed )
	{
		auto const *electron = dynamic_cast<meq::SplineProfile const *>(
			fixed.species[ 0 ].density.get() );
		BOOST_TEST_REQUIRE( electron != nullptr );
		BOOST_TEST_REQUIRE( electron->uniformKnots() );
		BOOST_TEST_REQUIRE( electron->knots().front().psi == 0.0 );
		BOOST_TEST_REQUIRE( electron->knots().back().psi == 1.0 );
		return static_cast<unsigned int>( electron->knots().size() - 1 );
	}

	/// The reference's plasma in the midplane gauge, on @a run's own crossings
	/// at its answer.
	MirrorProfiles transferToMidplane( MirrorRun &run, MirrorProfiles const &fixed )
	{
		return transferOnto( fixed, *answerCrossing( run, tableIntervals( fixed ) ),
		                     meq::RotatingSource::Gauge::Midplane );
	}

	/// The fixed-radius answer the midplane problem is transferred on, and
	/// that transfer: MR-4c's best run, level 2 at k = 4, whose own error
	/// M-263 bounds at 2.4e-08.
	struct Transferred
	{
		MirrorProfiles profiles;
		double largestTransfer = 0.0;  ///< max | ln( n_mid/n_fixed ) | over the knots
	};

	Transferred const &transferredOnMeq()
	{
		static Transferred const t = []
		{
			// The suite running claims the ladder's stem; the midplane suite
			// unless another, in its own process, already has.
			if ( ladderStemInUse == ladderStem )
				ladderStemInUse = midplaneStem;
			Options o;
			o.level = bestLevel;
			o.order = bestOrder;
			o.stem = ladderStemInUse;
			std::unique_ptr<MirrorRun> const run = runMirror( o );
			BOOST_TEST_REQUIRE( run->converged, run->failure );
			Transferred out;
			MirrorProfiles const fixed = mirrorProfiles();
			out.profiles = transferToMidplane( *run, fixed );
			// The same run is MR-4c's reference(), so it is filed rather than
			// solved again.
			sampledCache().emplace( std::make_tuple( bestLevel, bestOrder, true ),
			                        reduce( run, o ) );
			for ( double x = 0.0; x <= 1.0; x += 1.0/256.0 )
				out.largestTransfer = std::max( out.largestTransfer,
					std::abs( std::log( ( *out.profiles.species[ 0 ].density )( x )
					                    /( *fixed.species[ 0 ].density )( x ) ) ) );
			return out;
		}();
		return t;
	}

	void printHistoryLine( char const *label, Sampled const &s )
	{
		std::printf( "    %-24s %2d steps, %2d crossings tabulated:", label,
		             s.newtonSteps, s.midplaneRefreshes );
		for ( double r : s.history )
			std::printf( " %.2e", r );
		std::printf( "\n" );
	}

	// ------------------------------------------------------------------
	// The three-species plasma at M = 6
	// ------------------------------------------------------------------

	/// f( Psi ) = edge + ( peak - edge ) C psi_n^3 ( 1 - psi_n )^4 with
	/// psi_n = 1 - Psi, export_rotating.py's skew shape: flat to order 2 at the
	/// inner label and 3 at the outer, and exact at all three derivative
	/// levels. Clamped outside [ 0, 1 ], where the box gates the source anyway.
	class SkewProfile : public meq::Profile
	{
		public:
			SkewProfile( double peakIn, double edgeIn ) : peak( peakIn ), edge( edgeIn ) {}

			double operator()( double psi ) const override
			{
				double const x = 1.0 - clamp( psi );
				return edge + ( peak - edge )*shape*x*x*x*std::pow( 1.0 - x, 4 );
			}

			double prime( double psi ) const override
			{
				if ( psi < 0.0 || psi > 1.0 )
					return 0.0;
				double const x = 1.0 - psi;
				// d/dPsi = -d/dx.
				return -( peak - edge )*shape*x*x*std::pow( 1.0 - x, 3 )*( 3.0 - 7.0*x );
			}

			double doublePrime( double psi ) const override
			{
				if ( psi < 0.0 || psi > 1.0 )
					return 0.0;
				double const x = 1.0 - psi;
				double const y = 1.0 - x;
				double const u = 2.0*x*y*y*y*( 3.0 - 7.0*x ) - 3.0*x*x*y*y*( 3.0 - 7.0*x )
				                 - 7.0*x*x*y*y*y;
				return ( peak - edge )*shape*u;
			}

		private:
			static double clamp( double psi ) { return std::min( 1.0, std::max( 0.0, psi ) ); }
			double peak, edge;
			static constexpr double shape = 1.0/( ( 27.0/343.0 )*( 256.0/2401.0 ) );
	};

	double const elementaryCharge = 1.602176634e-19;
	double const electronMass = 9.1093837139e-31;
	double const deuteronMass = 3.3435837768e-27;
	/// C6+: twelve atomic mass units less six electrons.
	double const carbonMass = 12.0*1.66053906892e-27 - 6.0*electronMass;

	/// M-167's Mach number, M^2 = m_D omega^2 R^2/T_e, set to 6 at the
	/// profiles' peak on the radius MR-4's peak crosses the midplane at, and
	/// re-read on the solved state at psi_n = 0.5, where M-167 defines it.
	double const machTarget = 6.0;
	double const machRadius = 0.24;
	/// The drift of M-167's two-species omega, 1.30x across the plasma, so
	/// the labels rotate nearly as fast as the core -- which is what puts the
	/// fixed-radius exponents of section 7 at the labels.
	double const omegaDrift = 1.30;
	/// The midplane densities fall to 1e-6 of their peak at the labels. The
	/// box's precondition is F vanishing there, and under this gauge a flat n_s0
	/// still carries -n m_s omega^2 R_0 R_0'/T_s at a label, where omega does
	/// not vanish: so it is the density that does, to a millionth.
	double const densityFloor = 1.0e-6;
	/// The electrons' midplane peak, MEQ_MIDPLANE_DENSITY to sweep it.
	double const densityPeak = std::getenv( "MEQ_MIDPLANE_DENSITY" )
		? std::atof( std::getenv( "MEQ_MIDPLANE_DENSITY" ) ) : 5.0e18;

	/// At @a peak, the electrons' midplane peak, densityPeak unless given, in
	/// @a gauge, and with the densities falling to @a floor of it at the labels.
	MirrorProfiles threeSpecies( double peak = densityPeak,
		meq::RotatingSource::Gauge gauge = meq::RotatingSource::Gauge::Midplane,
		double floor = densityFloor )
	{
		double const tElectronPeak = 1.28e3*elementaryCharge;
		double const omegaPeak = machTarget*std::sqrt( tElectronPeak/deuteronMass )/machRadius;

		auto temperature = [ & ]( double peakEv )
		{
			return std::make_shared<SkewProfile const>( peakEv*elementaryCharge,
			                                            0.1*peakEv*elementaryCharge );
		};

		MirrorProfiles out;
		out.omega = std::make_shared<SkewProfile const>( omegaPeak, omegaPeak/omegaDrift );
		out.referenceRadius = 0.0;
		out.gauge = gauge;
		out.species.resize( 3 );
		out.species[ 0 ].mass = electronMass;
		out.species[ 0 ].charge = -1.0;
		out.species[ 0 ].temperature = temperature( 1.28e3 );
		out.species[ 0 ].density = std::make_shared<SkewProfile const>( peak, floor*peak );
		out.species[ 1 ].mass = deuteronMass;
		out.species[ 1 ].charge = 1.0;
		out.species[ 1 ].temperature = temperature( 1.53e3 );
		out.species[ 2 ].mass = carbonMass;
		out.species[ 2 ].charge = 6.0;
		out.species[ 2 ].temperature = temperature( 1.53e3 );
		out.species[ 2 ].density = std::make_shared<SkewProfile const>( 0.01*peak, floor*0.01*peak );
		out.species[ 1 ].density = meq::neutralisingDensity( out.species, 1 );
		return out;
	}

	char const *const threeNames[] = { "electron", "deuterium", "carbon" };
}

BOOST_AUTO_TEST_SUITE( midplane )

/*
 * 1. ONE EQUILIBRIUM, TWO GAUGES.
 *
 * geq's problem in the fixed-radius gauge is MR-4's. Transferred to the
 * midplane gauge on MEQ's best fixed-radius answer -- the closed form's shared
 * exponent at that answer's own crossings -- it is the SAME plasma, so a
 * midplane solve must land on the same field. It cannot land exactly: the
 * midplane solve tabulates its crossings on its OWN iterate, at its own
 * discretisation error, where the transfer used the best answer's. So the
 * difference is of the order of MEQ's own error at that level and falls with
 * it, which is the assertion; a wrong Jacobian term, a wrong sign on the
 * crossing's slope or a transfer the gauge did not mean would all leave a
 * difference that stays.
 *
 * k = 3 over levels 0..2, MR-4c's asymptotic rung.
 */
BOOST_AUTO_TEST_CASE( theMirrorPosedInBothGaugesIsOneEquilibrium )
{
	Transferred const &t = transferredOnMeq();
	std::printf( "\n  ONE EQUILIBRIUM, TWO GAUGES: the reference's densities transferred to the "
	             "midplane on MEQ's level %d, k = %d answer; | ln( n_mid/n_fixed ) | reaches %.4f\n",
	             bestLevel, bestOrder, t.largestTransfer );
	BOOST_TEST( t.largestTransfer > 0.1,
		"the transfer moves the densities by only " << t.largestTransfer << " in the log, so "
		"the two gauges' tables are nearly one table and this case cannot tell them apart" );

	int const order = 3;
	std::vector<Difference> gauge, own;
	std::vector<double> scale;
	for ( int level : { 0, 1, 2 } )
	{
		Sampled const &fixed = sampled( level, order, true );
		Options o;
		o.level = level;
		o.order = order;
		o.stem = midplaneStem;
		o.profiles = &t.profiles;
		Sampled const midplane = sample( o );

		std::printf( "\n    level %d, k = %d, %d elements\n", level, order, midplane.elements );
		printHistoryLine( "fixed radius", fixed );
		printHistoryLine( "midplane", midplane );
		BOOST_TEST( midplane.midplaneRefreshes == 1,
			"the crossings were tabulated " << midplane.midplaneRefreshes << " times over "
			<< midplane.newtonSteps << " Newton steps: once, where the knots start, is the "
			"contract -- they are unknowns after that" );
		BOOST_TEST( fixed.midplaneRefreshes == 0 );
		std::string why;
		if ( tailIsReadable( midplane.history ) )
			BOOST_TEST( quadraticTail( midplane.history, why ), "level " << level << ": " << why
				<< "the midplane gauge's Newton is not quadratic, so its Jacobian is not the residual's" );

		gauge.push_back( difference( midplane.psiStar, fixed.psiStar ) );
		own.push_back( difference( fixed.psiStar, reference().psiStar ) );
		scale.push_back( std::abs( fixed.psiAxis - fixed.psiBoundary ) );
		printDifference( "midplane - fixed radius, psi*", gauge.back() );
		printDifference( "fixed radius - best, psi*", own.back() );
		std::printf( "    labels: psi_ax %+.3e, psi_bnd %+.3e between the gauges\n",
		             midplane.psiAxis - fixed.psiAxis, midplane.psiBoundary - fixed.psiBoundary );
	}
	std::printf( "\n    rates of midplane - fixed radius, max: %.3f then %.3f; rms: %.3f then %.3f\n",
	             std::log2( gauge[ 0 ].max/gauge[ 1 ].max ), std::log2( gauge[ 1 ].max/gauge[ 2 ].max ),
	             std::log2( gauge[ 0 ].rms/gauge[ 1 ].rms ), std::log2( gauge[ 1 ].rms/gauge[ 2 ].rms ) );
	std::fflush( stdout );

	for ( std::size_t i = 0; i < gauge.size(); ++i )
		BOOST_TEST( gauge[ i ].max < own[ i ].max,
			"level " << i << ": the two gauges differ by " << gauge[ i ].max << ", more than "
			"the fixed-radius answer's own distance from the best, " << own[ i ].max
			<< " -- so the difference is not discretisation" );
	// In single the finer two differences sit under the format's floor --
	// 7.6e-08 and 8.2e-08 rms against 16 eps x |psi_bnd|, where double reads
	// 6.7e-09 and 1.0e-10 -- so whether they fall is round-off. A double-build
	// acceptance there, skipped loudly; the bound against the answer's own
	// error above still runs. M-272.
	double const fluxScale = std::abs( sampled( 2, order, true ).psiBoundary );
	if ( !meq::tests::aboveTheFormatsFloor( gauge.back().rms, fluxScale ) )
	{
		std::printf( "    *** SKIPPED: the gauges' difference ends at %.3e rms, under this "
		             "build's floor on psi*, so whether it falls is not measurable. M-272\n",
		             gauge.back().rms );
		return;
	}
	BOOST_TEST( ( gauge[ 0 ].rms > gauge[ 1 ].rms && gauge[ 1 ].rms > gauge[ 2 ].rms ),
		"the gauges' difference does not fall with the mesh: " << gauge[ 0 ].rms << ", "
		<< gauge[ 1 ].rms << ", " << gauge[ 2 ].rms );
}

/*
 * 2. MEQ AGAINST geq IN THE MIDPLANE GAUGE: MR-4d's matrix, on tables
 * export_rotating.py transfers on geq's own 513^2 crossings -- untracked,
 * since they are made from geq's answer, beside it in the runs directory.
 *
 * WHAT THE TRANSFER IS, AND IT IS NOT AN IDENTITY. MR-4a holds geq's N_s as
 * prescribed flux functions, which are geq's ZERO-POTENTIAL densities, not its
 * midplane ones, so the midplane tables are geq's own quasineutrality solved
 * at each surface's crossing on geq's solved state. They therefore carry
 * geq's 513^2 error in R_0 into the problem -- and the matrix is asserted on
 * MR-4d's own terms regardless: falling at every step, at geq's rate, and at
 * 513 below geq's own 257 - 513 difference. MEASURED, it lands at 0.39 of it,
 * against the fixed-radius gauge's 0.44 ( M-263 ).
 */
BOOST_AUTO_TEST_CASE( meqAgreesWithGeqInTheMidplaneGauge )
{
	std::vector<double> const g129 = readGeqPsi( "psi_direct_129.dat" );
	std::vector<double> const g257 = readGeqPsi( "psi_direct_257.dat" );
	std::vector<double> const g513 = readGeqPsi( "psi_direct_513.dat" );
	std::string const directory = rotatingAnswerDirectory();
	for ( char const *name : { "nmid_electron.dat", "nmid_deuterium.dat" } )
		meqtest::geq::openOutput( directory, name );

	MirrorProfiles geqMidplane = mirrorProfiles();
	geqMidplane.species[ 0 ].density = std::make_shared<meq::SplineProfile const>(
		meq::SplineProfile::fromFile( directory + "/nmid_electron.dat" ) );
	geqMidplane.species[ 1 ].density = std::make_shared<meq::SplineProfile const>(
		meq::SplineProfile::fromFile( directory + "/nmid_deuterium.dat" ) );
	geqMidplane.referenceRadius = 0.0;
	geqMidplane.gauge = meq::RotatingSource::Gauge::Midplane;

	// geq's transfer against MEQ's, on the same plasma: both are the physical
	// midplane density, read on two codes' crossings.
	Transferred const &onMeq = transferredOnMeq();
	double tableGap = 0.0;
	for ( double x = 0.0; x <= 1.0; x += 1.0/512.0 )
		tableGap = std::max( tableGap, std::abs( std::log(
			( *geqMidplane.species[ 0 ].density )( x )/( *onMeq.profiles.species[ 0 ].density )( x ) ) ) );
	std::printf( "\n  MR-4d IN THE MIDPLANE GAUGE: geq's midplane table against MEQ's, max "
	             "| ln ratio | %.3e\n", tableGap );

	Options o;
	o.level = bestLevel;
	o.order = bestOrder;
	o.stem = midplaneStem;
	o.profiles = &geqMidplane;
	Sampled const best = sample( o );
	printHistoryLine( "midplane, geq's tables", best );

	std::vector<double> const meq = withoutCorners( best.psiStar );
	Difference const d129 = difference( meq, g129 );
	Difference const d257 = difference( meq, g257 );
	Difference const d513 = difference( meq, g513 );
	Difference const geqCoarse = difference( withoutCorners( g129 ), g257 );
	Difference const geqFine = difference( withoutCorners( g257 ), g513 );
	Difference const fixedRadius = difference( meq, withoutCorners( reference().psiStar ) );
	printDifference( "MEQ midplane - geq direct 129", d129 );
	printDifference( "MEQ midplane - geq direct 257", d257 );
	printDifference( "MEQ midplane - geq direct 513", d513 );
	printDifference( "geq 257 - geq 513", geqFine );
	printDifference( "MEQ midplane - MEQ fixed radius", fixedRadius );
	std::printf( "    rates, max: MEQ - geq %.3f then %.3f; geq's own %.3f\n",
	             std::log2( d129.max/d257.max ), std::log2( d257.max/d513.max ),
	             std::log2( geqCoarse.max/geqFine.max ) );
	std::printf( "    rates, rms: MEQ - geq %.3f then %.3f; geq's own %.3f\n",
	             std::log2( d129.rms/d257.rms ), std::log2( d257.rms/d513.rms ),
	             std::log2( geqCoarse.rms/geqFine.rms ) );
	std::printf( "    at 513: %.3f ( max ) and %.3f ( rms ) of geq's own 257 - 513\n",
	             d513.max/geqFine.max, d513.rms/geqFine.rms );
	std::fflush( stdout );

	double const fluxScale = std::abs( best.psiBoundary );
	if ( !meq::tests::aboveTheFormatsFloor( geqFine.max, fluxScale ) )
	{
		std::printf( "    *** SKIPPED: geq's own 257 - 513 difference is under this build's "
		             "floor; the conversion bound alone is asserted. M-272\n" );
		BOOST_TEST( d513.max < meq::tests::exactAgreementFloor( fluxScale ) );
		return;
	}
	for ( auto const &norm : { std::make_tuple( "max", d129.max, d257.max, d513.max, geqFine.max,
	                                            std::log2( geqCoarse.max/geqFine.max ) ),
	                           std::make_tuple( "rms", d129.rms, d257.rms, d513.rms, geqFine.rms,
	                                            std::log2( geqCoarse.rms/geqFine.rms ) ) } )
	{
		char const *name = std::get<0>( norm );
		double const a = std::get<1>( norm ), b = std::get<2>( norm ), c = std::get<3>( norm );
		BOOST_TEST( ( a > b && b > c ),
			name << ": MEQ under the midplane gauge against geq reads " << a << ", " << b << ", "
			<< c << " at 129, 257, 513 and does not fall, so it is a conversion -- the transfer, "
			"the crossings, the slope of R_0 -- and not geq's resolution" );
		for ( auto const &step : { std::make_pair( a, b ), std::make_pair( b, c ) } )
			BOOST_TEST( std::log2( step.first/step.second ) > std::get<5>( norm ) - 0.30,
				name << ": MEQ under the midplane gauge falls at "
				<< std::log2( step.first/step.second ) << " against geq where geq converges at "
				<< std::get<5>( norm ) );
		BOOST_TEST( c < std::get<4>( norm ),
			name << ": at 513 MEQ under the midplane gauge is " << c << " from geq, above geq's "
			"own 257 - 513 difference " << std::get<4>( norm ) );
	}
}

/*
 * 3. THE CASE THE FIXED RADIUS CANNOT POSE: ELECTRONS, DEUTERIUM AND CARBON
 * AT M = 6 ON THE MIRROR.
 *
 * M-167 section 7 measured why it could not be: the fixed-radius transfer's
 * exponent is linear in the species mass, +187 to -310 for carbon at M = 6, and
 * a double ends at exp( 709.8 ). Those are the phi_0 = 0 bound; with phi_0 in,
 * on THIS solved plasma at the best single R_ref, carbon's fixed-radius table
 * spans about 1e117 -- inside a double, but a Hermite table holds it to
 * ( E'h )^4/384 relative, thousands of knots. Here the densities are given where they
 * are physical -- at each surface's midplane crossing -- and the solve is asked
 * three things:
 *
 *   it converges, from the vacuum, at every level -- at a midplane peak of
 *   5e18 m^-3. At constant tau = 1 the coarsest mesh stalled at 5e19, and that
 *   was the stabilisation's element-local fold, not the gauge; under the
 *   mirror's stabilisation 5e19 converges there, which the vacuum-gauge
 *   suite's case 5 asserts;
 *   its n_s0 range is the physical one: every n_s0 is its species' density
 *   where its surface crosses z = 0, to round-off, on the solved state; and the
 *   fixed-radius tables of the SAME solved plasma, at the best single R_ref,
 *   are recomputed and printed as M-167's table -- in the log, since they do
 *   not fit in a double;
 *   it self-converges: psi* at geq's nodes, level against level.
 */
BOOST_AUTO_TEST_CASE( threeSpeciesAtMachSixSolveInTheMidplaneGauge )
{
	MirrorProfiles const profiles = threeSpecies();
	int const order = 3;

	std::printf( "\n  THREE SPECIES AT M = %.1f, MIDPLANE GAUGE: omega %.4e .. %.4e rad/s, "
	             "carbon at 1%% of the electrons\n", machTarget, ( *profiles.omega )( 0.0 ),
	             ( *profiles.omega )( 4.0/7.0 ) );

	std::vector<Sampled> runs;
	std::unique_ptr<MirrorRun> finest;
	for ( int level : { 0, 1, 2 } )
	{
		Options o;
		o.level = level;
		o.order = order;
		o.stem = midplaneStem;
		o.profiles = &profiles;
		if ( level == 2 )
		{
			finest = runMirror( o );
			runs.push_back( reduce( finest, o ) );
		}
		else
			runs.push_back( sample( o ) );
		std::printf( "    level %d, %d elements:", level, runs.back().elements );
		printHistoryLine( "", runs.back() );
	}

	// SELF-CONVERGENCE.
	Difference const d01 = difference( runs[ 0 ].psiStar, runs[ 1 ].psiStar );
	Difference const d12 = difference( runs[ 1 ].psiStar, runs[ 2 ].psiStar );
	printDifference( "psi*, level 0 - level 1", d01 );
	printDifference( "psi*, level 1 - level 2", d12 );
	std::printf( "    rate: max %.3f, rms %.3f; labels psi_ax %.10e %.10e %.10e\n",
	             std::log2( d01.max/d12.max ), std::log2( d01.rms/d12.rms ),
	             runs[ 0 ].psiAxis, runs[ 1 ].psiAxis, runs[ 2 ].psiAxis );
	// In single the solves stop at the format's floor and the level 1 - 2
	// difference is round-off: 5.4e-08 rms against 16 eps x |psi_bnd|, and
	// 3.2e-07 in max where double reads 1.1e-08. Neither norm of a field under
	// the floor measures the discretisation, so the self-convergence is a
	// double-build acceptance there, skipped loudly; the n_s0 range below still
	// runs. M-272.
	if ( !meq::tests::aboveTheFormatsFloor( d12.rms, std::abs( runs[ 2 ].psiBoundary ) ) )
		std::printf( "    *** SKIPPED: the level 1 - 2 difference is %.3e rms, under this "
		             "build's floor on psi*, so the self-convergence is not measurable. "
		             "M-272\n", d12.rms );
	else
	{
		BOOST_TEST( d12.rms < 0.25*d01.rms,
			"psi* does not self-converge: level to level " << d01.rms << " then " << d12.rms );
		BOOST_TEST( d12.max < 0.25*d01.max,
			"psi* does not self-converge in max: " << d01.max << " then " << d12.max );
	}

	// THE n_s0 RANGE, ON THE SOLVED STATE.
	meq::NormalisedRotatingSource &source = *finest->source;
	meq::GradShafranovSolver &solver = *finest->solver;
	source.setNormalisation( solver.psiAxis(), solver.psiBoundary() );
	meq::RotatingSource const &inner = source.unnormalised();
	meq::Profile const *crossing = inner.midplaneCrossing();
	BOOST_TEST_REQUIRE( crossing != nullptr );

	// Every n_s0 is the density at its crossing: the transfer factor is ONE.
	double worstAtCrossing = 0.0;
	for ( std::size_t s = 0; s < 3; ++s )
		for ( double x = 0.0; x <= 1.0; x += 1.0/128.0 )
			worstAtCrossing = std::max( worstAtCrossing,
				std::abs( inner.densityExponent( s, ( *crossing )( x ), x ) ) );

	// M-167's Mach on the solved state, at psi_n = 0.5.
	double const machSolved = std::sqrt( deuteronMass*std::pow( ( *profiles.omega )( 0.5 ), 2 )
		*std::pow( ( *crossing )( 0.5 ), 2 )/( *profiles.species[ 0 ].temperature )( 0.5 ) );

	// The fixed-radius gauge's exponent at a radius: ln( n_s0^fixed/n_s0^mid ),
	// which is the midplane gauge's own exponent evaluated at R_ref. The best
	// R_ref minimises its largest magnitude, which M-167 section 7 shows is
	// already what the export picks. A radius at which the closure cannot be
	// evaluated on some surface reads as unusable rather than ending the case.
	auto exponent = [ & ]( std::size_t s, double rRef, double x )
	{
		try
		{
			return inner.densityExponent( s, rRef, x );
		}
		catch ( std::runtime_error const & )
		{
			return std::numeric_limits<double>::infinity();
		}
	};
	auto worstAt = [ & ]( double rRef )
	{
		double worst = 0.0;
		for ( std::size_t s = 0; s < 3; ++s )
			for ( double x = 0.0; x <= 1.0; x += 1.0/128.0 )
				worst = std::max( worst, std::abs( exponent( s, rRef, x ) ) );
		return worst;
	};
	double bestRef = 0.1, bestWorst = std::numeric_limits<double>::infinity();
	for ( double r = 0.10; r <= 0.35; r += 0.0005 )
		if ( worstAt( r ) < bestWorst )
		{
			bestWorst = worstAt( r );
			bestRef = r;
		}

	std::printf( "\n    M-167's M on the solved state at psi_n = 0.5: %.3f; R_0 runs %.4f .. %.4f m\n",
	             machSolved, ( *crossing )( 0.0 ), ( *crossing )( 1.0 ) );
	std::printf( "    ln( n_s( R_0 )/n_s0 ) under the midplane gauge, worst over species and "
	             "surfaces: %.3e\n", worstAtCrossing );
	std::printf( "    the SAME plasma's fixed-radius exponent ln( n_s0^fixed/n_s0^mid ) at the best "
	             "R_ref = %.4f m, at psi_n = 0.05 / 0.5 / 0.95:\n", bestRef );
	double carbonDecades = 0.0;
	for ( std::size_t s = 0; s < 3; ++s )
	{
		std::printf( "      %-10s", threeNames[ s ] );
		for ( double psiN : { 0.05, 0.5, 0.95 } )
			std::printf( " %+9.1f", exponent( s, bestRef, 1.0 - psiN ) );
		double lnMin = std::numeric_limits<double>::infinity();
		double lnMax = -lnMin, midMin = lnMin, midMax = -lnMin;
		for ( double x = 0.0; x <= 1.0; x += 1.0/128.0 )
		{
			double const lnMid = std::log( ( *inner.species()[ s ].density )( x ) );
			double const lnFixed = lnMid + exponent( s, bestRef, x );
			lnMin = std::min( lnMin, lnFixed );
			lnMax = std::max( lnMax, lnFixed );
			midMin = std::min( midMin, lnMid );
			midMax = std::max( midMax, lnMid );
		}
		double const decades = ( lnMax - lnMin )/std::log( 10.0 );
		std::printf( "   n_s0 spans 1e%.1f in the fixed-radius gauge and 1e%.1f in the midplane gauge\n",
		             decades, ( midMax - midMin )/std::log( 10.0 ) );
		if ( s == 2 )
			carbonDecades = decades;
		// The midplane gauge's table spans what was prescribed and no more:
		// the density floor's six decades.
		BOOST_TEST( ( midMax - midMin )/std::log( 10.0 ) < 6.0 + 1.0e-6 );
	}
	std::fflush( stdout );

	BOOST_TEST( worstAtCrossing < 1.0e-10,
		"a density under the midplane gauge is not its table's value where its surface "
		"crosses z = 0: the exponent there reads " << worstAtCrossing );
	BOOST_TEST( std::abs( machSolved - machTarget ) < 1.5,
		"the solved state sits at M = " << machSolved << ", not the operating point" );
	BOOST_TEST( carbonDecades > 100.0,
		"carbon's fixed-radius densities span only 1e" << carbonDecades << ", so this case "
		"does not exercise what the midplane gauge is for" );
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * THE VACUUM-MIDPLANE GAUGE ON THE MIRROR. MIRROR-PLAN.md section 5,
 * [M-281](MEASUREMENTS.md#m-281).
 *
 * meq::RotatingSource::Gauge::FixedCurve pins phi_0 on a curve R_g( Psi ) that
 * never moves, and the solver gives it the VACUUM's midplane crossing once,
 * before the first residual. F is then a pointwise function of ( R, z, psi )
 * and the Jacobian is exact with nothing frozen. One ctest,
 * RotatingMirrorVacuumGauge, asks:
 *
 *   1. geq's plasma posed on the vacuum's curve is the fixed-radius plasma,
 *      and Newton is QUADRATIC on it, analytic border columns or differenced;
 *   2. M-278 section 2's frozen vacuum table, again: three quadratic steps
 *      through either route, so a fixed table leaves nothing out of the
 *      Jacobian -- and why the vacuum's crossing is no estimate to check the
 *      midplane gauge's precondition on;
 *   3. three species at M = 6 at two densities, and how far the plasma
 *      moves its crossings from the vacuum's;
 *   4. the true-midplane densities relabelled off a vacuum-gauge answer,
 *      posed again in the midplane gauge, are the same plasma;
 *   5. the three-species ceiling on the coarsest mesh was constant tau's
 *      element-local fold, which the mirror's stabilisation lifts.
 */
namespace
{
	char const *const vacuumStem = "RotatingMirrorVacuumGauge";

	/// The vacuum's crossing R_0vac( Psi ) from the filaments alone: psi_c
	/// along z = 0 and its slope, normalised by psi_c at the two label points.
	/// What GradShafranovSolver::installVacuumGauge() gives a fixed-curve
	/// source, built here without a solver so that profiles can be stated on
	/// it before one exists; theVacuumCurveIsTheConductorsOwnCrossing holds
	/// the two together.
	std::shared_ptr<meq::SplineProfile const> vacuumCrossing( unsigned int intervals )
	{
		meqtest::geq::GeqCase const &p = problem();
		meq::ConductorField const filaments = geqFilaments();
		double const axis = filaments.psi( p.axisPointR, p.axisPointZ );
		double const boundary = filaments.psi( p.boundaryPointR, p.boundaryPointZ );
		meq::MidplaneField const field = [ & ]( double radius, double &x, double &dXdR )
		{
			double gradR = 0.0, gradZ = 0.0;
			filaments.gradPsi( radius, 0.0, gradR, gradZ );
			x = ( filaments.psi( radius, 0.0 ) - boundary )/( axis - boundary );
			dXdR = gradR/( axis - boundary );
			return true;
		};
		return meq::tabulateMidplaneCrossing( field, 0.0, rhoGamma, 0.0, 1.0, intervals );
	}

	/// The fixed-curve source @a profiles would build, carrying @a curve: what
	/// labelCheck() is asked of.
	meq::RotatingSource::LabelCheck labelCheckOn( MirrorProfiles const &profiles,
	                                              std::shared_ptr<meq::Profile const> curve )
	{
		meq::RotatingSource const source( profiles.species, profiles.omega,
			std::make_shared<meq::ConstantProfile const>( 0.0 ), 0.0, machine().mu0,
			meq::RotatingSource::Closure::Automatic,
			profiles.gauge == meq::RotatingSource::Gauge::FixedRadius
				? meq::RotatingSource::Gauge::FixedCurve : profiles.gauge );
		return source.labelCheck( std::move( curve ) );
	}

	void printLabelCheck( char const *label, meq::RotatingSource::LabelCheck const &c )
	{
		std::printf( "    %-44s F at a label %.3e of its peak ( Psi = %g, R = %.4f )\n", label,
		             c.relative, c.label, c.radius );
	}

	/// The last triple's order of a merit history, or NaN.
	double lastOrder( std::vector<double> const &history )
	{
		std::size_t const n = history.size();
		return n >= 4 ? meq::tests::newtonOrder( history[ n - 3 ], history[ n - 2 ], history[ n - 1 ] )
		              : std::numeric_limits<double>::quiet_NaN();
	}
}

BOOST_AUTO_TEST_SUITE( vacuum_gauge )

/*
 * 1. geq's PLASMA ON THE VACUUM'S CURVE IS THE FIXED-RADIUS PLASMA, AND NEWTON
 * IS QUADRATIC ON IT.
 *
 * Two fixed curves -- R = R_ref and R = R_0vac( Psi ) -- describe one plasma
 * with no solve between them: the transfer is the fixed-radius density read at
 * the vacuum's radius, exact at every knot. So the vacuum-gauge solve must land
 * on the fixed-radius answer to round-off rather than to the discretisation,
 * in the fixed radius's own three steps. The solver installs the curve once,
 * from the filaments, and refreshes nothing; and the border columns
 * DIFFERENCED -- which carry whatever moving the labels does to R_g( Psi ),
 * analytic or not -- give the same history, which is the FD check of the
 * columns on the assembled residual.
 *
 * k = 3 over levels 0..2, the midplane suite's rung.
 */
BOOST_AUTO_TEST_CASE( theMirrorOnTheVacuumCurveIsTheFixedRadiusMirrorWithAQuadraticNewton )
{
	ladderStemInUse = vacuumStem;
	MirrorProfiles const fixed = mirrorProfiles();
	std::shared_ptr<meq::SplineProfile const> const vacuum = vacuumCrossing( tableIntervals( fixed ) );
	MirrorProfiles const onVacuum = transferOnto( fixed, *vacuum, meq::RotatingSource::Gauge::FixedCurve );
	meq::RotatingSource::LabelCheck const check = labelCheckOn( onVacuum, vacuum );
	std::printf( "\n  geq'S PLASMA ON THE VACUUM'S CURVE: R_0vac runs %.6f .. %.6f m\n", ( *vacuum )( 0.0 ),
	             ( *vacuum )( 1.0 ) );
	printLabelCheck( "on the vacuum's curve", check );
	BOOST_TEST( check.relative <= meq::labelTolerance );

	int const order = 3;
	for ( int level : { 0, 1, 2 } )
	{
		Sampled const &radius = sampled( level, order, true );
		Options o;
		o.level = level;
		o.order = order;
		o.stem = vacuumStem;
		o.profiles = &onVacuum;
		std::unique_ptr<MirrorRun> run = runMirror( o );
		BOOST_TEST_REQUIRE( run->converged, run->failure );
		BOOST_TEST( run->solver->vacuumGaugeInstalls() == 1 );
		BOOST_TEST( run->solver->midplaneGaugeRefreshes() == 0 );
		meq::Profile const *installed = run->source->unnormalised().midplaneCrossing();
		BOOST_TEST_REQUIRE( installed != nullptr );
		double worstCurve = 0.0;
		for ( double x = 0.0; x <= 1.0; x += 1.0/1024.0 )
			worstCurve = std::max( worstCurve, std::abs( ( *installed )( x ) - ( *vacuum )( x ) ) );
		Sampled const v = reduce( run, o );
		Difference const gauges = difference( v.psiStar, radius.psiStar );
		double const scale = std::abs( radius.psiBoundary );

		std::printf( "\n    level %d, k = %d, %d elements; the installed curve against the filaments' own "
		             "%.3e m\n", level, order, v.elements, worstCurve );
		printHistoryLine( "fixed radius", radius );
		printHistoryLine( "vacuum curve", v );
		printDifference( "vacuum curve - fixed radius, psi*", gauges );
		BOOST_TEST( worstCurve < 1.0e-12, "the solver installed a curve " << worstCurve
			<< " m from the filaments' own crossing, so it is not the vacuum's" );
		BOOST_TEST( v.newtonSteps == radius.newtonSteps );
		std::string why;
		if ( tailIsReadable( v.history ) )
			BOOST_TEST( quadraticTail( v.history, why ), "level " << level << ": " << why
				<< "the vacuum-curve Newton is not quadratic, so its Jacobian is not the residual's" );
		BOOST_TEST( gauges.max <= std::max( 1.0e-9*scale, meq::tests::exactAgreementFloor( scale ) ),
			"level " << level << ": the one plasma on two fixed curves solves to answers "
			<< gauges.max << " apart, which is not round-off" );

		if ( level == 0 )
		{
			Options d = o;
			d.differencedColumns = true;
			Sampled const dv = sample( d );
			printHistoryLine( "vacuum curve, differenced", dv );
			Difference const columns = difference( dv.psiStar, v.psiStar );
			printDifference( "differenced - analytic columns, psi*", columns );
			BOOST_TEST_REQUIRE( dv.history.size() == v.history.size() );
			// A DIFFERENCED COLUMN IS GOOD TO ITS FORMAT'S DIFFERENCING FLOOR,
			// 4 sqrt( eps ) in single, and a step taken with it moves the next
			// merit by that fraction of the residual the step solved. Measured
			// in single, 1.0e-06 after a first merit of 0.548 and 2.5e-09 after
			// 2.4e-05: 2e-06 and 1e-04 of the residual each step solved, against
			// 1.4e-03. Zero in double, where the history agrees to 1e-6 as it
			// stands. M-219.
			double const differencedFloor = meq::tests::differencedAgreementFloor( 0.0 );
			for ( std::size_t i = 0; i < v.history.size(); ++i )
				BOOST_TEST( std::abs( dv.history[ i ] - v.history[ i ] )
				            <= 1.0e-6*v.history[ i ] + 1.0e-15*v.history.front()
				               + ( i > 0 ? differencedFloor*v.history[ i - 1 ] : 0.0 ),
					"step " << i << ": the differenced columns' merit " << dv.history[ i ]
					<< " against the analytic " << v.history[ i ] );
			BOOST_TEST( columns.max <= std::max( 1.0e-9*scale, meq::tests::exactAgreementFloor( scale ) ) );
		}
	}
	std::fflush( stdout );
}

/*
 * 2. M-278 SECTION 2's FROZEN TABLE, AGAIN: QUADRATIC, SO A FIXED TABLE LEAVES
 * NOTHING OUT OF THE JACOBIAN.
 *
 * M-278 says the midplane problem with the vacuum's table frozen throughout
 * converged in 6 steps where the fixed radius takes 3, which read as a term
 * still missing with the table fixed -- the table's dependence on the labels
 * through Psi, say. Reproduced at level 0, k = 3, on the same true-midplane
 * tables, by the fixed curve given the filaments' crossing: 3 QUADRATIC STEPS,
 * so the six are not reproduced and a fixed table leaves nothing out. The
 * retired per-iteration refresh was the whole of the old midplane gauge's
 * linear rate ([M-281](MEASUREMENTS.md#m-281)); the midplane gauge
 * beside it, its crossings border unknowns, is quadratic too.
 *
 * AND THE TABLES SAY WHY THE MIDPLANE GAUGE CANNOT BE CHECKED BEFORE ITS SOLVE.
 * F vanishes at their labels on their own solved crossings, below
 * meq::labelTolerance, and does not on the vacuum's: an estimate of R_0' from
 * the vacuum would refuse a problem that meets the precondition.
 */
BOOST_AUTO_TEST_CASE( aFrozenCrossingIsQuadraticSoAFixedTableLeavesNothingOut )
{
	ladderStemInUse = vacuumStem;
	Transferred const &t = transferredOnMeq();
	std::shared_ptr<meq::SplineProfile const> const vacuum = vacuumCrossing( 2048 );
	MirrorProfiles onCurve = t.profiles;
	onCurve.gauge = meq::RotatingSource::Gauge::FixedCurve;

	int const order = 3;
	Options base;
	base.level = 0;
	base.order = order;
	base.stem = vacuumStem;

	Options curve = base;
	curve.profiles = &onCurve;
	curve.gaugeCurve = vacuum;
	Sampled const c = sample( curve );

	// The midplane gauge's own answer, whose crossings the tables were
	// transferred to be read on: its knots are unknowns, so it is
	// quadratic too.
	Options midplane = base;
	midplane.profiles = &t.profiles;
	std::unique_ptr<MirrorRun> midplaneRun = runMirror( midplane );
	BOOST_TEST_REQUIRE( midplaneRun->converged, midplaneRun->failure );
	std::shared_ptr<meq::SplineProfile const> const ownCrossing = answerCrossing( *midplaneRun, 1024 );
	Sampled const r = reduce( midplaneRun, midplane );

	std::printf( "\n  M-278 SECTION 2, AT THIS COMMIT: the true-midplane tables, level 0, k = %d\n", order );
	printHistoryLine( "fixed curve, the vacuum's", c );
	printHistoryLine( "midplane, knots", r );
	std::string why;
	if ( tailIsReadable( c.history ) )
		BOOST_TEST( quadraticTail( c.history, why ), why );
	if ( tailIsReadable( r.history ) )
		BOOST_TEST( quadraticTail( r.history, why ), why );
	BOOST_TEST( c.newtonSteps <= 3,
		"a fixed table takes " << c.newtonSteps << " steps, so its Jacobian is not the whole one" );

	meq::RotatingSource::LabelCheck const onVacuum = labelCheckOn( t.profiles, vacuum );
	meq::RotatingSource::LabelCheck const onOwn = labelCheckOn( t.profiles, ownCrossing );
	printLabelCheck( "on the vacuum's crossing", onVacuum );
	printLabelCheck( "on the solved crossing", onOwn );
	BOOST_TEST( onOwn.relative <= meq::labelTolerance );
	BOOST_TEST( onVacuum.relative > meq::labelTolerance,
		"the vacuum's crossing reads these tables as meeting the box's precondition too, so it would "
		"serve the midplane gauge as an estimate after all" );
	std::fflush( stdout );
}

/*
 * 3. ELECTRONS, DEUTERIUM AND C6+ AT M = 6 ON THE VACUUM'S CURVE: the midplane
 * suite's profiles, read on R_0vac, at a midplane peak of 5e18 m^-3 and four
 * times it. The profile set meets the box's precondition there. At both
 * densities every level converges in three quadratic steps and the answer
 * self-converges. At constant tau = 1 level 0 at 2e19 did not, which was the
 * stabilisation's element-local fold and not the gauge (case 5,
 * [M-281](MEASUREMENTS.md#m-281), [M-291](MEASUREMENTS.md#m-291)); under the
 * mirror's default it is asserted like the rest
 * ([M-294](MEASUREMENTS.md#m-294)).
 *
 * AND WHAT THE GAUGE COSTS IN RANGE: the transfer to the solved crossing,
 * max | m_s omega^2 ( R_0vac^2 - R_0^2 )/2T_s | per species, against the same
 * solved plasma's fixed-radius exponent at the best single R_ref.
 */
BOOST_AUTO_TEST_CASE( threeSpeciesAtMachSixOnTheVacuumCurve )
{
	ladderStemInUse = vacuumStem;
	std::shared_ptr<meq::SplineProfile const> const vacuum = vacuumCrossing( 2048 );
	int const order = 3;
	for ( double peak : { 5.0e18, 2.0e19 } )
	{
		MirrorProfiles const three = threeSpecies( peak, meq::RotatingSource::Gauge::FixedCurve );
		meq::RotatingSource::LabelCheck const check = labelCheckOn( three, vacuum );
		std::printf( "\n  THREE SPECIES AT M = %.1f ON THE VACUUM'S CURVE, midplane peak %.1e m^-3\n",
		             machTarget, peak );
		printLabelCheck( "the profile set", check );
		BOOST_TEST( check.relative <= meq::labelTolerance );

		std::vector<Sampled> runs;
		for ( int level : { 0, 1, 2 } )
		{
			Options o;
			o.level = level;
			o.order = order;
			o.stem = vacuumStem;
			o.profiles = &three;
			std::unique_ptr<MirrorRun> run = runMirror( o );
			if ( !run->converged )
			{
				std::printf( "    level %d: NOT CONVERGED, merit %.3e after %zu steps: %s\n", level,
				             run->solver->newtonResiduals().back(),
				             run->solver->newtonResiduals().size() - 1, run->failure.c_str() );
				BOOST_TEST( run->converged,
					"three species at " << peak << " do not converge at level " << level << " under the "
					"mirror's stabilisation, which lifts constant tau's fold on level 0 (case 5): "
					<< run->failure );
				continue;
			}
			BOOST_TEST( run->solver->vacuumGaugeInstalls() == 1 );
			std::shared_ptr<meq::SplineProfile const> const own = answerCrossing( *run, 1024 );
			meq::RotatingSource const &inner = run->source->unnormalised();
			double shift = 0.0, fixedWorst = 0.0;
			std::vector<double> centrifugal( 3, 0.0 );
			for ( double x = 0.0; x <= 1.0; x += 1.0/1024.0 )
			{
				double const r0 = ( *own )( x ), rg = ( *vacuum )( x );
				shift = std::max( shift, std::abs( r0 - rg ) );
				double const w = std::pow( ( *three.omega )( x ), 2 );
				for ( std::size_t s = 0; s < 3; ++s )
					centrifugal[ s ] = std::max( centrifugal[ s ], std::abs( three.species[ s ].mass*w
						*( rg*rg - r0*r0 )/( 2.0*( *three.species[ s ].temperature )( x ) ) ) );
			}
			// The best single R_ref, as the midplane suite scans it: the fixed
			// radius's exponent ln( n_s( R_ref )/n_s( R_0 ) ) on this answer, a
			// radius the closure cannot be evaluated at reading as unusable.
			double bestWorst = std::numeric_limits<double>::infinity();
			for ( double rRef = 0.10; rRef <= 0.35; rRef += 0.0025 )
			{
				double worst = 0.0;
				try
				{
					for ( std::size_t s = 0; s < 3; ++s )
						for ( double x = 0.0; x <= 1.0; x += 1.0/128.0 )
							worst = std::max( worst, std::abs( inner.densityExponent( s, rRef, x )
							                                   - inner.densityExponent( s, ( *own )( x ), x ) ) );
				}
				catch ( std::runtime_error const & )
				{
					worst = std::numeric_limits<double>::infinity();
				}
				bestWorst = std::min( bestWorst, worst );
			}
			fixedWorst = bestWorst;
			runs.push_back( reduce( run, o ) );
			std::printf( "    level %d:", level );
			printHistoryLine( "", runs.back() );
			std::printf( "      max | R_0 - R_0vac | %.4e m; max | m omega^2 ( R_0vac^2 - R_0^2 )/2T |: "
			             "electron %.3e, deuterium %.3e, carbon %.3e; the best fixed radius's %.1f\n",
			             shift, centrifugal[ 0 ], centrifugal[ 1 ], centrifugal[ 2 ], fixedWorst );
			std::string why;
			if ( tailIsReadable( runs.back().history ) )
				BOOST_TEST( quadraticTail( runs.back().history, why ), "level " << level << ": " << why );
			BOOST_TEST( runs.back().newtonSteps <= 4 );
			BOOST_TEST( centrifugal[ 2 ] < 1.0e-1*fixedWorst,
				"carbon's transfer to the solved crossing is " << centrifugal[ 2 ] << " against the "
				"fixed radius's " << fixedWorst << ", so the vacuum's curve does not keep the range down" );
		}
		if ( runs.size() == 3 )
		{
			Difference const d01 = difference( runs[ 0 ].psiStar, runs[ 1 ].psiStar );
			Difference const d12 = difference( runs[ 1 ].psiStar, runs[ 2 ].psiStar );
			printDifference( "psi*, level 0 - level 1", d01 );
			printDifference( "psi*, level 1 - level 2", d12 );
			if ( meq::tests::aboveTheFormatsFloor( d12.rms, std::abs( runs[ 2 ].psiBoundary ) ) )
				BOOST_TEST( d12.rms < 0.25*d01.rms );
		}
	}
	std::fflush( stdout );
}

/*
 * 4. THE TRUE-MIDPLANE DENSITIES RELABELLED OFF A VACUUM-GAUGE ANSWER POSE THE
 * SAME PLASMA IN THE MIDPLANE GAUGE: meq::densityOnCurve() on the answer's own
 * crossing -- what the .nc writes -- and the tables solved again at the same
 * level, k = 3, under Gauge::Midplane. geq's two species and the three at
 * 5e18. The difference is the midplane gauge's knots': it pins the potential
 * on their interpolant, within a micron of the crossing the tables were read
 * at, and four times the knots takes it down by more than half; the control is
 * the SAME tables read as midplane densities, which is another plasma.
 * [M-282](MEASUREMENTS.md#m-282), [M-284](MEASUREMENTS.md#m-284).
 */
BOOST_AUTO_TEST_CASE( theRelabelledDensitiesPoseTheSamePlasmaInTheMidplaneGauge )
{
	ladderStemInUse = vacuumStem;
	int const order = 3;
	MirrorProfiles const fixed = mirrorProfiles();
	MirrorProfiles const two = transferOnto( fixed, *vacuumCrossing( tableIntervals( fixed ) ),
	                                         meq::RotatingSource::Gauge::FixedCurve );
	MirrorProfiles const three = threeSpecies( 5.0e18, meq::RotatingSource::Gauge::FixedCurve );
	for ( MirrorProfiles const *p : { &two, &three } )
	{
		Options o;
		o.level = 0;
		o.order = order;
		o.stem = vacuumStem;
		o.profiles = p;
		std::unique_ptr<MirrorRun> run = runMirror( o );
		BOOST_TEST_REQUIRE( run->converged, run->failure );
		meq::RotatingSource const inner = run->source->unnormalised();
		std::shared_ptr<meq::SplineProfile const> const shipped = answerCrossing( *run, 1024 );
		std::shared_ptr<meq::SplineProfile const> const finer = answerCrossing( *run, 4096 );
		Sampled const v = reduce( run, o );
		double scale = 0.0;
		for ( double x : v.plasma )
			if ( std::isfinite( x ) )
				scale = std::max( scale, std::abs( x ) );

		auto relabelled = [ & ]( meq::SplineProfile const &crossing )
		{
			MirrorProfiles out = *p;
			out.gauge = meq::RotatingSource::Gauge::Midplane;
			for ( std::size_t s = 0; s < out.species.size(); ++s )
				if ( s != 1 )
					out.species[ s ].density = meq::densityOnCurve( inner, s, crossing );
			out.species[ 1 ].density = meq::neutralisingDensity( out.species, 1 );
			return out;
		};
		MirrorProfiles const atShipped = relabelled( *shipped );
		MirrorProfiles const atFiner = relabelled( *finer );
		MirrorProfiles asGiven = *p;
		asGiven.gauge = meq::RotatingSource::Gauge::Midplane;

		Options m = o;
		m.profiles = &atShipped;
		Sampled const a = sample( m );
		m.profiles = &atFiner;
		m.gaugeKnots = 128;
		Sampled const b = sample( m );
		m.profiles = &asGiven;
		m.gaugeKnots = 0;
		Sampled const c = sample( m );
		Difference const dShipped = difference( a.psiStar, v.psiStar );
		Difference const dFiner = difference( b.psiStar, v.psiStar );
		Difference const dOther = difference( c.psiStar, v.psiStar );

		std::printf( "\n  RELABELLED AND POSED AGAIN, %zu species, level 0, k = %d: the plasma's own psi* "
		             "reaches %.3e\n", p->species.size(), order, scale );
		printHistoryLine( "vacuum curve", v );
		printHistoryLine( "midplane, relabelled", a );
		printDifference( "relabelled, 32 knots - vacuum curve, psi*", dShipped );
		printDifference( "relabelled, 128 knots - vacuum curve, psi*", dFiner );
		printDifference( "the tables as given - vacuum curve", dOther );
		// In single both knot counts leave the relabelled plasma under the
		// format's floor on psi* -- 2.7e-08 and 2.0e-08 for two species, where
		// 16 eps x | psi_bnd | is 1.3e-07 -- so whether more knots take it down
		// is round-off. Skipped loudly there; the control below still runs.
		// M-272.
		if ( meq::tests::aboveTheFormatsFloor( dFiner.max, v.psiBoundary ) )
			BOOST_TEST( dFiner.max < 0.5*dShipped.max,
				"128 knots leave the relabelled plasma " << dFiner.max << " from the one it was read off, "
				"against 32 knots' " << dShipped.max << ", so the difference is not the knots'" );
		else
			std::printf( "    *** SKIPPED: 128 knots leave the relabelled plasma %.3e from the one it was "
			             "read off, under this build's floor on psi*, so whether the knots set it is not "
			             "measurable. M-272\n", dFiner.max );
		BOOST_TEST( dOther.max > 1.0e2*dShipped.max,
			"the tables as given and the relabelled ones solve within " << dOther.max << " of each other, so "
			"this case cannot tell one plasma from another" );
	}
	std::fflush( stdout );
}

/*
 * 5. THE THREE-SPECIES CEILING ON THE COARSEST MESH WAS CONSTANT tau'S, NOT THE
 * PLASMA'S, AND THE MIRROR'S STABILISATION LIFTS IT.
 *
 * At a constant tau = 1 one element's local HDG block goes singular when the
 * reaction ( 1/R ) dF/dpsi reaches its lowest eigenvalue,
 * ( k + 1 ) tau |dK|/|K|: the discrete problem folds inside that element, near
 * 1.37e19 on level 0 in the vacuum-midplane gauge, so case 3 stalled at 2e19
 * there and converged on level 1, and the midplane gauge did the same.
 * tau = tau0/max( R, R_f ) makes the criterion dF/dpsi < ( k + 1 ) tau0 |dK|/|K|
 * wherever R >= R_f and moves the element's fold past what the plasma allows:
 * the vacuum-midplane problem's own fold, mesh converged near 4.16e19, and in
 * the midplane gauge tau's own, linear in tau0 and near 2.15e20 at the
 * mirror's tau0 = 3. So under the mirror's default -- every case's Options --
 * 2e19 on level 0 converges in BOTH midplane gauges, and 5e19, the density
 * the constant-tau ceiling stalled at on levels 0 and 1, converges on level 0
 * in the midplane gauge: quadratically, in the steps level 1 takes. The
 * control that makes this a statement about tau is constant tau = 1 at 2e19,
 * which must stall. 5e19 has no control, since constant tau stalls there
 * already at 2e19, and none on the vacuum's curve, where it is past the posed
 * problem's fold. [M-291](MEASUREMENTS.md#m-291),
 * [M-293](MEASUREMENTS.md#m-293),
 * [M-294](MEASUREMENTS.md#m-294).
 */
BOOST_AUTO_TEST_CASE( theCoarsestMeshsCeilingWasConstantTausAndTheMirrorsStabilisationLiftsIt )
{
	ladderStemInUse = vacuumStem;
	struct Case
	{
		char const *name;
		meq::RotatingSource::Gauge gauge;
		double peak;
		bool control;
	};
	for ( Case const &c : { Case{ "vacuum-midplane", meq::RotatingSource::Gauge::FixedCurve, 2.0e19, true },
	                        Case{ "midplane, 32 knots", meq::RotatingSource::Gauge::Midplane, 2.0e19, true },
	                        Case{ "midplane, 32 knots", meq::RotatingSource::Gauge::Midplane, 5.0e19, false } } )
	{
		MirrorProfiles const three = threeSpecies( c.peak, c.gauge );
		std::printf( "\n  THREE SPECIES AT M = %.1f, %.1e m^-3, %s gauge, level 0, k = 3\n",
		             machTarget, c.peak, c.name );

		Options o;
		o.level = 0;
		o.order = 3;
		o.stem = vacuumStem;
		o.profiles = &three;
		std::unique_ptr<MirrorRun> const scaled = runMirror( o );
		BOOST_TEST_REQUIRE( scaled->converged, c.name << ": the mirror's stabilisation, tau = " << o.tau
			<< "/max( R, " << o.tauFloorRadius << " m ), does not converge at " << c.peak
			<< " on level 0, so it does not lift constant tau's local fold: " << scaled->failure );
		std::vector<double> const &history = scaled->solver->newtonResiduals();
		std::printf( "    %-32s %2zu steps:", "the mirror's stabilisation", history.size() - 1 );
		for ( double r : history )
			std::printf( " %.2e", r );
		std::printf( "\n" );
		std::string why;
		if ( tailIsReadable( history ) )
			BOOST_TEST( quadraticTail( history, why ), c.name << ": " << why );
		BOOST_TEST( history.size() - 1 <= 4u, c.name << ": " << history.size() - 1 << " steps" );

		if ( !c.control )
			continue;
		Options d = o;
		d.tau = 1.0;
		d.diffusionScaledTau = false;
		std::unique_ptr<MirrorRun> const constant = runMirror( d );
		std::vector<double> const &control = constant->solver->newtonResiduals();
		std::printf( "    %-32s %s after %zu steps, merit %.2e\n", "constant tau = 1, the control",
		             constant->converged ? "CONVERGED" : "stalls", control.size() - 1, control.back() );
		BOOST_TEST( !constant->converged, c.name << ": constant tau = 1 converges at " << c.peak
			<< " on level 0, so this case no longer separates the stabilisation's fold from the "
			"plasma's; re-measure the ceiling (M-291) before moving the density" );
	}
	std::fflush( stdout );
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * THE MIDPLANE GAUGE WITH ITS CROSSINGS AS BORDER UNKNOWNS.
 * MIRROR-PLAN.md section 5; [M-283](MEASUREMENTS.md#m-283).
 *
 * GradShafranovSolver::setMidplaneGaugeKnots() makes R_0 at K knots in Psi
 * unknowns of the bordered Newton, the source's table the C^2 spline of
 * R_0^2/2 through them, and adds one row per knot: the flux along the midplane
 * at R_0k is its surface's. So F is pointwise in ( R, z, psi ) and the knots,
 * and the Jacobian is the whole one. One ctest, RotatingMirrorBorder, for:
 *
 *   1. THE NEW BLOCKS AGAINST THEIR OWN DERIVATIVES: the rows against q and
 *      psi -- exact, since the midplane's flux is linear in the state -- and
 *      against the knot; the rows' meaning at convergence; the columns against
 *      a central difference of the assembled source.
 *   2. QUADRATIC, on geq's two-species plasma, and to one answer from the
 *      vacuum, from that answer perturbed, and from the vacuum's table held;
 *      with the held-table control that says a FIXED table's Jacobian is
 *      complete.
 *   3. THREE SPECIES AT M = 6 AND 5e18, the same.
 *   4. geq's MIRROR UNDER IT, as MR-4d asks.
 *   5. THE REFUSALS.
 *
 * K = 32, the count the knot sweep of [M-283](MEASUREMENTS.md#m-283)
 * puts 170x under the discretisation error at level 2 and at the 1024-interval table's
 * own agreement with the answer.
 */
namespace
{
	char const *const borderStem = "RotatingMirrorBorder";
	int const borderKnots = 32;

	/// geq's two-species plasma posed in the midplane gauge, M-278 section 4:
	/// export_rotating.py's nmid_*.dat, untracked, beside geq's answer.
	MirrorProfiles geqMidplaneProfiles()
	{
		std::string const directory = rotatingAnswerDirectory();
		for ( char const *name : { "nmid_electron.dat", "nmid_deuterium.dat" } )
			meqtest::geq::openOutput( directory, name );
		MirrorProfiles out = mirrorProfiles();
		out.species[ 0 ].density = std::make_shared<meq::SplineProfile const>(
			meq::SplineProfile::fromFile( directory + "/nmid_electron.dat" ) );
		out.species[ 1 ].density = std::make_shared<meq::SplineProfile const>(
			meq::SplineProfile::fromFile( directory + "/nmid_deuterium.dat" ) );
		out.referenceRadius = 0.0;
		out.gauge = meq::RotatingSource::Gauge::Midplane;
		return out;
	}

	/// The table @a run's source holds, resampled as a Hermite table on
	/// tabulateMidplaneCrossing()'s own 1024 intervals.
	std::shared_ptr<meq::Profile const> tableOf( MirrorRun const &run )
	{
		meq::Profile const *table = run.source->unnormalised().midplaneCrossing();
		BOOST_TEST_REQUIRE( table != nullptr );
		return std::make_shared<meq::SplineProfile const>(
			[ table ]( double x ) { return ( *table )( x ); },
			[ table ]( double x ) { return table->prime( x ); }, 1024 );
	}

	/// @a run's answer as a start for another solver on the same mesh: the
	/// whole iterate, the labels, the coefficients and its crossings.
	Start startOf( MirrorRun const &run )
	{
		Start out;
		out.state = run.solver->state();
		out.psiAxis = run.solver->psiAxis();
		out.psiBoundary = run.solver->psiBoundary();
		out.coefficients = run.solver->exteriorCoefficients();
		out.crossing = tableOf( run );
		return out;
	}

	int stepsOf( MirrorRun const &run )
	{
		return static_cast<int>( run.solver->newtonResiduals().size() ) - 1;
	}

	void printMerits( char const *label, MirrorRun const &run )
	{
		std::printf( "    %-40s %2d steps:", label, stepsOf( run ) );
		for ( double r : run.solver->newtonResiduals() )
			std::printf( " %.2e", r );
		std::printf( "\n" );
		std::fflush( stdout );
	}

	/// The largest contraction e_{n+1}/e_n over the steps above 1e-13,
	/// against the next: a quadratic iteration's contraction falls with e_n,
	/// a linear one's stays put. Returned as the last two such ratios.
	std::pair<double, double> lastContractions( std::vector<double> const &h )
	{
		std::vector<double> ratios;
		for ( std::size_t i = 1; i < h.size(); ++i )
			if ( h[ i ] > 1.0e-13 )
				ratios.push_back( h[ i ]/h[ i - 1 ] );
		if ( ratios.size() < 2 )
			return { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
		return { ratios[ ratios.size() - 2 ], ratios.back() };
	}

	double maxDifference( MirrorRun const &a, MirrorRun const &b )
	{
		mfem::GridFunction d( a.solver->potential() );
		d -= b.solver->potential();
		return std::max( { static_cast<double>( d.Normlinf() ),
			std::abs( a.solver->psiAxis() - b.solver->psiAxis() ),
			std::abs( a.solver->psiBoundary() - b.solver->psiBoundary() ) } );
	}

	/// That a solve from another start landed on @a answer's: to 1e-11, where
	/// the two can be told apart from the format. Two whole solves, each
	/// stopped at ten eps of its own reference, meet only to eps times the
	/// trace system's conditioning -- exactAgreementFloor()'s 4096 eps, set
	/// by a solve read through two routes at 1270 eps. Measured in single,
	/// 6.1e-08 to 4.3e-07 on these mirrors, 7 to 55 eps of psi_bnd; under
	/// that floor the claim is skipped loudly, and in double it stands. M-272.
	void landsOnTheAnswer( double lands, MirrorRun const &answer, std::string const &label )
	{
		if ( meq::tests::aboveTheFormatsFloor( lands, answer.solver->psiBoundary(), 4096.0 ) )
			BOOST_TEST( lands < 1.0e-11, label << ": lands " << lands << " from the answer" );
		else
			std::printf( "    *** SKIPPED: %s lands %.3e from the answer, under this build's floor on "
			             "the flux. M-272\n", label.c_str(), lands );
	}
}

BOOST_AUTO_TEST_SUITE( border )

/*
 * 1. EVERY NEW BLOCK AGAINST ITS OWN DERIVATIVE, on a converged bordered solve
 * of geq's plasma at level 0, k = 3.
 *
 * ROWS AGAINST q AND psi. The midplane's flux is psi_h at the first run's start
 * plus R q_R integrated, so it is LINEAR in the state and its covector must
 * reproduce midplaneFieldOf() on any state, not only near the answer: on the
 * answer, and on the answer with every flux and potential entry perturbed.
 * AGAINST THE KNOT, a central difference of the covector's own sum in R against
 * the field's dpsi/dR there -- the corner. AGAINST THE LABELS, the rows'
 * meaning: at convergence each knot's Psi IS its knot's.
 *
 * COLUMNS AGAINST THE SOURCE: assembleGaugeColumns() against a central
 * difference, in each knot, of the assembled source term -- which is what
 * assembleCurrentColumn() returns at a current scale of one -- with the
 * interpolant through the moved knots installed.
 */
BOOST_AUTO_TEST_CASE( theCrossingRowsAndColumnsAreTheirOwnDerivatives )
{
	MirrorProfiles const geq = geqMidplaneProfiles();
	Options o;
	o.level = 0;
	o.order = 3;
	o.stem = borderStem;
	o.profiles = &geq;
	o.gaugeKnots = 16;
	std::unique_ptr<MirrorRun> const run = runMirror( o );
	BOOST_TEST_REQUIRE( run->converged, run->failure );
	meq::GradShafranovSolver &solver = *run->solver;
	std::vector<double> const crossings = solver.midplaneGaugeCrossings();
	BOOST_TEST_REQUIRE( crossings.size() == 16u );
	mfem::Vector const state = solver.state();
	auto const basis = std::make_shared<meq::CardinalSpline const>( 0.0, 1.0, crossings.size() );

	mfem::Vector perturbed( state );
	double const size = state.Normlinf();
	for ( int i = 0; i < perturbed.Size(); ++i )
		perturbed( i ) += 1.0e-3*size*std::sin( 1.7*i + 0.4 );
	meq::MidplaneField const field = solver.midplaneFieldOf( state );
	meq::MidplaneField const moved = solver.midplaneFieldOf( perturbed );

	double worstRow = 0.0, worstPerturbed = 0.0, worstCorner = 0.0, worstMeaning = 0.0;
	double const axis = solver.psiAxis(), boundary = solver.psiBoundary();
	for ( std::size_t k = 0; k < crossings.size(); ++k )
	{
		double const r = crossings[ k ];
		std::vector<int> dofs;
		std::vector<double> weights;
		BOOST_TEST_REQUIRE( solver.midplaneCovector( r, dofs, weights ) );
		auto dot = [ & ]( mfem::Vector const &v, double radius )
		{
			std::vector<int> d;
			std::vector<double> w;
			BOOST_TEST_REQUIRE( solver.midplaneCovector( radius, d, w ) );
			double total = 0.0;
			for ( std::size_t j = 0; j < d.size(); ++j )
				total += w[ j ]*v( d[ j ] );
			return total + run->filaments.psi( radius, 0.0 );
		};
		double psi = 0.0, slope = 0.0, psiMoved = 0.0, slopeMoved = 0.0;
		BOOST_TEST_REQUIRE( field( r, psi, slope ) );
		BOOST_TEST_REQUIRE( moved( r, psiMoved, slopeMoved ) );
		worstRow = std::max( worstRow, std::abs( dot( state, r ) - psi )/std::abs( psi ) );
		worstPerturbed = std::max( worstPerturbed, std::abs( dot( perturbed, r ) - psiMoved )/std::abs( psiMoved ) );
		double const h = 1.0e-6;
		double const difference = ( dot( state, r + h ) - dot( state, r - h ) )/( 2.0*h );
		worstCorner = std::max( worstCorner, std::abs( difference - slope )/std::abs( slope ) );
		worstMeaning = std::max( worstMeaning,
			std::abs( ( psi - boundary )/( axis - boundary ) - basis->knot( k ) ) );
	}

	// The columns, at three knots.
	meq::NormalisedRotatingSource &source = *run->source;
	source.setNormalisation( axis, boundary );
	std::vector<mfem::Vector> columns;
	BOOST_TEST_REQUIRE( solver.assembleGaugeColumns( state, *basis, crossings, columns ) );
	double worstColumn = 0.0;
	bool columnsMeasurable = true;
	for ( std::size_t k : { std::size_t( 0 ), crossings.size()/2, crossings.size() - 1 } )
	{
		double const h = 1.0e-6;
		std::vector<double> up = crossings, down = crossings;
		up[ k ] += h;
		down[ k ] -= h;
		mfem::Vector plus, minus;
		source.setMidplaneCrossing( meq::crossingThroughKnots( basis, up ) );
		solver.assembleCurrentColumn( state, plus );
		source.setMidplaneCrossing( meq::crossingThroughKnots( basis, down ) );
		solver.assembleCurrentColumn( state, minus );
		// THE DIFFERENCE IS READABLE ONLY WHILE WHAT IT DIFFERENCES MOVES BY
		// MORE THAN THE FORMAT'S FLOOR ON THE WHOLE SOURCE: the spline's
		// cardinal functions reach every element, so every entry carries its
		// own eps. In single the step moves the source by far less than that,
		// and the difference reads round-off. M-272.
		columnsMeasurable = columnsMeasurable && meq::tests::aboveTheFormatsFloor(
			2.0*h*static_cast<double>( columns[ k ].Norml2() ), static_cast<double>( plus.Norml2() ) );
		plus -= minus;
		plus /= 2.0*h;
		plus -= columns[ k ];
		worstColumn = std::max( worstColumn, double( plus.Norml2()/columns[ k ].Norml2() ) );
		BOOST_TEST( columns[ k ].Norml2() > 0.0 );
	}
	source.setMidplaneCrossing( meq::crossingThroughKnots( basis, crossings ) );

	std::printf( "\n  THE BORDERED GAUGE'S BLOCKS, level 0, k = 3, %zu knots\n", crossings.size() );
	std::printf( "    row against midplaneFieldOf(): %.3e on the answer, %.3e perturbed by 1e-3\n",
	             worstRow, worstPerturbed );
	std::printf( "    corner against a central difference in the knot: %.3e\n", worstCorner );
	std::printf( "    each knot's Psi against its knot at convergence: %.3e\n", worstMeaning );
	std::printf( "    columns against a central difference of the assembled source: %.3e\n", worstColumn );
	std::fflush( stdout );
	// IN SINGLE: the rows, the meaning and the columns are skipped where they
	// sit under the format's floor, and the corner, a central difference
	// against an assembled slope, is held to the differencing floor, 4 sqrt(
	// eps ). Measured in single: rows 1.2e-10 and 3.0e-10 and the meaning
	// 1.5e-08, against 16 eps = 1.9e-06; the corner 1.8e-04 against 1.4e-03;
	// the columns 3.9, since the knot's step moves the source by less than
	// the source's own round-off. Every bound stands in double. M-219, M-272.
	for ( auto const &check : { std::make_pair( "row on the answer", worstRow ),
	                            std::make_pair( "row perturbed", worstPerturbed ) } )
	{
		if ( meq::tests::aboveTheFormatsFloor( check.second, 1.0 ) )
			BOOST_TEST( check.second < 1.0e-13, check.first );
		else
			std::printf( "    *** SKIPPED: the %s, %.3e, is under this build's floor. M-272\n",
			             check.first, check.second );
	}
	BOOST_TEST( worstCorner < meq::tests::differencedAgreementFloor( 1.0e-7 ) );
	if ( meq::tests::aboveTheFormatsFloor( worstMeaning, 1.0 ) )
		BOOST_TEST( worstMeaning < 1.0e-10 );
	else
		std::printf( "    *** SKIPPED: each knot's Psi, %.3e from its knot, is under this build's "
		             "floor. M-272\n", worstMeaning );
	if ( columnsMeasurable )
		BOOST_TEST( worstColumn < 1.0e-6 );
	else
		std::printf( "    *** SKIPPED: the knot's step moves the assembled source by less than this "
		             "build's floor on it, so its central difference reads round-off. M-272\n" );
}

/*
 * 2. QUADRATIC, AND FROM ANY NEAR START TO THE SAME ANSWER.
 *
 * geq's two-species plasma at k = 3 on levels 0 and 1: from the vacuum, three
 * steps, as the fixed radius takes. At level 0, from the answer perturbed --
 * every entry, the labels and the crossings -- and from the plasma solved on
 * the vacuum's crossing HELD, it lands on the answer from the vacuum in two.
 * The held table converges in three itself, which is the control that says a
 * FIXED table's Jacobian is complete: the term the retired per-iteration
 * refresh left out is the table's pull by the field, and the knots carry it.
 */
BOOST_AUTO_TEST_CASE( theBorderedGaugeIsQuadraticWhereTheFrozenRefreshIsLinear )
{
	MirrorProfiles const geq = geqMidplaneProfiles();
	std::printf( "\n  geq's TWO SPECIES, k = 3, the bordered gauge, %d knots\n", borderKnots );
	for ( int level : { 0, 1 } )
	{
		Options b;
		b.level = level;
		b.order = 3;
		b.stem = borderStem;
		b.profiles = &geq;
		b.gaugeKnots = borderKnots;
		std::unique_ptr<MirrorRun> bordered = runMirror( b );
		BOOST_TEST_REQUIRE( bordered->converged, bordered->failure );
		std::printf( "\n    level %d\n", level );
		printMerits( "bordered, from the vacuum", *bordered );
		BOOST_TEST( stepsOf( *bordered ) <= 3 );
		std::string why;
		if ( tailIsReadable( bordered->solver->newtonResiduals() ) )
			BOOST_TEST( quadraticTail( bordered->solver->newtonResiduals(), why ), why );
		if ( level > 0 )
			continue;

		Options held = b;
		held.holdCrossing = true;
		std::unique_ptr<MirrorRun> vacuumHeld = runMirror( held );
		BOOST_TEST_REQUIRE( vacuumHeld->converged, vacuumHeld->failure );
		printMerits( "the vacuum's table held, from the vacuum", *vacuumHeld );
		BOOST_TEST( stepsOf( *vacuumHeld ) <= 3,
			"a held table took " << stepsOf( *vacuumHeld ) << " steps: a fixed table's "
			"Jacobian is not the whole one" );
		Start fromHeld = startOf( *vacuumHeld );
		fromHeld.crossing = answerCrossing( *vacuumHeld, 1024 );

		// The answer perturbed: every entry by 1e-3 of itself, the labels by
		// 1e-4, the crossings by 1e-4 m.
		Start perturbed = startOf( *bordered );
		for ( int i = 0; i < perturbed.state.Size(); ++i )
			perturbed.state( i ) *= 1.0 + 1.0e-3*std::sin( 1.7*i + 0.4 );
		perturbed.psiAxis *= 1.0 + 1.0e-4;
		perturbed.psiBoundary *= 1.0 - 1.0e-4;
		std::shared_ptr<meq::Profile const> const table = perturbed.crossing;
		perturbed.crossing = std::make_shared<meq::SplineProfile const>(
			[ table ]( double x ) { return ( *table )( x ) + 1.0e-4*std::sin( 7.0*x ); },
			[ table ]( double x ) { return table->prime( x ) + 7.0e-4*std::cos( 7.0*x ); }, 1024 );

		for ( auto const &start : { std::make_pair( "from the answer perturbed", static_cast<Start const *>( &perturbed ) ),
		                            std::make_pair( "from the vacuum's table held", static_cast<Start const *>( &fromHeld ) ) } )
		{
			// The vacuum start's target: a restart takes its cold reference at
			// its own labels, which moves the target enough to add a step.
			// meq::newtonTolerance() is the solver's own floor on it, ten eps.
			Options rb = b;
			rb.start = start.second;
			rb.relativeTolerance = 1.0e-30;
			rb.absoluteTolerance = meq::newtonTolerance( 1.0e-12 )*bordered->solver->newtonResiduals().front();
			std::unique_ptr<MirrorRun> const restarted = runMirror( rb );
			BOOST_TEST_REQUIRE( restarted->converged, restarted->failure );
			printMerits( ( std::string( "bordered, " ) + start.first ).c_str(), *restarted );
			double const lands = maxDifference( *restarted, *bordered );
			std::printf( "      lands on the bordered answer from the vacuum to %.3e\n", lands );
			BOOST_TEST( stepsOf( *restarted ) <= 2 );
			landsOnTheAnswer( lands, *bordered, start.first );
		}
	}
	std::fflush( stdout );
}

/*
 * 3. THREE SPECIES AT M = 6, at the midplane peak of 5e18 where the retired
 * frozen refresh took 7 linear steps: the bordered gauge converges in 3, and
 * from the plasma solved on the vacuum's crossing held in 2, to the same
 * answer.
 *
 * WHAT IT DOES NOT CARRY HIGHER: at constant tau = 1, 2e19 and 5e19 stalled on
 * this level from every start, as the frozen refresh did and as the held table
 * does -- so the left-out term was not what stopped them; the stabilisation's
 * element-local fold was ( [M-284](MEASUREMENTS.md#m-284),
 * [M-291](MEASUREMENTS.md#m-291) ). Under the mirror's stabilisation the
 * bordered gauge converges from the vacuum at both on this level, which the
 * vacuum-gauge suite's case 5 asserts; the held table at 5e19 is past the
 * vacuum-midplane problem's own fold and has no answer to start from.
 */
BOOST_AUTO_TEST_CASE( threeSpeciesAtMachSixAreQuadraticUnderTheBorderedGauge )
{
	MirrorProfiles const three = threeSpecies( 5.0e18 );
	Options b;
	b.level = 0;
	b.order = 3;
	b.stem = borderStem;
	b.profiles = &three;
	b.gaugeKnots = borderKnots;
	std::unique_ptr<MirrorRun> const bordered = runMirror( b );
	BOOST_TEST_REQUIRE( bordered->converged, bordered->failure );
	Options h = b;
	h.holdCrossing = true;
	std::unique_ptr<MirrorRun> held = runMirror( h );
	BOOST_TEST_REQUIRE( held->converged, held->failure );
	Start fromHeld = startOf( *held );
	fromHeld.crossing = answerCrossing( *held, 1024 );
	Options c = b;
	c.start = &fromHeld;
	c.relativeTolerance = 1.0e-30;
	c.absoluteTolerance = meq::newtonTolerance( 1.0e-12 )*bordered->solver->newtonResiduals().front();
	std::unique_ptr<MirrorRun> const continued = runMirror( c );
	BOOST_TEST_REQUIRE( continued->converged, continued->failure );

	std::printf( "\n  THREE SPECIES AT M = 6, 5e18, level 0, k = 3\n" );
	printMerits( "bordered, from the vacuum", *bordered );
	printMerits( "the vacuum's table held", *held );
	printMerits( "bordered, from the held answer", *continued );
	double const lands = maxDifference( *continued, *bordered );
	std::printf( "    the two bordered answers %.3e apart\n", lands );
	std::fflush( stdout );
	BOOST_TEST( stepsOf( *bordered ) <= 3 );
	BOOST_TEST( stepsOf( *continued ) <= 2 );
	landsOnTheAnswer( lands, *bordered, "from the held answer" );
}

/*
 * 4. geq's MIRROR UNDER THE BORDERED GAUGE: M-278 section 4's matrix -- MEQ at
 * level 2, k = 4, on geq's midplane tables, against geq's direct solves -- with
 * MR-4d's three asks: falling at every refinement of geq, at geq's own rate, and
 * at 513 under geq's own 257 - 513 difference.
 */
BOOST_AUTO_TEST_CASE( meqAgreesWithGeqUnderTheBorderedGauge )
{
	std::vector<double> const g129 = readGeqPsi( "psi_direct_129.dat" );
	std::vector<double> const g257 = readGeqPsi( "psi_direct_257.dat" );
	std::vector<double> const g513 = readGeqPsi( "psi_direct_513.dat" );
	MirrorProfiles const geq = geqMidplaneProfiles();
	ladderStemInUse = borderStem;
	Options o;
	o.level = bestLevel;
	o.order = bestOrder;
	o.stem = borderStem;
	o.profiles = &geq;
	o.gaugeKnots = borderKnots;
	std::unique_ptr<MirrorRun> const run = runMirror( o );
	BOOST_TEST_REQUIRE( run->converged, run->failure );
	std::printf( "\n  MR-4d UNDER THE BORDERED GAUGE, level %d, k = %d, %d knots\n", bestLevel,
	             bestOrder, borderKnots );
	printMerits( "bordered, from the vacuum", *run );
	Sampled const best = reduce( run, o );

	std::vector<double> const meq = withoutCorners( best.psiStar );
	Difference const d129 = difference( meq, g129 );
	Difference const d257 = difference( meq, g257 );
	Difference const d513 = difference( meq, g513 );
	Difference const geqCoarse = difference( withoutCorners( g129 ), g257 );
	Difference const geqFine = difference( withoutCorners( g257 ), g513 );
	Difference const fixedRadius = difference( meq, withoutCorners( reference().psiStar ) );
	printDifference( "MEQ bordered - geq direct 129", d129 );
	printDifference( "MEQ bordered - geq direct 257", d257 );
	printDifference( "MEQ bordered - geq direct 513", d513 );
	printDifference( "geq 257 - geq 513", geqFine );
	printDifference( "MEQ bordered - MEQ fixed radius", fixedRadius );
	std::printf( "    rates, max: MEQ - geq %.3f then %.3f; geq's own %.3f\n",
	             std::log2( d129.max/d257.max ), std::log2( d257.max/d513.max ),
	             std::log2( geqCoarse.max/geqFine.max ) );
	std::printf( "    at 513: %.3f ( max ) and %.3f ( rms ) of geq's own 257 - 513\n",
	             d513.max/geqFine.max, d513.rms/geqFine.rms );
	std::fflush( stdout );

	double const fluxScale = std::abs( best.psiBoundary );
	if ( !meq::tests::aboveTheFormatsFloor( geqFine.max, fluxScale ) )
	{
		std::printf( "    *** SKIPPED: geq's own 257 - 513 difference is under this build's "
		             "floor; the conversion bound alone is asserted. M-272\n" );
		BOOST_TEST( d513.max < meq::tests::exactAgreementFloor( fluxScale ) );
		return;
	}
	BOOST_TEST( ( d129.max > d257.max && d257.max > d513.max ) );
	for ( auto const &step : { std::make_pair( d129.max, d257.max ), std::make_pair( d257.max, d513.max ) } )
		BOOST_TEST( std::log2( step.first/step.second ) > std::log2( geqCoarse.max/geqFine.max ) - 0.30 );
	BOOST_TEST( d513.max < geqFine.max );
	BOOST_TEST( d513.rms < geqFine.rms );
}

/*
 * 5. WHAT THE BORDERED GAUGE REFUSES: a knot count not-a-knot cannot take, a
 * prescribed plasma current beside the knots, whose row reads the crossings and
 * has no corner against them, and a start that is not this machine's.
 */
BOOST_AUTO_TEST_CASE( theBorderedGaugeRefusesWhatItCannotCarry )
{
	MirrorProfiles const geq = geqMidplaneProfiles();
	Options o;
	o.level = 0;
	o.order = 2;
	o.stem = borderStem;
	o.profiles = &geq;
	o.maxIterations = 1;
	o.relativeTolerance = 1.0e-12;
	o.absoluteTolerance = 1.0e-14;
	std::unique_ptr<MirrorRun> const run = runMirror( o );
	meq::GradShafranovSolver &solver = *run->solver;
	for ( int bad : { -1, 0, 1, 2, 3 } )
		BOOST_CHECK_THROW( solver.setMidplaneGaugeKnots( bad ), std::invalid_argument );
	BOOST_CHECK_THROW( solver.setInitialState( solver.state(), std::vector<double>( 3, 0.0 ) ),
	                   std::invalid_argument );

	// A prescribed current beside the knots, on a second solve of the same
	// solver.
	solver.setMidplaneGaugeKnots( borderKnots );
	solver.setPlasmaCurrent( 1.0e-3 );
	try
	{
		solver.solve();
		BOOST_TEST( false, "a prescribed current beside the bordered gauge's knots was not refused" );
	}
	catch ( std::logic_error const &refused )
	{
		BOOST_TEST( std::string( refused.what() ).find( "prescribed plasma current" ) != std::string::npos,
		            refused.what() );
	}

	// A start from another mesh.
	Start wrong = startOf( *run );
	wrong.state.SetSize( wrong.state.Size() + 1 );
	wrong.state = 0.0;
	Options misfit = o;
	misfit.start = &wrong;
	std::unique_ptr<MirrorRun> const misfitRun = runMirror( misfit );
	BOOST_TEST( !misfitRun->converged );
	BOOST_TEST( misfitRun->failure.find( "another mesh" ) != std::string::npos, misfitRun->failure );
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * THE MIDPLANE GAUGE'S START, AND THE OUTER LOOP. MIRROR-PLAN.md section 5;
 * [M-284](MEASUREMENTS.md#m-284) and
 * [M-285](MEASUREMENTS.md#m-285).
 *
 * THE GUESS, MEASURED AND NOT TAKEN: that the same densities solved first on
 * the vacuum's crossing HELD -- the vacuum-midplane gauge, cheap and quadratic
 * -- is a start close enough that the bordered Newton then needs only a very
 * few of its dearer steps. It saves one bordered step of three at most, none
 * at 2e19, and its own three steps cost more than the one saved: the bordered
 * Newton from the vacuum is already quadratic, and its first step does what
 * the whole held solve does. One ctest, RotatingMirrorRoute:
 *
 *   1. the bordered gauge from the vacuum -- what ships -- against the held
 *      start then the bordered Newton, and the held start stopped where its
 *      merit meets the bordered problem's; and how far each start is from
 *      the answer, by pressure;
 *   2. the OUTER LOOP: whole held-table solves, the table each time the last
 *      answer's crossings -- Picard on the curve, with no knots. Linear, at
 *      the retired frozen refresh's own contraction, which grows with the
 *      pressure.
 */
namespace
{
	char const *const routeStem = "RotatingMirrorRoute";

	int dampedSteps( MirrorRun const &run )
	{
		int count = 0;
		for ( meq::GradShafranovSolver::BorderStep const &s : run.solver->borderSteps() )
			if ( s.damping > 0.0 && s.damping < 1.0 )
				++count;
		return count;
	}

	void printRun( char const *label, MirrorRun const &run, double seconds )
	{
		std::vector<double> const &h = run.solver->newtonResiduals();
		std::printf( "    %-36s %s %2d steps (damped %d), order %5.2f, %.3f s:", label,
		             run.converged ? "converged" : "FAILED   ", stepsOf( run ), dampedSteps( run ),
		             lastOrder( h ), seconds );
		for ( double v : h )
			std::printf( " %.2e", v );
		std::printf( "\n" );
		std::fflush( stdout );
	}

	double median( std::vector<double> v )
	{
		std::sort( v.begin(), v.end() );
		return v.empty() ? std::numeric_limits<double>::quiet_NaN() : v[ v.size()/2 ];
	}

	/// The worst | R_0 | apart of a table and the bordered answer's knots.
	double crossingGap( meq::Profile const &a, std::vector<double> const &knots )
	{
		double worst = 0.0;
		for ( std::size_t k = 0; k < knots.size(); ++k )
			worst = std::max( worst, std::abs( a( static_cast<double>( k )/( knots.size() - 1 ) ) - knots[ k ] ) );
		return worst;
	}

	/// The worst | R_0 | apart of two tables, over 1025 points.
	double crossingGap( meq::Profile const &a, meq::Profile const &b )
	{
		double worst = 0.0;
		for ( int k = 0; k <= 1024; ++k )
			worst = std::max( worst, std::abs( a( k/1024.0 ) - b( k/1024.0 ) ) );
		return worst;
	}

	/// psi* + psi_c of the vacuum a run started from, at geq's nodes.
	std::vector<double> vacuumStar( MirrorRun &run )
	{
		BOOST_TEST_REQUIRE( static_cast<bool>( run.vacuum ) );
		run.vacuum->solver->postProcess();
		Nodes const &nodes = geqNodes();
		return totalAt( run, run.vacuum->solver->postProcessedPotential(), nodes.radius, nodes.z );
	}

	double largest( std::vector<double> const &v )
	{
		double out = 0.0;
		for ( double x : v )
			if ( std::isfinite( x ) )
				out = std::max( out, std::abs( x ) );
		return out;
	}

	/// The cases, at the coarsest level each converges on: level 0 for all
	/// three under the mirror's stabilisation, where constant tau = 1 took
	/// three species at 2e19 to level 1 (vacuum-gauge suite, case 5).
	struct RouteCase
	{
		std::string name;
		MirrorProfiles profiles;
		int level;
	};

	std::vector<RouteCase> routeCases()
	{
		return { { "geq's two species", geqMidplaneProfiles(), 0 },
		         { "three species at 5e18", threeSpecies( 5.0e18 ), 0 },
		         { "three species at 2e19", threeSpecies( 2.0e19 ), 0 } };
	}
}

BOOST_AUTO_TEST_SUITE( route )

/*
 * 1. THE BORDERED GAUGE FROM THE VACUUM, AND THE HELD START IT WAS GUESSED TO
 * NEED. Per case, k = 3, the default knots, each solve three times for a
 * median wall clock (printed, never asserted: the machine is shared):
 *
 *   (a) the bordered Newton from the vacuum, its knots on the vacuum's
 *       crossings -- what solve() does;
 *   (b) the held table from the vacuum, then the bordered Newton from that
 *       answer, its knots on that answer's crossings;
 *   (b') the same with the held solve stopped once its merit is under the
 *       merit the bordered problem has at its answer, which (b) measures --
 *       the derived hand-over, below which held steps buy nothing.
 *
 * ASSERTED: (a) is quadratic in three steps; every route lands on one answer;
 * and the held start saves at most one bordered step, so its own steps -- at
 * two thirds of a bordered step's cost -- cannot pay for it.
 */
BOOST_AUTO_TEST_CASE( theHeldStartSavesAtMostOneBorderedStep )
{
	int const repeats = 3;
	std::shared_ptr<meq::SplineProfile const> const vacuum = vacuumCrossing( 1024 );
	for ( RouteCase const &c : routeCases() )
	{
		std::printf( "\n  %s, level %d, k = 3\n", c.name.c_str(), c.level );
		Options b;
		b.level = c.level;
		b.order = 3;
		b.stem = routeStem;
		b.profiles = &c.profiles;
		Options h = b;
		h.holdCrossing = true;
		h.heldCrossing = vacuum;

		std::vector<double> coldTimes, heldTimes, afterTimes, earlyTimes, afterEarlyTimes;
		std::unique_ptr<MirrorRun> cold, held, after, early, afterEarly;
		Start fromHeld, fromEarly;
		for ( int take = 0; take < repeats; ++take )
		{
			cold = runMirror( b );
			BOOST_TEST_REQUIRE( cold->converged, cold->failure );
			coldTimes.push_back( cold->plasmaSeconds );
			// ONE STOPPING TARGET FOR EVERY BORDERED SOLVE: (a)'s, which its
			// cold reference sets. A restart takes its reference at its own
			// labels, which moves the target by a factor near one -- enough to
			// add a step on a quadratic tail.
			double const target = meq::newtonTolerance( 1.0e-12 )*cold->solver->newtonResiduals().front();

			held = runMirror( h );
			BOOST_TEST_REQUIRE( held->converged, held->failure );
			heldTimes.push_back( held->plasmaSeconds );
			fromHeld = startOf( *held );
			fromHeld.crossing = answerCrossing( *held, 1024 );
			Options a = b;
			a.start = &fromHeld;
			a.relativeTolerance = 1.0e-30;
			a.absoluteTolerance = target;
			after = runMirror( a );
			BOOST_TEST_REQUIRE( after->converged, after->failure );
			afterTimes.push_back( after->plasmaSeconds );

			// THE HAND-OVER: the held solve stopped at the merit the bordered
			// problem reads at the held answer.
			Options he = h;
			he.relativeTolerance = 1.0e-30;
			he.absoluteTolerance = after->solver->newtonResiduals().front();
			early = runMirror( he );
			BOOST_TEST_REQUIRE( early->converged, early->failure );
			earlyTimes.push_back( early->plasmaSeconds );
			fromEarly = startOf( *early );
			fromEarly.crossing = answerCrossing( *early, 1024 );
			Options ae = a;
			ae.start = &fromEarly;
			afterEarly = runMirror( ae );
			BOOST_TEST_REQUIRE( afterEarly->converged, afterEarly->failure );
			afterEarlyTimes.push_back( afterEarly->plasmaSeconds );
		}

		double const tCold = median( coldTimes );
		double const tHeld = median( heldTimes ), tAfter = median( afterTimes );
		double const tEarly = median( earlyTimes ), tAfterEarly = median( afterEarlyTimes );
		printRun( "(a) bordered, from the vacuum", *cold, tCold );
		printRun( "(b) held, from the vacuum", *held, tHeld );
		printRun( "    then bordered", *after, tAfter );
		printRun( "(b') held to the hand-over", *early, tEarly );
		printRun( "    then bordered", *afterEarly, tAfterEarly );
		double const heldStep = tHeld/stepsOf( *held );
		double const borderedStep = tCold/stepsOf( *cold );
		std::printf( "    total: (a) %.3f s, (b) %.3f s, (b') %.3f s; a held step %.3f s, a bordered step "
		             "%.3f s, %.2fx\n", tCold, tHeld + tAfter, tEarly + tAfterEarly, heldStep,
		             borderedStep, borderedStep/heldStep );

		// How far each start is from the answer.
		std::vector<double> const knots = cold->solver->midplaneGaugeCrossings();
		Sampled const answer = reduce( cold, b );
		double const scale = largest( answer.plasma );
		Difference const dVacuum = difference( vacuumStar( *cold ), answer.psiStar );
		Sampled const heldAnswer = reduce( held, h );
		Difference const dHeld = difference( heldAnswer.psiStar, answer.psiStar );
		Sampled const earlyAnswer = reduce( early, h );
		Difference const dEarly = difference( earlyAnswer.psiStar, answer.psiStar );
		std::printf( "    from the answer: the vacuum psi* %.3e ( %.3f of the plasma's ), R_0 %.3e m\n",
		             dVacuum.max, dVacuum.max/scale, crossingGap( *vacuum, knots ) );
		std::printf( "                     the held answer %.3e ( %.3f ), R_0 %.3e m\n", dHeld.max,
		             dHeld.max/scale, crossingGap( *fromHeld.crossing, knots ) );
		std::printf( "                     the hand-over %.3e ( %.3f ), R_0 %.3e m\n", dEarly.max,
		             dEarly.max/scale, crossingGap( *fromEarly.crossing, knots ) );
		double const lands = std::max( maxDifference( *after, *cold ), maxDifference( *afterEarly, *cold ) );
		std::printf( "    the three answers within %.3e\n", lands );
		std::fflush( stdout );

		std::string why;
		BOOST_TEST( stepsOf( *cold ) <= 3, c.name );
		if ( tailIsReadable( cold->solver->newtonResiduals() ) )
			BOOST_TEST( quadraticTail( cold->solver->newtonResiduals(), why ), c.name << ": " << why );
		landsOnTheAnswer( lands, *cold, c.name );
		BOOST_TEST( stepsOf( *after ) >= stepsOf( *cold ) - 1,
			c.name << ": from the held answer the bordered Newton took " << stepsOf( *after )
			<< " steps against " << stepsOf( *cold ) << " from the vacuum, so the held start may now "
			"pay for itself -- measure the totals again; M-284" );
		BOOST_TEST( stepsOf( *afterEarly ) >= stepsOf( *cold ) - 1, c.name );
	}
}

/*
 * 2. THE OUTER LOOP ON THE CURVE: the held table solved whole, its table the
 * vacuum's crossing first and each time after the last answer's, warm-started
 * from that answer -- the frozen refresh with every inner solve converged. Its
 * fixed point is the midplane gauge's answer. Per case, how fast the table
 * stops moving, and what it costs against the bordered gauge on the same
 * case. The contraction is the frozen refresh's: about 0.016 for geq's two
 * species, 0.04 for three at 5e18, 0.18 at 2e19 -- the left-out pull of the
 * field on R_0 grows with the pressure.
 *
 * ASSERTED: it converges, it is LINEAR -- the last contraction within a factor
 * of two of the first -- and it lands on the bordered answer to that answer's
 * own knot interpolation: 1.5e-05 of the plasma's flux at most, for geq's two
 * species at level 0, where that level's discretisation error is 1.7e-02 of
 * it (M-278 section 3).
 */
BOOST_AUTO_TEST_CASE( theOuterLoopOnTheCurveIsLinear )
{
	std::shared_ptr<meq::SplineProfile const> const vacuum = vacuumCrossing( 1024 );
	for ( RouteCase const &c : routeCases() )
	{
		std::printf( "\n  OUTER LOOP, %s, level %d, k = 3\n", c.name.c_str(), c.level );
		Options b;
		b.level = c.level;
		b.order = 3;
		b.stem = routeStem;
		b.profiles = &c.profiles;
		std::unique_ptr<MirrorRun> bordered = runMirror( b );
		BOOST_TEST_REQUIRE( bordered->converged, bordered->failure );
		Sampled const answer = reduce( bordered, b );
		printRun( "bordered, from the vacuum", *bordered, bordered->plasmaSeconds );

		Options o = b;
		o.holdCrossing = true;
		std::shared_ptr<meq::Profile const> curve = vacuum;
		std::unique_ptr<Start> previous;
		std::vector<double> contractions;
		double lastGap = std::numeric_limits<double>::quiet_NaN();
		int inner = 0;
		double seconds = 0.0;
		bool settled = false;
		std::unique_ptr<MirrorRun> run;
		for ( int outer = 0; outer < 25 && !settled; ++outer )
		{
			Options step = o;
			step.heldCrossing = curve;
			step.start = previous.get();
			run = runMirror( step );
			BOOST_TEST_REQUIRE( run->converged, "outer step " << outer << ": " << run->failure );
			inner += stepsOf( *run );
			seconds += run->plasmaSeconds;
			std::shared_ptr<meq::Profile const> const next = answerCrossing( *run, 1024 );
			double const gap = crossingGap( *next, *curve );
			if ( std::isfinite( lastGap ) )
				contractions.push_back( gap/lastGap );
			std::printf( "    outer %2d: %d inner steps, the table moved %.3e m, contraction %.3e\n", outer,
			             stepsOf( *run ), gap, gap/lastGap );
			lastGap = gap;
			curve = next;
			previous = std::make_unique<Start>( startOf( *run ) );
			settled = gap < 1.0e-12;
		}
		Sampled const outerAnswer = reduce( run, o );
		Difference const gap = difference( outerAnswer.psiStar, answer.psiStar );
		std::printf( "    total: %d inner steps, %.2f s, against the bordered gauge's %d steps, %.2f s; "
		             "the two answers psi* %.3e ( %.2e of the plasma's )\n", inner, seconds,
		             stepsOf( *bordered ), bordered->plasmaSeconds, gap.max, gap.max/largest( answer.plasma ) );
		std::fflush( stdout );

		BOOST_TEST( settled, c.name << ": the outer loop's table was still moving by " << lastGap << " m" );
		BOOST_TEST_REQUIRE( contractions.size() >= 3u );
		// The last ratio above the floor, against the first.
		double const first = contractions.front();
		double const last = contractions[ contractions.size() - 2 ];
		BOOST_TEST( ( last > 0.5*first && last < 2.0*first ),
			c.name << ": the contraction went from " << first << " to " << last << ", so the loop is not linear" );
		BOOST_TEST( inner > stepsOf( *bordered ) );
		// IN SINGLE THE LOOP CAN READ SETTLED SHORT OF ITS FIXED POINT. Once
		// the table's pull on the field is under the inner solve's target, ten
		// eps, that solve takes no step, its crossings are its start's and the
		// table stops moving. Measured: at 2e19 the last moving step was
		// 3.2e-06 m at a contraction of 0.18, which leaves the answer 2.4e-04
		// of the plasma's flux from the bordered one. Where it stopped so, the
		// landing is the format's and is skipped loudly; double never takes a
		// zero step here and asserts it as it stands. M-272.
		bool stoppedOnTheTarget = false;
#ifdef MFEM_USE_SINGLE
		stoppedOnTheTarget = stepsOf( *run ) == 0;
#endif
		if ( stoppedOnTheTarget )
			std::printf( "    *** SKIPPED the landing: the last inner solve took no step, so the "
			             "loop stopped at this build's target rather than its fixed point. M-272\n" );
		else
			BOOST_TEST( gap.max < 1.0e-4*largest( answer.plasma ), c.name );
	}
}

BOOST_AUTO_TEST_SUITE_END()
