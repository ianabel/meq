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

		if ( o.coldStart )
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
				meq::RotatingSource::Closure::Automatic, profiles.gauge );
		run->source->setPlasmaSupportShape( meq::PlasmaSupportShape::Box, p.axialEdge );

		run->dtn = std::make_unique<meq::ExteriorDtN>( 0.0, rhoGamma, o.modes );
		run->solver = std::make_unique<meq::GradShafranovSolver>( *run->domain->sub,
		                                                          o.order );
		meq::GradShafranovSolver &solver = *run->solver;
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
		if ( o.relativeTolerance > 0.0 )
			solver.setNewtonControl( o.relativeTolerance, o.absoluteTolerance,
			                         o.maxIterations );

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
		double const relative = history[ n ]/history.front();
		if ( order > 1.5 || relative <= 2.0e-15 )
			return true;
		char buffer[ 256 ];
		std::snprintf( buffer, sizeof buffer, "the last triple's order is %.3f and "
		               "it ends at %.3e of the first residual, above round-off; ",
		               order, relative );
		why = buffer;
		return false;
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
		/// psi* + psi_c at the two label points, which psiAxis() and
		/// psiBoundary() read by psi_h.
		double starAxis = 0.0;
		double starBoundary = 0.0;
		/// The merit history, and how often the midplane gauge re-tabulated
		/// its crossings -- zero under the fixed radius.
		std::vector<double> history;
		int midplaneRefreshes = 0;
		int midplaneHolds = 0;  ///< iterations that kept the previous table
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
		s.midplaneHolds = run->solver->midplaneGaugeHolds();
		s.psiAxis = run->solver->psiAxis();
		s.psiBoundary = run->solver->psiBoundary();
		Nodes const &nodes = geqNodes();
		s.psi = totalAt( *run, run->solver->potential(), nodes.radius, nodes.z );
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
 * THREE SUITES, ONE ctest EACH, so that the solve's acceptance does not wait on
 * the resolution study: `solve` is MR-4b and its controls, `resolution` is
 * MR-4c and MR-4d, which share their runs through sampled(), and `driver` is
 * MR-5 -- the shipped example through the shipped binary, against the library.
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

	std::printf( "\n    %2s %5s %7s %6s | %10s %10s %6s | %10s %10s %6s | %10s %10s\n",
	             "k", "level", "NE", "s", "psi_h max", "rms", "rate", "psi* max", "rms",
	             "rate", "psi*( ax )", "psi*( bnd )" );
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<double> potential, post;
		std::vector<int> const levels = rungs( order );
		for ( int level : levels )
		{
			Sampled const &s = sampled( level, order, true );
			Difference const h = difference( s.psi, ref.psiStar );
			Difference const star = difference( s.psiStar, ref.psiStar );
			potential.push_back( h.rms );
			post.push_back( star.rms );
			double hRate = std::numeric_limits<double>::quiet_NaN();
			double starRate = std::numeric_limits<double>::quiet_NaN();
			std::size_t const n = potential.size();
			if ( n > 1 )
			{
				hRate = meq::tests::rate( potential[ n - 2 ], potential[ n - 1 ], 2.0 );
				starRate = meq::tests::rate( post[ n - 2 ], post[ n - 1 ], 2.0 );
			}
			std::printf( "    %2d %5d %7d %6.1f | %10.3e %10.3e %6.3f | %10.3e %10.3e "
			             "%6.3f | %10.2e %10.2e\n", order, level, s.elements, s.seconds,
			             h.max, h.rms, hRate, star.max, star.rms, starRate,
			             s.starAxis - ref.starAxis, s.starBoundary - ref.starBoundary );
			std::fflush( stdout );
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
	// In single both columns end under the format's floor -- aligned 5.5e-08
	// against 16 eps x psi_bnd = 1.25e-07 -- so the factor between them is
	// round-off. A double-build acceptance there, skipped loudly. M-272.
	if ( !meq::tests::aboveTheFormatsFloor( aligned.back(), ref.psiBoundary ) )
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
	/// under @a prefix, @a solverLines added under [solver], and each line
	/// starting with a key in @a lines replaced whole. TEXT SUBSTITUTION ON
	/// THE SHIPPED FILE, because a copy would drift and a drifted copy passing
	/// is worse than no test; shippedConfig() makes the tables it reads
	/// absolute, since the variant is written away from them.
	std::string configureMirror(
		std::string const &prefix, std::string const &solverLines,
		std::vector<std::pair<std::string, std::string>> const &lines = {} )
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
		for ( auto const &[ key, line ] : lines )
			replaceLine( key, line );
		if ( !solverLines.empty() )
		{
			std::size_t const at = text.find( "\n[solver]\n" );
			BOOST_TEST_REQUIRE( at != std::string::npos );
			text.insert( at + 10, solverLines );
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
			                        + " " + config, out.stem + ".log" );
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
 *   AND THE FILE'S PsiAxis IS THE FILAMENTS' OWN FLUX AT THE INNER POINT, so
 *   the axis row reads zero at the first step.
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
	BOOST_TEST( std::abs( source.psiAxisGuess() - filamentsAtAxis )
	            <= 1.0e-12*std::abs( filamentsAtAxis ),
		"[source] PsiAxis is " << source.psiAxisGuess() << " where the filaments' "
		"flux at the inner point is " << filamentsAtAxis << "; print it with "
		"meq::ConductorField over geqFilaments() at %.17g" );
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
	for ( char const *phrase : { "converged in ", "at the prescribed point",
	                             "psi_bnd starts at", "axial cutoff",
	                             "the plasma support is the BOX" } )
		BOOST_TEST( says( run.log, phrase ), "the driver did not say '" << phrase << "'" );
	for ( char const *phrase : { "no interior extremum", "wall-hugging",
	                             "NOT a magnetic axis", "NO interior extremum",
	                             "did not converge", "Retrying",
	                             "support is EMPTY" } )
		BOOST_TEST( !says( run.log, phrase ),
			"the driver said '" << phrase << "' of a mirror:\n" << run.log );
	int const driverSteps = static_cast<int>( numberAfter( run.log, "converged in " ) );
	double const driverAxis = numberAfter( run.log, "     psi_ax = " );

	// ---- the library, on geq's problem and the driver's mesh ---------------
	Options o;
	o.meshPath = mirrorMesh;
	o.level = 0;
	o.order = 2;
	o.modes = defaultModes;
	o.coldStart = true;
	o.neutraliseDeuterium = true;
	o.axisGuess = source.psiAxisGuess();
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
	std::printf( "    psi_ax printed %.6e, library %.12e\n", driverAxis, m.psiAxis );
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
		"the driver printed psi_ax = " << driverAxis << " against " << m.psiAxis );

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
 * A MIRROR WHOSE FIRST ATTEMPT FAILS IS NOT RETRIED: beside the box the
 * bordered Picard-then-Newton rung is skipped and the failure exits 2 as it
 * is. MEASUREMENTS.md M-RETRY.
 *
 * From [source] PsiAxis = 0.2 the bordered Newton fails -- in 30 steps, and
 * 300 do no better -- and the rung, run, converges in eight to a SECOND
 * DISCRETE ROOT: well converged, with one element of the box at -7.9e-02
 * where the shipped answer's largest |psi_h| is 6.6e-03. So exit 2 is the
 * assertion that matters: putting the rung back turns it into an exit 0
 * reporting that root.
 *
 * Of the rung's two triggers this is the failure. The annulus cannot fire on
 * a box run, since [boundary.xpoint] is refused beside it, and
 * theShippedMirrorSolvesThroughTheDriver asserts that on a converged one.
 */
BOOST_AUTO_TEST_CASE( aMirrorWhoseFirstAttemptFailsIsNotRetried )
{
	std::string const prefix = "mirror-driver-no-retry";
	std::string const config = configureMirror( prefix, "",
	                                            { { "PsiAxis = ", "PsiAxis = 0.2" } } );
	std::string const stem = std::string( driverScratch ) + "/" + prefix;
	std::remove( ( stem + ".nc" ).c_str() );
	int const status = runLogged( "MKL_NUM_THREADS=1 " + std::string( MEQ_RUN_PATH )
	                              + " " + config, stem + ".log" );
	std::string const text = readText( stem + ".log" );
	BOOST_TEST( status == 2,
		"meq-run exited " << status << " on a mirror whose bordered Newton cannot "
		"converge from PsiAxis = 0.2. An exit 0 is the retry rung run beside the "
		"box and reporting the root it reaches from there, which is not the "
		"shipped equilibrium:\n" << text );
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
	/// The two densities of the reference's two-species plasma, tabulated in
	/// the midplane gauge from @a run's own crossings, at the knots of the
	/// fixed-radius tables. Two species share one exponent, the closed form's
	///
	///     A = C( Psi )( R^2 - R_ref^2 )/2,
	///     C = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ),
	///
	/// so the density where a surface crosses z = 0 is n_s0 e^A at R = R_0, and
	/// its slope is the chain rule through C and R_0:
	///
	///     d/dPsi ( n e^A ) = e^A [ n' + n( C' ( R_0^2 - R_ref^2 )/2 + C R_0 R_0' ) ].
	///
	/// R_0 and R_0' are tabulateMidplaneCrossing()'s, on the solver's own
	/// midplane field at the answer -- the same reading the solver gives the
	/// source at every iterate.
	MirrorProfiles transferToMidplane( MirrorRun &run, MirrorProfiles const &fixed )
	{
		BOOST_TEST_REQUIRE( fixed.species.size() == 2u );
		BOOST_TEST_REQUIRE( ( fixed.gauge == meq::RotatingSource::Gauge::FixedRadius ) );
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

		auto const *electron = dynamic_cast<meq::SplineProfile const *>(
			fixed.species[ 0 ].density.get() );
		BOOST_TEST_REQUIRE( electron != nullptr );
		std::vector<meq::Knot> const &knots = electron->knots();
		BOOST_TEST_REQUIRE( electron->uniformKnots() );
		BOOST_TEST_REQUIRE( knots.front().psi == 0.0 );
		BOOST_TEST_REQUIRE( knots.back().psi == 1.0 );
		std::shared_ptr<meq::SplineProfile const> const crossing = meq::tabulateMidplaneCrossing(
			normalised, rMin, rMax, 0.0, 1.0, static_cast<unsigned int>( knots.size() - 1 ) );

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

			meq::Knot const &r = crossing->knots()[ i ];
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
		out.gauge = meq::RotatingSource::Gauge::Midplane;
		return out;
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
			ladderStemInUse = midplaneStem;
			Options o;
			o.level = bestLevel;
			o.order = bestOrder;
			o.stem = midplaneStem;
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
		std::printf( "    %-24s %2d steps, %2d crossings tabulated, %d held:", label,
		             s.newtonSteps, s.midplaneRefreshes, s.midplaneHolds );
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

	MirrorProfiles threeSpecies()
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
		out.gauge = meq::RotatingSource::Gauge::Midplane;
		out.species.resize( 3 );
		out.species[ 0 ].mass = electronMass;
		out.species[ 0 ].charge = -1.0;
		out.species[ 0 ].temperature = temperature( 1.28e3 );
		out.species[ 0 ].density = std::make_shared<SkewProfile const>( densityPeak, densityFloor*densityPeak );
		out.species[ 1 ].mass = deuteronMass;
		out.species[ 1 ].charge = 1.0;
		out.species[ 1 ].temperature = temperature( 1.53e3 );
		out.species[ 2 ].mass = carbonMass;
		out.species[ 2 ].charge = 6.0;
		out.species[ 2 ].temperature = temperature( 1.53e3 );
		out.species[ 2 ].density = std::make_shared<SkewProfile const>( 0.01*densityPeak, densityFloor*0.01*densityPeak );
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
		BOOST_TEST( midplane.midplaneRefreshes + midplane.midplaneHolds == midplane.newtonSteps + 1,
			"the crossings were tabulated " << midplane.midplaneRefreshes << " times and held "
			<< midplane.midplaneHolds << " over " << midplane.newtonSteps << " Newton steps: "
			"once before the first residual and once per iteration after it is the contract" );
		BOOST_TEST( fixed.midplaneRefreshes == 0 );

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
 *   5e18 m^-3. At 5e19 the frozen crossings' own dependence on the field is
 *   too strong for the iteration and it stalls near 8e-05, MEASURED: the
 *   coupling the freeze leaves out scales with the pressure;
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
	BOOST_TEST( d12.rms < 0.25*d01.rms,
		"psi* does not self-converge: level to level " << d01.rms << " then " << d12.rms );
	BOOST_TEST( d12.max < 0.25*d01.max,
		"psi* does not self-converge in max: " << d01.max << " then " << d12.max );

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
