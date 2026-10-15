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

	/// What the source is built from: the species, omega and R_ref.
	struct MirrorProfiles
	{
		std::vector<meq::Species> species;
		std::shared_ptr<meq::Profile const> omega;
		double referenceRadius = 0.0;
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

		MirrorProfiles profiles = mirrorProfiles();
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
				profiles.referenceRadius, axisStart, machine().mu0 );
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
	};

	/// Each ( level, k, aligned ) solved once per process, since MR-4d's best
	/// run is MR-4c's finest.
	Sampled const &sampled( int level, int order, bool aligned )
	{
		static std::map<std::tuple<int, int, bool>, Sampled> cache;
		auto const key = std::make_tuple( level, order, aligned );
		auto const found = cache.find( key );
		if ( found != cache.end() )
			return found->second;

		Options o;
		o.level = level;
		o.order = order;
		o.aligned = aligned;
		o.stem = ladderStem;
		std::unique_ptr<MirrorRun> const run = runMirror( o );
		BOOST_TEST_REQUIRE( run->converged, "k = " << order << ", level " << level
			<< ( aligned ? ", aligned" : ", cut" ) << ": " << run->failure );

		Sampled s;
		s.order = order;
		s.level = level;
		s.aligned = aligned;
		s.elements = run->domain->sub->GetNE();
		s.h = run->domain->h;
		s.newtonSteps = static_cast<int>( run->solver->newtonResiduals().size() ) - 1;
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
		return cache.emplace( key, std::move( s ) ).first->second;
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
 * FIVE SUITES, ONE ctest EACH, so that the solve's acceptance does not wait on
 * the resolution study: `solve` is MR-4b and its controls, `resolution` is
 * MR-4c and MR-4d, which share their runs through sampled(), and `driver` is
 * MR-5 -- the shipped example through the shipped binary, against the library.
 * `adaptive` judges the estimator on that example and `device` runs it on the
 * card; both need the binary, as `driver` does.
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
		                        + std::string( MEQ_RUN_PATH ) + " " + config
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
	                             "NOT a magnetic axis", "did not converge", "Retrying",
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
 * WHAT IT FOUND, TWICE, AND BOTH WERE THE SPLIT. MEASUREMENTS.md M-ADAPT.
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
			"modal sum and not the remainder's. MEASUREMENTS.md M-ADAPT" );
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
 * max |psi_h| and 5.7e-14 to 1.0e-13 of max |q_h| (M-DEVICE). So the bound
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
	double const hostAxis = numberAfter( host.log, "     psi_ax = " );
	double const cardAxis = numberAfter( card.log, "     psi_ax = " );

	std::printf( "\n  THE MIRROR ON THE CARD against --device cpu, %d elements, one host "
	             "thread each\n", onHost.mesh->GetNE() );
	std::printf( "    psi_h: %.3e over %d dofs, %.3e of max |psi_h| %.3e\n", worst,
	             onHost.psi->Size(), worst/largest, largest );
	std::printf( "    q_h:   %.3e over %d dofs, %.3e of max |q_h| %.3e\n", fluxWorst,
	             onHost.flux->Size(), fluxWorst/fluxLargest, fluxLargest );
	std::printf( "    Newton %d on the host, %d on the card; psi_ax printed %.6e and %.6e\n",
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
