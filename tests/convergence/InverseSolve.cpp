/*
 * IV-3 AND IV-4: COIL CURRENTS AS BORDER UNKNOWNS, ON ONE MACHINE, DRIVEN BY
 * ISOFLUX TARGETS AND THEN BY THE X-POINT. INVERSE-AND-RECONSTRUCTION-PLAN.md
 * section 9. The IV-4 case is the last in this file and has its own header.
 *
 * examples/machine-f-diiid.toml is the one machine that converges forward from
 * cold on every route (MEASUREMENTS.md M-103), so a failure here is the
 * inverse border's and not the machine's. The case:
 *
 *   1. solve the file forward, with its support sweeps, as apps/meq.cpp does;
 *   2. read isoflux targets OFF ITS OWN CONVERGED SEPARATRIX -- points where
 *      psi_h = psi_bnd, found by bisection along rays from the located axis;
 *   3. perturb some coil currents, make them unknowns, and solve again;
 *   4. assert the targets are met AND the currents come back.
 *
 * **TWO ASSERTIONS AND THE SECOND IS THE ONE WITH TEETH**, which is section
 * 1.4's point: a shape can be matched to plotting accuracy by a current set
 * nothing like the machine's, so a shape-only assertion is the test Lackner's
 * figures say will pass on a wrong answer.
 *
 * **ITS WEAKNESS, STATED**: the targets come from the answer, so the problem is
 * reachable by construction and this cannot detect a target the machine cannot
 * meet. IV-6 is the test that can. And it is posed so the currents ARE
 * determined -- four circuits against many more rows -- since recovery is only
 * a meaningful assertion where the rank says the targets fix the currents. The
 * rank is asserted for that reason rather than assumed.
 */
#define BOOST_TEST_MODULE InverseSolve
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <iterator>
#include <limits>
#include <sys/wait.h>

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/Config.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Inverse.hpp"
#include "meq/PointLocator.hpp"
#include "meq/Sampler.hpp"
#include "meq/SourceFactory.hpp"
#include "meq/WarmStart.hpp"

#include "RoundOffFloor.hpp"

#include "SourceTree.hpp"

namespace
{
	std::string const machineFile = meq::test::sourcePath( "examples/machine-f-diiid.toml" );

	/// Negative inside: the semicircle about the axis that is Gamma.
	mfem::PositionFunction semicircle( double radius, double centreZ )
	{
		return [ radius, centreZ ]( mfem::Vector const &x )
		{
			return std::hypot( x( 0 ), x( 1 ) - centreZ ) - radius;
		};
	}

	/// The machine, kept alive together; the solver aliases everything above
	/// it. The coil set is held WRITABLE, which is what an inverse solve needs
	/// and the forward fixtures in this tree do not.
	struct Machine
	{
		std::unique_ptr<meq::Configuration> config;
		std::unique_ptr<mfem::Mesh> background;
		std::unique_ptr<mfem::SubMesh> sub;
		std::unique_ptr<mfem::VertexConePath> path;
		mfem::Array<int> gammaHMarker;

		std::shared_ptr<meq::CoilSet> coils;
		std::shared_ptr<meq::NormalisedSource> plasma;
		std::shared_ptr<meq::CoilAugmentedNormalisedSource> source;
		std::unique_ptr<meq::ExteriorDtN> exterior;

		std::unique_ptr<mfem::Mesh> guessMesh;
		std::unique_ptr<mfem::GridFunction> guess;
		std::unique_ptr<mfem::GridFunction> carried;
		mfem::ConstantCoefficient zero{ 0.0 };

		std::unique_ptr<meq::GradShafranovSolver> solver;
	};

	/// A solver on the machine, set up as the driver sets one up, starting
	/// `psi_ax` and the X-point where it is told to. Every arm of a case that
	/// must start from the same place gets its own: setXPointBoundary() is
	/// one call per solver by design, so a solver cannot be put back.
	std::unique_ptr<meq::GradShafranovSolver> makeSolver( Machine &m,
	                                                      double psiAxisStart,
	                                                      double xPointR,
	                                                      double xPointZ )
	{
		meq::SourceConfig const &sourceConfig = m.config->getSource();
		auto solver = std::make_unique<meq::GradShafranovSolver>(
			*m.sub, m.config->getDiscretisation().polynomialDegree );
		solver->setNewtonControl(
			meq::newtonTolerance( m.config->getSolver().newtonRelativeTolerance ),
			meq::newtonTolerance( m.config->getSolver().newtonAbsoluteTolerance ),
			m.config->getSolver().newtonMaxIterations );
		solver->setExtension( *m.path, m.gammaHMarker );
		solver->setSource( *m.source, psiAxisStart );
		solver->setXPointBoundary( xPointR, xPointZ );
		solver->setPlasmaCurrent( sourceConfig.permeability()
		                          *sourceConfig.plasmaCurrent() );
		solver->setBoundaryData( m.zero );
		solver->setExteriorCoupling( *m.exterior );
		return solver;
	}

	/// The mesh is a build product of the file, so it is made the way the
	/// driver's users make it rather than assumed to be on disk.
	void ensureMesh( std::string const &file = machineFile )
	{
		std::string const command = std::string( MEQ_RUN_PATH )
		                            + " --mesh-only " + file
		                            + " > /dev/null 2>&1";
		BOOST_TEST_REQUIRE( std::system( command.c_str() ) == 0,
		                    "meq-run --mesh-only failed on " << file );
	}

	/// The assembly apps/meq.cpp performs for this file, in its order: cut D_h
	/// from the half disc at Gamma, wrap the plasma source around the coils and
	/// confine it ON THE WRAPPER, then the X-point border, the prescribed
	/// current and the exterior coupling, and the stored guess interpolated in.
	Machine buildMachine( std::string const &file = machineFile )
	{
		Machine m;
		m.config = std::make_unique<meq::Configuration>( file );
		meq::SourceConfig const &sourceConfig = m.config->getSource();

		m.background = std::make_unique<mfem::Mesh>(
			m.config->getMesh().file.c_str(), 1, 1 );
		double h = 0.0;
		for ( int e = 0; e < m.background->GetNE(); ++e )
			h = std::max( h, meq::elementDiameter( *m.background, e ) );

		meq::ExteriorConfig const &exteriorConfig = m.config->getBoundary().exterior;
		BOOST_TEST_REQUIRE( exteriorConfig.given );
		mfem::PositionFunction const levelSet =
			semicircle( exteriorConfig.radius, exteriorConfig.centreZ );

		std::vector<int> original( static_cast<std::size_t>( m.background->GetNE() ) );
		for ( int e = 0; e < m.background->GetNE(); ++e )
			original[ static_cast<std::size_t>( e ) ] = m.background->GetAttribute( e );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE( mfem::MarkLevelSetSubdomain(
			*m.background, levelSet, meq::markingOffset( 1.0, 1.5 ), marker, 1 ) > 0 );
		int const parentBoundaryMax = m.background->bdr_attributes.Size() > 0
			? m.background->bdr_attributes.Max() : 0;
		for ( int e = 0; e < m.background->GetNE(); ++e )
			m.background->SetAttribute( e, marker[ e ] ? 1 : 2 );
		m.background->SetAttributes();

		mfem::Array<int> domainAttribute( 1 );
		domainAttribute[ 0 ] = 1;
		m.sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( *m.background, domainAttribute ) );
		{
			mfem::Array<int> const &parent = m.sub->GetParentElementIDMap();
			for ( int e = 0; e < m.sub->GetNE(); ++e )
				m.sub->SetAttribute( e, original[ static_cast<std::size_t>( parent[ e ] ) ] );
			m.sub->SetAttributes();
			for ( int e = 0; e < m.background->GetNE(); ++e )
				m.background->SetAttribute( e, original[ static_cast<std::size_t>( e ) ] );
			m.background->SetAttributes();
		}
		int const gammaH = m.sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( gammaH > parentBoundaryMax );
		m.path = std::make_unique<mfem::VertexConePath>( *m.sub, gammaH, levelSet,
		                                                6.0*h );
		m.gammaHMarker.SetSize( gammaH );
		m.gammaHMarker = 0;
		m.gammaHMarker[ gammaH - 1 ] = 1;

		m.coils = meq::makeCoilSet( m.config->getCoils(),
		                            sourceConfig.permeability(), file );
		m.plasma = meq::makeNormalisedSource( sourceConfig, file );
		m.source = std::make_shared<meq::CoilAugmentedNormalisedSource>( m.plasma,
		                                                                 m.coils );
		// ON THE WRAPPER, AFTER IT EXISTS -- plan section 6.1 and M-82.
		m.source->setPlasmaSupport( sourceConfig.confinesToPlasma() );
		m.exterior = std::make_unique<meq::ExteriorDtN>( exteriorConfig.centreZ,
		                                                exteriorConfig.radius,
		                                                exteriorConfig.modes );

		meq::XPointConfig const &xpoint = m.config->getBoundary().xpoint;
		BOOST_TEST_REQUIRE( xpoint.given );

		m.solver = makeSolver( m, sourceConfig.psiAxisGuess(), xpoint.radius,
		                       xpoint.z );

		meq::InitialGuessConfig const &guessConfig = m.config->getInitialGuess();
		m.guessMesh = std::make_unique<mfem::Mesh>( guessConfig.meshFile.c_str(), 1, 1 );
		std::ifstream stream( guessConfig.file );
		BOOST_TEST_REQUIRE( stream.good() );
		m.guess = std::make_unique<mfem::GridFunction>( m.guessMesh.get(), stream );
		meq::FieldTransfer transfer( *m.guessMesh );
		m.carried = std::make_unique<mfem::GridFunction>(
			m.solver->potential().FESpace() );
		transfer.transfer( *m.guess, m.zero, *m.carried );
		m.solver->setInitialGuess( *m.carried );
		return m;
	}

	/// psi_h at a batch of points, `byNODES`. Every point must be in the mesh.
	std::vector<double> psiAt( meq::GradShafranovSolver const &solver,
	                           std::vector<double> const &radius,
	                           std::vector<double> const &z )
	{
		mfem::GridFunction const &potential = solver.potential();
		meq::PointLocator locator( *potential.FESpace()->GetMesh() );
		int const count = static_cast<int>( radius.size() );
		mfem::Vector positions( 2*count );
		for ( int i = 0; i < count; ++i )
		{
			positions( i ) = radius[ static_cast<std::size_t>( i ) ];
			positions( count + i ) = z[ static_cast<std::size_t>( i ) ];
		}
		locator.locate( positions );
		BOOST_TEST_REQUIRE( locator.missed() == 0 );
		mfem::Vector values;
		locator.interpolate( potential, values );
		return std::vector<double>( values.GetData(), values.GetData() + count );
	}

	/// Where psi_h = level along each ray from the axis, to round-off: the
	/// first downward crossing, bracketed by marching and closed by bisection.
	void contourPoints( meq::GradShafranovSolver const &solver, double level,
	                    std::vector<double> const &angles,
	                    std::vector<double> &radius, std::vector<double> &z )
	{
		double const r0 = solver.axisR();
		double const z0 = solver.axisZ();
		std::size_t const count = angles.size();
		std::vector<double> inside( count, 0.0 ), outside( count, -1.0 );

		for ( double reach = 0.02; reach < 1.5; reach += 0.02 )
		{
			std::vector<double> rr, zz;
			for ( double angle : angles )
			{
				rr.push_back( r0 + reach*std::cos( angle ) );
				zz.push_back( z0 + reach*std::sin( angle ) );
			}
			std::vector<double> const value = psiAt( solver, rr, zz );
			bool all = true;
			for ( std::size_t i = 0; i < count; ++i )
			{
				if ( outside[ i ] < 0.0 )
				{
					if ( value[ i ] < level )
						outside[ i ] = reach;
					else
						inside[ i ] = reach;
				}
				all = all && outside[ i ] > 0.0;
			}
			if ( all )
				break;
		}
		for ( std::size_t i = 0; i < count; ++i )
			BOOST_TEST_REQUIRE( outside[ i ] > 0.0, "ray " << i << " never left the plasma" );

		for ( int halving = 0; halving < 60; ++halving )
		{
			std::vector<double> rr, zz, middle( count );
			for ( std::size_t i = 0; i < count; ++i )
			{
				middle[ i ] = 0.5*( inside[ i ] + outside[ i ] );
				rr.push_back( r0 + middle[ i ]*std::cos( angles[ i ] ) );
				zz.push_back( z0 + middle[ i ]*std::sin( angles[ i ] ) );
			}
			std::vector<double> const value = psiAt( solver, rr, zz );
			for ( std::size_t i = 0; i < count; ++i )
				( value[ i ] < level ? outside[ i ] : inside[ i ] ) = middle[ i ];
		}
		radius.clear();
		z.clear();
		for ( std::size_t i = 0; i < count; ++i )
		{
			double const reach = 0.5*( inside[ i ] + outside[ i ] );
			radius.push_back( r0 + reach*std::cos( angles[ i ] ) );
			z.push_back( z0 + reach*std::sin( angles[ i ] ) );
		}
	}

	/// The driver's support loop: freeze the edge, refresh the component,
	/// solve, repeat until the support repeats. Returns whether it settled.
	bool solveWithSweeps( Machine &m, int &sweepsRun )
	{
		meq::GradShafranovSolver &solver = *m.solver;
		int const sweeps = m.config->getSolver().plasmaSupportSweeps;
		solver.prepare();
		mfem::GridFunction state( solver.potential() );
		double axis = m.config->getSource().psiAxisGuess();
		meq::XPointConfig const &xpoint = m.config->getBoundary().xpoint;
		double boundary = psiAt( solver, { xpoint.radius }, { xpoint.z } )[ 0 ];

		sweepsRun = 0;
		for ( int sweep = 0; sweep < sweeps; ++sweep )
		{
			int const before = solver.plasmaComponentElements();
			m.source->freezePlasmaEdge( axis, boundary );
			solver.setPlasmaSupportFrozen( false );
			solver.refreshPlasmaComponent( state );
			solver.setPlasmaSupportFrozen( true );
			if ( sweep > 0 )
				solver.setInitialGuess( state );
			solver.solve();
			++sweepsRun;
			state = solver.potential();
			axis = solver.psiAxis();
			boundary = solver.psiBoundary();
			if ( sweep > 0 && solver.plasmaComponentElements() == before )
				return true;
		}
		return false;
	}

	/// Points on the converged separatrix, rays at `count` equal angles
	/// about the axis, dropping any ray within 0.15 rad of the X-point: there
	/// the separatrix has a corner and the crossing is a saddle rather than a
	/// transversal one.
	void separatrixPoints( meq::GradShafranovSolver const &solver, int count,
	                       std::vector<double> &radius, std::vector<double> &z )
	{
		double const towardX = std::atan2( solver.xPointZ() - solver.axisZ(),
		                                   solver.xPointR() - solver.axisR() );
		std::vector<double> angles;
		for ( int i = 0; i < count; ++i )
		{
			double const angle = 2.0*M_PI*( i + 0.5 )/count;
			if ( std::abs( std::remainder( angle - towardX, 2.0*M_PI ) ) > 0.15 )
				angles.push_back( angle );
		}
		contourPoints( solver, solver.psiBoundary(), angles, radius, z );

		std::vector<double> const check = psiAt( solver, radius, z );
		double worst = 0.0;
		for ( double v : check )
			worst = std::max( worst, std::abs( v - solver.psiBoundary() ) );
		std::printf( "targets: %zu points on psi = psi_bnd, worst |psi - psi_bnd| %.3e\n",
		             radius.size(), worst );
		BOOST_TEST_REQUIRE( worst < meq::tests::scalarAgreementFloor(
		                                1.0e-12, solver.psiBoundary() ) );
	}

	/// What one inverse solve reports, printed and returned.
	struct InverseOutcome
	{
		bool solved = false;
		double worstTarget = 0.0;
		double worstCurrent = 0.0;
		double condition = 0.0;
	};

	/// Perturb the named coils by `factor`, free them one circuit each, ask
	/// for every separatrix point to share the first one's flux, and solve
	/// from the forward answer on its support -- which the ORIGINAL currents
	/// satisfy exactly, so the fixed point to come back to exists.
	InverseOutcome perturbFreeAndSolve( Machine &m,
	                                    std::vector<std::size_t> const &indices,
	                                    std::vector<double> const &factor,
	                                    std::vector<double> const &targetR,
	                                    std::vector<double> const &targetZ,
	                                    std::vector<double> &original )
	{
		meq::GradShafranovSolver &solver = *m.solver;
		std::vector<meq::GradShafranovSolver::FreeCircuit> circuits;
		original.clear();
		for ( std::size_t k = 0; k < indices.size(); ++k )
		{
			std::size_t const index = indices[ k ];
			double const current = m.coils->coil( index ).current();
			original.push_back( current );
			m.coils->setCurrent( index, factor[ k ]*current );
			circuits.push_back( { m.config->getCoils().coils[ index ].name,
			                      { index }, { 1.0 } } );
		}
		solver.clearInverse();
		solver.setFreeCurrents( m.coils, circuits );
		for ( std::size_t i = 1; i < targetR.size(); ++i )
			solver.addIsofluxTarget( targetR[ 0 ], targetZ[ 0 ],
			                         targetR[ i ], targetZ[ i ] );

		mfem::GridFunction const forward( solver.potential() );
		solver.setInitialGuess( forward );
		InverseOutcome out;
		out.solved = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const &failure )
		{
			out.solved = false;
			std::printf( "inverse solve threw: %s\n", failure.what() );
		}

		std::vector<double> const sigma = solver.inverseSingularValues();
		if ( !sigma.empty() && sigma.back() > 0.0 )
			out.condition = sigma.front()/sigma.back();
		std::printf( "inverse: %zu circuits, %zu rows, %d Newton iterations, "
		             "rank %d nullity %d, sigma %.3e .. %.3e\n",
		             indices.size(), targetR.size() - 1, solver.newtonIterations(),
		             solver.inverseRank(), solver.inverseNullity(),
		             sigma.empty() ? 0.0 : sigma.front(),
		             sigma.empty() ? 0.0 : sigma.back() );
		std::vector<double> const &history = solver.newtonResiduals();
		for ( std::size_t i = 0; i < history.size(); ++i )
			std::printf( "   it %2zu  %.6e\n", i, history[ i ] );

		for ( double v : solver.targetResiduals() )
			out.worstTarget = std::max( out.worstTarget, std::abs( v ) );
		std::vector<double> const recovered = solver.freeCircuitCurrents();
		for ( std::size_t k = 0; k < indices.size() && k < recovered.size(); ++k )
		{
			double const relative = std::abs( recovered[ k ] - original[ k ] )
			                        /std::abs( original[ k ] );
			out.worstCurrent = std::max( out.worstCurrent, relative );
			std::printf( "   %-4s  original %+.10e  perturbed %+.10e  recovered "
			             "%+.10e  rel %.3e\n", circuits[ k ].name.c_str(),
			             original[ k ], factor[ k ]*original[ k ], recovered[ k ],
			             relative );
			// THE SET THE SOURCE READS HOLDS THE ANSWER; there is no second copy.
			BOOST_TEST( m.coils->coil( indices[ k ] ).current() == recovered[ k ] );
		}
		std::printf( "worst target %.3e  worst current %.3e\n", out.worstTarget,
		             out.worstCurrent );
		return out;
	}

	/**
	 * Whether this build's format can resolve currents from targets at
	 * all, given the corner's condition number.
	 *
	 * **THE SPREAD OF THE SINGULAR VALUES MULTIPLIES A TARGET'S ROUND-OFF INTO
	 * A CURRENT**, so a recovered current is good to about `condition x eps`
	 * relative. Every coil free against 21 rows reads a condition of 7.96e+07:
	 * that is 1.8e-08 in double, where the recovery is a claim about the
	 * method, and 9.5 in single, where the currents are undetermined and F3B
	 * -- 482 A beside 1e+05 -- comes back with its sign flipped. The rank is
	 * full in both; what the format decides is only how WELL. The factor of
	 * 1000 leaves three digits for the method to be judged on.
	 */
	bool recoveryIsMeasurable( double condition )
	{
		return 1000.0*condition*static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() ) < 1.0;
	}

	/// The coil index of a named [[coils]] block.
	std::size_t coilNamed( Machine const &m, std::string const &name )
	{
		std::vector<meq::CoilParameters> const &list = m.config->getCoils().coils;
		for ( std::size_t i = 0; i < list.size(); ++i )
			if ( list[ i ].name == name )
				return i;
		BOOST_FAIL( "no coil named " << name );
		return 0;
	}
}

BOOST_AUTO_TEST_CASE( theInverseSolveRecoversAKnownCurrentSet )
{
	ensureMesh();
	Machine m = buildMachine();
	meq::GradShafranovSolver &solver = *m.solver;

	// 1. FORWARD, as the driver runs the file.
	int sweeps = 0;
	BOOST_TEST_REQUIRE( solveWithSweeps( m, sweeps ),
	                    "the forward support did not settle in " << sweeps
	                    << " sweeps" );
	BOOST_TEST_REQUIRE( solver.axisWasLocated() );
	double const forwardAxis = solver.psiAxis();
	std::printf( "forward: %d sweeps, psi_ax %.10e  psi_bnd %.10e  X ( %.6f, %.6f )\n",
	             sweeps, forwardAxis, solver.psiBoundary(), solver.xPointR(),
	             solver.xPointZ() );

	// 2. THE TARGETS, off the converged separatrix.
	std::vector<double> targetR, targetZ;
	separatrixPoints( solver, 24, targetR, targetZ );

	// 3. FOUR CIRCUITS, inboard and outboard, above and below, so the
	// perturbation moves the shape and not only its height.
	std::vector<std::size_t> const four = { coilNamed( m, "F2A" ),
	                                        coilNamed( m, "F2B" ),
	                                        coilNamed( m, "F6A" ),
	                                        coilNamed( m, "F6B" ) };
	std::vector<double> original;
	InverseOutcome const got = perturbFreeAndSolve( m, four,
		{ 1.03, 0.98, 1.02, 0.97 }, targetR, targetZ, original );
	BOOST_TEST_REQUIRE( got.solved );

	/*
	 * THE RANK SAYS THE TARGETS DETERMINE THE CURRENTS, without which the
	 * recovery below would be asserting a choice the minimum norm made.
	 *
	 * THE BOUNDS ARE THE MEASUREMENT WITH ROOM, NOT A GUESS. Measured: the
	 * targets close to 1.3e-13 and the currents return to 1.4e-11 relative from
	 * a two-to-three per cent perturbation, in two Newton steps, with psi_ax
	 * back on the forward answer to its tenth digit. Each bound is about a
	 * hundred times its reading.
	 *
	 * **THE NEWTON COUNT IS THE ASSERTION THAT SEES THE COLUMN.** The residual
	 * is exact whatever the coil column says, so a wrong column still reaches
	 * the answer -- linearly, at a rate set by how wrong it is. Two steps from
	 * 2.9e-02 to 4.5e-13 is the exact column; mutated by 10% it took seven,
	 * and with the currents not written at the trial points, five.
	 */
	BOOST_TEST( solver.inverseNullity() == 0 );
	// Targets are flux, so their scale is psi_ax's; currents are relative.
	BOOST_TEST( got.worstTarget
	            < meq::tests::scalarAgreementFloor( 1.0e-11, forwardAxis ) );
	BOOST_TEST( got.worstCurrent < meq::tests::scalarAgreementFloor( 1.0e-9, 1.0 ) );
	BOOST_TEST( std::abs( solver.psiAxis() - forwardAxis )
	            < meq::tests::scalarAgreementFloor( 1.0e-9*std::abs( forwardAxis ),
	                                                forwardAxis ) );
	BOOST_TEST( solver.newtonIterations() <= 3 );

	/*
	 * 4. AND EVERY COIL, on the same machine and the same 21 rows -- now at the
	 * original currents again, so the same separatrix. This is the case plan
	 * section 1.4 warns about and it is worth knowing which way it comes out:
	 * MEASURED, the corner stays of FULL rank, 33 of 33, and all eighteen
	 * currents return to 1.2e-07 relative or better. The worst is F3B, which
	 * carries 482 A where the others carry 1e+05; the lower set is
	 * consistently the less well determined, the X-point being below.
	 *
	 * So a separatrix sampled at 22 points DETERMINES DIII-D's eighteen
	 * currents in exact arithmetic. What it does not say is how well: that is
	 * the singular values' spread, which is printed, and which is what
	 * multiplies a measured target's error into a current's. The bound here is
	 * looser than the four-circuit one for the same reason.
	 */
	std::vector<std::size_t> every;
	std::vector<double> alternate;
	for ( std::size_t index = 0; index < m.coils->size(); ++index )
	{
		every.push_back( index );
		alternate.push_back( index % 2 == 0 ? 1.02 : 0.98 );
	}
	InverseOutcome const all = perturbFreeAndSolve( m, every, alternate, targetR,
	                                                targetZ, original );
	std::printf( "every coil: condition %.3e\n", all.condition );
	BOOST_TEST_REQUIRE( all.solved );
	BOOST_TEST( solver.inverseNullity() == 0 );
	BOOST_TEST( all.worstTarget
	            < meq::tests::scalarAgreementFloor( 1.0e-11, forwardAxis ) );
	if ( recoveryIsMeasurable( all.condition ) )
		BOOST_TEST( all.worstCurrent < 1.0e-5 );
	else
		BOOST_TEST_MESSAGE( "every coil: condition " << all.condition
		                    << " times this format's eps leaves no digit of the "
		                       "currents to recover; the rank is asserted and the "
		                       "recovery is not" );
}

/*
 * AND WITH TOO FEW TARGETS THE SHAPE IS MET BY CURRENTS THAT ARE NOT THE
 * MACHINE'S, WHICH IS LACKNER'S WARNING AND WHAT THE NULLITY IS FOR.
 *
 * Every coil free against five rows: the corner is 20 x 33 and at most rank 20,
 * so at least thirteen directions in the currents are chosen by the minimum
 * norm rather than by the targets. The run must CONVERGE -- the targets are
 * reachable -- and meet them; it must REPORT the nullity; and the currents must
 * NOT come back, which is the point. A solver that recovered them here would
 * be asserting information the targets do not carry.
 */
BOOST_AUTO_TEST_CASE( tooFewTargetsMeetTheShapeWithOtherCurrents )
{
	ensureMesh();
	Machine m = buildMachine();
	meq::GradShafranovSolver &solver = *m.solver;
	int sweeps = 0;
	BOOST_TEST_REQUIRE( solveWithSweeps( m, sweeps ) );

	std::vector<double> targetR, targetZ;
	separatrixPoints( solver, 6, targetR, targetZ );
	BOOST_TEST_REQUIRE( targetR.size() >= 5u );

	std::vector<std::size_t> every;
	std::vector<double> alternate;
	for ( std::size_t index = 0; index < m.coils->size(); ++index )
	{
		every.push_back( index );
		alternate.push_back( index % 2 == 0 ? 1.02 : 0.98 );
	}
	std::vector<double> original;
	InverseOutcome const got = perturbFreeAndSolve( m, every, alternate, targetR,
	                                                targetZ, original );
	BOOST_TEST_REQUIRE( got.solved );

	int const rows = static_cast<int>( targetR.size() ) - 1;
	BOOST_TEST( solver.inverseNullity() >= static_cast<int>( every.size() ) - rows );
	BOOST_TEST( got.worstTarget
	            < meq::tests::scalarAgreementFloor( 1.0e-11, solver.psiAxis() ) );
	BOOST_TEST( got.worstCurrent > 1.0e-4 );
}

/*
 * THE REFUSALS, each of which is a configuration that would otherwise converge
 * to something nobody asked for or move nothing at all.
 */
BOOST_AUTO_TEST_CASE( theInverseRefusesWhatItCannotPose )
{
	ensureMesh();
	Machine m = buildMachine();
	meq::GradShafranovSolver &solver = *m.solver;
	using Circuit = meq::GradShafranovSolver::FreeCircuit;

	// The setter's own: malformed circuits and weights.
	BOOST_CHECK_THROW( solver.setFreeCurrents( nullptr, { { "a", { 0 }, { 1.0 } } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, {} ), std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, { { "a", {}, {} } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, { { "a", { 0, 1 }, { 1.0 } } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, { { "a", { 99 }, { 1.0 } } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, { { "a", { 0 }, { 0.0 } } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils,
		{ Circuit{ "a", { 0 }, { 1.0 } }, Circuit{ "b", { 0 }, { 1.0 } } } ),
	                   std::invalid_argument );
	// F1A and F2A do not carry one current between them at multipliers 1, 1.
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils,
		{ { "a", { coilNamed( m, "F1A" ), coilNamed( m, "F2A" ) }, { 1.0, 1.0 } } } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.setFreeCurrents( m.coils, { { "a", { 0 }, { 1.0 } } },
	                                           { -1.0 } ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( solver.addIsofluxTarget( 1.0, 0.0, 1.0, 0.0 ),
	                   std::invalid_argument );
	BOOST_TEST( !solver.inverseIsPosed() );

	// solve()'s: the pairings only a whole configuration shows.
	solver.addIsofluxTarget( 1.2, 0.5, 1.2, -0.5 );
	BOOST_CHECK_THROW( solver.solve(), std::logic_error );

	solver.clearInverse();
	solver.setFreeCurrents( m.coils, { { "F1A", { 0 }, { 1.0 } } } );
	BOOST_CHECK_THROW( solver.solve(), std::logic_error );

	// A set the source does not carry moves nothing the residual reads.
	auto stranger = meq::makeCoilSet( m.config->getCoils(),
	                                  m.config->getSource().permeability() );
	solver.clearInverse();
	solver.setFreeCurrents( stranger, { { "F1A", { 0 }, { 1.0 } } } );
	solver.addIsofluxTarget( 1.2, 0.5, 1.2, -0.5 );
	BOOST_CHECK_THROW( solver.solve(), std::logic_error );

	solver.clearInverse();
	BOOST_TEST( !solver.inverseIsPosed() );
}

/*
 * IV-4: THE X-POINT ROWS, BOTH KINDS. INVERSE-AND-RECONSTRUCTION-PLAN.md
 * section 9.
 *
 * The forward answer's X-point target is moved by a stated distance, and two
 * lower coils are freed against its two rows -- a SQUARE problem, so the
 * currents that reach it are unique and the two kinds of row can be compared
 * on what they return:
 *
 *   * row four, addXPointTarget(), pins the saddle XP-3 is TRACKING;
 *   * row three, addCriticalPointTarget(), asks for a critical point there.
 *
 * Asserted: row four puts the tracked X-point on the target, and the
 * equilibrium survives RELEASE -- the same currents solved forward, which is
 * IV-6's chain on one machine. Row three, from the same start, reaches the
 * same currents. The currents themselves are REPORTED rather than asserted
 * against anything, since nothing independent says what they should be.
 *
 * AND THE CONTROL: a single free coil far from the null cannot move it to an
 * arbitrary point, so the same target with only that coil free must end at a
 * BEST FIT that says the target was not met.
 * Without it, a pass could be the X-point unknowns obeying their own rows
 * while the currents did nothing.
 */
BOOST_AUTO_TEST_CASE( theXPointTargetMovesTheTrackedSaddle )
{
	ensureMesh();
	Machine m = buildMachine();
	int sweeps = 0;
	BOOST_TEST_REQUIRE( solveWithSweeps( m, sweeps ) );

	// THE FORWARD SOLVER IS KEPT ALIVE: the saved field borrows its space.
	std::unique_ptr<meq::GradShafranovSolver> const forwardSolver =
		std::move( m.solver );
	mfem::GridFunction const forward( forwardSolver->potential() );
	double const forwardAxis = forwardSolver->psiAxis();
	double const forwardXR = forwardSolver->xPointR();
	double const forwardXZ = forwardSolver->xPointZ();
	double const targetR = forwardXR + 0.015;
	double const targetZ = forwardXZ - 0.015;
	std::printf( "forward X ( %.9f, %.9f ), target ( %.9f, %.9f ), moved %.4f m\n",
	             forwardXR, forwardXZ, targetR, targetZ,
	             std::hypot( targetR - forwardXR, targetZ - forwardXZ ) );

	std::vector<std::size_t> const pair = { coilNamed( m, "F4B" ),
	                                        coilNamed( m, "F8B" ) };
	std::vector<double> original;
	for ( std::size_t index : pair )
		original.push_back( m.coils->coil( index ).current() );

	/// A fresh solver at the forward answer, on its frozen support -- the
	/// source keeps the last sweep's edge -- with the original currents.
	auto restart = [ & ]() -> meq::GradShafranovSolver &
	{
		for ( std::size_t k = 0; k < pair.size(); ++k )
			m.coils->setCurrent( pair[ k ], original[ k ] );
		m.solver = makeSolver( m, forwardAxis, forwardXR, forwardXZ );
		m.solver->prepare();
		m.solver->setPlasmaSupportFrozen( false );
		m.solver->refreshPlasmaComponent( forward );
		m.solver->setPlasmaSupportFrozen( true );
		m.solver->setInitialGuess( forward );
		return *m.solver;
	};
	auto freePair = [ & ]( meq::GradShafranovSolver &solver )
	{
		std::vector<meq::GradShafranovSolver::FreeCircuit> circuits;
		for ( std::size_t index : pair )
			circuits.push_back( { m.config->getCoils().coils[ index ].name,
			                      { index }, { 1.0 } } );
		solver.setFreeCurrents( m.coils, circuits );
	};
	auto miss = [ & ]( meq::GradShafranovSolver const &solver )
	{
		return std::hypot( solver.xPointR() - targetR, solver.xPointZ() - targetZ );
	};
	auto report = [ & ]( meq::GradShafranovSolver const &solver, char const *label )
	{
		std::vector<double> const currents = solver.freeCircuitCurrents();
		std::printf( "%s: %d Newton iterations, rank %d nullity %d, X ( %.9f, %.9f )"
		             " off by %.3e m\n", label, solver.newtonIterations(),
		             solver.inverseRank(), solver.inverseNullity(),
		             solver.xPointR(), solver.xPointZ(), miss( solver ) );
		for ( std::size_t k = 0; k < currents.size(); ++k )
			std::printf( "   %-4s  %+.10e -> %+.10e  ( %+.3f%% )\n",
			             m.config->getCoils().coils[ pair[ k ] ].name.c_str(),
			             original[ k ], currents[ k ],
			             100.0*( currents[ k ] - original[ k ] )/original[ k ] );
		for ( double v : solver.newtonResiduals() )
			std::printf( "   %.6e\n", v );
		std::printf( "   dR %.3e dZ %.3e, target rows", solver.xPointR() - targetR,
		             solver.xPointZ() - targetZ );
		for ( double v : solver.targetResiduals() )
			std::printf( " %.3e", v );
		std::printf( "\n" );
	};

	// ROW FOUR.
	std::vector<double> rowFour;
	double rowFourMiss = 0.0, rowFourAxis = 0.0, releasedMiss = 0.0;
	int releasedIterations = 0;
	{
		meq::GradShafranovSolver &solver = restart();
		freePair( solver );
		solver.addXPointTarget( targetR, targetZ );
		solver.solve();
		report( solver, "row four" );
		rowFour = solver.freeCircuitCurrents();
		rowFourMiss = miss( solver );
		rowFourAxis = solver.psiAxis();

		// RELEASED: the solved currents are in the set, so a forward solve
		// from here is the same machine with nothing free.
		solver.clearInverse();
		solver.solve();
		releasedMiss = miss( solver );
		releasedIterations = solver.newtonIterations();
		std::printf( "released: %d Newton iterations, X off by %.3e m, psi_ax "
		             "%.10e against %.10e\n", releasedIterations, releasedMiss,
		             solver.psiAxis(), rowFourAxis );
	}

	// ROW THREE, from the same start.
	std::vector<double> rowThree;
	double rowThreeMiss = 0.0;
	{
		meq::GradShafranovSolver &solver = restart();
		freePair( solver );
		solver.addCriticalPointTarget( targetR, targetZ );
		solver.solve();
		report( solver, "row three" );
		rowThree = solver.freeCircuitCurrents();
		rowThreeMiss = miss( solver );
	}
	double currentsAgree = 0.0;
	for ( std::size_t k = 0; k < pair.size(); ++k )
		currentsAgree = std::max( currentsAgree,
		                          std::abs( rowThree[ k ] - rowFour[ k ] )
		                          /std::abs( rowFour[ k ] ) );
	std::printf( "row three against row four: currents agree to %.3e\n",
	             currentsAgree );

	// THE CONTROL: one distant coil cannot place the null.
	bool controlConverged = true;
	bool controlMet = true;
	double controlMiss = 0.0;
	{
		meq::GradShafranovSolver &solver = restart();
		std::size_t const far = coilNamed( m, "F6A" );
		solver.setFreeCurrents( m.coils, { { "F6A", { far }, { 1.0 } } } );
		solver.addXPointTarget( targetR, targetZ );
		solver.setNewtonControl(
			meq::newtonTolerance( m.config->getSolver().newtonRelativeTolerance ),
			meq::newtonTolerance( m.config->getSolver().newtonAbsoluteTolerance ),
			30 );
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			controlConverged = false;
		}
		controlMiss = miss( solver );
		controlMet = solver.inverseTargetsMet();
		std::printf( "control, F6A alone: %s after %d iterations, targets %s, X off "
		             "by %.3e m\n",
		             controlConverged ? "converged" : "did not converge",
		             solver.newtonIterations(), controlMet ? "MET" : "not met",
		             controlMiss );
	}

	/*
	 * MEASURED: a 2.1 cm move takes F4B down 13.4% and F8B up 58.3%, in five
	 * Newton steps on either row. Both rows put the tracked saddle on the
	 * target to the last bit -- the fixed point is the target itself, so the
	 * last update rounds onto it -- and agree on the currents to 5.7e-11.
	 * Released, the solve converges from the ORIGINAL forward field to the
	 * same saddle, 3.6e-12 m off, with psi_ax to 4e-11. The control
	 * converges to the BEST FIT one coil allows, its targets not met and the
	 * null 1.2 cm short -- the outcome IV-5's stationarity stop exists for.
	 *
	 * Mutated: H's sign flipped in row four's corner -- no convergence; row
	 * three reading q_r twice -- rank 16, nullity 1, the saddle 2.6 cm off.
	 * Negating a whole row with its value is NOT a mutation: -G = 0 is the
	 * same equation, and the test rightly cannot tell.
	 */
	// The saddle is ROOTED on q_h = 0 by the border's own rows, so the format
	// costs eps times the coordinates and not sqrt( eps ); the currents are
	// relative, so their scale is 1. Both floors are zero or inert in double.
	double const missBound = 1.0e-9 + meq::tests::rootedPositionFloor( targetR );
	BOOST_TEST( rowFourMiss < missBound );
	BOOST_TEST( releasedMiss < missBound );
	BOOST_TEST( rowThreeMiss < missBound );
	// The CURRENTS are another matter: two whole solves through different rows,
	// each reaching them through a corner of condition about 3.4e+08 (M-250).
	// Double's rows agree to 5.7e-11, which is 2.6e+05 eps; single's scatter
	// over 2.5e-06 to 3.5e-05 with MKL_NUM_THREADS, red at 1 and green at 4
	// under 64 eps. exactAgreementFloor()'s 4096 eps is a whole solve compared
	// two ways. A wrong row moves them by tens of per cent, and the floor is
	// zero in double. M-275.
	BOOST_TEST( currentsAgree < std::max( 1.0e-8, meq::tests::exactAgreementFloor( 1.0 ) ) );
	BOOST_TEST( controlConverged );
	BOOST_TEST( ( !controlMet && controlMiss > 1.0e-3 ),
	            "a single distant coil placed the null, so the pass above does "
	            "not show the currents doing the work" );
}

/*
 * IV-5: THROUGH THE DRIVER, WITH THE SUPPORT SWEEPS.
 * INVERSE-AND-RECONSTRUCTION-PLAN.md section 9.
 *
 * examples/machine-f-diiid-inverse.toml is the forward DIII-D file with four
 * currents Free at seeds 2 to 3 per cent off, and [inverse] targets read off
 * the forward equilibrium's separatrix and X-point. The driver runs it FROM
 * THE FILE'S COLD GUESS, inside the support sweeps, which IV-3 and IV-4 did
 * not: there the support was frozen at an answer.
 *
 * Asserted, in the manner of theDriverRunsTheAdaptiveLoop:
 *   * the driver exits 0, says the targets are met, and the support SETTLES --
 *     section 6.2's hazard, the support depending on currents that move;
 *   * the same file posed through the library by meq::poseInverse() and swept
 *     as the driver sweeps reaches the driver's psi_h dof for dof and the same
 *     sweep count;
 *   * the four currents come back to the forward file's.
 *
 * THE CASE THAT FORCED THE DESIGN IS THE FIRST SWEEP. Its support is frozen at
 * the cold guess's edge, so the targets -- read off the converged
 * equilibrium -- are not exactly reachable on it. A corner fitting every row
 * at once stalled there with the equilibrium's own border rows violated at
 * 1e-08 to 3e-07; the border rows are now constraints and the sweep stops at
 * the best fit, which the next sweep's support improves. M-252.
 */
namespace
{
	std::string const inverseFile = meq::test::sourcePath( "examples/machine-f-diiid-inverse.toml" );

	/// Run the shipped wrapper, capturing everything it prints.
	int runDriver( std::string const &file, std::string const &log )
	{
		// --gf: the case below reads the run's .mesh and _psi.gf back, and the
		// driver writes only the .nc unless asked.
		std::string const command = "MKL_NUM_THREADS=1 " + std::string( MEQ_RUN_PATH )
		                            + " " + file + " --gf > " + log + " 2>&1";
		int const status = std::system( command.c_str() );
		return WIFEXITED( status ) ? WEXITSTATUS( status ) : -1;
	}

	std::string readAll( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( std::istreambuf_iterator<char>( stream ),
		                    std::istreambuf_iterator<char>() );
	}

	/// One variable's values from `ncdump -p 9,17`, as the text between its
	/// `name =` and `;`, split on commas, with string quotes and blanks
	/// dropped.
	std::vector<std::string> ncdumpValues( std::string const &path,
	                                       std::string const &name )
	{
		std::string const scratch = "inverse-solve-ncdump.txt";
		std::string const command = "ncdump -p 9,17 -v " + name + " " + path
		                            + " > " + scratch + " 2>/dev/null";
		std::vector<std::string> out;
		int const status = std::system( command.c_str() );
		std::string const text = readAll( scratch );
		std::remove( scratch.c_str() );
		if ( status != 0 )
			return out;
		std::size_t at = text.find( "\ndata:" );
		if ( at == std::string::npos )
			return out;
		at = text.find( "\n " + name + " =", at );
		if ( at == std::string::npos )
			return out;
		std::size_t const end = text.find( ';', at );
		std::string const body = text.substr( at + name.size() + 4,
		                                      end - at - name.size() - 4 );
		std::string one;
		for ( char const c : body + "," )
		{
			if ( c == ',' )
			{
				out.push_back( one );
				one.clear();
			}
			else if ( c != '"' && c != ' ' && c != '\n' )
				one += c;
		}
		return out;
	}
}

BOOST_AUTO_TEST_CASE( theDriverSolvesAnInverseProblem )
{
	ensureMesh( inverseFile );
	std::string const log = "inverse-solve-driver.log";
	int const exit = runDriver( inverseFile, log );
	std::string const printed = readAll( log );
	std::remove( log.c_str() );
	BOOST_TEST_REQUIRE( exit == 0, "the driver exited " << exit << ":\n" << printed );
	BOOST_TEST( printed.find( "the targets are met" ) != std::string::npos );
	std::size_t const at = printed.find( "the plasma support took " );
	BOOST_TEST_REQUIRE( at != std::string::npos );
	int const driverSweeps = std::atoi( printed.c_str() + at + 24 );
	BOOST_TEST( printed.find( "sweeps and settled", at ) != std::string::npos,
	            "the driver's support did not settle" );

	// THE SAME FILE THROUGH THE LIBRARY, by the one translation the driver uses.
	Machine m = buildMachine( inverseFile );
	meq::poseInverse( *m.solver, m.coils, m.config->getCoils(),
	                  m.config->getInverse() );
	int librarySweeps = 0;
	BOOST_TEST_REQUIRE( solveWithSweeps( m, librarySweeps ) );
	BOOST_TEST_REQUIRE( m.solver->inverseTargetsMet() );
	BOOST_TEST( librarySweeps == driverSweeps );

	mfem::Mesh storedMesh( "machine-f-diiid-inverse.mesh", 1, 1 );
	std::ifstream stream( "machine-f-diiid-inverse_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction const stored( &storedMesh, stream );
	BOOST_TEST_REQUIRE( stored.Size() == m.solver->potential().Size() );
	double worst = 0.0;
	double largest = 0.0;
	for ( int i = 0; i < stored.Size(); ++i )
	{
		worst = std::max( worst, static_cast<double>(
			std::abs( stored( i ) - m.solver->potential()( i ) ) ) );
		largest = std::max( largest, static_cast<double>( std::abs( stored( i ) ) ) );
	}

	meq::Configuration const forward( machineFile );
	std::vector<meq::FreeCircuit> const circuits =
		meq::makeFreeCircuits( m.config->getCoils() );
	std::vector<double> const solved = m.solver->freeCircuitCurrents();
	double worstCurrent = 0.0;
	for ( std::size_t k = 0; k < circuits.size(); ++k )
	{
		double const original =
			forward.getCoils().coils[ circuits[ k ].coils.front() ].current;
		worstCurrent = std::max( worstCurrent,
		                         std::abs( solved[ k ] - original )/std::abs( original ) );
	}
	std::printf( "driver against library: %d and %d sweeps, psi_h to %.3e over %d "
	             "dofs ( max |psi_h| %.3e ), currents to the forward file's %.3e, "
	             "rank %d nullity %d\n", driverSweeps, librarySweeps, worst,
	             stored.Size(), largest, worstCurrent, m.solver->inverseRank(),
	             m.solver->inverseNullity() );

	BOOST_TEST( worst < meq::tests::scalarAgreementFloor( 1.0e-12, largest ) );
	BOOST_TEST( worstCurrent < meq::tests::scalarAgreementFloor( 1.0e-8, 1.0 ) );
	BOOST_TEST( m.solver->inverseNullity() == 0 );

	// THE FILE'S MATERIAL TABLE CARRIES THE SOLVED CURRENTS, not the
	// currents the run started from: each free coil's material_current is
	// its circuit's solved current times its multiplier.
	std::vector<std::string> const names =
		ncdumpValues( "machine-f-diiid-inverse.nc", "material_name" );
	std::vector<std::string> const written =
		ncdumpValues( "machine-f-diiid-inverse.nc", "material_current" );
	BOOST_TEST_REQUIRE( names.size() == written.size() );
	BOOST_TEST_REQUIRE( !names.empty(), "the .nc has no material table" );
	int checked = 0;
	for ( std::size_t k = 0; k < circuits.size(); ++k )
		for ( std::size_t j = 0; j < circuits[ k ].coils.size(); ++j )
		{
			std::string const &coil =
				m.config->getCoils().coils[ circuits[ k ].coils[ j ] ].name;
			auto const row = std::find( names.begin(), names.end(), coil );
			BOOST_TEST_REQUIRE( ( row != names.end() ),
			                    coil << " has no row in the material table" );
			double const current = std::strtod(
				written[ static_cast<std::size_t>( row - names.begin() ) ].c_str(),
				nullptr );
			double const expected = circuits[ k ].multipliers[ j ]*solved[ k ];
			BOOST_TEST( std::abs( current - expected )
			            < meq::tests::scalarAgreementFloor( 1.0e-8, std::abs( expected ) ),
			            coil << "'s material_current is " << current
			            << " and its solved current " << expected );
			++checked;
		}
	BOOST_TEST( checked == 4 );
}

/*
 * A TARGET SET THE FREE CURRENTS CANNOT MEET IS A BEST FIT, AND THE RUN SAYS SO
 * IN ITS EXIT STATUS. Plan section 4.5.
 *
 * The same targets with only F6A free -- the other three held at their
 * perturbed values -- cannot be met: one current cannot put a separatrix back
 * through 22 points. The run must converge anyway, to the best fit with every
 * other constraint of the equilibrium met, write everything, print that the
 * targets were not met, and exit 4 rather than 0.
 */
BOOST_AUTO_TEST_CASE( aTargetSetTheCurrentsCannotMeetIsABestFitAndExitsFour )
{
	ensureMesh( inverseFile );
	// Its tables and guess by absolute paths, since the variant is written into
	// the working directory and away from them.
	std::string text = meq::test::shippedConfig( "examples/machine-f-diiid-inverse.toml" );
	for ( char const *name : { "F2A", "F2B", "F6B" } )
	{
		std::size_t const block = text.find( std::string( "Name = \"" ) + name + "\"" );
		BOOST_TEST_REQUIRE( block != std::string::npos );
		std::size_t const free = text.find( "Free = true\n", block );
		BOOST_TEST_REQUIRE( free != std::string::npos );
		text.erase( free, std::string( "Free = true\n" ).size() );
	}
	std::size_t const prefix = text.find( "Prefix = \"machine-f-diiid-inverse\"" );
	BOOST_TEST_REQUIRE( prefix != std::string::npos );
	text.replace( prefix, std::string( "Prefix = \"machine-f-diiid-inverse\"" ).size(),
	              "Prefix = \"inverse-solve-best-fit\"" );
	std::string const file = "inverse-solve-best-fit.toml";
	{
		std::ofstream out( file );
		out << text;
	}

	std::string const log = "inverse-solve-best-fit.log";
	int const exit = runDriver( file, log );
	std::string const printed = readAll( log );
	std::remove( log.c_str() );
	std::size_t const at = printed.find( "inverse:" );
	std::printf( "best fit: exit %d\n%s\n", exit,
	             at == std::string::npos ? printed.c_str()
	                                     : printed.substr( at, 600 ).c_str() );
	BOOST_TEST( exit == 4 );
	BOOST_TEST( printed.find( "THE TARGETS WERE NOT MET" ) != std::string::npos );
	BOOST_TEST( printed.find( "inverse-solve-best-fit.nc" ) != std::string::npos,
	            "a best fit must still be written" );
	std::remove( file.c_str() );
}
