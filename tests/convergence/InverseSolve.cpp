/*
 * IV-3: COIL CURRENTS AS BORDER UNKNOWNS, ONE MACHINE, ISOFLUX ONLY.
 * INVERSE-AND-RECONSTRUCTION-PLAN.md section 9.
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

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/Config.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/Sampler.hpp"
#include "meq/SourceFactory.hpp"
#include "meq/WarmStart.hpp"

namespace
{
	char const *const machineFile = "examples/machine-f-diiid.toml";

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

	/// The mesh is a build product of the file, so it is made the way the
	/// driver's users make it rather than assumed to be on disk.
	void ensureMesh()
	{
		std::string const command = std::string( MEQ_RUN_PATH )
		                            + " --mesh-only " + machineFile
		                            + " > /dev/null 2>&1";
		BOOST_TEST_REQUIRE( std::system( command.c_str() ) == 0,
		                    "meq-run --mesh-only failed on " << machineFile );
	}

	/// The assembly apps/meq.cpp performs for this file, in its order: cut D_h
	/// from the half disc at Gamma, wrap the plasma source around the coils and
	/// confine it ON THE WRAPPER, then the X-point border, the prescribed
	/// current and the exterior coupling, and the stored guess interpolated in.
	Machine buildMachine()
	{
		Machine m;
		m.config = std::make_unique<meq::Configuration>( machineFile );
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
		                            sourceConfig.permeability(), machineFile );
		m.plasma = meq::makeNormalisedSource( sourceConfig, machineFile );
		m.source = std::make_shared<meq::CoilAugmentedNormalisedSource>( m.plasma,
		                                                                 m.coils );
		// ON THE WRAPPER, AFTER IT EXISTS -- plan section 6.1 and M-82.
		m.source->setPlasmaSupport( sourceConfig.confinesToPlasma() );
		m.exterior = std::make_unique<meq::ExteriorDtN>( exteriorConfig.centreZ,
		                                                exteriorConfig.radius,
		                                                exteriorConfig.modes );

		meq::XPointConfig const &xpoint = m.config->getBoundary().xpoint;
		BOOST_TEST_REQUIRE( xpoint.given );

		m.solver = std::make_unique<meq::GradShafranovSolver>(
			*m.sub, m.config->getDiscretisation().polynomialDegree );
		m.solver->setNewtonControl(
			meq::newtonTolerance( m.config->getSolver().newtonRelativeTolerance ),
			meq::newtonTolerance( m.config->getSolver().newtonAbsoluteTolerance ),
			m.config->getSolver().newtonMaxIterations );
		m.solver->setExtension( *m.path, m.gammaHMarker );
		m.solver->setSource( *m.source, sourceConfig.psiAxisGuess() );
		m.solver->setXPointBoundary( xpoint.radius, xpoint.z );
		m.solver->setPlasmaCurrent( sourceConfig.permeability()
		                            *sourceConfig.plasmaCurrent() );
		m.solver->setBoundaryData( m.zero );
		m.solver->setExteriorCoupling( *m.exterior );

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
		BOOST_TEST_REQUIRE( worst < 1.0e-12 );
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

		for ( double v : solver.isofluxResiduals() )
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
	BOOST_TEST( got.worstTarget < 1.0e-11 );
	BOOST_TEST( got.worstCurrent < 1.0e-9 );
	BOOST_TEST( std::abs( solver.psiAxis() - forwardAxis )
	            < 1.0e-9*std::abs( forwardAxis ) );
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
	BOOST_TEST( all.worstTarget < 1.0e-11 );
	BOOST_TEST( all.worstCurrent < 1.0e-5 );
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
	BOOST_TEST( got.worstTarget < 1.0e-11 );
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
