/*
 * XP-2: psi_bnd FROM THE LOCATED X-POINT, AS AN OUTER FIXED POINT.
 * FREE-BOUNDARY-PLAN.md section 10.6.
 *
 * FB-3's setBoundaryFluxPoint() pins psi_bnd = psi_h at a PRESCRIBED point,
 * which is right for a limiter -- the contact is a piece of hardware -- and
 * wrong for a divertor, whose X-point is a functional of the solution and moves
 * as Newton moves. XP-3 makes the X-point three more unknowns of the same
 * Newton; XP-2 is the halfway house that needs no new border at all:
 *
 *     solve with psi_bnd pinned at the CURRENT estimate of the X-point
 *     locate the saddle of the SOLVED q_h
 *     re-pin, re-solve
 *
 * an outer fixed point over ( r_X, z_X ), warm started from the previous sweep.
 *
 * THE DISCRETE CHOICE IS FROZEN WITHIN A SOLVE AND RE-DECIDED BETWEEN THEM,
 * which is section 10.5's answer to the one thing this stage cannot
 * differentiate: psi_bnd is a min over candidates and a min is not smooth, so
 * Newton is never asked to differentiate through it. Each solve is a smooth
 * problem with the bounding point FIXED.
 *
 *
 * WHY THIS RUNS ON THE MACHINE CASE AND NOT ON THE HALF-DISC FIXTURE, WHICH IS
 * A MEASURED FALSIFICATION AND NOT A PREFERENCE.
 *
 * FreeBoundaryCoupling.cpp's theTwoBorderSolveReportsATrueMagneticAxis is the
 * obvious host: it is already a physical free-boundary equilibrium with coils,
 * a prescribed current and ConfineToPlasma. It cannot be made diverted, and the
 * obstruction is structural rather than a matter of trying harder.
 *
 * psi vanishes IDENTICALLY on r = 0 -- the symmetry axis is fitted Dirichlet
 * boundary -- so the axis sits at normalised flux -psi_bnd/span. A NEGATIVE
 * psi_bnd therefore puts it at POSITIVE Psi, inside the plasma support, where
 * ConfineToPlasma does not switch gg' off; and F/r is mu_0 j_phi, so that is an
 * infinite toroidal current density on the axis. What the field then grows is
 * section 11.3's axis layer, and what CriticalPointFinder reports is the layer.
 * Measured on one such run: THIRTY critical points, every one at r < 0.1, a
 * ladder of maxima marching up the symmetry axis from z = -0.11 to z = +1.44,
 * and the "magnetic axis" at r = 0.022.
 *
 * That fixture reaches only psi_bnd = +9.5e-04 against psi_ax = 1.19e-02, a
 * margin of 8%, while a Shafranov vertical field alone contributes about
 * psi = B_v r^2/2 = -1.8e-02 at the limiter. The plasma's own flux outweighs
 * the coils' by a few per cent, so any conductor strong enough to make a null
 * tips psi_bnd negative first. Twelve configurations were run across three
 * routes -- an exterior filament swept in current, the same with the
 * equilibrium pair re-derived to hold the total vertical field, and a designed
 * three-conductor single null with the divertor meshed inside Omega -- and
 * every one of them ended at negative psi_bnd, at no off-axis saddle, or at a
 * saddle out in the vacuum ABOVE the plasma. Superposition predicted X-points
 * in three of six cases where the solve produced none, which is the other half
 * of it: the screen drops the plasma's response and the response is most of the
 * answer.
 *
 * examples/diverted-tokamak.toml IS THE MACHINE THIS RUNS ON INSTEAD, AND IT IS
 * BORROWED RATHER THAN DESIGNED. It is freegs4e's own A_testtokamak_classic:
 * an up-down asymmetric double null whose coil currents were SOLVED FOR by its
 * control system, which section 7.15 argues is what makes a free-boundary
 * Newton converge at all. Its psi_bndry IS the lower X-point's flux to every
 * digit and its psi_bnd/psi_ax is 0.39.
 *
 * A run of hand-designed single nulls on the LIMITED machine's geometry came
 * first and every one of them failed the same way -- negative psi_bnd and the
 * axis layer -- so the design machinery that used to live in this file is gone.
 * The finding it produced is kept: every design solve pinned psi_bnd at the
 * LIMITER, and a diverted machine's limiter is OUTSIDE its separatrix, so the
 * solve was being asked the wrong question rather than answering it badly.
 *
 *
 * THE INPUTS ARE CHECKED AND THE INNER SOLVE IS THE DEFECT. THIS CASE IS RED
 * AND WHAT FOLLOWS IS WHY.
 *
 * Both halves of the problem statement are verified independently of MEQ:
 *
 *   * the two profile tables reproduce the reference's own plasma current when
 *     integrated over the reference's own core -- 1.999667e+05 A against
 *     2.0e+05, so the span conversion is right and a healthy solve must report
 *     a profile scale of 1;
 *   * the Green's-function guess reproduces psi at the reference's magnetic
 *     axis to 4.1e-04 relative ( 8.268402e-02 against 8.271751e-02 ) and at its
 *     active X-point to 4.8e-04 ( 3.238847e-02 against 3.240413e-02 ), so the
 *     iteration starts on the physical branch.
 *
 * It does not stay there. Measured, one key changed at a time:
 *
 *   pin      k  guess    outcome
 *   X-point  2  33^2     99 steps to || r ||/|| r_0 || = 4.1e-11, psi_ax
 *                        8.532768e-02 -- and psi_h AT the reference axis reads
 *                        1.32e-02 against the reference's 8.27e-02. The located
 *                        axis is ( 2.2626, 0.5413 ), hard against Gamma
 *   X-point  3  33^2     14 steps, axis at ( 1.4691, 1.4401 )
 *   X-point  2  129^2    the residual falls to 3.4e-04 and STALLS there for 150
 *                        steps -- no runaway and no convergence -- with the
 *                        axis already at ( 2.2590, 0.3105 ) by step 13
 *   (1.68,0) 2  33^2     a LIMIT CYCLE: || r || descends to 7.1e-03 and jumps
 *                        back to 2.7e-02, six times over the 200-step cap. That
 *                        pin is inside the separatrix, where the reference
 *                        carries psi = 4.52e-02, so it is the easy question
 *
 * AND THE TRACE SAYS WHERE IT BREAKS. Printing the located axis at every
 * evaluation of the first run: the first fourteen are right -- ( 1.3636,
 * 0.0008 ) carrying psi = 8.077e-02, with the span climbing 4.14e-02, 6.20e-02,
 * 7.24e-02, ... , 8.268e-02 towards the reference's 8.2718e-02 -- and then
 * psi_ax collapses, 7.24e-02, 3.78e-02, 2.28e-02, with the axis following it
 * out to ( 2.2042, 0.2313 ). The solve is nearly there and then leaves.
 *
 * What that rules out: the profile conversion, the guess, the polynomial
 * degree, and the axis competition -- a conductor's O-point is excluded on both
 * paths now ( meq::Source::conductors ) and the runaway axis is in no coil.
 * What is left is the pair this stage exists to separate: a plasma support that
 * MOVES within the Newton, and a psi_bnd pinned exactly at a saddle, where the
 * enclosed area's derivative is not bounded. Section 10.5 says to freeze the
 * discrete choice within a solve, and the loop below does; what it does not yet
 * do is freeze the SUPPORT, and PlasmaConnectivity is a library control with no
 * TOML key, so that experiment belongs here rather than in an example.
 *
 * PER THE TESTING STANCE, THIS ASSERTS THE BEHAVIOUR WANTED AND FAILS UNTIL IT
 * IS THERE. It is not evidence that XP-2's loop is wrong -- the loop is the
 * eight lines below and nothing in the list above is about it -- it is evidence
 * that the loop has no equilibrium to iterate on yet.
 */
#define BOOST_TEST_MODULE XPointOuter
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/Config.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/SourceFactory.hpp"
#include "meq/WarmStart.hpp"

namespace
{
	/// The diverted machine, taken whole: mesh, profiles, coil currents,
	/// limiter contact, exterior expansion and guess. Nothing here is designed
	/// or tuned -- see the header for why that matters.
	char const *const machineFile = "examples/diverted-tokamak.toml";

	/// freegs4e's own answer for this machine, from
	/// tools/freegs4e-benchmark/A_testtokamak_classic.json. The ACTIVE null is
	/// the lower one; the upper sits 3.5e-03 further out.
	double const referenceXPointR = 1.093144118182931;
	double const referenceXPointZ = -0.6039650838688502;
	double const referencePsiBoundary = 3.240412550738516e-02;
	double const referencePsiAxis = 8.271751444840184e-02;
	double const referenceAxisR = 1.351273;
	double const referenceAxisZ = 0.062226;

	/// Negative inside. The semicircle about the axis, which is what makes the
	/// exterior expansion legal at all; the flat side IS r = 0 and stays
	/// ordinary fitted boundary.
	mfem::PositionFunction semicircle( double radius, double centreZ )
	{
		return [ radius, centreZ ]( mfem::Vector const &x )
		{
			return std::hypot( x( 0 ), x( 1 ) - centreZ ) - radius;
		};
	}

	/// Everything one solve needs, kept alive together. The solver aliases the
	/// submesh, the path and the source; the submesh aliases the background.
	/// Nothing here may outlive anything below it.
	struct Machine
	{
		std::unique_ptr<meq::Configuration> config;
		std::unique_ptr<mfem::Mesh> background;
		std::unique_ptr<mfem::SubMesh> sub;
		std::unique_ptr<mfem::VertexConePath> path;
		mfem::Array<int> gammaHMarker;
		int gammaH = 0;

		std::shared_ptr<meq::CoilSet const> coils;
		std::shared_ptr<meq::NormalisedSource> plasma;
		std::shared_ptr<meq::CoilAugmentedNormalisedSource> source;
		std::unique_ptr<meq::ExteriorDtN> exterior;

		std::unique_ptr<mfem::Mesh> guessMesh;
		std::unique_ptr<mfem::GridFunction> guess;
		std::unique_ptr<mfem::GridFunction> carried;
		mfem::ConstantCoefficient zero{ 0.0 };

		std::unique_ptr<meq::GradShafranovSolver> solver;
		bool converged = false;
		double mu0 = 0.0;
		double limiterR = 0.0;
		double limiterZ = 0.0;
	};

	/// The machine of machineFile: the gmsh half-disc with its conductors
	/// meshed to, the two profile tables, the coil currents, the limiter
	/// contact, the exterior expansion and the reference initial guess.
	///
	/// The sequence follows apps/meq.cpp, which is the only place this assembly
	/// exists; a test that assembled it differently would be measuring its own
	/// wiring.
	Machine buildMachine()
	{
		Machine m;
		m.config = std::make_unique<meq::Configuration>( machineFile );

		meq::SourceConfig const &sourceConfig = m.config->getSource();
		m.mu0 = sourceConfig.permeability();

		m.background = std::make_unique<mfem::Mesh>(
			m.config->getMesh().file.c_str(), 1, 1 );

		// A MESH READ FROM A FILE HAS TO BE MEASURED RATHER THAN COMPUTED, and
		// a zero search length is not a loose tolerance -- it is
		// mfem::VertexConePath aborting on the first vertex of Gamma_h.
		double h = 0.0;
		for ( int e = 0; e < m.background->GetNE(); ++e )
			h = std::max( h, meq::elementDiameter( *m.background, e ) );

		meq::ExteriorConfig const &exteriorConfig =
			m.config->getBoundary().exterior;
		BOOST_TEST_REQUIRE( exteriorConfig.given,
		                    "the machine file names no [boundary.exterior]" );
		mfem::PositionFunction const levelSet =
			semicircle( exteriorConfig.radius, exteriorConfig.centreZ );

		// The gmsh mesh carries MATERIAL attributes -- the four conductors and
		// the limiter region -- and MarkLevelSetSubdomain overwrites them to
		// select the subdomain, so they are saved and put back on both meshes.
		std::vector<int> original( static_cast<std::size_t>(
			m.background->GetNE() ) );
		for ( int e = 0; e < m.background->GetNE(); ++e )
			original[ static_cast<std::size_t>( e ) ] =
				m.background->GetAttribute( e );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE(
			mfem::MarkLevelSetSubdomain( *m.background, levelSet, 0.0, marker,
			                             1 ) > 0,
			"the subdomain inside Gamma is empty" );
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
				m.sub->SetAttribute(
					e, original[ static_cast<std::size_t>( parent[ e ] ) ] );
			m.sub->SetAttributes();
			for ( int e = 0; e < m.background->GetNE(); ++e )
				m.background->SetAttribute(
					e, original[ static_cast<std::size_t>( e ) ] );
			m.background->SetAttributes();
		}

		m.gammaH = m.sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( m.gammaH > parentBoundaryMax,
			"no boundary was generated when cutting the subdomain, so Gamma_h "
			"names a fitted edge" );
		m.path = std::make_unique<mfem::VertexConePath>( *m.sub, m.gammaH,
		                                                levelSet, 6.0*h );
		m.gammaHMarker.SetSize( m.gammaH );
		m.gammaHMarker = 0;
		m.gammaHMarker[ m.gammaH - 1 ] = 1;
		return m;
	}

	void completeMachine( Machine &m )
	{
		meq::SourceConfig const &sourceConfig = m.config->getSource();

		// THE FILE'S OWN [[coils]], THROUGH THE FACTORY THE DRIVER CALLS. The
		// currents are freegs4e's control system's, so this case borrows a
		// machine somebody else's code balanced rather than one built here --
		// which is the whole reason section 10 stopped designing single nulls
		// by hand.
		m.coils = meq::makeCoilSet( m.config->getCoils(),
		                            sourceConfig.permeability(), machineFile );

		m.plasma = meq::makeNormalisedSource( sourceConfig, machineFile );
		// SET ON THE PLASMA SOURCE, and the wrapper forwards it: a conductor
		// sits in the vacuum by construction, so confining the coil term to the
		// plasma would switch off every coil in the machine.
		m.plasma->setPlasmaSupport( true );
		m.source = std::make_shared<meq::CoilAugmentedNormalisedSource>(
			m.plasma, m.coils );

		meq::ExteriorConfig const &exteriorConfig =
			m.config->getBoundary().exterior;
		m.exterior = std::make_unique<meq::ExteriorDtN>( exteriorConfig.centreZ,
		                                                exteriorConfig.radius,
		                                                exteriorConfig.modes );

		meq::LimiterConfig const &limiterConfig = m.config->getBoundary().limiter;
		BOOST_TEST_REQUIRE( limiterConfig.given,
		                    "the machine file names no [boundary.limiter]" );
		m.limiterR = limiterConfig.r;
		m.limiterZ = limiterConfig.z;

		m.solver = std::make_unique<meq::GradShafranovSolver>(
			*m.sub, m.config->getDiscretisation().polynomialDegree );
		m.solver->setNewtonControl(
			m.config->getSolver().newtonRelativeTolerance,
			m.config->getSolver().newtonAbsoluteTolerance,
			m.config->getSolver().newtonMaxIterations );
		m.solver->setSource( *m.source, sourceConfig.psiAxisGuess() );
		m.solver->setPlasmaCurrent( m.mu0
		                            *sourceConfig.getMHD().plasmaCurrent );
		m.solver->setBoundaryData( m.zero );
		m.solver->setExtension( *m.path, m.gammaHMarker );
		m.solver->setExteriorCoupling( *m.exterior );

		/*
		 * THE GUESS IS PART OF THE PROBLEM STATEMENT AND NOT AN OPTIMISATION.
		 * examples/limited-tokamak.toml records that a cold bump start wanders
		 * for 200 iterations around || r || = 1.3 on this machine and never
		 * converges. The stored guess is freegs4e's own equilibrium
		 * reconstructed by Green's functions, on its own mesh, so it is the
		 * INTERPOLATING restart -- meq::FieldTransfer -- rather than the exact
		 * one.
		 */
		meq::InitialGuessConfig const &guessConfig = m.config->getInitialGuess();
		BOOST_TEST_REQUIRE( !guessConfig.file.empty(),
		                    "the machine file names no [initialguess]" );
		m.guessMesh = std::make_unique<mfem::Mesh>( guessConfig.meshFile.c_str(),
		                                            1, 1 );
		std::ifstream stream( guessConfig.file );
		BOOST_TEST_REQUIRE( stream.good(),
		                    "cannot open " << guessConfig.file );
		m.guess = std::make_unique<mfem::GridFunction>( m.guessMesh.get(),
		                                                stream );

		meq::FieldTransfer transfer( *m.guessMesh );
		m.carried = std::make_unique<mfem::GridFunction>(
			m.solver->potential().FESpace() );
		transfer.transfer( *m.guess, m.zero, *m.carried );
		m.solver->setInitialGuess( *m.carried );
	}

	/// Pin psi_bnd at ( r, z ) and solve. setBoundaryFluxPoint() clears the
	/// prepared flag, so the same solver may be re-pinned and re-solved -- which
	/// is the whole of XP-2's outer loop.
	///
	/// @param warm  psi_h of the previous sweep, or null to keep whatever guess
	///              is already set. CLAUDE.md's warm-start trap does not bite:
	///              solve() takes its convergence reference from the COLD
	///              iterate precisely so a good guess cannot make the target
	///              unreachable.
	bool solveAt( Machine &m, double r, double z,
	              mfem::GridFunction const *warm = nullptr )
	{
		m.solver->setBoundaryFluxPoint( r, z );
		if ( warm != nullptr )
			m.solver->setInitialGuess( *warm );
		try
		{
			m.solver->solve();
			m.converged = !m.solver->newtonResiduals().empty()
			              && m.solver->newtonResiduals().back() < 1.0e-8;
		}
		catch ( std::exception const & )
		{
			m.converged = false;
		}
		return m.converged;
	}

	/// How far inside its own element a located point sits, in reference
	/// coordinates: 0 on a face of the reference triangle, 1/3 at its centroid.
	///
	/// **A ROOT ON A MESH LINE IS NOT ONE ROOT, AND THIS IS THE DIAGNOSTIC.**
	/// q_h is discontinuous across a face, so a zero within the jump belongs to
	/// neither neighbour strictly and each side's polynomial carries its own.
	/// Measured on Soloviev::nstx() at k = 1, n = 16, where the X-point at
	/// r = 0.699700 sits 3.0e-04 from the mesh line r = 0.700000: elements 237
	/// and 238 each hold a root, BOTH with overshoot exactly zero, 6.9e-04
	/// apart. Neither is wrong and "the same root" is not defined there.
	///
	/// It matters to XP-2 rather than to XP-0 because the ELEMENT is what the
	/// border row is built from: a sub-h^(k+1) move of the field can flip which
	/// element is reported, and the outer fixed point then chatters between two
	/// answers instead of converging. The failing case above reads 0.004 here
	/// against 0.16 to 0.20 on a healthy one.
	double elementDepth( meq::CriticalPoint const &point )
	{
		return std::min( { point.referenceX, point.referenceY,
		                   1.0 - point.referenceX - point.referenceY } );
	}

	/// The saddle of q_h that is a genuine X-point rather than the symmetry
	/// axis.
	///
	/// psi vanishes identically on r = 0, so the axis carries near-zeros of q
	/// that sweep() reports as saddles. They are not X-points and no divertor
	/// put them there; anything at an r below a fraction of the plasma's minor
	/// radius is one of those.
	///
	/// Of the rest, the one bounding the plasma is the one carrying the LARGEST
	/// psi, psi having a maximum at the axis here: walking out from the core,
	/// the innermost saddle is met first.
	bool findXPoint( meq::CriticalPointFinder const &finder,
	                 meq::CriticalPoint &found )
	{
		std::vector<meq::CriticalPoint> const points = finder.sweep();
		bool any = false;
		for ( meq::CriticalPoint const &p : points )
		{
			if ( p.type != meq::CriticalPointType::Saddle || p.r < 0.30 )
				continue;
			if ( !any || p.psi > found.psi )
			{
				found = p;
				any = true;
			}
		}
		return any;
	}

	/// psi_h at an arbitrary point, for the one thing the solver does not
	/// expose: whether the LIMITER is inside the plasma or outside it once the
	/// boundary flux has been taken from the X-point instead.
	bool potentialAt( meq::GradShafranovSolver const &solver, double r,
	                  double z, double &value )
	{
		mfem::Mesh &mesh = *solver.potential().FESpace()->GetMesh();
		mfem::DenseMatrix points( 2, 1 );
		points( 0, 0 ) = r;
		points( 1, 0 ) = z;

		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> local;
		if ( mesh.FindPoints( points, elements, local ) < 1
		     || elements[ 0 ] < 0 )
			return false;

		value = solver.potential().GetValue( elements[ 0 ], local[ 0 ] );
		return true;
	}
}

BOOST_AUTO_TEST_CASE( theBoundaryFluxConvergesToTheLocatedXPoint )
{
	/*
	 * THE BOOTSTRAP PIN IS DELIBERATELY NOT THE ANSWER, and 0.05 m is the
	 * distance a machine's own drawings would get you to: it is 8% of the
	 * null's height below the midplane and about two thirds of an element on
	 * this mesh, so the loop has somewhere to travel and the seeded search
	 * still reaches from one sweep to the next.
	 *
	 * Pinning AT the reference X-point would make the first sweep's step the
	 * only thing measured, and it would measure zero.
	 */
	double const bootstrapR = referenceXPointR + 0.05;
	double const bootstrapZ = referenceXPointZ + 0.05;

	Machine m = buildMachine();
	completeMachine( m );

	// THE BOOTSTRAP: one solve pinned at the prior, and the saddle of its own
	// solved q_h.
	BOOST_TEST_REQUIRE( solveAt( m, bootstrapR, bootstrapZ ),
	                    "the diverted machine did not converge pinned at ( "
	                    << bootstrapR << ", " << bootstrapZ << " ), so the "
	                    "outer loop has no starting point" );

	meq::CriticalPoint x;
	{
		meq::CriticalPointFinder finder( *m.solver );
		BOOST_TEST_REQUIRE( findXPoint( finder, x ),
			"the limiter solve of the diverted machine carries no off-axis "
			"saddle, so this configuration is not diverted and XP-2 has nothing "
			"to iterate on" );
	}

	std::printf( "\n  XP-2: psi_bnd FROM THE LOCATED X-POINT, bootstrapped at "
	             "( %.3f, %.3f )\n", bootstrapR, bootstrapZ );
	std::printf( "    %-6s %5s %19s %11s %13s %13s %11s %7s %8s\n",
	             "sweep", "its", "X-point", "step", "psi_ax", "psi_bnd",
	             "psi_X - psi_h", "depth", "route" );
	std::printf( "    %-6s %5zu   (%6.3f,%7.3f) %11s %13.6e %13.6e %11s %7.3f "
	             "%8s\n", "boot", m.solver->newtonResiduals().size() - 1, x.r,
	             x.z, "-", m.solver->psiAxis(), m.solver->psiBoundary(), "-",
	             elementDepth( x ), "SWEEP" );
	std::fflush( stdout );

	struct Sweep
	{
		double r = 0.0;
		double z = 0.0;
		double step = 0.0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		double pinResidual = 0.0;
		double depth = 0.0;
		bool seeded = false;
		std::size_t iterations = 0;
	};
	std::vector<Sweep> sweeps;
	int seededSweeps = 0;

	// THE WARM START LIVES ACROSS THE WHOLE LOOP, the guess being BORROWED and
	// having to outlive the solve it seeds. Seeded from the bootstrap.
	mfem::GridFunction previous( m.solver->potential() );

	bool broke = false;
	for ( int outer = 0; outer < 8; ++outer )
	{
		double const pinnedR = x.r;
		double const pinnedZ = x.z;
		if ( !solveAt( m, pinnedR, pinnedZ, &previous ) )
		{
			broke = true;
			std::printf( "    %-6d %5s %19s\n", outer + 1, "-", "NO SOLVE" );
			break;
		}

		/*
		 * THE SEEDED SADDLE SEARCH, WHICH IS WHAT IT WAS BUILT FOR. The
		 * previous sweep's X-point is a prior worth having: the point moves by
		 * a fraction of an element, and CriticalPointFinder::sweep() costs one
		 * Newton per element where the seeded search costs the rings around
		 * one -- 32x to 205x on the XP-0 fixture, widening with refinement.
		 *
		 * ITS REACH IS ABOUT ONE AND A HALF ELEMENTS, so a step larger than
		 * that DECLINES rather than returning a wrong point, and the sweep is
		 * the fallback. How often that fires is a measurement rather than an
		 * implementation detail, so it is counted and printed.
		 */
		meq::CriticalPointFinder finder( *m.solver );
		meq::CriticalPoint next;
		bool const seeded = finder.tryFindCriticalPointFrom(
			pinnedR, pinnedZ, meq::AxisSense::Saddle, next );
		if ( seeded )
			++seededSweeps;
		if ( !seeded && !findXPoint( finder, next ) )
		{
			broke = true;
			std::printf( "    %-6d %5zu %19s\n", outer + 1,
			             m.solver->newtonResiduals().size() - 1, "NO SADDLE" );
			break;
		}

		Sweep sweep;
		sweep.r = next.r;
		sweep.z = next.z;
		sweep.step = std::hypot( next.r - pinnedR, next.z - pinnedZ );
		sweep.psiAxis = m.solver->psiAxis();
		sweep.psiBoundary = m.solver->psiBoundary();
		// THE FIXED POINT'S OWN DEFINING PROPERTY. The border makes psi_bnd
		// equal psi_h at the PINNED point exactly, whatever that point is; what
		// is not automatic is that the pinned point IS the saddle.
		sweep.pinResidual = next.psi - m.solver->psiBoundary();
		sweep.depth = elementDepth( next );
		sweep.seeded = seeded;
		sweep.iterations = m.solver->newtonResiduals().size() - 1;
		sweeps.push_back( sweep );

		std::printf( "    %-6d %5zu   (%6.3f,%7.3f) %11.3e %13.6e %13.6e "
		             "%11.2e %7.3f %8s\n", outer + 1, sweep.iterations, sweep.r,
		             sweep.z, sweep.step, sweep.psiAxis, sweep.psiBoundary,
		             sweep.pinResidual, sweep.depth,
		             sweep.seeded ? "seeded" : "SWEEP" );
		std::fflush( stdout );

		x = next;
		previous = m.solver->potential();
		if ( sweep.step < 1.0e-10 )
			break;
	}

	BOOST_TEST_REQUIRE( !broke, "the outer loop broke down" );
	BOOST_TEST_REQUIRE( sweeps.size() >= 3,
	                    "the outer loop ran " << sweeps.size() << " sweeps, "
	                    "too few to say anything about contraction" );

	Sweep const &last = sweeps.back();
	double const span = last.psiAxis - last.psiBoundary;

	/*
	 * THE PRECONDITION ON THE FIXTURE, CHECKED RATHER THAN TRUSTED. An X-point
	 * landing on a mesh line is held by two elements at once, each with its own
	 * root and both strictly contained, and which one is reported can flip for
	 * a sub-h^(k+1) move of the field -- so the border row would be built from a
	 * different element between sweeps and the loop would chatter rather than
	 * converge. elementDepth() says at length what was measured.
	 */
	BOOST_TEST( last.depth > 0.05,
		"the converged X-point sits " << last.depth << " into its element in "
		"reference coordinates, so it is on a mesh line and two elements hold "
		"it. That is a statement about where this machine's null falls on this "
		"mesh, not about the solver: move the designed target." );

	// AND THE SEEDED SEARCH IS WHAT THE LOOP RUNS ON, worth asserting rather
	// than only printing -- a loop silently falling back to a full sweep every
	// time is a different cost and would go unnoticed.
	std::printf( "    %d of %zu sweeps reached the X-point from the previous "
	             "one\n", seededSweeps, sweeps.size() );
	BOOST_TEST( seededSweeps + 1 >= static_cast<int>( sweeps.size() ),
		"only " << seededSweeps << " of " << sweeps.size() << " sweeps reached "
		"the X-point from the previous one. The seeded search reaches about one "
		"and a half elements, so repeated fallback is a statement about the "
		"outer map rather than about the search." );

	// ONE: IT IS A FIXED POINT, and it stops well inside an element -- a loop
	// merely wandering within one would look converged on a coarser measure.
	BOOST_TEST( last.step < 1.0e-8,
		"the outer iteration left the X-point moving by " << last.step
		<< " metres on its last sweep, so psi_bnd is whatever the loop happened "
		"to stop at rather than the flux at the saddle." );

	// TWO: IT CONTRACTS.
	BOOST_TEST( last.step < 1.0e-3*sweeps.front().step,
		"the outer step went from " << sweeps.front().step << " to "
		<< last.step << " over " << sweeps.size() << " sweeps. XP-2 is a "
		"fixed-point iteration and a slow one is a different algorithm from a "
		"fast one -- section 10.4's three-row border is the answer if this is "
		"what the map looks like." );

	// THREE: psi_bnd IS THE FLUX AT THE SADDLE.
	BOOST_TEST( std::abs( last.pinResidual ) < 1.0e-10*span,
		"psi at the located saddle differs from the psi_bnd the solve returned "
		"by " << last.pinResidual << ", which is "
		<< std::abs( last.pinResidual )/span << " of the span. The loop has "
		"converged to a point that is not the saddle." );

	/*
	 * FOUR: THE PLASMA IS DIVERTED AND THE NULL XP-2 FOUND IS THE ACTIVE ONE.
	 *
	 * This machine has TWO off-axis saddles -- freegs4e puts the lower at
	 * psi = 3.240412550738516e-02 and the upper at 2.891018784455272e-02, and
	 * that 3.5e-03 gap is the whole of what makes it a SINGLE null. The
	 * bounding surface is the null of LARGEST psi, psi having its maximum on
	 * the axis, so a loop that converged on the upper one would satisfy every
	 * assertion above and describe a different plasma.
	 *
	 * So the case sweeps the converged field for every off-axis saddle and
	 * asserts both halves: that there is more than one, which is the fixture's
	 * own premise, and that the one the loop stopped at carries the largest
	 * psi.
	 */
	std::vector<meq::CriticalPoint> nulls;
	{
		meq::CriticalPointFinder finder( *m.solver );
		for ( meq::CriticalPoint const &p : finder.sweep() )
			if ( p.type == meq::CriticalPointType::Saddle && p.r > 0.30 )
				nulls.push_back( p );
	}
	std::printf( "    %zu off-axis saddles in the converged field:", nulls.size() );
	for ( meq::CriticalPoint const &p : nulls )
		std::printf( "  ( %6.3f, %7.3f ) psi %11.4e", p.r, p.z, p.psi );
	std::printf( "\n" );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( nulls.size() >= 2,
		"the converged field carries " << nulls.size() << " off-axis saddle( s ). "
		"This machine is an up-down asymmetric DOUBLE null and the reference "
		"finds both, so one is a statement about the discretisation or about "
		"which branch the solve reached, and it makes the active-null check "
		"below vacuous." );

	double largest = -std::numeric_limits<double>::infinity();
	for ( meq::CriticalPoint const &p : nulls )
		largest = std::max( largest, p.psi );
	BOOST_TEST( std::abs( largest - last.psiBoundary ) < 1.0e-8*span,
		"the loop settled on a null carrying psi = " << last.psiBoundary
		<< " where the field's largest off-axis saddle carries " << largest
		<< ". The bounding surface is a MIN over candidates -- the null met "
		"first walking out from the axis -- and this one is not it." );

	/*
	 * AND AGAINST THE INDEPENDENT CODE, WHICH IS THE ONLY ASSERTION HERE THAT
	 * THE ANSWER IS RIGHT RATHER THAN SELF-CONSISTENT.
	 *
	 * The tolerance is per cent and the reason is recorded in
	 * examples/diverted-tokamak.toml: fgsref.py splines its own analytic
	 * profile before solving and the fit MOVES it -- 2.514e-05 of the amplitude
	 * in p' and 1.707e-02 in ff' -- while MEQ's tables are the analytic shape.
	 * So this is a cross-check at the level the two problem statements agree,
	 * not a convergence measurement.
	 */
	double const nullGap = std::hypot( last.r - referenceXPointR,
	                                   last.z - referenceXPointZ );
	std::printf( "    against freegs4e: X-point ( %.6f, %.6f ) vs ( %.6f, "
	             "%.6f ), %.3e m apart;  psi_bnd %.6e vs %.6e;  psi_ax %.6e vs "
	             "%.6e\n", last.r, last.z, referenceXPointR, referenceXPointZ,
	             nullGap, last.psiBoundary, referencePsiBoundary, last.psiAxis,
	             referencePsiAxis );
	std::fflush( stdout );

	BOOST_TEST( nullGap < 2.0e-2,
		"the located X-point is " << nullGap << " m from freegs4e's ( "
		<< referenceXPointR << ", " << referenceXPointZ << " )." );
	BOOST_TEST( std::abs( last.psiBoundary - referencePsiBoundary )
	            < 2.0e-2*std::abs( referencePsiBoundary ),
		"psi_bnd is " << last.psiBoundary << " against freegs4e's "
		<< referencePsiBoundary );
	BOOST_TEST( std::abs( last.psiAxis - referencePsiAxis )
	            < 2.0e-2*std::abs( referencePsiAxis ),
		"psi_ax is " << last.psiAxis << " against freegs4e's "
		<< referencePsiAxis );

	/*
	 * AND THE CORE IS WHERE THE REFERENCE PUTS IT, WHICH IS A DIFFERENT CLAIM
	 * FROM ANY OF THE ABOVE AND IS THE ONE THAT CATCHES A WRONG BRANCH.
	 *
	 * psi_ax, psi_bnd and I_p are all border unknowns, so a solve that has run
	 * away to a different equilibrium still reports three plausible numbers and
	 * a converged residual. What it cannot fake is psi_h AT the reference's own
	 * magnetic axis: on the branch this machine reached before XP-2 -- axis
	 * driven out to ( 2.2626, 0.5413 ), against Gamma -- psi at ( 1.3513,
	 * 0.0622 ) read 1.32e-02 where the reference carries 8.27e-02, while
	 * psi_ax itself read 8.53e-02 and looked healthy.
	 */
	double corePsi = 0.0;
	BOOST_TEST_REQUIRE( potentialAt( *m.solver, referenceAxisR, referenceAxisZ,
	                                 corePsi ),
	                    "the reference's magnetic axis is not in the mesh" );
	std::printf( "    psi_h at the reference axis ( %.4f, %.4f ) = %.6e "
	             "against psi_ax = %.6e\n", referenceAxisR, referenceAxisZ,
	             corePsi, last.psiAxis );
	std::fflush( stdout );

	BOOST_TEST( std::abs( corePsi - last.psiAxis ) < 5.0e-2*span,
		"psi_h at the reference equilibrium's magnetic axis is " << corePsi
		<< " where this solve's psi_ax is " << last.psiAxis << ", "
		<< std::abs( corePsi - last.psiAxis )/span << " of the span apart. The "
		"core of this solve is not where the reference's is, so the solve has "
		"converged to a different equilibrium and every border unknown above "
		"describes it rather than this machine." );

	// FIVE: AND THE EQUILIBRIUM IS THE ONE THE SHIPPED MACHINE CASE CALLS
	// HEALTHY, by the same three checks.
	meq::GradShafranovSolver::AxisSourceCheck const axisSource =
		m.solver->checkAxisSource();
	meq::CriticalPointFinder finder( *m.solver );
	meq::AxisAgreement const axis =
		finder.checkAxis( m.solver->psiAxis(), m.solver->psiBoundary() );

	std::printf( "    | F | on r = 0 is %.3e, psi_ax attained at ( %.3f, %.3f ),"
	             " Psi at the O-point %.4f over %d extrema and %d saddles, "
	             "mu0 I_p delivered %.6e\n", axisSource.worstOnAxis, axis.nodeR,
	             axis.nodeZ, axis.normalisedFlux, axis.extrema, axis.saddles,
	             m.solver->plasmaCurrent() );
	std::fflush( stdout );

	BOOST_TEST( axisSource.bounded,
		"| F | on the symmetry axis is " << axisSource.worstOnAxis
		<< ", so F/r = mu_0 j_phi is an unbounded toroidal current density on "
		"r = 0." );
	BOOST_TEST( axis.nodeR > 0.30,
		"psi_ax is attained at r = " << axis.nodeR << ", which is on or beside "
		"the symmetry axis -- section 11.3's axis layer." );
	BOOST_TEST( axis.agrees,
		"the O-point of q_h carries Psi = " << axis.normalisedFlux
		<< " against the 1 it must carry by definition." );

	double const wanted = m.mu0*m.config->getSource().getMHD().plasmaCurrent;
	BOOST_TEST( std::abs( m.solver->plasmaCurrent() - wanted )
	            < 1.0e-5*std::abs( wanted ),
		"the delivered current is " << m.solver->plasmaCurrent()
		<< " against the " << wanted << " prescribed" );
}
