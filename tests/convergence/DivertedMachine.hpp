#ifndef MEQ_TESTS_DIVERTED_MACHINE_HPP
#define MEQ_TESTS_DIVERTED_MACHINE_HPP

/*
 * THE DIVERTED MACHINE, ASSEMBLED ONCE FOR EVERY CASE THAT NEEDS ONE.
 * FREE-BOUNDARY-PLAN.md section 10.
 *
 * examples/diverted-tokamak.toml is freegs4e's own A_testtokamak_classic: an
 * up-down asymmetric double null whose coil currents were SOLVED FOR by its
 * control system, which section 7.15 argues is what makes a free-boundary
 * Newton converge at all. It is BORROWED rather than designed, and the header
 * of XPointOuter.cpp records the twelve hand-built single nulls that failed
 * before it and why -- the short form being that this tree's own geometry
 * cannot carry a null without driving psi_bnd negative, which puts the
 * symmetry axis inside the plasma support and grows section 11.3's axis layer.
 *
 * IT LIVES IN A HEADER BECAUSE XP-2 AND XP-3 SOLVE THE SAME MACHINE AND MUST
 * SOLVE THE SAME ONE. XP-3's acceptance is agreement with XP-2, so a second
 * copy of this assembly would make a disagreement between them ambiguous
 * between the two borders and the two fixtures -- which is the one thing that
 * comparison exists to rule out.
 *
 * The sequence follows apps/meq.cpp, which is the only place this assembly
 * exists; a fixture that assembled it differently would be measuring its own
 * wiring.
 *
 * Include it AFTER <boost/test/unit_test.hpp>: the builder reports a broken
 * fixture through BOOST_TEST_REQUIRE, so that a missing mesh or an empty
 * subdomain names itself rather than surfacing as a solver failure.
 */

#include <algorithm>
#include <cmath>
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

namespace meqtest
{
	/// The diverted machine, taken whole: mesh, profiles, coil currents,
	/// limiter contact, exterior expansion and guess. Nothing here is designed
	/// or tuned -- see the header for why that matters.
	inline char const *const machineFile = "examples/diverted-tokamak.toml";

	/// freegs4e's own answer for this machine, from
	/// tools/freegs4e-benchmark/A_testtokamak_classic.json. The ACTIVE null is
	/// the lower one; the upper sits 3.5e-03 further out.
	inline double const referenceXPointR = 1.093144118182931;
	inline double const referenceXPointZ = -0.6039650838688502;
	inline double const referencePsiBoundary = 3.240412550738516e-02;
	inline double const referencePsiAxis = 8.271751444840184e-02;
	inline double const referenceAxisR = 1.351273;
	inline double const referenceAxisZ = 0.062226;

	/// Negative inside. The semicircle about the axis, which is what makes the
	/// exterior expansion legal at all; the flat side IS r = 0 and stays
	/// ordinary fitted boundary.
	inline mfem::PositionFunction semicircle( double radius, double centreZ )
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
	/// @param refinements uniform refinements of the gmsh mesh before anything
	///        else happens, for a case that needs the SAME machine at more than
	///        one resolution. Zero is the shipped mesh and is what every case
	///        that does not ask gets. Refining here rather than re-running
	///        `halfdisc.py` at a smaller `--size` is deliberate: a re-meshed
	///        geometry moves which elements the conductors and the limiter
	///        region occupy, which is a second variable in a study about the
	///        first. Uniform refinement subdivides and carries the material
	///        attributes down, so the machine is unchanged and only `h` moves.
	inline Machine buildMachine( int refinements = 0 )
	{
		Machine m;
		m.config = std::make_unique<meq::Configuration>( machineFile );

		meq::SourceConfig const &sourceConfig = m.config->getSource();
		m.mu0 = sourceConfig.permeability();

		m.background = std::make_unique<mfem::Mesh>(
			m.config->getMesh().file.c_str(), 1, 1 );
		for ( int r = 0; r < refinements; ++r )
			m.background->UniformRefinement();

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

	inline void completeMachine( Machine &m )
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
		m.source = std::make_shared<meq::CoilAugmentedNormalisedSource>(
			m.plasma, m.coils );

		/*
		 * ON THE WRAPPER AND AFTER IT EXISTS, WHICH IS WHAT apps/meq.cpp DOES
		 * AND WHAT THIS FIXTURE DID NOT.
		 *
		 * setPlasmaSupport() is virtual so that a call on the WRAPPER reaches
		 * the source that evaluates the profiles; forwarding runs that way and
		 * only that way. Setting it on the inner source before wrapping leaves
		 * the wrapper's own flag false -- and the wrapper is what the solver
		 * holds, so GradShafranovSolver::plasmaComponentWanted() reads false and
		 * XP-1's flood fill never runs. F is still confined, because the inner
		 * source's pointwise test is live, so nothing fails loudly.
		 *
		 * WHAT THAT COSTS IS THE WHOLE OF XP-1 ON THE ONE DIVERTED CASE IN THIS
		 * TREE. Section 10.3: a diverted plasma's level set is DISCONNECTED
		 * across the X-point, so the pointwise test picks up the private flux
		 * region and gives it a current channel nobody asked for. That is the
		 * configuration XP-1 was built for and this fixture was running without
		 * it -- the fill reporting 0 of 0 elements over 0 components, which is
		 * what "no fill is live" prints and not what an empty plasma prints.
		 *
		 * The conductors stay outside the support either way: the wrapper's
		 * override confines the PLASMA term only, and answers fOutsidePlasma()
		 * with the coil term exactly, so a coil in the vacuum is untouched.
		 */
		m.source->setPlasmaSupport( true );

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
	inline bool solveAt( Machine &m, double r, double z,
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
	inline double elementDepth( meq::CriticalPoint const &point )
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
	inline bool findXPoint( meq::CriticalPointFinder const &finder,
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
	inline bool valueAt( mfem::GridFunction const &field, double r, double z,
	              double &value )
	{
		mfem::Mesh &mesh = *field.FESpace()->GetMesh();
		mfem::DenseMatrix points( 2, 1 );
		points( 0, 0 ) = r;
		points( 1, 0 ) = z;

		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> local;
		if ( mesh.FindPoints( points, elements, local ) < 1
		     || elements[ 0 ] < 0 )
			return false;

		value = field.GetValue( elements[ 0 ], local[ 0 ] );
		return true;
	}

	inline bool potentialAt( meq::GradShafranovSolver const &solver, double r,
	                  double z, double &value )
	{
		return valueAt( solver.potential(), r, z, value );
	}

	/// **FIX THE SUPPORT FOR THE NEXT SOLVE AT @a state, AND RE-DECIDE IT HERE
	/// RATHER THAN INSIDE NEWTON.** Section 10.5's prescription, applied to the
	/// support as well as to the bounding point: the threshold insidePlasma()
	/// tests against and the connected component the fill reaches are both held
	/// for the whole of one solve and both moved between solves.
	///
	/// BOTH HALVES, because either alone leaves the support moving. The
	/// threshold is the source's -- meq::NormalisedSource::freezePlasmaEdge --
	/// and the component is the solver's; the mask is element granular, so
	/// freezing it while the pointwise test drifts still moves the edge inside
	/// every element the fill reached.
	inline void freezeSupportAt( Machine &m, mfem::GridFunction const &state,
	                      double axis, double boundary )
	{
		m.source->freezePlasmaEdge( axis, boundary );
		// The PUBLIC refresh, which setPlasmaSupportFrozen() deliberately does
		// not suppress: this is the outer loop moving the support, and it is
		// the only thing that may.
		m.solver->setPlasmaSupportFrozen( false );
		m.solver->refreshPlasmaComponent( state );
		m.solver->setPlasmaSupportFrozen( true );
	}
}

#endif
