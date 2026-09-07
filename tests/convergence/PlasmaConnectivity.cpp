/*
 * THE PLASMA IS A CONNECTED SET AND `{ Psi > 0 }` IS NOT -- stage XP-1 of
 * FREE-BOUNDARY-PLAN.md section 10.6.
 *
 * meq::NormalisedSource::insidePlasma() is a POINTWISE test on the value,
 * ( psi - psi_bnd )*span > 0, with no connectivity in it whatever. Section 10.3
 * records what that costs on a DIVERTED plasma -- across an X-point the level
 * psi = psi_X cuts a neighbourhood into four sectors and two opposite ones
 * carry Psi > 0, so the source switches on in the private flux region as well
 * as in the plasma -- and calls it latent, on the grounds that no shipped
 * example sets ConfineToPlasma and MEQ has no diverted case.
 *
 * IT IS NOT LATENT, AND THE HALF-DISC CASE BELOW IS WHY. A LIMITER
 * configuration that ships and converges -- FreeBoundaryCoupling's own
 * half-disc, two borders, an exterior coupling, NO X-POINT ANYWHERE -- already
 * has a `{ psi > psi_bnd }` with more than one component. A level of a field
 * that is not monotone in radius cuts the domain into as many lobes as it
 * likes, and a saddle is only the most dramatic way to get one.
 *
 * WHAT THIS FILE MEASURES, IN THE ORDER THE ANSWERS DEPEND ON EACH OTHER:
 *
 *   1. theIntegratorHonoursTheElementItIsGiven -- the plumbing, directly.
 *      meq::SourceIntegrator reads ElementTransformation::ElementNo and
 *      nothing else identifies an element, so if that number were wrong the
 *      mask would apply to arbitrary elements and every measurement below
 *      would be noise. Checked against a mask holding ONE element.
 *
 *   2. theFillSeparatesThePrivateFluxRegionFromThePlasma -- section 10.6's
 *      sharp question, on Soloviev::iterExample2(), whose X-point is known in
 *      closed form. Does a FACE-neighbour fill need the explicit X-point
 *      blocking freegs4e's grid-based core_mask requires? IT DOES: the answer
 *      is no on section 10.3's reasoning and yes in fact, because two lobes
 *      meeting at a point are joined by the BAND of elements straddling the
 *      separatrix, which near a saddle is several elements wide. The two-rule
 *      fill is what fixes it and it locates nothing. Measured as an element
 *      count AND as the integral of |F|, because a count alone does not say
 *      the difference matters.
 *
 *   3. theLimiterCaseAlreadyHasMoreThanOneLobe -- the fixture above, solved,
 *      with the components counted at the converged answer.
 *
 *   4. theConnectivityChangeIsNotAJumpAtTheSizeThatMatters -- what it costs
 *      Newton. FB-4's discriminator was that the largest step between
 *      neighbouring samples DOES NOT SHRINK as the sampling is refined, and it
 *      is applied here to a connectivity change instead of a quadrature rule.
 *      IT SHRINKS: the control falls by 3.99 and 4.00 as the sampling is
 *      quartered and the connected column by 3.94 and 3.96, with the worst step
 *      at 0.0017% of int |F|. What changes hands as the level slides is a
 *      STRADDLING element, where Psi is near zero and the profile with it -- so
 *      the O( 1 ) jump the structure predicts is not exercised at a size this
 *      experiment can see, and the fill may live inside the Newton loop. The
 *      ring-depth rule this replaced DID jump, 1.68 then 1.20, which is what
 *      says the property belongs to the rule rather than to the problem.
 *
 *   5. theFillDoesNotMoveASingleLobeAnswer -- the control. An UNCONFINED solve
 *      must not see the fill at all, bit for bit; a CONFINED solve whose
 *      support is connected must reach the same answer to the tolerance it was
 *      given, in no more steps. Without this every case above is compatible
 *      with a fill that quietly perturbs the answer.
 *
 * AND ONE THING WAS FOUND THAT IS NOT ABOUT CONNECTIVITY AT ALL, in case 3:
 * psi_ax on that shipped converged answer is attained in an element TOUCHING
 * r = 0, and reads 2.5x the largest psi_h anywhere off the symmetry axis. That
 * is where the fill's seed rule comes from, and it is worth knowing about the
 * fixture independently of XP-1.
 */

#define BOOST_TEST_MODULE PlasmaConnectivity
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PlasmaComponent.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"

using meq::GradShafranovSolver;
using meq::PlasmaComponent;
using meq::analytic::SolovievEquilibrium;

namespace
{
	/// p'( Psi ) = amplitude * Psi^power, exact at all three derivative levels.
	///
	/// Copied from FreeBoundaryCoupling.cpp rather than shared, so that a change
	/// to either file's fixture cannot silently move the other's numbers. The
	/// power is a PRECONDITION and not a knob: FB-4 measured that p'( 0 ) != 0
	/// makes the residual discontinuous and Newton converges from nowhere.
	class PowerProfile : public meq::Profile
	{
		public:
			PowerProfile( double amplitudeIn, int powerIn )
				: amplitude( amplitudeIn ), power( powerIn ) {}

			double operator()( double psi ) const override
			{
				return amplitude*std::pow( psi, power );
			}
			double prime( double psi ) const override
			{
				return power < 1 ? 0.0
				       : amplitude*power*std::pow( psi, power - 1 );
			}
			double doublePrime( double psi ) const override
			{
				return power < 2 ? 0.0
				       : amplitude*power*( power - 1 )*std::pow( psi, power - 2 );
			}

		private:
			double amplitude;
			int power;
	};

	/// The coil term alone, as a meq::Source, so that what survives on a masked
	/// element can be compared against a closed form rather than against the
	/// same expression computed the same way one line up.
	class CoilOnly : public meq::Source
	{
		public:
			explicit CoilOnly( std::shared_ptr<meq::CoilSet const> coilsIn )
				: coils( std::move( coilsIn ) ) {}

			double f( double r, double z, double ) const override
			{
				return coils->f( r, z );
			}
			double dFdPsi( double, double, double ) const override
			{
				return 0.0;
			}

		private:
			std::shared_ptr<meq::CoilSet const> coils;
	};

	/// The half-disc of FreeBoundaryCoupling.cpp, transcribed rather than shared.
	///
	/// Gamma is the semicircle rho = 1.5 about the AXIS, cut out of a background
	/// box that reaches r = 0 exactly, with Gamma_h the generated arc and the
	/// flat side inherited fitted boundary -- the axis is not an approximation of
	/// anything and needs no transfer.
	///
	/// **THE BACKGROUND IS STATIC ON PURPOSE.** mfem::SubMesh keeps a pointer to
	/// its parent and mfem::VertexConePath reads the parent's edges, so a
	/// background built as a local and left behind is a dangling pointer that
	/// crashes inside MFEM with no MEQ frame in the trace -- CLAUDE.md records
	/// three fixtures that got exactly this wrong.
	double const halfDiscGamma = 1.5;
	double const halfDiscBox = 1.7;

	double halfDiscLevelSet( mfem::Vector const &x )
	{
		return std::hypot( x( 0 ), x( 1 ) ) - halfDiscGamma;
	}

	struct HalfDisc
	{
		mfem::SubMesh *sub = nullptr;
		mfem::VertexConePath *path = nullptr;
		mfem::Array<int> gammaHMarker;
		double h = 0.0;
	};

	HalfDisc makeHalfDisc( int n )
	{
		static std::vector<std::unique_ptr<mfem::Mesh>> backgrounds;
		static std::vector<std::unique_ptr<mfem::SubMesh>> subs;
		static std::vector<std::unique_ptr<mfem::VertexConePath>> paths;

		auto background = std::make_unique<mfem::Mesh>(
			mfem::Mesh::MakeCartesian2D( n, 2*n, mfem::Element::TRIANGLE, false,
			                             halfDiscBox, 2.0*halfDiscBox ) );
		background->Transform( [ ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 );
			out( 1 ) = in( 1 ) - halfDiscBox;
		} );

		HalfDisc out;
		out.h = halfDiscBox/static_cast<double>( n );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE( mfem::MarkLevelSetSubdomain(
			*background, halfDiscLevelSet, 0.0, marker, 1 ) > 0,
			"the half-disc is empty at n = " << n );
		for ( int e = 0; e < background->GetNE(); ++e )
			background->SetAttribute( e, marker[ e ] ? 1 : 2 );
		background->SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		auto sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( *background, domainAttr ) );

		int const gammaH = sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( sub->bdr_attributes.Size() >= 2,
			"D_h has only one boundary attribute at n = " << n
			<< ", so the axis was not inherited and Gamma_h has swallowed it" );

		out.gammaHMarker.SetSize( gammaH );
		out.gammaHMarker = 0;
		out.gammaHMarker[ gammaH - 1 ] = 1;

		auto path = std::make_unique<mfem::VertexConePath>(
			*sub, gammaH, halfDiscLevelSet, 6.0*out.h );

		out.sub = sub.get();
		out.path = path.get();

		backgrounds.push_back( std::move( background ) );
		subs.push_back( std::move( sub ) );
		paths.push_back( std::move( path ) );
		return out;
	}

	/// The midplane components of `{ psi > psi_bnd }`, by sampling: what a
	/// pointwise support test would switch the source on in.
	///
	/// A ONE-DIMENSIONAL count, and it is a lower bound on the two-dimensional
	/// one rather than the same number -- a lobe that misses the midplane is
	/// invisible to it. It is here because it is the cheapest possible statement
	/// of the finding and it needs no mesh: run a ray across the domain and
	/// count sign changes.
	int midplaneComponents( GradShafranovSolver &solver, double psiBnd,
	                        double span, double rMin, double rMax, int samples )
	{
		mfem::GridFunction const &potential = solver.potential();
		mfem::Mesh *mesh = potential.FESpace()->GetMesh();

		int components = 0;
		bool previous = false;
		for ( int i = 0; i <= samples; ++i )
		{
			double const r = rMin + ( rMax - rMin )*i/samples;

			mfem::DenseMatrix point( 2, 1 );
			point( 0, 0 ) = r;
			point( 1, 0 ) = 0.0;
			mfem::Array<int> elements;
			mfem::Array<mfem::IntegrationPoint> reference;
			mesh->FindPoints( point, elements, reference );

			bool now = false;
			if ( elements[ 0 ] >= 0 )
			{
				double const psi = potential.GetValue( elements[ 0 ], reference[ 0 ] );
				now = ( psi - psiBnd )*span > 0.0;
			}
			if ( now && !previous )
				++components;
			previous = now;
		}
		return components;
	}
}


/*
 * THE PLUMBING, AND IT RESTS ENTIRELY ON ONE INTEGER.
 *
 * meq::SourceIntegrator applies the mask by ElementTransformation::ElementNo.
 * Nothing else in an AssembleElementVector call identifies which element it is,
 * so if that number were not the element -- if MFEM's own workspaces left it
 * stale, say -- the mask would switch elements off at random and every
 * measurement in this file would be noise that happened to look like a result.
 *
 * So it is checked directly rather than inferred: a component holding exactly
 * ONE element, assembled on that element and on a neighbour, must give the full
 * source term on the first and fOutsidePlasma() on the second.
 *
 * AND THE COIL WRAPPER IS THE SECOND HALF, because the trap is the one this
 * tree has now met four times: a wrapper that forgets to forward. A coil sits
 * in the vacuum region by construction, so an element the fill did not reach
 * must still carry its current -- and the base class answers zero.
 */
BOOST_AUTO_TEST_CASE( theIntegratorHonoursTheElementItIsGiven )
{
	int const order = 2;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	for ( int v = 0; v < mesh.GetNV(); ++v )
		mesh.GetVertex( v )[ 0 ] += 1.0;          // r in [ 1, 2 ], off the axis.

	mfem::L2_FECollection collection( order, 2, mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace space( &mesh, &collection );

	auto pPrime = std::make_shared<PowerProfile const>( 0.6, 1 );
	auto ggPrime = std::make_shared<PowerProfile const>( 0.05, 1 );
	auto plasma = std::make_shared<meq::NormalisedMHDSource>( pPrime, ggPrime,
	                                                          0.1, 1.0 );
	plasma->setNormalisation( 0.1, 0.0 );

	// A coil well inside the box, so that every element sees some of it.
	auto coils = std::make_shared<meq::CoilSet>();
	coils->add( meq::Coil( 1.5, 0.5, 0.6, 0.6, 1.0e5 ) );
	meq::CoilAugmentedNormalisedSource augmented( plasma, coils );

	// A component holding element 0 alone: adjacency with no edges at all, so
	// the fill cannot leave the seed.
	PlasmaComponent single;
	std::vector<int> offsets( static_cast<std::size_t>( mesh.GetNE() ) + 1, 0 );
	single.setAdjacency( offsets, {} );
	std::vector<char> everywhere( static_cast<std::size_t>( mesh.GetNE() ), 1 );
	single.fill( everywhere, 0 );
	BOOST_TEST( single.componentNodes() == 1 );

	auto assemble = [ & ]( meq::Source const &source, PlasmaComponent const *mask,
	                       int element )
	{
		meq::SourceIntegrator integrator( source, 4 );
		integrator.setPlasmaComponent( mask );

		mfem::FiniteElement const &el = *space.GetFE( element );
		mfem::IsoparametricTransformation transformation;
		mesh.GetElementTransformation( element, &transformation );

		mfem::Array<int> dofs;
		space.GetElementDofs( element, dofs );
		mfem::Vector coefficients( dofs.Size() );
		// A constant psi of 0.05, which is inside the plasma at psi_ax = 0.1.
		coefficients = 0.05;

		mfem::Vector out;
		integrator.AssembleElementVector( el, transformation, coefficients, out );
		return out.Norml2();
	};

	std::printf( "\n  THE MASK IS APPLIED BY ElementTransformation::ElementNo\n" );
	std::printf( "    %-34s %14s %14s\n", "source", "element 0", "element 1" );

	double const plainSeed = assemble( *plasma, &single, 0 );
	double const plainOther = assemble( *plasma, &single, 1 );
	double const coilSeed = assemble( augmented, &single, 0 );
	double const coilOther = assemble( augmented, &single, 1 );
	double const unmaskedOther = assemble( *plasma, nullptr, 1 );
	CoilOnly const coilSource( coils );
	double const coilOnly = assemble( coilSource, nullptr, 1 );

	std::printf( "    %-34s %14.6e %14.6e\n", "plasma alone, masked",
	             plainSeed, plainOther );
	std::printf( "    %-34s %14.6e %14.6e\n", "plasma + coils, masked",
	             coilSeed, coilOther );
	std::printf( "    %-34s %14.6e %14.6e\n", "plasma alone, no mask",
	             plainSeed, unmaskedOther );
	std::printf( "    the coil term alone on element 1 reads %.6e\n", coilOnly );
	std::fflush( stdout );

	// THE SEED IS UNTOUCHED. If ElementNo were wrong this is the assertion that
	// would most likely still pass, which is why it is not the only one.
	BOOST_TEST( plainSeed > 0.0 );
	BOOST_TEST( std::abs( plainSeed - assemble( *plasma, nullptr, 0 ) )
	            <= 0.0,
	            "the mask changed the element it HOLDS, so it is not being read "
	            "as a per-element test at all" );

	// AND EVERY OTHER ELEMENT IS OFF. Exactly zero, not small: fOutsidePlasma()
	// returns a literal 0.0 for a source with no coils.
	BOOST_TEST( plainOther == 0.0,
	            "element 1 assembled " << plainOther << " with a mask that holds "
	            "only element 0. Either ElementNo is not the element or the mask "
	            "is not being consulted" );
	BOOST_TEST( unmaskedOther > 0.0,
	            "element 1 assembles nothing even WITHOUT a mask, so the case "
	            "above is vacuous" );

	// THE COIL SURVIVES THE MASK, which is the forwarding obligation. Without
	// CoilAugmentedNormalisedSource::fOutsidePlasma() this reads exactly zero
	// and every conductor in the machine is switched off outside the plasma.
	BOOST_TEST( coilOther > 0.0,
	            "the coil term vanished on a masked element, so "
	            "CoilAugmentedNormalisedSource::fOutsidePlasma() is not "
	            "forwarding -- the same trap setPlasmaSupport() and "
	            "normalisationDerivatives() already document" );
	BOOST_TEST( std::abs( coilOther - coilOnly ) <= 1.0e-14*coilOnly,
	            "the masked element carries " << coilOther << " where the coil "
	            "term alone is " << coilOnly << ": what survives outside the "
	            "plasma must be the coils and nothing else" );
}


/*
 * SECTION 10.6's SHARP QUESTION: DOES A FACE-NEIGHBOUR FILL LEAK THROUGH AN
 * X-POINT?
 *
 * `../freegs4e` fills on a uniform ( R, z ) grid and has to block a
 * neighbourhood of each X-point first, because on a grid the two lobes -- which
 * meet at a POINT -- are eight-neighbours. Section 10.3 predicts MEQ needs no
 * such blocking, "except in whichever element actually contains the X-point,
 * which is cut by both branches and is a face neighbour of both lobes. So the
 * expected leak is ONE ELEMENT WIDE".
 *
 * The fixture is Soloviev::iterExample2(), whose saddle sits at exactly
 * ( 0.88384, -0.704 ) with psi = 3.3e-52 -- so psi_bnd = 0 is the separatrix
 * value in closed form and there is nothing to locate. Its F is single-signed
 * NEGATIVE, so psi is a subsolution, the magnetic axis is an interior MINIMUM
 * and the span psi_ax - psi_bnd is negative; the product test in
 * NormalisedSource::insidePlasma() is written to be sign-agnostic for exactly
 * this reason and this fixture is what exercises it.
 *
 * MEASURED AS AN ELEMENT COUNT AND AS THE INTEGRAL OF |F|, because section 10.6
 * asks for both: a count alone does not say whether the difference matters, and
 * the private flux region is where the source would be switched on by mistake.
 */
BOOST_AUTO_TEST_CASE( theFillSeparatesThePrivateFluxRegionFromThePlasma )
{
	SolovievEquilibrium const equilibrium = SolovievEquilibrium::iterExample2();

	// The saddle, in closed form. Soloviev.hpp records it as exactly
	// x_sep = 1 - 1.1 delta eps and y_sep = -1.1 kappa eps.
	double const xPointR = 0.88384;
	double const xPointZ = -0.704;

	/// Everything one mesh has to say about the question.
	struct Reading
	{
		int elements = 0;
		int candidates = 0;
		int xPointElement = -1;
		double largestPsiN = 0.0;

		int pointwiseBelow = 0;
		double pointwiseMass = 0.0;

		int oneRuleBelow = 0;
		int blockedBelow = 0;
		double blockedMass = 0.0;

		int interiorComponents = 0;
		int twoRuleKept = 0;
		int twoRuleRing = 0;
		int twoRuleBelow = 0;
		double twoRuleMass = 0.0;
		int twoRuleLostAbove = 0;

		// The interior rule alone, with no band shared out: the control that
		// says what the watershed is buying.
		int interiorOnlyBelow = 0;
		int interiorOnlyLost = 0;

		int twoRuleLostBelow = 0;
	};

	auto measure = [ & ]( int n )
	{
		Reading out;

		// A box holding the plasma ( r in [ 0.68, 1.32 ], up to z = 0.64 ) and a
		// generous slab of the private flux region below the X-point.
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, 2*n, mfem::Element::TRIANGLE, false, 0.9, 1.8 );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			mesh.GetVertex( v )[ 0 ] += 0.55;
			mesh.GetVertex( v )[ 1 ] -= 1.10;
		}
		out.elements = mesh.GetNE();

		/*
		 * psi_ax, THE MAGNETIC AXIS, AND IT MUST BE SOUGHT INSIDE THE PLASMA
		 * RATHER THAN OVER THE WHOLE BOX -- WHICH IS A FINDING AND NOT A FIXTURE
		 * DETAIL.
		 *
		 * This equilibrium's F is single-signed NEGATIVE, so psi is a
		 * subsolution, the magnetic axis is an interior MINIMUM and psi = 0 is
		 * the separatrix -- CriticalPoints.hpp records the same sign trap for
		 * findAxis(), which seeds from BOTH nodal extremes for exactly this
		 * reason. Measured on this box, the most negative psi is -1.097e-01 at
		 * ( 0.879, -1.100 ), which is the bottom edge and is DEEP IN THE PRIVATE
		 * FLUX REGION; the magnetic axis is -3.930e-02 at ( 1.050, 0.025 ).
		 *
		 * So an argmax of Psi over the whole domain -- which is what
		 * GradShafranovSolver::refreshPlasmaComponent uses, and what psi_ax's own
		 * border constrains -- would seed the fill in the WRONG LOBE and name the
		 * private flux region as the plasma. It cannot fire on a free-boundary
		 * solve, where Gamma bounds the domain and the private flux inside it is
		 * a pocket rather than an unbounded sector, but the exposure is real.
		 */
		double psiAxis = 0.0;
		int axisElement = -1;
		for ( int e = 0; e < out.elements; ++e )
		{
			mfem::Vector centre;
			mesh.GetElementCenter( e, centre );
			double const psi = equilibrium.psi( centre( 0 ), centre( 1 ) );
			if ( std::abs( centre( 1 ) ) < 0.4 && psi < psiAxis )
			{
				psiAxis = psi;
				axisElement = e;
			}
		}
		double const psiBnd = 0.0;
		double const span = psiAxis - psiBnd;
		BOOST_TEST_REQUIRE( psiAxis < 0.0 );
		BOOST_TEST_REQUIRE( axisElement >= 0 );

		// The two rules GradShafranovSolver::refreshPlasmaComponent builds, on
		// the element's vertices rather than its dofs -- the analytic field is
		// exact everywhere, so this is the same question without a solve.
		std::vector<char> carriesPlasma( static_cast<std::size_t>( out.elements ), 0 );
		std::vector<char> interior( static_cast<std::size_t>( out.elements ), 1 );
		double bestPsiN = -std::numeric_limits<double>::infinity();

		for ( int e = 0; e < out.elements; ++e )
		{
			mfem::Array<int> vertices;
			mesh.GetElementVertices( e, vertices );
			for ( int i = 0; i < vertices.Size(); ++i )
			{
				double const *v = mesh.GetVertex( vertices[ i ] );
				double const psiN = ( equilibrium.psi( v[ 0 ], v[ 1 ] ) - psiBnd )/span;
				if ( psiN > 0.0 )
					carriesPlasma[ static_cast<std::size_t>( e ) ] = 1;
				else
					interior[ static_cast<std::size_t>( e ) ] = 0;
				bestPsiN = std::max( bestPsiN, psiN );
			}

			// The element containing the saddle, by the ordinary inverse map.
			mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( e, &transformation );
			mfem::Vector target( 2 );
			target( 0 ) = xPointR;
			target( 1 ) = xPointZ;
			mfem::IntegrationPoint reference;
			if ( transformation.TransformBack( target, reference ) ==
			     mfem::InverseElementTransformation::Inside )
				out.xPointElement = e;
		}
		out.largestPsiN = bestPsiN;

		for ( char c : carriesPlasma )
			out.candidates += c ? 1 : 0;

		auto belowTheXPoint = [ & ]( int e )
		{
			mfem::Vector centre;
			mesh.GetElementCenter( e, centre );
			return centre( 1 ) < xPointZ;
		};

		auto below = [ & ]( std::vector<char> const &keep )
		{
			int count = 0;
			for ( int e = 0; e < out.elements; ++e )
				if ( keep[ static_cast<std::size_t>( e ) ] && belowTheXPoint( e ) )
					++count;
			return count;
		};

		auto lostAbove = [ & ]( std::vector<char> const &keep )
		{
			int count = 0;
			for ( int e = 0; e < out.elements; ++e )
				if ( carriesPlasma[ static_cast<std::size_t>( e ) ]
				     && !keep[ static_cast<std::size_t>( e ) ]
				     && !belowTheXPoint( e ) )
					++count;
			return count;
		};

		// int |F| over a set of elements, with the pointwise test still applied
		// INSIDE each: the fill decides which elements the source may live in and
		// the value test decides where in them, exactly as the assembly does.
		auto sourceMass = [ & ]( std::vector<char> const &keep )
		{
			double total = 0.0;
			for ( int e = 0; e < out.elements; ++e )
			{
				if ( !keep[ static_cast<std::size_t>( e ) ] )
					continue;
				mfem::IsoparametricTransformation transformation;
				mesh.GetElementTransformation( e, &transformation );
				mfem::IntegrationRule const &rule =
					mfem::IntRules.Get( mesh.GetElementBaseGeometry( e ), 6 );
				for ( int i = 0; i < rule.GetNPoints(); ++i )
				{
					mfem::IntegrationPoint const &ip = rule.IntPoint( i );
					transformation.SetIntPoint( &ip );
					mfem::Vector point;
					transformation.Transform( ip, point );
					double const psi = equilibrium.psi( point( 0 ), point( 1 ) );
					if ( ( psi - psiBnd )*span <= 0.0 )
						continue;
					total += ip.weight*transformation.Weight()
					         *std::abs( equilibrium.f( point( 0 ), point( 1 ), psi ) );
				}
			}
			return total;
		};

		PlasmaComponent component;
		{
			mfem::Table const &neighbours = mesh.ElementToElementTable();
			std::vector<int> offsets{ 0 };
			std::vector<int> list;
			for ( int e = 0; e < out.elements; ++e )
			{
				int const *row = neighbours.GetRow( e );
				for ( int i = 0; i < neighbours.RowSize( e ); ++i )
					if ( row[ i ] >= 0 )
						list.push_back( row[ i ] );
				offsets.push_back( static_cast<int>( list.size() ) );
			}
			component.setAdjacency( std::move( offsets ), std::move( list ) );
		}

		auto keptBy = [ & ]( PlasmaComponent const &c )
		{
			std::vector<char> keep( static_cast<std::size_t>( out.elements ), 0 );
			for ( int e = 0; e < out.elements; ++e )
				if ( c.holds( e ) )
					keep[ static_cast<std::size_t>( e ) ] = 1;
			return keep;
		};

		out.pointwiseBelow = below( carriesPlasma );
		out.pointwiseMass = sourceMass( carriesPlasma );

		component.fill( carriesPlasma, axisElement );
		out.oneRuleBelow = below( keptBy( component ) );

		component.fill( carriesPlasma, axisElement, out.xPointElement );
		std::vector<char> const blocked = keptBy( component );
		out.blockedBelow = below( blocked );
		out.blockedMass = sourceMass( blocked );

		// The interior rule alone -- the fill with no band at all -- as the
		// control for what the watershed adds back.
		{
			std::vector<char> keep( static_cast<std::size_t>( out.elements ), 0 );
			for ( int e = 0; e < out.elements; ++e )
				keep[ static_cast<std::size_t>( e ) ] =
					interior[ static_cast<std::size_t>( e ) ];
			// Only the seed's own interior component, not every interior element.
			component.fill( interior, axisElement );
			keep = keptBy( component );
			out.interiorOnlyBelow = below( keep );
			out.interiorOnlyLost = lostAbove( keep );
		}

		// WHAT THE SOLVER RUNS: traverse the interior, share the band out by a
		// watershed. No depth, no X-point, no parameter.
		component.fill( interior, carriesPlasma, axisElement );
		std::vector<char> const twoRule = keptBy( component );
		out.interiorComponents = component.componentCount();
		out.twoRuleRing = component.ringNodes();
		out.twoRuleBelow = below( twoRule );
		out.twoRuleMass = sourceMass( twoRule );
		out.twoRuleLostAbove = lostAbove( twoRule );
		out.twoRuleLostBelow = below( twoRule );
		for ( int e = 0; e < out.elements; ++e )
			out.twoRuleKept += twoRule[ static_cast<std::size_t>( e ) ] ? 1 : 0;

		return out;
	};

	std::printf( "\n  A DIVERTED PLASMA: DOES A FACE-NEIGHBOUR FILL LEAK THROUGH "
	             "THE SADDLE?\n" );
	std::printf( "    Soloviev::iterExample2, X-point at ( %.5f, %.3f ) in closed "
	             "form\n", xPointR, xPointZ );

	std::vector<Reading> readings;
	for ( int n : { 24, 48, 96 } )
		readings.push_back( measure( n ) );

	/*
	 * THE BAND IS WHAT HAS TO BE SHARED OUT, AND SHARING IT IS NOT OPTIONAL.
	 *
	 * The interior rule alone separates the two lobes perfectly and drops an
	 * O( h ) band of the plasma edge with them: 179, 360, 720 elements over a
	 * fourfold refinement each time, which is the 1/h of a one-element-thick
	 * band around the boundary rather than anything that vanishes. FB-4's
	 * result is that psi* keeps k+2 exactly when k <= j, and an element-aligned
	 * support would put an O( h^(1+j) ) perturbation under it, so the band has
	 * to come back.
	 *
	 * The watershed brings it back and keeps the separation, at all three
	 * resolutions and with no parameter to choose.
	 */
	std::printf( "\n    THE BAND: ( elements below the saddle | plasma elements "
	             "dropped )\n" );
	std::printf( "    %-24s", "rule" );
	for ( Reading const &r : readings )
		std::printf( " %18d", r.elements );
	std::printf( "   <- mesh elements\n" );
	std::printf( "    %-24s", "interior alone" );
	for ( Reading const &r : readings )
		std::printf( " %8d | %7d", r.interiorOnlyBelow, r.interiorOnlyLost );
	std::printf( "\n    %-24s", "interior + watershed" );
	for ( Reading const &r : readings )
		std::printf( " %8d | %7d", r.twoRuleLostBelow, r.twoRuleLostAbove );
	std::printf( "\n" );

	Reading const &fine = readings.back();
	std::printf( "\n    ON THE FINEST, %d elements:\n", fine.elements );
	std::printf( "    %-32s %10s %14s %14s\n", "support", "elements",
	             "below X-point", "int |F|" );
	std::printf( "    %-32s %10d %14d %14.6e\n", "pointwise Psi > 0",
	             fine.candidates, fine.pointwiseBelow, fine.pointwiseMass );
	std::printf( "    %-32s %10d %14d %14s\n", "one-rule fill, nothing blocked",
	             fine.candidates, fine.oneRuleBelow, "( identical )" );
	std::printf( "    %-32s %10s %14d %14.6e\n", "one-rule, X-point element blocked",
	             "", fine.blockedBelow, fine.blockedMass );
	std::printf( "    %-32s %10d %14d %14.6e\n", "TWO-RULE, nothing located",
	             fine.twoRuleKept, fine.twoRuleBelow, fine.twoRuleMass );
	std::printf( "    the two-rule mask is %d interior elements and %d straddling "
	             "ones, and drops %d plasma element( s )\n",
	             fine.twoRuleKept - fine.twoRuleRing, fine.twoRuleRing,
	             fine.twoRuleLostAbove );
	std::printf( "    the largest Psi anywhere on the box is %.3f, and it is NOT "
	             "at the axis\n", fine.largestPsiN );
	std::fflush( stdout );

	for ( Reading const &r : readings )
	{
		BOOST_TEST_REQUIRE( r.xPointElement >= 0,
		                    "the X-point is not in the mesh at " << r.elements
		                    << " elements, so this case measures nothing" );

		// THE POINTWISE TEST IS WRONG, which is what makes the rest worth doing.
		BOOST_TEST( r.pointwiseBelow > 0,
		            "the pointwise test switches the source on in NO element "
		            "below the X-point, so the box needs to reach further down" );

		// AND SECTION 10.3's PREDICTION IS FALSE. It reads: "A fill over FACE
		// neighbours cannot cross a shared vertex, so it should be blocked at
		// the saddle automatically -- except in whichever element actually
		// contains the X-point ... So the expected leak is ONE ELEMENT WIDE."
		// The one-rule fill takes the WHOLE private flux region, because the
		// band of elements STRADDLING the separatrix is several elements wide
		// near a saddle and every one of them carries Psi > 0 at some vertex.
		BOOST_TEST( r.oneRuleBelow == r.pointwiseBelow,
		            "the one-rule fill reached " << r.oneRuleBelow << " of the "
		            "pointwise test's " << r.pointwiseBelow << " private-flux "
		            "elements. If it now blocks some of them the straddling band "
		            "has stopped bridging the lobes, and section 10.3's "
		            "prediction is worth re-reading rather than this assertion "
		            "being relaxed" );

		// THE TWO-RULE FILL IS THE ANSWER AND IT LOCATES NOTHING. Traversing the
		// strictly interior elements separates the lobes with no X-point
		// knowledge whatever -- which is what the prediction wanted, by a
		// different mechanism from the one it proposed.
		BOOST_TEST( r.interiorComponents >= 2,
		            "the interior rule leaves " << r.interiorComponents
		            << " component( s ), so the private flux region was never "
		            "separate and nothing is being separated" );
		BOOST_TEST( r.twoRuleBelow == 0,
		            "the two-rule fill still reaches " << r.twoRuleBelow
		            << " elements below the saddle" );

		// AND IT KEEPS THE PLASMA EDGE, which is the reason for the ring and is
		// not optional: FB-4's result is that psi* holds k+2 exactly when
		// k <= j, and a support switched off in an O( h ) band inside the plasma
		// would put an O( h^(1+j) ) perturbation under it.
		BOOST_TEST( r.twoRuleLostAbove <= 2,
		            "the two-rule fill dropped " << r.twoRuleLostAbove
		            << " elements that carry plasma and are NOT in the private "
		            "flux region. The rings are supposed to keep every "
		            "straddling element, so anything beyond the one or two "
		            "pinched at the saddle itself is the plasma edge being "
		            "eaten" );

		// THE DIFFERENCE MATTERS, in the source's own units rather than a count.
		BOOST_TEST( r.twoRuleMass < 0.95*r.pointwiseMass,
		            "the connectivity test removed "
		            << ( 1.0 - r.twoRuleMass/r.pointwiseMass )*100.0
		            << "% of int |F|, which is too little to be the private flux "
		            "region" );

		// AND IT IS NEVER WORSE THAN BLOCKING THE SADDLE'S OWN ELEMENT, which
		// is section 10.3's own proposed cure and is the control here.
		//
		// **THAT CURE IS RESOLUTION-DEPENDENT AND THE TWO-RULE FILL IS NOT**,
		// which is the second finding of this case and was not expected.
		// Blocking one element works on the two finer meshes and fails on the
		// coarsest, reaching 161 of 2304 elements below the saddle: the
		// straddling band there is wide enough relative to the plasma that
		// removing its middle element leaves the two lobes joined by its
		// neighbours. So the cure that needs an X-point finder -- stage XP-0,
		// not built -- is also the cure that stops working first.
		BOOST_TEST( r.twoRuleBelow <= r.blockedBelow,
		            "blocking the X-point's element left " << r.blockedBelow
		            << " elements below the saddle and the two-rule fill "
		            << r.twoRuleBelow << ". The two-rule fill is supposed to be "
		            "at least as good everywhere, since it needs no X-point at "
		            "all" );

		// THE BAND REALLY IS THE PLASMA EDGE, which is what says the watershed
		// is adding back something that matters rather than tidying. It is one
		// element thick around the boundary, so it must grow like 1/h.
		BOOST_TEST( r.interiorOnlyLost > 50,
		            "the interior rule alone drops only " << r.interiorOnlyLost
		            << " plasma elements, so there is no band to share out and "
		            "the watershed is buying nothing" );
		BOOST_TEST( r.interiorOnlyBelow == 0,
		            "the interior rule alone reaches " << r.interiorOnlyBelow
		            << " elements below the saddle, so the two lobes are not "
		            "separated by the traversal at all and everything below "
		            "rests on the band" );
	}

	// TWO ROUTES TO ONE MASK, where both routes work. On the finest mesh the
	// two-rule fill and the blocked one-rule fill remove the SAME region to
	// better than the ring's own width, which is what says neither answer is an
	// artefact of how it was computed.
	BOOST_TEST( std::abs( fine.blockedMass - fine.twoRuleMass )
	            < 0.02*fine.twoRuleMass,
	            "blocking the X-point's element keeps int |F| = "
	            << fine.blockedMass << " and the two-rule fill "
	            << fine.twoRuleMass << ". They are two ways of removing the same "
	            "region and should agree to about the ring's own width" );
	BOOST_TEST( fine.blockedBelow == 0 );

	// AND THE BAND SCALES LIKE THE BOUNDARY IT IS, which is what identifies it:
	// 179 / 360 / 720 over two fourfold refinements is 1/h and not h^0 or h^2.
	for ( std::size_t d = 0; d + 1 < readings.size(); ++d )
		BOOST_TEST( readings[ d + 1 ].interiorOnlyLost
		            > 1.5*readings[ d ].interiorOnlyLost,
		            "the plasma dropped by the interior rule went "
		            << readings[ d ].interiorOnlyLost << " -> "
		            << readings[ d + 1 ].interiorOnlyLost
		            << " over a fourfold refinement. It is a band one element "
		            "thick around the plasma edge, so it must grow like 1/h -- if "
		            "it does not, the watershed is covering something else" );
}


/*
 * AND IT IS LIVE ON A LIMITER CASE, WHICH IS WHAT MAKES XP-1 URGENT RATHER THAN
 * SCHEDULED.
 *
 * Section 10.3 says the defect "cannot fire today: no shipped example sets
 * ConfineToPlasma and MEQ has no diverted case". The first half is still true.
 * The second is beside the point: the half-disc below has NO X-POINT ANYWHERE,
 * it converges in six Newton steps, it is exactly the configuration
 * FreeBoundaryCoupling::theTwoBordersConvergeTogether drives at limiter
 * R = 1.20 -- psi_ax and psi_bnd reproduce that case to every printed digit --
 * and its converged `{ psi > psi_bnd }` is 705 elements in more than one piece:
 * two components of the interior on the mesh, and two lobes on the midplane
 * counted by a ray that uses no adjacency at all.
 *
 * The mechanism is more ordinary than a saddle and that is the point: psi_bnd
 * is a LEVEL, and a level of any field that is not monotone in radius cuts the
 * domain into as many pieces as it likes.
 *
 * **AND THE CASE FOUND SOMETHING ELSE ABOUT THAT FIXTURE, REPORTED HERE
 * BECAUSE NOTHING ELSE MEASURES IT.** psi_ax on this converged answer -- the
 * quantity the bordered Newton constrains, and the largest nodal psi_h by
 * definition -- is attained in an element TOUCHING r = 0, in the corner where
 * Gamma meets the symmetry axis. It reads 1.0916e-01 there against 4.4472e-02
 * as the largest value anywhere off the axis, a factor of 2.5. That corner is
 * where two different data meet (psi = 0 on the fitted axis, the transferred
 * exterior trace on Gamma_h) at the one place the lifting weight C = r
 * vanishes, and FB-A measured the flux mass ( r q, v ) giving those elements a
 * weight of order h. So the one thing the fill must not do is take psi_ax's
 * argmax for the magnetic axis, and refreshPlasmaComponent() does not: it seeds
 * off the symmetry axis, and this case is where that rule was measured into
 * existence.
 */
BOOST_AUTO_TEST_CASE( theLimiterCaseAlreadyHasMoreThanOneLobe )
{
	int const order = 2;
	int const n = 24;
	double const mu0 = 1.0;

	auto pPrime = std::make_shared<PowerProfile const>( 0.6, 1 );
	auto ggPrime = std::make_shared<PowerProfile const>( 0.05, 1 );

	HalfDisc d = makeHalfDisc( n );
	meq::ExteriorDtN const dtn( 0.0, halfDiscGamma, 4 );
	mfem::ConstantCoefficient zero( 0.0 );
	mfem::FunctionCoefficient guess( [ ]( mfem::Vector const &x )
	{
		double const dr = x( 0 ) - 0.75;
		double const dz = x( 1 );
		double const t = 1.0 - ( dr*dr + dz*dz )/( 0.40*0.40 );
		return t > 0.0 ? 0.1*t : 0.0;
	} );

	meq::NormalisedMHDSource source( pPrime, ggPrime, 0.1, mu0 );
	GradShafranovSolver solver( *d.sub, order );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-9, 1.0e-12, 150 );
	solver.setSource( source, 0.1 );
	solver.setBoundaryData( zero );
	solver.setExtension( *d.path, d.gammaHMarker );
	solver.setBoundaryFluxPoint( 1.20, 0.0 );
	solver.setExteriorCoupling( dtn );
	solver.solve();

	double const psiAxis = solver.psiAxis();
	double const psiBnd = solver.psiBoundary();
	double const span = psiAxis - psiBnd;

	std::printf( "\n  A LIMITER CASE WITH NO X-POINT, AND ITS SUPPORT IS NOT "
	             "CONNECTED\n" );
	std::printf( "    psi_ax %.9e   psi_bnd %.9e   newton %d   %d elements\n",
	             psiAxis, psiBnd, solver.newtonIterations(), d.sub->GetNE() );

	// THE PROFILE ITSELF, so that "more than one component" is visible rather
	// than only asserted. `+` is inside the pointwise support.
	std::printf( "    psi on z = 0, r = 0.05 .. 1.45:\n      " );
	int midplaneRuns = 0;
	bool previous = false;
	for ( int i = 0; i <= 56; ++i )
	{
		double const r = 0.05 + 1.40*i/56.0;
		mfem::DenseMatrix point( 2, 1 );
		point( 0, 0 ) = r;
		point( 1, 0 ) = 0.0;
		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> reference;
		d.sub->FindPoints( point, elements, reference );
		bool now = false;
		if ( elements[ 0 ] >= 0 )
		{
			double const psi = solver.potential().GetValue( elements[ 0 ],
			                                                reference[ 0 ] );
			now = ( psi - psiBnd )*span > 0.0;
		}
		if ( now && !previous )
			++midplaneRuns;
		previous = now;
		std::printf( "%c", elements[ 0 ] < 0 ? ' ' : ( now ? '+' : '.' ) );
	}
	std::printf( "\n" );

	// WHERE psi_ax IS, AND WHERE THE LARGEST psi_h OFF THE AXIS IS. If those two
	// agree the peak is a plasma; if they do not it is the corner degeneracy
	// FB-A measured, and the difference is what the seed rule exists for.
	mfem::GridFunction const &potential = solver.potential();
	mfem::Array<int> dofs;
	double onAxisPeak = -std::numeric_limits<double>::infinity();
	double offAxisPeak = -std::numeric_limits<double>::infinity();
	double onAxisR = 0.0, onAxisZ = 0.0, offAxisR = 0.0, offAxisZ = 0.0;

	for ( int e = 0; e < d.sub->GetNE(); ++e )
	{
		mfem::Array<int> vertices;
		d.sub->GetElementVertices( e, vertices );
		bool touchesAxis = false;
		for ( int i = 0; i < vertices.Size(); ++i )
			if ( d.sub->GetVertex( vertices[ i ] )[ 0 ] <= 0.0 )
				touchesAxis = true;

		potential.FESpace()->GetElementDofs( e, dofs );
		for ( int i = 0; i < dofs.Size(); ++i )
		{
			double const value = potential( dofs[ i ] );
			mfem::Vector centre;
			if ( touchesAxis && value > onAxisPeak )
			{
				onAxisPeak = value;
				d.sub->GetElementCenter( e, centre );
				onAxisR = centre( 0 );
				onAxisZ = centre( 1 );
			}
			if ( !touchesAxis && value > offAxisPeak )
			{
				offAxisPeak = value;
				d.sub->GetElementCenter( e, centre );
				offAxisR = centre( 0 );
				offAxisZ = centre( 1 );
			}
		}
	}
	std::printf( "    largest psi_h ON  the symmetry axis: %.6e at ( %.3f, %.3f )"
	             "   <- this is psi_ax\n", onAxisPeak, onAxisR, onAxisZ );
	std::printf( "    largest psi_h OFF the symmetry axis: %.6e at ( %.3f, %.3f )"
	             "   <- where the fill is seeded\n",
	             offAxisPeak, offAxisR, offAxisZ );

	// The fill AT THE CONVERGED ANSWER. The solve above ran with the support
	// OFF, so nothing here perturbed it; this is a report on the answer it
	// found, which is the honest way to ask what a confined run would have met.
	source.setPlasmaSupport( true );
	solver.refreshPlasmaComponent( solver.potential() );

	std::printf( "    %-38s %8d\n", "elements carrying Psi > 0 somewhere",
	             solver.plasmaCandidateElements() );
	std::printf( "    %-38s %8d\n", "components of the interior",
	             solver.plasmaComponentCount() );
	std::printf( "    %-38s %8d\n", "elements the fill reached",
	             solver.plasmaComponentElements() );
	std::printf( "    %-38s %8d\n", "midplane lobes, by sampling", midplaneRuns );
	std::fflush( stdout );

	// THE FIXTURE HAS SOMETHING TO FIX, ON THE MIDPLANE AND ON THE MESH. If
	// either ever reads 1 the case is vacuous and must say so rather than pass
	// quietly -- which is the whole reason plasmaComponentCount() is a public
	// number and not an implementation detail of the fill.
	BOOST_TEST( midplaneRuns >= 2,
	            "a ray across the midplane crosses " << midplaneRuns << " lobe( s )"
	            " of { psi > psi_bnd }. Section 10.3's claim that the pointwise "
	            "support test is a LATENT defect would then be correct after all, "
	            "and this case measures nothing" );
	BOOST_TEST( solver.plasmaComponentCount() > 1,
	            "the converged support has one component on the mesh against "
	            << midplaneRuns << " lobes on the midplane. The two routes "
	            "disagree about whether it is connected" );

	// AND THE FILL REMOVES IT.
	BOOST_TEST( solver.plasmaComponentElements() <
	            solver.plasmaCandidateElements(),
	            "the fill reached every candidate element, so it separated "
	            "nothing" );

	// THE SEED RULE IS DOING REAL WORK HERE, and this is the assertion that
	// records why it exists: psi_ax is 2.5x the largest value off the symmetry
	// axis, so seeding at psi_ax's own argmax would start the fill in a corner
	// element where the flux mass degenerates.
	BOOST_TEST( onAxisPeak > 1.5*offAxisPeak,
	            "psi_ax is " << onAxisPeak << " on the symmetry axis against "
	            << offAxisPeak << " off it. If those have come together the "
	            "corner degeneracy is gone, the off-axis seed rule is buying "
	            "nothing here, and the comment above it should say so" );
	/*
	 * AND psi_ax IS NO LONGER THAT PEAK, WHICH IS THE POINT OF OPTION 3.
	 *
	 * This asserted `onAxisPeak == psiAxis` -- "the same quantity by definition"
	 * -- and under AxisConstraint::NodalMaximum it was. It is not any more:
	 * psi_ax is the flux at the LOCATED magnetic axis, a zero of q_h, so on this
	 * fixture it reads 1.26e-01 at the plasma while the r = 0 layer sits at
	 * 3.08e-01. The border no longer follows the layer, which is precisely the
	 * defect FREE-BOUNDARY-PLAN.md section 11 opens with.
	 *
	 * So what is asserted is the SEPARATION, and it is the stronger statement:
	 * the layer is still there -- this fixture carries the 1/r pole of section
	 * 11.3 and nothing here repairs the FIELD -- and psi_ax has stopped
	 * reporting it. Both halves matter, and a test that only checked the second
	 * would pass on a fixture whose layer had quietly gone away.
	 */
	BOOST_TEST( onAxisPeak > 1.5*psiAxis,
	            "the symmetry-axis layer reads " << onAxisPeak << " and psi_ax "
	            "reads " << psiAxis << ". This case exists because those are far "
	            "apart: if they have come together either the layer is gone -- "
	            "check whether the fixture's gg' still fails to vanish at r = 0, "
	            "which is what puts it there -- or psi_ax is following the layer "
	            "again, which would mean the axis constraint has reverted to "
	            "AxisConstraint::NodalMaximum" );
}


/*
 * WHAT THE FILL COSTS NEWTON, AND IT IS THE QUESTION THAT DECIDES WHERE THE
 * FILL MAY LIVE.
 *
 * FB-4 measured that with a fixed quadrature rule and p'( 0 ) != 0 the
 * assembled residual is genuinely DISCONTINUOUS in the unknowns, and that
 * Newton then converges from nowhere -- not from the exact solution, not under
 * PicardThenNewton. The discriminator was not the size of the step but its
 * behaviour under refinement: a JUMP's largest step between neighbouring
 * samples does not shrink as the sampling is quartered, and a steep slope's
 * does.
 *
 * THE STRUCTURAL WORRY IS THAT A CONNECTIVITY TEST HAS NO CONTINUITY TO APPEAL
 * TO. When a lobe joins or leaves the seed's component its whole integral joins
 * or leaves the residual, O( 1 ), however smoothly the profile vanishes at the
 * plasma edge -- so the comfortable argument that j >= 1 makes everything C^1
 * would not apply. That is a prediction, and this case measures it rather than
 * repeating it: slide psi_bnd, which moves the level across the domain, and
 * watch the largest step between neighbouring samples as the sampling is
 * refined.
 *
 * **MEASURED, THE SHIPPED RULE DOES NOT PRODUCE A JUMP AT ALL, AND WHY IS THE
 * FINDING.** Both columns quarter as the sampling is quartered. The reason is
 * the watershed: what changes hands as psi_bnd slides is a STRADDLING element,
 * where Psi is near zero -- and at j >= 1 the profile vanishes there, so an
 * element transferring between components carries a vanishing integral with it.
 * The O( 1 ) jump the structure predicts needs a lobe with real mass to switch,
 * which needs the level to cross a local extremum of psi rather than sweep
 * past it.
 *
 * IT IS A PROPERTY OF THE RULE AND NOT OF THE PROBLEM, which is what makes it
 * worth asserting. An earlier version of the fill took a fixed number of RINGS
 * around the interior instead of sharing the band out, and on this same
 * experiment its column read 1.68 and then 1.20 against the control's 3.99 and
 * 4.00 -- settling on a floor at 1.6e-04 of int |F|, which is a jump. Rings
 * transfer whole lobes; a watershed transfers one element at a time.
 *
 * The integrand is `int |F|` over the confined support, which is what the
 * residual's source term is built from and is a scalar -- so the experiment is
 * PlasmaEdgeConvergence's, on the quantity XP-1 changes.
 */
BOOST_AUTO_TEST_CASE( theConnectivityChangeIsNotAJumpAtTheSizeThatMatters )
{
	int const order = 2;
	int const n = 24;

	auto pPrime = std::make_shared<PowerProfile const>( 0.6, 1 );
	auto ggPrime = std::make_shared<PowerProfile const>( 0.05, 1 );

	HalfDisc d = makeHalfDisc( n );
	meq::ExteriorDtN const dtn( 0.0, halfDiscGamma, 4 );
	mfem::ConstantCoefficient zero( 0.0 );
	mfem::FunctionCoefficient guess( [ ]( mfem::Vector const &x )
	{
		double const dr = x( 0 ) - 0.75;
		double const dz = x( 1 );
		double const t = 1.0 - ( dr*dr + dz*dz )/( 0.40*0.40 );
		return t > 0.0 ? 0.1*t : 0.0;
	} );

	meq::NormalisedMHDSource source( pPrime, ggPrime, 0.1, 1.0 );
	GradShafranovSolver solver( *d.sub, order );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-9, 1.0e-12, 150 );
	solver.setSource( source, 0.1 );
	solver.setBoundaryData( zero );
	solver.setExtension( *d.path, d.gammaHMarker );
	solver.setExteriorCoupling( dtn );
	solver.solve();

	double const psiAxis = solver.psiAxis();
	source.setPlasmaSupport( true );

	// The potential block alone, which refreshPlasmaComponent() accepts beside
	// the full NPC unknown -- it is the only block the fill reads.
	mfem::Vector const state( solver.potential() );

	/*
	 * `int |F|` over the confined support, as a function of psi_bnd.
	 *
	 * The element loop and the quadrature are meq::SourceIntegrator's own, so
	 * this is the residual's source term summed rather than a proxy for it, and
	 * the two connectivity settings differ in nothing but the mask.
	 */
	auto sourceMass = [ & ]( double psiBnd, bool connected )
	{
		source.setNormalisation( psiAxis, psiBnd );
		solver.setPlasmaConnectivity( connected
			? GradShafranovSolver::PlasmaConnectivity::Component
			: GradShafranovSolver::PlasmaConnectivity::Pointwise );
		solver.refreshPlasmaComponent( state );

		mfem::GridFunction const &potential = solver.potential();
		mfem::FiniteElementSpace const &space = *potential.FESpace();
		mfem::Mesh &mesh = *space.GetMesh();

		double total = 0.0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			if ( !solver.elementInPlasma( e ) )
				continue;

			mfem::FiniteElement const &el = *space.GetFE( e );
			mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( e, &transformation );
			mfem::IntegrationRule const &rule =
				mfem::IntRules.Get( el.GetGeomType(), 2*el.GetOrder() + 4 );

			for ( int i = 0; i < rule.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( i );
				transformation.SetIntPoint( &ip );
				mfem::Vector point;
				transformation.Transform( ip, point );
				double const psi = potential.GetValue( e, ip );
				total += ip.weight*transformation.Weight()
				         *std::abs( source.f( point( 0 ), point( 1 ), psi ) );
			}
		}
		return total;
	};

	// A window of psi_bnd wide enough that the level sweeps across the lobes.
	double const lower = 0.0;
	double const upper = 0.5*psiAxis;

	std::printf( "\n  IS A CONNECTIVITY CHANGE A JUMP? SLIDE psi_bnd AND REFINE "
	             "THE SAMPLING\n" );
	std::printf( "    the largest step between neighbouring samples, and the "
	             "factor it falls by\n" );
	std::printf( "    %-24s %14s %14s %8s %14s %8s\n", "support",
	             "401 samples", "1601", "x", "6401", "x" );

	double worstConnected[ 3 ] = { 0.0, 0.0, 0.0 };
	double worstPointwise[ 3 ] = { 0.0, 0.0, 0.0 };
	double scale = 0.0;

	int samples = 401;
	for ( int level = 0; level < 3; ++level )
	{
		double previousConnected = sourceMass( lower, true );
		double previousPointwise = sourceMass( lower, false );
		scale = std::max( scale, previousPointwise );

		for ( int i = 1; i <= samples; ++i )
		{
			double const psiBnd = lower + ( upper - lower )*i/samples;
			double const connected = sourceMass( psiBnd, true );
			double const pointwise = sourceMass( psiBnd, false );

			worstConnected[ level ] = std::max(
				worstConnected[ level ], std::abs( connected - previousConnected ) );
			worstPointwise[ level ] = std::max(
				worstPointwise[ level ], std::abs( pointwise - previousPointwise ) );

			previousConnected = connected;
			previousPointwise = pointwise;
			scale = std::max( scale, pointwise );
		}
		samples *= 4;
	}

	std::printf( "    %-24s %14.4e %14.4e %8.2f %14.4e %8.2f\n",
	             "pointwise Psi > 0", worstPointwise[ 0 ], worstPointwise[ 1 ],
	             worstPointwise[ 0 ]/worstPointwise[ 1 ], worstPointwise[ 2 ],
	             worstPointwise[ 1 ]/worstPointwise[ 2 ] );
	std::printf( "    %-24s %14.4e %14.4e %8.2f %14.4e %8.2f\n",
	             "connected component", worstConnected[ 0 ], worstConnected[ 1 ],
	             worstConnected[ 0 ]/worstConnected[ 1 ], worstConnected[ 2 ],
	             worstConnected[ 1 ]/worstConnected[ 2 ] );
	std::printf( "    int |F| itself is at most %.4e, so the connected column's "
	             "worst step is %.4f%% of it\n",
	             scale, 100.0*worstConnected[ 2 ]/scale );
	std::fflush( stdout );

	/*
	 * THE ANSWER: THE SHIPPED RULE IS AS CONTINUOUS AS THE POINTWISE ONE.
	 *
	 * The control quarters as the sampling is quartered, which is what a
	 * Lipschitz function does, and so does the connected column -- 3.94 and then
	 * 3.96 against 3.99 and 4.00. There is no floor to settle on, so on this
	 * configuration the connectivity test contributes no discontinuity that this
	 * experiment can see, and the fill may live inside the Newton loop rather
	 * than being frozen per solve and re-decided between them.
	 *
	 * That is not a proof that no configuration produces a jump. What is
	 * asserted is the property the shipped rule has to have: that the elements
	 * changing hands are the ones where the source is smallest, so that a
	 * connectivity change is a perturbation of the size of an element rather
	 * than the size of a lobe.
	 */
	bool const controlQuarters = worstPointwise[ 0 ] > 3.0*worstPointwise[ 1 ]
	                             && worstPointwise[ 1 ] > 3.0*worstPointwise[ 2 ];
	BOOST_TEST( controlQuarters,
	            "the pointwise column fell by "
	            << worstPointwise[ 0 ]/worstPointwise[ 1 ] << " and then "
	            << worstPointwise[ 1 ]/worstPointwise[ 2 ]
	            << " over quarterings of the sampling. It is supposed to fall by "
	            "four -- p' vanishes at the edge, so F -> 0 continuously as a "
	            "point crosses -- and if it does not then this experiment cannot "
	            "separate a jump from a slope and neither column means anything" );

	bool const connectedQuarters = worstConnected[ 0 ] > 3.0*worstConnected[ 1 ]
	                               && worstConnected[ 1 ] > 3.0*worstConnected[ 2 ];
	BOOST_TEST( connectedQuarters,
	            "the connected column fell by "
	            << worstConnected[ 0 ]/worstConnected[ 1 ] << " and then "
	            << worstConnected[ 1 ]/worstConnected[ 2 ]
	            << " where the control fell by four. A column that stops "
	            "shrinking is a JUMP, which is what a rule that transfers whole "
	            "lobes gives -- the ring-depth fill this replaced read 1.68 and "
	            "1.20 here. The watershed transfers one straddling element at a "
	            "time, where Psi is near zero and the profile with it, so if this "
	            "has started jumping the band is no longer what changes hands" );

	BOOST_TEST( worstConnected[ 2 ] < 1.0e-3*scale,
	            "the largest step from a connectivity change is "
	            << worstConnected[ 2 ] << " against an int |F| of " << scale
	            << ", which is " << 100.0*worstConnected[ 2 ]/scale << "%. A lobe "
	            "leaving the component takes its whole integral with it, so a "
	            "step that large means the fill is cutting the plasma rather "
	            "than trimming pockets off it -- and FB-4 measured what a "
	            "residual with an O( 1 ) discontinuity does to Newton" );

	/*
	 * WHAT THIS DOES NOT ESTABLISH, SAID HERE RATHER THAN LEFT TO BE ASSUMED.
	 *
	 * The natural next question is whether a CONFINED solve on this same
	 * half-disc converges with the fill live, and it cannot be answered on this
	 * fixture: run with ConfineToPlasma on, it reaches the iteration cap at
	 * 4.24e-02 with the connectivity test OFF and 6.11e-02 with it on. Neither
	 * converges, so the pair is not a measurement of what the fill costs -- it
	 * is a measurement of a configuration that does not converge confined, which
	 * CLAUDE.md already records for this geometry (the tabulated profile over
	 * [ 0, 1 ] "creeps at 0.99 a step", and ConfineToPlasma finds a different
	 * branch at four times the axis flux).
	 *
	 * So the end-to-end statement available today is the one
	 * theFillDoesNotMoveASingleLobeAnswer makes: on a confined solve whose
	 * converged support IS connected, the fill costs nothing -- the same answer
	 * to 9.5e-10 against a stopping tolerance of 1e-9, in 24 Newton steps
	 * against 47. A confined solve that converges AND whose converged support is
	 * disconnected is what would close the gap, and MEQ has no such case; the
	 * two properties have never yet been met together.
	 */
}


/*
 * THE CONTROL, AND WITHOUT IT EVERY CASE ABOVE IS COMPATIBLE WITH A FILL THAT
 * QUIETLY PERTURBS THE ANSWER.
 *
 * Two halves, and they are different statements.
 *
 *   AN UNCONFINED SOLVE MUST NOT SEE THE FILL AT ALL. With
 *   setPlasmaSupport() off there is no moving support, no mask is built, no
 *   adjacency is even computed, and every existing solve in this tree is
 *   untouched by XP-1 existing.
 *
 *   AND A CONFINED SOLVE WHOSE CONVERGED SUPPORT IS CONNECTED MUST REACH THE
 *   SAME ANSWER IN NO MORE STEPS. That is the sharper half, and it is a
 *   TOLERANCE rather than a bit: the fill refreshes before every residual, so
 *   at an intermediate iterate -- where psi is still a long way from the
 *   answer -- the level really can fragment and the fill really can drop a
 *   piece. The two runs then take different routes to the same root. What must
 *   not happen is that the route costs anything.
 *
 * The unconfined half IS bit for bit, and has to be: setPlasmaSupport() is off
 * there, so no mask is ever built and the assembly takes the identical branch.
 * A tolerance on that one would let a real perturbation through.
 */
BOOST_AUTO_TEST_CASE( theFillDoesNotMoveASingleLobeAnswer )
{
	int const order = 2;

	/*
	 * gg' = 0 HERE, WHERE EVERY OTHER CASE IN THIS FILE USES 0.05, AND IT IS A
	 * PHYSICS CHANGE RATHER THAN A TUNING ONE.
	 *
	 * This domain reaches r = 0, and F( 0, z, . ) is gg' and nothing else --
	 * p' is killed by its own r^2. So a non-zero gg' is a toroidal current
	 * density diverging like 1/r on the symmetry axis, which
	 * GradShafranovSolver::checkAxisSource() is the guard for and
	 * FREE-BOUNDARY-PLAN.md section 11.3 is the account of. With gg' = 0 the
	 * axis carries no current, which is what a plasma reaching r = 0 must
	 * satisfy -- a levitated dipole and a magnetic mirror both do, and both
	 * have no toroidal field for exactly this reason.
	 *
	 * THIS CASE IS A CONTROL ABOUT THE FILL AND NOT ABOUT FREE BOUNDARY, so it
	 * needs a solve whose support is connected and nothing more; it does not
	 * need to be a machine. The cases above keep gg' = 0.05 deliberately,
	 * because what they measure is connectivity on the fixture as it was
	 * published.
	 */
	auto pPrime = std::make_shared<PowerProfile const>( 0.6, 1 );
	auto ggPrime = std::make_shared<PowerProfile const>( 0.0, 1 );
	mfem::ConstantCoefficient zero( 0.0 );

	struct Result
	{
		double axis = 0.0;
		int iterations = 0;
		int candidates = 0;
		int reached = 0;
		int components = 0;
		std::vector<double> potential;
	};

	auto worstDifference = [ ]( Result const &a, Result const &b )
	{
		double worst = 0.0;
		BOOST_TEST_REQUIRE( a.potential.size() == b.potential.size() );
		for ( std::size_t i = 0; i < a.potential.size(); ++i )
			worst = std::max( worst,
			                  std::abs( a.potential[ i ] - b.potential[ i ] ) );
		return worst;
	};

	std::printf( "\n  THE CONTROL: THE FILL MUST NOT MOVE AN ANSWER IT DOES NOT "
	             "SEPARATE\n" );

	/*
	 * HALF ONE: the half-disc of the cases above, UNCONFINED.
	 */
	{
		HalfDisc d = makeHalfDisc( 24 );
		meq::ExteriorDtN const dtn( 0.0, halfDiscGamma, 4 );
		mfem::FunctionCoefficient guess( [ ]( mfem::Vector const &x )
		{
			double const dr = x( 0 ) - 0.75;
			double const dz = x( 1 );
			double const t = 1.0 - ( dr*dr + dz*dz )/( 0.40*0.40 );
			return t > 0.0 ? 0.1*t : 0.0;
		} );

		auto run = [ & ]( GradShafranovSolver::PlasmaConnectivity choice )
		{
			meq::NormalisedMHDSource source( pPrime, ggPrime, 0.1, 1.0 );
			GradShafranovSolver solver( *d.sub, order );
			solver.setPlasmaConnectivity( choice );
			solver.setInitialGuess( guess );
			solver.setNewtonControl( 1.0e-9, 1.0e-12, 150 );
			solver.setSource( source, 0.1 );
			solver.setBoundaryData( zero );
			solver.setExtension( *d.path, d.gammaHMarker );
			solver.setExteriorCoupling( dtn );
			solver.solve();

			Result out;
			out.axis = solver.psiAxis();
			out.iterations = solver.newtonIterations();
			mfem::GridFunction const &potential = solver.potential();
			out.potential.assign( potential.GetData(),
			                      potential.GetData() + potential.Size() );
			return out;
		};

		Result const pointwise =
			run( GradShafranovSolver::PlasmaConnectivity::Pointwise );
		Result const component =
			run( GradShafranovSolver::PlasmaConnectivity::Component );
		double const worst = worstDifference( pointwise, component );

		std::printf( "    UNCONFINED half-disc, %zu dofs\n",
		             pointwise.potential.size() );
		std::printf( "      pointwise  psi_ax %.9e  newton %d\n",
		             pointwise.axis, pointwise.iterations );
		std::printf( "      component  psi_ax %.9e  newton %d\n",
		             component.axis, component.iterations );
		std::printf( "      worst difference in psi_h: %.3e\n", worst );

		BOOST_TEST( worst == 0.0,
		            "switching PlasmaConnectivity moved an UNCONFINED solve by "
		            << worst << ". The fill is inert unless setPlasmaSupport() "
		            "is on, so this is the connectivity test reaching a solve "
		            "that never asked for it" );
		BOOST_TEST( pointwise.iterations == component.iterations );
		BOOST_TEST( pointwise.axis == component.axis );
	}

	/*
	 * HALF TWO: a CONFINED solve on a fitted rectangle, where the support is one
	 * lobe by construction.
	 *
	 * psi = 0 on the whole boundary and F > 0 inside, so psi > 0 = psi_bnd
	 * throughout the interior and `{ Psi > 0 }` is the domain itself. The fill
	 * therefore reaches every candidate element and holds() is true everywhere
	 * it was true before -- which is exactly the configuration in which a
	 * perturbation would be a defect rather than a difference.
	 */
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			16, 16, mfem::Element::TRIANGLE, false, 0.8, 1.6 );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			mesh.GetVertex( v )[ 0 ] += 0.6;
			mesh.GetVertex( v )[ 1 ] -= 0.8;
		}

		mfem::FunctionCoefficient guess( [ ]( mfem::Vector const &x )
		{
			double const dr = x( 0 ) - 1.0;
			double const dz = x( 1 );
			double const t = 1.0 - ( dr*dr + dz*dz )/( 0.40*0.40 );
			return t > 0.0 ? 0.1*t : 0.0;
		} );

		auto run = [ & ]( GradShafranovSolver::PlasmaConnectivity choice )
		{
			meq::NormalisedMHDSource source( pPrime, ggPrime, 0.1, 1.0 );
			source.setPlasmaSupport( true );

			GradShafranovSolver solver( mesh, order );
			solver.setPlasmaConnectivity( choice );
			solver.setInitialGuess( guess );
			solver.setNewtonControl( 1.0e-9, 1.0e-12, 60 );
			solver.setSource( source, 0.1 );
			solver.setBoundaryData( zero );
			solver.solve();

			Result out;
			out.axis = solver.psiAxis();
			out.iterations = solver.newtonIterations();
			out.candidates = solver.plasmaCandidateElements();
			out.reached = solver.plasmaComponentElements();
			out.components = solver.plasmaComponentCount();
			mfem::GridFunction const &potential = solver.potential();
			out.potential.assign( potential.GetData(),
			                      potential.GetData() + potential.Size() );
			return out;
		};

		Result const pointwise =
			run( GradShafranovSolver::PlasmaConnectivity::Pointwise );
		Result const component =
			run( GradShafranovSolver::PlasmaConnectivity::Component );
		double const worst = worstDifference( pointwise, component );

		std::printf( "    CONFINED fitted rectangle, %d elements, %zu dofs\n",
		             mesh.GetNE(), pointwise.potential.size() );
		std::printf( "      pointwise  psi_ax %.9e  newton %d\n",
		             pointwise.axis, pointwise.iterations );
		std::printf( "      component  psi_ax %.9e  newton %d   "
		             "%d component( s ), %d of %d elements\n",
		             component.axis, component.iterations, component.components,
		             component.reached, component.candidates );
		std::printf( "      worst difference in psi_h: %.3e\n", worst );
		std::fflush( stdout );

		// THE FILL RAN AND REACHED EVERYTHING, which is what makes the
		// bit-identity below a statement about the fill rather than about a
		// branch that was never taken.
		BOOST_TEST( component.candidates > 0,
		            "the confined solve has no candidate elements at all, so the "
		            "support test is not on and this half measures nothing" );
		BOOST_TEST( component.components == 1,
		            "the confined support has " << component.components
		            << " components on a fitted rectangle where psi = 0 on the "
		            "whole boundary and F > 0 inside. It is supposed to be one "
		            "lobe by construction, so this is not the control it claims "
		            "to be" );
		BOOST_TEST( component.reached == component.candidates,
		            "the fill reached " << component.reached << " of "
		            << component.candidates << " candidates on a single-lobe "
		            "support, so it dropped elements it had no reason to" );

		/*
		 * AND THE ANSWER DID NOT MOVE -- BUT TO THE SOLVER'S OWN TOLERANCE AND
		 * NOT BIT FOR BIT, AND THE DIFFERENCE IS INSTRUCTIVE RATHER THAN A
		 * WEAKENING.
		 *
		 * At CONVERGENCE the fill holds all 512 elements, so the assembled
		 * system is identical and the two are converged solutions of one
		 * problem. What differs is the PATH: the support is refreshed before
		 * every residual, and at some intermediate iterate -- where psi is still
		 * a long way from the answer -- the level psi = psi_bnd really did cut
		 * the domain into pieces, and the fill really did drop one. So the two
		 * runs take different routes to the same root and agree to 9.5e-10
		 * against a stopping tolerance of 1e-9, which is where two converged
		 * answers to one system are entitled to agree.
		 *
		 * AND THE ROUTE WITH THE FILL IS THE SHORTER ONE: 24 Newton steps
		 * against 47. That is not what this case was written to find and it is
		 * not asserted as a benefit -- one configuration is not a result -- but
		 * it is the opposite of the cost the fill was suspected of, and it is
		 * the same mechanism: a pointwise support switches the source on in
		 * pockets the iterate throws up, and the iterate then has to un-do them.
		 */
		BOOST_TEST( worst < 1.0e-8,
		            "the connectivity test moved a CONFINED single-lobe solve by "
		            << worst << ". At convergence the fill holds every element, "
		            "so the two are converged solutions of the SAME discrete "
		            "system and must agree to about the stopping tolerance" );
		BOOST_TEST( component.iterations <= pointwise.iterations,
		            "the confined solve took " << pointwise.iterations
		            << " Newton steps without the fill and "
		            << component.iterations << " with it. The fill is refreshed "
		            "before every residual, so a support that fragments and "
		            "re-joins mid-iteration is a discontinuity Newton has to "
		            "step over -- and if that is now costing steps rather than "
		            "saving them, the jump measured in "
		            "theConnectivityChangeIsNotAJumpAtTheSizeThatMatters has grown" );
		BOOST_TEST( std::abs( pointwise.axis - component.axis ) < 1.0e-8,
		            "psi_ax came out " << pointwise.axis << " without the fill "
		            "and " << component.axis << " with it" );
	}
}
