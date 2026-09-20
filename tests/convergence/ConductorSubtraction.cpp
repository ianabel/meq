// COIL-SUBTRACTION-PLAN.md CS-2: the split psi = psi_c + psi_p on a
// FIXED-boundary case, where nothing else moves.
//
// THE ACCEPTANCE IS AN IDENTITY RATHER THAN A RATE, AND THAT IS WHY IT IS
// SHARP. Put a filament inside Omega, give the solve no plasma at all, and hand
// it the filament's OWN field as the physical Dirichlet datum. Then
//
//     Delta* psi_p = -mu0 r J_plasma = 0        in Omega
//     psi_p|_Gamma = psi|_Gamma - psi_c|_Gamma = 0
//
// so psi_p is identically zero and the reported total psi_c + psi_p is the
// filament's field EXACTLY -- at every point, including arbitrarily close to
// the conductor, which no meshed solve of any degree can do. That is the plan's
// central claim, §2(b), reduced to something a zero-tolerance assertion can
// check.
//
// AND IT DISCRIMINATES, WHICH A "SOLVE AND LOOK REASONABLE" TEST WOULD NOT.
// Drop the subtraction and the same configuration imposes psi_c|_Gamma on the
// boundary of a source-free problem, so psi_p becomes the harmonic extension of
// that datum -- a perfectly convergent solve of the wrong problem, large, smooth
// and entirely plausible. The second case below measures exactly that, so the
// first one cannot pass by the shift being absent.

#define BOOST_TEST_MODULE ConductorSubtraction
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "mfem.hpp"

#include <limits>
#include <memory>

#include "meq/ConductorField.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "meq/FieldViews.hpp"
#include "meq/GradShafranov.hpp"

namespace
{
	// A vacuum: no plasma current anywhere. The whole right-hand side of the
	// remainder equation, and it is zero by construction rather than by a
	// profile that happens to vanish.
	class VacuumSource : public meq::Source
	{
		public:
			double f( double, double, double ) const override
			{
				return 0.0;
			}

			double dFdPsi( double, double, double ) const override
			{
				return 0.0;
			}
	};

	// A box well off the axis, so that 1/r is bounded and nothing here is about
	// the axis. 8 x 8 puts vertices on multiples of 0.1 in r and 0.1 in z.
	mfem::Mesh makeBox( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false, 0.8, 0.8 );
		// Shift into [ 0.6, 1.4 ] x [ -0.4, 0.4 ].
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			double *x = mesh.GetVertex( v );
			x[ 0 ] += 0.6;
			x[ 1 ] -= 0.4;
		}
		return mesh;
	}

	// DELIBERATELY NOT ON A VERTEX. The mesh's vertices are at multiples of 0.1
	// in both coordinates, so ( 1.03, 0.017 ) is on none of them -- which is
	// what setConductorField() requires and what the refusal case below checks
	// by moving it onto one.
	meq::ConductorField insideFilament()
	{
		meq::ConductorField field;
		field.add( meq::CurrentFilament( 1.03, 0.017, 2.5e5 ) );
		return field;
	}

	double maxAbs( mfem::GridFunction const &g )
	{
		double worst = 0.0;
		for ( int i = 0; i < g.Size(); ++i )
			worst = std::max( worst, std::fabs( g( i ) ) );
		return worst;
	}
}

// THE REMAINDER IS IDENTICALLY ZERO, AND THE TOTAL IS THE FILAMENT EXACTLY.
BOOST_AUTO_TEST_CASE( theRemainderOfAVacuumFilamentIsZero )
{
	meq::ConductorField const conductors = insideFilament();
	VacuumSource const source;

	mfem::FunctionCoefficient physicalDatum(
		[ &conductors ]( mfem::Vector const &x )
		{
			return conductors.psi( x( 0 ), x( 1 ) );
		} );

	std::printf( "\n  CS-2: a filament inside Omega, no plasma, the datum its "
	             "own field\n" );
	std::printf( "    %2s %10s %14s %14s\n", "k", "elements", "max |psi_p|",
	             "max |psi_c| on Gamma" );

	for ( int degree = 1; degree <= 3; ++degree )
	{
		mfem::Mesh mesh = makeBox( 8 );
		meq::GradShafranovSolver solver( mesh, degree );
		solver.setSource( source );
		solver.setBoundaryData( physicalDatum );
		solver.setConductorField( conductors );
		solver.solve();

		double const remainder = maxAbs( solver.potential() );

		// The scale the remainder is being compared against: what the datum
		// WOULD have imposed, and therefore how big a missing subtraction is.
		double scale = 0.0;
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			double const *x = mesh.GetVertex( v );
			scale = std::max( scale,
			                  std::fabs( conductors.psi( x[ 0 ], x[ 1 ] ) ) );
		}

		std::printf( "    %2d %10d %14.3e %14.3e\n", degree,
		             mesh.GetNE(), remainder, scale );

		// Not a tolerance on a converged answer: psi_p solves a homogeneous
		// problem with homogeneous data, so this is round-off in the assembly
		// and the trace solve and nothing else.
		BOOST_TEST( remainder < 1.0e-12*scale,
		            "psi_p is not zero: the remainder should solve Delta* "
		            "psi_p = 0 with psi_p|_Gamma = 0" );
	}
}

// AND WITHOUT THE SUBTRACTION THE SAME CONFIGURATION SOLVES THE WRONG PROBLEM,
// CONVERGENTLY. This is the control that gives the case above its teeth.
BOOST_AUTO_TEST_CASE( withoutTheSubtractionTheRemainderIsLargeAndPlausible )
{
	meq::ConductorField const conductors = insideFilament();
	VacuumSource const source;

	mfem::FunctionCoefficient physicalDatum(
		[ &conductors ]( mfem::Vector const &x )
		{
			return conductors.psi( x( 0 ), x( 1 ) );
		} );

	mfem::Mesh mesh = makeBox( 8 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( physicalDatum );
	// setConductorField() NOT called -- the only difference from the case above.
	solver.solve();

	double const unshifted = maxAbs( solver.potential() );
	std::printf( "\n  control, no setConductorField: max |psi_h| = %.6e\n",
	             unshifted );
	std::printf( "  ( it converged, and it is the harmonic extension of "
	             "psi_c|_Gamma rather than zero )\n" );

	BOOST_TEST( unshifted > 1.0e-3,
	            "the control is not large, so the case above proves nothing "
	            "about the datum shift" );
}

// THE REFUSAL, WHICH IS COIL-SUBTRACTION-PLAN.md SECTION 7.1's DECISION.
BOOST_AUTO_TEST_CASE( aFilamentOnAMeshVertexIsRefusedAndOneNearItIsNot )
{
	VacuumSource const source;
	mfem::ConstantCoefficient zero( 0.0 );

	mfem::Mesh mesh = makeBox( 8 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( zero );

	// ( 1.0, 0.0 ) IS a vertex of this mesh: the box is [ 0.6, 1.4 ] x
	// [ -0.4, 0.4 ] at 8 x 8, so vertices sit on multiples of 0.1 in both.
	meq::ConductorField onAVertex;
	onAVertex.add( meq::CurrentFilament( 1.0, 0.0, 2.5e5 ) );
	BOOST_CHECK_THROW( solver.setConductorField( onAVertex ),
	                   std::invalid_argument );

	// A MILLIMETRE AWAY IS FINE, and that is the half of the policy that says
	// this is a coincidence test and not a clearance test. psi_c is large there
	// and that is correct.
	meq::ConductorField nearAVertex;
	nearAVertex.add( meq::CurrentFilament( 1.0 + 1.0e-3, 0.0, 2.5e5 ) );
	BOOST_CHECK_NO_THROW( solver.setConductorField( nearAVertex ) );
	BOOST_TEST( solver.conductorField() == &nearAVertex );
}

// THE SOURCE AND ITS JACOBIAN SEE THE TOTAL, NOT THE REMAINDER.
//
// THE VACUUM CASE ABOVE CANNOT SEE THIS AND THAT IS WHY THIS ONE EXISTS. There
// f is identically zero, so it returns the same number whichever flux it is
// handed -- a test that cannot fail on the thing being changed. Under the split
// the solver hands SourceIntegrator the remainder psi_p while J_plasma is a
// function of the PHYSICAL flux, so f and dFdPsi must both be evaluated at
// psi_p + psi_c. Getting it wrong is silent: the solve converges, at the full
// rate, to a different equilibrium.
//
// THE COMPARISON USES ONLY PRE-EXISTING MACHINERY ON THE OTHER SIDE, which is
// what makes it a check rather than a restatement. The control does the shift
// in the CALLER -- a source wrapper that adds psi_c to the flux it is asked
// about, and a datum the caller shifts itself -- with no conductor field set at
// all, so it exercises none of the new code. If the solver's internal shift is
// missing, applied once too often, or applied to f but not to dFdPsi, the two
// arms disagree.
BOOST_AUTO_TEST_CASE( theSourceIsEvaluatedAtTheTotalFluxAndNotTheRemainder )
{
	meq::ConductorField const conductors = insideFilament();

	// Linear in psi, so dFdPsi is a nonzero constant and a Jacobian evaluated
	// at the wrong flux is a DIFFERENT number rather than the same one -- the
	// middle rung of CLAUDE.md's fixture ladder, chosen for exactly that.
	struct LinearSource : public meq::Source
	{
		double f( double r, double, double psi ) const override
		{
			return 0.8*r*psi;
		}

		double dFdPsi( double r, double, double ) const override
		{
			return 0.8*r;
		}
	};

	// The same source asked about the TOTAL, for the control arm.
	struct ShiftedSource : public meq::Source
	{
		explicit ShiftedSource( meq::ConductorField const &c ) : conductors( c )
		{
		}

		double f( double r, double z, double psi ) const override
		{
			return 0.8*r*( psi + conductors.psi( r, z ) );
		}

		double dFdPsi( double r, double, double ) const override
		{
			return 0.8*r;
		}

		meq::ConductorField const &conductors;
	};

	mfem::FunctionCoefficient physicalDatum(
		[ &conductors ]( mfem::Vector const &x )
		{
			return 0.05 + conductors.psi( x( 0 ), x( 1 ) );
		} );

	// The control's datum, shifted by the caller since no conductor field is
	// set: psi|_Gamma - psi_c|_Gamma, which is what the solver forms itself in
	// the arm under test.
	mfem::FunctionCoefficient shiftedDatum(
		[]( mfem::Vector const & )
		{
			return 0.05;
		} );

	LinearSource const plain;
	ShiftedSource const shifted( conductors );

	mfem::Mesh meshA = makeBox( 8 );
	meq::GradShafranovSolver split( meshA, 2 );
	split.setSource( plain );
	split.setBoundaryData( physicalDatum );
	split.setConductorField( conductors );
	split.solve();

	mfem::Mesh meshB = makeBox( 8 );
	meq::GradShafranovSolver control( meshB, 2 );
	control.setSource( shifted );
	control.setBoundaryData( shiftedDatum );
	control.solve();

	BOOST_TEST_REQUIRE( split.potential().Size() == control.potential().Size() );

	double worst = 0.0;
	double scale = 0.0;
	for ( int i = 0; i < split.potential().Size(); ++i )
	{
		worst = std::max( worst, std::fabs( split.potential()( i )
		                                    - control.potential()( i ) ) );
		scale = std::max( scale, std::fabs( control.potential()( i ) ) );
	}

	std::printf( "\n  the source sees the total: |split - control| = %.3e "
	             "against a scale of %.3e\n", worst, scale );
	std::printf( "  ( the control shifts in the CALLER and sets no conductor "
	             "field, so it runs none of the new code )\n" );

	BOOST_TEST( scale > 1.0e-3, "the control is trivial, so this proves nothing" );
	BOOST_TEST( worst < 1.0e-11*scale,
	            "the split's source or Jacobian is not being evaluated at "
	            "psi_c + psi_p" );
}

// THE EVALUATION ABSTRACTION: FREE WHEN THE SPLIT IS OFF, CORRECT WHEN IT IS ON.
//
// COIL-SUBTRACTION-PLAN.md §8.1 counts about sixty-eight places that read the
// solved field, and §3 names the hazard: every consumer that forgets psi_c is a
// silent wrong answer. meq::PotentialView is what makes the forgetting hard --
// it is built by NAME, physical() or remainder(), so the choice is made in a
// word at the call site rather than by knowing what the class does.
//
// BOTH HALVES ARE ASSERTED AT ZERO TOLERANCE, because both are identities.
// With no conductors the view must be the SAME NUMBER as the bare GridFunction,
// not a close one -- that is the property that lets sixty-eight call sites
// migrate without any of them moving an existing answer. With conductors it
// must be exactly the bare value plus psi_c at that point.
BOOST_AUTO_TEST_CASE( theViewIsFreeWithoutConductorsAndExactWithThem )
{
	meq::ConductorField const conductors = insideFilament();
	VacuumSource const source;

	mfem::FunctionCoefficient physicalDatum(
		[ &conductors ]( mfem::Vector const &x )
		{
			return conductors.psi( x( 0 ), x( 1 ) );
		} );

	mfem::Mesh mesh = makeBox( 8 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( physicalDatum );
	solver.solve();

	mfem::GridFunction const &solved = solver.potential();

	meq::PotentialView const bare
		= meq::PotentialView::remainder( solved );
	meq::PotentialView const nulled
		= meq::PotentialView::physical( solved, nullptr );
	meq::PotentialView const withConductors
		= meq::PotentialView::physical( solved, &conductors );

	BOOST_TEST( !bare.carriesConductors() );
	BOOST_TEST( !nulled.carriesConductors() );
	BOOST_TEST( withConductors.carriesConductors() );

	mfem::IntegrationPoint ip;
	double worstFree = 0.0;
	double worstShift = 0.0;
	double scale = 0.0;

	for ( int element = 0; element < mesh.GetNE(); ++element )
	{
		mfem::IntegrationRule const &rule = mfem::IntRules.Get(
			mesh.GetElementBaseGeometry( element ), 4 );

		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			ip = rule.IntPoint( q );
			double const raw = solved.GetValue( element, ip );

			// FREE: the null view is the bare value, bit for bit.
			worstFree = std::max( worstFree,
			                      std::fabs( bare.value( element, ip ) - raw ) );
			worstFree = std::max( worstFree,
			                      std::fabs( nulled.value( element, ip ) - raw ) );

			// EXACT: the physical view is that value plus psi_c right there.
			mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( element, &transformation );
			transformation.SetIntPoint( &ip );
			mfem::Vector point( 3 );
			transformation.Transform( ip, point );

			double const expected
				= raw + conductors.psi( point( 0 ), point( 1 ) );
			worstShift = std::max( worstShift,
			                       std::fabs( withConductors.value( element, ip )
			                                  - expected ) );
			scale = std::max( scale, std::fabs( expected ) );

			// AND THE TWO OVERLOADS AGREE, so a consumer that already holds a
			// transformation is not taking a different code path by accident.
			worstShift = std::max(
				worstShift,
				std::fabs( withConductors.value( element, ip )
				           - withConductors.value( transformation, ip ) ) );
		}
	}

	std::printf( "\n  PotentialView over %d elements: free-path error %.3e, "
	             "shifted-path error %.3e against a scale of %.3e\n",
	             mesh.GetNE(), worstFree, worstShift, scale );

	BOOST_TEST( worstFree == 0.0, boost::test_tools::tolerance( 0.0 ) );
	BOOST_TEST( worstShift == 0.0, boost::test_tools::tolerance( 0.0 ) );
	BOOST_TEST( scale > 1.0e-3, "the field is trivial, so this proves nothing" );
}

// THE TWO ROUTES AGREE, AND THE MESHED ONE CONVERGES TO THE SUBTRACTED ONE.
//
// COIL-SUBTRACTION-PLAN.md §0a-pre is a STANDING REQUIREMENT: MEQ must always be
// able to solve finite-sized coils accurately the old way, and the split is an
// option rather than a replacement. This case is what gives that teeth, and it
// is only possible because CS-1b put RECTANGLES in meq::ConductorField -- a
// rectangle is the conductor that can go either way, so it is the only one that
// can be the cross-check. A filament can only ever be subtracted.
//
// THE SUBTRACTED ARM IS THE TRUTH HERE, WHICH IS THE RIGHT WAY ROUND. With the
// coil's own field as the datum and no plasma, the remainder is identically
// zero and the total is psi_c EXACTLY -- §7.4's identity. The meshed arm solves
// Delta* psi = -mu0 r j_phi with the same datum and a source that is a top hat
// on the coil, so it approximates that same field to the mesh's order. So the
// difference between them IS the meshed route's discretisation error, and it
// must FALL under refinement. A difference that did not fall would say the two
// routes are solving different problems, which is exactly what §0b forbids.
BOOST_AUTO_TEST_CASE( theMeshedCoilConvergesToTheSubtractedOne )
{
	// A rectangle well inside the box, clear of its boundary, and ALIGNED TO
	// THE MESH AT EVERY LEVEL: the box is [ 0.6, 1.4 ] x [ -0.4, 0.4 ] and the
	// half-extents are 0.10, so the coil's edges sit at 0.90, 1.10, -0.10 and
	// +0.10 -- which are vertices at n = 8, 16 and 32 alike.
	//
	// THAT ALIGNMENT IS THE EXPERIMENT AND NOT A CONVENIENCE. The source is a
	// TOP HAT, so psi is not C^2 across the coil's edge; if the edge cuts
	// element interiors then WHICH elements it cuts is not a smooth function of
	// h, and the error is not even monotone. Measured, with half-extents of
	// 0.06 instead: 6.756e-04, 1.276e-03, 1.917e-04 over the same three levels
	// -- a rate of -0.92 and then +2.73. That is the unfitted-geometry
	// behaviour CLAUDE.md's "Unfitted convergence needs a two-tier rate
	// assertion" records for the extension path, met here for the same reason,
	// and it is a fact about cutting a discontinuity rather than about either
	// conductor route. MEQ meshes TO its coils -- tools/mesh/halfdisc.py exists
	// to do exactly that -- so the aligned case is the configuration this plan
	// is about.
	meq::Coil const conductor( 1.00, 0.00, 0.10, 0.10, 1.0e5 );

	meq::CoilSet meshed;
	meshed.add( conductor );

	meq::ConductorField subtracted;
	subtracted.add( conductor );

	// The coil as a DOMAIN SOURCE: F = mu0 r j_phi inside it and zero outside,
	// which is meq::CoilSet::f() and is what meq::CoilAugmentedSource wraps.
	struct MeshedCoilSource : public meq::Source
	{
		explicit MeshedCoilSource( meq::CoilSet const &c ) : coils( c ) {}

		double f( double r, double z, double ) const override
		{
			return coils.f( r, z );
		}

		double dFdPsi( double, double, double ) const override
		{
			return 0.0;
		}

		meq::CoilSet const &coils;
	};

	MeshedCoilSource const source( meshed );
	VacuumSource const vacuum;

	mfem::FunctionCoefficient datum(
		[ &meshed ]( mfem::Vector const &x )
		{
			return meshed.psi( x( 0 ), x( 1 ) );
		} );

	std::printf( "\n  §0a-pre: the same rectangle MESHED against SUBTRACTED\n" );
	std::printf( "    %8s %10s %14s %8s\n", "elements", "dofs",
	             "max |meshed - psi_c|", "rate" );

	double previous = 0.0;
	double worstRate = 1.0e30;

	for ( int n : { 8, 16, 32 } )
	{
		mfem::Mesh meshA = makeBox( n );
		meq::GradShafranovSolver meshedSolver( meshA, 2 );
		meshedSolver.setSource( source );
		meshedSolver.setBoundaryData( datum );
		meshedSolver.solve();

		// The subtracted arm on the same mesh, as the control that its own
		// remainder is still identically zero with a RECTANGLE rather than a
		// filament -- the identity must not depend on which kind it is.
		mfem::Mesh meshB = makeBox( n );
		meq::GradShafranovSolver splitSolver( meshB, 2 );
		splitSolver.setSource( vacuum );
		splitSolver.setBoundaryData( datum );
		splitSolver.setConductorField( subtracted );
		splitSolver.solve();

		BOOST_TEST( maxAbs( splitSolver.potential() ) < 1.0e-12,
		            "the rectangle's remainder is not zero" );

		// Compare the meshed solution against psi_c, at quadrature points
		// OUTSIDE the conductor: inside it the source is a top hat and the
		// meshed route's order is capped by the jump, which is a different
		// claim from the one being made here.
		double worst = 0.0;
		for ( int element = 0; element < meshA.GetNE(); ++element )
		{
			mfem::IntegrationRule const &rule = mfem::IntRules.Get(
				meshA.GetElementBaseGeometry( element ), 4 );
			mfem::IsoparametricTransformation transformation;
			meshA.GetElementTransformation( element, &transformation );

			for ( int q = 0; q < rule.GetNPoints(); ++q )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( q );
				transformation.SetIntPoint( &ip );
				mfem::Vector point( 3 );
				transformation.Transform( ip, point );

				// Outside the conductor, where both routes represent a smooth
				// field. Inside it the meshed one is resolving a top hat and
				// its order is capped by the jump, which is a different claim.
				if ( std::fabs( point( 0 ) - 1.00 ) < 0.13
				     && std::fabs( point( 1 ) - 0.00 ) < 0.13 )
					continue;

				double const got = meshedSolver.potential().GetValue( element,
				                                                      ip );
				worst = std::max( worst,
				                  std::fabs( got - subtracted.psi( point( 0 ),
				                                                   point( 1 ) ) ) );
			}
		}

		double const rate = previous > 0.0
			? std::log2( previous/worst ) : 0.0;
		if ( previous > 0.0 )
			worstRate = std::min( worstRate, rate );

		std::printf( "    %8d %10d %14.3e %8s", meshA.GetNE(),
		             meshedSolver.potential().Size(), worst,
		             previous > 0.0 ? "" : "   --\n" );
		if ( previous > 0.0 )
			std::printf( " %7.2f\n", rate );
		previous = worst;
	}

	// THE ASSERTION IS THAT IT FALLS, not that it reaches any particular order.
	// The source is a top hat, so psi is only C^1 across the coil's edge and no
	// polynomial degree recovers a clean k+1 globally -- which is the ORDER
	// argument and is not what §0a-pre needs. What it needs is that the two
	// routes converge to the same field.
	BOOST_TEST( worstRate > 0.9,
	            "the meshed coil is not converging to the subtracted one: "
	            "worst observed rate " << worstRate << ". The two routes are "
	            "solving different problems, which COIL-SUBTRACTION-PLAN.md "
	            "§0b forbids" );
}

// CS-4: A NORMALISED SOURCE UNDER THE SPLIT, WHICH THE REFUSAL USED TO FORBID.
//
// psi_ax IS AN UNKNOWN OF A BORDERED NEWTON AND A FUNCTIONAL OF THE PHYSICAL
// FLUX, so until every consumer of psi read psi_c + psi_p this combination was
// refused rather than approximated. The assertion here is the DEFINITION rather
// than a tolerance on an answer: the psi_ax the solver reports must be the total
// field at the axis it reports. If any consumer were still reading the
// remainder -- peakAt(), the fill, the limiter, the source -- the two would
// differ by psi_c, which on a machine is not a small number.
//
// AND THE CONTROL IS WHAT SAYS THE ASSERTION HAS TEETH: psi_c at that axis is
// reported beside it, so a reader can see the number the test would be wrong by.
BOOST_AUTO_TEST_CASE( aNormalisedSourceRunsUnderTheSplitAndReportsThePhysicalAxis )
{
	meq::ConductorField const conductors = insideFilament();

	mfem::Mesh mesh = makeBox( 16 );
	auto pPrime = std::make_shared< meq::ConstantProfile const >( 0.45 );
	auto ggPrime = std::make_shared< meq::ConstantProfile const >( 0.30 );
	meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );

	// THE PHYSICAL flux is zero on the box boundary; setBoundaryData() keeps
	// that meaning and the solver imposes psi|_Gamma - psi_c|_Gamma itself.
	mfem::ConstantCoefficient datum( 0.0 );

	// A separable sine bump, which is HighBetaConvergence's own guess for a
	// bordered psi_ax and is what puts a maximum inside the box for the axis
	// search to find. Without it the search has no interior extremum to locate
	// and reports ( 0, 0 ) -- which this case detected on its first run.
	mfem::FunctionCoefficient guess(
		[]( mfem::Vector const &x )
		{
			return 0.30*std::sin( M_PI*( x( 0 ) - 0.6 )/0.8 )
			       *std::sin( M_PI*( x( 1 ) + 0.4 )/0.8 );
		} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setBoundaryData( datum );
	solver.setInitialGuess( guess );

	// THE REFUSAL IS GONE: this call used to throw.
	BOOST_REQUIRE_NO_THROW( solver.setConductorField( conductors ) );
	BOOST_REQUIRE_NO_THROW( solver.setSource( source, 0.30 ) );

	solver.solve();

	double const reported = solver.psiAxis();

	/*
	 * psi_ax's CONSTRAINT IS psi_ax - max psi = 0, SO THE TEST IS THAT MAXIMUM.
	 * This configuration constrains against the NODAL PEAK rather than a
	 * located critical point -- axisR() is only filled by the located-axis
	 * mode, and reads ( 0, 0 ) here, which the first version of this case
	 * asserted against and correctly refused to pass on.
	 *
	 * So the peak is recomputed here from outside the solver, over the total
	 * and over the remainder alone, and psi_ax must match the FORMER. That is
	 * exactly what peakAt() was changed to do, and the two differ by psi_c --
	 * which the control below asserts is not a small number, or the case could
	 * not tell them apart.
	 */
	mfem::GridFunction const &remainder = solver.potential();
	mfem::FiniteElementSpace const *space = remainder.FESpace();

	double peakTotal = -std::numeric_limits< double >::infinity();
	double peakRemainder = -std::numeric_limits< double >::infinity();

	mfem::Array< int > dofs;
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		space->GetElementDofs( e, dofs );
		mfem::FiniteElement const *fe = space->GetFE( e );
		mfem::IntegrationRule const &nodes = fe->GetNodes();

		mfem::IsoparametricTransformation transformation;
		mesh.GetElementTransformation( e, &transformation );

		for ( int i = 0; i < dofs.Size(); ++i )
		{
			mfem::IntegrationPoint const &ip = nodes.IntPoint( i );
			transformation.SetIntPoint( &ip );
			mfem::Vector point( 3 );
			transformation.Transform( ip, point );

			int const dof = dofs[ i ] >= 0 ? dofs[ i ] : -1 - dofs[ i ];
			double const solved = remainder( dof );
			peakRemainder = std::max( peakRemainder, solved );
			peakTotal = std::max( peakTotal,
			                      solved + conductors.psi( point( 0 ),
			                                               point( 1 ) ) );
		}
	}

	std::printf( "\n  CS-4: a normalised source under the split, %d Newton steps\n",
	             solver.newtonIterations() );
	std::printf( "    psi_ax reported %12.6e\n", reported );
	std::printf( "    peak of psi_c + psi_p %12.6e   <- what it must be\n",
	             peakTotal );
	std::printf( "    peak of psi_p alone   %12.6e   <- what it would be if "
	             "peakAt() still read the remainder\n", peakRemainder );

	BOOST_TEST_REQUIRE( std::isfinite( reported ) );

	// The control: the two candidates must be far apart, or matching one of
	// them says nothing.
	BOOST_TEST( std::fabs( peakTotal - peakRemainder ) > 1.0e-3,
	            "the total and the remainder peak at nearly the same value, so "
	            "this case cannot tell which one psi_ax followed" );

	BOOST_TEST( std::fabs( reported - peakTotal ) < 1.0e-6,
	            "psi_ax is not the peak of the PHYSICAL flux" );
}

// psi_c IS EXACTLY ZERO WHEN NO CONDUCTOR FIELD IS SET, which is what lets
// CS-4's consumers add it unconditionally instead of branching.
BOOST_AUTO_TEST_CASE( theConductorPsiIsZeroWithoutASplit )
{
	mfem::Mesh mesh = makeBox( 4 );
	meq::GradShafranovSolver solver( mesh, 1 );

	BOOST_TEST( solver.conductorField() == nullptr );
	BOOST_TEST( solver.conductorPsi( 1.0, 0.1 ) == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );
}
