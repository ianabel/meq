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
#include <vector>

#include "meq/ConductorField.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/FieldViews.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/SurfaceAverage.hpp"

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

// CS-4 TRANCHE TWO: A FLUX SURFACE IS A LEVEL SET OF THE TOTAL, AND TRACING
// THE REMAINDER RETURNS A PERFECTLY GOOD CURVE THAT IS NOT ONE.
//
// This is the failure the stage exists to close and it has no symptom of its
// own. meq::ContourTracer roots the field it was given, so under the split it
// would root psi_p -- converge, close, report a small corrector residual and a
// clean turning number, and hand meq::surfaceAverages a curve that is not a
// flux surface. V', q, the metric and every column of _surfaces.nc would then
// be an average over the wrong curve, with nothing in the file saying so.
//
// THE TWO COLUMNS ARE THE WHOLE CASE. Along the traced contour:
//
//     spread of psi_p + psi_c    must be at the corrector's own tolerance
//     spread of psi_p alone      must be LARGE, or the case cannot tell them
//                                apart and proves nothing
//
// and the second is the control in exactly the sense CS-2's harmonic-extension
// arm is: it is what the answer would look like if the shift were absent.
//
// The evaluation is INDEPENDENT of the seam under test. A second tracer is
// built on the same two GridFunctions through the bare-field constructor, which
// cannot know about conductors, so it returns psi_p; the test adds psi_c itself
// from meq::ConductorField. What is shared is the element walk, which is not
// what changed.
BOOST_AUTO_TEST_CASE( theTracedSurfaceIsALevelSetOfTheTotalAndNotTheRemainder )
{
	meq::ConductorField const conductors = insideFilament();

	mfem::Mesh mesh = makeBox( 24 );
	auto pPrime = std::make_shared< meq::ConstantProfile const >( 0.45 );
	auto ggPrime = std::make_shared< meq::ConstantProfile const >( 0.30 );
	meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );

	// The PHYSICAL flux is zero on the box boundary; the solver subtracts
	// psi_c|_Gamma itself. Same configuration as the case above.
	mfem::ConstantCoefficient datum( 0.0 );
	mfem::FunctionCoefficient guess(
		[]( mfem::Vector const &x )
		{
			return 0.30*std::sin( M_PI*( x( 0 ) - 0.6 )/0.8 )
			       *std::sin( M_PI*( x( 1 ) + 0.4 )/0.8 );
		} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setBoundaryData( datum );
	solver.setInitialGuess( guess );
	solver.setConductorField( conductors );
	solver.setSource( source, 0.30 );
	solver.solve();
	solver.postProcess();

	// THE TRACER TAKES THE CONDUCTORS FROM THE SOLVER AND THE CALLER NEVER
	// ASKS. That is the point of doing it in the constructor: apps/meq.cpp
	// builds its tracer from the solver and needed no edit at all.
	meq::ContourTracer const tracer( solver );
	BOOST_TEST_REQUIRE( tracer.conductorField() == &conductors );

	// The control tracer: the same two fields through the door that cannot
	// know about a split, so it answers psi_p and q_p.
	meq::ContourTracer const bare( solver.postProcessedPotential(),
	                               solver.postProcessedFlux() );
	BOOST_TEST_REQUIRE( bare.conductorField() == nullptr );

	// THE SEAM IDENTITY, AT ZERO TOLERANCE, BEFORE ANY TRACING. Whatever the
	// geometry turns out to be, sampleAt() must be the bare answer plus psi_c
	// at that point and nothing else -- so a failure of the case below is a
	// failure of the tracing and not of the shift.
	double worstSeam = 0.0;
	for ( int i = 0; i < 7; ++i )
		for ( int j = 0; j < 7; ++j )
		{
			double const r = 0.70 + 0.10*i;
			double const z = -0.30 + 0.10*j;

			double psiTotal = 0.0, qRTotal = 0.0, qZTotal = 0.0;
			double psiBare = 0.0, qRBare = 0.0, qZBare = 0.0;
			if ( !tracer.sampleAt( r, z, psiTotal, qRTotal, qZTotal ) )
				continue;
			BOOST_TEST_REQUIRE( bare.sampleAt( r, z, psiBare, qRBare, qZBare ) );

			double bR = 0.0, bZ = 0.0;
			conductors.poloidalField( r, z, bR, bZ );

			worstSeam = std::max( worstSeam,
				std::fabs( psiTotal - ( psiBare + conductors.psi( r, z ) ) ) );
			worstSeam = std::max( worstSeam,
			                      std::fabs( qRTotal - ( qRBare + bZ ) ) );
			worstSeam = std::max( worstSeam,
			                      std::fabs( qZTotal - ( qZBare - bR ) ) );
		}

	/*
	 * THE LEVEL IS READ OFF THE FIELD RATHER THAN CHOSEN, because what value
	 * the physical flux takes at a given point is a property of this solve and
	 * not something a test may assume. The start point is on the far side of
	 * the box from the filament, so the contour through it encloses both the
	 * plasma's own maximum and the conductor.
	 */
	double level = 0.0, qRAt = 0.0, qZAt = 0.0;
	BOOST_TEST_REQUIRE( tracer.sampleAt( 0.78, 0.00, level, qRAt, qZAt ) );

	meq::Contour const contour = tracer.trace( level, 0.78, 0.00 );

	double worstTotal = 0.0;
	double lowRemainder = std::numeric_limits< double >::infinity();
	double highRemainder = -std::numeric_limits< double >::infinity();

	for ( std::size_t i = 0; i < contour.points.size(); ++i )
	{
		meq::ContourPoint const &point = contour.points[ i ];

		double psiBare = 0.0, qR = 0.0, qZ = 0.0;
		BOOST_TEST_REQUIRE( bare.sampleAt( point.r, point.z, psiBare, qR, qZ ) );

		double const total = psiBare + conductors.psi( point.r, point.z );
		worstTotal = std::max( worstTotal, std::fabs( total - level ) );
		lowRemainder = std::min( lowRemainder, psiBare );
		highRemainder = std::max( highRemainder, psiBare );
	}

	double const spreadRemainder = highRemainder - lowRemainder;

	/*
	 * THE FALSIFYING ARM, AND IT IS THE CURVE THIS STAGE REPLACES RATHER THAN
	 * AN ARTIFICIAL ONE.
	 *
	 * Ask the conductor-blind tracer for "the surface through ( 0.78, 0 )" --
	 * which is the trace of psi_p at the value psi_p takes there, and is
	 * exactly what apps/meq.cpp wrote into _surfaces.nc before this change.
	 * It closes, it reports a corrector residual at the same tolerance, and
	 * nothing about it says it is wrong. What says so is the PHYSICAL flux
	 * along it, measured here and reported beside the real surface's.
	 */
	double psiAtStart = 0.0, qRStart = 0.0, qZStart = 0.0;
	BOOST_TEST_REQUIRE( bare.sampleAt( 0.78, 0.00, psiAtStart, qRStart,
	                                   qZStart ) );
	meq::Contour const blind = bare.trace( psiAtStart, 0.78, 0.00 );

	double lowBlind = std::numeric_limits< double >::infinity();
	double highBlind = -std::numeric_limits< double >::infinity();
	for ( std::size_t i = 0; i < blind.points.size(); ++i )
	{
		meq::ContourPoint const &point = blind.points[ i ];

		double psiBare = 0.0, qR = 0.0, qZ = 0.0;
		BOOST_TEST_REQUIRE( bare.sampleAt( point.r, point.z, psiBare, qR,
		                                   qZ ) );

		double const total = psiBare + conductors.psi( point.r, point.z );
		lowBlind = std::min( lowBlind, total );
		highBlind = std::max( highBlind, total );
	}
	double const spreadBlind = blind.points.empty()
		? 0.0 : highBlind - lowBlind;

	std::printf( "\n  CS-4: the traced surface under the split\n" );
	std::printf( "    seam, sampleAt against psi_p + psi_c    %11.4e   (exact)\n",
	             worstSeam );
	std::printf( "    level traced                            %11.4e\n", level );
	std::printf( "    points %zu, closed %d, corrector target %11.4e\n",
	             contour.points.size(), contour.closed() ? 1 : 0,
	             contour.correctorTarget );
	std::printf( "    worst | psi_p + psi_c - level |         %11.4e   <- the "
	             "surface\n", worstTotal );
	std::printf( "    spread of psi_p alone along it          %11.4e   <- what "
	             "it is NOT a surface of\n", spreadRemainder );
	std::printf( "    the conductor-blind curve: %zu points, closed %d\n",
	             blind.points.size(), blind.closed() ? 1 : 0 );
	std::printf( "    spread of psi_p + psi_c along THAT      %11.4e   <- the "
	             "defect, had it stood\n", spreadBlind );

	BOOST_TEST( worstSeam == 0.0,
	            "sampleAt() is not exactly the solved field plus psi_c, so the "
	            "shift is being applied somewhere other than at the seam" );

	BOOST_TEST_REQUIRE( contour.points.size() > 20u );
	BOOST_TEST( contour.closed(),
	            "the contour did not close, so this case is measuring a trace "
	            "that failed rather than the field it traced" );

	// THE CONTROL FIRST: without a wide separation the assertion below is
	// empty, because psi_p would then be nearly constant on this curve anyway.
	BOOST_TEST( spreadRemainder > 1.0e-3,
	            "psi_p barely varies along the traced contour, so this case "
	            "cannot tell a level set of the total from one of the remainder" );

	BOOST_TEST( worstTotal < 1.0e-4*spreadRemainder,
	            "the traced curve is not a level set of psi_p + psi_c, which is "
	            "what a flux surface is under COIL-SUBTRACTION-PLAN.md's split" );

	// AND THE ARM THAT SAYS THE DEFECT WAS REAL. A curve that closed cleanly
	// on the remainder carries a physical flux that is not constant along it
	// at all, by four orders or more over what the real surface manages.
	BOOST_TEST_REQUIRE( blind.points.size() > 20u );
	BOOST_TEST( spreadBlind > 1.0e4*worstTotal,
	            "the conductor-blind curve is very nearly a flux surface after "
	            "all, so this fixture does not exhibit the defect the shift was "
	            "added to close and the case above proves nothing" );
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

/*
 * CS-3: THE SPLIT UNDER AN EXTERIOR COUPLING, AND IT IS FB-7's CASE TURNED
 * INSIDE OUT.
 *
 * FB-7 puts a conductor OUTSIDE Gamma and asks the coupling to carry it: the
 * exterior field is entirely the conductor's, so the continuous answer is
 * a = 0 and `aConductorOutsideGammaReachesTheCoupledSolve` asserts on that.
 * CS-3 puts one INSIDE Gamma and subtracts it, and every part of that statement
 * inverts:
 *
 *     the interior      psi_p solves Delta* psi_p = 0, so with no plasma the
 *                       continuous answer is psi_p IDENTICALLY ZERO -- psi is
 *                       psi_c and nothing else
 *     the exterior      `a` is NOT zero: it is psi_c's own Gegenbauer trace on
 *                       Gamma, which is the whole datum rather than an error
 *
 * TWO THINGS HAD TO MOVE TOGETHER FOR THAT TO HOLD, AND EITHER ONE ALONE IS A
 * PLAUSIBLE WRONG ANSWER -- COIL-SUBTRACTION-PLAN.md section 3 says so in
 * advance:
 *
 *     the DIRICHLET half   prepare() transfers g - psi_c inward, because what
 *                          lives on Gamma is the remainder
 *     the NEUMANN half     transmissionConstraint() ADDS the subtracted
 *                          conductors' moment, because row . state is now
 *                          int ( q_p . nu ) C_m and the condition is about
 *                          q_p + q_c
 *
 * AND THE ASSERTIONS ARE CHOSEN SO THAT EACH FAILURE IS CAUGHT SEPARATELY.
 * Leave the datum unshifted and psi_p is forced to carry the whole modal trace
 * on Gamma while its own equation says it is harmonic -- O( psi_c ), not small.
 * Get the moment's SIGN wrong and `a` is off by twice the conductor's flux,
 * which does not fall with h; the same disguise FB-7 records, where a wrong
 * sign fails to converge rather than diverging. Drop the moment entirely and
 * `a` is wrong by one conductor field. All three are caught by the first
 * assertion, and the second -- the modal sum reproducing psi_c ON Gamma -- is
 * the sharp one, because it reads the Dirichlet half directly rather than
 * through what the interior did with it.
 */
namespace
{
	/// Gamma at 1.5, the box at 1.7: FreeBoundaryCoupling.cpp's half-disc,
	/// copied rather than shared.
	///
	/// COPIED DELIBERATELY. That file is now the whole of the ctest suite --
	/// CLAUDE.md records it at 1149.60 s and 991.99 s on two runs of one tree,
	/// against `naming` at about 350 s -- so a case added there costs the suite
	/// its own runtime plus that binary's link. This one runs in seconds beside
	/// the rest of CS-2's cases, and the duplication is fifty lines of fixture.
	double const gammaRadius = 1.5;
	double const boxRadius = 1.7;

	double gammaLevelSet( mfem::Vector const &x )
	{
		return std::hypot( x( 0 ), x( 1 ) ) - gammaRadius;
	}

	/// The background, the cut subdomain and its transfer path, held together
	/// because a mfem::SubMesh KEEPS A POINTER TO ITS PARENT -- CLAUDE.md's
	/// trap, which cost three fixtures a segfault in mfem::VertexConePath's
	/// constructor. Returning the pair is what makes the parent outlive it.
	struct HalfDisc
	{
		std::unique_ptr<mfem::Mesh> background;
		std::unique_ptr<mfem::SubMesh> sub;
		std::unique_ptr<mfem::VertexConePath> path;
		mfem::Array<int> gammaHMarker;
		double h = 0.0;
	};

	HalfDisc makeHalfDisc( int n )
	{
		HalfDisc d;
		d.background = std::make_unique<mfem::Mesh>(
			mfem::Mesh::MakeCartesian2D( n, 2*n, mfem::Element::TRIANGLE, false,
			                             boxRadius, 2.0*boxRadius ) );
		d.background->Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 );
			out( 1 ) = in( 1 ) - boxRadius;
		} );
		d.h = boxRadius/static_cast<double>( n );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE( mfem::MarkLevelSetSubdomain(
			*d.background, gammaLevelSet, 0.0, marker, 1 ) > 0,
			"the half-disc is empty at n = " << n );
		for ( int e = 0; e < d.background->GetNE(); ++e )
			d.background->SetAttribute( e, marker[ e ] ? 1 : 2 );
		d.background->SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		d.sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( *d.background, domainAttr ) );

		// The arc is GENERATED by SubMesh and takes the new attribute; the flat
		// side is INHERITED from the box's r = 0 edge. So Gamma_h is the arc
		// alone and the axis stays ordinary fitted boundary, which is right --
		// the axis approximates nothing and needs no transfer.
		int const gammaH = d.sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( d.sub->bdr_attributes.Size() >= 2,
			"D_h has one boundary attribute at n = " << n
			<< ", so Gamma_h has swallowed the axis" );
		d.gammaHMarker.SetSize( gammaH );
		d.gammaHMarker = 0;
		d.gammaHMarker[ gammaH - 1 ] = 1;

		d.path = std::make_unique<mfem::VertexConePath>(
			*d.sub, gammaH, gammaLevelSet, 6.0*d.h );
		return d;
	}

	/// An up-down symmetric pair WELL inside Gamma, and the margin is the
	/// experiment's precision rather than tidiness.
	///
	/// The exterior is a Gegenbauer series in rho^( 1 - n ), truncated at
	/// `modes` terms, so what it cannot represent of psi_c falls like
	/// ( reach/gammaRadius )^n. The farthest point of these rectangles is
	/// hypot( 0.55, 0.25 ) = 0.604, which is 0.40 of Gamma -- so twelve modes
	/// leave about 1e-5 and the truncation is below the discretisation rather
	/// than beside it. Put the conductors at 0.9 instead and the ratio is 0.73,
	/// the residual is 15%, and psi_p would carry the SERIES rather than the
	/// solve.
	meq::ConductorField insideGammaPair()
	{
		meq::ConductorField field;
		field.add( meq::Coil( 0.50, +0.20, 0.05, 0.05, 6.0e5 ) );
		field.add( meq::Coil( 0.50, -0.20, 0.05, 0.05, 6.0e5 ) );
		return field;
	}
}

BOOST_AUTO_TEST_CASE( theSplitReachesTheExteriorCoupling )
{
	int const order = 2;
	int const modes = 12;

	meq::ConductorField const conductors = insideGammaPair();
	// RECTANGLES AND NOT FILAMENTS, so that the mesh-coincidence refusal is not
	// what this case is about: a rectangle has no line singularity for a node
	// to land on, which is CS-1b's contract.
	BOOST_TEST_REQUIRE( conductors.containment( 0.0, gammaRadius ) > 0.0 );

	VacuumSource noPlasma;

	// psi_c's own scale on Gamma -- what `a` must reproduce and what psi_p is
	// measured against. A bound relative to 1 would be a statement about the
	// coil current.
	double const datumScale = std::abs( conductors.psi( gammaRadius, 0.0 ) );
	BOOST_TEST_REQUIRE( datumScale > 0.0 );

	std::printf( "\n  CS-3: THE SPLIT UNDER AN EXTERIOR COUPLING\n" );
	std::printf( "    no plasma and the conductors INSIDE Gamma, so psi = psi_c "
	             "and the remainder is zero;\n"
	             "    psi_c on Gamma is %.4e and `a` must reproduce it\n\n",
	             datumScale );
	std::printf( "    %5s %7s %14s %14s\n",
	             "n", "newton", "max |psi_p|", "modal - psi_c" );

	std::vector<double> remainders;
	double worstTrace = 0.0;

	for ( int n : { 12, 24, 48 } )
	{
		HalfDisc d = makeHalfDisc( n );
		meq::ExteriorDtN const dtn( 0.0, gammaRadius, modes );
		mfem::ConstantCoefficient zero( 0.0 );

		meq::GradShafranovSolver solver( *d.sub, order );
		solver.setSource( noPlasma );
		solver.setBoundaryData( zero );
		solver.setExtension( *d.path, d.gammaHMarker );
		solver.setConductorField( conductors );
		solver.setExteriorCoupling( dtn );
		solver.solve();

		// (1) THE REMAINDER. psi_p is what the solver holds; the physical flux
		// is psi_p + psi_c, and the claim is that the second term is all of it.
		double const remainder = maxAbs( solver.potential() );

		// (2) THE MODAL SUM ON Gamma, WHICH IS THE SHARP ONE. It reads the
		// Dirichlet half directly: `a` is determined by psi_p = g - psi_c on
		// Gamma, so if the shift were absent the modal sum would reproduce
		// something else entirely. Sampled around the arc rather than at one
		// point, because a single point can agree by accident -- the modes are
		// orthogonal, not independent, in the value at one place.
		std::vector<double> const &a = solver.exteriorCoefficients();
		BOOST_TEST_REQUIRE( static_cast<int>( a.size() ) == modes );
		double trace = 0.0;
		double biggest = 0.0;
		for ( double const angle : { 0.05, 0.4, 0.9, 1.4, 1.9, 2.4, 3.0 } )
		{
			double const r = gammaRadius*std::sin( angle );
			double const z = gammaRadius*std::cos( angle );
			double modal = 0.0;
			for ( std::size_t m = 0; m < a.size(); ++m )
				modal += a[ m ]*dtn.basis(
					meq::ExteriorDtN::firstMode() + static_cast<int>( m ), r, z );
			trace = std::max( trace, std::abs( modal - conductors.psi( r, z ) ) );
			biggest = std::max( biggest, std::abs( modal ) );
		}
		BOOST_TEST_REQUIRE( biggest > 0.0,
			"the exterior coefficients are identically zero, so the coupling "
			"carried no datum at all -- which is FB-7's answer and not this "
			"case's: psi_c is INSIDE Gamma here and the exterior has to "
			"represent it" );
		worstTrace = std::max( worstTrace, trace );

		remainders.push_back( remainder );

		std::printf( "    %5d %7d %14.4e %14.4e\n",
		             n, solver.newtonIterations(), remainder, trace );
		std::fflush( stdout );

		// ONE NEWTON STEP. With no plasma the residual is affine in ( x, a ) and
		// psi_c enters as a constant, so more than one step would mean the
		// conductor had got into the Jacobian.
		BOOST_TEST( solver.newtonIterations() <= 1,
			"the coupled vacuum solve under the split took "
			<< solver.newtonIterations() << " Newton steps where the residual is "
			"affine" );
	}

	std::printf( "\n    worst modal-minus-psi_c on Gamma: %.4e   "
	             "( %.2e of the datum )\n\n",
	             worstTrace, worstTrace/datumScale );
	std::fflush( stdout );

	/*
	 * AND IT IS AN IDENTITY RATHER THAN A RATE, WHICH IS WHY THERE IS NO RATE
	 * COLUMN ABOVE AND WHY THE GATE IS 1e-5 RATHER THAN SOMETHING LOOSE.
	 *
	 * Measured: 8.80e-09, 1.43e-08, 1.47e-08 at n = 12, 24, 48 against a datum
	 * of 6.38e-02 -- FLAT, and 2.3e-07 of the field it is a remainder from.
	 * That is the right shape and not a disappointment: psi_p is identically
	 * zero in the continuum and the discrete problem has nothing to converge,
	 * because the remainder equation's right-hand side is zero and its boundary
	 * data is zero to the truncation. What is left is the twelve-mode series
	 * residual plus the linear solve's round-off, and neither depends on h --
	 * the `modal - psi_c` column reads 1.20e-08, 1.21e-08, 1.20e-08 and says
	 * so directly.
	 *
	 * SO A RATE ASSERTION HERE WOULD BE WRONG IN BOTH DIRECTIONS: it would fail
	 * on correct code (there is no rate) and it would pass on code that was
	 * merely converging towards the right answer from somewhere large. CS-2's
	 * theRemainderOfAVacuumFilamentIsZero makes the same choice for the same
	 * reason, and section 7.4 of the plan is the argument.
	 */
	/*
	 * THE REMAINDER IS SMALL AGAINST THE FIELD IT IS A REMAINDER FROM, AND
	 * BOTH HALVES OF CS-3 WERE FALSIFIED AGAINST IT RATHER THAN ARGUED FOR.
	 *
	 * Each was broken in turn, rebuilt, and this case re-run at n = 48 --
	 * MEASUREMENTS.md M-143 is the published table:
	 *
	 *     as built                             1.47e-08    2.3e-07 of the datum
	 *     interior moment's SIGN flipped       4.20e-02    6.6e-01 of it
	 *     datum shift removed from prepare()   2.10e-02    3.3e-01 of it
	 *
	 * -- and in both broken arms the modal sum on Gamma reads 4.12e-02 against
	 * psi_c, 65% of the datum, so the second assertion goes red too. That is
	 * what says the two halves are separately load bearing rather than one fix
	 * written twice, which is exactly the trap COIL-SUBTRACTION-PLAN.md section
	 * 3 predicts for this pair: "the two must move together or the datum is
	 * double counted".
	 *
	 * AND NEITHER BROKEN ARM DIVERGES OR FAILS TO CONVERGE. Both take one
	 * Newton step and report a smooth, plausible field -- FB-7 records the same
	 * disguise for the transmission row's own sign. A residual check could not
	 * tell; only comparing against the closed form can.
	 */
	BOOST_TEST( remainders.back() < 1.0e-5*datumScale,
		"the solved remainder is " << remainders.back() << " against a psi_c of "
		<< datumScale << " on Gamma -- " << remainders.back()/datumScale
		<< " of it, where the continuous answer is exactly zero. CS-3's two "
		"halves are the datum shift in prepare() and the interior conductor "
		"moment in transmissionConstraint(), and either one absent leaves this "
		"of order one" );

	// AND THE EXTERIOR REPRODUCES psi_c ON Gamma, which is the Dirichlet half
	// on its own. The bound is the series truncation at twelve modes, not the
	// mesh: see insideGammaPair() for the 0.40 ratio that sets it.
	BOOST_TEST( worstTrace < 1.0e-5*datumScale,
		"the modal sum differs from psi_c on Gamma by " << worstTrace
		<< " against a datum of " << datumScale << ". With no plasma the "
		"exterior field IS the conductors', so `a` is psi_c's Gegenbauer trace "
		"and nothing else" );
}

/*
 * AND THE GEOMETRY IS REFUSED IN BOTH ORDERS, WHICH IS THE ONLY WAY A
 * PRECONDITION ABOUT A PAIR CAN BE ENFORCED BY TWO INDEPENDENT SETTERS.
 *
 * A SUBTRACTED conductor must be strictly INSIDE Gamma and an exterior one
 * strictly OUTSIDE it, and the two refusals are mirror images with dual
 * reasons: the exterior expansion is a series in rho^( 1 - n ) standing for the
 * sources within Gamma, so a subtracted conductor beyond it puts a singularity
 * in the region that series describes -- and nothing downstream could say so.
 * The transmission would close, every border would reach machine zero, and the
 * answer would be a machine nobody described.
 */
BOOST_AUTO_TEST_CASE( aSubtractedConductorOutsideGammaIsRefusedEitherWayRound )
{
	HalfDisc d = makeHalfDisc( 12 );
	meq::ExteriorDtN const dtn( 0.0, gammaRadius, 4 );
	mfem::ConstantCoefficient zero( 0.0 );
	VacuumSource noPlasma;

	// BEYOND Gamma: hypot( 1.75, 0.10 ) against a radius of 1.5.
	meq::ConductorField beyond;
	beyond.add( meq::Coil( 1.70, 0.00, 0.05, 0.05, 1.0e5 ) );
	BOOST_TEST_REQUIRE( beyond.containment( 0.0, gammaRadius ) < 0.0 );

	// STRADDLING it, which is the case that belongs to neither route and is why
	// the measure is to the FARTHEST point rather than to the centre.
	meq::ConductorField straddling;
	straddling.add( meq::Coil( 1.45, 0.00, 0.10, 0.10, 1.0e5 ) );
	BOOST_TEST_REQUIRE( straddling.containment( 0.0, gammaRadius ) < 0.0 );

	int const order = 2;
	auto build = [ & ]()
	{
		auto solver = std::make_unique<meq::GradShafranovSolver>( *d.sub, order );
		solver->setSource( noPlasma );
		solver->setBoundaryData( zero );
		solver->setExtension( *d.path, d.gammaHMarker );
		return solver;
	};

	// COUPLING FIRST, THEN THE CONDUCTORS.
	for ( meq::ConductorField const *bad : { &beyond, &straddling } )
	{
		auto solver = build();
		solver->setExteriorCoupling( dtn );
		BOOST_CHECK_THROW( solver->setConductorField( *bad ),
		                   std::invalid_argument );
	}

	// CONDUCTORS FIRST, THEN THE COUPLING -- the other order, and the refusal
	// has to come from the other end.
	for ( meq::ConductorField const *bad : { &beyond, &straddling } )
	{
		auto solver = build();
		solver->setConductorField( *bad );
		BOOST_CHECK_THROW( solver->setExteriorCoupling( dtn ),
		                   std::invalid_argument );
	}

	// AND A CONDUCTOR INSIDE IS ACCEPTED IN BOTH ORDERS, or the refusals above
	// could be passing because the pair is refused outright.
	meq::ConductorField const good = insideGammaPair();
	{
		auto solver = build();
		solver->setExteriorCoupling( dtn );
		BOOST_CHECK_NO_THROW( solver->setConductorField( good ) );
	}
	{
		auto solver = build();
		solver->setConductorField( good );
		BOOST_CHECK_NO_THROW( solver->setExteriorCoupling( dtn ) );
	}

	// WITHOUT A COUPLING THERE IS NO Gamma AND NOTHING TO BE OUTSIDE OF, which
	// is the fixed-boundary case CS-2 is about: a subtracted conductor may sit
	// anywhere at all there, and refusing one would be a rule invented rather
	// than derived.
	{
		auto solver = build();
		BOOST_CHECK_NO_THROW( solver->setConductorField( beyond ) );
	}
}

// CS-5: MEQ'S OWN EXPENSIVE QUADRATURE IS THE REFERENCE, BECAUSE THE BENCHMARK
// CANNOT BE ONE.
//
// COIL-SUBTRACTION-PLAN.md §7.2 retires CS-5 as staged and says why:
// MEASUREMENTS.md M-139 puts the freegs4e comparison's floor at 7.7e-04 and it
// does NOT refine, because freegs4e's ShapedCoil caps at six points per
// triangle. A change whose entire claim is that the conductors are resolved
// EXACTLY cannot be accepted against a reference that resolves them to six
// points.
//
// The replacement is a self-comparison, and it is the right one here because
// the claim is about REPRESENTATION rather than about physics: MEQ owns the
// kernel -- meq::coilPsi() at whatever order is asked for, up to
// meq::maximumCoilQuadratureOrder -- so a reference field can be built far
// beyond anything a solve would use, on the SAME conductors, and the shipped
// rule measured against it.
//
// WHAT THIS CASE IS FOR, IN ONE LINE: to say how much of the split's error is
// the quadrature, so that everything else measured about it can be attributed
// elsewhere. It is the instrument that M-139 denied CS-5 as staged.
//
// AND IT IS SPLIT BY WHERE THE FIELD POINT IS, because the two regimes
// converge differently and merging them would hide the slower one. OFF the
// conductor the integrand is analytic and Gauss-Legendre is geometric in the
// separation; ON or INSIDE it the kernel carries a logarithm, the rule's cubic
// grading turns that into algebraic convergence, and no order in the shipped
// range reaches round-off. That is a property of the integrand and not a defect
// -- and it does not reach a subtracted solve, whose whole premise is that the
// conductors are OUTSIDE the elements being integrated.
BOOST_AUTO_TEST_CASE( theConductorQuadratureIsBoundedByAnExpensiveReference )
{
	// FAR BEYOND A SOLVE AND FAR BELOW THE CAP. meq::maximumCoilQuadratureOrder
	// is 256 and requireOrder()'s message records that the rule reaches
	// round-off at about 48, so 160 is a reference in the sense that matters:
	// raising it further changes nothing this case can see.
	int const reference = 160;
	int const shipped = meq::defaultCoilQuadratureOrder;

	// FOUR SHAPES AND THE ASPECT RATIO IS THE VARIABLE. A machine's conductors
	// are not all compact: MAST-U's solenoid is 12 mm by 3.18 m, which is 265
	// to 1, and it is the hardest case for any rule graded on a single length.
	struct Shape
	{
		char const *name;
		meq::Coil coil;
	};

	std::vector< Shape > const shapes = {
		{ "square 0.1",  meq::Coil( 1.0000, 0.000, 0.0500, 0.0500, 1.0e5 ) },
		{ "D3, 1.5:1",   meq::Coil( 0.8065, 1.982, 0.0431, 0.0284, 2.43e5 ) },
		{ "PX1, 16:1",   meq::Coil( 0.2415, 1.224, 0.0125, 0.2014, 1.55e5 ) },
		{ "solenoid",    meq::Coil( 0.1948, 0.000, 0.0060, 1.5900, 1.62e6 ) },
	};

	double const mu0 = meq::vacuumPermeability;

	std::printf( "\n  CS-5: the shipped order %d against an order %d reference\n",
	             shipped, reference );
	std::printf( "    %-12s %10s %12s %12s %12s\n", "coil", "gap/size",
	             "psi rel", "|grad| rel", "order 8 psi" );

	double worstClear = 0.0;
	double worstSurfacePsi = 0.0;
	double worstSurfaceGrad = 0.0;
	double worstInside = 0.0;
	double leastGain = 1.0e30;

	for ( Shape const &shape : shapes )
	{
		meq::Coil const &c = shape.coil;
		// THE COIL'S OWN SIZE, not its smaller half-extent. The far field is
		// set by the whole conductor, so this is the length the separation is
		// similar in -- measured across all four shapes, which is what four
		// shapes are here for.
		double const size = std::max( c.halfWidth(), c.halfHeight() );

		// TWO DIRECTIONS AND THE SECOND IS THE HARD ONE. Out past the CORNER
		// the point is far from most of the conductor whatever its aspect
		// ratio; out from the middle of the LONG FACE it is a coil-size away
		// from the near edge and a whole half-length from the far one, which
		// is where a rule graded on one length is weakest. On a 265:1 solenoid
		// the two differ by orders of magnitude and quoting only the corner
		// would understate the bound.
		for ( int direction = 0; direction < 2; ++direction )
		for ( double ratio : { -1.0, 0.0, 0.1, 0.5, 2.0, 10.0 } )
		{
			// ratio < 0 is the coil's CENTRE, which is the worst point of all:
			// the field point is in the middle of the domain being integrated.
			// It is the same point in both directions, so it is taken once.
			bool const inside = ratio < 0.0;
			if ( inside && direction == 1 )
				continue;

			double const r = inside ? c.centreR()
			               : direction == 0 ? c.rMax() + ratio*size
			                                : c.rMax() + ratio*size;
			double const z = inside ? c.centreZ()
			               : direction == 0 ? c.zMax() + ratio*size
			                                : c.centreZ();

			double const exact = meq::coilPsi( c, r, z, reference, mu0 );
			double eR = 0.0;
			double eZ = 0.0;
			meq::coilGradPsi( c, r, z, eR, eZ, reference, mu0 );
			double const gradScale = std::hypot( eR, eZ );

			double const got = meq::coilPsi( c, r, z, shipped, mu0 );
			double aR = 0.0;
			double aZ = 0.0;
			meq::coilGradPsi( c, r, z, aR, aZ, shipped, mu0 );

			double const coarse = meq::coilPsi( c, r, z, 8, mu0 );

			double const psiError = std::fabs( got - exact )
			                        /std::max( std::fabs( exact ), 1.0e-300 );
			double const gradError = std::hypot( aR - eR, aZ - eZ )
			                         /std::max( gradScale, 1.0e-300 );
			double const coarseError = std::fabs( coarse - exact )
			                           /std::max( std::fabs( exact ), 1.0e-300 );

			char where[ 32 ];
			std::snprintf( where, sizeof( where ), "%s%s",
			               inside ? "centre" :
			               ratio == 0.0 ? "0" :
			               ratio == 0.1 ? "0.1" :
			               ratio == 0.5 ? "0.5" :
			               ratio == 2.0 ? "2" : "10",
			               inside ? "" : direction == 0 ? " corner" : " face" );
			std::printf( "    %-12s %10s %12.3e %12.3e %12.3e\n", shape.name,
			             where, psiError, gradError, coarseError );

			double const worst = std::max( psiError, gradError );
			if ( inside )
				worstInside = std::max( worstInside, worst );
			else if ( ratio == 0.0 )
			{
				// SPLIT, BECAUSE psi AND ITS GRADIENT ARE NOT THE SAME
				// STATEMENT HERE and merging them loses the finding. On the
				// conductor's own surface psi_c is still at 1e-12 -- the
				// logarithm is integrable and the rule grades onto it -- while
				// grad psi_c is four to seven orders worse, because
				// differentiating under the integral puts the singularity one
				// power up and the grading is not steep enough for that.
				worstSurfacePsi = std::max( worstSurfacePsi, psiError );
				worstSurfaceGrad = std::max( worstSurfaceGrad, gradError );
			}
			else
			{
				worstClear = std::max( worstClear, worst );
				// REFINEMENT, WHICH IS THE PROPERTY M-139 SAYS THE BENCHMARK
				// LACKS. Only where there is room for it: past a few coil
				// sizes even order 8 is at round-off and the ratio of two
				// round-off figures is noise.
				if ( coarseError > 1.0e-10 )
					leastGain = std::min( leastGain, coarseError/std::max( psiError, 1.0e-16 ) );
			}
		}
	}

	std::printf( "    worst: clear of the coil %.3e; ON its surface, psi %.3e "
	             "and grad psi %.3e; at its centre %.3e\n",
	             worstClear, worstSurfacePsi, worstSurfaceGrad, worstInside );
	std::printf( "    least gain from order 8 to order %d, clear of the coil: "
	             "%.1fx\n", shipped, leastGain );

	// CLEAR OF THE CONDUCTOR THE SHIPPED RULE IS AT ROUND-OFF, and that is the
	// regime a subtracted solve lives in: Omega_h does not contain the
	// conductors -- that is the whole point of the split -- so every quadrature
	// point of every element is a point of this kind.
	BOOST_TEST( worstClear < 1.0e-8,
	            "the shipped conductor quadrature is only good to "
	            << worstClear << " relative CLEAR of the conductor, where the "
	            "integrand is analytic. A subtracted solve evaluates psi_c "
	            "nowhere else, so this bounds the split's quadrature error and "
	            "it has to be far below the discretisation" );

	// ON AND INSIDE IT, ALGEBRAIC AND LOOSER -- stated rather than hidden. The
	// only consumers are a MESHED run's own comparison and the guess builder,
	// neither of which is a claim about exactness.
	//
	// **AND THE GRADIENT IS THE WEAK ONE, BY SEVEN ORDERS.** The bound below
	// is set by the 265:1 solenoid at the middle of its long face: psi_c there
	// is 2.2e-12 and grad psi_c is 3.5e-05, on the same rule at the same
	// point. That is not a defect to fix at this order -- it is what
	// differentiating a logarithmic kernel under the integral does to a rule
	// graded for the kernel itself -- and it reaches nothing: a subtracted
	// solve has no quadrature point on a conductor's surface, because the
	// conductors are not in Omega_h. A MESHED run never evaluates psi_c at
	// all.
	BOOST_TEST( worstSurfacePsi < 1.0e-8,
	            "psi_c itself is only " << worstSurfacePsi << " relative ON "
	            "the conductor's surface, where the rule grades onto an "
	            "integrable logarithm and should still be near round-off" );
	BOOST_TEST( worstSurfaceGrad < 1.0e-4,
	            "grad psi_c is " << worstSurfaceGrad << " relative ON the "
	            "conductor's surface" );
	BOOST_TEST( worstInside < 1.0e-5,
	            "the shipped quadrature is " << worstInside << " relative at "
	            "the conductor's centre" );

	// AND IT REFINES, which is the whole of why this can be a reference and
	// freegs4e's six points per triangle cannot.
	BOOST_TEST( leastGain > 20.0,
	            "raising the order from 8 to " << shipped << " bought only "
	            << leastGain << "x at its weakest point clear of the "
	            "conductor. A reference that does not refine is M-139's "
	            "freegs4e floor, which is exactly what §7.2 says CS-5 cannot "
	            "be accepted against" );
}

// CS-5, THE SECOND HALF: THE SAME QUESTION ASKED OF A SOLVE RATHER THAN OF THE
// KERNEL.
//
// The case above bounds meq::coilPsi() against its own expensive form. That is
// necessary and it is not the claim: what CS-5 owes is a statement about the
// SPLIT, which is the kernel plus a boundary datum plus a solve. This puts the
// same rectangle through GradShafranovSolver at two quadrature orders -- the
// shipped one and one far beyond it -- with everything else identical, and
// measures how far apart the two reported equilibria are.
//
// THE COMPARISON IS AGAINST THE DISCRETISATION AND NOT AGAINST ZERO, which is
// the point §7.2 makes about M-139: a reference is only a reference if it
// resolves the thing being measured better than the thing being measured
// resolves itself. theMeshedCoilConvergesToTheSubtractedOne reads
// 1.6e-04 to 1.8e-05 on these meshes; a quadrature contribution of that size
// would make the split's acceptance a measurement of the quadrature rule.
BOOST_AUTO_TEST_CASE( theSolvedSplitIsInsensitiveToTheConductorQuadrature )
{
	meq::Coil const conductor( 1.00, 0.00, 0.10, 0.10, 1.0e5 );
	int const reference = 96;

	// THE DATUM IS THE EXPENSIVE FIELD IN BOTH ARMS, and that is what makes
	// this measure anything at all. Shift the datum with the quadrature and
	// the remainder is identically zero either way -- a beautiful, vacuous
	// agreement -- because the same error would be subtracted from both sides
	// of psi_p|_Gamma = g - psi_c|_Gamma. Holding g fixed at the reference
	// leaves the quadrature's own error as the only thing that moves.
	meq::CoilSet exact;
	exact.add( conductor );
	exact.setQuadratureOrder( reference );

	mfem::FunctionCoefficient datum(
		[ &exact ]( mfem::Vector const &x )
		{
			return exact.psi( x( 0 ), x( 1 ) );
		} );

	VacuumSource const vacuum;

	std::printf( "\n  CS-5: the SOLVED split at two conductor quadratures\n" );
	std::printf( "    %8s %16s %16s %10s\n", "elements", "max |psi diff|",
	             "max |psi|", "relative" );

	double worstRelative = 0.0;

	for ( int n : { 8, 16 } )
	{
		std::vector< double > total[ 2 ];
		double scale = 0.0;

		for ( int arm = 0; arm < 2; ++arm )
		{
			meq::ConductorField field;
			field.add( conductor );
			field.setQuadratureOrder( arm == 0 ? meq::defaultCoilQuadratureOrder
			                                   : reference );

			mfem::Mesh mesh = makeBox( n );
			meq::GradShafranovSolver solver( mesh, 2 );
			solver.setSource( vacuum );
			solver.setBoundaryData( datum );
			solver.setConductorField( field );
			solver.solve();

			// THE TOTAL AT EVERY QUADRATURE POINT, which is the field a
			// consumer sees -- psi_p alone would compare two remainders of two
			// different splits and mean nothing.
			for ( int element = 0; element < mesh.GetNE(); ++element )
			{
				mfem::IntegrationRule const &rule = mfem::IntRules.Get(
					mesh.GetElementBaseGeometry( element ), 4 );
				mfem::IsoparametricTransformation transformation;
				mesh.GetElementTransformation( element, &transformation );

				for ( int q = 0; q < rule.GetNPoints(); ++q )
				{
					mfem::IntegrationPoint const &ip = rule.IntPoint( q );
					transformation.SetIntPoint( &ip );
					mfem::Vector point( 3 );
					transformation.Transform( ip, point );

					double const value =
						solver.potential().GetValue( element, ip )
						+ field.psi( point( 0 ), point( 1 ) );
					total[ arm ].push_back( value );
					scale = std::max( scale, std::fabs( value ) );
				}
			}
		}

		BOOST_TEST_REQUIRE( total[ 0 ].size() == total[ 1 ].size() );

		double worst = 0.0;
		for ( std::size_t i = 0; i < total[ 0 ].size(); ++i )
			worst = std::max( worst, std::fabs( total[ 0 ][ i ]
			                                    - total[ 1 ][ i ] ) );

		double const relative = scale > 0.0 ? worst/scale : 0.0;
		worstRelative = std::max( worstRelative, relative );

		std::printf( "    %8d %16.3e %16.3e %10.2e\n",
		             static_cast< int >( total[ 0 ].size() ), worst, scale,
		             relative );
	}

	// FOUR ORDERS BELOW THE DISCRETISATION IS THE BAR, and it is set from the
	// sibling case's own table rather than picked: that reads 1.6e-04 at n = 8
	// and 1.8e-05 at n = 16 for the same conductor on the same meshes, so a
	// quadrature contribution at 1e-9 is a thousandth of the smaller of them
	// and the split's acceptance is a measurement of the DISCRETISATION.
	BOOST_TEST( worstRelative < 1.0e-9,
	            "the solved split moves by " << worstRelative << " relative "
	            "when the conductor quadrature is raised from the shipped "
	            "order to " << reference << ". COIL-SUBTRACTION-PLAN.md §7.2 "
	            "needs this far below the discretisation, or the split is "
	            "being accepted against its own quadrature rule" );
}

// §0b's CENTRAL CLAIM, ON THE ONE COMBINATION NO OTHER CASE COVERS:
// SPLIT + EXTERIOR COUPLING + PLASMA.
//
// THIS IS THE CELL NO OTHER ACCEPTANCE COVERS, WHICH IS §13.1's LESSON AND IS
// THE WHOLE REASON IT IS SEPARATE:
//
//   * theSplitReachesTheExteriorCoupling       -- coupling, NO plasma
//   * theMeshedCoilConvergesToTheSubtractedOne -- plasma, NO coupling
//   * DriverAcceptance's split case (M-144)    -- plasma, NO coupling
//
// Each is sharp and none has both, so until this case existed the combination
// every real machine is in had no acceptance at all -- and posing one found
// four defects, M-148.
//
// THE FIXTURE IS THE VACUUM CASE ABOVE WITH A PLASMA PUT IN IT, deliberately:
// same half-disc, same Gamma, same two rectangles inside it, same DtN. The only
// thing added is a source that carries current, so anything this case sees and
// its sibling does not is the plasma's doing.
//
// THE COMPARISON IS AGAINST THE MESHED ROUTE AND NOT AGAINST A CLOSED FORM,
// because there is no closed form for this and §0b's claim is a comparison
// anyway: the same conductors as a DOMAIN SOURCE, on the same mesh, must give
// the same equilibrium to the discretisation.
BOOST_AUTO_TEST_CASE( theSplitAndTheMeshedRouteAgreeUnderAnExteriorCoupling )
{
	int const order = 2;
	int const modes = 12;

	meq::ConductorField const subtracted = insideGammaPair();
	meq::CoilSet meshed;
	for ( std::size_t i = 0; i < subtracted.coils().size(); ++i )
		meshed.add( subtracted.coils().coil( i ) );

	// A PLASMA THAT IS A PLASMA: linear in psi, so the problem is genuinely
	// semi-linear and the plasma current is a functional of the answer rather
	// than a prescribed field. Kept small against the conductors' own scale so
	// that the equilibrium is one the conductors dominate -- which is the
	// regime the split exists for and the one where psi_c and psi_p are
	// comparable rather than one being noise on the other.
	struct LinearPlasma : public meq::Source
	{
		double f( double r, double, double psi ) const override
		{
			return 0.35*r*psi;
		}

		double dFdPsi( double r, double, double ) const override
		{
			return 0.35*r;
		}
	};

	// THE SAME SOURCE PLUS THE COILS, which is meq::CoilSet::f() added to the
	// plasma term -- the meshed route's whole content.
	struct MeshedPlasma : public meq::Source
	{
		explicit MeshedPlasma( meq::CoilSet const &c ) : coils( c ) {}

		double f( double r, double z, double psi ) const override
		{
			return 0.35*r*psi + coils.f( r, z );
		}

		double dFdPsi( double r, double, double ) const override
		{
			return 0.35*r;
		}

		meq::CoilSet const &coils;
	};

	LinearPlasma const plasma;
	MeshedPlasma const withCoils( meshed );

	std::printf( "\n  §0b UNDER AN EXTERIOR COUPLING, WITH A PLASMA\n" );
	std::printf( "    %5s %10s %16s %8s\n", "n", "elements",
	             "meshed - split", "rate" );

	double previous = 0.0;
	double first = 0.0;
	double last = 0.0;
	int firstN = 0;
	int lastN = 0;
	double lastRelative = 0.0;

	for ( int n : { 16, 24, 32 } )
	{
	std::vector< double > total[ 2 ];
	std::vector< std::pair< double, double > > at[ 2 ];
	double scale = 0.0;
	int elements = 0;

	for ( int arm = 0; arm < 2; ++arm )
	{
		HalfDisc d = makeHalfDisc( n );
		meq::ExteriorDtN const dtn( 0.0, gammaRadius, modes );
		mfem::ConstantCoefficient zero( 0.0 );

		meq::GradShafranovSolver solver( *d.sub, order );
		if ( arm == 0 )
			solver.setSource( withCoils );
		else
			solver.setSource( plasma );
		solver.setBoundaryData( zero );
		solver.setExtension( *d.path, d.gammaHMarker );
		if ( arm == 1 )
			solver.setConductorField( subtracted );
		solver.setExteriorCoupling( dtn );
		solver.solve();

		// THE PHYSICAL FLUX AT EVERY QUADRATURE POINT, which is the field a
		// consumer sees. Comparing psi_p against psi would compare a remainder
		// with a total and be large for the right reason.
		mfem::Mesh &mesh = *d.sub;
		for ( int element = 0; element < mesh.GetNE(); ++element )
		{
			mfem::IntegrationRule const &rule = mfem::IntRules.Get(
				mesh.GetElementBaseGeometry( element ), 4 );
			mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( element, &transformation );

			for ( int q = 0; q < rule.GetNPoints(); ++q )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( q );
				transformation.SetIntPoint( &ip );
				mfem::Vector point( 3 );
				transformation.Transform( ip, point );

				// NOT ON THE CONDUCTORS, and this exclusion is its sibling's
				// and carries its reason: inside and beside a rectangle the
				// MESHED route is resolving a top hat cut by element
				// interiors, so psi is only C^1 there and no polynomial
				// degree recovers it. theMeshedCoilConvergesToTheSubtractedOne
				// excludes the same band for the same reason and measures a
				// rate 2.69 outside it.
				//
				// **AND LEAVING IT IN MEASURES THE WRONG ARM, WHICH WAS
				// MEASURED BEFORE IT WAS EXCLUDED**: the worst point of the
				// whole comparison landed at ( 0.5602, 0.1995 ), one
				// centimetre outside the coil at ( 0.50, 0.20 ), with the
				// meshed arm reading 5.12e-02 against the subtracted arm's
				// 1.83e-01. The SUBTRACTED number is the accurate one there --
				// it is an exact quadrature of the same rectangle -- so a
				// comparison including that point reports the split's
				// ADVANTAGE as its error.
				bool nearConductor = false;
				for ( std::size_t c = 0; c < subtracted.coils().size(); ++c )
				{
					meq::Coil const &one = subtracted.coils().coil( c );
					if ( std::fabs( point( 0 ) - one.centreR() )
					         < 2.0*one.halfWidth()
					     && std::fabs( point( 1 ) - one.centreZ() )
					         < 2.0*one.halfHeight() )
						nearConductor = true;
				}
				if ( nearConductor )
					continue;

				double value = solver.potential().GetValue( element, ip );
				if ( arm == 1 )
					value += subtracted.psi( std::max( 0.0, point( 0 ) ),
					                         point( 1 ) );
				total[ arm ].push_back( value );
				at[ arm ].push_back( std::make_pair( point( 0 ),
				                                     point( 1 ) ) );
				scale = std::max( scale, std::fabs( value ) );
			}
		}

		elements = mesh.GetNE();
	}

	BOOST_TEST_REQUIRE( total[ 0 ].size() == total[ 1 ].size() );

	double worst = 0.0;
	std::size_t worstAt = 0;
	for ( std::size_t i = 0; i < total[ 0 ].size(); ++i )
	{
		double const gap = std::fabs( total[ 0 ][ i ] - total[ 1 ][ i ] );
		if ( gap > worst )
		{
			worst = gap;
			worstAt = i;
		}
	}
	double const relative = scale > 0.0 ? worst/scale : 0.0;

	// WHERE, AND NOT ONLY HOW MUCH. A disagreement concentrated at Gamma is a
	// boundary term; one concentrated on the conductors is the source; one
	// spread evenly is neither. That is the first question anyone reading a red
	// here will have -- and it is what identified the exclusion above.
	double const rate = previous > 0.0 ? std::log2( previous/worst ) : 0.0;
	if ( previous == 0.0 )
	{
		first = worst;
		firstN = n;
	}
	last = worst;
	lastN = n;
	previous = worst;
	lastRelative = relative;

	std::printf( "    %5d %10d %16.3e", n, elements, worst );
	if ( rate != 0.0 )
		std::printf( " %8.2f", rate );
	std::printf( "        worst at ( %.4f, %.4f ), |x| = %.4f\n",
	             at[ 0 ][ worstAt ].first, at[ 0 ][ worstAt ].second,
	             std::hypot( at[ 0 ][ worstAt ].first,
	                         at[ 0 ][ worstAt ].second ) );
	}

	// THE GATE IS THE MESHED ARM'S OWN DISCRETISATION, as M-144's is and for
	// the same reason: the coils are cut by element interiors, so the meshed
	// source carries a jump inside an element and its error is per cent level.
	// A split that agreed to round-off would be suspicious; one that differs by
	// more than the meshed arm can be wrong by is solving another problem.
	BOOST_TEST( lastRelative < 1.0e-1,
	            "the subtracted and meshed routes disagree by " << lastRelative
	            << " relative under an exterior coupling with a plasma. "
	            "COIL-SUBTRACTION-PLAN.md §0b requires the split to be a change "
	            "of REPRESENTATION and not of answer" );

	/*
	 * AND THE DIFFERENCE FALLS, WHICH IS THE HALF WITH TEETH -- BUT ACROSS THE
	 * WHOLE SEQUENCE AND NOT PER PAIR, AND THAT IS FORCED.
	 *
	 * A tolerance says the two routes are close on one mesh; a rate says they
	 * are converging to the same field, which is §0b's actual claim.
	 *
	 * **PER PAIR IS NOT AVAILABLE HERE AND ITS SIBLING SHOWS WHY.**
	 * theMeshedCoilConvergesToTheSubtractedOne can assert per pair only because
	 * its coil's edges sit on mesh VERTICES at every level -- it says so, and
	 * records 6.756e-04, 1.276e-03, 1.917e-04 with unaligned half-extents, a
	 * rate of -0.92 and then +2.73. This fixture is a half-disc with a curved
	 * Gamma and cannot align anything: the meshed source is a top hat whose
	 * edge cuts element interiors, and WHICH elements it cuts is not a smooth
	 * function of h. Measured here: 1.260e-01, 8.506e-03, 9.206e-03, so +3.89
	 * and then -0.11, with the worst point hopping between coils.
	 *
	 * That is the same two-tier rule CLAUDE.md states for the extension path
	 * under *Unfitted convergence needs a two-tier rate assertion*, met here
	 * for the same reason, and it is geometry rather than either conductor
	 * route.
	 */
	double const endToEnd =
		std::log2( first/last )
		/std::log2( static_cast< double >( lastN )/firstN );
	std::printf( "    across the sequence: %.2f\n", endToEnd );

	BOOST_TEST( endToEnd > 0.9,
	            "the two routes are not converging to each other: "
	            << endToEnd << " across the sequence under an exterior "
	            "coupling with a plasma" );
}
