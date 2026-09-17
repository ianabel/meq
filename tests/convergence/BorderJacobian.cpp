/*
 * THE BORDERED JACOBIAN'S ROWS, AGAINST A DIFFERENCE OF THE CONSTRAINTS.
 *
 * tests/convergence/NewtonConvergence.cpp checks the assembled FIELD Jacobian
 * against a central difference of the assembled residual, and CLAUDE.md names
 * that as the check a convergence table cannot replace: *A wrong Jacobian is
 * invisible to a convergence table*, because Newton reaches the same discrete
 * solution whatever carried it there. **The BORDER has never had the
 * equivalent**, and MEASUREMENTS.md M-118 records two defects that one would
 * have caught -- a fabricated current column and a deleted corner -- plus a
 * third that is still open and is what this case is built to settle.
 *
 *
 * THE OPEN QUESTION, AND WHY IT HAS A SHARP FORM.
 *
 * The axis constraint is
 *
 *     C_Ax( u ) = s - psi_h( x* ),    x* defined by q_h( x* ) = 0,
 *
 * and differentiating it gives two terms,
 *
 *     dC_Ax/du = -[ dpsi_h/du |_x*  +  grad psi_h( x* ) . d( x* )/du ].
 *
 * GradShafranov.cpp assembles the first and drops the second by the envelope
 * theorem. That step was justified by "grad_bar( psi ) = r q, so
 * grad( psi_h )( x* ) = 0 at a zero of q_h IDENTICALLY" -- and the relation
 * quoted is the CONTINUOUS one. The discrete flux equation gives, per element
 * and for all v in the flux space,
 *
 *     ( r q_h - grad_bar psi_h, v )_K = -< psi_h - psihat_h, v.n >_dK,
 *
 * so r q_h - grad_bar psi_h is the local lifting of the trace jump rather than
 * zero. cornerEntry()'s XP-3 arm says exactly this about the same identity --
 * "psi_h and q_h are separate solved fields whose identity is only weak" -- and
 * the two cannot both be right.
 *
 *
 * THE SHARP FORM IS THAT THE NEGLECTED TERM LIVES ENTIRELY ON THE FLUX DOFS,
 * WHERE MEQ'S ROW IS EXACTLY ZERO.
 *
 * x* is a root of q_h, so it moves when the FLUX moves and not when the
 * potential does. MEQ's row is supported on the axis element's POTENTIAL dofs,
 * so what it asserts about the flux block is
 *
 *     dC_Ax/d( flux dofs ) = 0, exactly.
 *
 * The truth is -grad psi_h( x* ) . d( x* )/d( flux ), and it is zero if and only if
 * grad psi_h( x* ) is. **So the whole question is decidable by a perturbation
 * that moves q_h and leaves psi_h alone**, which is what the second case below
 * does: it is not an approximation of the row, it is the row's own claim about
 * a block it does not touch.
 *
 * That framing is also why this case needs no new solver API. Everything it
 * reads -- flux(), potential(), and meq::CriticalPointFinder over bare fields --
 * is public, and the quantity in question is a property of the SOLVED FIELDS
 * rather than of the assembly code.
 *
 *
 * WHAT THIS CASE DOES NOT COVER, said plainly so the gap is not mistaken for
 * coverage. The border COLUMNS -- dR/dc -- are checked by
 * HighBetaConvergence.cpp's theAnalyticColumnAgreesWithTheDifferencedOne and
 * theAnalyticCurrentCornerAgreesWithTheDifferencedOne, which run the two routes
 * against each other through a whole solve. XP-3's four grad-q entries and the
 * psi_bnd row's two are exact arithmetic on an element polynomial and are
 * pinned by XPointBorder's observed order. What was uncovered, and is covered
 * here, is the ENVELOPE ASSUMPTION the axis row rests on.
 */
#define BOOST_TEST_MODULE BorderJacobian
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"
#include "meq/CriticalPoints.hpp"
#include "convergence/ConvergenceHarness.hpp"

namespace
{
	using meq::tests::standardBox;

	mfem::FunctionCoefficient bump( double height )
	{
		meq::tests::Rectangle const box = standardBox();
		double const rMin = box.rMin;
		double const zMin = box.zMin;
		double const width = box.width();
		double const depth = box.height();
		return mfem::FunctionCoefficient(
			[ height, rMin, zMin, width, depth ]( mfem::Vector const &x )
			{
				return height*std::sin( M_PI*( x( 0 ) - rMin )/width )
				       *std::sin( M_PI*( x( 1 ) - zMin )/depth );
			} );
	}

	/// One converged bordered solve, and everything read off it afterwards.
	/// Held together because the solver owns the spaces the fields live in, so
	/// a Machine that returned only the fields would hand out danglers.
	struct Machine
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::shared_ptr<meq::ConstantProfile const> pPrime;
		std::shared_ptr<meq::ConstantProfile const> ggPrime;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<mfem::ConstantCoefficient> zero;
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		double h = 0.0;
		bool converged = false;
	};

	std::unique_ptr<Machine> solveAt( int n, int order )
	{
		auto m = std::make_unique<Machine>();
		meq::tests::Rectangle const box = standardBox();
		m->mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( box, n ) );
		m->h = box.width()/n;

		m->pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		m->ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
		m->source = std::make_unique<meq::NormalisedMHDSource>( m->pPrime,
		                                                        m->ggPrime,
		                                                        1.0, 1.0 );
		m->zero = std::make_unique<mfem::ConstantCoefficient>( 0.0 );
		m->guess = std::make_unique<mfem::FunctionCoefficient>( bump( 0.30 ) );

		double const limiterR = box.rMin + 0.68*( box.rMax - box.rMin );
		double const limiterZ = box.zMin + 0.31*( box.zMax - box.zMin );

		m->solver = std::make_unique<meq::GradShafranovSolver>( *m->mesh, order );
		m->solver->setBoundaryFluxPoint( limiterR, limiterZ );
		m->solver->setSource( *m->source, 0.30 );
		m->solver->setBoundaryData( *m->zero );
		m->solver->setInitialGuess( *m->guess );
		m->solver->setNewtonControl( 1.0e-12, 1.0e-14, 40 );
		try
		{
			m->solver->solve();
			m->converged = true;
		}
		catch ( std::exception const &error )
		{
			std::printf( "    n = %d did not converge: %s\n", n, error.what() );
		}
		return m;
	}

	/// The value and gradient of a scalar field at a physical point, through the
	/// element that contains it. Returns false when the point is not in the
	/// mesh, which a located axis never is.
	bool valueAndGradientAt( mfem::GridFunction const &field, double r, double z,
	                         double &value, mfem::Vector &gradient )
	{
		mfem::Mesh *mesh = field.FESpace()->GetMesh();
		mfem::DenseMatrix point( 2, 1 );
		point( 0, 0 ) = r;
		point( 1, 0 ) = z;
		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> reference;
		if ( mesh->FindPoints( point, elements, reference ) < 1
		     || elements[ 0 ] < 0 )
			return false;

		// The caller-allocated overload: Mesh::GetElementTransformation( int )
		// hands out shared scratch, which CLAUDE.md records at six other sites.
		thread_local mfem::IsoparametricTransformation scratch;
		mesh->GetElementTransformation( elements[ 0 ], &scratch );
		scratch.SetIntPoint( &reference[ 0 ] );
		value = field.GetValue( scratch, reference[ 0 ] );
		field.GetGradient( scratch, gradient );
		return true;
	}
}

/*
 * THE IDENTITY IS WEAK, AND THE TWO SIDES CONVERGE AT DIFFERENT ORDERS.
 *
 * At the located axis q_h is zero by construction -- it is the root the finder
 * returned -- so the comment's claim reduces to grad psi_h( x* ) being zero
 * too. It is not: it is the local lifting of the trace jump, and it converges
 * one order slower than q_h does. That gap IS the mixed method, which is the
 * reason this is worth measuring rather than asserting either way.
 */
BOOST_AUTO_TEST_CASE( theFluxAndThePotentialGradientAgreeOnlyWeaklyAtTheAxis )
{
	int const order = 2;

	std::printf( "\n  grad psi_h AT A ZERO OF q_h, WHICH THE AXIS ROW CALLS ZERO\n" );
	std::printf( "    %5s %10s %14s %14s %14s %10s\n", "n", "h", "|q_h(x*)|",
	             "|grad psi_h|", "|grad psi*|", "star/h" );

	std::vector<double> sizes;
	std::vector<double> gradients;

	for ( int n : { 8, 12, 16, 24 } )
	{
		std::unique_ptr<Machine> m = solveAt( n, order );
		BOOST_TEST_REQUIRE( m->converged, "n = " << n << " must solve" );

		meq::CriticalPointFinder finder( *m->solver );
		meq::CriticalPoint axis;
		BOOST_TEST_REQUIRE( finder.tryFindAxis( axis ),
		                    "no axis located at n = " << n );

		double value = 0.0;
		mfem::Vector gradient;
		BOOST_TEST_REQUIRE( valueAndGradientAt( m->solver->potential(), axis.r,
		                                        axis.z, value, gradient ),
		                    "the located axis is not in the mesh at n = " << n );

		double const gradientNorm = gradient.Norml2();

		/*
		 * AND THE SAME AT THE POST-PROCESSED POTENTIAL, because psi* is the
		 * field whose gradient IS the solved flux -- that is what the local
		 * post-processing solves for. If the envelope argument is to hold
		 * anywhere it is here, and this column is what says by how much.
		 *
		 * It is NOT exact even so: psi* is the L2 projection of r q_h onto the
		 * gradients of P^(k+2) on the element, so grad psi* is the closest
		 * gradient field to r q_h rather than r q_h itself. What it should buy
		 * is an ORDER -- grad psi* converges at O( h^(k+1) ) where grad psi_h
		 * manages O( h^k ) -- and an order is the difference between a term
		 * worth assembling and one that is not.
		 */
		m->solver->postProcess();
		double starValue = 0.0;
		mfem::Vector starGradient;
		BOOST_TEST_REQUIRE( valueAndGradientAt( m->solver->postProcessedPotential(),
		                                        axis.r, axis.z, starValue,
		                                        starGradient ),
		                    "psi* cannot be evaluated at the located axis at n = "
		                    << n );
		double const starNorm = starGradient.Norml2();
		double rate = 0.0;
		if ( !sizes.empty() )
			rate = std::log( gradients.back()/gradientNorm )
			       /std::log( sizes.back()/m->h );
		std::printf( "    %5d %10.4f %14.4e %14.4e %14.4e %10.4f\n", n, m->h,
		             axis.fluxResidual, gradientNorm, starNorm,
		             starNorm/std::max( gradientNorm, 1.0e-300 ) );
		std::fflush( stdout );

		sizes.push_back( m->h );
		gradients.push_back( gradientNorm );

		// q_h IS ZERO THERE BY CONSTRUCTION and this is the control: if the
		// root finder had not converged, everything below would be measuring
		// its residual rather than the identity.
		BOOST_TEST( axis.fluxResidual < 1.0e-9,
		            "the located axis is not a root of q_h at n = " << n
		            << ": |q_h| = " << axis.fluxResidual );
	}

	/*
	 * THE ASSERTION IS THAT IT IS NOT ZERO, which is the opposite direction
	 * from most of this suite and is deliberate: the claim under test is an
	 * IDENTITY, so any single measurement above round-off falsifies it, and a
	 * tolerance would only decide how loudly.
	 */
	for ( std::size_t i = 0; i < gradients.size(); ++i )
		BOOST_TEST( gradients[ i ] > 1.0e-10,
		            "grad psi_h at the located axis reads " << gradients[ i ]
		            << ", which is round-off -- the envelope argument's "
		               "identity would then hold and the axis row needs no "
		               "correction. Re-read GradShafranov.cpp's axis row." );
}

/*
 * AND THE CONSEQUENCE, COMPUTED RATHER THAN DIFFERENCED: THE AXIS ROW NEGLECTS
 * A TERM ITS OWN MACHINERY COULD ASSEMBLE.
 *
 * THE FIRST VERSION OF THIS CASE DIFFERENCED IT AND THE INSTRUMENT WAS WRONG.
 * Perturbing the flux, relocating x* and differencing psi_h( x* ) gave
 * -5.760e-04, -1.158e-03, -2.321e-03 as the step was halved -- doubling each
 * time, which means the NUMERATOR was constant at 1.09e-05 whatever the step.
 * A constant numerator is the signature of a discontinuity and not of a
 * derivative: the relocated axis crosses a face, psi_h is L2 and jumps there by
 * O( h^(k+1) ), and the difference reports the jump however finely it is
 * probed. See the case below, which keeps that observation because it is a
 * finding in its own right.
 *
 * THE TERM IS AVAILABLE IN CLOSED FORM, SO NOTHING HAS TO BE DIFFERENCED.
 * MEQ's flux space is a scalar collection at vdim = 2 (`GradShafranov.cpp`'s
 * `fluxFes`), so with `q_h = sum_j c_(j,d) phi_j e_d`,
 *
 *     dq_d/dc_(j,d') = phi_j delta_(d d')
 *
 * and from q_h( x* ) = 0,
 *
 *     d( x* )/dc_(j,d') = -( grad q_h )^-1 e_d' phi_j( x* ).
 *
 * The neglected piece of the row is -grad psi_h . d( x* )/dc, so writing
 * w := ( grad q_h )^-T grad psi_h its entry on that dof is w_d' phi_j( x* ) and
 * its norm over the element's flux dofs is | w | | phi( x* ) |. The row MEQ
 * DOES assemble has norm | phi( x* ) | over the potential dofs, so when the two
 * spaces carry the same degree -- they do -- the ratio of neglected to kept is
 * simply
 *
 *     | ( grad q_h( x* ) )^-T grad psi_h( x* ) |.
 *
 * Every factor is already assembled elsewhere in the solver for XP-3's own
 * rows: `xFluxJacobian` is grad q_h and `xFluxShape` is dq_h/du. Writing the
 * correction is one 2x2 solve and a scatter.
 *
 * THIS CASE ASSERTS THE BEHAVIOUR THAT IS WANTED AND IS EXPECTED TO BE RED
 * UNTIL THE TERM IS WRITTEN, per CLAUDE.md's testing stance -- a test that
 * asserted the defect would pass while the defect stood, and that is exactly
 * what makes a green suite compatible with a broken solver. The bar is
 * NEGLIGIBLE-FOR-A-JACOBIAN rather than round-off: an exactly zero neglected
 * term is not achievable, and Newton does not need one.
 */
BOOST_AUTO_TEST_CASE( theAxisRowNeglectsATermItsOwnMachineryCouldAssemble )
{
	int const order = 2;

	std::printf( "\n  THE AXIS ROW'S NEGLECTED ENVELOPE TERM, IN CLOSED FORM\n" );
	std::printf( "    %5s %10s %14s %14s %14s %12s\n", "n", "h",
	             "|grad psi_h|", "|grad q_h|", "|neglected|", "of the row" );

	for ( int n : { 8, 12, 16, 24 } )
	{
		std::unique_ptr<Machine> m = solveAt( n, order );
		BOOST_TEST_REQUIRE( m->converged, "n = " << n << " must solve" );

		meq::CriticalPointFinder finder( *m->solver );
		meq::CriticalPoint axis;
		BOOST_TEST_REQUIRE( finder.tryFindAxis( axis ), "no axis at n = " << n );

		mfem::GridFunction const &potential = m->solver->potential();
		mfem::GridFunction const &flux = m->solver->flux();

		// THE AXIS ELEMENT, located once and used for every evaluation below, so
		// that grad psi_h, grad q_h and both shape vectors come from ONE
		// polynomial. Reading them from different elements is the trap this
		// whole case is about.
		mfem::Mesh *mesh = potential.FESpace()->GetMesh();
		mfem::DenseMatrix point( 2, 1 );
		point( 0, 0 ) = axis.r;
		point( 1, 0 ) = axis.z;
		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> local;
		BOOST_TEST_REQUIRE( mesh->FindPoints( point, elements, local ) == 1,
		                    "the located axis is not in the mesh at n = " << n );

		thread_local mfem::IsoparametricTransformation scratch;
		mesh->GetElementTransformation( elements[ 0 ], &scratch );
		scratch.SetIntPoint( &local[ 0 ] );

		mfem::Vector gradPsi;
		potential.GetGradient( scratch, gradPsi );

		mfem::DenseMatrix gradQ;
		flux.GetVectorGradient( scratch, gradQ );

		// w = ( grad q_h )^-T grad psi_h, solved rather than inverted.
		mfem::DenseMatrix transposed( gradQ, 't' );
		mfem::Vector w( 2 );
		{
			mfem::DenseMatrixInverse inverse( transposed );
			inverse.Mult( gradPsi, w );
		}

		// The two shape vectors, so the ratio is reported honestly rather than
		// assumed equal. They are the same degree today and the assertion does
		// not depend on that staying true.
		mfem::FiniteElement const &potentialEl =
			*potential.FESpace()->GetFE( elements[ 0 ] );
		mfem::FiniteElement const &fluxEl =
			*flux.FESpace()->GetFE( elements[ 0 ] );
		mfem::Vector potentialShape( potentialEl.GetDof() );
		mfem::Vector fluxShape( fluxEl.GetDof() );
		potentialEl.CalcShape( local[ 0 ], potentialShape );
		fluxEl.CalcShape( local[ 0 ], fluxShape );

		double const kept = potentialShape.Norml2();
		double const neglected = w.Norml2()*fluxShape.Norml2();
		double const ratio = neglected/kept;

		std::printf( "    %5d %10.4f %14.4e %14.4e %14.4e %12.4e\n", n, m->h,
		             gradPsi.Norml2(), gradQ.FNorm(), neglected, ratio );
		std::fflush( stdout );

		/*
		 * THE BAR. A Jacobian entry a millionth of the row it sits beside
		 * cannot move a Newton step; anything larger can, and the observed
		 * order is the only thing that would ever show it -- which is exactly
		 * the defect CLAUDE.md's *A wrong Jacobian is invisible to a
		 * convergence table* says no error norm can see.
		 */
		bool const acceptable =
			meq::GradShafranovSolver::axisRowCarriesEnvelopeTerm
			|| ratio <= 1.0e-6;
		BOOST_TEST( acceptable,
		            "at n = " << n << " the axis row NEGLECTS a term of norm "
		            << neglected << " beside a kept row of norm " << kept
		            << " ( ratio " << ratio << " ). It is "
		               "grad psi_h( x* )^T ( grad q_h )^-1 ( flux shape at x* ) "
		               "on the axis element's FLUX dofs, where the row is "
		               "currently zero. The envelope theorem drops it only if "
		               "grad psi_h( x* ) = 0, which the case above measures as "
		               "false. xFluxJacobian and xFluxShape already hold every "
		               "factor. THE FIX is to assemble it in "
		               "solveWithNormalisation()'s constraintLocated branch and "
		               "flip GradShafranovSolver::axisRowCarriesEnvelopeTerm, "
		               "which is the constant this assertion reads." );
	}
}

/*
 * AND THE AXIS CONSTRAINT IS DISCONTINUOUS IN THE FLUX, WHICH IS A SECOND
 * FINDING AND NOT THE SAME ONE.
 *
 * MEASUREMENTS.md M-115's surviving hypothesis for the line-search plateau is
 * that `psi_bnd - psi_h( x_X )` is a point evaluation of an L2 field and jumps
 * when the X-point crosses a face. M-117 records that the AXIS row carries 0.75
 * of MAST's merit and that there is no instrument for it at all -- BorderStep
 * carries trialXR, trialXZ, trialXElement and trialBoundaryConstraint for the
 * X-point and nothing for the axis.
 *
 * This is that instrument, in its cheapest form. Perturbing the flux moves x*,
 * and when the move crosses a face the value of psi_h there jumps by the L2
 * discontinuity -- so the difference quotient reports a CONSTANT numerator
 * however finely the step is refined, which is what the first version of the
 * case above ran into. The element index is printed because it is the evidence:
 * a constant numerator with an unchanged element would be something else.
 *
 * NOTHING IS ASSERTED ABOUT THE SIZE. What is asserted is the shape -- that the
 * numerator does not fall with the step -- because that is what distinguishes a
 * discontinuity from a derivative and is the claim M-115 needs. A tolerance on
 * the jump would be a tolerance on h and on where the axis happens to sit.
 */
BOOST_AUTO_TEST_CASE( theAxisConstraintJumpsWhenTheAxisCrossesAFace )
{
	int const order = 2;
	int const n = 16;

	std::unique_ptr<Machine> m = solveAt( n, order );
	BOOST_TEST_REQUIRE( m->converged, "the reference solve must converge" );

	mfem::GridFunction const &flux = m->solver->flux();
	mfem::GridFunction const &potential = m->solver->potential();

	// A FIXED SEED, so a failure is reproducible; random, because a hand-chosen
	// direction can be accidentally orthogonal to what is under test.
	std::mt19937 generator( 20260917u );
	std::normal_distribution<double> normal( 0.0, 1.0 );
	mfem::Vector direction( flux.Size() );
	for ( int i = 0; i < direction.Size(); ++i )
		direction( i ) = normal( generator );
	direction /= direction.Norml2();

	double const scale = std::max( flux.Norml2(), 1.0e-30 );

	std::printf( "\n  THE AXIS CONSTRAINT UNDER A FLUX PERTURBATION\n" );
	std::printf( "    %12s %16s %16s %10s %10s\n", "step", "numerator",
	             "|x* moved|", "elem -", "elem +" );

	auto axisUnder = [ & ]( double signedStep, int &element, double &r,
	                        double &z )
	{
		mfem::GridFunction perturbed( flux );
		perturbed.Add( signedStep, direction );

		// The bare-field finder, which exists for exactly this: the root is
		// taken on the PERTURBED q and the value on the UNPERTURBED psi, so
		// what moves is the point alone.
		meq::CriticalPointFinder moved( perturbed, potential );
		meq::CriticalPoint found;
		if ( !moved.tryFindAxis( found ) )
			return std::numeric_limits<double>::quiet_NaN();
		element = found.element;
		r = found.r;
		z = found.z;

		double value = 0.0;
		mfem::Vector gradient;
		if ( !valueAndGradientAt( potential, found.r, found.z, value, gradient ) )
			return std::numeric_limits<double>::quiet_NaN();
		return value;
	};

	std::vector<double> numerators;
	for ( double const relative : { 1.0e-4, 5.0e-5, 2.5e-5 } )
	{
		double const step = relative*scale;

		int elementPlus = -1, elementMinus = -1;
		double rPlus = 0.0, zPlus = 0.0, rMinus = 0.0, zMinus = 0.0;
		double const plus = axisUnder( step, elementPlus, rPlus, zPlus );
		double const minus = axisUnder( -step, elementMinus, rMinus, zMinus );
		bool const bothFinite = std::isfinite( plus ) && std::isfinite( minus );
		BOOST_TEST_REQUIRE( bothFinite,
		                    "the axis could not be relocated at step " << step );

		numerators.push_back( plus - minus );
		std::printf( "    %12.4e %16.6e %16.4e %10d %10d\n", step, plus - minus,
		             std::hypot( rPlus - rMinus, zPlus - zMinus ), elementMinus,
		             elementPlus );
		std::fflush( stdout );
	}

	/*
	 * A DERIVATIVE'S NUMERATOR IS PROPORTIONAL TO THE STEP AND WOULD FALL BY
	 * FOUR OVER THESE THREE. A jump's does not move at all. The bar is set at
	 * "fell by less than half" so that the finding is the SHAPE and not a
	 * particular constant.
	 */
	double const first = std::abs( numerators.front() );
	double const last = std::abs( numerators.back() );
	BOOST_TEST( last >= 0.5*first,
	            "the numerator fell from " << first << " to " << last
	            << " as the step was quartered, so the axis constraint is "
	               "responding smoothly here after all and this case's premise "
	               "-- M-115's face-crossing hypothesis, on the axis row -- "
	               "needs re-reading" );
}
