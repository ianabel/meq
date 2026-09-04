#define BOOST_TEST_MODULE SurfaceAverageConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/SurfaceAverage.hpp"

#include "analytic/FluxSurfaceReference.hpp"
#include "analytic/Soloviev.hpp"
#include "ConvergenceHarness.hpp"

/*
 * The acceptance test for INVERSION-PLAN.md stage IN-2: flux-surface averages,
 * ON THE FITTED PATH, where Gamma_h is Gamma and there is no band -- so every
 * rate here is attributable to the discretisation and the quadrature and to
 * nothing else.
 *
 * FOUR ACCEPTANCES, AND THEY ARE INDEPENDENT OF EACH OTHER ON PURPOSE.
 *
 * 1. AGAINST A CONVERGED REFERENCE ON THE EXACT FIELD.
 *    tests/analytic/FluxSurfaceReference.hpp computes the same averages from the
 *    analytic psi and grad psi by rays plus the periodic trapezoid, and reaches
 *    round-off. IT NEVER TOUCHES psi_h, so a comparison against it is
 *    INVERSION-PLAN.md section 2's error ( a ) -- the discretisation -- with
 *    ( b ) and ( c ) removed by construction. It is a reference VALUE and not a
 *    closed form, and that distinction is stated wherever the number is printed:
 *    there IS no closed form, because psi is elementary and an integral over a
 *    contour of it is not.
 *
 *    Two axes are measured separately, because they are separate errors:
 *    SPECTRAL CONVERGENCE IN THE NUMBER OF ANGULAR POINTS at fixed mesh, and
 *    convergence IN h at fixed angular resolution.
 *
 * 2. THE AVERAGED GRAD-SHAFRANOV IDENTITY, WHICH NEEDS NO REFERENCE VALUE AT
 *    ALL. Three averages check each other with nothing but the equation. The
 *    right-hand side is written with the F THE SOLVER IS FED and never as a
 *    hand-derived -mu_0 p' - g g' < R^-2 >, for the reason
 *    SolovievEquilibrium::deltaStarFD() exists: an independent quantity is the
 *    only thing that catches a misread formula, and a right-hand side
 *    re-derived by hand is not independent of the hand that derived it.
 *
 * 3. TWO INDEPENDENT EXTRACTIONS MUST AGREE. The predictor-corrector trace with
 *    its cubic Hermite, and the ray-based angle parametrisation, are genuinely
 *    different routes to the same surface -- different nodes, different metric,
 *    different rule. Agreement is worth more than either being plausible. The
 *    THIRD leg, section 3.3's implicit quadrature, is explicitly deferred and
 *    the header of src/meq/SurfaceAverage.hpp says so; read this as two of
 *    three.
 *
 * 4. THE BAND FLAG IS CARRIED. Zero here, and asserted to be zero rather than
 *    assumed: the fitted path has no band, so a false is the correct answer
 *    rather than an absent one, and an edit that starts reporting one on this
 *    path should fail loudly.
 *
 * THE RESULT THAT WAS NOT EXPECTED, AND IT IS IN THE FIRST TABLE: psi* DOES NOT
 * BUY AN ORDER HERE. The tracer's default pairing roots psi*, which converges at
 * k+2, so the natural expectation is k+2 in the averages. Every quantity in this
 * file converges at k+1 with EITHER pairing, because the weight divides by
 * | grad psi | = r | q | and q* converges at k+1 like q_h -- the reconstruction
 * buys its order in the potential and there is no k+2 flux to be had. What psi*
 * buys is a constant, measured at a factor of 1.3 to 7.5 on the finest mesh. The
 * table asserts k+1 on both columns and prints the ratio at every degree.
 *
 * A RATE IS TAKEN OVER THE ROOT MEAN SQUARE OF FOUR SURFACES, NOT OVER ONE, and
 * that is a measurement rather than a convenience. The error of an average is a
 * SIGNED quantity whose magnitude is not monotone in h: on a single surface the
 * same sequence gives pairwise rates of 0.79 and then 5.93, because the
 * integrand's error changes sign around the contour and the cancellation is not
 * the same at every mesh. Four surfaces are smooth enough to carry a rate, and
 * a profile is what the consumer reads anyway.
 *
 * WHY nstx() AND WHY A BOX THAT IS NOT standardBox(). The reference values are
 * recorded on SolovievEquilibrium::nstx(), whose axis sits at r = 1.318, and its
 * surfaces are ELONGATED: the Psi_N = 0.50 surface reaches r in [ 0.81, 1.66 ]
 * and z in [ -0.87, 0.90 ], which does not fit inside [0.6,1.4]x[-0.6,0.6] at
 * all. The box below is chosen so that every surface measured -- Psi_N = 0.15,
 * 0.25, 0.35 and 0.50 -- sits well inside it: the smallest psi anywhere on its
 * boundary is Psi_N = 0.679, so all four levels are strictly interior closed
 * curves.
 *
 * THE THIRD TABULATED SURFACE, Psi_N = 0.75, IS NOT MEASURED HERE AND THAT IS A
 * PROPERTY OF THE FIXTURE RATHER THAN OF THE METHOD. It reaches z = +/- 1.18 and
 * r = 1.73, so enclosing it needs a box whose margin at the coarsest mesh of any
 * dyadic sweep is under one cell -- at which point what would be measured is the
 * contour's distance to the mesh boundary. Its reference value is still asserted
 * below, because the reference is a function of the level and costs no solve.
 */

namespace
{

	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	double const twoPi = 6.283185307179586476925286766559;

	/// How far below k+1 a rate taken ACROSS THE WHOLE SEQUENCE may fall. Pairs
	/// are printed and never asserted on, which is the two-tier pattern of
	/// ExtensionConvergence.cpp and is here for the same reason: which elements
	/// a contour passes through is not a smooth function of h, so the constant
	/// in front of the power is not the same constant at every level. A genuine
	/// loss of an order reads about 1 below the target and fails this
	/// comfortably.
	double const sequenceSlack = 0.40;

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}

	/// The box. See the header for why it is not standardBox().
	Rectangle nstxBox()
	{
		return Rectangle{ 0.60, 1.90, -1.10, 1.10 };
	}

	/// The magnetic axis of the CLOSED FORM, to round-off. Newton on
	/// grad( psi ) = 0 with the Hessian by central differences; the Hessian's
	/// accuracy does not reach the answer, the fixed point being where the
	/// analytic gradient vanishes whatever steered it there. Lifted from
	/// FluxSurfaceConvergence.cpp, which is where the argument for it is.
	struct ExactAxis
	{
		double r;
		double z;
		double psi;
	};

	ExactAxis exactAxis( Equilibrium const &eq, double rGuess, double zGuess )
	{
		double r = rGuess;
		double z = zGuess;
		double const step = 1.0e-5;

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );

			double a0 = 0.0;
			double a1 = 0.0;
			double b0 = 0.0;
			double b1 = 0.0;
			double hessian[ 2 ][ 2 ];
			eq.gradPsi( r + step, z, a0, b0 );
			eq.gradPsi( r - step, z, a1, b1 );
			hessian[ 0 ][ 0 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 0 ] = ( b0 - b1 )/( 2.0*step );
			eq.gradPsi( r, z + step, a0, b0 );
			eq.gradPsi( r, z - step, a1, b1 );
			hessian[ 0 ][ 1 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 1 ] = ( b0 - b1 )/( 2.0*step );

			double const det = hessian[ 0 ][ 0 ]*hessian[ 1 ][ 1 ]
			                   - hessian[ 0 ][ 1 ]*hessian[ 1 ][ 0 ];
			r += -(  hessian[ 1 ][ 1 ]*gr - hessian[ 0 ][ 1 ]*gz )/det;
			z += -( -hessian[ 1 ][ 0 ]*gr + hessian[ 0 ][ 0 ]*gz )/det;
		}

		ExactAxis axis;
		axis.r = r;
		axis.z = z;
		axis.psi = eq.psi( r, z );
		return axis;
	}

	/// psi at normalised flux @a fraction, taken from the CLOSED FORM so that it
	/// is the same level at every mesh and every degree.
	///
	/// THE NORMALISATION IS THE REFERENCE HEADER'S: psi = 0 on the separatrix
	/// and psi_ax on the axis, so psi( Psi_N ) = psi_ax ( 1 - Psi_N ). It is NOT
	/// FluxSurfaceConvergence.cpp's, which interpolates between the axis and the
	/// smallest psi on its box boundary -- that one moves with the box and could
	/// not be compared against a tabulated reference at all.
	double levelAt( ExactAxis const &axis, double fraction )
	{
		return axis.psi*( 1.0 - fraction );
	}

	/// The reference: the same averages on the EXACT field, at enough angles to
	/// be converged to round-off. Never touches psi_h.
	meq::analytic::SurfaceQuadrature reference( Equilibrium const &eq,
	                                            ExactAxis const &axis,
	                                            double level, int angles = 1024 )
	{
		return meq::analytic::surfaceQuadrature( eq, axis.r, axis.z, level,
		                                         angles, 2.0 );
	}

	double referenceInverseRSquared( meq::analytic::SurfaceQuadrature const &s )
	{
		return s.average( []( meq::analytic::SurfacePoint const &p )
		{
			return 1.0/( p.r*p.r );
		} );
	}

	/// One solve, kept alive, with the post-processing done so that both
	/// candidate potentials are available. The harness's measure() destroys its
	/// solver on the way out, which is right for an error norm and no use here:
	/// the tracer borrows the potential and the flux and needs them to outlive
	/// the measurement. Member order is load bearing and the class is
	/// non-copyable because the coefficients capture `this`. Lifted from
	/// FluxSurfaceConvergence.cpp.
	class SolvedEquilibrium
	{
		public:
			SolvedEquilibrium( Equilibrium const &eqIn, Rectangle const &boxIn,
			                   int orderIn, int n )
				: eq( eqIn ),
				  mesh( meq::tests::makeMesh( boxIn, n ) ),
				  sourceCoeff( [ this ]( mfem::Vector const &x )
				  {
					  return eq.f( x( 0 ), x( 1 ), 0.0 );
				  } ),
				  psiCoeff( [ this ]( mfem::Vector const &x )
				  {
					  return eq.psi( x( 0 ), x( 1 ) );
				  } ),
				  solver( mesh, orderIn )
			{
				solver.setSource( sourceCoeff );
				solver.setBoundaryData( psiCoeff );
				solver.solve();
				solver.postProcess();
			}

			SolvedEquilibrium( SolvedEquilibrium const & ) = delete;
			SolvedEquilibrium &operator=( SolvedEquilibrium const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			meq::CriticalPoint axis() const
			{
				meq::CriticalPointFinder finder( solver );
				return finder.findAxis();
			}

		private:
			Equilibrium eq;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			meq::GradShafranovSolver solver;
	};

	/// A traced contour and its angle fit.
	///
	/// THERE USED TO BE A THIRD MEMBER HERE, the relative tolerance a ladder had
	/// settled on, and it is gone because the thing it worked around is fixed.
	/// AngleParametrisation::stalledRays and worstResidual are the honest
	/// signals and they are already in the fit.
	struct Fitted
	{
		meq::Contour contour;
		meq::AngleParametrisation fit;
	};

	/*
	 * TRACE, THEN FIT AT THE TIGHTEST TOLERANCE THE FIELD WILL SUPPORT.
	 *
	 * THIS LADDER IS A WORKAROUND AND IT IS LABELLED AS ONE. ContourTracer::
	 * fitByAngle() runs a safeguarded Newton along each ray and THROWS when it
	 * cannot reach | psi_h - c | <= tolerance x scale. On a discontinuous field
	 * that is not always attainable: { psi_h = c } is a union of per-element
	 * arcs offset by the DG jump, so a ray crossing a face where c falls inside
	 * the jump has NO point on it with psi_h = c at all, and no tolerance
	 * tighter than the jump can be met. The tracer's own corrector already
	 * handles exactly this -- FluxSurfaces.hpp documents keeping the best
	 * iterate and accepting after four non-improving steps -- and the ray
	 * Newton does not.
	 *
	 * MEASURED, so that this is a number rather than a worry: at k = 1 on the
	 * raw pairing the fit needs a tolerance of 1e-2 at n = 12 falling to 1e-5 at
	 * n = 96, tracking its own DG jump; at k >= 2, and at k = 1 with psi*, the
	 * tightest 1e-12 works at every mesh in this file. It also gets likelier
	 * with the ANGLE COUNT, simply because more rays are more chances: the same
	 * contour that fits at N = 256 can fail at N = 4096.
	 *
	 * The ladder asks for the tightest tolerance that works and reports it, so
	 * a table always says what precision the level set was actually located to.
	 * The right fix is in the tracer and is reported rather than done here.
	 */
	/// Trace one closed surface and fit it by angle, at the tracer's own
	/// tolerance and with no retry.
	///
	/// THIS USED TO BE A TOLERANCE LADDER AND THE LADDER IS GONE. fitByAngle()
	/// once threw where its ray Newton could not meet a tolerance that, on a
	/// discontinuous field, is sometimes UNATTAINABLE -- a ray crossing a face
	/// where c falls inside the jump has no point on it with psi_h = c at all.
	/// This function worked round that by loosening the tolerance a decade at a
	/// time until the fit succeeded, and reporting which rung it landed on.
	///
	/// The tracer now does what its own contour corrector always did: it keeps
	/// the best iterate the ray reached, accepts it, and counts it in
	/// AngleParametrisation::stalledRays with worstResidual saying how close it
	/// got. So the first rung always succeeds and the ladder is dead.
	///
	/// REMOVING IT IS NOT ONLY TIDYING. The ladder caught std::runtime_error,
	/// which is what fitByAngle throws for every OTHER reason too -- a ray on
	/// which every evaluation left the field, a vanishing flux, a surface that
	/// is not star-shaped. Those were being swallowed, retried at eight
	/// loosening tolerances, and finally reported as the ladder's own guess at
	/// what went wrong. A specific diagnosis was being converted into a vague
	/// one. Now they propagate.
	Fitted fittedSurface( meq::ContourTracer &tracer,
	                      meq::CriticalPoint const &axis, double level,
	                      std::size_t angles )
	{
		Fitted out;
		out.contour = tracer.traceFromAxis( level, axis );
		if ( !out.contour.closed() )
			throw std::runtime_error(
				"fittedSurface: the contour did not close, so there is no closed "
				"surface to average over" );

		out.fit = tracer.fitByAngle( out.contour, axis, angles );
		return out;
	}

	/// Trace one surface and build its quadrature by the equispaced-angle
	/// route, at the tightest tolerance the field supports.
	meq::SurfaceAverages averagesAt( meq::ContourTracer &tracer,
	                                 meq::CriticalPoint const &axis,
	                                 double level, std::size_t angles,
	                                 double *worstResidualOut = nullptr )
	{
		Fitted const fitted = fittedSurface( tracer, axis, level, angles );
		if ( worstResidualOut != nullptr )
			*worstResidualOut = fitted.fit.worstResidual;
		return meq::surfaceAverages( tracer, fitted.fit );
	}

	double relative( double value, double target )
	{
		return std::abs( value - target )/std::abs( target );
	}

}

/*
 * THE REFERENCE HAS NOT MOVED.
 *
 * No solve, no mesh, no tracer -- this is the exact field alone, and it exists
 * because every rate below is measured against these three numbers. If the
 * fixture's coefficients, the ray solve or the metric identity in
 * tests/analytic/FluxSurfaceReference.hpp ever change, this fails HERE rather
 * than turning up as a mysterious floor in a convergence table twenty lines
 * further down.
 *
 * The values are the ones recorded in that file's own header, converged to the
 * last printed digit at 256 angles. They are a reference VALUE and not a closed
 * form: there is no closed form for these quantities on this fixture, because
 * psi is elementary and an integral over a contour of it is not.
 */
BOOST_AUTO_TEST_CASE( theConvergedReferenceReproducesItsRecordedValues )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.0, 0.0 );

	std::printf( "\n  The reference, on the EXACT field of nstx(). A converged "
	             "VALUE, not a closed form.\n" );
	std::printf( "  axis ( %.12f, %.12f ), psi_ax = %.12e\n",
	             axis.r, axis.z, axis.psi );

	// The axis itself, first: every level below is built from it.
	BOOST_TEST( std::abs( axis.r - 1.318167937714 ) < 1.0e-9 );
	BOOST_TEST( std::abs( axis.z - 0.011088725858 ) < 1.0e-9 );
	BOOST_TEST( std::abs( axis.psi + 2.662896051834e-01 ) < 1.0e-12 );

	struct Recorded
	{
		double fraction;
		double vPrime;
		double inverseRSquared;
	};

	// tests/analytic/FluxSurfaceReference.hpp, its own header table.
	std::vector<Recorded> const recorded = {
		{ 0.25, 6.71413847786385e+01, 6.84662650462071e-01 },
		{ 0.50, 7.42523454479159e+01, 8.71200178854720e-01 },
		{ 0.75, 8.81259886962482e+01, 1.29438587807901e+00 }
	};

	std::printf( "  %6s %24s %14s %24s %14s %10s\n", "Psi_N", "V'", "vs recorded",
	             "< R^-2 >", "vs recorded", "transv" );

	for ( Recorded const &want : recorded )
	{
		double const level = levelAt( axis, want.fraction );
		meq::analytic::SurfaceQuadrature const s = reference( eq, axis, level );
		double const inverse = referenceInverseRSquared( s );

		std::printf( "  %6.2f %24.14e %14.2e %24.14e %14.2e %10.3f\n",
		             want.fraction, s.vPrime, relative( s.vPrime, want.vPrime ),
		             inverse, relative( inverse, want.inverseRSquared ),
		             s.transversality );

		BOOST_TEST( relative( s.vPrime, want.vPrime ) < 1.0e-12,
		            "V' at Psi_N = " << want.fraction << " reads " << s.vPrime
		            << " against the " << want.vPrime << " recorded in "
		            << "FluxSurfaceReference.hpp. Every rate in this file is "
		            << "measured against that number, so it is pinned here "
		            << "rather than trusted" );
		BOOST_TEST( relative( inverse, want.inverseRSquared ) < 1.0e-12,
		            "< R^-2 > at Psi_N = " << want.fraction << " reads "
		            << inverse << " against the recorded "
		            << want.inverseRSquared );
	}
	std::fflush( stdout );
}

/*
 * THE IDENTITY, ON THE EXACT FIELD, AND THE INSTRUMENT THAT WOULD HIDE IT.
 *
 * ( 1 / V' ) d/dpsi ( V' < | grad psi |^2 / R^2 > ) = - < F / R^2 >, with F the
 * source the solver is fed. Here it is evaluated on the analytic psi, where
 * every average is exact to round-off, so what is left is the d/dpsi ALONE.
 *
 * THE POINT OF THE TEST IS THE SECOND COLUMN. A plain central difference
 * carries its own O( step^2 ) truncation, and that floors the agreement however
 * exact everything else is -- the floor is the INSTRUMENT, not the identity, and
 * an identity checked at 1e-5 would pass with a real defect underneath it. This
 * is the third place in the tree where that applies; src/meq/Zernike.hpp's
 * derivative test is the first and INVERSION-PLAN.md IN-3 the second.
 */
BOOST_AUTO_TEST_CASE( theExactIdentityNeedsRichardsonAndNotACentralDifference )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.0, 0.0 );
	double const level = levelAt( axis, 0.25 );
	int const angles = 256;

	auto weighted = [ & ]( double c )
	{
		meq::analytic::SurfaceQuadrature const s = reference( eq, axis, c, angles );
		return s.vPrime*s.average( []( meq::analytic::SurfacePoint const &p )
		{
			return ( p.psiR*p.psiR + p.psiZ*p.psiZ )/( p.r*p.r );
		} );
	};

	meq::analytic::SurfaceQuadrature const here = reference( eq, axis, level,
	                                                         angles );
	double const rightHandSide = here.average(
		[ & ]( meq::analytic::SurfacePoint const &p )
		{
			return -eq.f( p.r, p.z, level )/( p.r*p.r );
		} );

	std::printf( "\n  The averaged Grad-Shafranov identity on the EXACT field, "
	             "Psi_N = 0.25, %d angles\n", angles );
	std::printf( "  - < F / R^2 > = %.14e, so the residuals below are relative "
	             "to an O(1) quantity\n", rightHandSide );
	std::printf( "  %12s %16s %16s %10s\n", "step", "Richardson", "central",
	             "ratio" );

	std::vector<double> richardsonColumn;
	std::vector<double> centralColumn;
	std::vector<double> steps;

	for ( double fraction : { 0.05, 0.025, 0.0125, 0.00625 } )
	{
		double const step = fraction*std::abs( axis.psi );
		double const coarse = ( weighted( level + step )
		                        - weighted( level - step ) )/( 2.0*step );
		double const fine = ( weighted( level + 0.5*step )
		                      - weighted( level - 0.5*step ) )/step;

		double const richardson = ( 4.0*fine - coarse )/3.0/here.vPrime
		                          - rightHandSide;
		double const central = coarse/here.vPrime - rightHandSide;

		steps.push_back( step );
		richardsonColumn.push_back( std::abs( richardson ) );
		centralColumn.push_back( std::abs( central ) );

		std::printf( "  %12.6f %16.4e %16.4e %10.1f\n", step,
		             std::abs( richardson ), std::abs( central ),
		             std::abs( central )/std::abs( richardson ) );
	}
	std::fflush( stdout );

	// The central difference is second order in the step and the extrapolated
	// one is fourth. Both are asserted, because the SEPARATION is the finding
	// and a separation with only one side measured is not one.
	double const ratio = steps.front()/steps.back();
	double const centralRate = rate( centralColumn.front(), centralColumn.back(),
	                                 ratio );
	double const richardsonRate = rate( richardsonColumn.front(),
	                                    richardsonColumn.back(), ratio );

	std::printf( "  central difference converges at %.3f in the step, "
	             "Richardson at %.3f\n", centralRate, richardsonRate );

	bool const centralIsSecondOrder = centralRate > 2.0 - 0.2
	                                  && centralRate < 2.0 + 0.3;
	BOOST_TEST( centralIsSecondOrder,
	            "the plain central difference converges at " << centralRate
	            << " rather than the second order it is. IT IS THE CONTROL: if "
	            << "it ever converges as fast as the extrapolation there is "
	            << "nothing to extrapolate and this test is worthless" );

	BOOST_TEST( richardsonRate > centralRate + 1.0,
	            "Richardson extrapolation converges at " << richardsonRate
	            << " against the central difference's " << centralRate
	            << ", which is not the order it buys" );

	for ( std::size_t i = 0; i < steps.size(); ++i )
		BOOST_TEST( richardsonColumn[ i ] < 0.05*centralColumn[ i ],
		            "at step " << steps[ i ] << " the extrapolated residual is "
		            << richardsonColumn[ i ] << " and the plain one "
		            << centralColumn[ i ] << ". The identity is exact on this "
		            << "field, so the difference between those two columns is "
		            << "the INSTRUMENT and nothing else" );
}

/*
 * SPECTRAL IN THE ANGLE COUNT, AND THE DIFFERENCED METRIC KEPT AS A LIVE
 * CONTROL.
 *
 * INVERSION-PLAN.md section 3.2 warns that it is easy to build a spectrally
 * accurate rule and then feed it a SECOND-ORDER Jacobian obtained by
 * differencing neighbouring node positions, at which point the whole scheme is
 * second order and nothing in its output says so. IN-1 measured that on an ARC
 * LENGTH. AN AVERAGE IS A DIFFERENT INTEGRAND AND IS ENTITLED TO ITS OWN
 * MEASUREMENT, which is the whole reason this table has four columns rather than
 * inheriting IN-1's two.
 *
 * The interesting half is < R^-2 >, because there is a plausible argument that
 * the trap should not bite it: an average is a RATIO of two integrals over the
 * same metric, so a bad metric ought to cancel. Measured, it does not -- what
 * the ratio buys is a smaller constant and nothing whatever in the order.
 *
 * The reference here is the SAME discrete surface at a large angle count, not
 * the analytic one: this table measures the quadrature and not the
 * discretisation, and mixing the two would attribute the field error to the
 * rule.
 */
BOOST_AUTO_TEST_CASE( theQuadratureIsOnlyAsGoodAsItsMetric )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = nstxBox();
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );
	double const level = levelAt( exact, 0.25 );

	SolvedEquilibrium solved( eq, box, 3, 32 );
	meq::CriticalPoint const axis = solved.axis();

	meq::ContourTracer tracer( solved.theSolver() );
	tracer.setTargetTurn( 0.05 );

	// ONE TOLERANCE FOR THE WHOLE SWEEP, chosen by the ladder at the LARGEST
	// angle count -- which is where the ray fit is most likely to meet the DG
	// jump, since more rays are more chances. Every row below is then located to
	// the same precision, which is what makes the columns comparable at all: a
	// tolerance that moved with N would put the ladder into the rate.
	std::size_t const settled = 2048;
	Fitted const settledFit = fittedSurface( tracer, axis, level, settled );
	meq::Contour const &contour = settledFit.contour;
	BOOST_TEST( contour.closed(), "the contour to be averaged did not close" );

	meq::SurfaceAverages const converged = meq::surfaceAverages( tracer,
	                                                             settledFit.fit );

	std::printf( "\n  Flux-surface averages of { psi* = c } against the angle "
	             "count, k = 3, n = 32, Psi_N = 0.25\n" );
	std::printf( "  settled at N = %d: V' = %.12e, < R^-2 > = %.12e\n",
	             static_cast<int>( settled ), converged.vPrime,
	             converged.inverseRSquared() );
	std::printf( "  the DG jump on this contour is %.2e in psi; that is where "
	             "the pointwise columns floor.\n  The settled fit stalled on %d "
	             "of %d rays, worst residual %.2e\n", contour.worstFaceJump,
	             settledFit.fit.stalledRays,
	             static_cast<int>( settled ), settledFit.fit.worstResidual );
	std::printf( "  %6s %14s %7s %14s %7s | %14s %7s %14s %7s\n",
	             "N", "V' from q", "rate", "V' diffed", "rate",
	             "<R^-2> from q", "rate", "<R^-2> diffed", "rate" );

	std::vector<std::size_t> const counts = { 16, 32, 64, 128 };
	std::vector<double> vFlux;
	std::vector<double> vDiff;
	std::vector<double> aFlux;
	std::vector<double> aDiff;

	double const settledInverse = converged.inverseRSquared();

	auto inverseSquare = []( meq::SurfaceNode const &node )
	{
		return 1.0/( node.r*node.r );
	};

	for ( std::size_t c = 0; c < counts.size(); ++c )
	{
		meq::SurfaceAverages const surface =
			meq::surfaceAverages( tracer, tracer.fitByAngle( contour, axis,
			                                                 counts[ c ] ) );

		BOOST_TEST( surface.differencedAvailable(),
		            "the equispaced-angle route must carry its differenced "
		            "control, or the comparison this table is built on has no "
		            "second column" );

		double const a = relative( surface.vPrime, converged.vPrime );
		double const b = relative( surface.vPrimeDifferenced, converged.vPrime );
		double const d = relative( surface.average( inverseSquare ),
		                           settledInverse );
		double const e = relative( surface.averageDifferenced( inverseSquare ),
		                           settledInverse );

		vFlux.push_back( a );
		vDiff.push_back( b );
		aFlux.push_back( d );
		aDiff.push_back( e );

		if ( c == 0 )
			std::printf( "  %6d %14.4e %7s %14.4e %7s | %14.4e %7s %14.4e %7s\n",
			             static_cast<int>( counts[ c ] ), a, "-", b, "-", d, "-",
			             e, "-" );
		else
			std::printf( "  %6d %14.4e %7.3f %14.4e %7.3f | %14.4e %7.3f "
			             "%14.4e %7.3f\n",
			             static_cast<int>( counts[ c ] ), a,
			             rate( vFlux[ c - 1 ], a, 2.0 ), b,
			             rate( vDiff[ c - 1 ], b, 2.0 ), d,
			             rate( aFlux[ c - 1 ], d, 2.0 ), e,
			             rate( aDiff[ c - 1 ], e, 2.0 ) );
	}
	std::fflush( stdout );

	double const ratio = static_cast<double>( counts.back() )/counts.front();

	struct Column
	{
		char const *name;
		std::vector<double> const &pointwise;
		std::vector<double> const &differenced;
	};

	std::vector<Column> const columns = {
		{ "V'", vFlux, vDiff },
		{ "< R^-2 >", aFlux, aDiff }
	};

	for ( Column const &column : columns )
	{
		double const pointwiseRate = rate( column.pointwise.front(),
		                                   column.pointwise.back(), ratio );
		double const differencedRate = rate( column.differenced.front(),
		                                     column.differenced.back(), ratio );

		std::printf( "  %-10s pointwise metric %.3f, differenced metric %.3f, "
		             "separation at N = %d is %.1e\n", column.name,
		             pointwiseRate, differencedRate,
		             static_cast<int>( counts.back() ),
		             column.differenced.back()/column.pointwise.back() );

		BOOST_TEST( differencedRate < 3.0,
		            column.name << ": the differenced Jacobian converges at "
		            << differencedRate << ". THAT COLUMN IS THE CONTROL -- it is "
		            << "the trap of INVERSION-PLAN.md section 3.2, a spectrally "
		            << "accurate rule fed a second-order metric -- and if it "
		            << "ever converges as fast as the pointwise one the "
		            << "comparison is empty and this test is worthless" );

		BOOST_TEST( differencedRate > 2.0 - 0.3,
		            column.name << ": the differenced Jacobian converges at "
		            << differencedRate << " rather than the second order a "
		            << "central difference gives" );

		BOOST_TEST( pointwiseRate > differencedRate + 1.0,
		            column.name << ": the pointwise metric from q converges at "
		            << pointwiseRate << " against the differenced one's "
		            << differencedRate << ". The whole content of the stage is "
		            << "that taking the Jacobian from q rather than from a "
		            << "difference is worth more than a constant" );

		BOOST_TEST( column.pointwise.back() < 1.0e-3*column.differenced.back(),
		            column.name << ": at N = " << counts.back()
		            << " the pointwise metric gives " << column.pointwise.back()
		            << " against the differenced one's "
		            << column.differenced.back() << ", which is not the "
		            << "separation the table is supposed to demonstrate" );
	}
}


/*
 * IN h, AGAINST THE CONVERGED REFERENCE ON THE EXACT FIELD -- ERROR ( a )
 * ALONE, AND THE ONE RESULT OF THIS STAGE THAT WAS NOT EXPECTED.
 *
 * The angle count is fixed and large, so the quadrature error of the previous
 * table is out of the way and what is left is the discretisation: how far
 * { psi_h = c } is from { psi = c }, and how far r | q_h | is from | grad psi |.
 * BOTH enter, and by different routes -- the first sets where the nodes are and
 * the second sets what they weigh.
 *
 * TWO COLUMNS, AND THE POINT IS THAT THEY HAVE THE SAME ORDER. ContourTracer's
 * default is Potential::PostProcessed, which roots psi* and takes its flux from
 * the reconstruction; Potential::Raw roots psi_h with q_h. psi* converges at
 * k+2 where psi_h converges at k+1, so the natural expectation -- and the one
 * this test was written to check rather than to assume -- is that the
 * post-processed pairing buys an order in the averages too.
 *
 * IT DOES NOT, AND THE REASON IS STRUCTURAL. The weight is
 * 2 pi R dl / | grad psi | and | grad psi | is r | q |: a FLUX, not a potential.
 * q_h converges at k+1 and so does q*, the reconstruction buying its extra
 * order in the POTENTIAL and not in the field that divides every weight in this
 * file. So the level set improves by an order and the weight does not, and the
 * average inherits the worse of the two. What psi* buys here is a CONSTANT --
 * measured below, a factor of a few -- and not an order. A constant is not an
 * order and the two must not be confused; this is the same shape of finding as
 * the band continuation of B in CLAUDE.md, where the flux carries the
 * potential's order and its own gradient cannot.
 *
 * SO k+1 IS ASSERTED ON BOTH COLUMNS, and Potential::Raw is what pins it.
 *
 * THE ERROR IS A SIGNED QUANTITY AND ITS MAGNITUDE IS NOT MONOTONE IN h, which
 * is measured rather than assumed: a single surface gives rates that dip and
 * rebound -- one pair reading 0.79 and the next 5.93 on the same sequence --
 * because the integrand's error changes sign around the contour and the
 * cancellation is not the same at every mesh. The rate is therefore taken over
 * the ROOT-MEAN-SQUARE of four surfaces rather than over one, which is also the
 * quantity a consumer of profiles actually cares about, and it is taken across
 * the whole sequence rather than pairwise. That is the two-tier pattern of
 * ExtensionConvergence.cpp and it is here for the same reason: which elements a
 * contour passes through is not a smooth function of h.
 *
 * THE LEVELS COME FROM THE CLOSED FORM and are therefore the same at every mesh
 * and every degree. A level read off the discrete solution would move with h,
 * and the reference it is compared against would move with it, so the rate would
 * be measured against a moving target and would mean nothing.
 */
BOOST_AUTO_TEST_CASE( theAveragesConvergeInTheMeshAgainstTheExactField )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = nstxBox();
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );

	std::vector<int> const sizes = { 12, 24, 48, 96 };
	std::vector<double> const fractions = { 0.15, 0.25, 0.35, 0.50 };
	std::size_t const angles = 256;

	// The reference values, once: they do not depend on the mesh.
	std::vector<double> levels;
	std::vector<double> wantVPrime;
	std::vector<double> wantInverse;
	for ( double fraction : fractions )
	{
		double const level = levelAt( exact, fraction );
		meq::analytic::SurfaceQuadrature const want = reference( eq, exact, level );
		levels.push_back( level );
		wantVPrime.push_back( want.vPrime );
		wantInverse.push_back( referenceInverseRSquared( want ) );
	}

	std::printf( "\n  Flux-surface averages against the CONVERGED REFERENCE on "
	             "the exact field.\n  A reference VALUE, not a closed form: "
	             "these quantities have none on this fixture.\n"
	             "  nstx() on [%.2f,%.2f]x[%.2f,%.2f], %d angles, RMS over "
	             "Psi_N =", box.rMin, box.rMax, box.zMin, box.zMax,
	             static_cast<int>( angles ) );
	for ( double fraction : fractions )
		std::printf( " %.2f", fraction );
	std::printf( "\n" );

	struct Choice
	{
		char const *name;
		meq::Potential which;
	};

	std::vector<Choice> const choices = {
		{ "raw psi_h, q_h", meq::Potential::Raw },
		{ "psi*, q*      ", meq::Potential::PostProcessed }
	};

	for ( int order = 1; order <= 3; ++order )
	{
		double const expected = order + 1.0;

		std::vector<std::vector<double>> vErrors( choices.size() );
		std::vector<std::vector<double>> aErrors( choices.size() );

		std::printf( "\n  k = %d, expecting k+1 = %.0f in both columns\n",
		             order, expected );
		std::printf( "  %6s %-16s %13s %7s %13s %7s %9s\n", "n", "pairing",
		             "rms V' error", "rate", "rms <R^-2>", "rate", "worst res" );

		for ( std::size_t i = 0; i < sizes.size(); ++i )
		{
			// ONE SOLVE PER MESH, shared by both pairings: they differ in which
			// fields the tracer reads and in nothing else, so solving twice
			// would measure the same discretisation twice.
			SolvedEquilibrium solved( eq, box, order, sizes[ i ] );
			meq::CriticalPoint const axis = solved.axis();

			for ( std::size_t c = 0; c < choices.size(); ++c )
			{
				meq::ContourTracer tracer( solved.theSolver(), choices[ c ].which );

				double sumV = 0.0;
				double sumA = 0.0;
				double worstTolerance = 0.0;

				for ( std::size_t j = 0; j < levels.size(); ++j )
				{
					double tolerance = 0.0;
					meq::SurfaceAverages const surface =
						averagesAt( tracer, axis, levels[ j ], angles, &tolerance );
					worstTolerance = std::max( worstTolerance, tolerance );

					double const v = relative( surface.vPrime, wantVPrime[ j ] );
					double const a = relative( surface.inverseRSquared(),
					                           wantInverse[ j ] );
					sumV += v*v;
					sumA += a*a;
				}

				double const v = std::sqrt( sumV/levels.size() );
				double const a = std::sqrt( sumA/levels.size() );
				vErrors[ c ].push_back( v );
				aErrors[ c ].push_back( a );

				if ( i == 0 )
					std::printf( "  %6d %-16s %13.4e %7s %13.4e %7s %9.0e\n",
					             sizes[ i ], choices[ c ].name, v, "-", a, "-",
					             worstTolerance );
				else
					std::printf( "  %6d %-16s %13.4e %7.3f %13.4e %7.3f %9.0e\n",
					             sizes[ i ], choices[ c ].name, v,
					             rate( vErrors[ c ][ i - 1 ], v, 2.0 ), a,
					             rate( aErrors[ c ][ i - 1 ], a, 2.0 ),
					             worstTolerance );
			}
		}
		std::fflush( stdout );

		double const ratio = static_cast<double>( sizes.back() )/sizes.front();

		for ( std::size_t c = 0; c < choices.size(); ++c )
		{
			double const vRate = rate( vErrors[ c ].front(), vErrors[ c ].back(),
			                           ratio );
			double const aRate = rate( aErrors[ c ].front(), aErrors[ c ].back(),
			                           ratio );

			std::printf( "  k = %d, %s over the sequence: V' %.3f, "
			             "< R^-2 > %.3f\n", order, choices[ c ].name, vRate,
			             aRate );

			BOOST_TEST( vRate > expected - sequenceSlack,
			            "k = " << order << ", " << choices[ c ].name
			            << ": V' converges at " << vRate << " rather than the "
			            << expected << " the flux's own order gives. Both the "
			            << "level set and the weight enter, and the weight is "
			            << "1/( r | q | ) -- so a lost order here is a lost "
			            << "order in q" );

			BOOST_TEST( aRate > expected - sequenceSlack,
			            "k = " << order << ", " << choices[ c ].name
			            << ": < R^-2 > converges at " << aRate
			            << " rather than " << expected );

			// A RATE IS BLIND TO A CONSTANT FACTOR, which is the failure
			// CLAUDE.md's testing stance is organised around, so the absolute
			// error on the finest mesh is asserted as well.
			BOOST_TEST( vErrors[ c ].back() < 1.0e-4,
			            "k = " << order << ", " << choices[ c ].name
			            << ": V' on the finest mesh is wrong by "
			            << vErrors[ c ].back()
			            << " relative, which is a rate about the wrong number" );
			BOOST_TEST( aErrors[ c ].back() < 1.0e-4,
			            "k = " << order << ", " << choices[ c ].name
			            << ": < R^-2 > on the finest mesh is wrong by "
			            << aErrors[ c ].back() << " relative" );
		}

		// WHAT psi* BUYS, AS A NUMBER. It is a constant and not an order; see
		// the comment on this test. Printed at every degree so that a change of
		// character -- the post-processed column pulling an order ahead, which
		// is what would happen if a k+2 flux ever became available -- shows up
		// as a trend rather than being invisible.
		std::printf( "  k = %d, psi* over psi_h on the finest mesh: V' x%.2f, "
		             "< R^-2 > x%.2f\n", order,
		             vErrors[ 0 ].back()/vErrors[ 1 ].back(),
		             aErrors[ 0 ].back()/aErrors[ 1 ].back() );
	}
}

/*
 * THE AVERAGED GRAD-SHAFRANOV IDENTITY ON THE DISCRETE FIELD, WHICH NEEDS NO
 * REFERENCE VALUE AT ALL.
 *
 *     ( 1 / V' ) d/dpsi ( V' < | grad psi |^2 / R^2 > )  =  - < F / R^2 >
 *
 * Three averages check each other with nothing but the equation, and the
 * right-hand side is written with THE F THE SOLVER IS FED rather than as a
 * hand-derived -mu_0 p' - g g' < R^-2 >. The two are equal for Solov'ev, but the
 * second is re-derived by hand and so is not independent of the hand that
 * derived it; the first applies to every fixture in tests/analytic/. It is the
 * discipline SolovievEquilibrium::deltaStarFD() follows and it is why the F is a
 * callable in the facility's signature.
 *
 * THE PROPERTY TO ASSERT IS CONVERGENCE, NOT SMALLNESS. The identity is exact
 * for psi and holds only to the discretisation error for psi_h, so a residual
 * that is merely small says nothing; a residual that falls at the field's own
 * order says the three averages agree with each other in the limit for the
 * right reason.
 *
 * AND THE SECOND COLUMN IS THE INSTRUMENT, MADE VISIBLE. The d/dpsi is taken by
 * Richardson extrapolation in the answer column and by a plain central
 * difference in the control. The control's own O( step^2 ) truncation FLOORS it
 * -- flat to three figures across a sixteen-fold refinement at k = 3 -- while
 * the extrapolated column keeps falling. An identity checked with that
 * instrument would pass at 1e-5 with a real defect underneath it, which is the
 * whole reason INVERSION-PLAN.md puts the requirement in capitals.
 *
 * THE STEP IS A REAL CHOICE AND IT IS NOT MONOTONE. Too large and the
 * extrapolation's own truncation dominates; too small and the difference of two
 * nearly equal integrals divides the surfaces' own jump noise by the step and
 * amplifies it. Measured at k = 2, n = 96: 1.5e-08 at 2% of | psi_ax | and
 * 3.4e-07 at 1%, the smaller step being twenty times WORSE. Two percent is what
 * is used, and the k = 3 column then meets the extrapolation's own floor rather
 * than the discretisation -- which is why the rate is asserted at k = 2 and the
 * k = 3 row is printed with that explanation beside it.
 */
BOOST_AUTO_TEST_CASE( theAveragedIdentityHoldsOnTheDiscreteFieldAndConverges )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = nstxBox();
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );
	double const level = levelAt( exact, 0.25 );
	double const step = 0.02*std::abs( exact.psi );
	std::size_t const angles = 256;

	std::vector<int> const sizes = { 12, 24, 48, 96 };

	// THE F THE SOLVER IS FED, forwarded rather than re-derived: SolvedEquilibrium
	// hands the solver eq.f( r, z, 0 ) as its source coefficient, and this is the
	// same function.
	auto source = [ &eq ]( double r, double z, double psi )
	{
		return eq.f( r, z, psi );
	};

	std::printf( "\n  The averaged Grad-Shafranov identity on { psi* = c }, "
	             "Psi_N = 0.25, %d angles, step %.5f\n",
	             static_cast<int>( angles ), step );

	for ( int order = 2; order <= 3; ++order )
	{
		std::vector<double> richardson;
		std::vector<double> central;
		double scale = 0.0;

		std::printf( "  k = %d  %6s %15s %7s %15s %7s %11s\n", order, "n",
		             "Richardson", "rate", "central diff", "rate", "- <F/R^2>" );

		for ( std::size_t i = 0; i < sizes.size(); ++i )
		{
			SolvedEquilibrium solved( eq, box, order, sizes[ i ] );
			meq::CriticalPoint const axis = solved.axis();
			meq::ContourTracer tracer( solved.theSolver() );

			// The identity traces five surfaces of its own, so the ladder is run
			// on the base level first and the tolerance it settles on is left in
			// place for all of them.
			double tolerance = 0.0;
			averagesAt( tracer, axis, level, angles, &tolerance );
			tracer.setTolerance( tolerance );

			meq::AveragedEquationResidual const extrapolated =
				meq::averagedGradShafranovResidual( tracer, axis, level, angles,
				                                    step, source );
			meq::AveragedEquationResidual const plain =
				meq::averagedGradShafranovResidual(
					tracer, axis, level, angles, step, source,
					meq::FluxDerivative::CentralDifference );
			tracer.setTolerance( 1.0e-12 );

			BOOST_TEST( !extrapolated.extended,
			            "the fitted path has no band, so no surface in this "
			            "identity may report one" );

			richardson.push_back( std::abs( extrapolated.residual ) );
			central.push_back( std::abs( plain.residual ) );
			scale = std::abs( extrapolated.rightHandSide );

			if ( i == 0 )
				std::printf( "  %8s %6d %15.4e %7s %15.4e %7s %11.5f\n", "",
				             sizes[ i ], richardson.back(), "-", central.back(),
				             "-", extrapolated.rightHandSide );
			else
				std::printf( "  %8s %6d %15.4e %7.3f %15.4e %7.3f %11.5f\n", "",
				             sizes[ i ], richardson.back(),
				             rate( richardson[ i - 1 ], richardson.back(), 2.0 ),
				             central.back(),
				             rate( central[ i - 1 ], central.back(), 2.0 ),
				             extrapolated.rightHandSide );
		}
		std::fflush( stdout );

		double const ratio = static_cast<double>( sizes.back() )/sizes.front();
		double const extrapolatedRate = rate( richardson.front(),
		                                      richardson.back(), ratio );
		double const plainRate = rate( central.front(), central.back(), ratio );

		std::printf( "  k = %d over the sequence: Richardson %.3f, central "
		             "%.3f; the right-hand side is %.4f\n", order,
		             extrapolatedRate, plainRate, scale );

		// THE CONTROL, FIRST: a central difference must STOP converging, or
		// there is nothing for the extrapolation to buy and this table is
		// worthless.
		BOOST_TEST( central.back() > 0.25*central[ sizes.size() - 2 ],
		            "k = " << order << ": the plain central difference went from "
		            << central[ sizes.size() - 2 ] << " to " << central.back()
		            << " over the last refinement, so it has NOT floored. THAT "
		            << "COLUMN IS THE CONTROL: it is supposed to be limited by "
		            << "its own O( step^2 ) truncation, and if it converges with "
		            << "the mesh instead then the instrument is not the limit "
		            << "and this comparison is empty" );

		BOOST_TEST( richardson.back() < 0.02*central.back(),
		            "k = " << order << ": on the finest mesh the extrapolated "
		            << "residual is " << richardson.back() << " and the plain "
		            << "one " << central.back() << ". The identity is the same "
		            << "identity in both, so the gap between them is the "
		            << "INSTRUMENT and nothing else" );

		BOOST_TEST( richardson.back() < 1.0e-7*scale,
		            "k = " << order << ": the identity residual on the finest "
		            << "mesh is " << richardson.back() << " against a "
		            << "right-hand side of " << scale << ". Three averages are "
		            << "supposed to check each other here with nothing but the "
		            << "equation" );

		if ( order == 2 )
		{
			// ASSERTED AT k = 2 ONLY, and the reason is measured rather than
			// convenient: at k = 3 the discrete residual on the finest mesh
			// falls to 2.8e-09, which is where the RICHARDSON extrapolation's
			// own O( step^4 ) truncation sits at this step -- 2.4e-09,
			// interpolated from the exact-field sweep in
			// theExactIdentityNeedsRichardsonAndNotACentralDifference. A rate
			// asserted there would be a rate about the instrument again, one
			// level up.
			BOOST_TEST( extrapolatedRate > order + 1.0 - sequenceSlack,
			            "k = " << order << ": the identity residual converges at "
			            << extrapolatedRate << " rather than at the field's own "
			            << order + 1.0 << ". A residual that is merely small "
			            << "says nothing; one that falls at the discretisation's "
			            << "order says the three averages agree for the right "
			            << "reason" );
		}
	}
}

/*
 * TWO INDEPENDENT EXTRACTIONS OF THE SAME SURFACE MUST AGREE -- AND THIS IS TWO
 * LEGS OF THREE.
 *
 * The predictor-corrector trace with its cubic Hermite, and the ray-based angle
 * parametrisation, share the field, the level and the axis and NOTHING else:
 *
 *   nodes    the tracer's own accepted points        |   equispaced in poloidal
 *                                                   |   angle about the axis
 *   metric   the interpolant's own | dx/dt |        |   sqrt( rho'^2 + rho^2 )
 *                                                   |   from IN-1's identity
 *   rule     composite Gauss-Legendre per segment   |   periodic trapezoid
 *   hypothesis  none                                |   star-shapedness
 *
 * Agreement is worth more than either being plausible, and disagreement would
 * not say which was wrong -- which is precisely why both are kept rather than
 * one being chosen. THE THIRD LEG IS NOT BUILT: INVERSION-PLAN.md section 3.3's
 * implicit quadrature, a rule constructed on the level set with no curve
 * extracted at all, is deferred deliberately (Saye is for hyperrectangles and
 * meq is on triangles, so it would be Fries-Omerovic or moment fitting). Read
 * what follows as two of three and not as a closed case.
 *
 * HOW WELL THEY CAN AGREE IS SET BY THE DG JUMP, AND THAT IS THE FINDING OF
 * THIS TEST RATHER THAN A CAVEAT ON IT. { psi_h = c } is NOT one curve: it is a
 * union of per-element arcs offset from each other by the jump in psi_h across
 * every face. Two routes that place their nodes differently -- equispaced in
 * angle against the trace's own spacing along the curve -- therefore sample
 * different arcs at every face crossing, and no amount of quadrature closes
 * that. Measured at k = 2, n = 24: the two differ by 3.2e-07 in V' where the ray
 * fit is 3.5e-07 from the exact field. THEY AGREE TO ABOUT THEIR OWN ERROR, and
 * they cannot do better.
 *
 * WHICH IS WHY THE ASSERTION IS ON THE RATE. A gap the size of the
 * discretisation error is what two correct discretisations of the same integral
 * must produce; a gap that does NOT converge, or converges at a different order,
 * or is much larger than either route's own error, is a systematic difference of
 * definition -- a missing 2 pi R, a metric taken about the wrong point, a weight
 * with the gradient on the wrong side. Those are exactly the defects a
 * single-route convergence table cannot see, because they converge too.
 *
 * AND THE GAP IS TAKEN OVER FOUR SURFACES, for the reason the h table above
 * gives: the difference of two SIGNED errors cancels by different amounts at
 * different meshes, and on one surface it reads 0.99 for one pair and 2.89 for
 * the next on the same sequence. The root mean square over four surfaces is
 * smooth enough to carry a rate and is also the quantity a consumer of profiles
 * cares about.
 *
 * AND A THIRD, CHEAP CHECK THAT IS NOT A THIRD EXTRACTION. An integral over a
 * curve does not depend on how the curve is parametrised, so fitting the SAME
 * contour about a deliberately displaced ray origin must give the same V'. It is
 * a real check on the metric -- it is the invariance the identity
 * | dx/dtheta | = sqrt( rho'^2 + rho^2 ) exists to satisfy -- and it costs one
 * more fit. It does not make a third leg: the field, the curve and the rule are
 * all shared. Measured at 4e-10 and better over displacements up to a third of
 * the minor radius.
 */
BOOST_AUTO_TEST_CASE( theTwoExtractionsAgree )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = nstxBox();
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );

	std::vector<int> const sizes = { 12, 24, 48 };
	std::vector<double> const fractions = { 0.15, 0.25, 0.35, 0.50 };

	std::vector<double> levels;
	std::vector<double> wantVPrime;
	std::vector<double> wantInverse;
	for ( double fraction : fractions )
	{
		double const level = levelAt( exact, fraction );
		meq::analytic::SurfaceQuadrature const want = reference( eq, exact, level );
		levels.push_back( level );
		wantVPrime.push_back( want.vPrime );
		wantInverse.push_back( referenceInverseRSquared( want ) );
	}

	std::printf( "\n  Two extractions of { psi* = c }: the ray fit and the "
	             "traced Hermite. The third leg\n  -- implicit quadrature on "
	             "the level set, section 3.3 -- is NOT built, so this is two\n"
	             "  of three. Delta_s = h/16, and RMS over Psi_N =" );
	for ( double fraction : fractions )
		std::printf( " %.2f", fraction );
	std::printf( "\n" );

	for ( int order = 2; order <= 3; ++order )
	{
		std::printf( "  k = %d %6s %8s %13s %13s %13s %7s %13s %7s\n", order,
		             "n", "points", "angle vs ref", "Herm vs ref", "angle-Herm",
		             "rate", "<R^-2> gap", "rate" );

		std::vector<double> gaps;
		std::vector<double> inverseGaps;

		for ( std::size_t i = 0; i < sizes.size(); ++i )
		{
			SolvedEquilibrium solved( eq, box, order, sizes[ i ] );
			meq::CriticalPoint const axis = solved.axis();

			double sumAngle = 0.0;
			double sumHermite = 0.0;
			double sumGap = 0.0;
			double sumInverse = 0.0;
			int points = 0;

			for ( std::size_t j = 0; j < levels.size(); ++j )
			{
				meq::ContourTracer tracer( solved.theSolver() );

				// THE TRACE STEP IS TIED TO h, WHICH IS NOT A DETAIL.
				// INVERSION-PLAN.md section 8's fourth risk is "Delta_s chosen
				// once and left": the whole section 2 argument rests on the
				// representation error being tuned to the mesh, and a fixed step
				// reintroduces exactly the floor the method exists to avoid --
				// which then looks like the method failing rather than the
				// parameter being wrong. Measured: at Delta_s = h/4 the Hermite
				// route is 174 times further from the exact field than the ray
				// fit at k = 3, n = 48, entirely because of its own
				// O( Delta_s^4 ); at h/16 it is closer than the ray fit. Section
				// 2 predicts Delta_s <~ h^((k+1)/4) for the interpolation to be
				// subdominant and that is what this is, with the constant
				// measured rather than assumed. Curvature control and the
				// element ceiling are both off so that the step really is the
				// one asked for.
				double const spacing = box.width()/( 16.0*sizes[ i ] );
				tracer.setTargetTurn( 0.0 );
				tracer.setLocalStepCeiling( 0.0 );
				tracer.setStep( spacing );

				Fitted const fitted = fittedSurface( tracer, axis,
				                                              levels[ j ], 1024 );

				meq::SurfaceAverages const angle =
					meq::surfaceAverages( tracer, fitted.fit );
				meq::SurfaceAverages const hermite =
					meq::surfaceAverages( tracer, fitted.contour, 10 );
				tracer.setTolerance( 1.0e-12 );

				double const angleError = relative( angle.vPrime,
				                                    wantVPrime[ j ] );
				double const hermiteError = relative( hermite.vPrime,
				                                      wantVPrime[ j ] );
				double const gap = relative( hermite.vPrime, angle.vPrime );
				double const inverseGap = relative( hermite.inverseRSquared(),
				                                    angle.inverseRSquared() );

				sumAngle += angleError*angleError;
				sumHermite += hermiteError*hermiteError;
				sumGap += gap*gap;
				sumInverse += inverseGap*inverseGap;
				points = static_cast<int>( fitted.contour.points.size() );

				BOOST_TEST( relative( angle.arcLength(), hermite.arcLength() )
				            < 1.0e-4,
				            "k = " << order << ", n = " << sizes[ i ]
				            << ", Psi_N = " << fractions[ j ] << ": the two "
				            << "routes disagree about the arc length of the "
				            << "same curve by "
				            << relative( angle.arcLength(),
				                         hermite.arcLength() ) );
			}

			double const count = static_cast<double>( levels.size() );
			double const angleError = std::sqrt( sumAngle/count );
			double const hermiteError = std::sqrt( sumHermite/count );
			double const gap = std::sqrt( sumGap/count );
			double const inverseGap = std::sqrt( sumInverse/count );

			gaps.push_back( gap );
			inverseGaps.push_back( inverseGap );

			if ( i == 0 )
				std::printf( "  %8s %6d %8d %13.4e %13.4e %13.4e %7s %13.4e "
				             "%7s\n", "", sizes[ i ], points, angleError,
				             hermiteError, gap, "-", inverseGap, "-" );
			else
				std::printf( "  %8s %6d %8d %13.4e %13.4e %13.4e %7.3f %13.4e "
				             "%7.3f\n", "", sizes[ i ], points, angleError,
				             hermiteError, gap,
				             rate( gaps[ i - 1 ], gap, 2.0 ), inverseGap,
				             rate( inverseGaps[ i - 1 ], inverseGap, 2.0 ) );

			// NEITHER ROUTE MAY BE FURTHER FROM THE OTHER THAN FROM THE TRUTH,
			// by more than a small factor. That is what says the gap is the
			// discretisation the two share rather than a difference of
			// definition -- a missing 2 pi R, a metric taken about the wrong
			// point, a gradient on the wrong side of the division. Those are
			// exactly the defects a single-route convergence table cannot see,
			// because they converge too.
			double const worst = std::max( angleError, hermiteError );
			BOOST_TEST( gap < 5.0*worst,
			            "k = " << order << ", n = " << sizes[ i ]
			            << ": the two extractions differ by " << gap << " in V' "
			            << "where the worse of them is " << worst << " from the "
			            << "exact field" );

			// AND THE HERMITE ROUTE MUST BE AN ANSWER IN ITS OWN RIGHT, not
			// merely consistent with the other one.
			BOOST_TEST( hermiteError < 20.0*angleError,
			            "k = " << order << ", n = " << sizes[ i ]
			            << ": the Hermite route is " << hermiteError
			            << " from the exact field against the ray fit's "
			            << angleError << ". It is an independent extraction and "
			            << "is expected to be comparable, not merely "
			            << "consistent" );
		}
		std::fflush( stdout );

		double const ratio = static_cast<double>( sizes.back() )/sizes.front();
		double const gapRate = rate( gaps.front(), gaps.back(), ratio );
		double const inverseRate = rate( inverseGaps.front(), inverseGaps.back(),
		                                 ratio );

		std::printf( "  k = %d: the disagreement converges at %.3f in V' and "
		             "%.3f in < R^-2 >\n", order, gapRate, inverseRate );

		BOOST_TEST( gapRate > order + 1.0 - sequenceSlack,
		            "k = " << order << ": the disagreement between the two "
		            << "extractions converges at " << gapRate << " rather than "
		            << "at the field's own " << order + 1.0 << ". A gap that "
		            << "does not converge at the discretisation's order is not "
		            << "the discretisation, and two routes that share only the "
		            << "field and the level have nothing else to disagree "
		            << "about" );

		BOOST_TEST( inverseRate > order + 1.0 - sequenceSlack,
		            "k = " << order << ": the disagreement in < R^-2 > converges "
		            << "at " << inverseRate );
	}

	// PARAMETRISATION INVARIANCE. Same contour, same field, same rule, rays
	// drawn from somewhere else.
	{
		SolvedEquilibrium solved( eq, box, 3, 32 );
		meq::CriticalPoint const axis = solved.axis();
		meq::ContourTracer tracer( solved.theSolver() );
		tracer.setTargetTurn( 0.05 );

		Fitted const home = fittedSurface(
			tracer, axis, levelAt( exact, 0.25 ), 512 );
		meq::SurfaceAverages const centred = meq::surfaceAverages( tracer,
		                                                           home.fit );

		std::printf( "\n  The same contour re-parametrised about a displaced "
		             "origin. An integral over a curve\n  does not know how the "
		             "curve was parametrised, so these must be zero.\n" );
		std::printf( "  %8s %14s %14s %10s\n", "offset", "dV'/V'", "d<R^-2>",
		             "transv" );

		for ( double offset : { 0.05, 0.10, 0.15 } )
		{
			meq::CriticalPoint displaced = axis;
			displaced.r += offset;
			displaced.z -= 0.7*offset;

			meq::AngleParametrisation const fit =
				tracer.fitByAngle( home.contour, displaced, 512 );
			meq::SurfaceAverages const moved = meq::surfaceAverages( tracer, fit );

			double const dv = relative( moved.vPrime, centred.vPrime );
			double const da = relative( moved.inverseRSquared(),
			                            centred.inverseRSquared() );

			std::printf( "  %8.2f %14.4e %14.4e %10.3f\n", offset, dv, da,
			             fit.transversality );

			BOOST_TEST( dv < 1.0e-8,
			            "moving the ray origin by " << offset << " moved V' by "
			            << dv << " relative. The integral is over a curve and "
			            << "cannot depend on the rays used to sample it; a "
			            << "dependence is the metric identity being wrong" );
			BOOST_TEST( da < 1.0e-8,
			            "moving the ray origin by " << offset
			            << " moved < R^-2 > by " << da << " relative" );
		}
		tracer.setTolerance( 1.0e-12 );
		std::fflush( stdout );
	}
}

/*
 * THE CONTRACT: THE BAND FLAG, THE DIFFERENCED CONTROL, AND THE NAMED WRAPPERS.
 *
 * Small, cheap, and none of it is decoration.
 *
 * THE BAND FLAG READS FALSE HERE AND THAT IS ASSERTED RATHER THAN ASSUMED. The
 * fitted path has no band -- Gamma_h IS Gamma -- so false is the CORRECT answer
 * and not an absent one, and an edit that starts reporting a band on this path
 * should fail loudly. INVERSION-PLAN.md section 4.3 requires that a surface
 * crossing the band be flagged rather than counted, so that IN-2 reports the two
 * populations separately instead of quoting one rate over both; CLAUDE.md
 * records the count-not-mask version of the same mistake in the .nc as half of a
 * real defect.
 *
 * THE DIFFERENCED CONTROL MUST REFUSE WHERE IT DOES NOT EXIST. The contour route
 * has no equispaced neighbour to difference against, and returning the pointwise
 * answer there would make the control agree with the answer for the worst
 * possible reason -- which is exactly how a control quietly stops being one.
 *
 * AND THE NAMED QUANTITIES MUST BE THIN WRAPPERS AND NOT SECOND
 * IMPLEMENTATIONS. That is the whole shape of the facility: MANTA-COUPLING.md
 * section 5's list of geometry slots "is negotiated with the transport physics
 * case, not fixed by MaNTA" and is being revised, so each named quantity is one
 * line over one primitive rather than a function of its own. If a wrapper ever
 * drifts from the integrand it is supposed to be, the convergence tables above
 * would carry the drift and say nothing about it.
 */
BOOST_AUTO_TEST_CASE( theFacilityCarriesItsFlagsAndItsWrappersAreThin )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = nstxBox();
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );
	double const level = levelAt( exact, 0.25 );

	SolvedEquilibrium solved( eq, box, 2, 24 );
	meq::CriticalPoint const axis = solved.axis();
	meq::ContourTracer tracer( solved.theSolver() );

	Fitted const fitted = fittedSurface( tracer, axis, level, 256 );

	meq::SurfaceAverages const angle = meq::surfaceAverages( tracer, fitted.fit );
	meq::SurfaceAverages const hermite = meq::surfaceAverages( tracer,
	                                                           fitted.contour, 8 );
	tracer.setTolerance( 1.0e-12 );

	std::printf( "\n  The contract, k = 2, n = 24, Psi_N = 0.25: %d angle nodes, "
	             "%d Hermite nodes\n", static_cast<int>( angle.nodes.size() ),
	             static_cast<int>( hermite.nodes.size() ) );

	// THE BAND. False, zero, and None -- on the fitted path all three are the
	// right answer and all three are asserted.
	std::printf( "  band: angle route extended %d over %d nodes, deepest %.1e, "
	             "extension %s\n", angle.extendedNodes,
	             static_cast<int>( angle.nodes.size() ), angle.deepestBandNode,
	             meq::bandExtensionName( angle.bandExtension ) );

	BOOST_TEST( !angle.extended,
	            "the fitted path has no band between Gamma_h and Gamma, so no "
	            "node of a surface inside it may be flagged as extended. A true "
	            "here means either the tracer is extending where it need not, or "
	            "the flag is being set by something other than the tracer" );
	BOOST_TEST( angle.extendedNodes == 0 );
	BOOST_TEST( angle.deepestBandNode == 0.0 );
	bool const noExtensionConfigured =
		angle.bandExtension == meq::BandExtension::None;
	BOOST_TEST( noExtensionConfigured,
	            "the fitted path configures no band extension, so a surface "
	            "built on it must not claim one as its provenance" );
	BOOST_TEST( !hermite.extended );
	BOOST_TEST( hermite.extendedNodes == 0 );

	for ( meq::SurfaceNode const &node : angle.nodes )
		BOOST_TEST( !node.extended,
		            "a node of a surface wholly inside Omega_h is flagged as "
		            "extended" );

	// THE DIFFERENCED CONTROL. Present where the trap lives, refused where it
	// does not.
	BOOST_TEST( angle.differencedAvailable(),
	            "the equispaced-angle route must carry the differenced metric "
	            "control: it is the trap of INVERSION-PLAN.md section 3.2 and "
	            "there is nothing to compare a spectral claim against without "
	            "it" );
	BOOST_TEST( !hermite.differencedAvailable(),
	            "the contour route has no equispaced neighbour to difference "
	            "against, so it must not claim a differenced control" );

	bool refused = false;
	try
	{
		hermite.averageDifferenced( []( meq::SurfaceNode const &node )
		{
			return 1.0/( node.r*node.r );
		} );
	}
	catch ( std::runtime_error const & )
	{
		refused = true;
	}
	BOOST_TEST( refused,
	            "averageDifferenced() returned a number for a surface that "
	            "carries no differenced control. Falling back on the pointwise "
	            "metric there would make the control agree with the answer for "
	            "the worst possible reason" );

	// THE WRAPPERS. Bit-identical to the facility calls they stand for -- these
	// are the same arithmetic in the same order, so anything but equality is a
	// wrapper that has drifted.
	double const inverse = angle.average( []( meq::SurfaceNode const &node )
	{
		return 1.0/( node.r*node.r );
	} );
	double const gradient = angle.average( []( meq::SurfaceNode const &node )
	{
		return node.gradient*node.gradient/( node.r*node.r );
	} );

	BOOST_TEST( angle.inverseRSquared() == inverse,
	            "inverseRSquared() is not average( 1/R^2 ). It is supposed to be "
	            "one line over the facility, not a second implementation" );
	BOOST_TEST( angle.gradPsiSquaredOverRSquared() == gradient,
	            "gradPsiSquaredOverRSquared() is not "
	            "average( |grad psi|^2 / R^2 )" );

	// < | grad psi |^2 / R^2 > IS < | q |^2 >, exactly, since grad psi = r q.
	// Checked because it is stated in the header and because a reader meeting
	// the long form is entitled to know the short one is the same number.
	double const fluxSquared = angle.average( []( meq::SurfaceNode const &node )
	{
		return node.qR*node.qR + node.qZ*node.qZ;
	} );
	BOOST_TEST( relative( gradient, fluxSquared ) < 1.0e-14,
	            "< |grad psi|^2 / R^2 > reads " << gradient << " against "
	            << "< |q|^2 > = " << fluxSquared << ". grad_bar( psi ) = r q, so "
	            << "these are the same quantity and the header says so" );

	// THE SAFETY FACTOR IS NEVER CALLED q. It is V' g < R^-2 > / 4 pi^2, RoPP
	// (142), with g( psi ) = R B_toroidal supplied by the caller because a
	// meq::Source does not carry one. Pinned against the product it stands for,
	// which fixes the 4 pi^2 -- an algebraic combination is exactly as easy to
	// get wrong by a constant as anything else, and nothing else here would see
	// it.
	double const g = 1.7;
	// 4 pi^2 is ( 2 pi )^2, written that way so the constant is visible.
	double const expected = angle.vPrime*g*angle.inverseRSquared()
	                        /( twoPi*twoPi );
	std::printf( "  safety factor at g = %.1f: %.12e\n", g,
	             angle.safetyFactor( g ) );
	BOOST_TEST( relative( angle.safetyFactor( g ), expected ) < 1.0e-15,
	            "safetyFactor( g ) is not V' g < R^-2 > / 4 pi^2" );

	// arcLength() IS ALGEBRAICALLY THE PRODUCING OBJECT'S OWN LENGTH, which is
	// not an independent measurement and is not quoted as one -- see the header.
	// What it checks is that every weight really is 2 pi R | dx/ds | ds /
	// | grad psi |: a missing r, a wrong dtheta or a gradient on the wrong side
	// of the division all break the identity.
	std::printf( "  arc length: facility %.14e, IN-1's fit %.14e, Hermite "
	             "%.14e\n", angle.arcLength(), fitted.fit.length(),
	             hermite.arcLength() );

	BOOST_TEST( relative( angle.arcLength(), fitted.fit.length() ) < 1.0e-13,
	            "the arc length through the facility is " << angle.arcLength()
	            << " and AngleParametrisation::length() is "
	            << fitted.fit.length() << ". They are algebraically the same "
	            << "number, so a difference is a weight that is not built the "
	            << "way SurfaceAverage.hpp says it is" );
	BOOST_TEST( relative( hermite.arcLength(),
	                      fitted.contour.hermiteLength( 10 ) ) < 1.0e-13,
	            "the arc length through the facility is " << hermite.arcLength()
	            << " and Contour::hermiteLength() is "
	            << fitted.contour.hermiteLength( 10 ) );

	std::fflush( stdout );
}
