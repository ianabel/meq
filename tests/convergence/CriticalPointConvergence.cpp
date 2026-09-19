#define BOOST_TEST_MODULE CriticalPointConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/GradShafranov.hpp"

// omp_get_max_threads() and omp_set_num_threads(), for
// theAxisSweepDoesNotDependOnTheThreadCount. Guarded on the same pair
// src/meq/Threading.hpp is: with either absent meq's sweep is serial by
// construction and the case has nothing to compare.
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
#include <omp.h>
#endif

#include "analytic/Soloviev.hpp"
#include "ConvergenceHarness.hpp"

/*
 * The acceptance test for INVERSION-PLAN.md stage IN-A: the magnetic axis as a
 * genuine critical point of psi_h, and the Poincare-Hopf audit over the whole
 * domain.
 *
 * TWO MEASUREMENTS, AND THEY ARE INDEPENDENT.
 *
 * The first is a rate. The located axis is where q_h vanishes; the true axis is
 * where q vanishes; so the distance between them is the pointwise error of q_h
 * divided by the local |dq/dx|, by the implicit function theorem applied to
 * q( x ) = 0. q converges at k+1 -- SolovievConvergence.cpp is that measurement
 * and this file reproduces its L2 column alongside, so that the two are read
 * together -- and the divisor is an O( 1 ) constant fixed by the equilibrium
 * rather than by the mesh. So THE AXIS POSITION IS EXPECTED AT k+1, and that is
 * what is asserted.
 *
 * It is worth being clear about what would happen if the axis were instead
 * found from the potential, which is the obvious alternative and is what a code
 * without a solved flux has to do. Newton on grad( psi_h ) = 0 has a residual
 * that is a DIFFERENTIATED L2 field, converging at k rather than k+1, and the
 * located point would then converge one order more slowly. That is the whole
 * content of "the residual is a solved field" in CriticalPoints.hpp, and it is
 * the reason this stage is cheap to make accurate.
 *
 * The second is a count, and it is not a rate at all. audit() walks the
 * boundary and returns the topological degree of q there, which by the
 * Poincare index theorem is the SUM of the indices of the zeros inside. It must
 * come out an integer, and the two things worth asserting are that it comes out
 * the RIGHT integer and that it goes on doing so where a spurious critical
 * point is most likely -- on a mesh far too coarse for the equilibrium. That is
 * the control, and it is the reason the under-resolved cases below are run at
 * all.
 *
 * THE SIGN OF THE SOLOV'EV AXIS, WHICH IS NOT WHAT THE PLAN ASSUMES.
 *
 * INVERSION-PLAN.md section 6 argues from the maximum principle that with
 * single-signed F >= 0 the interior critical point is a MAXIMUM. Every
 * Solov'ev fixture in tests/analytic has F = -( ( 1 - A ) r^2 + A ), which on
 * these domains is single-signed NEGATIVE -- so psi is a subsolution, its
 * maximum is on the boundary, and the axis is an interior MINIMUM. Measured
 * here, and it is a minimum on all four configurations.
 *
 * The plan's conclusion is untouched: a maximum and a minimum both have index
 * +1 in two dimensions, so #extrema - #saddles = 1 either way and one axis
 * still means no interior saddle. What it does invalidate is the seeding rule
 * "start from the largest nodal value of psi": that finds a corner of the
 * benchmark rectangle here. CriticalPointFinder::findAxis() seeds from both
 * nodal extremes for exactly this reason.
 *
 * WHEN chi IS A THEOREM, MEASURED RATHER THAN ASSUMED.
 *
 * Poincare-Hopf needs q TRANSVERSE to the boundary, which is a condition on
 * q . n and not on the boundary being a level set. An earlier draft of this
 * file asserted the opposite -- that a rectangle cut out of a larger
 * equilibrium could not be transverse -- and the audit contradicted it: on the
 * standard box q . n keeps one sign the whole way round with
 * min |q . n|/|q| = 0.15, so the hypothesis holds and winding == chi == 1 is a
 * theorem there rather than a coincidence. It is asserted on that basis below.
 *
 * The wide box of theWindingNumberIsASumOfIndicesAndNotACount is the other
 * case: it reaches past an X-point, q . n changes sign, transversality reads
 * 0.00 and the degree is 0 against chi = 1. That is not a contradiction and not
 * a defect -- the hypothesis is simply not satisfied -- and it is why
 * IndexAudit reports transversality instead of leaving the reader to guess.
 *
 * The degree equals the sum of the interior indices in BOTH cases, because that
 * half of the statement needs no transversality at all.
 *
 * AND THE SAME TWO MEASUREMENTS AGAIN AT THE X-POINT, WHICH IS XP-0.
 *
 * FREE-BOUNDARY-PLAN.md section 10.6's first stage is this file's axis study
 * with the saddle in place of the axis: the position of the X-point against a
 * closed form at a measured rate, and audit() reading +1 and -1 over a box
 * holding both critical points. It needs no free boundary and no normalisation,
 * only the same fixed-boundary Solov'ev solve the axis study runs. Three things
 * about it are worth knowing before reading the cases at the foot of this
 * file:
 *
 *   * the finder reaches the saddle with NO seed. findAxis() cannot -- it looks
 *     for an extremum -- and sweep() can, returning exactly one saddle on a box
 *     that holds exactly one, at every order and every mesh tried;
 *   * the reference has to be checked first and one of the four Solov'ev
 *     fixtures fails that check, which is why the closed-form case is separate
 *     and carries that fixture as its control;
 *   * a box holding an axis and an X-point CANNOT read audit.consistent(), and
 *     that is arithmetic rather than a discretisation failure -- see the case.
 */

namespace
{

	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	/// A critical point of the CLOSED FORM, to round-off: the magnetic axis when
	/// seeded near the axis, the X-point when seeded near the X-point.
	///
	/// Newton on grad( psi ) = 0 with the Hessian by central differences. The
	/// Hessian's accuracy does not reach the answer -- the fixed point of this
	/// iteration is where the analytic gradient vanishes, whatever steered it
	/// there -- which is the same observation CriticalPoints.hpp makes about its
	/// own Jacobian, made here in a place where it can be checked independently.
	struct ExactCritical
	{
		double r;
		double z;
		double psi;
		double gradient;
		double determinant;
		double trace;
	};

	ExactCritical exactCriticalPoint( Equilibrium const &eq,
	                                  double rGuess, double zGuess )
	{
		double r = rGuess;
		double z = zGuess;
		double const step = 1.0e-5;

		double hessian[ 2 ][ 2 ] = { { 0.0, 0.0 }, { 0.0, 0.0 } };

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );

			double a0 = 0.0;
			double a1 = 0.0;
			double b0 = 0.0;
			double b1 = 0.0;
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
			double const dr = -(  hessian[ 1 ][ 1 ]*gr - hessian[ 0 ][ 1 ]*gz )/det;
			double const dz = -( -hessian[ 1 ][ 0 ]*gr + hessian[ 0 ][ 0 ]*gz )/det;

			r += dr;
			z += dz;

			if ( std::abs( dr ) + std::abs( dz ) < 1.0e-15 )
				break;
		}

		ExactCritical axis;
		axis.r = r;
		axis.z = z;
		axis.psi = eq.psi( r, z );

		double gr = 0.0;
		double gz = 0.0;
		eq.gradPsi( r, z, gr, gz );
		axis.gradient = std::sqrt( gr*gr + gz*gz );
		axis.determinant = hessian[ 0 ][ 0 ]*hessian[ 1 ][ 1 ]
		                   - hessian[ 0 ][ 1 ]*hessian[ 1 ][ 0 ];
		axis.trace = hessian[ 0 ][ 0 ] + hessian[ 1 ][ 1 ];
		return axis;
	}

	/// One solve, kept alive.
	///
	/// The harness's measure() destroys its solver on the way out, which is
	/// exactly right for an error norm and no use at all here: the critical
	/// point finder borrows the flux and the potential and needs them to outlive
	/// the measurement. Member order is load bearing -- the mesh and the
	/// coefficients are referenced by the solver and must be constructed before
	/// it and destroyed after it -- and the class is non-copyable because the
	/// coefficients capture `this`.
	///
	/// The LINEAR path, as SolovievConvergence.cpp uses: the Solov'ev source does
	/// not depend on psi, so there is nothing for Newton to do and a
	/// FunctionCoefficient is the honest way to say so.
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
				  fluxCoeff( 2, [ this ]( mfem::Vector const &x, mfem::Vector &value )
				  {
					  eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
				  } ),
				  solver( mesh, orderIn )
			{
				solver.setSource( sourceCoeff );
				solver.setBoundaryData( psiCoeff );
				solver.solve();
				h = boxIn.width()/static_cast<double>( n );
			}

			SolvedEquilibrium( SolvedEquilibrium const & ) = delete;
			SolvedEquilibrium &operator=( SolvedEquilibrium const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			double errorFlux()
			{
				return solver.fluxError( fluxCoeff );
			}

			double errorPsi()
			{
				return solver.potentialError( psiCoeff );
			}

			double meshSize() const
			{
				return h;
			}

			mfem::Mesh &theMesh()
			{
				return mesh;
			}

		private:
			Equilibrium eq;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			mfem::VectorFunctionCoefficient fluxCoeff;
			meq::GradShafranovSolver solver;
			double h = 0.0;
	};

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}

	/// Slack on the rate asserted across the WHOLE mesh sequence -- see the
	/// comment above theMagneticAxisConvergesAtTheFluxesOwnOrder for why it is
	/// asserted there and not pair by pair. 0.25 rather than
	/// SolovievConvergence.cpp's 0.15 because a pointwise quantity is noisier
	/// than an integrated one: the located axis sits somewhere different inside
	/// its element at every refinement, so the constant in front of h^(k+1) is
	/// not the same constant at every level. Measured, the three sequences come
	/// out at 2.340, 3.484 and 4.447 against design orders of 2, 3 and 4 -- they
	/// clear k+1 itself, so the slack is not what is holding the assertion up.
	double const rateSlack = 0.25;

	struct AxisMeasurement
	{
		double h;
		int traceDofs;
		double distance;
		double errorFlux;
		double pointwiseFlux;
		double residual;
		double overshoot;
		meq::CriticalPointType type;
	};

	void printCriticalTable( char const *what, int order, ExactCritical const &exact,
	                         std::vector<AxisMeasurement> const &points )
	{
		std::printf( "\n  %s as a zero of q, Solov'ev NSTX, k = %d\n", what, order );
		std::printf( "  exact point ( %.12f, %.12f ), |grad psi| = %.2e, %s\n",
		             exact.r, exact.z, exact.gradient,
		             exact.determinant > 0.0
		               ? ( exact.trace > 0.0 ? "minimum" : "maximum" ) : "saddle" );
		std::printf( "  %8s %9s %14s %7s %14s %7s %14s %7s %7s %9s\n",
		             "h", "trace", "|x - x_ax|", "rate", "L2(q)", "rate",
		             "|q_h-q|(x_ax)", "rate", "ratio", "over" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			AxisMeasurement const &p = points[ i ];
			double const ratio = p.distance/p.pointwiseFlux;
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %14.6e %7s %14.6e %7s %14.6e %7s %7.2f %9.1e\n",
				             p.h, p.traceDofs, p.distance, "-", p.errorFlux, "-",
				             p.pointwiseFlux, "-", ratio, p.overshoot );
			}
			else
			{
				double const refinement = points[ i - 1 ].h/p.h;
				std::printf( "  %8.5f %9d %14.6e %7.3f %14.6e %7.3f %14.6e %7.3f %7.2f %9.1e\n",
				             p.h, p.traceDofs, p.distance,
				             rate( points[ i - 1 ].distance, p.distance, refinement ),
				             p.errorFlux,
				             rate( points[ i - 1 ].errorFlux, p.errorFlux, refinement ),
				             p.pointwiseFlux,
				             rate( points[ i - 1 ].pointwiseFlux, p.pointwiseFlux,
				                   refinement ),
				             ratio, p.overshoot );
			}
		}
		std::fflush( stdout );
	}

	void printAudit( char const *label, meq::IndexAudit const &audit,
	                 std::vector<meq::CriticalPoint> const &found )
	{
		int sum = 0;
		for ( std::size_t i = 0; i < found.size(); ++i )
			sum += found[ i ].index;

		std::printf( "  %-28s winding %+2d (defect %.1e)  chi %+2d  loops %d  "
		             "found %2d sum %+2d  worst turn %.3f  min|q| %.2e  "
		             "transverse %s (%.2f)\n",
		             label, audit.windingNumber, audit.windingDefect,
		             audit.eulerCharacteristic, audit.boundaryLoops,
		             static_cast<int>( found.size() ), sum, audit.worstTurn,
		             audit.smallestFlux, audit.transverse ? "yes" : "no",
		             audit.transversality );
		for ( std::size_t i = 0; i < found.size(); ++i )
			std::printf( "  %-28s   %-10s at ( %.6f, %.6f )  psi %+.6e  "
			             "index %+d  |q| %.2e  over %.1e\n",
			             "", meq::criticalPointName( found[ i ].type ),
			             found[ i ].r, found[ i ].z, found[ i ].psi,
			             found[ i ].index, found[ i ].fluxResidual,
			             found[ i ].overshoot );
		std::fflush( stdout );
	}

	/// | q_h - q | at a given physical point, evaluated inside whichever element
	/// holds it. This is the quantity the located axis inherits, and it is
	/// measured rather than inferred so that the two columns can be read against
	/// each other.
	double pointwiseFluxError( mfem::Mesh &mesh, mfem::GridFunction const &flux,
	                           Equilibrium const &eq, double r, double z )
	{
		mfem::DenseMatrix point( 2, 1 );
		point( 0, 0 ) = r;
		point( 1, 0 ) = z;

		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> ips;
		mesh.FindPoints( point, elements, ips );
		if ( elements[ 0 ] < 0 )
			return std::numeric_limits<double>::quiet_NaN();

		mfem::Vector value( 2 );
		flux.GetVectorValue( elements[ 0 ], ips[ 0 ], value );

		double exactR = 0.0;
		double exactZ = 0.0;
		eq.flux( r, z, exactR, exactZ );
		return std::sqrt( ( value( 0 ) - exactR )*( value( 0 ) - exactR )
		                  + ( value( 1 ) - exactZ )*( value( 1 ) - exactZ ) );
	}

	std::vector<int> const axisMeshes = { 4, 8, 16, 32 };

	/*
	 * XP-0'S REFERENCE, WHICH IS A GEOMETRIC RULE AND NOT A TRANSCRIBED DECIMAL.
	 *
	 * Cerfon & Freidberg close the twelve-coefficient up-down asymmetric form by
	 * putting the X-point at x_sep = 1 - 1.1 delta eps, y_sep = -1.1 kappa eps,
	 * and Soloviev.hpp's nstx() is the solve of those twelve conditions at
	 * eps = 0.78, kappa = 2, delta = 0.35. So the reference below is the rule
	 * that DEFINED the fixture rather than a number read off a printout of it,
	 * and theClosedFormXPointIsASaddleOfTheClosedForm is what checks the two
	 * agree.
	 */
	double const nstxEpsilon = 0.78;
	double const nstxKappa = 2.0;
	double const nstxDelta = 0.35;

	double const nstxXPointR = 1.0 - 1.1*nstxDelta*nstxEpsilon;   // 0.699700
	double const nstxXPointZ = -1.1*nstxKappa*nstxEpsilon;        // -1.716000

	/*
	 * THE X-POINT BOX, AND WHY IT IS NOT standardBox().
	 *
	 * The benchmark rectangle stops at z = -0.6 and the X-point is at
	 * z = -1.716, so the standard box does not reach it -- which is why the
	 * axis study above has never had a saddle in it. This one is the same 0.8
	 * across, so its h column is directly comparable with the axis table's, and
	 * it is centred on the X-point at ( 0.437, 0.480 ) of its own extent. The
	 * closed form has exactly TWO critical points over the whole region
	 * [ 0.3, 1.9 ] x [ -2.1, 1.9 ] -- measured -- so this box contains the
	 * saddle and nothing else, and "the saddle sweep() found" needs no
	 * tie-break to be well defined.
	 */
	inline Rectangle xPointBox()
	{
		return Rectangle{ 0.35, 1.15, -2.10, -1.30 };
	}

	/*
	 * AND THE BOX THAT HOLDS BOTH, for the audit. Tall rather than square,
	 * because the two critical points are 1.83 apart in z and 0.62 in r; the
	 * cells are 2:1 and the mesh family is still dyadic and shape regular, which
	 * is all a rate or a degree needs. theWindingNumberIsASumOfIndicesAndNotACount
	 * already uses a 2:1 box for the same reason.
	 */
	inline Rectangle divertedBox()
	{
		return Rectangle{ 0.35, 1.55, -2.10, 0.30 };
	}

	/*
	 * FIVE LEVELS RATHER THAN THE AXIS STUDY'S FOUR, AND THE EXTRA ONE IS THE
	 * INSTRUMENT RATHER THAN THE ANSWER.
	 *
	 * A pointwise error carries a constant that is wherever in its element the
	 * point happens to fall, so a rate taken across a short sequence measures
	 * the ratio of two of those constants as much as it measures an order.
	 * MEASURED, the whole-sequence rate of the X-point position:
	 *
	 *     levels                  k = 1    k = 2    k = 3
	 *     { 4, 8, 16, 32 }         1.685    2.796    4.500
	 *     { 4, 8, 16, 32, 64 }     1.814    3.566    4.223
	 *
	 * -- against targets of 2, 3 and 4 less rateSlack. Four levels puts k = 1
	 * UNDER the target and five puts it over; carried on to n = 128 it reads
	 * 1.916 from n = 4 and 2.017 from n = 8. The pointwise error of q at the
	 * same point does the same thing in the same places, 1.696 over four levels
	 * against 1.862 over five, so what the short sequence was short of was
	 * levels and not order.
	 *
	 * THE MARGIN AT k = 1 IS 0.064 AND IT DOES NOT GET BETTER BY ADDING A SIXTH
	 * LEVEL, which is worth knowing before anybody tries. n = 128 lifts k = 1 to
	 * 1.916 and drops k = 2 to 2.884, because the X-point at n = 64 happens to
	 * fall where q_h is unusually good and the next refinement has nothing left
	 * to gain there. There is no sequence that is comfortable at all three
	 * orders at once; that is what a pointwise quantity is like, and it is why
	 * the rate is read against the pointwise error of the FIELD as well as
	 * against k+1.
	 *
	 * Fifteen solves cost about six seconds, the coarse ones being free, so the
	 * longer lever arm is bought rather than argued for.
	 */
	std::vector<int> const xPointMeshes = { 4, 8, 16, 32, 64 };

	/*
	 * A SECOND X-POINT BOX, FOR THE SEEDED SEARCH, AND THE REASON IS A
	 * MEASUREMENT RATHER THAN A PREFERENCE.
	 *
	 * theSeededSaddleSearchCostsLessThanASweep asserts that the seeded search
	 * and sweep() return the same root TO ROUND-OFF, which is only a
	 * well-posed thing to ask where there IS one root to return. q_h is
	 * discontinuous across a face, so where the true X-point sits within the
	 * jump of that face BOTH neighbours hold a clean root of their own
	 * polynomial, O( h^(k+1) ) apart, and "the same root" stops being defined:
	 * sweep() merges the two and keeps whichever it saw first, while a seeded
	 * search stops at the first ring that yields one and keeps whichever its
	 * seed reached. Neither is wrong and they need not agree.
	 *
	 * xPointBox() PUTS THE X-POINT THERE AT ONE MESH IN THREE, WHICH IS AN
	 * ACCIDENT OF ITS CORNER. nstx()'s X-point is at r = 0.699700, and
	 * [ 0.35, 1.15 ] divided sixteen ways has a mesh line at r = 0.700000 --
	 * 3.0e-04 away, 0.6% of a cell. Measured at k = 1, n = 16 on that box:
	 * element 237 holds the root at r = 0.699826 and element 238 holds its own
	 * at r = 0.700027, both strictly inside, 6.9e-04 apart, and the two routes
	 * return different ones of them. At k = 2 and 3 the flux is accurate enough
	 * that both sides put the root on the same side of the line and the
	 * question does not arise, which is why this is a coarse-mesh effect rather
	 * than a defect.
	 *
	 * SO THE BOX IS CHOSEN TO KEEP THE X-POINT OFF THE MESH LINES -- vertical,
	 * horizontal AND the diagonal of the split -- at every mesh in the sweep.
	 * The worst margin over n = 8, 16, 32 is 0.136 of a cell here against
	 * xPointBox()'s 0.006, a factor of 23, and it is the best available over
	 * boxes of the same 0.8 width on a 0.005 grid of corners. IT IS ALSO
	 * CHECKED RATHER THAN TRUSTED: the case asserts the located saddle's own
	 * distance to the boundary of its element, so a fixture that drifts back
	 * into the degenerate case says so instead of failing the agreement
	 * assertion with a message about round-off.
	 *
	 * It holds exactly one critical point, for xPointBox()'s reason: the closed
	 * form has two over [ 0.3, 1.9 ] x [ -2.1, 1.9 ], this box lies inside that
	 * region, and the magnetic axis -- measured at ( 1.3182, 0.0111 ) -- is
	 * nowhere near z in [ -2.095, -1.295 ]. The case REQUIREs the sweep to
	 * return exactly one saddle rather than trusting that.
	 */
	inline Rectangle seededXPointBox()
	{
		return Rectangle{ 0.34, 1.14, -2.095, -1.295 };
	}

	/// How far a located point sits inside its OWN reference element, in
	/// reference units: min( x, y, 1 - x - y ) on a triangle, so 1/3 at the
	/// centroid and zero on a face.
	///
	/// It is the complement of CriticalPoint::overshoot, which measures the
	/// same thing on the other side of the face and reads zero for everything
	/// inside. A root at depth zero is one the neighbour may hold just as
	/// well, and that is the situation seededXPointBox() exists to avoid.
	double referenceDepth( meq::CriticalPoint const &point )
	{
		return std::min( std::min( point.referenceX, point.referenceY ),
		                 1.0 - point.referenceX - point.referenceY );
	}

}

/*
 * The reference the rate is measured against, checked before it is used.
 *
 * A convergence table against a wrong reference converges beautifully to the
 * wrong place -- the standing hazard this project's testing stance is organised
 * around, and the one that produced two wrong sets of Solov'ev coefficients. So
 * the closed-form axis is verified to be a critical point of the closed form
 * before anything is compared with it, and its Hessian is reported so that the
 * "it is a minimum, not a maximum" claim in the file comment above is a
 * measurement.
 */
BOOST_AUTO_TEST_CASE( theClosedFormAxisIsACriticalPointOfTheClosedForm )
{
	struct Case
	{
		char const *name;
		Equilibrium eq;
	};

	std::vector<Case> cases = {
		{ "nstx",            Equilibrium::nstx() },
		{ "nstxAsPublished", Equilibrium::nstxAsPublished() },
		{ "iterExample2",    Equilibrium::iterExample2() },
		{ "nstxExample3",    Equilibrium::nstxExample3() }
	};

	std::printf( "\n  The closed-form magnetic axis of each Solov'ev fixture\n" );
	for ( std::size_t i = 0; i < cases.size(); ++i )
	{
		ExactCritical const axis = exactCriticalPoint( cases[ i ].eq, 1.0, 0.0 );
		std::printf( "  %-18s ( %.12f, %.12f )  psi %+.6e  |grad psi| %.2e  "
		             "det %+.4f  tr %+.4f  %s\n",
		             cases[ i ].name, axis.r, axis.z, axis.psi, axis.gradient,
		             axis.determinant, axis.trace,
		             axis.determinant > 0.0
		               ? ( axis.trace > 0.0 ? "minimum" : "maximum" ) : "saddle" );

		BOOST_TEST( axis.gradient < 1.0e-12,
		            cases[ i ].name << ": the reference axis is not a critical point, "
		            << "|grad psi| = " << axis.gradient );

		// Positive determinant AND positive trace is a positive definite
		// Hessian, which is a minimum. Both are asserted, because det > 0 alone
		// admits a maximum and this file's account of the sign turns on which it
		// is.
		BOOST_TEST( axis.determinant > 0.0,
		            cases[ i ].name << ": the axis is a saddle, det = "
		            << axis.determinant );
		BOOST_TEST( axis.trace > 0.0,
		            cases[ i ].name << ": the axis is a maximum, not the minimum "
		            << "this file's sign argument says it is; trace = " << axis.trace );
	}
	std::fflush( stdout );
}

/*
 * IN-A'S ACCEPTANCE RATE.
 *
 * The axis is a zero of q_h, so its distance from the true axis is the
 * POINTWISE error of q_h there, divided by dq/dx. That is a different quantity
 * from the L2 error of q, and the difference is what shapes this test.
 *
 * BOTH ARE MEASURED AND PRINTED SIDE BY SIDE, which is what turns the argument
 * into a measurement: the table carries the L2 error of q, the pointwise error
 * of q at the exact axis, and the position error, and the last two agree to
 * within a factor of four at every point of every sequence. The factor is not
 * free either -- dq/dx at the axis has eigenvalues 1.30 and 0.306, so by the
 * implicit function theorem the position error is between 0.77 and 3.27 times
 * the pointwise flux error, and the measured range over all twelve points is
 * 0.77 to 3.37. So the root finder adds nothing to the error of the field it is
 * rooting, and the residual it roots really is the solved one.
 *
 * WHY THE RATE IS ASSERTED OVER THE WHOLE SEQUENCE RATHER THAN PAIR BY PAIR.
 *
 * A pointwise error is not an L2 error and does not fall smoothly. The measured
 * pointwise flux error at the axis, k = 1, over h = 0.2 to 0.025, converges at
 * 2.20, 3.47 and 1.18 on successive pairs -- averaging 2 and visiting neither
 * neighbour of it -- because the axis sits somewhere different inside its
 * element at every refinement and the DG error oscillates within an element.
 * The position error inherits exactly that raggedness: 4.17, 1.54, 1.31. At
 * k = 3 the same columns are cleaner but still not smooth, 3.87 / 4.06 / 4.44
 * and 3.79 / 5.94 / 3.61.
 *
 * So the assertion is the pattern ExtensionConvergence.cpp already uses for a
 * quantity whose per-pair rate is not a rate: monotone decrease at every
 * refinement, and the rate ACROSS THE WHOLE SEQUENCE at k+1. Measured, that
 * comes out 2.34, 3.48 and 4.45 for k = 1, 2, 3 against targets of 2, 3 and 4.
 *
 * A per-pair assertion loose enough to admit 1.31 at k = 1 would be 0.7 of
 * slack, which is not an assertion about anything.
 */
BOOST_AUTO_TEST_CASE( theMagneticAxisConvergesAtTheFluxesOwnOrder )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();
	ExactCritical const exact = exactCriticalPoint( eq, 1.0, 0.0 );

	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<AxisMeasurement> points;

		for ( std::size_t m = 0; m < axisMeshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, axisMeshes[ m ] );

			meq::CriticalPointFinder finder( run.theSolver() );
			meq::CriticalPoint const axis = finder.findAxis();

			AxisMeasurement point;
			point.h = run.meshSize();
			point.traceDofs = run.theSolver().numTraceDofs();
			point.distance = std::sqrt( ( axis.r - exact.r )*( axis.r - exact.r )
			                            + ( axis.z - exact.z )*( axis.z - exact.z ) );
			point.errorFlux = run.errorFlux();
			point.pointwiseFlux = pointwiseFluxError( run.theMesh(),
			                                          run.theSolver().flux(),
			                                          eq, exact.r, exact.z );
			point.residual = axis.fluxResidual;
			point.overshoot = axis.overshoot;
			point.type = axis.type;
			points.push_back( point );
		}

		printCriticalTable( "Magnetic axis", order, exact, points );

		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			// The residual is what says the returned point is a root at all. A
			// point reported at |q| = 1e-3 is a place Newton gave up, and every
			// rate below it would be meaningless.
			BOOST_TEST( points[ i ].residual < 1.0e-10,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the located axis has |q_h| = " << points[ i ].residual
			            << ", which is not a root" );

			// The sign argument in the file comment, asserted where it is used:
			// F < 0 on this box, so the axis is an interior MINIMUM. If this ever
			// reads maximum, the flux has been handed over with the wrong sign --
			// which the winding number cannot see, since index( -v ) = index( v )
			// in two dimensions.
			bool const isMinimum
				= ( points[ i ].type == meq::CriticalPointType::Minimum );
			BOOST_TEST( isMinimum,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the Solov'ev axis came out a "
			            << meq::criticalPointName( points[ i ].type )
			            << " rather than a minimum" );

			// THE STRUCTURAL ASSERTION, and the one that says what IN-A is for:
			// the located point is no worse than the pointwise error of the field
			// it was rooted in, up to the conditioning of dq/dx. 10 against a
			// linearised bound of 3.27 and a measured worst of 3.37.
			double const ratio = points[ i ].distance/points[ i ].pointwiseFlux;
			BOOST_TEST( ratio < 10.0,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the located axis is " << ratio
			            << " times further from the true axis than the pointwise "
			            << "error of q_h there. The root finder is supposed to add "
			            << "nothing to the error of the field it roots" );
		}

		// Monotone decrease. A pointwise quantity does not converge smoothly, but
		// it must converge.
		for ( std::size_t i = 1; i < points.size(); ++i )
			BOOST_TEST( points[ i ].distance < points[ i - 1 ].distance,
			            "k = " << order << ": refining from h = " << points[ i - 1 ].h
			            << " to " << points[ i ].h << " moved the axis error from "
			            << points[ i - 1 ].distance << " to " << points[ i ].distance );

		double const refinement = points.front().h/points.back().h;
		double const measured = rate( points.front().distance,
		                              points.back().distance, refinement );
		double const expected = order + 1.0 - rateSlack;

		std::printf( "  over the whole sequence: the axis position converges at "
		             "%.3f, wanted %.2f\n", measured, expected );
		std::fflush( stdout );

		BOOST_TEST( measured >= expected,
		            "k = " << order
		            << ": the axis position converged at " << measured
		            << " across the whole sequence, wanted " << expected
		            << ". The axis is a zero of q, so it inherits q's order; a rate "
		            << "near " << order
		            << " would say the residual being rooted is a differentiated "
		            << "field rather than a solved one" );
	}
}

/*
 * The audit, on every Solov'ev configuration in tests/analytic that has its axis
 * inside the benchmark box.
 *
 * frcExample1 is excluded and the reason is geometric rather than numerical:
 * its axis sits at r = 1.4071, just outside the box's rMax = 1.4, so there is no
 * critical point to find and the winding number would be zero -- correctly.
 */
BOOST_AUTO_TEST_CASE( theWindingNumberAuditReturnsOneOnTheAnalyticFixtures )
{
	struct Case
	{
		char const *name;
		Equilibrium eq;
	};

	std::vector<Case> cases = {
		{ "nstx",            Equilibrium::nstx() },
		{ "nstxAsPublished", Equilibrium::nstxAsPublished() },
		{ "iterExample2",    Equilibrium::iterExample2() },
		{ "nstxExample3",    Equilibrium::nstxExample3() }
	};

	Rectangle const box = meq::tests::standardBox();

	std::printf( "\n  Poincare-Hopf audit, k = 2, n = 16\n" );

	for ( std::size_t i = 0; i < cases.size(); ++i )
	{
		SolvedEquilibrium run( cases[ i ].eq, box, 2, 16 );

		meq::CriticalPointFinder finder( run.theSolver() );
		meq::IndexAudit const audit = finder.audit();
		std::vector<meq::CriticalPoint> const found = finder.sweep();

		printAudit( cases[ i ].name, audit, found );

		BOOST_TEST( audit.windingNumber == 1,
		            cases[ i ].name << ": the degree of q on the boundary is "
		            << audit.windingNumber << ", wanted 1" );

		// A degree is an integer. This is not a discretisation error that shrinks
		// with h -- it is either at round-off or the walk failed to resolve a
		// rotation, in which case worstTurn is at pi and the number above is
		// arbitrary.
		BOOST_TEST( audit.windingDefect < 1.0e-9,
		            cases[ i ].name << ": the accumulated turning is "
		            << audit.turning << " of a full turn, which is not an integer; "
		            << "worst single turn " << audit.worstTurn );

		BOOST_TEST( audit.worstTurn < 1.5,
		            cases[ i ].name << ": a single sample turned by "
		            << audit.worstTurn << " radians, so the boundary walk is "
		            << "under-sampled and the winding number is not to be believed" );

		BOOST_TEST( audit.eulerCharacteristic == 1,
		            cases[ i ].name << ": chi of a triangulated rectangle came out "
		            << audit.eulerCharacteristic << ", not 1" );

		BOOST_TEST( audit.boundaryLoops == 1,
		            cases[ i ].name << ": the boundary threaded into "
		            << audit.boundaryLoops << " loops, not 1" );

		// q points consistently outward across this boundary -- measured, see the
		// file comment -- so Poincare-Hopf applies and the degree is not merely
		// equal to chi, it is required to be. Both halves are asserted, because
		// asserting the conclusion without the hypothesis would make it a
		// coincidence dressed as a theorem.
		BOOST_TEST( audit.transverse,
		            cases[ i ].name << ": q is not transverse to the boundary, "
		            << "min |q.n|/|q| = " << audit.transversality
		            << ". The comparison against chi below is then not entitled" );

		BOOST_TEST( audit.consistent(),
		            cases[ i ].name << ": the degree is " << audit.windingNumber
		            << " and chi is " << audit.eulerCharacteristic
		            << " on a boundary q IS transverse to, where Poincare-Hopf says "
		            << "they must agree" );

		// The sweep is not exhaustive, so it cannot prove the winding number
		// right. What it can do is fail to agree with it, which IS evidence -- of
		// a missed root, a misclassified one, or a duplicate.
		int sum = 0;
		for ( std::size_t j = 0; j < found.size(); ++j )
			sum += found[ j ].index;
		BOOST_TEST( sum == audit.windingNumber,
		            cases[ i ].name << ": the sweep found " << found.size()
		            << " critical points whose indices sum to " << sum
		            << ", against a boundary degree of " << audit.windingNumber
		            << ". The sweep is not exhaustive, so this disagreement means a "
		            << "root was missed or misclassified, not that the degree is wrong" );
	}
}

/*
 * THE CONTROL, and it is what makes the test above mean anything.
 *
 * Spurious critical points come from noise in q_h, so they are likeliest where
 * q_h is worst -- on a mesh far too coarse for the equilibrium. n = 2 puts four
 * triangles across a domain whose flux varies by a factor of several, at k = 1,
 * which is as under-resolved as this benchmark gets while still solving.
 *
 * THE WINDING NUMBER IS 1 AT EVERY ONE OF THEM, which is IN-A's acceptance
 * criterion for this control. If it ever comes out otherwise that is a FINDING
 * and not a test to relax: it would say the discrete field has acquired a
 * critical point the continuous one does not have, and the right response is to
 * record the mesh it happened on.
 *
 * AND ONE THING THAT DID COME OUT, WHICH IS WORTH RECORDING RATHER THAN TIDYING
 * AWAY. At h = 0.4, k = 1 the boundary degree reads 1 and the sweep finds NO
 * zero of q_h anywhere in the mesh. That is not a failure of the root finder
 * and it is not a contradiction: q_h is DISCONTINUOUS, and the Poincare index
 * theorem is about continuous fields. Each element's polynomial puts its own
 * zero a little way into a neighbour's territory, the neighbour's polynomial
 * does the same in the other direction, and at a jump of O( h^(k+1) ) against
 * an element of size h there is a window in which the zero belongs to neither.
 * Measured: with CriticalPointFinder's containment allowance raised to 0.2 of a
 * reference element the sweep does find a candidate, 0.102 outside its own
 * element and 0.076 from the true axis -- an answer as coarse as the mesh that
 * produced it. The default allowance of 0.10 refuses it, and refusing is right.
 *
 * The window closes as the mesh refines, which is why every finer mesh here
 * finds the axis and why this is a property of h rather than a defect. The
 * assertion below pins it to n = 2 exactly: a finer mesh losing the axis fails.
 *
 * Note what the audit can and cannot see. A spurious maximum and a spurious
 * saddle sum to zero, so they leave the winding number at 1 and pass -- which is
 * why sweep()'s count is printed beside it.
 */
BOOST_AUTO_TEST_CASE( theAuditSurvivesADeliberatelyUnderResolvedMesh )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();

	std::vector<int> const coarse = { 2, 3, 4, 6 };

	std::printf( "\n  The under-resolved control, k = 1\n" );

	for ( std::size_t i = 0; i < coarse.size(); ++i )
	{
		SolvedEquilibrium run( eq, box, 1, coarse[ i ] );

		meq::CriticalPointFinder finder( run.theSolver() );
		meq::IndexAudit const audit = finder.audit();
		std::vector<meq::CriticalPoint> const found = finder.sweep();

		char label[ 64 ];
		std::snprintf( label, sizeof( label ), "k = 1, n = %2d, h = %.3f",
		               coarse[ i ], run.meshSize() );
		printAudit( label, audit, found );

		BOOST_TEST( audit.windingNumber == 1,
		            label << ": the degree of q on the boundary is "
		            << audit.windingNumber
		            << " on an under-resolved mesh. That is a finding rather than a "
		            << "test to relax -- record the mesh it happened on" );

		BOOST_TEST( audit.windingDefect < 1.0e-9,
		            label << ": the turning is " << audit.turning
		            << " of a full turn, which is not an integer" );

		// Whatever the sweep reports must be a genuine root, whether or not it
		// found all of them.
		for ( std::size_t j = 0; j < found.size(); ++j )
			BOOST_TEST( found[ j ].fluxResidual < 1.0e-10,
			            label << ": a reported critical point has |q_h| = "
			            << found[ j ].fluxResidual << ", which is not a root" );

		if ( found.empty() )
		{
			// The recorded finding above, pinned to the one mesh it happens on. A
			// finer mesh reaching this branch is a regression; n = 2 reaching the
			// other branch is an improvement and passes.
			BOOST_TEST( coarse[ i ] == 2,
			            label << ": the sweep found no zero of q_h at all. That is "
			            << "expected only at h = 0.4, where the jump in q_h is wide "
			            << "enough that the zero belongs to no element -- see the "
			            << "comment above this test" );
			continue;
		}

		int sum = 0;
		for ( std::size_t j = 0; j < found.size(); ++j )
			sum += found[ j ].index;
		BOOST_TEST( sum == audit.windingNumber,
		            label << ": the sweep's indices sum to " << sum
		            << " against a boundary degree of " << audit.windingNumber );
	}
}

/*
 * A DEGREE IS A SUM OF INDICES AND NEVER A COUNT, demonstrated rather than
 * asserted.
 *
 * INVERSION-PLAN.md section 5 says this in capitals because it was got wrong
 * once during the survey that produced the plan, and a warning in a header is
 * worth less than a live case. iterExample2 has its axis at ( 1.051, 0.024 )
 * and an X-point at ( 0.88384, -0.704 ); a box containing both contains an
 * extremum of index +1 -- a minimum here, since this fixture's F is negative --
 * and a saddle of index -1, and the degree of q on its boundary is ZERO.
 *
 * So a caller reading "degree 0" as "no critical points here" would be wrong by
 * two critical points, one of which is an X-point. The audit is a
 * certification, never an exclusion test.
 */
BOOST_AUTO_TEST_CASE( theWindingNumberIsASumOfIndicesAndNotACount )
{
	Equilibrium const eq = Equilibrium::iterExample2();

	// Deliberately not the standard box: this one reaches down past the X-point
	// at z = -0.704, which the standard box's zMin = -0.6 stops short of.
	Rectangle const box = { 0.7, 1.4, -0.9, 0.6 };

	SolvedEquilibrium run( eq, box, 2, 24 );

	meq::CriticalPointFinder finder( run.theSolver() );
	meq::IndexAudit const audit = finder.audit();
	std::vector<meq::CriticalPoint> const found = finder.sweep();

	std::printf( "\n  Degree is a sum, not a count: iterExample2 on "
	             "[%.1f,%.1f]x[%.1f,%.1f], k = 2, n = 24\n",
	             box.rMin, box.rMax, box.zMin, box.zMax );
	printAudit( "axis and X-point", audit, found );

	BOOST_TEST( audit.windingNumber == 0,
	            "the degree of q on a boundary enclosing one extremum and one "
	            << "saddle is " << audit.windingNumber << ", wanted 0" );

	BOOST_TEST( audit.eulerCharacteristic == 1,
	            "chi of a triangulated rectangle came out "
	            << audit.eulerCharacteristic << ", not 1" );

	// And this is the point: the degree does not agree with chi here, and the
	// domain is emphatically not empty.
	BOOST_TEST( !audit.consistent(),
	            "the degree and chi agree on a domain containing a saddle, which "
	            << "would make this test's demonstration empty" );

	// Why that is not a contradiction: Poincare-Hopf's hypothesis fails. q . n
	// changes sign on this boundary, because reaching past the X-point means
	// reaching past the flux surface that turns round there. Asserting this is
	// what stops the row above being read as a defect.
	BOOST_TEST( !audit.transverse,
	            "q is transverse to this boundary, with min |q.n|/|q| = "
	            << audit.transversality
	            << ". If that is so then Poincare-Hopf applies and a degree of "
	            << audit.windingNumber << " against chi = "
	            << audit.eulerCharacteristic << " is a genuine contradiction "
	            << "rather than an inapplicable hypothesis" );

	int extrema = 0;
	int saddles = 0;
	for ( std::size_t i = 0; i < found.size(); ++i )
	{
		if ( found[ i ].type == meq::CriticalPointType::Saddle )
			++saddles;
		if ( found[ i ].type == meq::CriticalPointType::Maximum
		     || found[ i ].type == meq::CriticalPointType::Minimum )
			++extrema;
	}

	BOOST_TEST( extrema == 1,
	            "the sweep found " << extrema << " extrema, wanted 1" );
	BOOST_TEST( saddles == 1,
	            "the sweep found " << saddles << " saddles, wanted 1. Without one "
	            << "the zero degree above is not the cancellation this test is "
	            << "about" );

	// The title claim, spelled out: +1 and -1 sum to the degree, which is the
	// half of Poincare-Hopf that needs no transversality and holds here.
	int sum = 0;
	for ( std::size_t i = 0; i < found.size(); ++i )
		sum += found[ i ].index;
	BOOST_TEST( sum == audit.windingNumber,
	            "the two critical points have indices summing to " << sum
	            << " against a boundary degree of " << audit.windingNumber );

	// The X-point of iterExample2 is published: Soloviev.hpp records it at
	// exactly ( x_sep, y_sep ) = ( 0.88384, -0.704 ), which is where the twelve
	// constraints put it by construction. Checking the located saddle against
	// that is a second, independent statement that this is the X-point and not
	// some artefact of an under-resolved corner.
	for ( std::size_t i = 0; i < found.size(); ++i )
	{
		if ( found[ i ].type != meq::CriticalPointType::Saddle )
			continue;
		double const dr = found[ i ].r - 0.88384;
		double const dz = found[ i ].z + 0.704;
		double const distance = std::sqrt( dr*dr + dz*dz );
		std::printf( "  located saddle is %.3e from the published X-point\n",
		             distance );
		BOOST_TEST( distance < 1.0e-3,
		            "the located saddle is " << distance
		            << " from the published X-point at ( 0.88384, -0.704 )" );
	}

	// findAxis() must refuse here rather than return the saddle or guess between
	// two candidates -- it looks for an interior EXTREMUM, and there is exactly
	// one, so it succeeds and returns that one.
	meq::CriticalPoint axis;
	BOOST_TEST( finder.tryFindAxis( axis ),
	            "findAxis() failed on a domain with exactly one interior extremum" );
	bool const axisIsMinimum = ( axis.type == meq::CriticalPointType::Minimum );
	BOOST_TEST( axisIsMinimum,
	            "findAxis() returned a " << meq::criticalPointName( axis.type )
	            << " where the sweep found a minimum" );
}

/*
 * IN-A'S AXIS IS NOT GradShafranovSolver::psiAxis(), AND THIS MEASURES THE GAP.
 *
 * The solver's psi_ax is the extreme NODAL value of psi_h, chosen because the
 * bordered Newton of the normalised-profile path needs a constraint it can
 * differentiate -- one nodal value is one entry of the discrete unknown, so the
 * border row is sparse. IN-A's axis is the critical point. INVERSION-PLAN.md
 * lists conflating them as risk 7 and CriticalPoints.hpp says at length not to,
 * so here is the number that makes it concrete.
 *
 * WHAT IS MEASURED, and the second column is the sharp one:
 *
 *   position   the distance from the critical point to the extreme nodal point.
 *              That is the distance to the nearest Gauss-Lobatto node, which is
 *              O( h ) at every k and does NOT converge at the discretisation's
 *              rate. Measured, it is between 0.10 h and 0.38 h across the
 *              sequence at k = 2: first order in h with a constant that is
 *              wherever in its element the extremum happens to fall.
 *
 *   value      psi_h at the critical point against the extreme nodal value.
 *              psi_h is smooth and quadratic about its own extremum, so a node
 *              O( h ) away is O( h^2 ) off in value -- second order at EVERY k,
 *              where psi_h's own L2 error is k+1.
 *
 * AND THE ASSERTION IS ON THE RATIO, NOT ON A RATE, because neither column
 * above converges smoothly: both inherit the wandering constant in the position
 * gap, and squaring it in the value column makes that worse. The measured
 * per-pair rates of the value gap at k = 2 are 0.06, 1.54 and 0.78 -- averaging
 * about 2 and visiting neither neighbour of it. What IS stable is that the gap
 * is far larger than psi_h's own L2 error and grows relative to it under
 * refinement, which is exactly the statement "these are different quantities
 * and refinement separates them rather than reconciling them":
 *
 *       k        coarsest        finest
 *       2            2.05         202
 *       3          293           4204
 *
 * -- the value gap as a multiple of L2( psi_h ) on the same mesh.
 *
 * Note the value gap is NOT the error in psi_ax as an approximation to the true
 * axis flux, and it is not asserted to be. It is the disagreement between two
 * ways of reading the same discrete solution.
 */
BOOST_AUTO_TEST_CASE( theCriticalPointIsNotTheExtremeNodalValue )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();

	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<double> hs;
		std::vector<double> positionGaps;
		std::vector<double> valueGaps;
		std::vector<double> psiErrors;

		for ( std::size_t m = 0; m < axisMeshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, axisMeshes[ m ] );

			meq::CriticalPointFinder finder( run.theSolver() );
			meq::CriticalPoint const axis = finder.findAxis();

			// The extreme nodal value, and where its node is. Read the same way
			// GradShafranovSolver's bordered Newton reads it -- over the dofs of
			// the potential space -- and located by asking the element for its
			// nodal points, which is what makes "the distance to the nearest node"
			// a thing that can be measured at all.
			mfem::GridFunction const &psiH = run.theSolver().potential();
			mfem::FiniteElementSpace const &space = *psiH.FESpace();

			double best = 0.0;
			double bestR = 0.0;
			double bestZ = 0.0;
			bool haveBest = false;

			mfem::Array<int> dofs;
			for ( int element = 0; element < run.theMesh().GetNE(); ++element )
			{
				space.GetElementDofs( element, dofs );
				mfem::FiniteElement const *fe = space.GetFE( element );
				mfem::IntegrationRule const &nodes = fe->GetNodes();
				if ( nodes.Size() != dofs.Size() )
					continue;

				mfem::ElementTransformation *trans
					= run.theMesh().GetElementTransformation( element );

				for ( int i = 0; i < dofs.Size(); ++i )
				{
					double const value = psiH( dofs[ i ] );
					// The Solov'ev axis is a MINIMUM on this box, so the extreme
					// nodal value in the sense that matters is the smallest one.
					// The high-beta source, whose F has the other sign, wants the
					// largest -- see the file comment.
					if ( haveBest && value >= best )
						continue;

					mfem::Vector physical( 2 );
					trans->Transform( nodes[ i ], physical );

					best = value;
					bestR = physical( 0 );
					bestZ = physical( 1 );
					haveBest = true;
				}
			}

			BOOST_TEST_REQUIRE( haveBest );

			hs.push_back( run.meshSize() );
			positionGaps.push_back( std::sqrt( ( axis.r - bestR )*( axis.r - bestR )
			                                   + ( axis.z - bestZ )*( axis.z - bestZ ) ) );
			valueGaps.push_back( std::abs( axis.psi - best ) );
			psiErrors.push_back( run.errorPsi() );
		}

		std::printf( "\n  The critical point against the extreme nodal value, k = %d\n",
		             order );
		std::printf( "  %8s %14s %7s %14s %7s %14s %7s\n",
		             "h", "|dx| position", "rate", "|dpsi| value", "rate",
		             "L2(psi_h)", "rate" );
		for ( std::size_t i = 0; i < hs.size(); ++i )
		{
			if ( i == 0 )
			{
				std::printf( "  %8.5f %14.6e %7s %14.6e %7s %14.6e %7s\n",
				             hs[ i ], positionGaps[ i ], "-", valueGaps[ i ], "-",
				             psiErrors[ i ], "-" );
			}
			else
			{
				double const ratio = hs[ i - 1 ]/hs[ i ];
				std::printf( "  %8.5f %14.6e %7.3f %14.6e %7.3f %14.6e %7.3f\n",
				             hs[ i ],
				             positionGaps[ i ],
				             rate( positionGaps[ i - 1 ], positionGaps[ i ], ratio ),
				             valueGaps[ i ],
				             rate( valueGaps[ i - 1 ], valueGaps[ i ], ratio ),
				             psiErrors[ i ],
				             rate( psiErrors[ i - 1 ], psiErrors[ i ], ratio ) );
			}
		}
		std::fflush( stdout );

		// The two are different quantities, so the gap between them must not be
		// at round-off. Scaled by h so that "different" means "different by an
		// amount the mesh sets", not by an absolute number that would have to be
		// retuned at every order.
		for ( std::size_t i = 0; i < hs.size(); ++i )
			BOOST_TEST( positionGaps[ i ] > 1.0e-3*hs[ i ],
			            "k = " << order << ", h = " << hs[ i ]
			            << ": the critical point and the extreme nodal point are "
			            << positionGaps[ i ] << " apart, which is round-off. They "
			            << "are supposed to be different quantities" );

		// And at k >= 2 the sharp statement, made on the ratio rather than on a
		// rate for the reason set out above this test: the disagreement between
		// the two readings dwarfs psi_h's own L2 error, and refinement makes it
		// dwarf it by more. No amount of refinement makes them interchangeable.
		if ( order >= 2 )
		{
			double const coarsestRatio = valueGaps.front()/psiErrors.front();
			double const finestRatio = valueGaps.back()/psiErrors.back();

			std::printf( "  the nodal gap is %.4g times L2( psi_h ) on the coarsest "
			             "mesh and %.4g on the finest\n",
			             coarsestRatio, finestRatio );
			std::fflush( stdout );

			BOOST_TEST( finestRatio > 10.0,
			            "k = " << order << ": on the finest mesh the disagreement "
			            << "between the critical point and the extreme nodal value "
			            << "is only " << finestRatio
			            << " times psi_h's own L2 error. They are supposed to be "
			            << "different quantities by far more than the solution error" );

			BOOST_TEST( finestRatio > 10.0*coarsestRatio,
			            "k = " << order << ": the disagreement went from "
			            << coarsestRatio << " to " << finestRatio
			            << " times L2( psi_h ) under a factor of eight in h. It is "
			            << "supposed to grow, because the gap is second order in h "
			            << "at every k while psi_h converges at k+1" );
		}
	}
}

/*
 * THE REFERENCE COORDINATES ARE THE ONES THE ROOT WAS FOUND AT, AND THE TEST
 * THAT SAYS SO PUSHES THEM BACK THROUGH THE ELEMENT MAP.
 *
 * CriticalPoint carries ( referenceX, referenceY ) so that a caller can call
 * CalcShape() on the axis element at the axis -- FREE-BOUNDARY-PLAN.md section
 * 11.5's option 3, where those shape functions ARE a row of a bordered
 * Jacobian. The alternative is to hand back ( r, z ) alone and let the caller
 * invert the map with TransformBack, and CLAUDE.md records why that is a trap:
 * a CLAMPED inverse returns a point on the element boundary rather than
 * failing, and a field that is constant over the element cannot tell the
 * difference.
 *
 * SO WHAT IS ASSERTED IS THE ROUND TRIP, which is the only thing that can catch
 * the plumbing being wrong: transform( referencePoint() ) must be ( r, z ). It
 * is required at round-off rather than at a tolerance, because both sides are
 * the SAME transformation applied to the SAME point -- rootInElement() sets
 * found.r from exactly this call. A tolerance here would pass on a stale or a
 * defaulted reference point whenever the element happened to be small.
 *
 * AND psi AT THAT POINT IS CHECKED TOO, because a caller of option 3 evaluates
 * BOTH: the shape functions give the Jacobian row and their contraction with
 * the element's dofs gives the residual. If those two disagreed the border
 * would be differentiating a different quantity from the one it constrains.
 */
BOOST_AUTO_TEST_CASE( theReferenceCoordinatesReproduceTheLocatedPoint )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();

	std::printf( "\n  THE REFERENCE POINT, PUSHED BACK THROUGH THE ELEMENT MAP\n" );
	std::printf( "    %5s %5s %8s %11s %11s %13s %13s\n",
	             "k", "n", "element", "xi", "eta", "round trip", "psi round trip" );

	for ( int order = 1; order <= 3; ++order )
	{
		for ( std::size_t m = 0; m < axisMeshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, axisMeshes[ m ] );

			meq::CriticalPointFinder finder( run.theSolver() );
			meq::CriticalPoint const axis = finder.findAxis();

			BOOST_TEST_REQUIRE( axis.element >= 0,
				"the axis was located but carries no element, so there is nothing "
				"to evaluate a shape function on" );

			mfem::IntegrationPoint const ip = axis.referencePoint();

			// A LOCAL transformation, never Mesh::GetElementTransformation( int ):
			// that hands out shared scratch and "calling this function resets
			// pointers obtained from previous calls". CriticalPoints.cpp uses a
			// thread_local for the same reason and CLAUDE.md records six call
			// sites that had to be repaired for it.
			mfem::IsoparametricTransformation transformation;
			run.theMesh().GetElementTransformation( axis.element, &transformation );

			mfem::Vector physical( 2 );
			transformation.Transform( ip, physical );

			double const dr = physical( 0 ) - axis.r;
			double const dz = physical( 1 ) - axis.z;
			double const roundTrip = std::sqrt( dr*dr + dz*dz );

			// psi from the reference point, against the psi the finder reported.
			double const psiThere =
				run.theSolver().potential().GetValue( axis.element, ip );
			double const psiGap = std::abs( psiThere - axis.psi );

			std::printf( "    %5d %5d %8d %11.4e %11.4e %13.3e %13.3e\n",
			             order, axisMeshes[ m ], axis.element,
			             axis.referenceX, axis.referenceY, roundTrip, psiGap );
			std::fflush( stdout );

			// ROUND-OFF, NOT A TOLERANCE. Both sides are the same map applied to
			// the same reference point, so anything above the last few bits of
			// the coordinates means the reference point is not the one the root
			// was found at. Scaled by the position so that this reads the same on
			// a domain of any size.
			double const scale = std::max( 1.0, std::abs( axis.r )
			                                    + std::abs( axis.z ) );
			BOOST_TEST( roundTrip <= 1.0e-13*scale,
				"the reference point ( " << axis.referenceX << ", "
				<< axis.referenceY << " ) of element " << axis.element
				<< " transforms to ( " << physical( 0 ) << ", " << physical( 1 )
				<< " ) where the located axis is ( " << axis.r << ", " << axis.z
				<< " ). These are the same map applied to the same point, so a gap "
				"means the reference coordinates are stale, defaulted, or from "
				"another element." );

			BOOST_TEST( psiGap <= 1.0e-12*std::max( 1.0, std::abs( axis.psi ) ),
				"psi_h at the reference point is " << psiThere
				<< " where the finder reported " << axis.psi
				<< ". A caller of option 3 evaluates both the shape functions and "
				"their contraction with the element dofs at this point, so these "
				"disagreeing means the border would differentiate one quantity and "
				"constrain another." );
		}
	}
}

/*
 * THE WARM-START ENTRY POINT: SAME ANSWER, A FRACTION OF THE WORK.
 *
 * tryFindAxisFrom() exists so that a bordered Newton can re-locate the axis
 * once per Jacobian without paying for a sweep of every element. That is only
 * worth having if it gives the SAME point, so the first assertion is agreement
 * with findAxis() to round-off -- not to a tolerance, because both routes root
 * the same element's polynomial with the same Newton and the same stopping
 * rule, so where they agree at all they agree to the last bits. A tolerance
 * would hide the case where the seeded search finds a NEIGHBOURING element's
 * version of the same root, which is a real thing q_h's discontinuity permits
 * and is exactly what the containment machinery exists to arbitrate.
 *
 * SEEDED FROM A DISTANCE, because that is the case it is for. In a continuation
 * the seed is the previous iterate's axis, which is a step away rather than on
 * top of the answer, so seeding exactly at the answer would test nothing about
 * the widening. Half an element is a step much larger than a converging Newton
 * takes near the end.
 *
 * AND THE COST IS MEASURED IN FIELD EVALUATIONS RATHER THAN IN SECONDS. This
 * project's standing rule is that a timing on this machine is a measurement
 * about the machine; the honest cost here is how many ELEMENTS are searched,
 * which is a property of the algorithm and reproduces anywhere. sweep() roots
 * every element; this roots the seed and two rings.
 */
BOOST_AUTO_TEST_CASE( theSeededSearchFindsTheSameAxisForLessWork )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();

	std::printf( "\n  THE SEEDED SEARCH AGAINST findAxis()\n" );
	std::printf( "    %5s %5s %9s %13s %13s %11s\n",
	             "k", "n", "elements", "seeded from", "gap to sweep", "reached" );

	for ( int order = 1; order <= 3; ++order )
	{
		for ( std::size_t m = 1; m < axisMeshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, axisMeshes[ m ] );

			meq::CriticalPointFinder finder( run.theSolver() );
			meq::CriticalPoint const reference = finder.findAxis();

			// The sense the fixture actually has, rather than a guess. Every
			// Solov'ev fixture here has F single-signed NEGATIVE, so psi is a
			// subsolution and the axis is an interior MINIMUM -- the file comment
			// above measures that. Passing the wrong one would make this case
			// return false everywhere and look like a defect in the search.
			meq::AxisSense const sense =
				reference.type == meq::CriticalPointType::Maximum
					? meq::AxisSense::Maximum : meq::AxisSense::Minimum;

			// HALF AN ELEMENT AWAY, DIAGONALLY, so that the seed is neither the
			// answer nor aligned with the mesh.
			double const h = run.meshSize();
			double const seedR = reference.r + 0.5*h;
			double const seedZ = reference.z - 0.5*h;

			meq::CriticalPoint seeded;
			bool const found =
				finder.tryFindAxisFrom( seedR, seedZ, sense, seeded );

			double const dr = seeded.r - reference.r;
			double const dz = seeded.z - reference.z;
			double const gap = found ? std::sqrt( dr*dr + dz*dz )
			                         : std::numeric_limits<double>::infinity();

			std::printf( "    %5d %5d %9d %13.3e %13.3e %11s\n",
			             order, axisMeshes[ m ], run.theMesh().GetNE(),
			             0.5*h*std::sqrt( 2.0 ), gap,
			             found ? "yes" : "NO" );
			std::fflush( stdout );

			BOOST_TEST_REQUIRE( found,
				"the seeded search found no " << ( sense == meq::AxisSense::Maximum
				                                   ? "maximum" : "minimum" )
				<< " half an element from the axis at k = " << order << ", n = "
				<< axisMeshes[ m ] << ", where findAxis() found one at ( "
				<< reference.r << ", " << reference.z << " )." );

			// ROUND-OFF. Both routes root the same polynomial with the same
			// Newton; they either reach the same root or they reach a different
			// element's version of it, and the second is what this is here to
			// catch.
			BOOST_TEST( gap <= 1.0e-12*std::max( 1.0, std::abs( reference.r ) ),
				"the seeded search returned ( " << seeded.r << ", " << seeded.z
				<< " ) where findAxis() returns ( " << reference.r << ", "
				<< reference.z << " ) -- " << gap << " apart. Seeded from half an "
				"element away, these should be the same root to the last bits; a "
				"gap of order h^(k+1) means it found a NEIGHBOURING element's "
				"version of it, which setContainment() is what arbitrates." );

			BOOST_TEST( seeded.element == reference.element,
				"the seeded search credits the axis to element " << seeded.element
				<< " where findAxis() credits it to " << reference.element
				<< ". The position agreed, so this is the containment tie-break "
				"landing differently -- harmless for a position and NOT harmless "
				"for a caller building a Jacobian row on that element's dofs." );
		}
	}
}

/*
 * IT RETURNS false RATHER THAN LYING, WHICH IS THE HALF THAT MAKES THE OTHER
 * HALF USABLE.
 *
 * tryFindAxisFrom() takes the extremum NEAREST its seed where the seed region
 * offers several, which is right for a caller following one critical point and
 * would be badly wrong as a way to DISCOVER one. So the contract that matters
 * is what it does when there is nothing to find: it must decline, not return
 * whatever the widening happened to reach.
 *
 * TWO WAYS OF HAVING NOTHING TO FIND, and they are different failures:
 *
 *   * seeded in a CORNER of the domain, far from the axis, where the Solov'ev
 *     psi is monotone and q does not vanish at all within reach;
 *   * seeded ON the axis but asking for the WRONG SENSE. Every fixture here has
 *     an interior MINIMUM, so asking for a Maximum must fail even though a
 *     critical point is sitting under the seed. That is the sense trap of the
 *     file header made into an assertion: it is also what a caller who handed
 *     the raw flux block -- which holds -q, turning every Minimum into a
 *     Maximum -- would see, so a search that answered anyway would make that
 *     mistake invisible.
 *
 * The corner case is deliberately NOT asserted to fail for a specific reason:
 * with seedRings widened far enough it would eventually reach the axis, and the
 * assertion is about the DEFAULT radius. Raising the rings is a documented way
 * to search further, so the case pins the shipped default rather than a
 * property of the algorithm.
 */
BOOST_AUTO_TEST_CASE( theSeededSearchDeclinesRatherThanGuessing )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = meq::tests::standardBox();
	int const order = 2;
	int const n = 16;

	SolvedEquilibrium run( eq, box, order, n );
	meq::CriticalPointFinder finder( run.theSolver() );
	meq::CriticalPoint const reference = finder.findAxis();

	meq::AxisSense const right = reference.type == meq::CriticalPointType::Maximum
		? meq::AxisSense::Maximum : meq::AxisSense::Minimum;
	meq::AxisSense const wrong = right == meq::AxisSense::Maximum
		? meq::AxisSense::Minimum : meq::AxisSense::Maximum;

	std::printf( "\n  WHEN THERE IS NOTHING TO FIND ( k = %d, n = %d )\n",
	             order, n );
	std::printf( "    the axis is a %s at ( %.4f, %.4f )\n",
	             meq::criticalPointName( reference.type ),
	             reference.r, reference.z );

	// (1) FAR AWAY. The bottom-left corner of the box, which on standardBox()
	// against nstx() is most of the domain away from the axis.
	meq::CriticalPoint far;
	bool const foundFar = finder.tryFindAxisFrom( box.rMin, box.zMin, right, far );
	std::printf( "    seeded at the ( %.2f, %.2f ) corner, right sense: %s\n",
	             box.rMin, box.zMin, foundFar ? "FOUND SOMETHING" : "declined" );

	BOOST_TEST( !foundFar,
		"seeded in the corner at ( " << box.rMin << ", " << box.zMin
		<< " ) the search returned a " << meq::criticalPointName( far.type )
		<< " at ( " << far.r << ", " << far.z << " ), " << far.element
		<< ". Two rings of face neighbours do not reach the axis from there, so "
		"whatever this is it is not the axis -- and a caller following a branch "
		"would have been handed it silently. If the search radius has been "
		"widened deliberately, this case pins the SHIPPED default and is what "
		"has to move." );

	// (2) ON the axis, wrong sense. There IS a critical point under the seed.
	meq::CriticalPoint mistyped;
	bool const foundWrong =
		finder.tryFindAxisFrom( reference.r, reference.z, wrong, mistyped );
	std::printf( "    seeded ON the axis, asking for a %s: %s\n",
	             wrong == meq::AxisSense::Maximum ? "maximum" : "minimum",
	             foundWrong ? "FOUND SOMETHING" : "declined" );
	std::fflush( stdout );

	BOOST_TEST( !foundWrong,
		"seeded exactly on a " << meq::criticalPointName( reference.type )
		<< " and asked for the opposite sense, the search returned a "
		<< meq::criticalPointName( mistyped.type ) << " at ( " << mistyped.r
		<< ", " << mistyped.z << " ). The sense filter is what stands between a "
		"caller and the raw flux block, which holds -q and turns every Minimum "
		"into a Maximum with every winding number unchanged." );

	// AND THE CONTROL: the same seed with the RIGHT sense must succeed, or the
	// two refusals above are not evidence about the sense filter at all -- they
	// would be evidence that the search never finds anything.
	meq::CriticalPoint proper;
	bool const foundRight =
		finder.tryFindAxisFrom( reference.r, reference.z, right, proper );
	BOOST_TEST( foundRight,
		"the control failed: seeded on the axis with the CORRECT sense the "
		"search found nothing, so the two refusals above say nothing about the "
		"sense filter." );
}

/*
 * XP-0, FREE-BOUNDARY-PLAN.md section 10.6: THE X-POINT AGAINST A CLOSED FORM.
 *
 * Everything above XP-1 in that plan assumes MEQ can locate an X-point at the
 * order its flux converges at, and these three cases are the measurement that
 * says it can. They need no free boundary, no normalisation and no coupling: a
 * Solov'ev source, an exact Dirichlet datum, and a saddle whose position is
 * fixed by the twelve constraints that built the fixture.
 *
 * THE FIRST THING TO ESTABLISH IS THAT THE REFERENCE IS A REFERENCE, and
 * Soloviev.hpp says in as many words that for one of its four fixtures it is
 * not: nstxAsPublished()'s prescribed X-point is not where its psi has a saddle
 * at all. A rate study against that point would converge, at a clean rate, to a
 * place the finder is right not to be -- the standing hazard of this file. So
 * the closed-form check is a separate case with its own control, run before
 * anything is compared against it.
 *
 * WHY nstx() AND NOT ONE OF THE OTHER TWO. All three of nstx(), iterExample2()
 * and nstxExample3() pass the closed-form check below, so any of them would
 * serve. nstx() is the one section 10.2 names, it is the fixture the axis study
 * above already runs, and its X-point is the furthest of the three from its own
 * axis -- 1.83 against nstxExample3's 1.60 and iterExample2's 0.75 -- so a box
 * that holds the saddle and nothing else is roomy rather than contrived.
 */
BOOST_AUTO_TEST_CASE( theClosedFormXPointIsASaddleOfTheClosedForm )
{
	struct Case
	{
		char const *name;
		Equilibrium eq;
		double r;
		double z;

		/// Whether the point above is where this fixture's psi actually has its
		/// saddle. False for exactly one of them, and that one is the control.
		bool prescribedIsTheSaddle;
	};

	std::vector<Case> cases = {
		{ "nstx",            Equilibrium::nstx(),
		  nstxXPointR, nstxXPointZ, true },
		{ "iterExample2",    Equilibrium::iterExample2(),
		  0.88384,     -0.704,      true },
		{ "nstxExample3",    Equilibrium::nstxExample3(),
		  0.71257,     -1.4586,     true },
		{ "nstxAsPublished", Equilibrium::nstxAsPublished(),
		  nstxXPointR, nstxXPointZ, false }
	};

	std::printf( "\n  THE PRESCRIBED X-POINT AGAINST THE SADDLE OF THE "
	             "CLOSED FORM\n" );
	std::printf( "  %-17s %-24s %11s   %-24s %11s %10s %9s\n",
	             "", "prescribed", "|grad psi|", "saddle of psi", "|grad psi|",
	             "det", "moved by" );

	for ( std::size_t i = 0; i < cases.size(); ++i )
	{
		double gr = 0.0;
		double gz = 0.0;
		cases[ i ].eq.gradPsi( cases[ i ].r, cases[ i ].z, gr, gz );
		double const gradientThere = std::sqrt( gr*gr + gz*gz );

		ExactCritical const saddle =
			exactCriticalPoint( cases[ i ].eq, cases[ i ].r, cases[ i ].z );

		double const dr = saddle.r - cases[ i ].r;
		double const dz = saddle.z - cases[ i ].z;
		double const moved = std::sqrt( dr*dr + dz*dz );

		std::printf( "  %-17s ( %.6f, %.6f ) %11.2e   ( %.6f, %.6f ) %11.2e %10.4f %9.2e\n",
		             cases[ i ].name, cases[ i ].r, cases[ i ].z, gradientThere,
		             saddle.r, saddle.z, saddle.gradient, saddle.determinant,
		             moved );
		std::fflush( stdout );

		// Whatever Newton converged to is a critical point of the closed form,
		// on every fixture. That much is about the iteration, not about the
		// fixture, and it has to hold before either branch below means anything.
		BOOST_TEST( saddle.gradient < 1.0e-12,
		            cases[ i ].name << ": the reference saddle is not a critical "
		            << "point, |grad psi| = " << saddle.gradient );

		// A NEGATIVE DETERMINANT IS WHAT MAKES IT A SADDLE rather than an
		// extremum, and it is the half that distinguishes an X-point from a
		// second magnetic axis. Note the Hessian's determinant and dq/dx's
		// differ by r^2 > 0, so the sign is the same quantity either way.
		BOOST_TEST( saddle.determinant < 0.0,
		            cases[ i ].name << ": the point near the prescribed X-point "
		            << "is an extremum, det = " << saddle.determinant
		            << ". An X-point is an indefinite Hessian" );

		if ( cases[ i ].prescribedIsTheSaddle )
		{
			// The twelve constraints put psi_x = psi_y = 0 AT the prescribed
			// point, so this is round-off and not a tolerance.
			BOOST_TEST( gradientThere < 1.0e-12,
			            cases[ i ].name << ": grad psi at the prescribed X-point ( "
			            << cases[ i ].r << ", " << cases[ i ].z << " ) is "
			            << gradientThere << ", so the twelve constraints that "
			            << "built this fixture are not satisfied at it" );

			BOOST_TEST( moved < 1.0e-12,
			            cases[ i ].name << ": Newton walked " << moved
			            << " away from the prescribed X-point to find the saddle. "
			            << "The prescribed point IS the saddle for this fixture, "
			            << "so anything above round-off means the coefficients "
			            << "and the geometry have come apart" );
		}
		else
		{
			// THE CONTROL. Soloviev.hpp records that the printed coefficients
			// satisfy none of the twelve conditions and that the true saddle
			// sits at ( 0.6958, -1.8069 ) with psi = -8.7e-3. Asserting that
			// here is what makes "nstx() is the right fixture" a measurement
			// rather than a preference: if this ever passed the two checks
			// above, the two sets would be interchangeable and the choice would
			// carry no information.
			BOOST_TEST( gradientThere > 1.0e-3,
			            cases[ i ].name << ": grad psi at the prescribed X-point "
			            << "is " << gradientThere << ", i.e. this fixture DOES "
			            << "satisfy the X-point conditions after all. It is here "
			            << "as the counter-example, so that is a finding about "
			            << "Soloviev.hpp rather than about this test" );

			BOOST_TEST( moved > 1.0e-2,
			            cases[ i ].name << ": its saddle is only " << moved
			            << " from the prescribed X-point, where Soloviev.hpp "
			            << "records 9.1e-2. The control has stopped being a "
			            << "control" );
		}
	}
}

/*
 * XP-0'S ACCEPTANCE RATE: THE LOCATED X-POINT AGAINST THE EXACT ONE.
 *
 * Structurally this is the axis study one rung along, and deliberately so --
 * the acceptance in FREE-BOUNDARY-PLAN.md section 10.6 is "the rate findAxis()
 * reaches for the axis", so the measurement has to be the same measurement. The
 * saddle is a zero of q_h exactly as the axis is, its distance from the true
 * saddle is the POINTWISE error of q_h there divided by dq/dx, and both columns
 * are printed side by side so that claim can be read rather than believed.
 *
 * NOTHING IS SEEDED, WHICH IS THE FIRST THING WORTH RECORDING ABOUT THE FINDER.
 * findAxis() is no use here -- it looks for an interior EXTREMUM and refuses a
 * saddle by construction -- and tryFindAxisFrom() takes an AxisSense for the
 * same reason. What reaches the X-point is sweep(), which roots every element
 * from its own centre and its own quietest flux node and merges what it finds.
 * On this box it returns exactly ONE saddle at every k and every n, so "the
 * saddle sweep() found" is well defined with no prior, no tie-break and no
 * hand-placed guess. That is a stronger statement than the acceptance asked
 * for: it wanted a located point, and what it gets is a located point that no
 * knowledge of the answer went into.
 *
 * The table's first row is the proof of that, incidentally: at k = 1, n = 4 the
 * located saddle is 3.5e-3 from the reference. If any part of the answer had
 * come from the reference it would not be.
 *
 * THE CONDITIONING WINDOW IS COMPUTED, NOT QUOTED, AND IT IS TWO SIDED.
 *
 * q_h( x_h ) = 0 and q( x* ) = 0 give x_h - x* = -J^-1 ( q_h - q )( x* ) to
 * first order, with J = dq/dx = Hess( psi )/r at the saddle -- SYMMETRIC, since
 * a Hessian is, so its singular values are |lambda_1| and |lambda_2| and the
 * ratio of the position error to the pointwise flux error is trapped between
 * 1/|lambda|_max and 1/|lambda|_min. For nstx()'s X-point those are 0.899102
 * and 0.578734, so the window is [ 1.112, 1.728 ] -- and the measured ratio
 * over all fifteen points lands in it at both ends, 1.11 and 1.73. The window
 * is derived below from the fixture's own Hessian rather than written down, so
 * a change of fixture moves it correctly.
 *
 * That is a sharper statement than the axis study can make, where the measured
 * range 0.77 to 3.37 sits inside a window of 0.77 to 3.27 without touching it,
 * and it is the whole of what "the root finder adds nothing to the error of the
 * field it roots" means.
 *
 * THE RATE IS ASSERTED THE AXIS STUDY'S WAY, FOR THE AXIS STUDY'S REASON:
 * monotone decrease at every refinement, and the rate ACROSS THE WHOLE SEQUENCE
 * at k+1 less rateSlack. The per-pair rates are not rates -- at k = 2 they read
 * 3.19, 3.72, 1.48, 5.88 over five refinements, and carried to n = 128 the next
 * one is 0.16, because at n = 64 the X-point happens to fall where q_h is
 * unusually good and there is nothing left to gain. Any per-pair assertion loose
 * enough to admit 0.16 would be an assertion about nothing.
 */
BOOST_AUTO_TEST_CASE( theXPointConvergesAtTheFluxesOwnOrder )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = xPointBox();
	ExactCritical const exact =
		exactCriticalPoint( eq, nstxXPointR, nstxXPointZ );

	BOOST_TEST_REQUIRE( exact.determinant < 0.0,
		"the reference for this rate study is not a saddle, so there is nothing "
		"to measure an X-point finder against" );

	// The conditioning of the root, from the closed form's own Hessian.
	// dq/dx = Hess( psi )/r at a critical point, exactly, because the term that
	// differentiates the 1/r carries grad psi and grad psi vanishes here.
	double const jacobianDet = exact.determinant/( exact.r*exact.r );
	double const jacobianTrace = exact.trace/exact.r;
	double const discriminant =
		std::sqrt( jacobianTrace*jacobianTrace - 4.0*jacobianDet );
	double const lambdaPlus = 0.5*( jacobianTrace + discriminant );
	double const lambdaMinus = 0.5*( jacobianTrace - discriminant );
	double const biggest = std::max( std::abs( lambdaPlus ),
	                                 std::abs( lambdaMinus ) );
	double const smallest = std::min( std::abs( lambdaPlus ),
	                                  std::abs( lambdaMinus ) );

	// 15% either side of the linearised window, which is what the neglected
	// quadratic term and the difference between J and its discrete counterpart
	// are worth. Measured, the fifteen points fill 1.11 to 1.73 against a
	// window of 1.112 to 1.728, so the headroom is 14% below and 15% above.
	double const ratioFloor = 1.0/( 1.15*biggest );
	double const ratioCeiling = 1.15/smallest;

	std::printf( "\n  X-POINT: dq/dx has eigenvalues %+.6f and %+.6f, so the "
	             "position error\n  must be between %.3f and %.3f times the "
	             "pointwise flux error\n",
	             lambdaPlus, lambdaMinus, 1.0/biggest, 1.0/smallest );

	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<AxisMeasurement> points;

		for ( std::size_t m = 0; m < xPointMeshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, xPointMeshes[ m ] );

			meq::CriticalPointFinder finder( run.theSolver() );
			std::vector<meq::CriticalPoint> const found = finder.sweep();

			int saddles = 0;
			meq::CriticalPoint saddle;
			for ( std::size_t j = 0; j < found.size(); ++j )
			{
				if ( found[ j ].type != meq::CriticalPointType::Saddle )
					continue;
				++saddles;
				saddle = found[ j ];
			}

			// UNSEEDED, SO THE COUNT IS THE CONTRACT. This box holds exactly one
			// saddle of the closed form, so a sweep reporting two has invented
			// one and a sweep reporting none has lost it -- and in either case
			// "the located X-point" below would be a choice rather than an
			// answer. Required rather than merely checked, because the rest of
			// the loop would then be measuring the wrong point.
			BOOST_TEST_REQUIRE( saddles == 1,
				"k = " << order << ", n = " << xPointMeshes[ m ]
				<< ": the sweep found " << saddles << " saddles on a box holding "
				"exactly one. Nothing is seeded here, so which one is the X-point "
				"is not a question this test is entitled to answer." );

			AxisMeasurement point;
			point.h = run.meshSize();
			point.traceDofs = run.theSolver().numTraceDofs();
			double const dr = saddle.r - exact.r;
			double const dz = saddle.z - exact.z;
			point.distance = std::sqrt( dr*dr + dz*dz );
			point.errorFlux = run.errorFlux();
			point.pointwiseFlux = pointwiseFluxError( run.theMesh(),
			                                          run.theSolver().flux(),
			                                          eq, exact.r, exact.z );
			point.residual = saddle.fluxResidual;
			point.overshoot = saddle.overshoot;
			point.type = saddle.type;
			points.push_back( point );
		}

		printCriticalTable( "X-point", order, exact, points );

		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			BOOST_TEST( points[ i ].residual < 1.0e-10,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the located X-point has |q_h| = "
			            << points[ i ].residual << ", which is not a root" );

			double const ratio = points[ i ].distance/points[ i ].pointwiseFlux;

			BOOST_TEST( ratio < ratioCeiling,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the located X-point is " << ratio
			            << " times further from the true saddle than the pointwise "
			            << "error of q_h there, against a linearised ceiling of "
			            << 1.0/smallest
			            << ". The root finder is supposed to add nothing to the "
			            << "error of the field it roots" );

			// AND THE FLOOR, which the axis study does not assert and which is
			// available here because dq/dx is better conditioned at this saddle
			// than at that axis. It cannot be beaten: J is symmetric, so no
			// direction of flux error maps to a position error smaller than
			// 1/|lambda|_max times it. A ratio below this is not a better answer,
			// it is a measurement that has stopped meaning what it says --
			// typically a located point that is not the root of the field being
			// differenced against.
			BOOST_TEST( ratio > ratioFloor,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": the position error is only " << ratio
			            << " times the pointwise flux error, below the linearised "
			            << "floor of " << 1.0/biggest
			            << ". A symmetric dq/dx cannot do that" );
		}

		for ( std::size_t i = 1; i < points.size(); ++i )
			BOOST_TEST( points[ i ].distance < points[ i - 1 ].distance,
			            "k = " << order << ": refining from h = "
			            << points[ i - 1 ].h << " to " << points[ i ].h
			            << " moved the X-point error from "
			            << points[ i - 1 ].distance << " to "
			            << points[ i ].distance );

		double const refinement = points.front().h/points.back().h;
		double const measured = rate( points.front().distance,
		                              points.back().distance, refinement );
		double const alsoTheField = rate( points.front().pointwiseFlux,
		                                  points.back().pointwiseFlux,
		                                  refinement );
		double const expected = order + 1.0 - rateSlack;

		std::printf( "  over the whole sequence: the X-point position converges "
		             "at %.3f and the pointwise\n  flux error at the exact "
		             "X-point at %.3f, wanted %.2f\n",
		             measured, alsoTheField, expected );
		std::fflush( stdout );

		// AND THE ROBUST HALF OF THE SAME STATEMENT, which is the one that
		// separates "the finder is wrong" from "the field is coarse here": the
		// position rate must be the pointwise flux error's own rate, whatever
		// that happens to be. Measured, the two differ by 0.048, 0.022 and 0.057
		// at k = 1, 2, 3 -- so this holds at a third of the slack the rate
		// itself is asserted at, and it holds even where the rate is thin.
		BOOST_TEST( std::abs( measured - alsoTheField ) < 0.3,
		            "k = " << order << ": the X-point position converges at "
		            << measured << " and the pointwise error of q at that same "
		            << "point at " << alsoTheField
		            << ". The located point is a root of q_h, so those are the "
		            << "same number up to the conditioning of dq/dx; a gap means "
		            << "the finder has stopped tracking the field it roots" );

		BOOST_TEST( measured >= expected,
		            "k = " << order
		            << ": the X-point position converged at " << measured
		            << " across the whole sequence, wanted " << expected
		            << ". The X-point is a zero of q exactly as the axis is, so "
		            << "it inherits q's order; a rate near " << order
		            << " would say the residual being rooted is a differentiated "
		            << "field rather than a solved one. The pointwise error of q "
		            << "at the same point reads " << alsoTheField
		            << ", which is where to look first -- if the two agree the "
		            << "finder is tracking its field and the shortfall is the "
		            << "field's" );
	}
}

/*
 * THE AUDIT OVER A BOX HOLDING BOTH: +1 AT THE AXIS, -1 AT THE X-POINT.
 *
 * This is the second half of XP-0's acceptance and it is the check a diverted
 * configuration will lean on, because it is the only one that does not need to
 * know where anything is. theWindingNumberIsASumOfIndicesAndNotACount already
 * demonstrates the arithmetic on iterExample2; what this adds is the fixture
 * XP-0 is measured on, at three orders, with each of the two points identified
 * against its own closed form rather than merely counted.
 *
 * THE DEGREE IS 0 AND chi IS 1, AND THAT IS FORCED RATHER THAN OBSERVED.
 *
 * A reader coming to this from the acceptance criterion may expect
 * audit.consistent(). It cannot hold, and the argument is one line:
 * Poincare-Hopf says that if q points outward everywhere on the boundary then
 * the sum of the interior indices is chi = 1; the interior indices here are +1
 * and -1 and sum to 0; therefore q does NOT point outward everywhere, the
 * hypothesis fails, and the degree is entitled to disagree with chi. There is
 * no box containing exactly one axis and one X-point on which the audit reads
 * consistent, whatever the mesh -- so asserting consistent() here would be
 * asserting something no discretisation could deliver.
 *
 * What IS asserted is the half of the theorem that needs no hypothesis at all:
 * the degree equals the sum of the indices of the zeros inside. That holds
 * here, at every order, and it is the half a diverted free-boundary solve would
 * use -- "the sweep has found everything the boundary says is in there".
 */
BOOST_AUTO_TEST_CASE( theAuditReadsPlusOneAtTheAxisAndMinusOneAtTheXPoint )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = divertedBox();

	ExactCritical const exactSaddle =
		exactCriticalPoint( eq, nstxXPointR, nstxXPointZ );
	ExactCritical const exactAxisPoint = exactCriticalPoint( eq, 1.0, 0.0 );

	std::printf( "\n  Poincare-Hopf over a box holding the axis AND the X-point, "
	             "[%.2f,%.2f]x[%.2f,%.2f], n = 16\n",
	             box.rMin, box.rMax, box.zMin, box.zMax );

	for ( int order = 1; order <= 3; ++order )
	{
		SolvedEquilibrium run( eq, box, order, 16 );

		meq::CriticalPointFinder finder( run.theSolver() );
		meq::IndexAudit const audit = finder.audit();
		std::vector<meq::CriticalPoint> const found = finder.sweep();

		char label[ 64 ];
		std::snprintf( label, sizeof( label ), "nstx diverted k=%d h=%.4f",
		               order, run.meshSize() );
		printAudit( label, audit, found );

		int extrema = 0;
		int saddles = 0;
		int sum = 0;
		double saddleGap = -1.0;
		double axisGap = -1.0;
		for ( std::size_t i = 0; i < found.size(); ++i )
		{
			sum += found[ i ].index;

			BOOST_TEST( found[ i ].fluxResidual < 1.0e-10,
			            label << ": a reported critical point has |q_h| = "
			            << found[ i ].fluxResidual << ", which is not a root" );

			if ( found[ i ].type == meq::CriticalPointType::Saddle )
			{
				++saddles;
				double const dr = found[ i ].r - exactSaddle.r;
				double const dz = found[ i ].z - exactSaddle.z;
				saddleGap = std::sqrt( dr*dr + dz*dz );
				BOOST_TEST( found[ i ].index == -1,
				            label << ": the saddle carries index "
				            << found[ i ].index << ", not -1" );
			}
			if ( found[ i ].type == meq::CriticalPointType::Maximum
			     || found[ i ].type == meq::CriticalPointType::Minimum )
			{
				++extrema;
				double const dr = found[ i ].r - exactAxisPoint.r;
				double const dz = found[ i ].z - exactAxisPoint.z;
				axisGap = std::sqrt( dr*dr + dz*dz );
				BOOST_TEST( found[ i ].index == +1,
				            label << ": the extremum carries index "
				            << found[ i ].index << ", not +1" );
			}
		}

		BOOST_TEST( extrema == 1,
		            label << ": the sweep found " << extrema
		            << " extrema on a box holding one magnetic axis" );
		BOOST_TEST( saddles == 1,
		            label << ": the sweep found " << saddles
		            << " saddles on a box holding one X-point" );

		// AND THEY ARE THE RIGHT TWO POINTS, not merely the right two indices.
		// A tenth of a cell, which is loose on purpose -- the sharp positional
		// statement is theXPointConvergesAtTheFluxesOwnOrder and this is here so
		// that "one +1 and one -1" cannot be satisfied by a pair of artefacts.
		// Measured at k = 1, the worst of the six gaps is 3.2e-3 against a
		// threshold of 7.5e-3.
		double const near = 0.1*run.meshSize();
		if ( saddles == 1 )
			BOOST_TEST( saddleGap < near,
			            label << ": the located saddle is " << saddleGap
			            << " from the closed form's X-point at ( "
			            << exactSaddle.r << ", " << exactSaddle.z << " )" );
		if ( extrema == 1 )
			BOOST_TEST( axisGap < near,
			            label << ": the located extremum is " << axisGap
			            << " from the closed form's axis at ( "
			            << exactAxisPoint.r << ", " << exactAxisPoint.z << " )" );

		// The half of Poincare-Hopf that needs no transversality, and the one a
		// diverted solve would actually use.
		BOOST_TEST( sum == audit.windingNumber,
		            label << ": the indices sum to " << sum
		            << " against a boundary degree of " << audit.windingNumber
		            << ". The sweep is not exhaustive, so this disagreement means "
		            << "a root was missed or misclassified, not that the degree "
		            << "is wrong" );

		BOOST_TEST( audit.windingNumber == 0,
		            label << ": the degree of q on a boundary enclosing one axis "
		            << "and one X-point is " << audit.windingNumber
		            << ", wanted 0. A degree is a SUM of indices" );

		BOOST_TEST( audit.windingDefect < 1.0e-9,
		            label << ": the accumulated turning is " << audit.turning
		            << " of a full turn, which is not an integer; worst single "
		            << "turn " << audit.worstTurn );

		BOOST_TEST( audit.worstTurn < 1.5,
		            label << ": a single sample turned by " << audit.worstTurn
		            << " radians, so the boundary walk is under-sampled and the "
		            << "degree is not to be believed" );

		BOOST_TEST( audit.eulerCharacteristic == 1,
		            label << ": chi of a triangulated rectangle came out "
		            << audit.eulerCharacteristic << ", not 1" );

		BOOST_TEST( audit.boundaryLoops == 1,
		            label << ": the boundary threaded into " << audit.boundaryLoops
		            << " loops, not 1" );

		// The forced pair, asserted so that a reader meeting degree 0 against
		// chi 1 is told which hypothesis is missing rather than left to wonder
		// whether the audit is broken. If either of these ever flipped, the
		// other one must have flipped too -- and a transverse boundary with a
		// degree of 0 and chi of 1 would be a genuine contradiction.
		BOOST_TEST( !audit.transverse,
		            label << ": q IS transverse to this boundary, min |q.n|/|q| = "
		            << audit.transversality
		            << ". Then Poincare-Hopf applies and a degree of "
		            << audit.windingNumber << " against chi = "
		            << audit.eulerCharacteristic
		            << " is a contradiction rather than an inapplicable "
		            << "hypothesis" );

		BOOST_TEST( !audit.consistent(),
		            label << ": the degree and chi agree on a box containing a "
		            << "saddle. One +1 and one -1 sum to 0 and chi is 1, so they "
		            << "cannot -- something has been missed" );
	}
}

/*
 * THE SEEDED SADDLE SEARCH: THE SAME X-POINT sweep() FINDS, FOR A THIRTIETH TO
 * A TWO-HUNDREDTH OF THE WORK.
 *
 * XP-0 above locates the X-point with sweep(), which is the only UNSEEDED
 * route to a saddle: findAxis() looks for an interior extremum and refuses one
 * by construction. That is the right instrument for XP-0, whose whole point is
 * that nothing is seeded -- but it costs one Newton per element per seed, and
 * an outer iteration that relocates the X-point once per Jacobian would pay
 * that on every step. tryFindCriticalPointFrom() with AxisSense::Saddle is the
 * same trade tryFindAxisFrom() makes for the axis, and this case is the same
 * measurement one critical point along.
 *
 * MOST OF WHAT IS ASSERTED HERE WOULD PASS ON A SWEEP. TWO THINGS WOULD NOT.
 *
 * That the seeded search finds a saddle, that it is the SAME saddle to
 * round-off, that it is credited to the same element and that it is classified
 * a saddle would all pass on an implementation that simply called sweep() and
 * filtered -- which is exactly what this entry point exists not to do. The two
 * that cannot are the COUNT and its TREND with refinement, and they are counts
 * rather than wall clocks for this project's standing reason: a timing on one
 * machine is a measurement about that machine, and the ratio of element-local
 * Newtons is a property of the algorithm.
 *
 * newtonSolves() counts every ( element, seed ) pair the finder attempted,
 * converged or not. sweep() attempts every element; this attempts the seed
 * element and however many rings it takes.
 *
 * SEEDED OFF THE CLOSED FORM'S X-POINT AND NOT OFF THE LOCATED ONE, half an
 * element away in each coordinate. A seed taken from the answer would make the
 * search look better than it is; the closed form's X-point is a prior that
 * exists before any of these meshes do, and 0.7 h is a step far larger than a
 * converging outer iteration takes near the end. It is also the offset
 * theSeededSearchFindsTheSameAxisForLessWork uses for the axis, so the two cost
 * tables mean the same thing.
 */
BOOST_AUTO_TEST_CASE( theSeededSaddleSearchCostsLessThanASweep )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = seededXPointBox();
	std::vector<int> const meshes = { 8, 16, 32 };

	std::printf( "\n  THE SEEDED SADDLE SEARCH AGAINST sweep()\n" );
	std::printf( "    %3s %4s %9s %12s %12s %10s %9s %8s %12s\n",
	             "k", "n", "elements", "sweep solves", "seed solves",
	             "seed elts", "saving", "depth", "gap" );

	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<double> savings;

		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			SolvedEquilibrium run( eq, box, order, meshes[ m ] );
			meq::CriticalPointFinder finder( run.theSolver() );

			// THE REFERENCE, AND ITS COST. sweep() is what XP-0 uses and what the
			// saving is measured against, so both come from the same call.
			finder.resetCounters();
			std::vector<meq::CriticalPoint> const all = finder.sweep();
			long const sweepSolves = finder.newtonSolves();
			long const sweepElements = finder.elementsRooted();

			int saddles = 0;
			meq::CriticalPoint reference;
			for ( std::size_t j = 0; j < all.size(); ++j )
			{
				if ( all[ j ].type != meq::CriticalPointType::Saddle )
					continue;
				++saddles;
				reference = all[ j ];
			}

			BOOST_TEST_REQUIRE( saddles == 1,
				"k = " << order << ", n = " << meshes[ m ] << ": the sweep found "
				<< saddles << " saddles on a box holding exactly one, so there is "
				"no reference for the seeded search to agree with." );

			// AND THE ROOT MUST BE ONE ROOT, which is a property of the box
			// rather than of the search. See seededXPointBox(): where the
			// X-point falls within q_h's jump across a face, both neighbours
			// hold a clean root of their own polynomial and the two routes are
			// entitled to return different ones. The threshold is a twentieth
			// of a reference element against a measured worst of 0.19 here, and
			// against 0.004 on xPointBox() at k = 1, n = 16, which is the
			// configuration that made this a checked precondition.
			double const depth = referenceDepth( reference );
			BOOST_TEST_REQUIRE( depth > 0.05,
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the located X-point sits " << depth
				<< " of a reference element from the boundary of the element "
				"holding it, so a face neighbour's polynomial holds its own "
				"version of the same root and \"the same root to round-off\" is "
				"not a defined question. That is a property of where the box's "
				"mesh lines fall and not of the search -- see "
				"seededXPointBox()." );

			// HALF AN ELEMENT AWAY, DIAGONALLY, from the CLOSED FORM's X-point
			// -- the same offset theSeededSearchFindsTheSameAxisForLessWork uses
			// for the axis, so the two cost tables are comparable.
			double const h = run.meshSize();
			double const seedR = nstxXPointR + 0.5*h;
			double const seedZ = nstxXPointZ - 0.5*h;

			finder.resetCounters();
			meq::CriticalPoint seeded;
			bool const found = finder.tryFindCriticalPointFrom(
				seedR, seedZ, meq::AxisSense::Saddle, seeded );
			long const seededSolves = finder.newtonSolves();
			long const seededElements = finder.elementsRooted();

			double const dr = seeded.r - reference.r;
			double const dz = seeded.z - reference.z;
			double const gap = found ? std::sqrt( dr*dr + dz*dz )
			                         : std::numeric_limits<double>::infinity();

			std::printf( "    %3d %4d %9d %12ld %12ld %10ld %8.1fx %8.3f %12.3e\n",
			             order, meshes[ m ], run.theMesh().GetNE(),
			             sweepSolves, seededSolves, seededElements,
			             seededSolves > 0
			               ? static_cast<double>( sweepSolves )/
			                 static_cast<double>( seededSolves ) : 0.0,
			             depth, gap );
			std::fflush( stdout );

			BOOST_TEST_REQUIRE( found,
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded search found no saddle half an element from the "
				"closed form's X-point at ( " << nstxXPointR << ", "
				<< nstxXPointZ << " ), where the sweep found one at ( "
				<< reference.r << ", " << reference.z << " )." );

			// ROUND-OFF, for the axis case's reason: both routes root the same
			// element's polynomial with the same Newton and the same stopping
			// rule, so they either reach the same root or they reach a
			// NEIGHBOURING element's version of it -- which is O( h^(k+1) ) away
			// and is what the containment machinery arbitrates.
			BOOST_TEST( gap <= 1.0e-12*std::max( 1.0, std::abs( reference.r ) ),
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded search returned ( " << seeded.r << ", "
				<< seeded.z << " ) where the sweep returns ( " << reference.r
				<< ", " << reference.z << " ) -- " << gap << " apart. Seeded half "
				"an element away these are the same root to the last bits; a gap "
				"of order h^(k+1) means a neighbouring element's version of it, "
				"and "
				"the depth check above is what says the box did not ask for "
				"that." );

			BOOST_TEST( seeded.element == reference.element,
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded search credits the X-point to element "
				<< seeded.element << " where the sweep credits it to "
				<< reference.element << ". The position agreed, so this is the "
				"containment tie-break landing differently -- harmless for a "
				"position and NOT harmless for a caller building a bordered "
				"Jacobian row on that element's dofs, which is what XP-3 does." );

			BOOST_TEST( ( seeded.type == meq::CriticalPointType::Saddle ),
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded search returned a "
				<< meq::criticalPointName( seeded.type )
				<< " for AxisSense::Saddle." );

			// THE ASSERTION WITH TEETH. Everything above passes on an
			// implementation that calls sweep() and filters; this does not. An
			// ORDER OF MAGNITUDE rather than merely "fewer", because "fewer"
			// would be met by a sweep that skipped one element: the worst
			// measured here is 32x, at the coarsest mesh, so ten leaves three
			// times the headroom and still says something.
			BOOST_TEST( seededSolves*10 < sweepSolves,
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded saddle search ran " << seededSolves
				<< " element-local Newtons against the sweep's " << sweepSolves
				<< " over " << run.theMesh().GetNE()
				<< " elements. The entry point exists for the cost and nothing "
				"else -- if it is not an order of magnitude cheaper than a sweep "
				"it is a sweep." );

			BOOST_TEST( seededElements < sweepElements,
				"k = " << order << ", n = " << meshes[ m ]
				<< ": the seeded search entered " << seededElements
				<< " elements against the sweep's " << sweepElements << "." );

			savings.push_back( seededSolves > 0
			                   ? static_cast<double>( sweepSolves )/
			                     static_cast<double>( seededSolves ) : 0.0 );
		}

		// AND THE SAVING IMPROVES WITH REFINEMENT, which is the property an
		// outer iteration needs and is a different statement from any single
		// row of the table. The search is bounded by rings -- an absolute
		// element count set by how far the seed is out, not by how many
		// elements there are -- while a sweep is the whole mesh. So the
		// per-Jacobian cost of relocating the X-point does not grow with the
		// problem, and a case that only ever measured one mesh could not say
		// so. Measured 32x, 128x, 205x over n = 8, 16, 32.
		BOOST_TEST( savings.back() > savings.front(),
			"k = " << order << ": the seeded search saves " << savings.front()
			<< "x on the coarsest mesh and " << savings.back()
			<< "x on the finest. It is bounded by rings and the sweep is bounded "
			"by the mesh, so refining must widen the gap; a saving that shrinks "
			"means the ring growth has started tracking the element count." );
	}
}

/*
 * AND IT DECLINES WHERE THERE IS NO SADDLE, WHICH IS WHAT MAKES THE ABOVE
 * USABLE.
 *
 * tryFindCriticalPointFrom() takes the point of the requested sense NEAREST its
 * seed where the seed region offers several -- right for a caller following one
 * critical point and badly wrong as a way to DISCOVER one -- so the contract
 * that matters is what it does when there is nothing of that sense to find.
 *
 * THREE WAYS OF HAVING NOTHING TO FIND, on a box that holds an axis AND an
 * X-point so that the search is never merely staring at an empty field:
 *
 *   * seeded ON the magnetic axis. There is a critical point directly under the
 *     seed and it is an extremum, and the X-point is 1.73 away in z -- far
 *     outside six rings. This is the mirror of
 *     theSeededSearchDeclinesRatherThanGuessing's sense trap, and it is the one
 *     that matters most for an outer iteration: an X-point that has wandered
 *     off must be reported missing, not replaced by the axis;
 *   * seeded in a CORNER, where q vanishes nowhere within reach;
 *   * seeded on the X-point asking for an EXTREMUM, which must fail for the
 *     same reason read the other way.
 *
 * A Degenerate is never an acceptable answer to any of them: it is an admission
 * that the determinant is at round-off and no classification is entitled, and a
 * caller handed one as an X-point would follow it.
 *
 * AND THE REFUSALS AT THE OTHER TWO ENTRY POINTS ARE ASSERTED HERE TOO, because
 * they are refusals of a DIFFERENT kind and mixing them up is the trap. A seed
 * with no saddle near it is a fact about the equilibrium and returns false; an
 * AxisSense::Saddle handed to findAxis() is a mistake in the call -- its seeds
 * are the extreme nodal values of psi_h, which is where an axis is and is not
 * where an X-point is -- and answering "not found" to it would report the
 * absence of a saddle on the strength of never having looked.
 */
BOOST_AUTO_TEST_CASE( theSeededSaddleSearchDeclinesRatherThanGuessing )
{
	Equilibrium const eq = Equilibrium::nstx();
	Rectangle const box = divertedBox();
	int const order = 2;
	int const n = 16;

	SolvedEquilibrium run( eq, box, order, n );
	meq::CriticalPointFinder finder( run.theSolver() );

	ExactCritical const exactSaddle =
		exactCriticalPoint( eq, nstxXPointR, nstxXPointZ );
	ExactCritical const exactAxis = exactCriticalPoint( eq, 1.0, 0.0 );

	std::printf( "\n  THE SEEDED SADDLE SEARCH WHEN THERE IS NO SADDLE TO FIND "
	             "( k = %d, n = %d )\n", order, n );
	std::printf( "    the axis is at ( %.4f, %.4f ), the X-point at "
	             "( %.4f, %.4f ), %.2f apart in z\n",
	             exactAxis.r, exactAxis.z, exactSaddle.r, exactSaddle.z,
	             std::abs( exactAxis.z - exactSaddle.z ) );

	// THE CONTROL FIRST, because without it every refusal below is equally well
	// explained by a search that never finds anything.
	meq::CriticalPoint onTarget;
	bool const foundTarget = finder.tryFindCriticalPointFrom(
		exactSaddle.r, exactSaddle.z, meq::AxisSense::Saddle, onTarget );
	std::printf( "    seeded ON the X-point, asking for a saddle: %s\n",
	             foundTarget ? "found" : "DECLINED" );

	BOOST_TEST_REQUIRE( foundTarget,
		"the control failed: seeded on the X-point and asked for a saddle the "
		"search found nothing, so the three refusals below say nothing about the "
		"sense filter or the search radius." );

	double const controlDr = onTarget.r - exactSaddle.r;
	double const controlDz = onTarget.z - exactSaddle.z;
	BOOST_TEST( std::sqrt( controlDr*controlDr + controlDz*controlDz )
	            < 0.1*run.meshSize(),
		"the control found a saddle at ( " << onTarget.r << ", " << onTarget.z
		<< " ) rather than the X-point at ( " << exactSaddle.r << ", "
		<< exactSaddle.z << " )." );

	// (1) ON THE AXIS, ASKING FOR A SADDLE. There is a critical point under the
	// seed; it is the wrong kind, and the right kind is 1.73 away in z.
	meq::CriticalPoint atAxis;
	bool const foundAtAxis = finder.tryFindCriticalPointFrom(
		exactAxis.r, exactAxis.z, meq::AxisSense::Saddle, atAxis );
	std::printf( "    seeded ON the axis, asking for a saddle: %s\n",
	             foundAtAxis ? "FOUND SOMETHING" : "declined" );

	BOOST_TEST( !foundAtAxis,
		"seeded exactly on the magnetic axis and asked for a saddle, the search "
		"returned a " << meq::criticalPointName( atAxis.type ) << " at ( "
		<< atAxis.r << ", " << atAxis.z << " ), element " << atAxis.element
		<< ". The X-point is at ( " << exactSaddle.r << ", " << exactSaddle.z
		<< " ), far outside the seed rings, so whatever this is it is not it -- "
		"and an outer iteration following an X-point would have been handed it "
		"silently." );

	// AND NOTHING WAS WRITTEN INTO THE RESULT, which is the contract's other
	// half and is what stops a declining search from handing back a Degenerate.
	// A default CriticalPoint carries element -1 and type Degenerate; if the
	// search ever wrote a candidate it rejected, this is what would see it.
	BOOST_TEST( atAxis.element == -1,
		"the search declined and still wrote element " << atAxis.element
		<< " into the result, a " << meq::criticalPointName( atAxis.type )
		<< " at ( " << atAxis.r << ", " << atAxis.z
		<< " ). A caller who does not check the return value would follow it." );

	// (2) FAR AWAY, in the corner of the box.
	meq::CriticalPoint corner;
	bool const foundCorner = finder.tryFindCriticalPointFrom(
		box.rMin, box.zMax, meq::AxisSense::Saddle, corner );
	std::printf( "    seeded at the ( %.2f, %.2f ) corner, asking for a saddle: "
	             "%s\n", box.rMin, box.zMax,
	             foundCorner ? "FOUND SOMETHING" : "declined" );

	BOOST_TEST( !foundCorner,
		"seeded in the corner at ( " << box.rMin << ", " << box.zMax
		<< " ) the search returned a " << meq::criticalPointName( corner.type )
		<< " at ( " << corner.r << ", " << corner.z
		<< " ). Six rings of face neighbours do not reach the X-point from "
		"there. If setSeedRings()' default has been widened deliberately, this "
		"case pins the SHIPPED default and is what has to move." );

	// (3) ON THE X-POINT, ASKING FOR AN EXTREMUM. The mirror of (1).
	meq::CriticalPoint mistyped;
	bool const foundMistyped = finder.tryFindCriticalPointFrom(
		exactSaddle.r, exactSaddle.z, meq::AxisSense::Either, mistyped );
	std::printf( "    seeded ON the X-point, asking for an extremum: %s\n",
	             foundMistyped ? "FOUND SOMETHING" : "declined" );
	std::fflush( stdout );

	BOOST_TEST( !foundMistyped,
		"seeded exactly on the X-point and asked for an extremum, the search "
		"returned a " << meq::criticalPointName( mistyped.type ) << " at ( "
		<< mistyped.r << ", " << mistyped.z
		<< " ). AxisSense::Either is EITHER EXTREMUM and never a saddle; a "
		"default that had quietly grown to mean 'any critical point' would make "
		"findAxis() return the X-point of a diverted equilibrium as the "
		"magnetic axis, at the right rate, with every table intact." );

	// AND THE OTHER KIND OF REFUSAL, which is about the CALL rather than about
	// the field. Both of these seed from the extreme nodal values of psi_h, so
	// a saddle is not something they could look for.
	BOOST_CHECK_THROW( finder.findAxis( meq::AxisSense::Saddle ),
	                   std::invalid_argument );

	meq::CriticalPoint unused;
	BOOST_CHECK_THROW( finder.tryFindAxis( unused, meq::AxisSense::Saddle ),
	                   std::invalid_argument );

	// tryFindAxisFrom() COULD serve it -- it is tryFindCriticalPointFrom() --
	// and refuses anyway, so that "axis" keeps meaning axis at the one entry
	// point a caller with a prior reaches for.
	BOOST_CHECK_THROW( finder.tryFindAxisFrom( exactSaddle.r, exactSaddle.z,
	                                           meq::AxisSense::Saddle, unused ),
	                   std::invalid_argument );
}

/*
 * ===========================================================================
 * THE SWEEP IS THREADED AND THE ANSWER IS NOT ALLOWED TO KNOW
 * ===========================================================================
 *
 * sweep() roots every element of the mesh, which is embarrassingly parallel --
 * and the deduplication that follows is not, because the merge rule below
 * `points[ j ].type == point.type && distance < reach` keeps "the one least
 * outside its own element" and compares each candidate against the list built
 * SO FAR. That list is in element order. Deduplicating inside the parallel
 * region would make it thread-arrival order, and the measurement recorded on
 * that rule says exactly what that costs: on the Solov'ev benchmark at
 * k = 1, n = 6 the candidate strictly inside its element is 2.7e-3 from the
 * true axis while the neighbour sitting 8.5e-2 outside is 6.1e-3, and the
 * neighbour has the LOWER element index. So the two candidates are 3.4e-3
 * apart in the answer and which one survives is decided by the order they are
 * seen in.
 *
 * **THAT IS THE FAILURE THIS CASE EXISTS FOR AND IT IS THE QUIET KIND.** The
 * axis would move at the 1e-3 level with OMP_NUM_THREADS, which on a coarse
 * mesh reads as discretisation error rather than as a bug, and every
 * convergence table in this file would still pass.
 *
 * So the assertion is EQUALITY AT 0.000e+00, entry for entry, in every field
 * of CriticalPoint -- not a tolerance. The design claims the threaded sweep
 * reproduces the serial one exactly, and an assertion should say what the
 * design claims.
 *
 * k = 1, n = 6 IS THE MESH NAMED IN THAT MEASUREMENT, and it is used here for
 * that reason: it is a configuration where the sweep is KNOWN to find
 * duplicates the merge has to choose between. Two finer meshes ride along so
 * that a case which found no duplicates at all -- and so could not fail --
 * would show up as a suspiciously small candidate count rather than as a pass.
 */
BOOST_AUTO_TEST_CASE( theAxisSweepDoesNotDependOnTheThreadCount )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	BOOST_TEST_MESSAGE( "MFEM is built without MFEM_USE_OPENMP and "
	                    "MFEM_THREAD_SAFE together, so meq's sweep is serial "
	                    "by construction and there is nothing to compare" );
#else
	Equilibrium const eq = Equilibrium::nstx();

	int const ambient = omp_get_max_threads();
	std::printf( "\n  THE SWEEP AT ONE THREAD AND AT %d\n", ambient );
	std::printf( "    %-14s %-4s %-4s %10s %12s %12s\n", "box", "k", "n",
	             "candidates", "newton", "elements" );
	std::fflush( stdout );

	/*
	 * THE FIRST ROW IS THE ONE WITH TEETH AND THE FILE ALREADY SAYS WHY.
	 *
	 * seededXPointBox()'s own comment records the measurement: on xPointBox()
	 * at k = 1, n = 16 the X-point at r = 0.699700 sits 3.0e-04 from the mesh
	 * line at r = 0.700000, element 237 holds a root at r = 0.699826 and
	 * element 238 holds its own at r = 0.700027, BOTH STRICTLY INSIDE and
	 * 6.9e-04 apart -- so both have `overshoot` exactly zero, the merge rule's
	 * `point.overshoot < points[ j ].overshoot` is FALSE either way, and
	 * "sweep() merges the two and keeps whichever it saw first". That is the
	 * one configuration in this file where the surviving candidate is decided
	 * by the ORDER the candidates are seen in and by nothing else.
	 *
	 * **WITHOUT IT THIS CASE CANNOT FAIL.** Checked rather than assumed: with
	 * the concatenation deliberately put into thread-arrival order the axis
	 * rows below stay green, because there the two candidates have DIFFERENT
	 * overshoots and the merge rule is decisive -- which is what the rule is
	 * for. A test that cannot fail is worse than no test, so the fixture that
	 * makes it fail is the first row and the others ride along as controls.
	 */
	struct Case
	{
		char const *name;
		Rectangle box;
		int order;
		int n;
	};
	std::vector<Case> const cases = {
		{ "xPoint", xPointBox(), 1, 16 },
		{ "xPoint", xPointBox(), 1, 8 },
		{ "standard", meq::tests::standardBox(), 1, 6 },
		{ "standard", meq::tests::standardBox(), 2, 10 },
	};

	for ( Case const &one : cases )
	{
			SolvedEquilibrium run( eq, one.box, one.order, one.n );
			meq::CriticalPointFinder finder( run.theSolver() );

			omp_set_num_threads( 1 );
			finder.resetCounters();
			std::vector<meq::CriticalPoint> const serial = finder.sweep();
			long const serialNewton = finder.newtonSolves();
			long const serialElements = finder.elementsRooted();

			omp_set_num_threads( ambient );
			finder.resetCounters();
			std::vector<meq::CriticalPoint> const threaded = finder.sweep();
			long const threadedNewton = finder.newtonSolves();
			long const threadedElements = finder.elementsRooted();

			std::printf( "    %-14s %-4d %-4d %10zu %12ld %12ld\n", one.name,
			             one.order, one.n, serial.size(), serialNewton,
			             serialElements );
			std::fflush( stdout );

			BOOST_TEST_REQUIRE( serial.size() == threaded.size(),
				one.name << " k = " << one.order << ", n = " << one.n
				<< ": the sweep returned " << serial.size()
				<< " candidates at one thread and " << threaded.size() << " at "
				<< ambient
				<< ". The deduplication must run serially over the candidates "
				"concatenated in ELEMENT order; a different count means it did "
				"not." );

			// THE COUNTERS TOO, because they are what a caller measures a
			// seeded search against a sweep by, and an atomic increment is
			// exactly as order independent as a sum of ones.
			BOOST_TEST( serialNewton == threadedNewton,
				one.name << " k = " << one.order << ", n = " << one.n << ": "
				<< serialNewton << " element-local Newton solves at one thread "
				"against " << threadedNewton << " at " << ambient );
			BOOST_TEST( serialElements == threadedElements,
				one.name << " k = " << one.order << ", n = " << one.n << ": "
				<< serialElements << " elements entered at one thread against "
				<< threadedElements << " at " << ambient );

			for ( std::size_t i = 0; i < serial.size(); ++i )
			{
				meq::CriticalPoint const &a = serial[ i ];
				meq::CriticalPoint const &b = threaded[ i ];

				// EVERY FIELD, AND AT ZERO. Position is what a reader looks at,
				// but `element` and `overshoot` are what the merge rule DECIDES
				// on -- so a dedup that picked the other candidate would be
				// caught here even where the two sit at the same place to
				// printing precision.
				BOOST_TEST( a.r - b.r == 0.0,
					one.name << " k = " << one.order << ", n = " << one.n
					<< ", candidate " << i << ": r moved by " << ( a.r - b.r )
					<< " between one thread and " << ambient );
				BOOST_TEST( a.z - b.z == 0.0,
					one.name << " k = " << one.order << ", n = " << one.n
					<< ", candidate " << i << ": z moved by " << ( a.z - b.z ) );
				BOOST_TEST( a.psi - b.psi == 0.0,
					one.name << " k = " << one.order << ", n = " << one.n
					<< ", candidate " << i << ": psi moved by "
					<< ( a.psi - b.psi ) );
				BOOST_TEST( a.fluxResidual - b.fluxResidual == 0.0 );
				BOOST_TEST( a.determinant - b.determinant == 0.0 );
				BOOST_TEST( a.trace - b.trace == 0.0 );
				BOOST_TEST( a.overshoot - b.overshoot == 0.0,
					one.name << " k = " << one.order << ", n = " << one.n
					<< ", candidate " << i
					<< ": overshoot moved, so the merge kept a DIFFERENT "
					"candidate at the two thread counts. That is the 1e-3 "
					"failure this case exists for." );
				BOOST_TEST( a.referenceX - b.referenceX == 0.0 );
				BOOST_TEST( a.referenceY - b.referenceY == 0.0 );
				BOOST_TEST( a.element == b.element,
					one.name << " k = " << one.order << ", n = " << one.n
					<< ", candidate " << i << ": element " << a.element
					<< " at one thread and " << b.element << " at " << ambient
					<< ". The merge chose differently." );
				bool const sameType = ( a.type == b.type );
				BOOST_TEST( sameType );
				BOOST_TEST( a.index == b.index );
			}
	}

	omp_set_num_threads( ambient );
#endif
}
