#define BOOST_TEST_MODULE CriticalPointConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/GradShafranov.hpp"

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
 */

namespace
{

	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	/// The magnetic axis of the CLOSED FORM, to round-off.
	///
	/// Newton on grad( psi ) = 0 with the Hessian by central differences. The
	/// Hessian's accuracy does not reach the answer -- the fixed point of this
	/// iteration is where the analytic gradient vanishes, whatever steered it
	/// there -- which is the same observation CriticalPoints.hpp makes about its
	/// own Jacobian, made here in a place where it can be checked independently.
	struct ExactAxis
	{
		double r;
		double z;
		double psi;
		double gradient;
		double determinant;
		double trace;
	};

	ExactAxis exactAxis( Equilibrium const &eq, double rGuess, double zGuess )
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

		ExactAxis axis;
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

	void printAxisTable( int order, ExactAxis const &exact,
	                     std::vector<AxisMeasurement> const &points )
	{
		std::printf( "\n  Magnetic axis as a zero of q, Solov'ev NSTX, k = %d\n", order );
		std::printf( "  exact axis ( %.12f, %.12f ), |grad psi| = %.2e, %s\n",
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
		ExactAxis const axis = exactAxis( cases[ i ].eq, 1.0, 0.0 );
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
	ExactAxis const exact = exactAxis( eq, 1.0, 0.0 );

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

		printAxisTable( order, exact, points );

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
