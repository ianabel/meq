#define BOOST_TEST_MODULE SurfaceFitConvergence
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
#include "meq/SurfaceFit.hpp"

#include "analytic/FluxSurfaceReference.hpp"
#include "analytic/Soloviev.hpp"
#include "ConvergenceHarness.hpp"

/*
 * The acceptance test for INVERSION-PLAN.md stage IN-3: the representation.
 * R( discRadius, theta ) and z( discRadius, theta ) as truncated Zernike
 * expansions on the disc, fitted to traced surface POINTS.
 *
 * IT RUNS ON AN ANALYTIC EQUILIBRIUM AND THAT IS THE POINT RATHER THAN A
 * CONVENIENCE. How fast a REPRESENTATION converges is a property of the
 * representation and not of the discretisation, so the surfaces here come from
 * tests/analytic/FluxSurfaceReference.hpp -- rays plus a safeguarded root solve
 * on the exact psi, converged to 1e-15 -- and the mesh is sidestepped entirely.
 * It also sidesteps a real limitation rather than a hypothetical one: IN-2
 * found that a fitted rectangle cannot hold Psi_N = 0.75 on nstx() at all, so a
 * discrete study is boxed in where an analytic one is not. The discrete leg is
 * the LAST case in this file and it exists to show the library works on the
 * consumer's own path, not to measure a rate.
 *
 * ---------------------------------------------------------------------------
 * THE HEADLINE FINDING, AND IT WAS NOT WHAT THIS STAGE WAS BRIEFED TO LOOK FOR.
 * ---------------------------------------------------------------------------
 *
 * THE ANGULAR LABEL DECIDES WHETHER THE BASIS CONVERGES AT ALL, AND THE OBVIOUS
 * LABEL -- THE GEOMETRIC POLOIDAL ANGLE ABOUT THE MAGNETIC AXIS, WHICH IS WHAT
 * BOTH OF MEQ'S SURFACE SAMPLERS PRODUCE -- IS THE WRONG ONE.
 *
 * A function smooth in Cartesian ( x, y ) has exactly ONE angular harmonic
 * multiplying discRadius^1, because its expansion begins f( 0 ) + a x + b y.
 * Parametrise a family of nested ELLIPSES by the geometric angle and the
 * coefficient of discRadius^1 carries cos( theta ), cos( 3 theta ),
 * cos( 5 theta ) and so on for ever -- and every one of those beyond the first
 * is a mode the Zernike index constraint excludes at radial degree one,
 * precisely because it is not smooth at the origin. So the parametrisation puts
 * content exactly where the basis refuses to look.
 *
 * theAngleLabelDecidesWhetherTheBasisConvergesAtAll is that statement with the
 * answer known exactly: the same ellipse points fit to 3e-15 at degree TWO
 * under their own parameter and are wrong by 1e-2 at degree TWENTY under the
 * geometric angle. Same points, same basis, same solver.
 *
 * meq::relabelByAxisShape() is the repair and it needs no field: the innermost
 * traced surface IS the axis ellipse to O( Psi_N ), so a three-parameter linear
 * fit of a quadratic form to it gives the tilt and the flattening, and the
 * relabelling is then an exact reparametrisation of the circle.
 *
 * AND IT IS EXACT AT ONE ORDER AND NO BETTER, WHICH IS THE SECOND FINDING. A
 * relabelling that does not depend on discRadius has one function's worth of
 * freedom and can therefore fix one order. The discRadius^2 and higher
 * coefficients keep whatever harmonics the surface shaping puts there, so the
 * fit stays limited near the axis; measured, the envelope decays geometrically
 * with the innermost surface at Psi_N = 0.10 and leaves an algebraic tail when
 * it is pulled in to Psi_N = 0.02. That is INVERSION-PLAN.md section 4.4's
 * tension -- axis regularity against everything else -- as a table, and it is
 * IN-4's to resolve rather than this stage's.
 *
 * ---------------------------------------------------------------------------
 * WHAT IS ASSERTED, AND WHY EACH CONTROL IS LIVE RATHER THAN QUOTED.
 * ---------------------------------------------------------------------------
 *
 * 1. THE ANGLE, above. The control is the identical point cloud relabelled.
 *
 * 2. THE RADIAL COORDINATE. Zernike.hpp argues that parametrising by Psi_N
 *    rather than by sqrt( Psi_N ) puts a square-root branch point at the axis
 *    and that every basis then converges ALGEBRAICALLY against it, "with
 *    nothing in a convergence table to say why". That was an argument and this
 *    is the measurement. FitRadialCoordinate::NormalisedFlux is the losing
 *    branch and it is a first-class option, not test scaffolding.
 *
 * 3. THE PARITY CONSTRAINT. FitBasis::TensorProduct admits any power of the
 *    radius times any Fourier mode, so it admits cos( theta ) with no radial
 *    factor at all -- which has no limit at the origin. IT FITS THE SAMPLE
 *    CLOUD BETTER THAN ZERNIKE DOES, on more modes, and it is useless: its
 *    design matrix condition number is 1e16 and the point it puts at the axis
 *    is a CURVE rather than a point. That is the shape CLAUDE.md keeps
 *    recording -- a defect no residual column can see.
 *
 * 4. THE DERIVATIVE, WHICH IS THE ONE THAT MATTERS. MANTA-COUPLING.md needs
 *    dGeometry/dpsi and that is the reason a representation exists at all, so
 *    the fit is differentiated and compared against a difference of surfaces
 *    traced INDEPENDENTLY -- not against a difference of the fit's own output,
 *    which would check nothing.
 *
 *    IT IS STATED IN discRadius AND NOT IN Psi_N. INVERSION-PLAN.md section 4.4
 *    records that d/dPsi_N diverges like 1/( 2 sqrt( Psi_N ) ) at the axis in
 *    EVERY equilibrium -- the m = 1 coefficient of R is the minor radius, which
 *    grows like sqrt( Psi ) -- so an acceptance written over the whole domain
 *    against dGeometry/dPsi would be comparing infinities and would fail for a
 *    reason with nothing to do with the fit. d/d( discRadius ) is finite
 *    everywhere and is what is asserted; the conversion is checked separately,
 *    on a flux range bounded away from the axis, and its GROWTH as that bound
 *    is lowered is reported because that is the number a future coupling reads
 *    to decide where its innermost node may sit.
 *
 *    AND THE DIFFERENCE IS RICHARDSON-EXTRAPOLATED, ( 4 D( h/2 ) - D( h ) )/3.
 *    A plain central difference carries its own O( h^2 ) truncation, so the
 *    comparison is floored by the INSTRUMENT rather than by the derivative.
 *    This tree has measured that three times already -- Zernike's own
 *    derivative test at 1.3e-07 plain against 1.4e-11 extrapolated, and IN-2's
 *    identity control flat to three figures across a sixteenfold refinement --
 *    and the plain column is kept live here so it is a fourth measurement and
 *    not a fourth citation.
 *
 * 5. THE METRIC, cheap and with real teeth: the fit's own | dx/dtheta | against
 *    the pointwise identity | dx/dtheta | = sqrt( rho'^2 + rho^2 ) built from
 *    the field's gradient, which is IN-1's quantity. Two independent routes to
 *    the same number, one of them geometric and one of them the field.
 *
 * 6. THE AXIS. Every mode with m != 0 vanishes at discRadius = 0, so a Zernike
 *    fit puts every surface's collapse point at ONE place by construction, to
 *    round-off, whatever its coefficients are. That is asserted as an exact
 *    zero rather than a tolerance. What is NOT free is whether that point is
 *    the magnetic axis -- it is an extrapolation into the hole the sample set
 *    leaves -- and that error is measured against the closed form.
 *
 * 7. CONDITIONING, WHERE THE PLAUSIBLE ANSWER IS WRONG. Zernike is orthogonal
 *    on the disc under discRadius d( discRadius ) = d( Psi_N )/2, so
 *    Gauss-Legendre LEVELS IN Psi_N with equispaced angles ought to make the
 *    discrete inner product the continuous one and the fit a projection.
 *    MEASURED, THE LAYOUT MOVES NOTHING -- three layouts agree to within 25% at
 *    every inner limit tried, and Gauss is sometimes the worst of them, because
 *    the orthogonality argument needs the nodes to span the whole disc and a
 *    sample set has a HOLE in the middle. What does move it, by four orders and
 *    by eleven thousand respectively, is how big that hole is and whether the
 *    disc edge is rescaled onto the outermost surface. The wrong answer is
 *    asserted DEAD, not merely omitted: if a layout ever moves the conditioning
 *    by an order, the paragraph in SurfaceFit.hpp saying it does not is stale.
 *
 * 8. AND THE WHOLE CHAIN ON A SOLVED FIELD -- solve, root the axis, trace, ray
 *    fit, least squares -- against the exact surfaces, so that nothing above is
 *    a property of the analytic shortcut alone.
 *
 * 9. THE TRACER REPAIR IN-3 NEEDED BEFORE IT COULD START, which is not a test
 *    of the fit at all and says so where it stands.
 */

namespace
{

	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	double const twoPi = 6.283185307179586476925286766559;

	/// The magnetic axis of the CLOSED FORM, to round-off. Newton on
	/// grad( psi ) = 0 with the Hessian by central differences; the Hessian's
	/// accuracy does not reach the answer, the fixed point being where the
	/// analytic gradient vanishes whatever steered it there. Lifted from
	/// SurfaceAverageConvergence.cpp, which lifted it from
	/// FluxSurfaceConvergence.cpp, which is where the argument for it is.
	struct ExactAxis
	{
		double r;
		double z;
		double psi;
	};

	void hessianOf( Equilibrium const &eq, double r, double z, double h[ 2 ][ 2 ] )
	{
		double const step = 1.0e-5;
		double a0 = 0.0;
		double a1 = 0.0;
		double b0 = 0.0;
		double b1 = 0.0;

		eq.gradPsi( r + step, z, a0, b0 );
		eq.gradPsi( r - step, z, a1, b1 );
		h[ 0 ][ 0 ] = ( a0 - a1 )/( 2.0*step );
		h[ 1 ][ 0 ] = ( b0 - b1 )/( 2.0*step );
		eq.gradPsi( r, z + step, a0, b0 );
		eq.gradPsi( r, z - step, a1, b1 );
		h[ 0 ][ 1 ] = ( a0 - a1 )/( 2.0*step );
		h[ 1 ][ 1 ] = ( b0 - b1 )/( 2.0*step );
	}

	ExactAxis exactAxis( Equilibrium const &eq, double rGuess, double zGuess )
	{
		double r = rGuess;
		double z = zGuess;

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );

			double hessian[ 2 ][ 2 ];
			hessianOf( eq, r, z, hessian );

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

	/// THE NORMALISATION IS FluxSurfaceReference.hpp's: psi = 0 on the
	/// separatrix and psi_ax on the axis, so psi( Psi_N ) = psi_ax( 1 - Psi_N ).
	double levelAt( ExactAxis const &axis, double normalisedFlux )
	{
		return axis.psi*( 1.0 - normalisedFlux );
	}

	/// One surface on the EXACT field, at equispaced geometric poloidal angle.
	meq::analytic::SurfaceQuadrature exactSurface( Equilibrium const &eq,
	                                              ExactAxis const &axis,
	                                              double normalisedFlux,
	                                              int angles )
	{
		// reach 3.0 and four thousand march steps rather than the fixture's
		// defaults: the surfaces wanted here run from Psi_N = 0.02, whose minor
		// radius is about 0.04, up to Psi_N = 0.9, which reaches past r = 1.7.
		// One reach has to bracket both, so the march has to be fine enough for
		// the small one and long enough for the large one. The march only ever
		// BRACKETS; every digit comes from the safeguarded Newton beneath it.
		return meq::analytic::surfaceQuadrature( eq, axis.r, axis.z,
		                                         levelAt( axis, normalisedFlux ),
		                                         angles, 3.0, 4000 );
	}

	/// ONE point of the surface psi = level along the ray at a GIVEN geometric
	/// angle, which the fixture cannot supply because it samples equispaced
	/// angles only. Same safeguarded Newton, same bracket-by-marching, so the
	/// point carries the same 1e-15.
	///
	/// This is what makes the derivative acceptance a comparison against
	/// INDEPENDENTLY traced surfaces rather than against the fit's own output:
	/// the four surfaces the Richardson difference needs are each solved from
	/// scratch and the fit never enters.
	void rayPoint( Equilibrium const &eq, ExactAxis const &axis,
	               double normalisedFlux, double theta, double &r, double &z )
	{
		double const level = levelAt( axis, normalisedFlux );
		double const cosine = std::cos( theta );
		double const sine = std::sin( theta );

		auto along = [ & ]( double distance )
		{
			return eq.psi( axis.r + distance*cosine, axis.z + distance*sine )
			       - level;
		};

		double lower = 0.0;
		double lowerValue = along( 0.0 );
		double upper = -1.0;

		for ( int step = 1; step <= 4000; ++step )
		{
			double const distance = 3.0*step/4000.0;
			double const value = along( distance );
			if ( lowerValue*value <= 0.0 )
			{
				upper = distance;
				break;
			}
			lower = distance;
			lowerValue = value;
		}

		if ( !( upper > 0.0 ) )
			throw std::runtime_error( "rayPoint: the ray did not bracket the level" );

		double distance = 0.5*( lower + upper );

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double const rr = axis.r + distance*cosine;
			double const zz = axis.z + distance*sine;
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( rr, zz, gr, gz );

			double const value = eq.psi( rr, zz ) - level;
			double const slope = gr*cosine + gz*sine;

			if ( value*lowerValue > 0.0 )
			{
				lower = distance;
				lowerValue = value;
			}
			else
			{
				upper = distance;
			}

			double next = distance;
			if ( std::abs( slope ) > 0.0 )
				next = distance - value/slope;
			if ( !( next > lower && next < upper ) )
				next = 0.5*( lower + upper );

			double const move = std::abs( next - distance );
			distance = next;
			if ( move < 1.0e-16*( 1.0 + std::abs( distance ) ) )
				break;
		}

		r = axis.r + distance*cosine;
		z = axis.z + distance*sine;
	}

	/// How the levels of a sample set are laid out in Psi_N. See case 7: this is
	/// a conditioning choice and nothing else.
	enum class Layout
	{
		EquispacedFlux,   ///< equispaced in Psi_N. The obvious one.
		EquispacedRadius, ///< equispaced in sqrt( Psi_N ), the disc coordinate.
		GaussFlux         ///< Gauss-Legendre in Psi_N, which is the disc measure.
	};

	char const *layoutName( Layout layout )
	{
		switch ( layout )
		{
			case Layout::EquispacedFlux:   return "equispaced in Psi_N";
			case Layout::EquispacedRadius: return "equispaced in discRadius";
			case Layout::GaussFlux:        return "Gauss-Legendre in Psi_N";
		}
		return "unknown";
	}

	/// Gauss-Legendre nodes on [ -1, 1 ], by Newton on the Legendre polynomial
	/// with its own three-term recurrence supplying value and derivative. Only
	/// the nodes are wanted; the weights would be a least-squares weighting and
	/// case 7 measures the layout alone.
	std::vector<double> gaussNodes( int count )
	{
		std::vector<double> nodes( static_cast<std::size_t>( count ), 0.0 );

		for ( int i = 0; i < count; ++i )
		{
			double x = std::cos( M_PI*( i + 0.75 )/( count + 0.5 ) );

			for ( int iteration = 0; iteration < 100; ++iteration )
			{
				double previous = 1.0;
				double current = x;
				for ( int n = 2; n <= count; ++n )
				{
					double const next = ( ( 2*n - 1 )*x*current
					                      - ( n - 1 )*previous )/n;
					previous = current;
					current = next;
				}

				double const derivative = count*( x*current - previous )
				                          /( x*x - 1.0 );
				double const step = current/derivative;
				x -= step;
				if ( std::abs( step ) < 1.0e-15 )
					break;
			}

			nodes[ static_cast<std::size_t>( i ) ] = x;
		}

		return nodes;
	}

	std::vector<double> levelsFor( Layout layout, double smallest, double largest,
	                               int count )
	{
		std::vector<double> levels;

		if ( layout == Layout::GaussFlux )
		{
			std::vector<double> const nodes = gaussNodes( count );
			for ( double node : nodes )
				levels.push_back( 0.5*( smallest + largest )
				                  + 0.5*( largest - smallest )*node );
			std::sort( levels.begin(), levels.end() );
			return levels;
		}

		for ( int i = 0; i < count; ++i )
		{
			double const fraction = static_cast<double>( i )/( count - 1.0 );

			if ( layout == Layout::EquispacedFlux )
			{
				levels.push_back( smallest + ( largest - smallest )*fraction );
			}
			else
			{
				double const radius = std::sqrt( smallest )
					+ ( std::sqrt( largest ) - std::sqrt( smallest ) )*fraction;
				levels.push_back( radius*radius );
			}
		}

		return levels;
	}

	/// The sample cloud on the exact field, at the geometric poloidal angle.
	/// Relabelling is a separate step, deliberately: every case that compares
	/// the two labels has to be comparing the SAME points.
	std::vector<meq::SurfaceSample> exactSamples( Equilibrium const &eq,
	                                              ExactAxis const &axis,
	                                              std::vector<double> const &levels,
	                                              int angles )
	{
		std::vector<meq::SurfaceSample> samples;

		for ( double normalisedFlux : levels )
		{
			meq::analytic::SurfaceQuadrature const surface =
				exactSurface( eq, axis, normalisedFlux, angles );

			for ( meq::analytic::SurfacePoint const &point : surface.points )
			{
				meq::SurfaceSample sample;
				sample.normalisedFlux = normalisedFlux;
				sample.theta = point.theta;
				sample.r = point.r;
				sample.z = point.z;
				samples.push_back( sample );
			}
		}

		return samples;
	}

	/// The worst distance between the fit and the points it was fitted to.
	double worstFitError( meq::SurfaceFit const &fit,
	                      std::vector<meq::SurfaceSample> const &samples )
	{
		double worst = 0.0;

		for ( meq::SurfaceSample const &sample : samples )
		{
			double r = 0.0;
			double z = 0.0;
			fit.position( sample.normalisedFlux, sample.theta, r, z );
			worst = std::max( worst, std::hypot( r - sample.r, z - sample.z ) );
		}

		return worst;
	}

	/// The envelope over a WINDOW of two degrees rather than one, and the reason
	/// is a real trap: a surface family with an up-down or left-right symmetry
	/// puts all of its content on one parity of m, so the envelope at every
	/// degree of the other parity is an exact zero and a per-degree ratio
	/// alternates between 0 and infinity. Measured, the synthetic ellipse family
	/// of case 1 reads 1e-15 at every even degree while its fit error is 1e-2.
	double envelopeWindow( meq::SurfaceFit const &fit, int degree )
	{
		double envelope = fit.coefficientEnvelope( degree );
		if ( degree > 0 )
			envelope = std::max( envelope, fit.coefficientEnvelope( degree - 1 ) );

		return envelope;
	}

	/// The benchmark box for the discrete leg. NOT standardBox(): nstx()'s
	/// surfaces are elongated and Psi_N = 0.50 already reaches r in
	/// [ 0.81, 1.66 ] and z in [ -0.87, 0.90 ]. Lifted from
	/// SurfaceAverageConvergence.cpp with its reasoning.
	Rectangle nstxBox()
	{
		return Rectangle{ 0.60, 1.90, -1.10, 1.10 };
	}

	/// One solve, kept alive, with the post-processing done so that both
	/// candidate potentials are available. Lifted from
	/// SurfaceAverageConvergence.cpp; member order is load bearing and the class
	/// is non-copyable because the coefficients capture `this`.
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

			meq::CriticalPoint axis()
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

}

/*
 * ===========================================================================
 * 1. THE ANGLE LABEL, ON A FAMILY WHOSE ANSWER IS KNOWN EXACTLY.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theAngleLabelDecidesWhetherTheBasisConvergesAtAll )
{
	// Nested similar ellipses about ( 1.3, 0 ), semi-axes ( 0.5, 1.0 ) at
	// Psi_N = 1, with discRadius = sqrt( Psi_N ). Under the ellipse's own
	// parameter this family IS a polynomial of degree one in ( x, y ):
	// R = 1.3 + 0.5 x, z = 1.0 y. So a Zernike fit of ANY degree at or above one
	// must reproduce it to round-off, and there is no approximation theory in
	// the answer to argue about.
	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.02;
	double const largest = 0.60;
	double const centre = 1.3;
	double const shortAxis = 0.5;
	double const longAxis = 1.0;

	std::vector<meq::SurfaceSample> circles;
	std::vector<meq::SurfaceSample> byParameter;
	std::vector<meq::SurfaceSample> byGeometricAngle;

	for ( int i = 0; i < surfaces; ++i )
	{
		double const normalisedFlux = smallest
			+ ( largest - smallest )*i/( surfaces - 1.0 );
		double const radius = std::sqrt( normalisedFlux );

		for ( int j = 0; j < angles; ++j )
		{
			double const theta = twoPi*j/angles;

			circles.push_back( meq::SurfaceSample{
				normalisedFlux, theta, centre + shortAxis*radius*std::cos( theta ),
				shortAxis*radius*std::sin( theta ), 1.0 } );

			byParameter.push_back( meq::SurfaceSample{
				normalisedFlux, theta, centre + shortAxis*radius*std::cos( theta ),
				longAxis*radius*std::sin( theta ), 1.0 } );

			// The point of THIS ellipse that subtends the geometric angle theta:
			// tan( parameter ) = ( short/long ) tan( theta ).
			double const parameter = std::atan2( shortAxis*std::sin( theta ),
			                                     longAxis*std::cos( theta ) );
			byGeometricAngle.push_back( meq::SurfaceSample{
				normalisedFlux, theta,
				centre + shortAxis*radius*std::cos( parameter ),
				longAxis*radius*std::sin( parameter ), 1.0 } );
		}
	}

	meq::SurfaceFitOptions options;
	options.discEdge = largest;

	std::printf( "\n=== IN-3 case 1: the angular label, on exact ellipses ===\n" );
	std::printf( "    semi-axes ( %.2f, %.2f ) at Psi_N = 1, %d surfaces x %d"
	             " angles\n\n", shortAxis, longAxis, surfaces, angles );
	std::printf( "  %-4s  %-12s  %-12s  %-12s  %-12s\n", "L", "circles",
	             "ellipse par.", "geometric", "relabelled" );

	meq::AxisShape const shape = meq::axisShapeFromSamples( byGeometricAngle,
	                                                        centre, 0.0 );
	std::vector<meq::SurfaceSample> const relabelled =
		meq::relabelByAxisShape( byGeometricAngle, shape );

	double worstCircle = 0.0;
	double worstParameter = 0.0;
	double bestGeometric = 1.0;
	double worstRelabelled = 0.0;

	for ( int degree = 2; degree <= 20; degree += 2 )
	{
		meq::SurfaceFit const circleFit( degree, circles, options );
		meq::SurfaceFit const parameterFit( degree, byParameter, options );
		meq::SurfaceFit const geometricFit( degree, byGeometricAngle, options );
		meq::SurfaceFit const relabelledFit( degree, relabelled, options );

		double const circleError = worstFitError( circleFit, circles );
		double const parameterError = worstFitError( parameterFit, byParameter );
		double const geometricError = worstFitError( geometricFit,
		                                             byGeometricAngle );
		double const relabelledError = worstFitError( relabelledFit, relabelled );

		std::printf( "  %-4d  %.6e  %.6e  %.6e  %.6e\n", degree, circleError,
		             parameterError, geometricError, relabelledError );

		worstCircle = std::max( worstCircle, circleError );
		worstParameter = std::max( worstParameter, parameterError );
		bestGeometric = std::min( bestGeometric, geometricError );
		worstRelabelled = std::max( worstRelabelled, relabelledError );
	}

	std::printf( "\n  axis ellipse recovered from the innermost surface alone:"
	             " tilt %.3e rad, short/long %.6f (exact %.6f)\n",
	             shape.tilt, shape.semiAxisRatio, shortAxis/longAxis );

	// The relabelling must be the identity on a circular axis, which is the
	// cheapest check that it is the right way up rather than its own inverse.
	meq::AxisShape const round = meq::axisShapeFromSamples( circles, centre, 0.0 );
	BOOST_TEST( std::abs( round.semiAxisRatio - 1.0 ) < 1.0e-10 );
	BOOST_TEST( std::abs( meq::shapedPoloidalAngle( round, 1.0 ) - 1.0 ) < 1.0e-12 );

	// The recovered ellipse IS the family's own, because the family is exactly
	// elliptical -- so this is a check on axisShapeFromSamples() and not an
	// approximation.
	BOOST_TEST( std::abs( shape.semiAxisRatio - shortAxis/longAxis ) < 1.0e-10 );
	BOOST_TEST( std::abs( shape.tilt ) < 1.0e-10 );

	// The parametrisation that makes the family a polynomial reproduces it
	// exactly at every degree; so does the circular family, whose geometric
	// angle IS its parameter.
	BOOST_TEST( worstCircle < 1.0e-12 );
	BOOST_TEST( worstParameter < 1.0e-12 );
	BOOST_TEST( worstRelabelled < 1.0e-12 );

	// And the SAME POINTS under the geometric angle are not merely worse, they
	// do not converge: the best any degree up to twenty reaches is still ten
	// thousand times the round-off the other three sit at. If this ever falls
	// below 1e-3 the control is empty and every claim in SurfaceFit.hpp's angle
	// section needs re-measuring rather than re-wording.
	BOOST_TEST( bestGeometric > 1.0e-3 );
	BOOST_TEST( bestGeometric/std::max( worstRelabelled, 1.0e-16 ) > 1.0e6 );
}

/*
 * ===========================================================================
 * 2. THE COEFFICIENT ENVELOPE, AND THE RADIAL COORDINATE AS ITS CONTROL.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theEnvelopeDecaysGeometricallyAndTheFluxCoordinateDoesNot )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.10;
	double const largest = 0.60;
	int const degree = 20;

	std::vector<meq::SurfaceSample> const geometric = exactSamples(
		eq, axis, levelsFor( Layout::EquispacedFlux, smallest, largest, surfaces ),
		angles );
	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                       axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	meq::SurfaceFitOptions radial;
	radial.discEdge = largest;
	meq::SurfaceFitOptions flux = radial;
	flux.coordinate = meq::FitRadialCoordinate::NormalisedFlux;

	meq::SurfaceFit const radialFit( degree, samples, radial );
	meq::SurfaceFit const fluxFit( degree, samples, flux );

	std::printf( "\n=== IN-3 case 2: the coefficient envelope, nstx(), Psi_N in"
	             " [ %.2f, %.2f ] ===\n", smallest, largest );
	std::printf( "    %d surfaces x %d angles = %zu samples, %zu modes, axis"
	             " ellipse short/long %.4f\n", surfaces, angles,
	             radialFit.diagnostics().samples, radialFit.diagnostics().modes,
	             shape.semiAxisRatio );
	std::printf( "    condition numbers: discRadius %.3e, normalised flux"
	             " %.3e\n\n", radialFit.diagnostics().conditionNumber,
	             fluxFit.diagnostics().conditionNumber );
	std::printf( "  %-4s  %-14s  %-14s  %-10s\n", "l", "rho = sqrt(Psi)",
	             "Psi (CONTROL)", "ratio" );

	std::vector<double> radialEnvelope;
	std::vector<double> fluxEnvelope;

	for ( int l = 0; l <= degree; ++l )
	{
		double const here = envelopeWindow( radialFit, l );
		double const there = envelopeWindow( fluxFit, l );
		radialEnvelope.push_back( here );
		fluxEnvelope.push_back( there );
		std::printf( "  %-4d  %.6e    %.6e    %8.2f\n", l, here, there,
		             there/here );
	}

	double const radialError = worstFitError( radialFit, samples );
	double const fluxError = worstFitError( fluxFit, samples );
	std::printf( "\n  worst fit error: discRadius %.6e, normalised flux %.6e"
	             " -- a factor of %.1f\n", radialError, fluxError,
	             fluxError/radialError );

	// GEOMETRIC DECAY, asserted on the ratio of envelopes TWO degrees apart.
	// Two and not one, because envelopeWindow() takes a maximum over a pair and
	// a one-degree ratio of overlapping windows is exactly 1 whenever the same
	// coefficient is the maximum of both -- which for a family with any symmetry
	// is most of the time. Two degrees is also the natural step for this basis,
	// since l and l - 2 carry the same angular parities.
	double worstRatio = 0.0;
	for ( int l = degree/2; l <= degree; ++l )
		worstRatio = std::max( worstRatio,
		                       radialEnvelope[ l ]/radialEnvelope[ l - 2 ] );

	double worstControlRatio = 0.0;
	for ( int l = degree/2; l <= degree; ++l )
		worstControlRatio = std::max( worstControlRatio,
		                              fluxEnvelope[ l ]/fluxEnvelope[ l - 2 ] );

	double const radialDrop = radialEnvelope[ degree ]
	                          /radialEnvelope[ degree/2 ];
	double const controlDrop = fluxEnvelope[ degree ]/fluxEnvelope[ degree/2 ];

	std::printf( "  worst two-degree envelope ratio over the tail: discRadius"
	             " %.4f, normalised flux %.4f\n", worstRatio, worstControlRatio );
	std::printf( "  and over the whole tail, l = %d to %d: discRadius falls by"
	             " %.1f, the control by %.1f\n", degree/2, degree,
	             1.0/radialDrop, 1.0/controlDrop );

	BOOST_TEST( worstRatio < 0.75 );
	BOOST_TEST( radialDrop < 0.10 );

	// AND THE CONTROL DOES NOT. Algebraic decay is still decay, so the assertion
	// cannot be that the control stops falling; what makes Zernike.hpp's
	// argument a measurement is that the GAP GROWS. It is asserted three ways --
	// the control's own tail drop, the ratio at the top degree, and the growth
	// of that ratio -- because any one of them alone could be met by a control
	// that merely trails.
	BOOST_TEST( controlDrop > 5.0*radialDrop );
	BOOST_TEST( fluxEnvelope[ degree ]/radialEnvelope[ degree ] > 50.0 );
	BOOST_TEST( ( fluxEnvelope[ degree ]/radialEnvelope[ degree ] )
	            /( fluxEnvelope[ degree/2 ]/radialEnvelope[ degree/2 ] ) > 3.0 );
	BOOST_TEST( fluxError/radialError > 30.0 );
}

/*
 * ===========================================================================
 * 3. THE PARITY CONSTRAINT, AND THE AXIS.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theParityConstraintKeepsTheAxisAPoint )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.10;
	double const largest = 0.60;

	std::vector<meq::SurfaceSample> const geometric = exactSamples(
		eq, axis, levelsFor( Layout::EquispacedFlux, smallest, largest, surfaces ),
		angles );
	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                        axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	std::printf( "\n=== IN-3 case 3: the parity constraint at the axis ===\n" );
	std::printf( "    exact axis ( %.9f, %.9f )\n\n", axis.r, axis.z );
	std::printf( "  %-4s  %-8s  %-11s  %-11s  %-11s  %-11s\n", "L", "basis",
	             "modes", "cond", "fit error", "axis error" );

	double zernikeSpread = 0.0;
	double tensorSpread = 0.0;
	double zernikeAxisError = 0.0;
	double tensorAxisError = 0.0;
	double zernikeCondition = 0.0;
	double tensorCondition = 0.0;
	double zernikeFitError = 0.0;
	double tensorFitError = 0.0;
	std::vector<double> zernikeAxisHistory;

	for ( int degree = 4; degree <= 20; degree += 4 )
	{
		for ( int which = 0; which < 2; ++which )
		{
			meq::SurfaceFitOptions options;
			options.discEdge = largest;
			options.basis = which == 0 ? meq::FitBasis::Zernike
			                           : meq::FitBasis::TensorProduct;

			meq::SurfaceFit const fit( degree, samples, options );

			// The point the fit puts at discRadius = 0, swept over the angle. For
			// Zernike this is one point by construction; for the control it is a
			// curve, and its extent is the measurement.
			double minR = 1.0e30;
			double maxR = -1.0e30;
			double minZ = 1.0e30;
			double maxZ = -1.0e30;
			for ( int j = 0; j < 64; ++j )
			{
				double r = 0.0;
				double z = 0.0;
				fit.axisAtAngle( twoPi*j/64, r, z );
				minR = std::min( minR, r );
				maxR = std::max( maxR, r );
				minZ = std::min( minZ, z );
				maxZ = std::max( maxZ, z );
			}

			double const spread = std::max( maxR - minR, maxZ - minZ );
			double const axisError = std::hypot( 0.5*( minR + maxR ) - axis.r,
			                                     0.5*( minZ + maxZ ) - axis.z );
			double const error = worstFitError( fit, samples );

			std::printf( "  %-4d  %-8s  %-11zu  %.5e  %.5e  %.5e  spread %.3e\n",
			             degree, which == 0 ? "Zernike" : "tensor",
			             fit.diagnostics().modes,
			             fit.diagnostics().conditionNumber, error, axisError,
			             spread );

			if ( which == 0 )
			{
				zernikeSpread = std::max( zernikeSpread, spread );
				zernikeAxisError = axisError;
				zernikeCondition = fit.diagnostics().conditionNumber;
				zernikeFitError = error;
				zernikeAxisHistory.push_back( axisError );
			}
			else
			{
				tensorSpread = std::max( tensorSpread, spread );
				tensorAxisError = axisError;
				tensorCondition = fit.diagnostics().conditionNumber;
				tensorFitError = error;
			}
		}
	}

	std::printf( "\n  THE CONTROL FITS THE SAMPLE CLOUD BETTER AND IS USELESS:"
	             " %.3e against %.3e on the data,\n  at a condition number of"
	             " %.2e against %.2e, and it puts the axis on a curve %.2e"
	             " wide.\n", tensorFitError, zernikeFitError, tensorCondition,
	             zernikeCondition, tensorSpread );

	// EXACT, not a tolerance: every Zernike mode with m != 0 contains only
	// powers discRadius^|m| and above, so all of them vanish at the centre and
	// what is left has no theta in it at all. Round-off in the summation is the
	// only reason this is not identically zero.
	BOOST_TEST( zernikeSpread < 1.0e-14 );

	// The control's does not, because it admits cos( theta ) with no radial
	// factor -- a mode with no limit at the origin.
	BOOST_TEST( tensorSpread > 1.0e-4 );

	// The axis is an EXTRAPOLATION into the hole the sample set leaves, so it is
	// not free; what is free is that it is a point. It improves with degree and
	// lands far closer than the control's centroid does.
	BOOST_TEST( zernikeAxisError < 1.0e-5 );
	BOOST_TEST( tensorAxisError/zernikeAxisError > 30.0 );
	BOOST_TEST( zernikeAxisHistory.back() < 0.05*zernikeAxisHistory.front() );

	// And the conditioning, which is the other half of what parity buys.
	BOOST_TEST( tensorCondition/zernikeCondition > 1.0e6 );
}

/*
 * ===========================================================================
 * 4. THE DERIVATIVE. THE ONE THAT MATTERS.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theDerivativeAgreesWithIndependentlyTracedSurfaces )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.10;
	double const largest = 0.60;
	double const at = 0.35;
	double const step = 0.02;

	std::vector<meq::SurfaceSample> const geometric = exactSamples(
		eq, axis, levelsFor( Layout::EquispacedFlux, smallest, largest, surfaces ),
		angles );
	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                        axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	meq::SurfaceFitOptions options;
	options.discEdge = largest;

	double const centreRadius = std::sqrt( at/largest );

	// d( r, z )/d( discRadius ) at a geometric angle, from surfaces solved from
	// scratch at four neighbouring levels. The fit is not consulted anywhere in
	// here, which is what makes this an independent instrument.
	auto difference = [ & ]( double theta, double width, double &dR, double &dZ )
	{
		double rPlus = 0.0;
		double zPlus = 0.0;
		double rMinus = 0.0;
		double zMinus = 0.0;
		double const outer = ( centreRadius + width )*( centreRadius + width )
		                     *largest;
		double const inner = ( centreRadius - width )*( centreRadius - width )
		                     *largest;
		rayPoint( eq, axis, outer, theta, rPlus, zPlus );
		rayPoint( eq, axis, inner, theta, rMinus, zMinus );
		dR = ( rPlus - rMinus )/( 2.0*width );
		dZ = ( zPlus - zMinus )/( 2.0*width );
	};

	std::printf( "\n=== IN-3 case 4: d( r, z )/d( discRadius ) at Psi_N = %.2f"
	             " ===\n", at );
	std::printf( "    against surfaces traced independently at four"
	             " neighbouring levels, Richardson step %.3f in discRadius\n\n",
	             step );
	std::printf( "  %-4s  %-13s  %-13s  %-13s\n", "L", "fit error",
	             "d/d(discRad)", "relative" );

	std::vector<double> richardson;

	for ( int degree = 4; degree <= 20; degree += 4 )
	{
		meq::SurfaceFit const fit( degree, samples, options );

		double worstRichardson = 0.0;
		double worstRelative = 0.0;

		for ( int j = 0; j < 16; ++j )
		{
			double const theta = twoPi*j/16;
			double const label = meq::shapedPoloidalAngle( shape, theta );

			double fitR = 0.0;
			double fitZ = 0.0;
			fit.radialDerivative( at, label, fitR, fitZ );

			double coarseR = 0.0;
			double coarseZ = 0.0;
			double fineR = 0.0;
			double fineZ = 0.0;
			difference( theta, step, coarseR, coarseZ );
			difference( theta, 0.5*step, fineR, fineZ );

			double const exactR = ( 4.0*fineR - coarseR )/3.0;
			double const exactZ = ( 4.0*fineZ - coarseZ )/3.0;
			double const magnitude = std::hypot( exactR, exactZ );
			double const error = std::hypot( fitR - exactR, fitZ - exactZ );

			worstRichardson = std::max( worstRichardson, error );
			worstRelative = std::max( worstRelative, error/magnitude );
		}

		std::printf( "  %-4d  %.6e  %.6e  %.6e\n", degree,
		             worstFitError( fit, samples ), worstRichardson,
		             worstRelative );

		richardson.push_back( worstRichardson );
	}

	// The derivative of the fit converges. It is not asserted to reach the value
	// error -- differentiating a truncated expansion costs about a factor of the
	// degree, and it does -- but it must fall by orders across the sweep.
	std::printf( "\n  the derivative error falls by %.1f from L = 4 to L = 20\n",
	             richardson.front()/richardson.back() );
	BOOST_TEST( richardson.front()/richardson.back() > 100.0 );
	BOOST_TEST( richardson.back() < 1.0e-3 );

	// ===================================================================
	// AND THE INSTRUMENT, SWEPT, WHICH IS THE FOURTH TIME THIS TREE HAS HAD TO
	// MEASURE IT.
	//
	// The comparison above is only a statement about the FIT if the difference
	// it is compared against is better than the fit is. A plain central
	// difference carries its own O( h^2 ) truncation, so it is not: swept in the
	// step it converges at two toward the fit's own error and stops there, while
	// the extrapolated one is already there at the largest step and does not
	// move. THE FLAT COLUMN IS THE ANSWER AND THE CONVERGING ONE IS THE
	// INSTRUMENT, which is exactly the reading IN-2 records for the same
	// measurement on its averaged-equation identity.
	// ===================================================================
	meq::SurfaceFit const finest( 20, samples, options );

	std::printf( "\n  the instrument, at L = 20 where the fit's own derivative"
	             " error is %.6e:\n", richardson.back() );
	std::printf( "  %-9s  %-13s  %-13s  %-9s\n", "step", "plain (CTRL)",
	             "Richardson", "ratio" );

	double worstSeparation = 0.0;
	double plainAtCoarsest = 0.0;
	double plainAtFinest = 0.0;
	double richardsonSpread = 0.0;
	double richardsonFirst = 0.0;

	for ( double width : { 0.16, 0.08, 0.04, 0.02 } )
	{
		double worstPlain = 0.0;
		double worstRichardson = 0.0;

		for ( int j = 0; j < 16; ++j )
		{
			double const theta = twoPi*j/16;
			double const label = meq::shapedPoloidalAngle( shape, theta );

			double fitR = 0.0;
			double fitZ = 0.0;
			finest.radialDerivative( at, label, fitR, fitZ );

			double coarseR = 0.0;
			double coarseZ = 0.0;
			double fineR = 0.0;
			double fineZ = 0.0;
			difference( theta, width, coarseR, coarseZ );
			difference( theta, 0.5*width, fineR, fineZ );

			worstPlain = std::max( worstPlain,
			                       std::hypot( fitR - coarseR, fitZ - coarseZ ) );
			worstRichardson = std::max(
				worstRichardson,
				std::hypot( fitR - ( 4.0*fineR - coarseR )/3.0,
				            fitZ - ( 4.0*fineZ - coarseZ )/3.0 ) );
		}

		std::printf( "  %-9.3f  %.6e  %.6e  %8.1f\n", width, worstPlain,
		             worstRichardson, worstPlain/worstRichardson );

		worstSeparation = std::max( worstSeparation,
		                            worstPlain/worstRichardson );
		if ( plainAtCoarsest == 0.0 )
		{
			plainAtCoarsest = worstPlain;
			richardsonFirst = worstRichardson;
		}
		plainAtFinest = worstPlain;
		richardsonSpread = std::max( richardsonSpread,
		                             std::abs( worstRichardson/richardsonFirst
		                                       - 1.0 ) );
	}

	std::printf( "\n  the plain column falls by %.1f across an eightfold"
	             " refinement of the step and the extrapolated one moves by"
	             " %.1f%%\n", plainAtCoarsest/plainAtFinest,
	             100.0*richardsonSpread );

	// The plain difference is dominated by its own truncation at the coarse end
	// -- if it ever were not, every derivative comparison in this suite would
	// need re-reading -- and the extrapolated one is flat, because it has
	// already reached the fit's error and is measuring that rather than itself.
	BOOST_TEST( worstSeparation > 10.0 );
	BOOST_TEST( plainAtCoarsest/plainAtFinest > 10.0 );
	BOOST_TEST( richardsonSpread < 0.25 );
}

/*
 * ===========================================================================
 * 5. THE 1/( 2 discRadius ) CONVERSION, AND HOW IT GROWS.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theFluxDerivativeConversionGrowsLikeTheInverseRootFlux )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.10;
	double const largest = 0.90;
	int const degree = 20;

	std::vector<meq::SurfaceSample> const geometric = exactSamples(
		eq, axis, levelsFor( Layout::EquispacedFlux, smallest, largest, surfaces ),
		angles );
	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                        axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	meq::SurfaceFitOptions options;
	options.discEdge = largest;
	meq::SurfaceFit const fit( degree, samples, options );

	// FIRST: the conversion is a CHAIN RULE and not a modelling choice, so on a
	// range bounded away from the axis it must agree with a difference in Psi_N
	// to the same accuracy the discRadius derivative does. Psi_N in [ 0.1, 0.9 ]
	// per INVERSION-PLAN.md IN-3.
	std::printf( "\n=== IN-3 case 5: d/dPsi_N, and what it costs near the axis"
	             " ===\n\n" );
	std::printf( "  the chain rule, on Psi_N in [ 0.1, 0.9 ]:\n" );
	std::printf( "  %-9s  %-13s  %-13s\n", "Psi_N", "worst rel err", "|d/dPsi|" );

	double worstChain = 0.0;

	for ( double at : { 0.10, 0.30, 0.50, 0.70, 0.90 } )
	{
		double worstHere = 0.0;
		double magnitude = 0.0;

		for ( int j = 0; j < 16; ++j )
		{
			double const theta = twoPi*j/16;
			double const label = meq::shapedPoloidalAngle( shape, theta );

			double fitR = 0.0;
			double fitZ = 0.0;
			fit.fluxDerivative( at, label, fitR, fitZ );

			// A Richardson difference of the FIT in Psi_N. This checks the
			// conversion factor alone -- the surfaces are not re-solved, because
			// what is under test here is 1/( 2 discRadius ) and not the fit.
			auto quotient = [ & ]( double width, double &dR, double &dZ )
			{
				double rPlus = 0.0;
				double zPlus = 0.0;
				double rMinus = 0.0;
				double zMinus = 0.0;
				fit.position( at + width, label, rPlus, zPlus );
				fit.position( at - width, label, rMinus, zMinus );
				dR = ( rPlus - rMinus )/( 2.0*width );
				dZ = ( zPlus - zMinus )/( 2.0*width );
			};

			double coarseR = 0.0;
			double coarseZ = 0.0;
			double fineR = 0.0;
			double fineZ = 0.0;
			// A small step and not a convenient one: R varies like sqrt( Psi_N )
			// so its Psi-derivatives grow toward the axis, and at Psi_N = 0.10 a
			// step of 0.01 leaves the extrapolation's own truncation at 1.6e-06
			// -- which would be read as a fault in the conversion factor. This
			// is the SAME reading as the instrument sweep in case 4, met from
			// the other end.
			quotient( 0.004, coarseR, coarseZ );
			quotient( 0.002, fineR, fineZ );
			double const exactR = ( 4.0*fineR - coarseR )/3.0;
			double const exactZ = ( 4.0*fineZ - coarseZ )/3.0;

			double const size = std::hypot( exactR, exactZ );
			worstHere = std::max( worstHere,
			                      std::hypot( fitR - exactR, fitZ - exactZ )/size );
			magnitude = std::max( magnitude, std::hypot( fitR, fitZ ) );
		}

		std::printf( "  %-9.2f  %.6e  %.6e\n", at, worstHere, magnitude );
		worstChain = std::max( worstChain, worstHere );
	}

	BOOST_TEST( worstChain < 1.0e-6 );

	// SECOND: how the conversion grows as the innermost node is lowered, which
	// is the number a future coupling reads. INVERSION-PLAN.md section 4.4:
	// MaNTA enforces vanishing fluxes on axis and can keep its nodes off
	// Psi = 0, so what is left is a conditioning question and this is it.
	//
	// THE PRODUCT IS THE POINT, NOT THE MAGNITUDE. If the growth really is the
	// coordinate's 1/( 2 sqrt( Psi_N ) ) then 2 sqrt( Psi_N ) | d/dPsi_N | tends
	// to a constant -- and that constant is the m = 1 amplitude, which is the
	// minor radius scale. A caller putting its innermost node at Psi_1 sees a
	// geometry block of magnitude about that constant over 2 sqrt( Psi_1 ).
	std::printf( "\n  the growth of the conversion toward the axis"
	             " (Psi_N below %.2f is an EXTRAPOLATION):\n", smallest );
	std::printf( "  %-9s  %-13s  %-13s  %-13s\n", "Psi_N", "1/(2 sqrt)",
	             "max |d/dPsi|", "product" );

	std::vector<double> products;

	for ( double at : { 0.50, 0.30, 0.20, 0.10, 0.05, 0.02, 0.01, 0.005 } )
	{
		double magnitude = 0.0;
		for ( int j = 0; j < 64; ++j )
		{
			double dR = 0.0;
			double dZ = 0.0;
			fit.fluxDerivative( at, twoPi*j/64, dR, dZ );
			magnitude = std::max( magnitude, std::hypot( dR, dZ ) );
		}

		double const product = 2.0*std::sqrt( at )*magnitude;
		products.push_back( product );
		std::printf( "  %-9.3f  %13.4f  %13.4f  %13.4f\n", at,
		             1.0/( 2.0*std::sqrt( at ) ), magnitude, product );
	}

	std::printf( "\n  the product settles at %.4f, so a node at Psi_1 carries a"
	             " geometry derivative of about %.3f / sqrt( Psi_1 )\n",
	             products.back(), 0.5*products.back() );

	// The product must be bounded and settling, which is the statement that the
	// growth is the coordinate's and not the fit's. A fit whose own error blew
	// up in the hole would send this column climbing instead.
	BOOST_TEST( products.back() < 2.0*products.front() );
	BOOST_TEST( std::abs( products.back()/products[ products.size() - 2 ] - 1.0 )
	            < 0.05 );
}

/*
 * ===========================================================================
 * 6. THE METRIC: THE FIT'S OWN | dx/dtheta | AGAINST THE FIELD, POINTWISE.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theAngularMetricAgreesWithTheFieldPointwise )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const smallest = 0.10;
	double const largest = 0.60;

	std::vector<double> const levels = levelsFor( Layout::EquispacedFlux,
	                                              smallest, largest, surfaces );
	std::vector<meq::SurfaceSample> const geometric = exactSamples( eq, axis,
	                                                                levels,
	                                                                angles );
	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                        axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	meq::SurfaceFitOptions options;
	options.discEdge = largest;

	std::printf( "\n=== IN-3 case 6: | dx/dtheta | from the fit against the"
	             " field's own pointwise identity ===\n" );
	std::printf( "    the field route is IN-1's sqrt( rho'^2 + rho^2 ) with"
	             " rho' = -rho ( grad psi . u' )/( grad psi . u ),\n"
	             "    which shares no arithmetic with a Zernike expansion\n\n" );
	std::printf( "  %-4s  %-13s  %-13s\n", "L", "fit error", "worst rel metric" );

	std::vector<double> worstMetric;

	for ( int degree = 4; degree <= 20; degree += 4 )
	{
		meq::SurfaceFit const fit( degree, samples, options );
		double worst = 0.0;

		for ( std::size_t i = 0; i < levels.size(); i += 6 )
		{
			meq::analytic::SurfaceQuadrature const surface =
				exactSurface( eq, axis, levels[ i ], angles );

			for ( std::size_t j = 0; j < surface.points.size(); j += 7 )
			{
				meq::analytic::SurfacePoint const &point = surface.points[ j ];
				double const label = meq::shapedPoloidalAngle( shape, point.theta );

				// The fit is parametrised by the LABEL, the field's identity by
				// the geometric angle, so the chain factor d( label )/d( theta )
				// is not optional -- and it is the factor a reader is most likely
				// to drop, because both quantities are called | dx/dtheta |.
				double const width = 1.0e-6;
				double const slope =
					( meq::shapedPoloidalAngle( shape, point.theta + width )
					  - meq::shapedPoloidalAngle( shape, point.theta - width ) )
					/( 2.0*width );

				double const metric = fit.angularSpeed( levels[ i ], label )
				                      *std::abs( slope );
				worst = std::max( worst,
				                  std::abs( metric - point.metric )/point.metric );
			}
		}

		std::printf( "  %-4d  %.6e  %.6e\n", degree,
		             worstFitError( fit, samples ), worst );
		worstMetric.push_back( worst );
	}

	BOOST_TEST( worstMetric.back() < 1.0e-3 );
	BOOST_TEST( worstMetric.front()/worstMetric.back() > 50.0 );
}

/*
 * ===========================================================================
 * 7. CONDITIONING: THE HOLE AND THE DISC EDGE, AND NOT THE LAYOUT.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theHoleAndTheDiscEdgeDecideTheConditioningAndNotTheLayout )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const axis = exactAxis( eq, 1.3, 0.0 );

	int const angles = 96;
	int const surfaces = 24;
	double const largest = 0.60;
	int const degree = 16;

	std::printf( "\n=== IN-3 case 7: conditioning ===\n" );
	std::printf( "    THE PLAUSIBLE STORY IS THAT THE LEVELS SHOULD BE"
	             " GAUSS-LEGENDRE IN Psi_N. Zernike is orthogonal on the disc\n"
	             "    under discRadius d( discRadius ) = d( Psi_N )/2, so Gauss"
	             " levels with equispaced angles ought to make the\n"
	             "    discrete inner product the continuous one and the fit a"
	             " projection. THEY DO NOT, because the orthogonality\n"
	             "    argument needs the nodes to span the whole disc and a"
	             " sample set has a HOLE in the middle -- nobody\n"
	             "    traces the surface at Psi_N = 0, because it is a"
	             " point.\n\n" );
	std::printf( "  %-9s  %-26s  %-9s  %-12s  %-12s\n", "Psi_min", "levels",
	             "discEdge", "condition", "fit error" );

	struct Reading
	{
		double condition;
		double error;
	};

	std::vector<double> holeConditions;
	double worstLayoutSpread = 0.0;
	double worstEdgeGain = 0.0;
	double worstEdgeDisagreement = 0.0;

	for ( double smallest : { 0.02, 0.05, 0.10, 0.25 } )
	{
		std::vector<Reading> rescaled;
		std::vector<Reading> partial;

		for ( Layout layout : { Layout::EquispacedFlux, Layout::EquispacedRadius,
		                        Layout::GaussFlux } )
		{
			std::vector<double> const levels = levelsFor( layout, smallest,
			                                              largest, surfaces );
			std::vector<meq::SurfaceSample> const geometric =
				exactSamples( eq, axis, levels, angles );
			meq::AxisShape const shape = meq::axisShapeFromSamples( geometric,
			                                                        axis.r,
			                                                        axis.z );
			std::vector<meq::SurfaceSample> const samples =
				meq::relabelByAxisShape( geometric, shape );

			for ( double edge : { 1.0, largest } )
			{
				meq::SurfaceFitOptions options;
				options.discEdge = edge;
				meq::SurfaceFit const fit( degree, samples, options );
				Reading const reading{ fit.diagnostics().conditionNumber,
				                       worstFitError( fit, samples ) };

				std::printf( "  %-9.2f  %-26s  %-9.2f  %.6e  %.6e\n", smallest,
				             layoutName( layout ), edge, reading.condition,
				             reading.error );

				if ( edge == largest )
					rescaled.push_back( reading );
				else
					partial.push_back( reading );
			}
		}

		std::printf( "\n" );

		// The three layouts, at the same hole and the same edge.
		double best = rescaled.front().condition;
		double worst = rescaled.front().condition;
		for ( Reading const &reading : rescaled )
		{
			best = std::min( best, reading.condition );
			worst = std::max( worst, reading.condition );
		}
		worstLayoutSpread = std::max( worstLayoutSpread, worst/best );

		// The edge, layout by layout.
		for ( std::size_t i = 0; i < rescaled.size(); ++i )
		{
			worstEdgeGain = std::max( worstEdgeGain, partial[ i ].condition
			                          /rescaled[ i ].condition );
			worstEdgeDisagreement = std::max(
				worstEdgeDisagreement,
				std::abs( partial[ i ].error - rescaled[ i ].error )
				/std::max( partial[ i ].error, rescaled[ i ].error ) );
		}

		holeConditions.push_back( rescaled.front().condition );
	}

	std::printf( "  THE HOLE, at the rescaled edge and equispaced levels:"
	             " %.3e, %.3e, %.3e, %.3e as the innermost surface\n"
	             "  is pulled out from Psi_N = 0.02 to 0.25 -- %.0f orders for a"
	             " hole of %.2f to %.2f of the radius.\n",
	             holeConditions[ 0 ], holeConditions[ 1 ], holeConditions[ 2 ],
	             holeConditions[ 3 ],
	             std::log10( holeConditions[ 3 ]/holeConditions[ 0 ] ),
	             std::sqrt( 0.02/largest ), std::sqrt( 0.25/largest ) );
	std::printf( "  THE EDGE: leaving it at Psi_N = 1 costs up to %.0fx.\n"
	             "  THE LAYOUT: the three differ by at most %.2fx, and"
	             " Gauss-Legendre is sometimes the WORST of them.\n",
	             worstEdgeGain, worstLayoutSpread );

	// THE HOLE IS THE LEVER, by four orders.
	BOOST_TEST( holeConditions[ 3 ]/holeConditions[ 0 ] > 1.0e3 );

	// SO IS THE EDGE, and the two extents fit the SAME FUNCTION -- a Zernike
	// expansion of degree L spans the polynomials of degree L in ( x, y ), which
	// is closed under scaling ( x, y ). Asserting the agreement is what says the
	// rescaling is a change of coordinate and not a change of model, and it is
	// what makes the conditioning gain free rather than a trade.
	BOOST_TEST( worstEdgeGain > 100.0 );
	BOOST_TEST( worstEdgeDisagreement < 1.0e-6 );

	// AND THE LAYOUT IS NOT. This is the assertion that keeps the wrong story
	// dead: if a layout ever moves the conditioning by an order, the paragraph
	// in SurfaceFit.hpp saying it does not is stale and needs re-measuring.
	BOOST_TEST( worstLayoutSpread < 2.0 );
}

/*
 * ===========================================================================
 * 8. THE DISCRETE LEG: THE SAME FIT ON A SOLVED FIELD, THROUGH THE TRACER.
 * ===========================================================================
 */
BOOST_AUTO_TEST_CASE( theFitReachesTheSolvedFieldThroughTheTracer )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const exact = exactAxis( eq, 1.3, 0.0 );

	int const order = 2;
	int const n = 48;
	int const angles = 64;
	int const surfaces = 12;
	double const smallest = 0.15;
	double const largest = 0.50;
	int const degree = 12;

	SolvedEquilibrium solved( eq, nstxBox(), order, n );
	meq::CriticalPoint const axis = solved.axis();
	meq::ContourTracer tracer( solved.theSolver() );

	std::printf( "\n=== IN-3 case 8: the same fit on a SOLVED field, k = %d,"
	             " n = %d ===\n", order, n );
	std::printf( "    tracer axis ( %.9f, %.9f ) against the closed form's"
	             " ( %.9f, %.9f )\n\n", axis.r, axis.z, exact.r, exact.z );

	std::vector<meq::SurfaceSample> geometric;
	std::vector<std::vector<double>> fieldSpeeds;
	int stalledRays = 0;
	double worstRayResidual = 0.0;

	for ( int i = 0; i < surfaces; ++i )
	{
		double const normalisedFlux = smallest
			+ ( largest - smallest )*i/( surfaces - 1.0 );
		double const level = levelAt( exact, normalisedFlux );

		meq::Contour const contour = tracer.traceFromAxis( level, axis );
		BOOST_TEST_REQUIRE( contour.closed() );

		// AT THE TRACER'S OWN DEFAULT TOLERANCE, 1e-12, AND WITH NO LADDER.
		// SurfaceAverageConvergence.cpp carries a loop that loosens the tolerance
		// decade by decade until fitByAngle() stops throwing, because the ray
		// Newton used to demand a residual that on a discontinuous field is
		// sometimes UNATTAINABLE -- a ray crossing a face where the level falls
		// inside the DG jump has no point on it at that level at all. It now
		// keeps its best iterate and accepts, exactly as the tracer's own
		// corrector already did, and counts what it accepted.
		meq::AngleParametrisation const fit =
			tracer.fitByAngle( contour, axis, static_cast<std::size_t>( angles ) );

		stalledRays += fit.stalledRays;
		worstRayResidual = std::max( worstRayResidual, fit.worstResidual );

		for ( std::size_t j = 0; j < fit.count(); ++j )
		{
			meq::SurfaceSample sample;
			sample.normalisedFlux = normalisedFlux;
			sample.theta = twoPi*static_cast<double>( j )/fit.count();
			sample.r = fit.pointR[ j ];
			sample.z = fit.pointZ[ j ];
			geometric.push_back( sample );
		}

		// THE METRIC AGAINST THE SOLVED FLUX, which is case 6's check with q_h
		// in place of the analytic gradient -- and it costs nothing extra now
		// that AngleParametrisation carries the q it computed rather than
		// discarding it. AngleParametrisation::speed IS sqrt( rho'^2 + rho^2 )
		// built pointwise from that q, so comparing it against the fit's own
		// | dx/dtheta | is the representation against the solved field, with no
		// arithmetic shared between them.
		BOOST_TEST_REQUIRE( fit.fluxR.size() == fit.count() );
		for ( std::size_t j = 0; j < fit.count(); ++j )
		{
			double const magnitude = std::hypot( fit.fluxR[ j ], fit.fluxZ[ j ] );
			BOOST_TEST_REQUIRE( magnitude > 0.0 );
		}

		fieldSpeeds.push_back( fit.speed );
	}

	std::printf( "  %d surfaces traced and fitted at the default tolerance:"
	             " %d stalled rays, worst ray residual %.3e\n",
	             surfaces, stalledRays, worstRayResidual );

	meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
	                                                        axis.z );
	std::vector<meq::SurfaceSample> const samples =
		meq::relabelByAxisShape( geometric, shape );

	meq::SurfaceFitOptions options;
	options.discEdge = largest;
	meq::SurfaceFit const fit( degree, samples, options );

	// The fit against the EXACT surfaces, which is what says the whole chain --
	// solve, root, trace, ray fit, least squares -- lands on the right geometry
	// rather than merely on a self-consistent one.
	double worstAgainstExact = 0.0;
	for ( int i = 0; i < surfaces; ++i )
	{
		double const normalisedFlux = smallest
			+ ( largest - smallest )*i/( surfaces - 1.0 );

		for ( int j = 0; j < 32; ++j )
		{
			double const theta = twoPi*j/32;
			double exactR = 0.0;
			double exactZ = 0.0;
			rayPoint( eq, exact, normalisedFlux, theta, exactR, exactZ );

			double r = 0.0;
			double z = 0.0;
			fit.position( normalisedFlux,
			              meq::shapedPoloidalAngle( shape, theta ), r, z );
			worstAgainstExact = std::max( worstAgainstExact,
			                              std::hypot( r - exactR, z - exactZ ) );
		}
	}

	// The metric, against the solved flux, at the fit's own nodes.
	double worstMetric = 0.0;
	for ( int i = 0; i < surfaces; ++i )
	{
		double const normalisedFlux = smallest
			+ ( largest - smallest )*i/( surfaces - 1.0 );

		for ( std::size_t j = 0; j < fieldSpeeds[ i ].size(); ++j )
		{
			double const theta = twoPi*static_cast<double>( j )
			                     /fieldSpeeds[ i ].size();
			double const width = 1.0e-6;
			double const slope =
				( meq::shapedPoloidalAngle( shape, theta + width )
				  - meq::shapedPoloidalAngle( shape, theta - width ) )
				/( 2.0*width );

			double const metric =
				fit.angularSpeed( normalisedFlux,
				                  meq::shapedPoloidalAngle( shape, theta ) )
				*std::abs( slope );
			worstMetric = std::max( worstMetric,
			                        std::abs( metric - fieldSpeeds[ i ][ j ] )
			                        /fieldSpeeds[ i ][ j ] );
		}
	}

	double const axisSpread = [ & ]
	{
		double minR = 1.0e30;
		double maxR = -1.0e30;
		double minZ = 1.0e30;
		double maxZ = -1.0e30;
		for ( int j = 0; j < 64; ++j )
		{
			double r = 0.0;
			double z = 0.0;
			fit.axisAtAngle( twoPi*j/64, r, z );
			minR = std::min( minR, r );
			maxR = std::max( maxR, r );
			minZ = std::min( minZ, z );
			maxZ = std::max( maxZ, z );
		}
		return std::max( maxR - minR, maxZ - minZ );
	}();

	double fitAxisR = 0.0;
	double fitAxisZ = 0.0;
	fit.axis( fitAxisR, fitAxisZ );

	std::printf( "  L = %d over %zu samples: condition %.3e, fit residual"
	             " %.3e, worst distance from the EXACT surfaces %.3e\n",
	             degree, fit.diagnostics().samples,
	             fit.diagnostics().conditionNumber,
	             worstFitError( fit, samples ), worstAgainstExact );
	std::printf( "  the fit's own axis ( %.9f, %.9f ) is %.3e from the closed"
	             " form's, and is theta-independent to %.3e\n", fitAxisR,
	             fitAxisZ, std::hypot( fitAxisR - exact.r, fitAxisZ - exact.z ),
	             axisSpread );
	std::printf( "  | dx/dtheta | from the fit against the SOLVED flux at every"
	             " node: worst relative %.3e\n", worstMetric );

	// THE FIX: the ray fit must not throw at the tracer's own tolerance. That is
	// the property, and it is what SurfaceAverageConvergence.cpp's ladder was
	// working around; if this ever throws again, that ladder is the record of
	// what it used to cost.
	BOOST_TEST( worstRayResidual < 1.0e-9 );

	// The chain lands on the right geometry. The bound is the discretisation's,
	// not the representation's: at k = 2 on n = 48 the traced contour is
	// O( h^(k+2) ) from the true one and the fit cannot be closer than the
	// points it was given.
	BOOST_TEST( worstAgainstExact < 5.0e-4 );

	// And the axis is still one point, on a field that has never heard of the
	// parity constraint.
	BOOST_TEST( axisSpread < 1.0e-14 );

	// The metric agrees with the solved flux to about what the fit itself is
	// worth on this mesh and degree. It is a LOOSER bound than case 6's, and
	// deliberately: there q comes from a closed form and here it is q_h, so the
	// discretisation is in the comparison as well as the truncation.
	BOOST_TEST( worstMetric < 5.0e-2 );
}

/*
 * ===========================================================================
 * 9. THE RAY FIT MUST ACCEPT WHERE THE CORRECTOR WOULD, RATHER THAN THROWING.
 * ===========================================================================
 *
 * NOT A TEST OF THE FIT. It is here because IN-3 could not have been built
 * without it: ContourTracer::fitByAngle() is the only route from a solved field
 * to a meq::SurfaceSample, and on a discontinuous field it used to THROW rather
 * than answer.
 *
 * The mechanism is the one FluxSurfaces.hpp already documents for the tracer's
 * own corrector. { psi_h = c } is a union of per-element arcs offset by the DG
 * jump, so a RAY crossing a face where c falls inside that jump has no point on
 * it with psi_h = c AT ALL, and no tolerance tighter than the jump is
 * attainable there. The corrector keeps its best iterate and accepts after four
 * non-improving steps; the ray Newton demanded the tolerance and threw.
 * Measured before the repair, at k = 1 on the raw pairing it failed at 1e-12 on
 * every mesh from n = 12 to 32, and the probability rises with the angle count
 * because more rays are more chances -- which is why this case asks for a lot of
 * them. tests/convergence/SurfaceAverageConvergence.cpp carries the tolerance
 * ladder that was written to work around it, and that ladder is now the record
 * of what it used to cost.
 *
 * IT ASSERTS THAT THE MECHANISM FIRES, not merely that nothing threw. A run
 * reading zero stalled rays would pass a no-throw assertion while telling
 * nothing whatever about the repair -- so if this ever reads zero, the case has
 * stopped exercising what it is for and needs a coarser mesh or more rays
 * rather than a quiet pass.
 */
BOOST_AUTO_TEST_CASE( theRayFitAcceptsWhereTheCorrectorWouldRatherThanThrowing )
{
	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const exact = exactAxis( eq, 1.3, 0.0 );

	int const order = 1;
	int const n = 16;
	std::size_t const angles = 512;

	SolvedEquilibrium solved( eq, nstxBox(), order, n );
	meq::CriticalPoint const axis = solved.axis();

	// THE RAW PAIRING, deliberately: psi_h jumps at O( h^(k+1) ) across a face
	// where psi* jumps at O( h^(k+2) ), so the raw pairing is where the
	// unattainable tolerance actually lives. Potential::PostProcessed is the
	// tracer's default and meets 1e-12 on this mesh, which is why it would make
	// a silent and useless version of this test.
	meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );

	std::printf( "\n=== IN-3 case 9: the ray fit on a field whose own tolerance"
	             " is unattainable ===\n" );
	std::printf( "    k = %d, n = %d, raw pairing, %zu rays, tracer default"
	             " tolerance 1e-12\n\n", order, n, angles );
	std::printf( "  %-8s  %-8s  %-9s  %-13s  %-13s\n", "Psi_N", "points",
	             "stalled", "worst resid", "target" );

	int totalStalled = 0;
	double worstResidual = 0.0;
	double worstJump = 0.0;

	for ( double normalisedFlux : { 0.20, 0.30, 0.40 } )
	{
		double const level = levelAt( exact, normalisedFlux );
		meq::Contour const contour = tracer.traceFromAxis( level, axis );
		BOOST_TEST_REQUIRE( contour.closed() );

		meq::AngleParametrisation const fit = tracer.fitByAngle( contour, axis,
		                                                         angles );

		std::printf( "  %-8.2f  %-8zu  %-9d  %.6e  %.6e  jump %.6e\n",
		             normalisedFlux, contour.points.size(), fit.stalledRays,
		             fit.worstResidual, contour.correctorTarget,
		             contour.worstFaceJump );

		totalStalled += fit.stalledRays;
		worstResidual = std::max( worstResidual, fit.worstResidual );
		worstJump = std::max( worstJump, contour.worstFaceJump );

		// The flux is carried on the fit now rather than recomputed by the
		// consumer, and it is the SAME q the identity below was built from.
		BOOST_TEST_REQUIRE( fit.fluxR.size() == fit.count() );
		for ( std::size_t j = 0; j < fit.count(); ++j )
		{
			double const magnitude = std::hypot( fit.fluxR[ j ], fit.fluxZ[ j ] );
			BOOST_TEST_REQUIRE( magnitude > 0.0 );
		}
	}

	std::printf( "\n  %d rays of %zu accepted at their best reachable residual"
	             " rather than at the tolerance\n", totalStalled, 3*angles );

	// THE MECHANISM FIRED. See the header of this case: a zero here is not a
	// pass, it is the case having stopped testing the thing it exists for.
	BOOST_TEST( totalStalled > 0 );

	// AND WHAT IT ACCEPTED IS AS CLOSE AS THE FIELD ALLOWS, WHICH IS A
	// SELF-CALIBRATING STATEMENT AND NOT A CHOSEN NUMBER. A stalled ray is one
	// that could not beat the DG jump of psi_h at the face it crossed, so its
	// residual is bounded by that jump -- and the tracer MEASURES the jump
	// independently, at every face the CONTOUR crosses, and reports it as
	// Contour::worstFaceJump. Comparing the two is comparing a ray's failure
	// against the field's own ambiguity, so it stays right at every mesh and
	// every degree, where a fixed tolerance would have to be re-tuned at each.
	// If this ever fails, the acceptance has stopped being "as close as the
	// field allows" and has become "wherever it happened to stop".
	std::printf( "  the worst residual anywhere is %.3e, against a measured DG"
	             " jump of %.3e -- a ratio of %.2f\n", worstResidual, worstJump,
	             worstResidual/worstJump );
	BOOST_TEST( worstResidual < 5.0*worstJump );
}
