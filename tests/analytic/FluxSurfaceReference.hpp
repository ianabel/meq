#ifndef MEQ_ANALYTIC_FLUXSURFACEREFERENCE_HPP
#define MEQ_ANALYTIC_FLUXSURFACEREFERENCE_HPP

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

/*
 * Flux-surface averages on an EXACT equilibrium, to converged precision.
 *
 * This is the reference INVERSION-PLAN.md IN-2 measures against, and it is a
 * reference VALUE rather than a closed form. There is no closed form: psi is
 * elementary but V'( psi ), < R^-2 > and the safety factor are integrals over a
 * CONTOUR of it, and the contours of the Cerfon-Freidberg twelve-coefficient
 * family are not curves whose arc length is elementary. What is available
 * instead is a quadrature that converges geometrically, which on a smooth
 * closed surface reaches round-off and is therefore as good as a closed form
 * for the purpose of testing something else against it.
 *
 * IT NEVER TOUCHES psi_h, AND THAT IS THE ENTIRE REASON IT EXISTS. Every
 * quantity here is computed from the analytic psi and grad psi of a fixture in
 * this directory. So a comparison against it measures INVERSION-PLAN.md
 * section 2's error (a) -- the discretisation -- with (b) and (c) removed by
 * construction, rather than measuring all three together and attributing the
 * total to whichever one is being worked on.
 *
 * ONE FACILITY, NOT A LIST OF QUANTITIES.
 *
 * average() takes a callable, so V', < R^-2 >, < |grad psi|^2 / R^2 > and
 * anything else are one line each rather than one function each. That is
 * deliberate and it is a response to what the consumer looks like:
 * MANTA-COUPLING.md section 5 says the geometry slots are "a set of scalar
 * functions of the flux label" whose concrete membership "is negotiated with
 * the transport physics case, not fixed by MaNTA" -- a list that is still
 * moving must not be baked in as a list of functions. The same shape is what
 * IN-2 should build in src/meq.
 *
 * THE RADIUS IS SOLVED, NOT MARCHED, AND rho' IS NOT DIFFERENCED.
 *
 * Each surface point is found by a safeguarded root solve along a ray from the
 * magnetic axis: psi( a + rho u( theta ) ) = c, with the derivative along the
 * ray being grad psi . u exactly. So the point sits on the level set to solver
 * tolerance and carries no marching error.
 *
 * The METRIC is the part that is easy to get wrong and is the whole subject of
 * INVERSION-PLAN.md section 3.2's block quote: it is easy to build a
 * spectrally accurate rule and then feed it a SECOND-ORDER Jacobian obtained by
 * differencing neighbouring node positions, at which point the whole scheme is
 * second order and nothing in its output says so. Differentiating the defining
 * identity psi( a + rho( theta ) u( theta ) ) = c gives
 *
 *     grad psi . ( rho' u + rho u' ) = 0    so    rho' = -rho ( grad psi . u' )
 *                                                       / ( grad psi . u ),
 *
 * which is POINTWISE from the gradient with nothing differenced, and then
 * | dx/dtheta | = sqrt( rho'^2 + rho^2 ). differencedMetric() computes the same
 * quantity the wrong way and is kept as the control, because "spectral" means
 * nothing without one.
 *
 * STAR-SHAPEDNESS IS A HYPOTHESIS AND IS REPORTED, NEVER ASSUMED.
 *
 * The denominator above is grad psi . u, which vanishes exactly when the ray
 * from the axis is TANGENT to the surface -- that is, when the surface is not
 * star-shaped about the axis. SurfaceQuadrature::transversality is
 * min | grad psi . u | / | grad psi | over the walk and is the direct analogue
 * of IndexAudit::transversality in src/meq/CriticalPoints.hpp, kept for the
 * same reason: a hypothesis that is measured is a hypothesis that cannot be
 * silently violated. INVERSION-PLAN.md section 3.4 records that ray methods
 * fail on indented cross-sections, at and beyond the separatrix, and near the
 * axis where the bracket degenerates; this is that failure made visible.
 *
 * WHY A RAY IS ACCEPTABLE HERE AND NOT AS MEQ'S PRIMARY ROUTE. Section 3.4
 * rejects ray bisection as the tracer because it assumes star-shapedness. It is
 * the right choice for a REFERENCE on an analytic equilibrium, where the
 * surfaces are known to be star-shaped and the hypothesis is checked every
 * call, and its independence from the predictor-corrector route is exactly what
 * makes agreement between the two worth something.
 *
 * ---------------------------------------------------------------------------
 * WHAT IT MEASURES, ON SolovievEquilibrium::nstx(), AXIS AT
 * ( 1.318167937714, 0.011088725858 ), psi_ax = -2.662896051834e-01
 * ---------------------------------------------------------------------------
 *
 * The axis is an interior MINIMUM here -- F is single-signed negative on this
 * fixture, so psi is a subsolution -- and the levels below are
 * psi_ax ( 1 - Psi_N ) with psi = 0 on the separatrix.
 *
 *   Psi_N   V'                        < R^-2 >                  arc length
 *   0.25    6.71413847786385e+01      6.84662650462071e-01      2.86563934
 *   0.50    7.42523454479159e+01      8.71200178854720e-01      4.25101127
 *   0.75    8.81259886962482e+01      1.29438587807901e+00      5.57924227
 *
 * Every one of those is converged to the last printed digit by 256 angles, and
 * transversality reads 0.75, 0.73 and 0.70 -- comfortably star-shaped -- with a
 * worst root-solve residual of 1.1e-15.
 *
 * THE ARC LENGTH IS WHERE THE POINT IS MADE. Successive differences at
 * Psi_N = 0.25, doubling the number of angles from 16:
 *
 *   from q, pointwise      1.08e-02   4.90e-05   2.98e-08   3.11e-15   1.78e-15
 *   differenced positions  1.37e-01   4.07e-02   1.07e-02   2.73e-03   6.84e-04
 *   chord sum              4.56e-02   1.07e-02   2.73e-03   6.84e-04   1.71e-04
 *
 * The first column is GEOMETRIC and reaches round-off at 256 angles. The other
 * two fall by a factor of four per doubling, which is O( N^-2 ), and they
 * converge to the SAME limit -- so nothing about their output says they are
 * three orders of magnitude worse. At 256 angles the pointwise answer is exact
 * and the differenced one is wrong in the fourth decimal.
 *
 * THAT MIDDLE ROW IS THE WHOLE OF SECTION 3.2's WARNING IN ONE LINE: the RULE
 * is spectral in both of the first two rows -- it is the same trapezoid over the
 * same points -- and the second is second order because its JACOBIAN is. A
 * spectrally accurate quadrature fed a differenced metric is a second-order
 * scheme wearing a spectral name.
 *
 * ONE PRACTICAL NOTE ON reach. It is a search radius and it is not optional to
 * think about: nstx is elongated, and at Psi_N = 0.75 the surface extends past
 * rho = 0.8 vertically while reaching it by rho = 0.4 outboard. A reach of 0.8
 * throws -- correctly, and with a message that says which of the two things
 * went wrong -- where 2.0 succeeds. The throw is the designed behaviour and is
 * better than the prior art's silent partial curve; see section 4.3.
 */

namespace meq
{

	namespace analytic
	{

		/// One point of the quadrature. Everything a surface average could want
		/// is here, so that an integrand is a function of this and of nothing
		/// else -- in particular no integrand needs to call back into the
		/// equilibrium and risk evaluating it at a different point.
		struct SurfacePoint
		{
			double theta = 0.0;    ///< angle about the magnetic axis
			double rho = 0.0;      ///< distance from the axis along that ray
			double rhoPrime = 0.0; ///< d rho / d theta, pointwise from grad psi
			double r = 0.0;
			double z = 0.0;
			double psiR = 0.0;     ///< d psi / d r at the point
			double psiZ = 0.0;
			double gradient = 0.0; ///< | grad psi |
			double metric = 0.0;   ///< | dx / d theta | = sqrt( rho'^2 + rho^2 )
			double residual = 0.0; ///< | psi( x ) - c |, the root solve's own error
		};

		/// A closed surface, sampled at equispaced angle about the axis.
		struct SurfaceQuadrature
		{
			std::vector<SurfacePoint> points;

			double level = 0.0;

			/// V'( psi ) = closed-integral 2 pi R dl / | grad psi |. THE 2 pi R
			/// IS PART OF THE DEFINITION and the flux-surface-averaged
			/// Grad-Shafranov identity below is stated for this one; a
			/// per-unit-length convention changes it.
			double vPrime = 0.0;

			/// min | grad psi . u | / | grad psi |: how far the surface is from
			/// failing to be star-shaped about the axis. See the header.
			double transversality = 0.0;

			/// The worst root-solve residual over the walk, in units of psi.
			double worstResidual = 0.0;

			/// < f >_psi = ( 1 / V' ) closed-integral 2 pi R f dl / | grad psi |.
			double average( std::function<double( SurfacePoint const & )> const &f ) const
			{
				double weighted = 0.0;
				double total = 0.0;

				for ( SurfacePoint const &p : points )
				{
					double const w = p.r*p.metric/p.gradient;
					weighted += w*f( p );
					total += w;
				}

				if ( !( total > 0.0 ) )
				{
					throw std::runtime_error(
						"SurfaceQuadrature::average: the surface has no weight" );
				}

				return weighted/total;
			}

			/// The same integral WITHOUT dividing by V', which is what an
			/// integrand that is already a density wants.
			double integrate( std::function<double( SurfacePoint const & )> const &f ) const
			{
				double sum = 0.0;
				double const dTheta = 2.0*M_PI/static_cast<double>( points.size() );

				for ( SurfacePoint const &p : points )
				{
					sum += 2.0*M_PI*p.r*f( p )*p.metric/p.gradient;
				}

				return sum*dTheta;
			}

			/// | dx / d theta | at point j by a central difference of the
			/// NEIGHBOURING POSITIONS -- the second-order control of section
			/// 3.2, kept so that "spectral" has something to be spectral
			/// against. Never used by average(); a caller asks for it
			/// deliberately.
			double differencedMetric( std::size_t j ) const
			{
				std::size_t const n = points.size();
				std::size_t const next = ( j + 1 )%n;
				std::size_t const previous = ( j + n - 1 )%n;
				double const dTheta = 2.0*M_PI/static_cast<double>( n );

				double const dr = ( points[ next ].r - points[ previous ].r )/( 2.0*dTheta );
				double const dz = ( points[ next ].z - points[ previous ].z )/( 2.0*dTheta );

				return std::sqrt( dr*dr + dz*dz );
			}
		};

		/**
		 * Sample the surface psi = level at @a angles equispaced angles about
		 * ( rAxis, zAxis ).
		 *
		 * @param eq        anything with psi( r, z ) and gradPsi( r, z, gr, gz ).
		 * @param reach     how far along a ray to search, in the same units as r.
		 *                  A bracket is found by marching and then closed by a
		 *                  safeguarded Newton; the march is only ever used to
		 *                  BRACKET, never to locate.
		 *
		 * @throws std::runtime_error if a ray fails to bracket the level, which
		 *         is what a non-star-shaped surface or a level outside the
		 *         plasma looks like from here. It throws rather than returning a
		 *         partial curve, per section 4.3's complaint against the prior
		 *         art.
		 */
		template<class Equilibrium>
		SurfaceQuadrature surfaceQuadrature( Equilibrium const &eq,
		                                     double rAxis, double zAxis,
		                                     double level, int angles,
		                                     double reach,
		                                     int marchSteps = 400 )
		{
			if ( angles < 8 )
			{
				throw std::invalid_argument(
					"surfaceQuadrature: need at least eight angles" );
			}

			SurfaceQuadrature surface;
			surface.level = level;
			surface.points.resize( static_cast<std::size_t>( angles ) );
			surface.transversality = 1.0;

			for ( int j = 0; j < angles; ++j )
			{
				double const theta = 2.0*M_PI*static_cast<double>( j )
				                     /static_cast<double>( angles );
				double const cosine = std::cos( theta );
				double const sine = std::sin( theta );

				auto along = [ & ]( double rho )
				{
					return eq.psi( rAxis + rho*cosine, zAxis + rho*sine ) - level;
				};

				// Bracket. The march exists only to find a sign change; every
				// digit of the answer comes from the solve below it.
				double lower = 0.0;
				double lowerValue = along( 0.0 );
				double upper = -1.0;

				for ( int step = 1; step <= marchSteps; ++step )
				{
					double const rho = reach*static_cast<double>( step )
					                   /static_cast<double>( marchSteps );
					double const value = along( rho );

					if ( lowerValue*value <= 0.0 )
					{
						upper = rho;
						break;
					}

					lower = rho;
					lowerValue = value;
				}

				if ( !( upper > 0.0 ) )
				{
					throw std::runtime_error(
						"surfaceQuadrature: a ray did not bracket the level -- the"
						" surface is not star-shaped about the given axis, or the"
						" level lies outside the plasma" );
				}

				// Safeguarded Newton: bisection keeps the bracket, Newton
				// supplies the order, and the derivative along the ray is
				// grad psi . u exactly.
				double rho = 0.5*( lower + upper );

				for ( int iteration = 0; iteration < 100; ++iteration )
				{
					double const r = rAxis + rho*cosine;
					double const z = zAxis + rho*sine;
					double gr = 0.0;
					double gz = 0.0;
					eq.gradPsi( r, z, gr, gz );

					double const value = eq.psi( r, z ) - level;
					double const slope = gr*cosine + gz*sine;

					if ( value*lowerValue > 0.0 )
					{
						lower = rho;
						lowerValue = value;
					}
					else
					{
						upper = rho;
					}

					double next = rho;

					if ( std::abs( slope ) > 0.0 )
					{
						next = rho - value/slope;
					}

					if ( !( next > lower && next < upper ) )
					{
						next = 0.5*( lower + upper );
					}

					double const move = std::abs( next - rho );
					rho = next;

					if ( move < 1.0e-15*( 1.0 + std::abs( rho ) ) )
					{
						break;
					}
				}

				SurfacePoint &p = surface.points[ static_cast<std::size_t>( j ) ];
				p.theta = theta;
				p.rho = rho;
				p.r = rAxis + rho*cosine;
				p.z = zAxis + rho*sine;
				eq.gradPsi( p.r, p.z, p.psiR, p.psiZ );
				p.gradient = std::sqrt( p.psiR*p.psiR + p.psiZ*p.psiZ );
				p.residual = std::abs( eq.psi( p.r, p.z ) - level );

				// rho' pointwise, from the gradient alone. u' = ( -sin, cos ).
				double const alongRay = p.psiR*cosine + p.psiZ*sine;
				double const acrossRay = -p.psiR*sine + p.psiZ*cosine;

				if ( !( std::abs( alongRay ) > 0.0 ) )
				{
					throw std::runtime_error(
						"surfaceQuadrature: the ray is tangent to the surface --"
						" grad psi . u vanishes, so the surface is not"
						" star-shaped about the given axis" );
				}

				p.rhoPrime = -rho*acrossRay/alongRay;
				p.metric = std::sqrt( p.rhoPrime*p.rhoPrime + rho*rho );

				surface.transversality = std::min(
					surface.transversality, std::abs( alongRay )/p.gradient );
				surface.worstResidual = std::max( surface.worstResidual, p.residual );
			}

			surface.vPrime = surface.integrate(
				[]( SurfacePoint const & ) { return 1.0; } );

			return surface;
		}

		/**
		 * The flux-surface average of the Grad-Shafranov equation itself, as a
		 * residual that a correct set of averages drives to zero.
		 *
		 *     ( 1 / V' ) d/dpsi ( V' < | grad psi |^2 / R^2 > )  =  < Delta* psi / R^2 >
		 *
		 * from < div G > = ( 1 / V' ) d/dpsi ( V' < G . grad psi > ) with
		 * G = grad psi / R^2, using div( grad psi / R^2 ) = Delta* psi / R^2.
		 * The right-hand side is -< F / R^2 >, and F is the source the solver is
		 * actually fed -- so this checks THREE averages against each other with
		 * nothing but the equation, and needs no reference value at all. It is
		 * the same discipline as SolovievEquilibrium::deltaStarFD(): an
		 * independent quantity is the only thing that catches a misread formula.
		 *
		 * THE d/dpsi MUST BE RICHARDSON-EXTRAPOLATED. A plain central difference
		 * carries its own O( step^2 ) truncation, which floors the agreement
		 * around 1e-6 however exact everything else is -- the floor is the
		 * INSTRUMENT, not the identity, and an identity checked at 1e-6 would
		 * pass with a real defect underneath it. Measured on a case where every
		 * piece is elementary, plain reads 8.1e-07 where ( 4 D( h/2 ) - D( h ) )
		 * / 3 reads 2.7e-13. This is the third place in the tree where that
		 * applies; see src/meq/Zernike.hpp's derivative test and
		 * INVERSION-PLAN.md IN-3.
		 *
		 * @returns left-hand side minus right-hand side, in the units of F / R^2.
		 */
		template<class Equilibrium>
		double averagedGradShafranovResidual( Equilibrium const &eq,
		                                      double rAxis, double zAxis,
		                                      double level, int angles,
		                                      double reach, double step )
		{
			auto weighted = [ & ]( double c )
			{
				SurfaceQuadrature const s = surfaceQuadrature(
					eq, rAxis, zAxis, c, angles, reach );
				return s.vPrime*s.average(
					[]( SurfacePoint const &p )
					{
						return ( p.psiR*p.psiR + p.psiZ*p.psiZ )/( p.r*p.r );
					} );
			};

			double const coarse = ( weighted( level + step )
			                        - weighted( level - step ) )/( 2.0*step );
			double const fine = ( weighted( level + 0.5*step )
			                      - weighted( level - 0.5*step ) )/step;
			double const derivative = ( 4.0*fine - coarse )/3.0;

			SurfaceQuadrature const here = surfaceQuadrature(
				eq, rAxis, zAxis, level, angles, reach );

			double const leftHandSide = derivative/here.vPrime;
			double const rightHandSide = here.average(
				[ & ]( SurfacePoint const &p )
				{
					// Delta* psi = -F, and F is the source the solver is fed.
					return -eq.f( p.r, p.z, level )/( p.r*p.r );
				} );

			return leftHandSide - rightHandSide;
		}

	}

}

#endif
