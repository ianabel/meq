#include "Coils.hpp"

#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

/*
 * WHAT IS DECIDED IN THIS FILE, AND WHY EACH WAY
 *
 * Coils.hpp carries the physics -- the derivation of F = mu0 r j_phi, the
 * cross-section integral, and what the quadrature achieves. This comment is
 * about the three numerical choices underneath it, each of which was measured
 * before it was believed.
 *
 *
 * 1. THE KERNEL IS CARLSON'S, IN THE COMPLEMENTARY MODULUS
 *
 * The textbook filament flux is
 *
 *     psi = ( mu0 I / 2 pi ) d [ ( 1 - k^2/2 ) K( k ) - E( k ) ],
 *     d^2 = ( a + r )^2 + ( z - z0 )^2,    k^2 = 4 a r / d^2,
 *
 * and it is unusable here for one reason: the quadrature has to evaluate it at
 * points arbitrarily close to the source, where k -> 1. Formed as
 * 2 sqrt( a r )/d, k rounds to exactly 1.0 once the field point is within about
 * 1e-8 of a source filament -- at which std::comp_ellint_1( 1.0 ) is NaN and
 * std::comp_ellint_2( 1.0 ) is a perfectly ordinary 1.0, so the failure is a
 * NaN propagating out of an otherwise plausible calculation.
 *
 * The complementary modulus has NO cancellation at all:
 *
 *     k'^2 = 1 - k^2 = [ ( a + r )^2 + dz^2 - 4 a r ] / d^2
 *                    = [ ( a - r )^2 + dz^2 ] / d^2,
 *
 * which is the squared distance from the field point to the source, over d^2 --
 * nothing is subtracted from one anywhere. Carlson's symmetric forms take
 * exactly that argument (DLMF 19.25.1):
 *
 *     K = R_F( 0, k'^2, 1 ),      E = K - ( k^2/3 ) R_D( 0, k'^2, 1 ),
 *
 * and the bracket collapses to k^2 [ R_D/3 - R_F/2 ], which is what
 * geometricKernel() below evaluates.
 *
 * MEASURED. Against std::comp_ellint at k = 1e-8, 1e-3, 0.1, 0.5, 0.9, 0.99,
 * 0.999999 and 1 - 1e-12, the two routes agree to a worst 7.0e-13 in K and
 * 4.3e-13 in E, both at the k -> 1 end where the standard library is itself
 * losing digits. Beyond it Carlson simply keeps going: at k'^2 = 1e-16, 1e-24,
 * 1e-40 and 1e-300 it returns 19.807, 29.017, 47.438 and 346.774, agreeing with
 * ln( 4/k' ) to every digit printed, while the textbook form has been NaN since
 * 1e-16. Boost is the implementation for the reason CLAUDE.md gives under
 * *Prefer a maintained library to a hand-rolled algorithm*.
 *
 * PROMOTION IS TURNED OFF, and that is a measurement rather than a habit. Boost
 * evaluates a double-precision special function in long double by default. Over
 * k'^2 from 1e-30 to 1, promoted and unpromoted R_F and R_D agree to a worst
 * 4.9e-16 -- round-off -- and the unpromoted pair costs 40.9 ns against 88.7,
 * so promotion is buying nothing here and charging 2.2x for it. Since the
 * kernel is evaluated 4 n^2 times per point that factor is the whole cost of
 * this class.
 *
 *
 * 2. THE GAUSS NODES ARE COMPUTED, NOT TABULATED, AND THE REASON IS THE API
 *
 * Boost's boost::math::quadrature::gauss<Real, N> takes N as a TEMPLATE
 * parameter, so a runtime-selectable order cannot use it without a switch over
 * a handful of instantiations -- which is what ExteriorDtN.cpp does, and it
 * can, because its order is fixed at 30 by an argument about polynomial degree.
 * Here the order is the knob a caller turns to trade accuracy against time, and
 * the measured table in Coils.hpp is a sweep over it, so it has to be a
 * run-time value.
 *
 * So the nodes come from Newton on Boost's own legendre_p and legendre_p_prime
 * -- the special function is still the library's, and only the root finding is
 * local, which is the same division of labour ExteriorDtN.cpp makes when it
 * uses legendre_p_prime rather than a recurrence. Checked against Boost's
 * tabulated gauss<double, 30>: abscissae and weights agree to a worst 1.5e-16
 * in ABSOLUTE value, and the rule integrates x^m over [ -1, 1 ] to a worst
 * absolute error of 4.4e-16 for every m up to 2n - 1 = 59, which is exactness
 * at the round-off floor.
 *
 *
 * 3. PANELS AND GRADING
 *
 * See Coils.hpp for what these buy. The mechanics:
 *
 *   * The rectangle is cut at the field point's coordinates CLAMPED to it, so a
 *     point inside gives four panels, a point outside in one coordinate gives
 *     two, and a point outside in both gives one. The clamp is what makes the
 *     three cases one piece of code.
 *   * Each panel is graded toward the corner nearest the field point --
 *     t = corner +- L u^3 -- which for an interior point is the singularity and
 *     for an exterior one is the nearest approach.
 *
 * The grading exponent is a compile-time constant here, not a parameter.
 * Coils.hpp records why: at exponent 4 the innermost node's offset underflows
 * against the coordinate it is added to, the node lands on the field point, and
 * the kernel is infinite. Cubic keeps a margin at every accepted order, and the
 * order is capped so that it keeps it.
 *
 * THE UNDERFLOW GUARD IS KEPT ANYWAY. A node whose distance to the field point
 * rounds to exactly zero is DROPPED rather than evaluated -- Boost would
 * otherwise throw a domain error from inside a quadrature loop, which is a
 * confusing place to meet a geometry problem. It cannot happen at the shipped
 * exponent and accepted orders; it is there so that if it ever does, the answer
 * is slightly wrong rather than absent. What is dropped is one node whose
 * quadrature weight carries a factor L u^2 of order 1e-15, so the omission is
 * far below round-off in the integral.
 */

namespace meq
{
	namespace
	{
		/// pi, to the precision the rest of this tree uses.
		constexpr double pi = 3.14159265358979323846;

		/// The grading exponent: t = corner + L u^gradingExponent. Fixed at 3
		/// for the floating-point reason in Coils.hpp, not because 3 is
		/// mathematically special -- the observed convergence rate is about
		/// 4 * gradingExponent and 4 would be better if it worked.
		constexpr int gradingExponent = 3;

		/// Boost's elliptic integrals without long-double promotion. Measured
		/// at 4.9e-16 agreement and 2.2x the speed; see the file comment.
		using EllipticPolicy = boost::math::policies::policy<
			boost::math::policies::promote_double<false> >;

		/// A Gauss-Legendre rule on [ -1, 1 ].
		struct GaussRule
		{
			std::vector<double> abscissa;
			std::vector<double> weight;
		};

		/// The rule of a given order, computed once per order per thread.
		///
		/// thread_local rather than static, because meq threads over surfaces
		/// and rays elsewhere in this tree and a shared mutable cache is
		/// exactly the kind of thing CLAUDE.md records going wrong quietly. A
		/// rule is a few hundred bytes and there will be one or two orders in
		/// play.
		GaussRule const &gaussRule( int order )
		{
			static thread_local std::map<int, GaussRule> cache;

			auto const found = cache.find( order );
			if ( found != cache.end() )
				return found->second;

			GaussRule rule;
			rule.abscissa.assign( static_cast<std::size_t>( order ), 0.0 );
			rule.weight.assign( static_cast<std::size_t>( order ), 0.0 );

			for ( int i = 0; i < order; ++i )
			{
				// The standard Chebyshev-like starting guess; Newton on
				// P_n( t ) then converges in a handful of steps for every
				// order this class accepts.
				double t = std::cos( pi*( i + 0.75 )/( order + 0.5 ) );
				for ( int step = 0; step < 100; ++step )
				{
					double const value = boost::math::legendre_p( order, t );
					double const slope =
						boost::math::legendre_p_prime( order, t );
					double const correction = -value/slope;
					t += correction;
					if ( std::abs( correction ) <= 1.0e-16 )
						break;
				}

				double const slope = boost::math::legendre_p_prime( order, t );
				rule.abscissa[ static_cast<std::size_t>( i ) ] = t;
				rule.weight[ static_cast<std::size_t>( i ) ] =
					2.0/( ( 1.0 - t*t )*slope*slope );
			}

			return cache.emplace( order, std::move( rule ) ).first->second;
		}

		/// The three geometric quantities the filament kernel is built from:
		/// d, k^2 and the COMPLEMENTARY k'^2, the last computed as the squared
		/// distance to the source over d^2 so that nothing cancels.
		void loopGeometry( double r, double z, double loopRadius,
		                   double loopHeight, double &d, double &kSquared,
		                   double &complementary )
		{
			double const dz = z - loopHeight;
			double const sum = loopRadius + r;
			double const difference = loopRadius - r;
			double const dSquared = sum*sum + dz*dz;

			d = std::sqrt( dSquared );
			kSquared = 4.0*loopRadius*r/dSquared;
			complementary = ( difference*difference + dz*dz )/dSquared;
		}

		/// psi of a filament, per unit current and per unit mu0:
		///
		///     G = ( 1/2 pi ) d k^2 [ R_D( 0, k'^2, 1 )/3 - R_F( 0, k'^2, 1 )/2 ]
		///
		/// which is identically ( 1/2 pi ) d [ ( 1 - k^2/2 ) K - E ]. The
		/// caller supplies the geometry, because it also needs k'^2 to decide
		/// whether the point is on the loop at all.
		double geometricKernel( double d, double kSquared, double complementary )
		{
			double const rf =
				boost::math::ellint_rf( 0.0, complementary, 1.0,
				                        EllipticPolicy() );
			double const rd =
				boost::math::ellint_rd( 0.0, complementary, 1.0,
				                        EllipticPolicy() );
			return d*kSquared*( rd/3.0 - rf/2.0 )/( 2.0*pi );
		}

		/// The two pieces the derivative kernel is built from, per unit
		/// current and per unit mu0:
		///
		///     A = k^2 ( R_D/3 - R_F/2 )          the psi bracket itself
		///     B = k^2 ( R_F - R_D/3 ) / k'^2     which is E/k'^2 - K
		///
		/// BOTH CARRY k^2 AS AN EXPLICIT FACTOR, which is what makes them
		/// exactly 0.0 on the axis rather than a cancellation of two numbers
		/// near pi/2. B written as E/k'^2 - K loses its figures there and in
		/// the far field alike; the identity
		///
		///     E/k'^2 - K = [ R_F - ( k^2/3 )R_D - k'^2 R_F ]/k'^2
		///                = k^2( R_F - R_D/3 )/k'^2
		///
		/// uses k'^2 = 1 - k^2 once and never subtracts anything from one.
		void kernelBrackets( double kSquared, double complementary,
		                     double &a, double &b )
		{
			double const rf =
				boost::math::ellint_rf( 0.0, complementary, 1.0,
				                        EllipticPolicy() );
			double const rd =
				boost::math::ellint_rd( 0.0, complementary, 1.0,
				                        EllipticPolicy() );
			a = kSquared*( rd/3.0 - rf/2.0 );
			b = kSquared*( rf - rd/3.0 )/complementary;
		}

		/// grad_bar of the filament kernel, per unit current and per unit mu0.
		///
		///     d_r G = ( 1/2 pi )[ ( ( a + r )/d ) A
		///                         + ( a( a^2 - r^2 + dz^2 )/d^3 ) B ]
		///     d_z G = ( 1/2 pi )[ ( dz/d ) A - ( 2 a r dz/d^3 ) B ]
		///
		/// See Coils.hpp for the derivation and for why the coefficient of B is
		/// written a( a^2 - r^2 + dz^2 )/d^3 rather than the a/d - 2ar( a+r )/d^3
		/// the chain rule hands you: they are identical, and the second is a
		/// difference of products that cancels to zero AT the loop, which is
		/// exactly where B is largest.
		void filamentGradKernel( double r, double z, double loopRadius,
		                         double loopHeight, double d, double kSquared,
		                         double complementary, double &dR, double &dZ )
		{
			double bracketA = 0.0;
			double bracketB = 0.0;
			kernelBrackets( kSquared, complementary, bracketA, bracketB );

			double const dz = z - loopHeight;
			double const sum = loopRadius + r;
			double const dCubed = d*d*d;

			// a^2 - r^2 + dz^2, one subtraction rather than a difference of
			// products, and exactly zero at the loop.
			double const radial =
				loopRadius*( loopRadius*loopRadius - r*r + dz*dz )/dCubed;

			dR = ( ( sum/d )*bracketA + radial*bracketB )/( 2.0*pi );
			dZ = ( ( dz/d )*bracketA
			       - ( 2.0*loopRadius*r*dz/dCubed )*bracketB )/( 2.0*pi );
		}

		void requireFinite( double value, char const *what, char const *where )
		{
			if ( !std::isfinite( value ) )
			{
				std::ostringstream message;
				message << where << ": " << what << " must be finite, but is "
				        << value;
				throw std::invalid_argument( message.str() );
			}
		}

		/// The refusal for a field point off the half-plane, WITH THE POINT IN
		/// IT. A caller who reaches this has extrapolated an evaluation past
		/// the axis, and the one thing they need to know is by how much: a
		/// radius at round-off is a different defect from one at element
		/// scale, and a message naming neither sends the reader to a debugger
		/// to recover a number the throw site already had.
		[[noreturn]] void refuseNegativeRadius( char const *where, double r,
		                                        double z )
		{
			std::ostringstream message;
			message << where << ": the field point radius must not be "
			           "negative; r = 0 is allowed and gives exactly zero. "
			           "The point is ( r, z ) = ( " << r << ", " << z << " )";
			throw std::invalid_argument( message.str() );
		}

		void requireOrder( int order, char const *where )
		{
			if ( order < 2 || order > maximumCoilQuadratureOrder )
			{
				std::ostringstream message;
				message << where << ": the quadrature order must be between 2 "
				           "and " << maximumCoilQuadratureOrder << ", but is "
				        << order
				        << ". The upper bound is a floating-point limit and not "
				           "a taste: the cubically graded nodes cluster as "
				           "n^-6 on the field point, and past it the innermost "
				           "one lands exactly on it. The rule reaches round-off "
				           "at about 48, so nothing useful is refused";
				throw std::invalid_argument( message.str() );
			}
		}

		/// The cross-section integral of the unit kernel over one ellipse for
		/// a field point INSIDE it, swept as chords along rays from the field
		/// point itself.
		///
		/// The area element rho d rho d theta carries the factor rho that
		/// meets the kernel's log( 1/rho ); cubic grading from rho = 0 then
		/// leaves t^5 log t, which Gauss takes to round-off. No node can land
		/// on the field point, every one being at rho > 0.
		double ellipseChordSweep( double r, double z, double centreR,
		                          double centreZ, double semiR, double semiZ,
		                          int order )
		{
			GaussRule const &rule = gaussRule( order );

			// Four angles per chord node. The theta rule is the accuracy
			// limit outside the ellipse and costs nothing on a ray that
			// misses, so it is the cheaper of the two axes to spend on.
			int const angles = 4*order;
			double const step = 2.0*pi/angles;

			// The field point in the ellipse's own coordinates. cSquared is
			// negative inside, zero on the boundary and positive outside, and
			// it is the constant term of the chord quadratic.
			double const u = ( r - centreR )/semiR;
			double const v = ( z - centreZ )/semiZ;
			double const outside = u*u + v*v - 1.0;

			double total = 0.0;
			for ( int m = 0; m < angles; ++m )
			{
				// The MIDPOINT rule, which is the trapezoid for a periodic
				// integrand and puts no node on theta = 0, where a field
				// point level with the centre would have its chord endpoints.
				double const theta = ( m + 0.5 )*step;
				double const cosine = std::cos( theta );
				double const sine = std::sin( theta );

				// Where the ray meets the ellipse: substituting
				// ( r + rho cos, z + rho sin ) into ( . /semiR )^2 +
				// ( . /semiZ )^2 = 1 gives a quadratic in rho.
				double const quadratic = cosine*cosine/( semiR*semiR )
				                         + sine*sine/( semiZ*semiZ );
				double const linear = 2.0*( u*cosine/semiR + v*sine/semiZ );
				double const discriminant =
					linear*linear - 4.0*quadratic*outside;

				// The ray misses the ellipse entirely, which is most of them
				// for a distant field point and is why one is cheap.
				if ( !( discriminant > 0.0 ) )
					continue;

				double const root = std::sqrt( discriminant );
				double const far = ( -linear + root )/( 2.0*quadratic );
				if ( !( far > 0.0 ) )
					continue;

				double const nearRoot = ( -linear - root )/( 2.0*quadratic );
				// Clipped at the field point itself: inside the ellipse the
				// near root is behind the ray and the chord starts at rho = 0,
				// which is the singular case the grading below is for.
				double const near = std::max( 0.0, nearRoot );
				double const length = far - near;
				if ( !( length > 0.0 ) )
					continue;

				double chord = 0.0;
				for ( int i = 0; i < order; ++i )
				{
					std::size_t const iu = static_cast<std::size_t>( i );
					double const t = 0.5*( rule.abscissa[ iu ] + 1.0 );
					double const tCubed = t*t*t;
					double const jacobian = 3.0*t*t;

					// Graded toward the near end, which is the singularity
					// when the point is inside and the nearest approach when
					// it is not.
					double const rho = near + length*tCubed;
					double const weight =
						0.5*rule.weight[ iu ]*length*jacobian;

					double const source = r + rho*cosine;
					double const height = z + rho*sine;

					double d = 0.0;
					double kSquared = 0.0;
					double complementary = 0.0;
					loopGeometry( r, z, source, height, d, kSquared,
					              complementary );

					// rho > 0 at every node, so this cannot fire on the field
					// point. It guards the node rounding onto the SOURCE ring
					// through the axis, which the refusal in ellipsePsi()
					// already makes unreachable.
					if ( !( complementary > 0.0 ) )
						continue;

					chord += weight*rho
					         *geometricKernel( d, kSquared, complementary );
				}
				total += step*chord;
			}

			return total;
		}

		/// The same integral for a field point OUTSIDE the ellipse, swept in
		/// the ELLIPSE'S OWN polar coordinates: every node lands in the
		/// ellipse and none of them is near the field point.
		///
		/// **THE CHORD SWEEP CANNOT DO THIS AND THE FAILURE IS SEVERE RATHER
		/// THAN GRADUAL.** Seen from outside, the ellipse subtends a CONE, and
		/// a uniform rule in theta over the whole circle spends its nodes
		/// mostly on rays that miss. A distant or small ellipse subtends a few
		/// degrees, so a handful of nodes -- sometimes none -- carry the entire
		/// integral. Measured before this branch existed: the far field of a
		/// shrinking column stopped converging to its filament at 1e-04 and
		/// then went BACKWARDS, and Delta* outside the column read order the
		/// interior value with an erratic sign.
		///
		/// Here the integrand is analytic over the whole disc -- the field
		/// point is off it -- so Gauss in s and the midpoint rule in phi are
		/// both spectral, and the Jacobian s makes the polar origin ordinary.
		/// What is left is a field point CLOSE to the boundary, where the
		/// kernel is nearly singular just outside the domain of integration;
		/// that is the one place this is three figures rather than machine
		/// precision, and it is a guess's worth. See Coils.hpp.
		double ellipseDiscSweep( double r, double z, double centreR,
		                         double centreZ, double semiR, double semiZ,
		                         int order )
		{
			GaussRule const &rule = gaussRule( order );

			// FOUR PER CHORD NODE RATHER THAN TWO, AND IT IS THE BOUNDARY
			// BAND THAT PAYS FOR IT. Halving this halves the cost -- 66 us a
			// point to 34 -- and leaves the far field and Delta* where they
			// were, because both are spectral out there. What it costs is the
			// one place this rule is weakest: the step across the boundary
			// between the two sweeps goes from 9.2e-04 to 3.5e-03, which
			// the_two_elliptical_sweeps_agree_across_their_own_seam refuses.
			int const angles = 4*order;
			double const step = 2.0*pi/angles;

			double total = 0.0;
			for ( int m = 0; m < angles; ++m )
			{
				double const phi = ( m + 0.5 )*step;
				double const cosine = std::cos( phi );
				double const sine = std::sin( phi );

				double ray = 0.0;
				for ( int i = 0; i < order; ++i )
				{
					std::size_t const iu = static_cast<std::size_t>( i );
					double const s = 0.5*( rule.abscissa[ iu ] + 1.0 );
					double const weight = 0.5*rule.weight[ iu ];

					double const source = centreR + semiR*s*cosine;
					double const height = centreZ + semiZ*s*sine;

					double d = 0.0;
					double kSquared = 0.0;
					double complementary = 0.0;
					loopGeometry( r, z, source, height, d, kSquared,
					              complementary );

					if ( !( complementary > 0.0 ) )
						continue;

					ray += weight*s
					       *geometricKernel( d, kSquared, complementary );
				}
				total += step*semiR*semiZ*ray;
			}

			return total;
		}

		/// The cross-section integral of the unit kernel over one ellipse.
		/// Multiply by mu0 and by the current density to get psi.
		///
		/// TWO RULES, SPLIT ON WHERE THE FIELD POINT IS, because the
		/// singularity is in the domain for one of them and outside it for the
		/// other and no single rule is good at both. See each sweep for what
		/// it is for and what it is bad at.
		double ellipseIntegral( double r, double z, double centreR,
		                        double centreZ, double semiR, double semiZ,
		                        int order )
		{
			double const u = ( r - centreR )/semiR;
			double const v = ( z - centreZ )/semiZ;

			// The boundary itself goes to the chord sweep, which grades onto
			// it; the disc sweep would have the singularity ON its own edge.
			if ( u*u + v*v <= 1.0 )
				return ellipseChordSweep( r, z, centreR, centreZ, semiR,
				                          semiZ, order );

			return ellipseDiscSweep( r, z, centreR, centreZ, semiR, semiZ,
			                         order );
		}

		/// The cross-section integral of the unit kernel over one coil, panelled
		/// and graded as the file comment describes. Multiply by mu0 and by the
		/// current density to get psi.
		double crossSectionIntegral( Coil const &coil, double r, double z,
		                             int order )
		{
			GaussRule const &rule = gaussRule( order );

			double const rLow = coil.rMin();
			double const rHigh = coil.rMax();
			double const zLow = coil.zMin();
			double const zHigh = coil.zMax();

			// The field point clamped into the rectangle. This one line is what
			// makes the interior, edge-on and exterior cases one piece of code:
			// it is the singularity when the point is inside and the nearest
			// point of the coil when it is not.
			double const rSplit = std::min( std::max( r, rLow ), rHigh );
			double const zSplit = std::min( std::max( z, zLow ), zHigh );

			double rEdge[ 3 ] = { rLow, rHigh, rHigh };
			int rPanels = 1;
			if ( rSplit > rLow && rSplit < rHigh )
			{
				rEdge[ 1 ] = rSplit;
				rEdge[ 2 ] = rHigh;
				rPanels = 2;
			}

			double zEdge[ 3 ] = { zLow, zHigh, zHigh };
			int zPanels = 1;
			if ( zSplit > zLow && zSplit < zHigh )
			{
				zEdge[ 1 ] = zSplit;
				zEdge[ 2 ] = zHigh;
				zPanels = 2;
			}

			double total = 0.0;
			for ( int pr = 0; pr < rPanels; ++pr )
			{
				double const rA = rEdge[ pr ];
				double const rB = rEdge[ pr + 1 ];
				double const rLength = rB - rA;
				// Grade toward whichever end of this panel is nearer the field
				// point. For a split panel that is the shared corner; for an
				// unsplit one it is the nearer edge of the coil.
				bool const rFromLow =
					( std::abs( rA - rSplit ) <= std::abs( rB - rSplit ) );

				for ( int pz = 0; pz < zPanels; ++pz )
				{
					double const zA = zEdge[ pz ];
					double const zB = zEdge[ pz + 1 ];
					double const zLength = zB - zA;
					bool const zFromLow =
						( std::abs( zA - zSplit ) <= std::abs( zB - zSplit ) );

					double panel = 0.0;
					for ( int i = 0; i < order; ++i )
					{
						std::size_t const iu = static_cast<std::size_t>( i );
						double const u = 0.5*( rule.abscissa[ iu ] + 1.0 );
						double const uCubed = u*u*u;
						double const uJacobian = 3.0*u*u;

						double const source = rFromLow
							? rA + rLength*uCubed
							: rB - rLength*uCubed;
						double const rWeight =
							0.5*rule.weight[ iu ]*rLength*uJacobian;

						for ( int j = 0; j < order; ++j )
						{
							std::size_t const jv =
								static_cast<std::size_t>( j );
							double const v =
								0.5*( rule.abscissa[ jv ] + 1.0 );
							double const vCubed = v*v*v;
							double const vJacobian = 3.0*v*v;

							double const height = zFromLow
								? zA + zLength*vCubed
								: zB - zLength*vCubed;
							double const zWeight =
								0.5*rule.weight[ jv ]*zLength*vJacobian;

							double d = 0.0;
							double kSquared = 0.0;
							double complementary = 0.0;
							loopGeometry( r, z, source, height, d, kSquared,
							              complementary );

							// The node has rounded onto the field point. See
							// the file comment: dropped, not evaluated, and
							// unreachable at the shipped grading exponent.
							if ( !( complementary > 0.0 ) )
								continue;

							panel += rWeight*zWeight
							         *geometricKernel( d, kSquared,
							                           complementary );
						}
					}
					total += panel;
				}
			}

			return total;
		}

		/// The cross-section integral of the GRADIENT kernel, panelled and
		/// graded exactly as crossSectionIntegral() is. Differentiation under
		/// the integral sign: the coil is fixed and the field point is the
		/// variable, so d/dr passes through.
		///
		/// The integrand is 1/distance where psi's is a logarithm, so inside
		/// the coil this is harder for the graded rule than psi is. Outside it
		/// -- FB-7's case, the conductor being beyond Gamma -- both are
		/// analytic and both are spectral.
		void crossSectionGradient( Coil const &coil, double r, double z,
		                           int order, double &outR, double &outZ )
		{
			GaussRule const &rule = gaussRule( order );

			double const rLow = coil.rMin();
			double const rHigh = coil.rMax();
			double const zLow = coil.zMin();
			double const zHigh = coil.zMax();

			double const rSplit = std::min( std::max( r, rLow ), rHigh );
			double const zSplit = std::min( std::max( z, zLow ), zHigh );

			double rEdge[ 3 ] = { rLow, rHigh, rHigh };
			int rPanels = 1;
			if ( rSplit > rLow && rSplit < rHigh )
			{
				rEdge[ 1 ] = rSplit;
				rEdge[ 2 ] = rHigh;
				rPanels = 2;
			}

			double zEdge[ 3 ] = { zLow, zHigh, zHigh };
			int zPanels = 1;
			if ( zSplit > zLow && zSplit < zHigh )
			{
				zEdge[ 1 ] = zSplit;
				zEdge[ 2 ] = zHigh;
				zPanels = 2;
			}

			double totalR = 0.0;
			double totalZ = 0.0;
			for ( int pr = 0; pr < rPanels; ++pr )
			{
				double const rA = rEdge[ pr ];
				double const rB = rEdge[ pr + 1 ];
				double const rLength = rB - rA;
				bool const rFromLow =
					( std::abs( rA - rSplit ) <= std::abs( rB - rSplit ) );

				for ( int pz = 0; pz < zPanels; ++pz )
				{
					double const zA = zEdge[ pz ];
					double const zB = zEdge[ pz + 1 ];
					double const zLength = zB - zA;
					bool const zFromLow =
						( std::abs( zA - zSplit ) <= std::abs( zB - zSplit ) );

					for ( int i = 0; i < order; ++i )
					{
						std::size_t const iu = static_cast<std::size_t>( i );
						double const u = 0.5*( rule.abscissa[ iu ] + 1.0 );
						double const uCubed = u*u*u;
						double const uJacobian = 3.0*u*u;

						double const source = rFromLow
							? rA + rLength*uCubed
							: rB - rLength*uCubed;
						double const rWeight =
							0.5*rule.weight[ iu ]*rLength*uJacobian;

						for ( int j = 0; j < order; ++j )
						{
							std::size_t const jv =
								static_cast<std::size_t>( j );
							double const v =
								0.5*( rule.abscissa[ jv ] + 1.0 );
							double const vCubed = v*v*v;
							double const vJacobian = 3.0*v*v;

							double const height = zFromLow
								? zA + zLength*vCubed
								: zB - zLength*vCubed;
							double const zWeight =
								0.5*rule.weight[ jv ]*zLength*vJacobian;

							double d = 0.0;
							double kSquared = 0.0;
							double complementary = 0.0;
							loopGeometry( r, z, source, height, d, kSquared,
							              complementary );

							// The same underflow guard crossSectionIntegral
							// carries, and for the same reason.
							if ( !( complementary > 0.0 ) )
								continue;

							double gradR = 0.0;
							double gradZ = 0.0;
							filamentGradKernel( r, z, source, height, d,
							                    kSquared, complementary,
							                    gradR, gradZ );

							double const weight = rWeight*zWeight;
							totalR += weight*gradR;
							totalZ += weight*gradZ;
						}
					}
				}
			}

			outR = totalR;
			outZ = totalZ;
		}
	}

	Coil::Coil( double centreRIn, double centreZIn, double halfWidthIn,
	            double halfHeightIn, double currentIn )
		: centreRValue( centreRIn ),
		  centreZValue( centreZIn ),
		  halfWidthValue( halfWidthIn ),
		  halfHeightValue( halfHeightIn ),
		  currentValue( currentIn )
	{
		requireFinite( centreRIn, "the centre radius", "meq::Coil" );
		requireFinite( centreZIn, "the centre height", "meq::Coil" );
		requireFinite( halfWidthIn, "the half-width", "meq::Coil" );
		requireFinite( halfHeightIn, "the half-height", "meq::Coil" );
		requireFinite( currentIn, "the current", "meq::Coil" );

		if ( !( halfWidthIn > 0.0 ) )
			throw std::invalid_argument(
				"meq::Coil: the half-width must be positive. A coil of zero "
				"radial extent is a filament, which has an infinite current "
				"density and no cross-section integral; filamentPsi() is the "
				"function for that case" );
		if ( !( halfHeightIn > 0.0 ) )
			throw std::invalid_argument(
				"meq::Coil: the half-height must be positive, for the reason "
				"the half-width must be" );

		// The same refusal meq::BoundaryShape makes, and for the same reason:
		// the Grad-Shafranov operator carries a 1/r which is not integrable
		// through the axis, so a coil reaching it is not merely unusual, it is
		// unsolvable. The filament kernel would give out there too -- psi is
		// exactly zero on the axis, so a coil straddling it would be
		// integrating through its own zero.
		if ( !( centreRIn - halfWidthIn > 0.0 ) )
		{
			std::ostringstream message;
			message << "meq::Coil: the coil reaches r = "
			        << centreRIn - halfWidthIn
			        << ", which is on or beyond the axis; the Grad-Shafranov "
			           "operator's 1/r is not integrable there";
			throw std::invalid_argument( message.str() );
		}
	}

	double Coil::centreR() const
	{
		return centreRValue;
	}

	double Coil::centreZ() const
	{
		return centreZValue;
	}

	double Coil::halfWidth() const
	{
		return halfWidthValue;
	}

	double Coil::halfHeight() const
	{
		return halfHeightValue;
	}

	double Coil::current() const
	{
		return currentValue;
	}

	double Coil::rMin() const
	{
		return centreRValue - halfWidthValue;
	}

	double Coil::rMax() const
	{
		return centreRValue + halfWidthValue;
	}

	double Coil::zMin() const
	{
		return centreZValue - halfHeightValue;
	}

	double Coil::zMax() const
	{
		return centreZValue + halfHeightValue;
	}

	double Coil::area() const
	{
		return 4.0*halfWidthValue*halfHeightValue;
	}

	double Coil::currentDensity() const
	{
		return currentValue/area();
	}

	bool Coil::contains( double r, double z ) const
	{
		// Closed, edges included; see the header for why that is stated rather
		// than merely chosen.
		return r >= rMin() && r <= rMax() && z >= zMin() && z <= zMax();
	}

	double filamentPsi( double r, double z, double loopRadius,
	                    double loopHeight, double current, double mu0 )
	{
		requireFinite( r, "the field point radius", "meq::filamentPsi" );
		requireFinite( z, "the field point height", "meq::filamentPsi" );
		requireFinite( loopRadius, "the loop radius", "meq::filamentPsi" );
		requireFinite( loopHeight, "the loop height", "meq::filamentPsi" );
		requireFinite( current, "the current", "meq::filamentPsi" );
		requireFinite( mu0, "mu0", "meq::filamentPsi" );

		if ( !( loopRadius > 0.0 ) )
			throw std::invalid_argument(
				"meq::filamentPsi: the loop radius must be positive" );
		if ( r < 0.0 )
			refuseNegativeRadius( "meq::filamentPsi", r, z );

		double d = 0.0;
		double kSquared = 0.0;
		double complementary = 0.0;
		loopGeometry( r, z, loopRadius, loopHeight, d, kSquared, complementary );

		if ( !( complementary > 0.0 ) )
			throw std::invalid_argument(
				"meq::filamentPsi: the field point is ON the loop, where psi "
				"is genuinely infinite -- this is a line current and the flux "
				"diverges logarithmically at it. A caller who has reached here "
				"has a geometry error rather than a large number" );

		// k^2 is an EXACT factor, so this is bit-exactly 0.0 on the axis, which
		// is the boundary condition the free-boundary problem imposes there.
		return mu0*current*geometricKernel( d, kSquared, complementary );
	}

	void filamentGradPsi( double r, double z, double loopRadius,
	                      double loopHeight, double current,
	                      double &dPsiDr, double &dPsiDz, double mu0 )
	{
		requireFinite( r, "the field point radius", "meq::filamentGradPsi" );
		requireFinite( z, "the field point height", "meq::filamentGradPsi" );
		requireFinite( loopRadius, "the loop radius", "meq::filamentGradPsi" );
		requireFinite( loopHeight, "the loop height", "meq::filamentGradPsi" );
		requireFinite( current, "the current", "meq::filamentGradPsi" );
		requireFinite( mu0, "mu0", "meq::filamentGradPsi" );

		if ( !( loopRadius > 0.0 ) )
			throw std::invalid_argument(
				"meq::filamentGradPsi: the loop radius must be positive" );
		if ( r < 0.0 )
			throw std::invalid_argument(
				"meq::filamentGradPsi: the field point radius must not be "
				"negative; r = 0 is allowed and gives exactly zero" );

		double d = 0.0;
		double kSquared = 0.0;
		double complementary = 0.0;
		loopGeometry( r, z, loopRadius, loopHeight, d, kSquared, complementary );

		if ( !( complementary > 0.0 ) )
			throw std::invalid_argument(
				"meq::filamentGradPsi: the field point is ON the loop, where "
				"the gradient is genuinely infinite -- this is a line current "
				"and grad psi diverges like 1/distance at it. There is no "
				"self-field and no self-force for a filament; meq::Coil, which "
				"has a finite cross-section, is the class for that" );

		double gradR = 0.0;
		double gradZ = 0.0;
		filamentGradKernel( r, z, loopRadius, loopHeight, d, kSquared,
		                    complementary, gradR, gradZ );

		// Both brackets carry k^2 as an explicit factor and k^2 is exactly zero
		// on the axis, so this is bit-exactly ( 0, 0 ) there -- the limit the
		// literal chain-rule form cannot reach.
		dPsiDr = mu0*current*gradR;
		dPsiDz = mu0*current*gradZ;
	}

	void filamentFlux( double r, double z, double loopRadius,
	                   double loopHeight, double current,
	                   double &qR, double &qZ, double mu0 )
	{
		filamentGradPsi( r, z, loopRadius, loopHeight, current, qR, qZ, mu0 );
		// 0/0 on the axis, deliberately not special-cased: the limit is finite
		// but this expression does not reach it, and a caller reading NaN has
		// learned something a quietly substituted limit would have hidden.
		qR /= r;
		qZ /= r;
	}

	double coilPsi( Coil const &coil, double r, double z, int order,
	                double mu0 )
	{
		requireFinite( r, "the field point radius", "meq::coilPsi" );
		requireFinite( z, "the field point height", "meq::coilPsi" );
		requireFinite( mu0, "mu0", "meq::coilPsi" );
		requireOrder( order, "meq::coilPsi" );

		if ( r < 0.0 )
			refuseNegativeRadius( "meq::coilPsi", r, z );

		return mu0*coil.currentDensity()
		       *crossSectionIntegral( coil, r, z, order );
	}

	void coilGradPsi( Coil const &coil, double r, double z,
	                  double &dPsiDr, double &dPsiDz, int order, double mu0 )
	{
		requireFinite( r, "the field point radius", "meq::coilGradPsi" );
		requireFinite( z, "the field point height", "meq::coilGradPsi" );
		requireFinite( mu0, "mu0", "meq::coilGradPsi" );
		requireOrder( order, "meq::coilGradPsi" );

		if ( r < 0.0 )
			refuseNegativeRadius( "meq::coilGradPsi", r, z );

		double gradR = 0.0;
		double gradZ = 0.0;
		crossSectionGradient( coil, r, z, order, gradR, gradZ );

		double const scale = mu0*coil.currentDensity();
		dPsiDr = scale*gradR;
		dPsiDz = scale*gradZ;
	}

	void coilFlux( Coil const &coil, double r, double z,
	               double &qR, double &qZ, int order, double mu0 )
	{
		coilGradPsi( coil, r, z, qR, qZ, order, mu0 );
		qR /= r;
		qZ /= r;
	}

	double filamentAxisFlux( double z, double loopRadius, double loopHeight,
	                         double current, double mu0 )
	{
		requireFinite( z, "the field point height", "meq::filamentAxisFlux" );
		requireFinite( loopRadius, "the loop radius", "meq::filamentAxisFlux" );
		requireFinite( loopHeight, "the loop height", "meq::filamentAxisFlux" );
		requireFinite( current, "the current", "meq::filamentAxisFlux" );
		requireFinite( mu0, "mu0", "meq::filamentAxisFlux" );

		if ( !( loopRadius > 0.0 ) )
			throw std::invalid_argument(
				"meq::filamentAxisFlux: the loop radius must be positive" );

		// d is the distance from the field point to the ring, which on the axis
		// is the same for every point of it -- that is the whole reason this
		// case is closed form where the general one needs elliptic integrals.
		// It is strictly positive because the ring never meets the axis.
		double const d = std::hypot( loopRadius, z - loopHeight );
		return 0.5*mu0*current*loopRadius*loopRadius/( d*d*d );
	}

	double coilAxisFlux( Coil const &coil, double z, int order, double mu0 )
	{
		requireFinite( z, "the field point height", "meq::coilAxisFlux" );
		requireFinite( mu0, "mu0", "meq::coilAxisFlux" );
		requireOrder( order, "meq::coilAxisFlux" );

		// A PLAIN TENSOR RULE, and the header says why: the axis is outside
		// every meq::Coil by that class's own refusal, so nothing here is
		// singular and the panelling and cubic grading crossSectionIntegral()
		// carries would buy nothing.
		GaussRule const &rule = gaussRule( order );

		double const rLow = coil.rMin();
		double const rHigh = coil.rMax();
		double const zLow = coil.zMin();
		double const zHigh = coil.zMax();
		double const rLength = rHigh - rLow;
		double const zLength = zHigh - zLow;

		double total = 0.0;
		for ( int i = 0; i < order; ++i )
		{
			std::size_t const iu = static_cast<std::size_t>( i );
			double const source = rLow + rLength*0.5*( rule.abscissa[ iu ] + 1.0 );
			double const rWeight = 0.5*rule.weight[ iu ]*rLength;

			for ( int j = 0; j < order; ++j )
			{
				std::size_t const jv = static_cast<std::size_t>( j );
				double const height =
					zLow + zLength*0.5*( rule.abscissa[ jv ] + 1.0 );
				double const zWeight = 0.5*rule.weight[ jv ]*zLength;

				double const d = std::hypot( source, z - height );
				total += rWeight*zWeight
				         *0.5*source*source/( d*d*d );
			}
		}

		// The current DENSITY, as every other coil kernel in this file does it:
		// the cross-section integral above carries the geometry and the density
		// carries the amperes.
		return mu0*coil.currentDensity()*total;
	}

	double ellipsePsi( double r, double z, double centreR, double centreZ,
	                   double semiR, double semiZ, double current, int order,
	                   double mu0 )
	{
		requireFinite( r, "the field point radius", "meq::ellipsePsi" );
		requireFinite( z, "the field point height", "meq::ellipsePsi" );
		requireFinite( centreR, "the centre radius", "meq::ellipsePsi" );
		requireFinite( centreZ, "the centre height", "meq::ellipsePsi" );
		requireFinite( semiR, "the radial semi-axis", "meq::ellipsePsi" );
		requireFinite( semiZ, "the vertical semi-axis", "meq::ellipsePsi" );
		requireFinite( current, "the current", "meq::ellipsePsi" );
		requireFinite( mu0, "mu0", "meq::ellipsePsi" );
		requireOrder( order, "meq::ellipsePsi" );

		if ( r < 0.0 )
			throw std::invalid_argument(
				"meq::ellipsePsi: the field point radius must not be "
				"negative" );

		if ( !( semiR > 0.0 ) || !( semiZ > 0.0 ) )
			throw std::invalid_argument(
				"meq::ellipsePsi: both semi-axes must be positive" );

		// The same refusal meq::Coil and meq::CurrentFilament make, and for
		// the same reason: the operator's 1/r is not integrable through r = 0
		// and psi is identically zero there, so a conductor reaching the axis
		// would be sitting in its own zero.
		if ( !( centreR - semiR > 0.0 ) )
		{
			std::ostringstream message;
			message << "meq::ellipsePsi: the ellipse reaches or crosses the "
			           "axis: centreR - semiR = " << centreR - semiR
			        << " and must be positive";
			throw std::invalid_argument( message.str() );
		}

		// The current density of a uniform ellipse, I/( pi a b ).
		double const density = current/( pi*semiR*semiZ );

		return mu0*density
		       *ellipseIntegral( r, z, centreR, centreZ, semiR, semiZ, order );
	}

	CurrentFilament::CurrentFilament( double radiusIn, double heightIn,
	                                  double currentIn )
		: radiusValue( radiusIn ),
		  heightValue( heightIn ),
		  currentValue( currentIn )
	{
		requireFinite( radiusIn, "the radius", "meq::CurrentFilament" );
		requireFinite( heightIn, "the height", "meq::CurrentFilament" );
		requireFinite( currentIn, "the current", "meq::CurrentFilament" );

		// The same refusal meq::Coil makes for a coil reaching the axis, and
		// for the same reason: the operator's 1/r is not integrable through
		// r = 0, and psi is identically zero on the axis so a conductor there
		// would be sitting in its own zero.
		if ( !( radiusIn > 0.0 ) )
		{
			std::ostringstream message;
			message << "meq::CurrentFilament: the radius must be positive, but "
			           "is " << radiusIn
			        << ". A filament on or beyond the axis is the same refusal "
			           "meq::Coil makes for a coil reaching it: the "
			           "Grad-Shafranov operator's 1/r is not integrable there";
			throw std::invalid_argument( message.str() );
		}
	}

	double CurrentFilament::radius() const
	{
		return radiusValue;
	}

	double CurrentFilament::height() const
	{
		return heightValue;
	}

	double CurrentFilament::current() const
	{
		return currentValue;
	}

	double filamentPsi( CurrentFilament const &filament, double r, double z,
	                    double mu0 )
	{
		return filamentPsi( r, z, filament.radius(), filament.height(),
		                    filament.current(), mu0 );
	}

	void filamentGradPsi( CurrentFilament const &filament, double r, double z,
	                      double &dPsiDr, double &dPsiDz, double mu0 )
	{
		filamentGradPsi( r, z, filament.radius(), filament.height(),
		                 filament.current(), dPsiDr, dPsiDz, mu0 );
	}

	void filamentFlux( CurrentFilament const &filament, double r, double z,
	                   double &qR, double &qZ, double mu0 )
	{
		filamentFlux( r, z, filament.radius(), filament.height(),
		              filament.current(), qR, qZ, mu0 );
	}

	double filamentAxisFlux( CurrentFilament const &filament, double z,
	                         double mu0 )
	{
		return filamentAxisFlux( z, filament.radius(), filament.height(),
		                         filament.current(), mu0 );
	}

	CoilSet::CoilSet( double mu0In )
		: permeability( mu0In ),
		  quadratureOrderValue( defaultCoilQuadratureOrder )
	{
		requireFinite( mu0In, "mu0", "meq::CoilSet" );
		if ( !( mu0In > 0.0 ) )
			throw std::invalid_argument(
				"meq::CoilSet: mu0 must be positive. Zero would leave every "
				"coil silently inert -- f() would be identically zero and the "
				"solve would converge beautifully to a vacuum -- which is a "
				"worse outcome than an error. Normalised units want 1" );
	}

	void CoilSet::add( Coil const &coil )
	{
		coilList.push_back( coil );
	}

	std::size_t CoilSet::size() const
	{
		return coilList.size();
	}

	bool CoilSet::empty() const
	{
		return coilList.empty();
	}

	Coil const &CoilSet::coil( std::size_t index ) const
	{
		if ( index >= coilList.size() )
			throw std::out_of_range(
				"meq::CoilSet::coil: index " + std::to_string( index )
				+ " is outside a set of " + std::to_string( coilList.size() )
				+ " coils" );
		return coilList[ index ];
	}

	std::vector<Coil> const &CoilSet::coils() const
	{
		return coilList;
	}

	double CoilSet::f( double r, double z ) const
	{
		// F = mu0 r j_phi, derived in the header. Summed rather than
		// short-circuited on the first hit, because overlapping coils really do
		// add their current densities.
		double density = 0.0;
		for ( Coil const &c : coilList )
			if ( c.contains( r, z ) )
				density += c.currentDensity();

		return permeability*r*density;
	}

	double CoilSet::totalCurrent() const
	{
		double sum = 0.0;
		for ( Coil const &c : coilList )
			sum += c.current();
		return sum;
	}

	int CoilSet::indexContaining( double r, double z ) const
	{
		for ( std::size_t i = 0; i < coilList.size(); ++i )
			if ( coilList[ i ].contains( r, z ) )
				return static_cast<int>( i );
		return -1;
	}

	double CoilSet::psi( double r, double z ) const
	{
		double sum = 0.0;
		for ( Coil const &c : coilList )
			sum += coilPsi( c, r, z, quadratureOrderValue, permeability );
		return sum;
	}

	double CoilSet::psiOf( std::size_t index, double r, double z ) const
	{
		return coilPsi( coil( index ), r, z, quadratureOrderValue,
		                permeability );
	}

	void CoilSet::gradPsi( double r, double z,
	                       double &dPsiDr, double &dPsiDz ) const
	{
		// Delta* is linear, so the gradient of the sum is the sum of the
		// gradients, exactly as psi() sums the fields. An empty set gives
		// ( 0, 0 ) without evaluating anything.
		double totalR = 0.0;
		double totalZ = 0.0;
		for ( Coil const &one : coilList )
		{
			double gradR = 0.0;
			double gradZ = 0.0;
			coilGradPsi( one, r, z, gradR, gradZ, quadratureOrderValue,
			             permeability );
			totalR += gradR;
			totalZ += gradZ;
		}
		dPsiDr = totalR;
		dPsiDz = totalZ;
	}

	void CoilSet::gradPsiOf( std::size_t index, double r, double z,
	                         double &dPsiDr, double &dPsiDz ) const
	{
		coilGradPsi( coil( index ), r, z, dPsiDr, dPsiDz, quadratureOrderValue,
		             permeability );
	}

	void CoilSet::flux( double r, double z, double &qR, double &qZ ) const
	{
		gradPsi( r, z, qR, qZ );
		// NaN on the axis in both components. gradPsi() is exactly zero there
		// and this divides it by zero; the limit is finite and is not reached.
		// See the header -- a semicircle centred on the axis MEETS it.
		qR /= r;
		qZ /= r;
	}

	void CoilSet::setQuadratureOrder( int order )
	{
		requireOrder( order, "meq::CoilSet::setQuadratureOrder" );
		quadratureOrderValue = order;
	}

	int CoilSet::quadratureOrder() const
	{
		return quadratureOrderValue;
	}

	double CoilSet::mu0() const
	{
		return permeability;
	}

	ExteriorCoilSet::ExteriorCoilSet( double mu0In )
		: permeability( mu0In ),
		  quadratureOrderValue( defaultCoilQuadratureOrder )
	{
		requireFinite( mu0In, "mu0", "meq::ExteriorCoilSet" );
		if ( !( mu0In > 0.0 ) )
			throw std::invalid_argument(
				"meq::ExteriorCoilSet: mu0 must be positive. Zero would leave "
				"every conductor silently inert -- psi() would be identically "
				"zero on Gamma and the coupled solve would converge beautifully "
				"to a machine with no conductors in it -- which is a worse "
				"outcome than an error. Normalised units want 1" );
	}

	void ExteriorCoilSet::add( Coil const &coil )
	{
		coilList.push_back( coil );
	}

	void ExteriorCoilSet::add( CurrentFilament const &filament )
	{
		filamentList.push_back( filament );
	}

	std::size_t ExteriorCoilSet::size() const
	{
		return coilList.size() + filamentList.size();
	}

	bool ExteriorCoilSet::empty() const
	{
		return coilList.empty() && filamentList.empty();
	}

	std::size_t ExteriorCoilSet::coilCount() const
	{
		return coilList.size();
	}

	std::size_t ExteriorCoilSet::filamentCount() const
	{
		return filamentList.size();
	}

	Coil const &ExteriorCoilSet::coil( std::size_t index ) const
	{
		if ( index >= coilList.size() )
			throw std::out_of_range(
				"meq::ExteriorCoilSet::coil: index " + std::to_string( index )
				+ " is outside a set of " + std::to_string( coilList.size() )
				+ " rectangular coils. The filaments have an index space of "
				"their own -- see filament()" );
		return coilList[ index ];
	}

	CurrentFilament const &ExteriorCoilSet::filament( std::size_t index ) const
	{
		if ( index >= filamentList.size() )
			throw std::out_of_range(
				"meq::ExteriorCoilSet::filament: index "
				+ std::to_string( index ) + " is outside a set of "
				+ std::to_string( filamentList.size() )
				+ " filaments. The rectangular coils have an index space of "
				"their own -- see coil()" );
		return filamentList[ index ];
	}

	std::vector<Coil> const &ExteriorCoilSet::coils() const
	{
		return coilList;
	}

	std::vector<CurrentFilament> const &ExteriorCoilSet::filaments() const
	{
		return filamentList;
	}

	double ExteriorCoilSet::totalCurrent() const
	{
		double sum = 0.0;
		for ( Coil const &c : coilList )
			sum += c.current();
		for ( CurrentFilament const &one : filamentList )
			sum += one.current();
		return sum;
	}

	double ExteriorCoilSet::psi( double r, double z ) const
	{
		double sum = 0.0;
		for ( Coil const &c : coilList )
			sum += coilPsi( c, r, z, quadratureOrderValue, permeability );
		for ( CurrentFilament const &one : filamentList )
			sum += filamentPsi( one, r, z, permeability );
		return sum;
	}

	void ExteriorCoilSet::gradPsi( double r, double z,
	                               double &dPsiDr, double &dPsiDz ) const
	{
		// Delta* is linear, so the gradient of the sum is the sum of the
		// gradients, exactly as psi() sums the fields. An empty set gives
		// ( 0, 0 ) without evaluating anything.
		double totalR = 0.0;
		double totalZ = 0.0;

		for ( Coil const &one : coilList )
		{
			double gradR = 0.0;
			double gradZ = 0.0;
			coilGradPsi( one, r, z, gradR, gradZ, quadratureOrderValue,
			             permeability );
			totalR += gradR;
			totalZ += gradZ;
		}

		for ( CurrentFilament const &one : filamentList )
		{
			double gradR = 0.0;
			double gradZ = 0.0;
			filamentGradPsi( one, r, z, gradR, gradZ, permeability );
			totalR += gradR;
			totalZ += gradZ;
		}

		dPsiDr = totalR;
		dPsiDz = totalZ;
	}

	void ExteriorCoilSet::flux( double r, double z,
	                            double &qR, double &qZ ) const
	{
		gradPsi( r, z, qR, qZ );
		// NaN on the axis in both components. gradPsi() is exactly zero there
		// and this divides it by zero; the limit is finite and is not reached.
		// See the header -- a semicircle centred on the axis MEETS it, and
		// GradShafranovSolver::conductorNormalFlux() is where the rule lives.
		qR /= r;
		qZ /= r;
	}

	double ExteriorCoilSet::clearance( double centreZ, double rhoGamma ) const
	{
		requireFinite( centreZ, "the centre height",
		               "meq::ExteriorCoilSet::clearance" );
		requireFinite( rhoGamma, "the radius of Gamma",
		               "meq::ExteriorCoilSet::clearance" );
		if ( !( rhoGamma > 0.0 ) )
			throw std::invalid_argument(
				"meq::ExteriorCoilSet::clearance: the radius of Gamma must be "
				"positive; got " + std::to_string( rhoGamma ) );

		// An empty set clears everything, and says so with an infinity rather
		// than with a large number a caller might mistake for a measurement.
		double least = std::numeric_limits<double>::infinity();

		for ( Coil const &one : coilList )
		{
			// The nearest point of the CLOSED rectangle to ( 0, centreZ ). The
			// centre is on the axis and Coil refuses rMin <= 0, so the nearest
			// radius is always rMin; only the height needs clamping. Taking the
			// nearest point rather than the centre is what makes a coil
			// straddling Gamma report a negative clearance.
			double const nearZ = std::min( std::max( centreZ, one.zMin() ),
			                               one.zMax() );
			least = std::min( least,
			                  std::hypot( one.rMin(), nearZ - centreZ ) );
		}

		for ( CurrentFilament const &one : filamentList )
			least = std::min( least, std::hypot( one.radius(),
			                                     one.height() - centreZ ) );

		return least - rhoGamma;
	}

	void ExteriorCoilSet::setQuadratureOrder( int order )
	{
		requireOrder( order, "meq::ExteriorCoilSet::setQuadratureOrder" );
		quadratureOrderValue = order;
	}

	int ExteriorCoilSet::quadratureOrder() const
	{
		return quadratureOrderValue;
	}

	double ExteriorCoilSet::mu0() const
	{
		return permeability;
	}



	CoilAugmentedSource::CoilAugmentedSource( std::shared_ptr<Source const> plasmaIn,
	                                          std::shared_ptr<CoilSet const> coilsIn )
		: plasmaSource( std::move( plasmaIn ) ), coilSet( std::move( coilsIn ) )
	{
		if ( !plasmaSource )
			throw std::invalid_argument( "meq::CoilAugmentedSource: the plasma source is null" );
		if ( !coilSet )
			throw std::invalid_argument( "meq::CoilAugmentedSource: the coil set is null" );
	}

	double CoilAugmentedSource::f( double r, double z, double psi ) const
	{
		return plasmaSource->f( r, z, psi ) + coilSet->f( r, z );
	}

	double CoilAugmentedSource::dFdPsi( double r, double z, double psi ) const
	{
		// The coils contribute nothing: their current is data, not a function
		// of psi. Adding a zero here would be harmless and is left out so that
		// the asymmetry with f() is visible rather than buried.
		return plasmaSource->dFdPsi( r, z, psi );
	}

	Source const & CoilAugmentedSource::plasma() const { return *plasmaSource; }
	CoilSet const & CoilAugmentedSource::coils() const { return *coilSet; }

	CoilSet const * CoilAugmentedSource::conductors() const
	{
		return coilSet.get();
	}

	CoilAugmentedNormalisedSource::CoilAugmentedNormalisedSource(
		std::shared_ptr<NormalisedSource> plasmaIn,
		std::shared_ptr<CoilSet const> coilsIn )
		: plasmaSource( std::move( plasmaIn ) ), coilSet( std::move( coilsIn ) )
	{
		if ( !plasmaSource )
			throw std::invalid_argument( "meq::CoilAugmentedNormalisedSource: the plasma source is null" );
		if ( !coilSet )
			throw std::invalid_argument( "meq::CoilAugmentedNormalisedSource: the coil set is null" );
	}

	double CoilAugmentedNormalisedSource::f( double r, double z, double psi ) const
	{
		return plasmaSource->f( r, z, psi ) + coilSet->f( r, z );
	}

	double CoilAugmentedNormalisedSource::dFdPsi( double r, double z, double psi ) const
	{
		return plasmaSource->dFdPsi( r, z, psi );
	}

	void CoilAugmentedNormalisedSource::setNormalisation( double psiAxis, double psiBoundary )
	{
		plasmaSource->setNormalisation( psiAxis, psiBoundary );
	}

	double CoilAugmentedNormalisedSource::normalisation() const
	{
		return plasmaSource->normalisation();
	}

	double CoilAugmentedNormalisedSource::boundaryNormalisation() const
	{
		return plasmaSource->boundaryNormalisation();
	}

	void CoilAugmentedNormalisedSource::setCurrentScale( double scale )
	{
		plasmaSource->setCurrentScale( scale );
		NormalisedSource::setCurrentScale( scale );
	}

	double CoilAugmentedNormalisedSource::currentScale() const
	{
		return plasmaSource->currentScale();
	}

	double CoilAugmentedNormalisedSource::scaledF( double r, double z,
	                                               double psi ) const
	{
		return plasmaSource->scaledF( r, z, psi );
	}

	double CoilAugmentedNormalisedSource::scaledDFdPsi( double r, double z,
	                                                    double psi ) const
	{
		return plasmaSource->scaledDFdPsi( r, z, psi );
	}

	bool CoilAugmentedNormalisedSource::normalisationDerivatives(
		double r, double z, double psi, double &dFdAxis,
		double &dFdBoundary ) const
	{
		// The coil term is independent of both normalisations, so the wrapped
		// source's answer IS the sum's answer.
		return plasmaSource->normalisationDerivatives( r, z, psi, dFdAxis,
		                                               dFdBoundary );
	}

	void CoilAugmentedNormalisedSource::setPlasmaSupport( bool confined )
	{
		// ON THE WRAPPED SOURCE, not on this one. NormalisedSource::insidePlasma
		// reads whichever object's flag, and the plasma term is evaluated
		// through plasmaSource -- so setting it here would leave the source that
		// actually computes F unconfined, and the support would silently do
		// nothing.
		plasmaSource->setPlasmaSupport( confined );
		NormalisedSource::setPlasmaSupport( confined );
	}

	void CoilAugmentedNormalisedSource::freezePlasmaEdge( double axis,
	                                                     double boundary )
	{
		// ON THE WRAPPED SOURCE FIRST, for the reason setPlasmaSupport() is:
		// the plasma term is evaluated through plasmaSource, so it is that
		// object's insidePlasma() that decides where F is switched off. The
		// base call after it keeps this wrapper's own plasmaEdgeIsFrozen() and
		// supportAxis() honest, which is what the solver's connected-component
		// fill reads.
		plasmaSource->freezePlasmaEdge( axis, boundary );
		NormalisedSource::freezePlasmaEdge( axis, boundary );
	}

	void CoilAugmentedNormalisedSource::thawPlasmaEdge()
	{
		plasmaSource->thawPlasmaEdge();
		NormalisedSource::thawPlasmaEdge();
	}

	double CoilAugmentedNormalisedSource::fOutsidePlasma( double r, double z ) const
	{
		// The coil term and nothing else, and it is the SAME expression f() adds
		// -- so on an element the fill did not reach this class contributes
		// exactly what it would have contributed with the plasma term absent,
		// bit for bit rather than to round-off.
		return coilSet->f( r, z );
	}

	CoilSet const * CoilAugmentedNormalisedSource::conductors() const
	{
		return coilSet.get();
	}

	NormalisedSource & CoilAugmentedNormalisedSource::plasma() const { return *plasmaSource; }
	CoilSet const & CoilAugmentedNormalisedSource::coils() const { return *coilSet; }

}
