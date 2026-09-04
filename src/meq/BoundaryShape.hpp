#ifndef MEQ_BOUNDARYSHAPE_HPP
#define MEQ_BOUNDARYSHAPE_HPP

/*
 * The plasma boundary as a parametrised closed curve, and as a level set.
 *
 * MEQ's fixed-boundary problem takes Gamma to be a known curve. On the fitted
 * path that curve is the mesh boundary; on the extension path of stage 5 it is
 * a level set, and everything downstream -- MarkLevelSetSubdomain, the SubMesh,
 * VertexConePath -- reaches the geometry only through a function that is
 * negative inside. So the whole job of this file is: parameters in, level set
 * out.
 *
 * NO MFEM HERE, deliberately, as in Profiles.hpp and Source.hpp. Plain doubles,
 * so the geometry is unit-testable without the finite element library and CI --
 * which cannot obtain MEQ's MFEM branch -- can test it. The mfem::Coefficient
 * adapter lives with the assembly that needs it.
 *
 * ONE PARAMETRISATION, NOT TWO. Miller and MXH are the same formula:
 *
 *   refs/MXH.pdf eqs (1)-(3), Arbon, Candy & Belli:
 *       R( theta ) = R0 + r cos( theta_R )
 *       Z( theta ) = Z0 + kappa r sin( theta )
 *       theta_R    = theta + c0 + sum_{n=1}^{N} [ c_n cos( n theta )
 *                                               + s_n sin( n theta ) ]
 *
 *   refs/Miller.pdf eq (34), Miller, Chu, Greene, Lin-Liu & Waltz:
 *       R = R0 + r cos[ theta + ( arcsin delta ) sin theta ]
 *       Z = kappa r sin theta
 *
 * and MXH eq (4) states the reduction outright: Turnbull-Miller is recovered by
 * keeping s_1 = arcsin( delta ) and s_2 = -zeta, zeta being the squareness. So
 * Miller is MXH truncated at N = 1, miller() is sugar over the general
 * constructor, and there is one evaluator to get wrong. Two implementations of
 * one curve is how Example 6's p-versus-F inconsistency happened.
 *
 * A TRANSCRIPTION NOTE. refs/HDG-GradShafranov.pdf Example 6 prints
 * arcsin( delta sin t ) where Miller and Cerfon & Freidberg both print
 * arcsin( delta ) sin t. The two agree at t = 0, +-pi/2 and pi and differ by at
 * most 6e-4 in r between; tests/analytic/MillerDShape.hpp implements Example 6's
 * form because it reproduces that paper's table, and this file implements the
 * other because two independent sources and the physical meaning of delta agree
 * on it. If a number here disagrees with MillerDShape by ~1e-4, that is why.
 *
 * PHYSICAL READINGS of the harmonics, from MXH section 2: s_1 is the
 * triangularity, -s_2 the squareness, c_0 the tilt, c_1 the ovality, and the
 * n-th harmonic induces an (n-2)-sided polygonal deformation.
 */

#include <stdexcept>
#include <string>
#include <vector>

namespace meq
{
	/// Thrown when a shape cannot be used: degenerate parameters, a curve that
	/// crosses the axis, or one that is not star shaped about its own centre.
	class ShapeError : public std::runtime_error
	{
		public:
			explicit ShapeError( std::string const &what )
				: std::runtime_error( what )
			{
			}
	};

	class BoundaryShape
	{
		public:
			/// The general MXH surface.
			///
			/// @param r0In        major radius of the centre, metres. > 0.
			/// @param z0In        height of the centre, metres.
			/// @param minorIn     r, the minor radius. > 0, and < r0In so that the
			///                    surface does not reach the axis, where the
			///                    operator's 1/r is not integrable.
			/// @param elongationIn kappa. > 0.
			/// @param cosIn       c_0, c_1, ... c_N. **Starts at c_0**, the tilt.
			///                    May be empty, which means no cosine harmonics.
			/// @param sinIn       s_1, s_2, ... s_N. **Starts at s_1**; there is no
			///                    s_0, because sin( 0 theta ) vanishes identically
			///                    and a slot for it would be a place to put a
			///                    number that silently does nothing.
			///
			/// @throws ShapeError if the parameters are degenerate, if the surface
			///         reaches r <= 0, or if it is not star shaped about
			///         ( r0In, z0In ) -- see levelSet() for why that last one is
			///         not a technicality.
			BoundaryShape( double r0In, double z0In, double minorIn,
			               double elongationIn,
			               std::vector<double> cosIn = {},
			               std::vector<double> sinIn = {} );

			/// The Miller D-shape, which is MXH at N <= 2.
			///
			/// @param deltaIn      triangularity, |delta| < 1. Enters as
			///                     s_1 = arcsin( delta ) -- NOT as delta itself.
			///                     That distinction cost this project a day once
			///                     already; see the Solov'ev coefficients in
			///                     CLAUDE.md.
			/// @param squarenessIn zeta, Turnbull-Miller's fourth parameter,
			///                     entering as s_2 = -zeta. Zero gives the
			///                     original three-parameter Miller shape.
			static BoundaryShape miller( double r0In, double z0In, double minorIn,
			                             double elongationIn, double deltaIn,
			                             double squarenessIn = 0.0 );

			/// The point on the curve at parameter @a theta, which is MXH's
			/// poloidal angle and not the polar angle about the centre.
			void point( double theta, double &r, double &z ) const;

			/// Negative inside Gamma, zero on it, positive outside.
			///
			/// This is the RADIAL GAP -- |( r, z ) - centre| minus the distance
			/// from the centre to the curve at the same polar angle -- and not a
			/// signed distance, which is larger. Nothing downstream needs a
			/// distance: MarkLevelSetSubdomain wants a sign at the vertices and
			/// VertexConePath wants a root along a ray, and both are unaffected.
			///
			/// The curve radius at a given polar angle is found by bisecting on
			/// theta, which is legitimate exactly when the polar angle is strictly
			/// increasing in theta -- that is, when the curve is star shaped about
			/// its centre. The constructor checks it, because a curve that is not
			/// would give bisection several roots to choose between and it would
			/// return one of them without complaint.
			double levelSet( double r, double z ) const;

			/// The polar angle about the centre of the curve point at @a theta,
			/// in [ 0, 2 pi ). Strictly increasing in theta for a star shaped
			/// curve, which is what the constructor verifies.
			double polarAngle( double theta ) const;

			double majorRadius() const { return r0; }
			double centreHeight() const { return z0; }
			double minorRadius() const { return minor; }
			double elongation() const { return elongationValue; }

			/// The bounding box, for checking that the background mesh contains
			/// the surface with room for the transfer paths to reach past it.
			void boundingBox( double &rMin, double &rMax,
			                  double &zMin, double &zMax ) const;

		private:
			/// theta_R, MXH eq (3).
			double shiftedAngle( double theta ) const;

			/// Bisect for the theta whose point has polar angle @a target.
			double parameterAtPolarAngle( double target ) const;

			/// Throws unless the polar angle increases strictly with theta.
			void requireStarShaped() const;

			double r0, z0, minor, elongationValue;
			std::vector<double> cosCoefficients;   ///< c_0 ... c_N
			std::vector<double> sinCoefficients;   ///< s_1 ... s_N
	};
}

#endif // MEQ_BOUNDARYSHAPE_HPP
