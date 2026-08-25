#ifndef MEQ_TESTS_SOLOVIEV_HPP
#define MEQ_TESTS_SOLOVIEV_HPP

/*
 * Solov'ev equilibria: exact solutions of the Grad-Shafranov equation for a
 * source that is constant in psi.
 *
 * This is meq's primary correctness benchmark. It is used ahead of any
 * convergence study, because a wrong sign or weight convention converges at
 * the right rate to the wrong function -- an order-of-accuracy study cannot
 * detect that, and only a comparison against a closed form can.
 *
 * Sources:
 *   refs/HDG-GradShafranov.pdf           eqs (10)-(12)  -- the expansion
 *   refs/HDG-GradShafranov-Adaptive.pdf  eqs (21)-(22)  -- the NSTX coefficients
 *   refs/Refs.md                         doi 10.1063/1.3328818 (Cerfon & Freidberg)
 *                                        for the geometric parametrisation
 */

#include <array>
#include <cmath>

namespace meq
{
namespace analytic
{

/*
 * SIGN CONVENTION -- read this before using the class.
 *
 * The toroidal elliptic operator is
 *
 *     Delta*(psi) := r d_r( (1/r) d_r psi ) + d_zz psi,
 *
 * and the Grad-Shafranov equation, in the form both papers state it, is
 *
 *     -Delta*(psi) = F(r,z,psi),
 *
 * equivalently  -div_bar( (1/r) grad_bar psi ) = F/r.
 *
 * For the Solov'ev profiles, with the flux normalised so that A + C = 1,
 *
 *     F(r,z,psi) = -( (1 - A) r^2 + A ).                               (*)
 *
 * THE TWO PAPERS DISAGREE ON THIS SIGN, and (*) is the one that is right.
 *
 *   HDG-GradShafranov.pdf eq (10)          F = -( (1-A) r^2 + A )   <-- correct
 *   HDG-GradShafranov-Adaptive.pdf eq (21) F =  ( (1-A) r^2 + A )   <-- sign slip
 *
 * The second contradicts its own eq (1). Applying Delta* directly to the
 * particular solution both papers publish settles it:
 *
 *     Delta*( r^4/8 )              =  r^2
 *     Delta*( (A/2) r^2 ln r )     =  A
 *     Delta*( -A r^4/8 )           = -A r^2
 *     ------------------------------------------
 *     Delta*( psi_P )              =  (1 - A) r^2 + A
 *
 * so F = -Delta*(psi_P) = -( (1-A) r^2 + A ), which is eq (10). The homogeneous
 * terms psi_1..psi_12 contribute nothing, being Delta*-harmonic by construction.
 *
 * deltaStarFD() below recomputes this numerically, and the test suite asserts
 * it, so the claim is checked rather than merely asserted here.
 */

/// An exact Solov'ev equilibrium: a particular solution plus a linear
/// combination of the twelve Delta*-harmonic functions of
/// refs/HDG-GradShafranov.pdf eq (12).
///
/// All lengths are normalised to the major radius; r must be strictly
/// positive, since the expansion contains log(r).
class SolovievEquilibrium
{
	public:
		/// @param aIn  the papers' parameter A: fixes the ratio of plasma
		///             pressure to magnetic pressure. The flux is normalised
		///             so that A + C = 1.
		/// @param cIn  the twelve coefficients c_1 ... c_12, fixed by imposing
		///             geometric constraints (Cerfon & Freidberg's procedure).
		SolovievEquilibrium( double aIn, std::array<double, 12> const &cIn )
			: a( aIn ), c( cIn )
		{
		}

		/// The NSTX-like high-beta case of HDG-GradShafranov-Adaptive.pdf
		/// section 4.1, eq (22c). This is the configuration whose exact
		/// solution meq's stage-3 benchmark is measured against.
		///
		/// c_10 IS NOT THE PRINTED VALUE. See the note below the coefficient
		/// list for the derivation; the paper's eq (22c) prints c_10 identical
		/// to c_7, which is a typesetting duplicate.
		static SolovievEquilibrium nstx()
		{
			return SolovievEquilibrium( -0.52, {
				-0.00147796615575325,
				-0.366568333204813,
				 0.002409406149732,
				-0.023957517168316,
				 0.000692888519765,
				-0.001768712177298,
				-0.000044132956899,
				 0.000433522611526,
				 0.008286849573230,
				-0.0012801846791351433,   // c_10, corrected -- see below
				-0.001299619729855,
				 0.000072050578303
			} );
		}

		/*
		 * WHY c_10 IS NOT THE PUBLISHED NUMBER.
		 *
		 * refs/HDG-GradShafranov-Adaptive.pdf eq (22c) prints
		 *
		 *     c_7  = -0.000044132956899
		 *     c_10 = -0.000044132956899
		 *
		 * identically. That is a typesetting duplicate, and using it gives an
		 * equilibrium that is not normalised: a Solov'ev boundary is the level
		 * set psi = 0, and for a single-null configuration that surface passes
		 * through the X-point, so the X-point must have psi = 0. With the
		 * printed value it does not -- measured, the saddle sits at
		 * (0.6958, -1.8069) with psi = -8.7e-3, so the zero level set is not
		 * the separatrix at all.
		 *
		 * The value here is the root of psi( X( c_10 ) ) = 0, the X-point being
		 * continued along as c_10 varies. It is 29.0 times the printed
		 * magnitude, and it puts the X-point at (0.6795, -1.8443) with
		 * psi = -3e-17.
		 *
		 * CORROBORATION, which is the reason to trust this rather than merely
		 * to prefer it. The configuration is meant to be NSTX at high beta,
		 * elongation kappa = 2.0. Measuring flux surfaces just inside the
		 * separatrix gives kappa = 2.09, 2.06, 2.03 at psi = -0.02, -0.01,
		 * -0.005 -- converging on 2.0 from outside. The printed coefficient
		 * produces no closed configuration to measure, and a least-squares fit
		 * over Cerfon-Freidberg's conditions that was also tried gives
		 * c_10 = -2.87e-3, for which the X-point lands at psi = +1.2e-2 -- i.e.
		 * OUTSIDE the psi = 0 surface, which contradicts a single null -- and
		 * kappa = 1.75 with a degenerate eps = 0.98.
		 *
		 * WHAT THIS DOES NOT AFFECT. psi_1 ... psi_12 are Delta*-harmonic, so
		 * no choice of c_i changes F, changes Delta*(psi), or changes any
		 * convergence rate. That is exactly why the error was invisible: the
		 * only thing a wrong coefficient moves is the absolute error against
		 * the exact solution, and the ceilings in the convergence tests were
		 * calibrated to the wrong solution. They are recalibrated in the same
		 * commit as this change.
		 *
		 * WHAT IS STILL NOT VERIFIED. The other eleven coefficients. Nothing
		 * here constrains them, and the same argument says nothing in the suite
		 * can. Reproducing the full set needs Cerfon-Freidberg's twelve
		 * geometric conditions, and that paper is not in refs/ -- doi
		 * 10.1063/1.3328818 if it is ever wanted.
		 */

		/// The poloidal flux function.
		double psi( double r, double z ) const
		{
			double value = psiParticular( r, z );
			for ( int i = 0; i < 12; ++i )
			{
				value += c[ i ] * homogeneous( i, r, z );
			}
			return value;
		}

		/// grad_bar(psi) = ( d_r psi, d_z psi ). Not the HDG flux: that is
		/// this divided by r, see flux().
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			gradPsiParticular( r, z, dPsiDr, dPsiDz );
			for ( int i = 0; i < 12; ++i )
			{
				double gr, gz;
				gradHomogeneous( i, r, z, gr, gz );
				dPsiDr += c[ i ] * gr;
				dPsiDz += c[ i ] * gz;
			}
		}

		/// The HDG flux q = grad_bar(psi) / r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The Grad-Shafranov source, F = -( (1-A) r^2 + A ). Independent of
		/// psi, which is what makes this equilibrium linear.
		double f( double r, double /*z*/, double /*psiValue*/ ) const
		{
			return -( ( 1.0 - a ) * r * r + a );
		}

		/// dF/dpsi. Identically zero: the Solov'ev source does not depend on
		/// psi, so a Newton solve on this problem must converge in a single
		/// step. That is a useful property to assert.
		double dFdPsi( double, double, double ) const
		{
			return 0.0;
		}

		/// Delta*(psi), by central differences of psi(). Used to verify the
		/// sign convention documented above, and to catch a mistyped term in
		/// the expansion, without depending on the hand-derived gradients.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			// r d_r( (1/r) d_r psi ) as a second difference of the inner
			// quantity, plus d_zz psi.
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) ) / ( 2.0 * h ) / rr;
			};

			double dRInner = ( innerR( r + h ) - innerR( r - h ) ) / ( 2.0 * h );
			double dZZ = ( psi( r, z + h ) - 2.0 * psi( r, z ) + psi( r, z - h ) )
			             / ( h * h );

			return r * dRInner + dZZ;
		}

		/// The papers' parameter A.
		double getA() const
		{
			return a;
		}

	private:
		double a;
		std::array<double, 12> c;

		/// psi_P = r^4/8 + A( (1/2) r^2 ln r - r^4/8 ),
		/// refs/HDG-GradShafranov.pdf eq (11).
		double psiParticular( double r, double z ) const
		{
			(void)z;
			double r2 = r * r;
			double r4 = r2 * r2;
			return r4 / 8.0 + a * ( 0.5 * r2 * std::log( r ) - r4 / 8.0 );
		}

		void gradPsiParticular( double r, double z, double &gr, double &gz ) const
		{
			(void)z;
			double r2 = r * r;
			double r3 = r2 * r;
			double logR = std::log( r );
			gr = r3 / 2.0 + a * ( r * logR + r / 2.0 - r3 / 2.0 );
			gz = 0.0;
		}

		/// psi_1 ... psi_12 of refs/HDG-GradShafranov.pdf eq (12), indexed
		/// from zero. Each satisfies Delta*(psi_i) = 0.
		static double homogeneous( int i, double r, double z )
		{
			double r2 = r * r,  r4 = r2 * r2,  r6 = r4 * r2;
			double z2 = z * z,  z3 = z2 * z,   z4 = z2 * z2;
			double z5 = z4 * z, z6 = z4 * z2;
			double logR = std::log( r );

			switch ( i )
			{
				case 0:  return 1.0;
				case 1:  return r2;
				case 2:  return z2 - r2 * logR;
				case 3:  return r4 - 4.0 * r2 * z2;
				case 4:  return 2.0 * z4 - 9.0 * z2 * r2 + 3.0 * r4 * logR
					            - 12.0 * r2 * z2 * logR;
				case 5:  return r6 - 12.0 * r4 * z2 + 8.0 * r2 * z4;
				case 6:  return 8.0 * z6 - 140.0 * z4 * r2 + 75.0 * z2 * r4
					            - 15.0 * r6 * logR + 180.0 * r4 * z2 * logR
					            - 120.0 * r2 * z4 * logR;
				case 7:  return z;
				case 8:  return z * r2;
				case 9:  return z3 - 3.0 * z * r2 * logR;
				case 10: return 3.0 * z * r4 - 4.0 * z3 * r2;
				case 11: return 8.0 * z5 - 45.0 * z * r4 - 80.0 * z3 * r2 * logR
					            + 60.0 * z * r4 * logR;
				default: return 0.0;
			}
		}

		static void gradHomogeneous( int i, double r, double z,
		                             double &gr, double &gz )
		{
			double r2 = r * r,  r3 = r2 * r,  r4 = r2 * r2, r5 = r4 * r;
			double z2 = z * z,  z3 = z2 * z,  z4 = z2 * z2, z5 = z4 * z;
			double logR = std::log( r );

			switch ( i )
			{
				case 0:
					gr = 0.0;
					gz = 0.0;
					break;
				case 1:
					gr = 2.0 * r;
					gz = 0.0;
					break;
				case 2:
					gr = -2.0 * r * logR - r;
					gz = 2.0 * z;
					break;
				case 3:
					gr = 4.0 * r3 - 8.0 * r * z2;
					gz = -8.0 * r2 * z;
					break;
				case 4:
					gr = -18.0 * z2 * r + 12.0 * r3 * logR + 3.0 * r3
					     - 24.0 * r * z2 * logR - 12.0 * r * z2;
					gz = 8.0 * z3 - 18.0 * z * r2 - 24.0 * r2 * z * logR;
					break;
				case 5:
					gr = 6.0 * r5 - 48.0 * r3 * z2 + 16.0 * r * z4;
					gz = -24.0 * r4 * z + 32.0 * r2 * z3;
					break;
				case 6:
					gr = -280.0 * z4 * r + 300.0 * z2 * r3
					     - 90.0 * r5 * logR - 15.0 * r5
					     + 720.0 * r3 * z2 * logR + 180.0 * r3 * z2
					     - 240.0 * r * z4 * logR - 120.0 * r * z4;
					gz = 48.0 * z5 - 560.0 * z3 * r2 + 150.0 * z * r4
					     + 360.0 * r4 * z * logR - 480.0 * r2 * z3 * logR;
					break;
				case 7:
					gr = 0.0;
					gz = 1.0;
					break;
				case 8:
					gr = 2.0 * r * z;
					gz = r2;
					break;
				case 9:
					gr = -6.0 * z * r * logR - 3.0 * z * r;
					gz = 3.0 * z2 - 3.0 * r2 * logR;
					break;
				case 10:
					gr = 12.0 * z * r3 - 8.0 * z3 * r;
					gz = 3.0 * r4 - 12.0 * z2 * r2;
					break;
				case 11:
					gr = -180.0 * z * r3 - 160.0 * z3 * r * logR - 80.0 * z3 * r
					     + 240.0 * z * r3 * logR + 60.0 * z * r3;
					gz = 40.0 * z4 - 45.0 * r4 - 240.0 * z2 * r2 * logR
					     + 60.0 * r4 * logR;
					break;
				default:
					gr = 0.0;
					gz = 0.0;
					break;
			}
		}
};

}
}

#endif // MEQ_TESTS_SOLOVIEV_HPP
