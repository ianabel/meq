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
		/// The NSTX-like single-null case of HDG-GradShafranov-Adaptive.pdf
		/// section 4.1, at A = -0.52.
		///
		/// THESE ARE NOT THAT PAPER'S PRINTED COEFFICIENTS. They are solved
		/// from Cerfon & Freidberg's twelve constraints, because the printed
		/// set does not satisfy them -- see the note below, and
		/// nstxAsPublished() if the printed numbers are wanted.
		static SolovievEquilibrium nstx()
		{
			return SolovievEquilibrium( -0.52, {
				-0.0016207703709346553,
				-0.3627285631694282,
				-0.0032978198208932852,
				-0.025453326073564573,
				 0.0018935982298685441,
				-0.0028994822166496844,
				-9.9473997747192579e-5,
				 0.0013798134064081668,
				 0.0079226303825184995,
				-0.0077012853250729403,
				-0.0016424391448371881,
				 0.00018410580253187123
			} );
		}

		/// The coefficients exactly as printed in
		/// HDG-GradShafranov-Adaptive.pdf eq (22c).
		///
		/// Kept because it is still a perfectly good exact solution -- every
		/// psi_i is Delta*-harmonic, so ANY coefficients solve the equation --
		/// and because it is the only way to reproduce that paper's numbers.
		/// What it is NOT is the NSTX geometry it is described as: see below.
		static SolovievEquilibrium nstxAsPublished()
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
				-0.000044132956899,
				-0.001299619729855,
				 0.000072050578303
			} );
		}

		/*
		 * WHY THE PUBLISHED COEFFICIENTS ARE NOT USED.
		 *
		 * refs/CerfonFreidberg.pdf section IX gives the twelve constraints that
		 * determine c_1 ... c_12 for an up-down asymmetric single-null
		 * configuration -- its eq (28). For the NSTX parameters
		 * eps = 0.78, kappa = 2, delta = 0.35, with the X-point at
		 * x_sep = 1 - 1.1 delta eps = 0.6997, y_sep = -1.1 kappa eps = -1.716,
		 * they are twelve linear conditions on the c_i:
		 *
		 *     psi = 0             at the outer and inner equatorial points,
		 *                         the upper high point, and the X-point
		 *     psi_y = 0           at both equatorial points
		 *     psi_x = 0           at the high point
		 *     psi_x = psi_y = 0   at the X-point
		 *     three curvature conditions with C&F's N_1, N_2, N_3
		 *
		 * Solving that 12x12 system at A = -0.52 gives the set in nstx().
		 * MEASURED: it satisfies all twelve to ~1e-17, puts the X-point exactly
		 * at (0.699700, -1.716000) with psi = -2.6e-18, and reproduces
		 * Delta*(psi) = -F to 9.9e-15.
		 *
		 * The printed set satisfies none of them. Evaluated at the same four
		 * points it gives
		 *
		 *     outer equatorial   psi = -7.54e-3
		 *     inner equatorial   psi = +3.92e-4
		 *     upper high point   psi = +2.92e-2     <-- 11% of the axis flux
		 *     X-point            psi = -9.76e-3
		 *
		 * where all four should be zero. So the surface psi = 0 for the printed
		 * coefficients is not the NSTX boundary it is described as, and the
		 * true saddle sits at (0.6958, -1.8069) with psi = -8.7e-3 rather than
		 * at the prescribed X-point with psi = 0.
		 *
		 * A PREVIOUS ATTEMPT AT THIS, recorded because it was wrong in an
		 * instructive way. The paper prints c_10 identical to c_7, which looked
		 * like a lone typesetting duplicate, and c_10 was "corrected" on its
		 * own by imposing psi = 0 at the X-point -- one of the twelve
		 * conditions. That gave -1.28e-3 and did make the separatrix the zero
		 * level set. But with the other eleven coefficients also wrong it was
		 * fitting one condition of twelve, and the agreement it produced
		 * (kappa near 2 on interior surfaces) was coincidence. Fixing one
		 * coefficient cannot repair a set that is wrong throughout.
		 *
		 * WHAT THIS DOES NOT CHANGE. Every psi_i is Delta*-harmonic, so no
		 * choice of c_i changes F, changes Delta*(psi), or changes any
		 * convergence rate -- both sets are exact solutions and both give k+1.
		 * The coefficients only decide which domain the benchmark is posed on,
		 * and how large the absolute error is. Which is exactly why this went
		 * unnoticed: the only check that can see it is the absolute-error
		 * ceiling, and before this it had never been compared against anything
		 * but itself.
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
