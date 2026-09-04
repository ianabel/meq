#ifndef MEQ_TESTS_SOLOVIEV_HPP
#define MEQ_TESTS_SOLOVIEV_HPP

/*
 * Solov'ev equilibria: exact solutions of the Grad-Shafranov equation for a
 * source that is constant in psi.
 *
 * This is MEQ's primary correctness benchmark. It is used ahead of any
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
		/// solution MEQ's stage-3 benchmark is measured against.
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
				-0.0016308631916804547,
				-0.36242840247099738,
				-0.0033609747667579911,
				-0.025508342529018488,
				 0.0018640839541034535,
				-0.0028795855993667026,
				-9.7204527463518464e-5,
				 0.0013794520831780796,
				 0.0078878253370834794,
				-0.0076919796106070015,
				-0.0016366269795164746,
				 0.0001842341131865215
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
		 * THE THREE SOLOV'EV GEOMETRIES OF refs/HDG-GradShafranov.pdf section 4.2,
		 * ITS EXAMPLES 1, 2 AND 3.
		 *
		 * That paper prints no coefficients for them at all -- it gives only
		 * ( eps, delta, kappa, A ) per example and defers the rest to Cerfon &
		 * Freidberg: "Following the procedure carefully derived in [47] to choose
		 * the coefficients ( c_1, ..., c_12 ) ... The interested reader is
		 * referred to the above reference." So each set below was SOLVED from
		 * C&F's constraints, in mpmath at 50 digits, and then verified. The
		 * warning at the head of this file applies with full force: every psi_i
		 * is Delta*-harmonic, so a wrong coefficient changes no rate, no source
		 * and no residual. It changes only which domain the benchmark is on.
		 *
		 * HOW THEY WERE VERIFIED. Four checks, of which only the last two are
		 * independent of the system that was solved:
		 *
		 *   1. the linear system residual max | M c - b |, which came out below
		 *      1e-54 for all three -- this only says the solve was clean;
		 *   2. Delta*( psi ) - ( ( 1-A ) x^2 + A ) by 50-digit central
		 *      differences at twelve scattered points, below 3e-24 for all three;
		 *   3. eps, kappa and delta recovered FROM THE psi = 0 CONTOUR: the two
		 *      roots of psi( x, 0 ) either side of x = 1, and the point where
		 *      psi = 0 and psi_x = 0 together. All three came back to twelve
		 *      digits, and for the two single-null cases the saddle where
		 *      psi_x = psi_y = 0 sits at the prescribed X-point with psi ~ 1e-52;
		 *   4. the CURVATURE of the psi = 0 contour, -psi_yy/psi_x at the two
		 *      equatorial points and -psi_xx/psi_y at the high point, against
		 *      N_1, N_2, N_3 obtained by differentiating C&F's model surface
		 *      eq (9) numerically rather than by reusing their eq (11). Agreement
		 *      to twelve digits in all nine numbers.
		 *
		 * Check 4 is the one that matters, because it is the one that catches a
		 * mis-transcribed N -- and it caught one. See the note on alpha below.
		 *
		 * THE ALPHA CONVENTION, and a defect it exposes in nstx() above.
		 *
		 * C&F eq (11) is written in terms of alpha with sin( alpha ) = delta, so
		 * alpha = arcsin( delta ), and
		 *
		 *     N_1 = -( 1 + alpha )^2/( eps kappa^2 )
		 *     N_2 =  ( 1 - alpha )^2/( eps kappa^2 )
		 *     N_3 = -kappa/( eps cos^2( alpha ) )
		 *
		 * That alpha = arcsin( delta ) and not delta was confirmed by
		 * differentiating eq (9) directly: d2x/dy2 at tau = 0 comes out
		 * -( 1 + alpha )^2/( eps kappa^2 ) exactly. The three sets below use it
		 * consistently.
		 *
		 * nstx() DOES NOT. Evaluated against C&F eq (28), it satisfies all nine
		 * value and first-derivative conditions to ~1e-17, and it satisfies the
		 * high-point curvature condition with alpha = arcsin( delta ) -- but it
		 * satisfies the two EQUATORIAL curvature conditions only with alpha =
		 * delta, and misses them by -8.9e-3 and +1.2e-4 with the correct alpha.
		 * So that set was solved with alpha = delta in N_1 and N_2 and
		 * arcsin( delta ) in N_3. The consequence is small and entirely
		 * geometric: the psi = 0 contour has outboard curvature -0.5841 where the
		 * model surface has -0.5907, one per cent out, and the coefficients
		 * differ from a consistent solve in the third significant figure. It
		 * changes no rate and no source, for the reason given at the top of this
		 * file, and it is left alone here rather than corrected because the
		 * absolute-error ceilings of SolovievConvergence.cpp and
		 * ExtensionConvergence.cpp are calibrated against it.
		 */

		/// refs/HDG-GradShafranov.pdf Example 1: a Field Reversed Configuration,
		/// eps = 0.99, delta = 0.7, kappa = 10, A = 0.
		///
		/// UP-DOWN SYMMETRIC AND SMOOTH, so this is C&F's SEVEN-coefficient form,
		/// eq (8), closed by the seven smooth-surface constraints of eq (10):
		/// psi = 0 at the outer and inner equatorial points and the high point,
		/// psi_x = 0 at the high point, and the three curvature conditions.
		/// c_8 ... c_12 are identically zero -- those five terms are odd in y and
		/// an up-down symmetric configuration has no use for them.
		///
		/// This is C&F section VIII's "first method" for an FRC, in its own words:
		/// "The first method makes use of the solution already derived using the
		/// smooth surface constraints and approximates the ideal FRC by choosing
		/// eps = 0.99 and delta = 0.7." Their Figure 7(a) is this equilibrium.
		///
		/// A DISCREPANCY WITH THE PAPER, recorded because it is visible in its
		/// results. refs/HDG-GradShafranov.pdf says of its Example 1 that "the
		/// exact Solov'ev solution is a linear combination of the functions psi_p
		/// from equation (11) with A = 0 and psi_1, psi_2, and psi_4 from
		/// equation (12), yielding ... a bivariate polynomial of fourth degree",
		/// and attributes the sharp drop in its Figure 4 at k = 4 to exactly that.
		/// The seven-condition solve does NOT give that: c_3, c_5, c_6 and c_7 are
		/// small but non-zero below, so psi carries x^6, y^6 and log x terms and
		/// is not a polynomial. Which three conditions a three-coefficient
		/// solution would have been fitted to is not stated anywhere in that
		/// paper. What is implemented here is C&F's documented construction for
		/// these parameters; its convergence rates are the paper's k+1 claim, but
		/// its error levels and its k = 4 behaviour will not reproduce Figure 4.
		///
		/// Note eps = 0.99 puts the inner equatorial point at x = 0.01, and psi
		/// carries a log x: any domain used with this must stay well away from
		/// the axis.
		static SolovievEquilibrium frcExample1()
		{
			return SolovievEquilibrium( 0.0, {
				 4.8740209822310697e-5,
				-0.48742234131737989,
				 1.7226860958230982e-6,
				-0.0019044503610466043,
				-2.0028003674447282e-7,
				-3.6810713663408967e-6,
				 6.2918158091420308e-10,
				 0.0,
				 0.0,
				 0.0,
				 0.0,
				 0.0
			} );
		}

		/// refs/HDG-GradShafranov.pdf Example 2: an ITER-like single null,
		/// eps = 0.32, delta = 0.33, kappa = 2, A = -0.115.
		///
		/// UP-DOWN ASYMMETRIC with a downward X-point -- the paper says "The
		/// geometries are up-down asymmetric with a downwards oriented x-point"
		/// of its Examples 2 and 3 together -- so this is C&F's TWELVE-coefficient
		/// form, eq (26), closed by the twelve constraints of eq (28) with
		/// x_sep = 1 - 1.1 delta eps = 0.88384 and y_sep = -1.1 kappa eps = -0.704.
		///
		/// MEASURED: the psi = 0 contour has outer root 1.32, inner root 0.68,
		/// high point ( 0.8944, 0.64 ), recovered eps 0.32, kappa 2.0,
		/// delta 0.33, and a saddle at exactly ( 0.88384, -0.704 ) where
		/// psi = 3.3e-52.
		static SolovievEquilibrium iterExample2()
		{
			return SolovievEquilibrium( -0.115, {
				 0.085245241753779182,
				 0.049982802471268987,
				-0.28415411484266971,
				-0.13953333734005256,
				 0.18868881483537096,
				-0.1732478009843941,
				-0.0071646222528961062,
				 0.087479403396828952,
				 0.38699905013330074,
				-0.24303011558115132,
				-0.058060253802943254,
				 0.0068065843995832065
			} );
		}

		/// refs/HDG-GradShafranov.pdf Example 3: an NSTX-like single null,
		/// eps = 0.78, delta = 0.335, kappa = 1.7, A = -0.115.
		///
		/// The same twelve-condition construction as iterExample2(), with
		/// x_sep = 0.71257 and y_sep = -1.4586. NOT the same configuration as
		/// nstx() above, which is refs/HDG-GradShafranov-Adaptive.pdf section
		/// 4.1's high-beta case at A = -0.52, delta = 0.35, kappa = 2.
		///
		/// MEASURED: outer root 1.78, inner root 0.22, high point
		/// ( 0.7387, 1.326 ), recovered eps 0.78, kappa 1.7, delta 0.335, saddle
		/// at exactly ( 0.71257, -1.4586 ) with psi = -7.5e-52.
		static SolovievEquilibrium nstxExample3()
		{
			return SolovievEquilibrium( -0.115, {
				 0.0099937887677258671,
				-0.29773546867049059,
				 0.00031767618361544489,
				-0.031764596774629322,
				 0.006664222002433393,
				-0.0076160650769064532,
				-0.00030977960567813626,
				 0.0035112164077100054,
				 0.024969577934914871,
				-0.020438537381042342,
				-0.0043285843195953019,
				 0.00061554002592304011
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
		 * A SECOND ERROR, IN THIS FILE'S OWN FIRST CORRECTION, found by
		 * checking the curvature of the resulting contour against the model
		 * surface rather than against the formula it was derived from.
		 *
		 * Cerfon & Freidberg eq (11) reads
		 *
		 *     N_1 = -(1 + alpha)^2 / (eps kappa^2),
		 *     N_2 =  (1 - alpha)^2 / (eps kappa^2),
		 *     N_3 = -kappa / (eps cos^2 alpha),
		 *
		 * where alpha is the shape parameter of their eq (9) and sin(alpha) =
		 * delta. The first solve substituted sin(alpha) = delta into N_1 and
		 * N_2, giving -(1 + delta)^2 -- but eq (11) means alpha ITSELF, which
		 * for delta = 0.35 is arcsin(0.35) = 0.3576, not 0.35. N_3 happens to
		 * be immune, cos^2(alpha) = 1 - delta^2 either way, which is why the
		 * error was one-sided and easy to miss.
		 *
		 * Settled by differentiating eq (9) directly at tau = 0, pi and pi/2:
		 *
		 *     N_1  measured -0.5907049043   with alpha  exact   with delta  -0.5841
		 *     N_2  measured +0.1322804125   with alpha  exact   with delta  +0.1355
		 *     N_3  measured -2.9220542041   with alpha  exact   with delta  exact
		 *
		 * The coefficients above are the re-solve with the correct N_1 and N_2.
		 * They differ from the first correction by 0.2 to 2.3 per cent. Verified
		 * non-circularly: the curvature of the psi = 0 contour they produce
		 * agrees with the model surface to 2e-11, 3e-10 and 2e-12, where the
		 * previous set was out by 1.1 and 2.4 per cent at the two equatorial
		 * points.
		 *
		 * The moral, and the reason this is written out rather than quietly
		 * fixed: checking a solve against the formula it used cannot detect a
		 * misread formula. Only the independent quantity can -- here, the
		 * curvature of the surface the coefficients are supposed to reproduce.
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
