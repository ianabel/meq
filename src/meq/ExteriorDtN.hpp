#ifndef MEQ_EXTERIORDTN_HPP
#define MEQ_EXTERIORDTN_HPP

#include <cstddef>
#include <functional>
#include <vector>

/*
 * The exterior of a semicircle, as an operator on its own trace: stage FB-0 of
 * FREE-BOUNDARY-PLAN.md.
 *
 * WHAT PROBLEM THIS SOLVES. Free boundary makes the domain unbounded -- the
 * vacuum region runs from the plasma out to infinity, with psi -> 0 there --
 * and a finite element mesh cannot represent that. The standard remedy is an
 * artificial boundary Gamma with an exact operator on it that stands in for
 * everything outside: a Dirichlet-to-Neumann map, which given psi on Gamma
 * returns the normal derivative the exterior solution would have. Couple that
 * to the interior discretisation and the truncation costs nothing.
 *
 * FOR THE AXISYMMETRIC OPERATOR ON A SEMICIRCLE THAT MAP IS DIAGONAL, and this
 * class is that statement. No layer potentials, no singular quadrature, no
 * elliptic integrals, no O( N^2 ) kernel evaluations -- one number per mode.
 * The derivation is FREE-BOUNDARY-PLAN.md section 3 and is summarised below
 * because a reader of this header should not have to fetch it.
 *
 * THE SEPARATION. In spherical coordinates r = rho sin( theta ),
 * z = zCentre + rho cos( theta ), mu = cos( theta ), the Grad-Shafranov
 * operator is
 *
 *     Delta* = d_rhorho + ( 1/rho^2 )( d_thetatheta - cot( theta ) d_theta )
 *
 * -- the ( 2/rho ) d_rho of the Laplacian is cancelled exactly by the
 * -( 1/r ) d_r that distinguishes Delta* from it, which is the accident this
 * whole class rests on. Separating psi = rho^alpha f( mu ) gives
 *
 *     ( 1 - mu^2 ) f'' + alpha( alpha - 1 ) f = 0,
 *
 * the GEGENBAUER EQUATION OF ORDER -1/2. So alpha = n or alpha = 1 - n, and the
 * angular functions are
 *
 *     C_n( mu ) = ( P_{n-2}( mu ) - P_n( mu ) ) / ( 2n - 1 ),   n >= 2,
 *
 * with P the Legendre polynomials.
 *
 * THREE CONSEQUENCES, AND EACH REMOVES SOMETHING.
 *
 *   * C_n( +-1 ) = 0 for every n >= 2, so EVERY MODE VANISHES ON THE AXIS
 *     IDENTICALLY. The flat side of the half-disc needs no separate treatment
 *     in the exterior at all -- it is satisfied by construction rather than
 *     imposed.
 *   * The exterior modes rho^( 1 - n ) all DECAY, and there is no admissible
 *     constant mode. So the coupling paper's undetermined constant u_infinity
 *     and its compatibility condition -- its two fiddliest pieces -- do not
 *     arise here. The axisymmetric problem is CLEANER than the Laplace one it
 *     is modelled on, which is not the usual direction.
 *   * The map is diagonal: if psi|_Gamma = sum a_n C_n( mu ) then
 *     psi_ext = sum a_n ( rho/rhoGamma )^( 1 - n ) C_n( mu ) and
 *
 *         d psi/d rho |_Gamma = sum a_n ( 1 - n )/rhoGamma C_n( mu ).
 *
 * AND THE ORTHOGONALITY WEIGHT IS MEQ'S OWN, which is the second accident.
 * Gegenbauer functions of order -1/2 are orthogonal in the weight
 * ( 1 - mu^2 )^-1, and on a semicircle centred on the axis
 *
 *     dGamma / r = rhoGamma d theta / ( rhoGamma sin theta ) = dmu/( 1 - mu^2 )
 *
 * -- exactly the weight the Grad-Shafranov weak form already carries, with
 * rhoGamma cancelling out of it. So
 *
 *     integral_Gamma C_m C_n dGamma/r = delta_mn h_n,
 *     h_n = 2 / ( n( n - 1 )( 2n - 1 ) ),
 *
 * and the exterior contributes a diagonal block with entry
 * ( n - 1 ) h_n / rhoGamma.
 *
 * MFEM-FREE, DELIBERATELY, like Profiles, Source, Zernike and SurfaceFit: it is
 * special functions and geometry, plain doubles in and out. That is what lets
 * CI -- which cannot obtain the MFEM branch MEQ needs -- build and test it, and
 * it is why the acceptance for FB-0 is a unit test rather than a convergence
 * study.
 *
 * WHAT THIS CLASS IS NOT. It is not a solver and it does not know about a mesh.
 * It supplies the four numbers a coupled system needs -- the basis on Gamma,
 * the symbol, the mass, and the field at an exterior point -- and FB-1 is where
 * they are assembled into one. It also assumes Gamma is a SEMICIRCLE centred on
 * the axis; the diagonality is a property of that geometry and of no other, so
 * an elongated or offset boundary is a different (and dense) operator.
 */

namespace meq
{

	/// The exterior Dirichlet-to-Neumann map of a semicircle centred on the
	/// axis, in the Gegenbauer basis where it is diagonal.
	///
	/// MODES ARE INDEXED BY THEIR DEGREE n, STARTING AT 2, NOT BY POSITION.
	/// n = 0 and n = 1 are inadmissible rather than merely omitted: h_n divides
	/// by n( n - 1 ), and the functions they would name do not vanish on the
	/// axis. Every method here therefore takes an n in [ 2, 2 + modeCount() ),
	/// and a caller looping over modes should loop over degrees. Getting this
	/// wrong is an off-by-two rather than an off-by-one, which is the kind that
	/// produces a plausible wrong answer rather than an obvious one.
	///
	/// AND THE FORMULA ITSELF IS ONLY VALID FOR n >= 2, which is a separate
	/// reason not to extend the loop downward "to see". Verified: the true
	/// Gegenbauer functions are C_0 = 1 and C_1 = -mu, while
	/// ( P_{n-2} - P_n )/( 2n - 1 ) returns 1 - mu for BOTH -- so the formula
	/// does not merely divide by zero there, it silently returns the wrong
	/// function. What n = 0 and n = 1 would name is a four-dimensional
	/// degenerate sector -- 1, z/rho, rho and z, all Delta*-harmonic -- and NOT
	/// ONE OF THEM VANISHES ON THE AXIS. That is why they are inadmissible here
	/// rather than merely inconvenient, and it is a stronger statement than the
	/// plan's: the two that also decay, 1 and z/rho, have INFINITE norm in the
	/// weight dGamma/r, so there is nothing for a compatibility condition to be
	/// imposed on. Nothing admissible is lost by starting at 2 -- { C_n }_{n>=2}
	/// is complete in L^2( dmu/( 1 - mu^2 ) ).
	class ExteriorDtN
	{
		public:
			/// @param zCentreIn   the axial position of the semicircle's centre.
			/// @param rhoGammaIn  its radius. Must be positive.
			/// @param modesIn     how many modes, so degrees 2 .. modesIn + 1.
			///                    Must be at least one.
			///
			/// Throws std::invalid_argument on a non-positive radius or a
			/// non-positive mode count, rather than producing a basis that is
			/// quietly empty.
			ExteriorDtN( double zCentreIn, double rhoGammaIn, int modesIn );

			/// C_n( cos theta ) at the point ( r, z ), where theta is measured
			/// from the centre.
			///
			/// The point need NOT lie on Gamma: only the direction matters,
			/// since C_n is a function of mu alone. That is deliberate -- it is
			/// what lets exterior() evaluate the field anywhere outside without
			/// a second angular routine -- but it means an accidental call with
			/// an interior point returns a number rather than an error.
			double basis( int n, double r, double z ) const;

			/// dC_n/dmu at the same point, for a caller assembling a tangential
			/// derivative on Gamma.
			double basisDerivative( int n, double r, double z ) const;

			/// ( 1 - n )/rhoGamma: the Dirichlet-to-Neumann symbol.
			///
			/// THE SIGN IS THE OUTWARD RADIAL DERIVATIVE, d psi/d rho, with rho
			/// increasing AWAY from the centre and so out of the computational
			/// domain. It is negative for every admissible n, which is the
			/// statement that an exterior mode decays. A weak form that wants
			/// the derivative along its own outward normal wants this number as
			/// it stands; one written with an inward normal wants it negated,
			/// and FREE-BOUNDARY-PLAN.md section 7 is explicit that this sign
			/// will be got wrong at least once and that FB-1 is what settles it.
			double symbol( int n ) const;

			/// h_n = 2/( n( n - 1 )( 2n - 1 ) ): the mass of mode n in the
			/// weight dGamma/r, which is the weight the weak form carries.
			///
			/// Note it does NOT depend on rhoGamma. That cancellation is not a
			/// simplification made here -- it is a property of the semicircle,
			/// and it is what makes the exterior block a function of the mode
			/// alone.
			double mass( int n ) const;

			/// ( n - 1 ) h_n / rhoGamma: the diagonal entry the coupled system
			/// actually receives, which is -symbol( n ) * mass( n ).
			///
			/// Supplied as its own method rather than left to the caller
			/// because the sign is where the mistake goes: this one is POSITIVE
			/// for every admissible n.
			double blockEntry( int n ) const;

			/// psi at an exterior point, from the coefficients.
			///
			/// @param a  one coefficient per mode, in degree order from n = 2.
			///           Its size must be modeCount().
			///
			/// For output, and for the field at the coils. Evaluating it INSIDE
			/// Gamma is meaningless -- the expansion is of the exterior solution
			/// and the interior one is what the mesh is for -- so a point closer
			/// to the centre than rhoGamma throws rather than extrapolating a
			/// decaying mode inward, where rho^( 1 - n ) grows without bound.
			double exterior( double r, double z,
			                 std::vector<double> const &a ) const;

			/// The coefficients of a trace, by projection in the weight
			/// dGamma/r.
			///
			/// @param trace  psi on Gamma, as a function of ( r, z ).
			///
			/// a_n = ( 1/h_n ) integral_Gamma psi C_n dGamma/r, which is exact
			/// for a trace that is itself a combination of the modes and is a
			/// best fit otherwise. THE QUADRATURE IS GAUSS-LEGENDRE IN mu AND
			/// THAT IS NOT AN APPROXIMATION FOR THE ORTHOGONALITY ITSELF: the
			/// singular weight cancels against the factor ( 1 - mu^2 ) that
			/// every C_n carries, leaving a polynomial, so the mass integrals
			/// come out exact. It is an approximation for a general trace, and
			/// quadraturePoints() says how many points are being spent.
			std::vector<double> coefficients(
				std::function<double( double, double )> const &trace ) const;

			/// How many modes: degrees 2 .. modeCount() + 1.
			int modeCount() const;
			/// The lowest admissible degree, which is 2. A named constant
			/// rather than a literal at every call site.
			static int firstMode();
			/// The highest degree this object carries.
			int lastMode() const;

			double zCentre() const;
			double rhoGamma() const;

			/// How many Gauss-Legendre points coefficients() spends. Fixed at
			/// construction from the mode count, since the integrand's degree
			/// grows with it.
			int quadraturePoints() const;

		private:
			/// Throws unless firstMode() <= n <= lastMode(), naming the range.
			void requireMode( int n ) const;

			/// mu = cos( theta ) = ( z - zCentre )/rho at ( r, z ), and rho.
			/// One place, because the sign of ( z - zCentre ) and the choice of
			/// which coordinate is the polar axis are exactly what a second
			/// copy would get wrong.
			void direction( double r, double z, double &mu, double &rho ) const;

			double zCentreValue;
			double rhoGammaValue;
			int modeCountValue;
			int quadraturePointCount;
	};

}

#endif // MEQ_EXTERIORDTN_HPP
