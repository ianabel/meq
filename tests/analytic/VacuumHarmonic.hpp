#ifndef MEQ_TESTS_VACUUMHARMONIC_HPP
#define MEQ_TESTS_VACUUMHARMONIC_HPP

/*
 * Delta*-harmonic flux functions that VANISH ON THE AXIS, for FB-A.
 *
 * FREE-BOUNDARY-PLAN.md section 7 puts FB-A first -- "the axis. A vacuum solve
 * on a mesh touching R = 0. No free boundary, no coupling" -- and ROADMAP.md
 * calls it the one thing measurable today. Its acceptance was written as "the
 * trace condition number bounded under refinement", and that is weaker than
 * this project accepts anywhere else. It does not have to be: a vacuum solve
 * needs Delta* psi = 0, and there are POLYNOMIAL Delta*-harmonic functions that
 * are exactly zero on the axis, so FB-A gets a rate against a closed form like
 * every other stage.
 *
 * MEQ's operator is Delta* psi = d_rr psi - ( 1/R ) d_r psi + d_zz psi, and
 * -Delta* psi = F, so a vacuum solve is F = 0. The three functions here:
 *
 *     psi = R^2                R^2 z                R^4 - 4 R^2 z^2
 *
 * Each is Delta*-harmonic and each is identically zero at R = 0. Checked
 * symbolically before being written down, and checked again at run time by
 * deltaStarFD() -- see the note on the fourth candidate below.
 *
 * THE CONTROL, AND IT IS THE INTERESTING ONE. Cerfon & Freidberg's twelve-term
 * basis also contains R^2 ln R - z^2, which is Delta*-harmonic and BOUNDED at
 * the axis but NOT zero there -- it tends to -z^2. axisNonZero() is that
 * function. A solve that handles R^2 and fails on this one has found something
 * about the axis rather than about the mesh, which is the discrimination
 * FREE-BOUNDARY-PLAN.md section 8's three non-measurements are missing. Note
 * that its GRADIENT carries a ln R and a 1/R, so q is unbounded at the axis and
 * only the potential is: it is a control, not a benchmark, and a mesh reaching
 * exactly R = 0 cannot integrate it.
 *
 * A FOURTH CANDIDATE WAS GUESSED AND WAS WRONG, which is why every one of these
 * is checked rather than asserted. R^2 ( R^2 - 4 z^2 ) z looks like the natural
 * odd-in-z partner of R^4 - 4 R^2 z^2 and is not harmonic at all --
 * Delta* of it is -16 R^2 z. The same species of error as the Solov'ev
 * coefficients, and caught the same way: by an independent check rather than by
 * re-reading the algebra.
 *
 * WHY NO SOURCE TERM. These are vacuum fields, so f() and dFdPsi() are
 * identically zero and a Newton solve on them is affine -- one step, exactly,
 * like Soloviev.hpp. That is deliberate: FB-A is about the operator's
 * degeneracy at R = 0 and nothing else, and a non-linearity would put a second
 * variable into a measurement that exists to isolate one.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/// A Delta*-harmonic flux function vanishing on the axis, as a weighted sum of
/// the three admissible polynomials.
///
/// R may be zero for psi() and for the gradient -- that is the entire point --
/// but NOT for a source divided by R, which is why f() returns a bare zero
/// rather than computing anything.
class VacuumHarmonic
{
	public:
		/// @param aIn  the weight on R^2.
		/// @param bIn  the weight on R^2 z.
		/// @param cIn  the weight on R^4 - 4 R^2 z^2.
		VacuumHarmonic( double aIn, double bIn, double cIn )
			: aValue( aIn ), bValue( bIn ), cValue( cIn )
		{
		}

		/// All three modes present, with weights of the same order so that no
		/// one of them dominates the error. A single mode would let a defect in
		/// the treatment of the others hide.
		static VacuumHarmonic mixed()
		{
			return VacuumHarmonic( 1.0, 0.8, 0.5 );
		}

		/// R^2 alone: the simplest field that reaches the axis, and the one to
		/// look at first when the mixed case fails.
		static VacuumHarmonic quadratic()
		{
			return VacuumHarmonic( 1.0, 0.0, 0.0 );
		}

		/// The poloidal flux function. Zero at R = 0 for every weight.
		double psi( double radius, double z ) const
		{
			double const r2 = radius*radius;
			return aValue*r2 + bValue*r2*z + cValue*( r2*r2 - 4.0*r2*z*z );
		}

		/// d psi / d R.
		double dPsiDr( double radius, double z ) const
		{
			return 2.0*aValue*radius + 2.0*bValue*radius*z
			       + cValue*( 4.0*radius*radius*radius - 8.0*radius*z*z );
		}

		/// d psi / d z.
		double dPsiDz( double radius, double z ) const
		{
			return bValue*radius*radius - 8.0*cValue*radius*radius*z;
		}

		/// The flux q = ( 1/R ) grad-bar psi, R component.
		///
		/// BOUNDED AT THE AXIS, AND THAT IS THE POINT. Every term of psi carries
		/// R^2, so every term of d_r psi carries R and the division is exact:
		/// q_r -> 2a as R -> 0 rather than diverging. This is the "q is bounded
		/// at the axis because psi ~ R^2" that FREE-BOUNDARY-PLAN.md section 8
		/// lists as one of its three non-measurements, written down so that it
		/// can be one.
		double fluxR( double radius, double z ) const
		{
			return 2.0*aValue + 2.0*bValue*z
			       + cValue*( 4.0*radius*radius - 8.0*z*z );
		}

		/// The flux q, z component. It VANISHES at the axis rather than merely
		/// being bounded: d_z psi = R^2 ( b - 8 c z ), so q_z = R ( b - 8 c z ).
		double fluxZ( double radius, double z ) const
		{
			return radius*( bValue - 8.0*cValue*z );
		}

		/// The HDG flux q = grad_bar( psi ) / R, in the harness's own signature.
		///
		/// NOT written as gradPsi() then divided by R, which is how every other
		/// fixture here does it -- because at R = 0 that division is 0/0 and the
		/// whole point of this fixture is that the LIMIT exists. fluxR() and
		/// fluxZ() carry the cancellation already done, so this is exact at the
		/// axis rather than a NaN.
		void flux( double radius, double z, double &qR, double &qZ ) const
		{
			qR = fluxR( radius, z );
			qZ = fluxZ( radius, z );
		}

		/// grad_bar( psi ), for a caller that wants it unweighted.
		void gradPsi( double radius, double z, double &dR, double &dZ ) const
		{
			dR = dPsiDr( radius, z );
			dZ = dPsiDz( radius, z );
		}

		/// The source. Identically zero: this is a vacuum field.
		double f( double /*R*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// And so is its derivative, which makes the solve affine.
		double dFdPsi( double /*R*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// Delta*( psi ), by central differences of psi(), in exactly the
		/// arrangement Soloviev.hpp and ManufacturedNonlinear.hpp use.
		///
		/// It must come out zero, and the suite asserts that rather than
		/// trusting the algebra -- see the fourth-candidate note in the file
		/// comment. Keep R well away from zero when calling it: the finite
		/// difference divides by R, so it is a check on the interior and not at
		/// the axis, where the operator's own 1/R is what is being studied.
		double deltaStarFD( double radius, double z, double h = 1.0e-4 ) const
		{
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( radius + h ) - innerR( radius - h ) )/( 2.0*h );
			double const dZZ = ( psi( radius, z + h ) - 2.0*psi( radius, z ) + psi( radius, z - h ) )
			                   /( h*h );

			return radius*dRInner + dZZ;
		}

		/// THE CONTROL: R^2 ln R - z^2. Delta*-harmonic, bounded at the axis,
		/// and NOT zero there -- it tends to -z^2 as R -> 0.
		///
		/// One of Cerfon & Freidberg's twelve terms, and the one the polynomial
		/// ones are usually paired with. Its gradient carries a ln R and its
		/// flux a 1/R, so it is a statement about the POTENTIAL only. Use it on
		/// a mesh whose inner edge is small but positive; a mesh reaching R = 0
		/// exactly cannot integrate it, which is itself the discrimination.
		static double axisNonZero( double radius, double z )
		{
			return radius*radius*std::log( radius ) - z*z;
		}

		/// Delta* of the control, by the same central differences. Also zero.
		static double axisNonZeroDeltaStarFD( double radius, double z,
		                                      double h = 1.0e-4 )
		{
			auto innerR = [ & ]( double rr )
			{
				return ( axisNonZero( rr + h, z ) - axisNonZero( rr - h, z ) )
				       /( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( radius + h ) - innerR( radius - h ) )/( 2.0*h );
			double const dZZ = ( axisNonZero( radius, z + h ) - 2.0*axisNonZero( radius, z )
			                     + axisNonZero( radius, z - h ) )/( h*h );

			return radius*dRInner + dZZ;
		}

		double a() const { return aValue; }
		double b() const { return bValue; }
		double c() const { return cValue; }

	private:
		double aValue;
		double bValue;
		double cValue;
};

}
}

#endif // MEQ_TESTS_VACUUMHARMONIC_HPP
