#ifndef MEQ_TESTS_VARYINGCENTRIFUGAL_HPP
#define MEQ_TESTS_VARYINGCENTRIFUGAL_HPP

/*
 * A rotating plasma whose CENTRIFUGAL EXPONENT VARIES ACROSS FLUX SURFACES, and
 * an independent closed form for the source it produces.
 *
 * Source of the equations, and nothing else is used:
 *   refs/RotatingGK.pdf  (96), the poloidal density variation; (97), the
 *                  quasineutrality condition that closes it; (136), the
 *                  generalised Grad-Shafranov equation -- Abel, Plunk, Wang,
 *                  Barnes, Cowley, Dorland & Schekochihin, Rep. Prog. Phys. 76
 *                  (2013) 116201
 *   docs/rotation.rst  the derivation, the gauge, and the units
 *
 * WHAT THIS IS FOR. At two species the exponent of (96) is the SAME for both,
 *
 *     n_s( r, psi ) = n_s0( psi ) exp[ C( psi ) h( r ) ],
 *     h( r ) := ( r^2 - rRef^2 )/2,
 *     C( psi ) = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ),
 *
 * and C is a flux function. C'( psi ) is non-zero exactly when omega^2/T varies
 * from surface to surface -- and it is the term Li & Zhu (Comput. Phys. Commun.
 * 260 (2021) 107264) print with the wrong sign in their (9), on both the
 * dOmega/dpsi and the dT/dpsi corrections. Neither of their own benchmarks can
 * see that, because both hold C constant.
 *
 * NEITHER CAN ANY OTHER CLOSED FORM IN THIS DIRECTORY, AND THAT IS STRUCTURAL
 * RATHER THAN AN OVERSIGHT. RotatingSoloviev.hpp holds T and Omega constant.
 * MaschkePerrin.hpp varies five profiles, but its (4.7) constrains precisely the
 * ratio omega^2/( Rbar T ), so T tracks omega^2 and C comes out CONSTANT -- and
 * C constant is what collapses that paper's (4.6) to its (4.8) and makes the
 * equation solvable at all. A varying C is what makes it unsolvable in closed
 * form. So the route to a benchmark here is a MANUFACTURED one and there is no
 * point looking for another paper.
 *
 * WHAT IS ALREADY COVERED, MEASURED RATHER THAN ASSUMED, because the gap is
 * narrower than "nothing touches C'" and a fixture should say what it adds.
 * RotatingSourceTests.cpp runs at C' != 0 and its
 * thePressureMatchesTheIsothermalClosedForm checks meq::RotatingSource::pressure
 * and densityExponent against an independent closed form there, to 1e-14. So
 * C( psi ) ITSELF is independently pinned. What is not is C' and C'', which
 * enter nothing but dp/dpsi and d2p/dpsi2 -- that is, nothing but f() and
 * dFdPsi(). Those are checked today only by central differences of quantities
 * MEQ itself computes, at 1e-6 to 1e-7, and a central difference of f() cannot
 * see a term missing from both f() and dFdPsi(). This file supplies the closed
 * form those two want, at round-off.
 *
 * RotatingNewtonConvergence.cpp also runs at C' != 0 -- its profiles are
 * unrelated polynomials, so C varies -- but its manufactured source subtracts
 * meq::RotatingSource's own f() at the exact solution, so the exact solution is
 * exact WHATEVER f() computes. It measures dFdPsi against f(); it cannot measure
 * either against the equations.
 *
 * TWO ROUTES TO THE SAME PRESSURE LIVE IN THIS FILE, AND THAT IS THE POINT.
 *
 *   route A   potentialByBisection() solves (97) numerically, by bisection on
 *             Sum_s Z_s n_s = 0, and pressureByBisection() sums n_s T_s at the
 *             root. C APPEARS NOWHERE IN IT.
 *   route B   pressure(), dPressureDPsi() and d2PressureDPsi2() evaluate the
 *             hand-derived closed form through C, C' and C''.
 *
 * Route B is what the driver compares against meq::RotatingSource. Route A is
 * what says route B's algebra is a consequence of (96) and (97) rather than a
 * transcription that happens to agree with MEQ's transcription of the same
 * algebra -- differentiating route A numerically reaches route B's derivatives,
 * and a shared slip in the C-chain would not survive that.
 *
 * THE CONSTRUCTION RUNS BACKWARDS FROM C, WHICH IS THE ONE DESIGN DECISION.
 * Choosing omega and the temperatures and then reading off C leaves C a ratio
 * of whatever they happen to be, and its two derivatives correspondingly
 * awkward to state. Instead C( psi ) is PRESCRIBED as a quadratic and omega is
 * derived from it,
 *
 *     omega( psi ) = sqrt[ C( psi ) ( Z_1 T_2 - Z_2 T_1 )/( Z_1 m_2 - Z_2 m_1 ) ],
 *
 * which is (97) read the other way. C, C' and C'' are then two lines each, and
 * the temperatures stay free -- so the fixture varies T_1, T_2, n_10, n_20,
 * omega AND gg', six profiles, with C' and C'' both non-zero throughout.
 *
 * EVERY SHAPE FUNCTION IS POSITIVE FOR EVERY REAL psi, WHICH IS A REQUIREMENT
 * AND NOT TIDINESS. A Newton iterate is under no obligation to stay in the
 * solution's range, and this fixture is driven through a solve; the closure
 * divides by Z_1 T_2 - Z_2 T_1 and omega is a square root, so a shape function
 * that turns over gives a throw or a NaN from inside a quadrature loop. Each
 * quadratic below has a negative discriminant, which is asserted in the driver
 * rather than left to the reader.
 *
 * UNITS: mu0 = 1, as everywhere else in this directory. The whole fixture is
 * dimensionless; RoPP's own units and the Gaussian-to-SI conversion are
 * meq::RotatingSource's business and are covered by RotatingSourceTests.cpp,
 * which runs in keV and m^-3.
 *
 * GAUGE: phi_0( rRef, psi ) = 0, which is the local gauge meq::RotatingSource
 * takes. n_s0 is therefore the physical density on r = rRef. Comparing tabulated
 * n_s0 against a code in a different gauge is meaningless; see
 * meq::RotatingSource's header.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace meq
{
namespace analytic
{

/**
 * Two species in sonic toroidal rotation, with a centrifugal exponent that
 * varies from flux surface to flux surface.
 *
 * The species are a helium-like ion at Z = +2 and electrons at Z = -1. THE
 * CHARGE ASYMMETRY IS DELIBERATE: every other rotating fixture in this tree runs
 * at Z = +-1, where Z_1 T_2 - Z_2 T_1 collapses to T_1 + T_2 and
 * Z_1 m_2 - Z_2 m_1 to m_1 + m_2, so a code that assumed the sum rather than the
 * charge-weighted combination would pass all of them. Here the combinations are
 * 2 T_2 + T_1 and 2 m_2 + m_1 and the two are told apart.
 *
 * WHAT IS NOT FREE, and it is worth knowing before adding a third shape
 * function. At two species quasineutrality on r = rRef forces
 * n_20 = -( Z_1/Z_2 ) n_10, so the two reference densities are ONE function and
 * the reference pressure is p_0 = n_10( T_1 - ( Z_1/Z_2 ) T_2 ), which is
 * n_10 ( Z_1 T_2 - Z_2 T_1 )/( -Z_2 ). p_0 is therefore proportional to the
 * closure's own denominator whatever profiles are chosen -- a property of two
 * species, not of this fixture, and not something a different choice removes.
 */
class VaryingCentrifugalPlasma
{
	public:
		/// What to do with C' and C'' where they enter dp/dpsi and d2p/dpsi2.
		///
		/// A MUTATION SWITCH IN THE FIXTURE RATHER THAN IN THE PRODUCTION CODE,
		/// so that "breaking the construction moves the answer" is a test case
		/// that runs on every build rather than an experiment somebody once did
		/// with a patched src/meq. The driver asserts that None agrees with
		/// meq::RotatingSource and that the other two do not.
		enum class Mutation
		{
			/// C' and C'' as derived. The only setting a measurement may use.
			None,

			/// C' -> 0 and C'' -> 0: the source a code gets if it differentiates
			/// p_0 and forgets that the exponent is a flux function too. This is
			/// the term that is absent from every published rotating benchmark,
			/// so it is the mistake that survives them all.
			DropExponentDrift,

			/// C' -> -C' and C'' -> -C'': Li & Zhu's (9), whose dOmega/dpsi and
			/// dT/dpsi corrections carry reversed signs. It moves every term
			/// LINEAR in C' or C'' by twice its size and leaves the C'^2 term of
			/// d2p/dpsi2 exactly where it was, that one being quadratic -- so it
			/// is twice DropExponentDrift's perturbation in F, where the only
			/// drift term is linear, and comparable to it in dF/dpsi, where the
			/// quadratic term is a third of the drift. Measured: 1.58 against
			/// 0.79 in F, and 1.02 against 0.87 in dF/dpsi. It is the more
			/// interesting of the two, being a mistake somebody actually made in
			/// print.
			FlipExponentDrift
		};

		/**
		 * @param exponentAmplitude  C_0, the size of the centrifugal exponent.
		 *                           C( psi ) = C_0 kappa( psi ) with kappa a
		 *                           fixed positive quadratic, so this scales the
		 *                           rotation and nothing else. Zero is legal and
		 *                           gives a static plasma.
		 * @param mutationIn         see Mutation. Defaults to None.
		 * @throws std::invalid_argument if @a exponentAmplitude is negative or
		 *         not finite -- a negative C is an inward centrifugal force.
		 */
		explicit VaryingCentrifugalPlasma( double exponentAmplitude = 2.0,
		                                   Mutation mutationIn = Mutation::None )
			: amplitude( exponentAmplitude ), mutationChoice( mutationIn )
		{
			if ( !( amplitude >= 0.0 ) || !std::isfinite( amplitude ) )
			{
				throw std::invalid_argument( "the exponent amplitude must be finite and non-negative" );
			}
		}

		/*
		 * THE CONFIGURATION, AND HOW ITS CONSTANTS WERE FIXED.
		 *
		 * C_0 = 2 puts the shared exponent C h over [ -1.41, +2.11 ] on the
		 * standard benchmark box with rRef at its middle, so the densities vary
		 * by a factor of 34 across it -- sonic, which is the regime RoPP (136) is
		 * about, and not a perturbation a dropped term could hide inside.
		 *
		 * MEASURED on the standard box [0.6,1.4] x [-0.6,0.6] with psi swept over
		 * [ -0.2, 1.0 ], which is the range the manufactured solve visits:
		 *
		 *   C                      1.76 to  4.40
		 *   C'                     1.00 to  3.40      never zero, zero at -0.7
		 *   C''                    2.00                exactly, kappa being quadratic
		 *   C h                   -1.41 to  2.11      so n varies by a factor of 33
		 *
		 * and the C'-carrying share is measured by the driver and printed: at the
		 * outboard edge it is 63 to 68% of dp/dpsi and 69 to 89% of d2p/dpsi2
		 * over psi in [ 0, 1 ]. Those are the whole argument for the fixture -- a
		 * term worth most of the Jacobian, checked today only by differencing the
		 * code that computes it.
		 */

		/// The standard configuration: C_0 = 2, sonic.
		static VaryingCentrifugalPlasma standard()
		{
			return VaryingCentrifugalPlasma( 2.0 );
		}

		/// C_0 = 0, so omega is identically zero and the exponent with it. THE
		/// CONTROL FOR THE WHOLE FIXTURE: p collapses to p_0( psi ), F to
		/// mu0 r^2 p_0' + g g', and anything the rotating case shows that this
		/// one shows too belongs to the profiles rather than to the rotation.
		static VaryingCentrifugalPlasma stationary()
		{
			return VaryingCentrifugalPlasma( 0.0 );
		}

		/// The same plasma with C' and C'' mutated; see Mutation.
		VaryingCentrifugalPlasma mutated( Mutation how ) const
		{
			return VaryingCentrifugalPlasma( amplitude, how );
		}

		// ---- constants ----

		static std::size_t speciesCount()
		{
			return 2;
		}

		/// m_s. Helium-4 in proton masses, and the electron/proton ratio: KEPT
		/// rather than sent to zero, because meq::RotatingSource keeps it and a
		/// fixture should exercise what the code keeps.
		double mass( std::size_t s ) const
		{
			require( s );
			return s == 0 ? 4.0 : 1.0/1836.152673;
		}

		/// Z_s, signed and in units of e.
		double charge( std::size_t s ) const
		{
			require( s );
			return s == 0 ? 2.0 : -1.0;
		}

		/// rRef, where phi_0 vanishes and n_s0 is the physical density. The
		/// middle of the standard box, so that h( r ) changes sign inside it --
		/// which matters, because a term proportional to h is invisible to a
		/// check made only where h > 0.
		static double referenceRadius()
		{
			return 1.0;
		}

		/// mu0. One: this fixture is dimensionless.
		static double mu0()
		{
			return 1.0;
		}

		/// C_0.
		double exponentAmplitude() const
		{
			return amplitude;
		}

		Mutation mutation() const
		{
			return mutationChoice;
		}

		// ---- the six flux functions, each exact at two derivative levels ----

		/// T_s( psi ). T_1 is affine and rising; T_2 is a quadratic that falls
		/// and turns over, so the two are not proportional and
		/// Z_1 T_2 - Z_2 T_1 is not a multiple of either. Both are positive for
		/// every real psi except T_1 below -1/0.35 = -2.857, which no iterate
		/// reaches.
		double temperature( std::size_t s, double psi ) const
		{
			require( s );
			return s == 0 ? 0.9*( 1.0 + 0.35*psi )
			              : 0.6*( 1.0 - 0.5*psi + 0.2*psi*psi );
		}

		double temperaturePrime( std::size_t s, double psi ) const
		{
			require( s );
			return s == 0 ? 0.9*0.35 : 0.6*( -0.5 + 0.4*psi );
		}

		double temperatureDoublePrime( std::size_t s, double ) const
		{
			require( s );
			return s == 0 ? 0.0 : 0.6*0.4;
		}

		/// n_0, the density amplitude. IT IS THE REACTION-RATIO KNOB AND NOTHING
		/// ELSE, which is why it is a named constant rather than folded into the
		/// shape below.
		///
		/// C is PRESCRIBED here and omega follows from it, so scaling the
		/// densities scales p_0, hence F and dF/dpsi, and leaves C, C' and C''
		/// exactly where they were -- the whole centrifugal structure this
		/// fixture exists to test is untouched by it. What it does move is
		/// max| dF/dpsi |/lambda_1, the diagnostic CLAUDE_HDGGS.md uses for how hard a
		/// source is on Newton, with lambda_1 = pi^2( 1/w^2 + 1/h^2 ) = 22.28 on
		/// the standard box.
		///
		/// MEASURED: at n_0 = 1 the ratio is 29.3, which is deep into the
		/// indefinite regime -- the linearised operator -Delta* - dF/dpsi has
		/// many negative eigenvalues there and a manufactured solve on it is a
		/// study of Newton's basin rather than of its Jacobian. At n_0 = 0.02 it
		/// is 0.597 over the range the solve visits, beside
		/// RotatingNewtonConvergence.cpp's 0.36, so Newton takes a handful of
		/// steps and what is measured is its ORDER.
		static double densityAmplitude()
		{
			return 0.02;
		}

		/// n_s0( psi ), the density on r = rRef. The shape is a quadratic with a
		/// negative discriminant, 0.64 - 3.6, so it is positive for every real
		/// psi; n_20 is -( Z_1/Z_2 ) n_10 = 2 n_10, which is quasineutrality on
		/// the reference curve, exactly and at every derivative level.
		double referenceDensity( std::size_t s, double psi ) const
		{
			return densityWeight( s )*densityAmplitude()*( 1.5 + 0.8*psi + 0.6*psi*psi );
		}

		double referenceDensityPrime( std::size_t s, double psi ) const
		{
			return densityWeight( s )*densityAmplitude()*( 0.8 + 1.2*psi );
		}

		double referenceDensityDoublePrime( std::size_t s, double ) const
		{
			return densityWeight( s )*densityAmplitude()*1.2;
		}

		/// omega( psi ), DERIVED from the prescribed C by (97) read backwards:
		/// omega^2 = C D/K with D = Z_1 T_2 - Z_2 T_1 and K = Z_1 m_2 - Z_2 m_1.
		/// Both C and D are positive for every real psi, so the square root is
		/// always real and no branch is needed.
		double omega( double psi ) const
		{
			return std::sqrt( omegaSquared( psi ) );
		}

		double omegaPrime( double psi ) const
		{
			double const w = omegaSquared( psi );
			if ( w == 0.0 )
			{
				// C_0 = 0. omega is identically zero, so its derivative is zero
				// and the 1/( 2 sqrt w ) below would be 0/0. This is not a
				// removable singularity dodged -- it is the exact derivative of
				// the constant zero function.
				return 0.0;
			}
			return 0.5*omegaSquaredPrime( psi )/std::sqrt( w );
		}

		double omegaDoublePrime( double psi ) const
		{
			double const w = omegaSquared( psi );
			if ( w == 0.0 )
			{
				return 0.0;
			}
			double const wPrime = omegaSquaredPrime( psi );
			return 0.5*omegaSquaredDoublePrime( psi )/std::sqrt( w )
			       - 0.25*wPrime*wPrime/( w*std::sqrt( w ) );
		}

		/// g dg/dpsi, affine. NON-CONSTANT deliberately: RotatingSoloviev.hpp,
		/// MaschkePerrin.hpp and RotatingNewtonConvergence.cpp all carry a
		/// constant g g', so ( g g' )' -- the one term of dF/dpsi that is not the
		/// pressure chain rule -- is exercised by no rotating fixture in the
		/// tree.
		///
		/// WHERE IT DOMINATES AND WHERE IT DOES NOT, because "small" would be a
		/// half-truth. The pressure term carries mu0 r^2 exp( C h ), which spans
		/// four orders across the box: at r = 1.4, psi = 1 it is 13.4 against
		/// this profile's ( g g' )' = 0.25, so g g' is 1.8% of dF/dpsi and the
		/// drift share quoted above is undiluted; at r = 0.6, psi = 0 the
		/// exponential damps it to 1.0e-3 and g g' is essentially all of dF/dpsi.
		/// That is the equation's own shape rather than a choice -- a centrifugal
		/// term is small on the inboard side -- and it is why the drift shares
		/// are reported at the outboard edge and the pointwise agreement is
		/// asserted everywhere.
		static double ggPrime( double psi )
		{
			return -0.4 + 0.25*psi;
		}

		static double ggPrimeDerivative( double )
		{
			return 0.25;
		}

		// ---- the closure of (96) and (97), in closed form ----

		/// h( r ) = ( r^2 - rRef^2 )/2, the radial factor every exponent carries.
		static double radialFactor( double r )
		{
			return 0.5*( r*r - referenceRadius()*referenceRadius() );
		}

		/// D( psi ) = Z_1 T_2 - Z_2 T_1, the closure's denominator. Positive for
		/// every real psi: 2.1 - 0.285 psi + 0.24 psi^2, discriminant -1.93.
		double closureDenominator( double psi ) const
		{
			return charge( 0 )*temperature( 1, psi ) - charge( 1 )*temperature( 0, psi );
		}

		double closureDenominatorPrime( double psi ) const
		{
			return charge( 0 )*temperaturePrime( 1, psi )
			       - charge( 1 )*temperaturePrime( 0, psi );
		}

		double closureDenominatorDoublePrime( double psi ) const
		{
			return charge( 0 )*temperatureDoublePrime( 1, psi )
			       - charge( 1 )*temperatureDoublePrime( 0, psi );
		}

		/// K = Z_1 m_2 - Z_2 m_1, the closure's numerator constant.
		double closureMassFactor() const
		{
			return charge( 0 )*mass( 1 ) - charge( 1 )*mass( 0 );
		}

		/// C( psi ) = C_0 kappa( psi ), kappa = 1 + 0.7 psi + 0.5 psi^2.
		///
		/// PRESCRIBED, not derived: omega is what follows from it. kappa's
		/// discriminant is 0.49 - 2 = -1.51, so kappa is positive for every real
		/// psi with a minimum of 0.755, which is what makes omega real
		/// everywhere.
		double exponentCoefficient( double psi ) const
		{
			return amplitude*( 1.0 + 0.7*psi + 0.5*psi*psi );
		}

		/// C'( psi ). THE TERM THIS FIXTURE EXISTS FOR.
		double exponentCoefficientPrime( double psi ) const
		{
			return amplitude*( 0.7 + 1.0*psi );
		}

		/// C''( psi ). Constant and non-zero, kappa being a genuine quadratic --
		/// a linear kappa would leave C'' = 0 and one of the three terms of
		/// d2p/dpsi2 untested.
		double exponentCoefficientDoublePrime( double ) const
		{
			return amplitude*1.0;
		}

		/// The shared exponent A( r, psi ) = C( psi ) h( r ) of (96). At two
		/// species BOTH species carry it, which is what makes Sum_s Z_s n_s
		/// vanish at every r once it vanishes at rRef.
		double densityExponent( double r, double psi ) const
		{
			return exponentCoefficient( psi )*radialFactor( r );
		}

		/// e phi_0( r, psi ), the potential that holds (97).
		///
		/// Derived here rather than quoted: eliminating the two exponents of (96)
		/// against Z_1 n_10 = -Z_2 n_20 gives
		/// omega^2 h ( m_1 T_2 - m_2 T_1 ) = y ( Z_1 T_2 - Z_2 T_1 ), which is
		/// linear in y. Zero at r = rRef exactly, by h( rRef ) = 0, which is the
		/// gauge.
		double potential( double r, double psi ) const
		{
			double const num = mass( 0 )*temperature( 1, psi )
			                   - mass( 1 )*temperature( 0, psi );
			return omegaSquared( psi )*radialFactor( r )*num/closureDenominator( psi );
		}

		/// d( e phi_0 )/dpsi at fixed r. Exposed because meq::RotatingSource
		/// exposes dPotentialDPsi() and the driver checks it: it is a second,
		/// independent place where a wrong chain rule through omega and T shows
		/// up, and unlike f() it isolates the potential from the pressure.
		double potentialPrime( double r, double psi ) const
		{
			double const num = mass( 0 )*temperature( 1, psi )
			                   - mass( 1 )*temperature( 0, psi );
			double const numPrime = mass( 0 )*temperaturePrime( 1, psi )
			                        - mass( 1 )*temperaturePrime( 0, psi );
			double const d = closureDenominator( psi );
			double const dPrime = closureDenominatorPrime( psi );

			double const ratio = num/d;
			double const ratioPrime = numPrime/d - num*dPrime/( d*d );

			return radialFactor( r )*( omegaSquaredPrime( psi )*ratio
			                           + omegaSquared( psi )*ratioPrime );
		}

		/// n_s( r, psi ) from (96).
		double density( std::size_t s, double r, double psi ) const
		{
			require( s );
			return referenceDensity( s, psi )*std::exp( densityExponent( r, psi ) );
		}

		// ---- the pressure and the source: route B ----

		/// p_0( psi ) = Sum_s n_s0 T_s, the pressure on the reference curve.
		double referencePressure( double psi ) const
		{
			double sum = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				sum += referenceDensity( s, psi )*temperature( s, psi );
			}
			return sum;
		}

		double referencePressurePrime( double psi ) const
		{
			double sum = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				sum += referenceDensityPrime( s, psi )*temperature( s, psi )
				       + referenceDensity( s, psi )*temperaturePrime( s, psi );
			}
			return sum;
		}

		double referencePressureDoublePrime( double psi ) const
		{
			double sum = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				sum += referenceDensityDoublePrime( s, psi )*temperature( s, psi )
				       + 2.0*referenceDensityPrime( s, psi )*temperaturePrime( s, psi )
				       + referenceDensity( s, psi )*temperatureDoublePrime( s, psi );
			}
			return sum;
		}

		/// p( r, psi ) = p_0( psi ) exp[ C( psi ) h( r ) ]. Both species carry
		/// the same exponent, so it comes outside the sum -- which is the only
		/// reason a two-species rotating pressure has a closed form at all.
		double pressure( double r, double psi ) const
		{
			return referencePressure( psi )*std::exp( densityExponent( r, psi ) );
		}

		/// dp/dpsi at FIXED r:
		///
		///     [ p_0' + p_0 C' h ] exp( C h ).
		///
		/// The second term is the whole of what a varying centrifugal exponent
		/// contributes to F, and it vanishes identically on r = rRef, where
		/// h = 0. That is the trap CLAUDE_FLOW.md records: the one radius where
		/// the gauge is exact is the one radius where a rotating source is
		/// indistinguishable from a static one, so a check placed there measures
		/// nothing.
		double dPressureDPsi( double r, double psi ) const
		{
			double const h = radialFactor( r );
			double const cPrime = effectiveExponentPrime( psi );

			return ( referencePressurePrime( psi )
			         + referencePressure( psi )*cPrime*h )
			       *std::exp( densityExponent( r, psi ) );
		}

		/// d2p/dpsi2 at FIXED r:
		///
		///     [ p_0'' + 2 p_0' C' h + p_0 ( C'' h + C'^2 h^2 ) ] exp( C h ).
		///
		/// THREE of the four terms carry C' or C''. That is the ratio the
		/// mutation controls exploit and it is why dropping the drift is not a
		/// small perturbation of the Jacobian.
		double d2PressureDPsi2( double r, double psi ) const
		{
			double const h = radialFactor( r );
			double const cPrime = effectiveExponentPrime( psi );
			double const cDoublePrime = effectiveExponentDoublePrime( psi );

			double const p0 = referencePressure( psi );

			return ( referencePressureDoublePrime( psi )
			         + 2.0*referencePressurePrime( psi )*cPrime*h
			         + p0*( cDoublePrime*h + cPrime*cPrime*h*h ) )
			       *std::exp( densityExponent( r, psi ) );
		}

		/// F = mu0 r^2 dp/dpsi|_r + g g', which is (136) collapsed. Returns F,
		/// NOT F/r: the 1/r belongs to the weak form, as everywhere else in this
		/// directory and in meq::Source.
		double f( double r, double /*z*/, double psi ) const
		{
			return mu0()*r*r*dPressureDPsi( r, psi ) + ggPrime( psi );
		}

		/// dF/dpsi = mu0 r^2 d2p/dpsi2 + ( g g' )'.
		double dFdPsi( double r, double /*z*/, double psi ) const
		{
			return mu0()*r*r*d2PressureDPsi2( r, psi ) + ggPrimeDerivative( psi );
		}

		// ---- (97) solved numerically: route A ----

		/// e phi_0 by BISECTION on Sum_s Z_s n_s = 0, with no C anywhere.
		///
		/// WHY BISECTION AND NOT NEWTON. meq::RotatingSource's general closure is
		/// a safeguarded Newton with an implicit-differentiation derivative; the
		/// point of this routine is to be a different algorithm reaching the same
		/// number, so it is the plainest bracketing method there is and it takes
		/// no derivative at all. Sum_s Z_s n_s is strictly decreasing in y --
		/// every term of its derivative carries Z_s^2 -- and runs from +infinity
		/// to -infinity when the charges have both signs, so a bracket always
		/// exists and bisection always converges.
		///
		/// One hundred halvings takes any starting bracket below the granularity
		/// of a double, so the loop bound is a fixed count rather than a
		/// tolerance: a deterministic number of flops, and no branch that could
		/// behave differently on two builds.
		double potentialByBisection( double r, double psi ) const
		{
			double scale = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				scale = std::max( scale, std::abs( temperature( s, psi ) ) );
			}
			if ( !( scale > 0.0 ) )
			{
				throw std::invalid_argument( "every temperature vanishes, so (97) has no scale" );
			}

			// The residual is decreasing, so it is positive at the low end.
			double lo = -scale;
			double hi = scale;
			for ( int i = 0; i < 200 && !( neutralityResidual( lo, r, psi ) > 0.0 ); ++i )
			{
				lo *= 2.0;
			}
			for ( int i = 0; i < 200 && !( neutralityResidual( hi, r, psi ) < 0.0 ); ++i )
			{
				hi *= 2.0;
			}

			for ( int i = 0; i < 100; ++i )
			{
				double const mid = 0.5*( lo + hi );
				if ( neutralityResidual( mid, r, psi ) > 0.0 )
				{
					lo = mid;
				}
				else
				{
					hi = mid;
				}
			}

			return 0.5*( lo + hi );
		}

		/// p = Sum_s n_s T_s with n_s built from (96) at the bisected potential.
		/// The independent anchor for everything above: it never forms C, never
		/// forms p_0, and never takes a psi-derivative.
		double pressureByBisection( double r, double psi ) const
		{
			double const y = potentialByBisection( r, psi );
			double const h = radialFactor( r );
			double const w = omegaSquared( psi );

			double sum = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				double const t = temperature( s, psi );
				double const exponent = mass( s )*w*h/t - charge( s )*y/t;
				sum += referenceDensity( s, psi )*std::exp( exponent )*t;
			}
			return sum;
		}

		/// Sum_s Z_s n_s at a trial potential, which is (97)'s left hand side.
		/// Exposed so the driver can assert it really does vanish at the bisected
		/// root rather than merely that the bisection returned.
		double neutralityResidual( double y, double r, double psi ) const
		{
			double const h = radialFactor( r );
			double const w = omegaSquared( psi );

			// Factored by the largest exponent, so that a high Mach number cannot
			// overflow the residual. The common factor is positive and cancels
			// out of the sign test, which is all bisection uses.
			double top = -std::numeric_limits<double>::infinity();
			double exponents[ 2 ] = { 0.0, 0.0 };
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				double const t = temperature( s, psi );
				exponents[ s ] = mass( s )*w*h/t - charge( s )*y/t;
				top = std::max( top, exponents[ s ] );
			}

			double sum = 0.0;
			for ( std::size_t s = 0; s < speciesCount(); ++s )
			{
				sum += charge( s )*referenceDensity( s, psi )
				       *std::exp( exponents[ s ] - top );
			}
			return sum;
		}

	private:
		static void require( std::size_t s )
		{
			if ( s >= speciesCount() )
			{
				throw std::out_of_range( "no such species" );
			}
		}

		/// The charge weight that makes quasineutrality on r = rRef exact:
		/// n_20 = -( Z_1/Z_2 ) n_10, so the weights are 1 and 2 here.
		double densityWeight( std::size_t s ) const
		{
			require( s );
			return s == 0 ? 1.0 : -charge( 0 )/charge( 1 );
		}

		/// omega^2 = C D/K, which is what makes C the prescribed function.
		double omegaSquared( double psi ) const
		{
			return exponentCoefficient( psi )*closureDenominator( psi )/closureMassFactor();
		}

		double omegaSquaredPrime( double psi ) const
		{
			return ( exponentCoefficientPrime( psi )*closureDenominator( psi )
			         + exponentCoefficient( psi )*closureDenominatorPrime( psi ) )
			       /closureMassFactor();
		}

		double omegaSquaredDoublePrime( double psi ) const
		{
			return ( exponentCoefficientDoublePrime( psi )*closureDenominator( psi )
			         + 2.0*exponentCoefficientPrime( psi )*closureDenominatorPrime( psi )
			         + exponentCoefficient( psi )*closureDenominatorDoublePrime( psi ) )
			       /closureMassFactor();
		}

		/// C' as the pressure derivatives see it, which is where the mutation
		/// lives. pressure() itself is NOT mutated: a mutation that moved p as
		/// well would be caught by the pressure check before it reached the
		/// Jacobian, and the control is supposed to isolate the drift term.
		double effectiveExponentPrime( double psi ) const
		{
			return driftScale()*exponentCoefficientPrime( psi );
		}

		double effectiveExponentDoublePrime( double psi ) const
		{
			return driftScale()*exponentCoefficientDoublePrime( psi );
		}

		double driftScale() const
		{
			switch ( mutationChoice )
			{
				case Mutation::DropExponentDrift:
					return 0.0;
				case Mutation::FlipExponentDrift:
					return -1.0;
				case Mutation::None:
				default:
					return 1.0;
			}
		}

		double amplitude;
		Mutation mutationChoice;
};

}
}

#endif // MEQ_TESTS_VARYINGCENTRIFUGAL_HPP
