#define BOOST_TEST_MODULE RotatingNewtonConvergence

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * FL-5 OF FLOW-PLAN.md: THE MANUFACTURED NON-LINEAR ROTATING CASE, AND THE
 * FIRST THING IN THE FLOW WORK THAT CAN SEE A WRONG JACOBIAN.
 *
 * Everything FL-0 to FL-4 measured is blind to meq::RotatingSource::dFdPsi.
 * Li & Zhu's rotating Solov'ev has T and Omega constant and p linear in psi, so
 * dF/dpsi is identically zero; so does Maschke & Perrin's, for the same reason
 * -- their (4.7) forces the exponent coefficient C to be a constant. FLOW-PLAN
 * §6.3 says it plainly: neither published rotating benchmark can see the
 * C'( psi ) term, and C' is precisely where the chain rule through omega and
 * the temperatures lives. A wrong dFdPsi leaves every error and every rate in
 * those studies unchanged to six figures -- CLAUDE.md measures exactly that for
 * a Jacobian 5% out -- and only the OBSERVED NEWTON ORDER moves, to 1.
 *
 * THE CONSTRUCTION, which is the whole point of this file.
 *
 * A rotating source's F depends on ( r, psi ) only, so it cannot in general
 * equal -Delta* psi for a psi somebody chose. The manufactured trick is to add
 * a remainder that carries NO psi dependence at all:
 *
 *     F_total( r, z, psi ) = rot.f( r, z, psi ) + h( r, z )
 *     h( r, z )            = -Delta*( psiExact )( r, z )
 *                            - rot.f( r, z, psiExact( r, z ) )
 *
 * At psi = psiExact the two rotating terms cancel and F_total = -Delta* psiExact,
 * so psiExact IS an exact solution of -Delta* psi = F_total. And because h is a
 * function of ( r, z ) alone,
 *
 *     dF_total/dpsi = rot.dFdPsi
 *
 * EXACTLY. THE ENTIRE JACOBIAN IS THE ROTATING SOURCE'S. That is what makes
 * this a test of meq::RotatingSource::dFdPsi rather than of a manufactured
 * expression somebody wrote beside it: there is nothing else in the Jacobian to
 * be right. It is the same device Example 5 uses -- its bracket
 * ( sc^2 - psi^2 + exp( -sc ) - exp( -psi ) ) vanishes at the exact solution --
 * except that here the surviving psi-dependence is not written out at all, it
 * is whatever meq::RotatingSource computes.
 *
 * psiExact is ManufacturedNonlinear's shape, which keeps this comparable with
 * NewtonConvergence.cpp's Example 5 study on the same box:
 *
 *     psi( r, z ) = sin( kr ( r + r0 ) ) cos( kz z ),
 *     r0 = -0.5,  kr = 1.15 pi,  kz = 1.15
 *
 * whose Delta* is closed form. Differentiating,
 *
 *     d_rr psi = -kr^2 psi,   d_zz psi = -kz^2 psi,
 *     ( 1/r ) d_r psi = ( kr/r ) cos( kr ( r + r0 ) ) cos( kz z ),
 *
 * and Delta* := d_rr - ( 1/r ) d_r + d_zz, so
 *
 *     Delta* psi = -( kr^2 + kz^2 ) psi
 *                  - ( kr/r ) cos( kr ( r + r0 ) ) cos( kz z ).
 *
 * That is checked against a central difference in the FIRST test below rather
 * than trusted, because a mistyped Delta* would make every rate in this file
 * measure the wrong equation and no rate could say so. Measured worst
 * disagreement over the box: 4.786e-07, which is the O( h^2 ) truncation floor
 * of the difference and not a discrepancy.
 *
 * THE PROFILES CARRY GENUINE psi-DEPENDENCE IN ALL THREE OF n_s0, T_s AND
 * omega, and that is a requirement rather than a flourish: with any one of them
 * constant the chain rule this file exists to check loses a term, and with all
 * three constant dF/dpsi is zero and the file tests nothing. The second test
 * case is the control that says so.
 *
 * TWO CONSTRAINTS ON CHOOSING THEM, both of which bite because psiExact runs
 * over roughly [ -0.11, +1.0 ] on the standard box and NOT over [ 0, 1 ]:
 *
 *   * every temperature must stay strictly POSITIVE over the whole visited
 *     range. The two-species closure divides by Z_1 T_2 - Z_2 T_1, which for
 *     opposite charges is T_1 + T_2 and goes singular when a temperature turns
 *     over; meq::RotatingSource throws from inside a quadrature loop when it
 *     does, which is a poor place to find out. The ranges are printed and the
 *     positivity asserted, so a later change to the profiles fails loudly.
 *
 *   * charge neutrality on r = rRef must hold at EVERY psi the solve reaches.
 *     meq::neutralisingDensity() makes that automatic and exact at all three
 *     derivative levels, which is why the electron density is built with it
 *     rather than hand-balanced -- the constructor only samples [ 0, 1 ].
 *
 * NORMALISED UNITS, mu0 = 1, as RotatingSourceConvergence.cpp does. The
 * benchmark is dimensionless; putting SI constants in and dividing them out
 * again would only add round-off.
 *
 * The reaction ratio max| dF/dpsi |/lambda_1 -- the diagnostic CLAUDE.md uses
 * for how hard a source is on Newton, with lambda_1 = pi^2( 1/w^2 + 1/h^2 ) =
 * 22.3 on this box -- comes out at 0.36 over the solution's actual range and
 * 0.76 over the whole ( r, psi ) box. Real, and mild enough that what gets
 * measured is Newton's ORDER rather than the edge of its basin: it takes three
 * steps on every mesh and every degree studied below.
 *
 * AND THE FILE WAS MUTATION TESTED, because a test of a Jacobian that cannot
 * see a wrong one is worth nothing and this suite's own stance says to check
 * rather than assume. Returning 1.05*rot.dFdPsi() from dFdPsi() below, and
 * changing nothing else:
 *
 *   every L2 error and every rate     UNCHANGED, to all seven digits printed,
 *   in all three tables               at k = 1, 2 and 3
 *   Newton iterations                 3 -> 6 on every mesh
 *   observed Newton order             1.980 -> 1.055, settling to 1.000
 *   the Jacobian check                3.5e-11 -> 2.1e-04
 *
 * That is CLAUDE.md's +5% experiment reproduced on the rotating source, and it
 * is the whole argument for this file existing beside FL-4's rate study: the
 * rate study is the part that does not move.
 */

namespace
{

	using meq::tests::Measurement;
	using meq::tests::bestNewtonOrder;
	using meq::tests::newtonOrder;
	using meq::tests::standardBox;

	// ManufacturedNonlinear::example5()'s parameters, so that this study sits on
	// the same flux shape and the same box as NewtonConvergence.cpp's. r0 is an
	// offset that places the arch of the sine inside the domain and is NOT a
	// major radius.
	double const pi = 3.14159265358979323846;
	double const psiR0 = -0.5;
	double const psiKr = 1.15*pi;
	double const psiKz = 1.15;

	// The gauge curve. phi_0 vanishes here and each n_s0 is the physical density
	// here; the box is r in [ 0.6, 1.4 ], so this is its middle and the exponent
	// ( r^2 - rRef^2 )/2 runs over [ -0.32, +0.48 ].
	double const referenceRadius = 1.0;

	/// A polynomial in psi, exact at all three derivative levels.
	///
	/// meq::Profile grew doublePrime() for FLOW-PLAN §5.3's reason: p is not a
	/// flux function, F is already dp/dpsi, and the Jacobian spends a second
	/// derivative of every input. A polynomial is used here rather than
	/// meq::HermiteCubicSpline deliberately -- a Hermite cubic is C^1, so its
	/// second derivative is piecewise linear and JUMPS at every interior knot,
	/// and a Newton step landing on a knot would see a one-sided Jacobian. That
	/// is a real effect nobody has measured; it has no business being mixed into
	/// a measurement of the chain rule.
	class PolynomialProfile : public meq::Profile
	{
		public:
			explicit PolynomialProfile( std::vector<double> coefficients )
				: c( std::move( coefficients ) )
			{
			}

			double operator()( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = c.size(); i > 0; --i )
					sum = sum*psi + c[ i - 1 ];

				return sum;
			}

			double prime( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = c.size(); i > 1; --i )
					sum = sum*psi + static_cast<double>( i - 1 )*c[ i - 1 ];

				return sum;
			}

			double doublePrime( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = c.size(); i > 2; --i )
					sum = sum*psi + static_cast<double>( ( i - 1 )*( i - 2 ) )*c[ i - 1 ];

				return sum;
			}

		private:
			std::vector<double> c;
	};

	std::shared_ptr<meq::Profile const> polynomial( std::vector<double> coefficients )
	{
		return std::make_shared<PolynomialProfile const>( std::move( coefficients ) );
	}

	/*
	 * THE PROFILE SET. Every one of the four flux functions below varies with
	 * psi, and each is chosen so that its own contribution to the chain rule is
	 * non-zero over the range the solve visits:
	 *
	 *   T_i   = 1.0 + 0.30 psi      positive for psi > -3.33
	 *   T_e   = 0.8 - 0.20 psi      positive for psi <  4.00
	 *   n_i0  = 2.0 + 1.00 psi + 0.50 psi^2
	 *   n_e0  = neutralisingDensity, hence exactly n_i0 for Z = +-1
	 *   omega = 1.2 + 0.30 psi
	 *
	 * n_i0 is 0.5( psi + 1 )^2 + 1.5 and so is positive for every real psi,
	 * which matters because a Newton iterate is under no obligation to stay in
	 * the solution's range. The temperatures have the margins printed by the
	 * control test below.
	 *
	 * The masses are hydrogenic: m_i = 1 and m_e = 1e-4, kept rather than sent
	 * to zero because meq::RotatingSource keeps them and this file should
	 * exercise what it keeps. Only the combinations Z_1 m_2 - Z_2 m_1 and
	 * m_1 T_2 - m_2 T_1 ever appear.
	 *
	 * g g' IS CONSTANT, AND THAT IS DELIBERATE. It contributes to F, so the term
	 * is exercised, but ( g g' )' is then exactly zero -- which means EVERY LAST
	 * BIT of dF_total/dpsi comes from mu0 r^2 d2p/dpsi2, the rotating pressure
	 * chain rule. meq::MHDSource already covers a psi-dependent g g', and mixing
	 * it in here would only give the Jacobian somewhere else to be right.
	 */
	double const ionMass = 1.0;
	double const electronMass = 1.0e-4;
	double const ggPrimeConstant = -0.4;

	std::vector<meq::Species> rotatingSpecies()
	{
		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = ionMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = polynomial( { 1.0, 0.30 } );
		species[ 0 ].density = polynomial( { 2.0, 1.00, 0.50 } );

		species[ 1 ].mass = electronMass;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = polynomial( { 0.8, -0.20 } );
		// Charge neutrality on r = rRef, closed exactly rather than balanced by
		// hand: the constructor samples psi over [ 0, 1 ] only, and this solve
		// visits negative psi.
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return species;
	}

	meq::RotatingSource rotatingSource()
	{
		return meq::RotatingSource( rotatingSpecies(),
			polynomial( { 1.2, 0.30 } ),
			polynomial( { ggPrimeConstant } ),
			referenceRadius,
			1.0 );
	}

	/**
	 * The manufactured non-linear rotating benchmark.
	 *
	 * It is a meq::Source, so it can be handed straight to
	 * GradShafranovSolver::setSource(); and it carries psi(), flux() and
	 * deltaStarFD() as the fixtures in tests/analytic do, so
	 * meq::tests::checkOrder() can measure it against its own exact solution.
	 *
	 * The whole psi-dependence of f() -- and therefore the whole of dFdPsi() --
	 * belongs to the meq::RotatingSource held by value below. The remainder
	 * added to it is a function of ( r, z ) alone; see the file header.
	 */
	class ManufacturedRotatingSource : public meq::Source
	{
		public:
			ManufacturedRotatingSource()
				: rot( rotatingSource() )
			{
			}

			/// The imposed poloidal flux.
			double psi( double r, double z ) const
			{
				return std::sin( psiKr*( r + psiR0 ) )*std::cos( psiKz*z );
			}

			/// grad_bar( psi ). Not the HDG flux: that is this divided by r.
			void gradPsi( double r, double z, double & dPsiDr, double & dPsiDz ) const
			{
				dPsiDr =  psiKr*std::cos( psiKr*( r + psiR0 ) )*std::cos( psiKz*z );
				dPsiDz = -psiKz*std::sin( psiKr*( r + psiR0 ) )*std::sin( psiKz*z );
			}

			/// The HDG flux q = grad_bar( psi )/r.
			void flux( double r, double z, double & qR, double & qZ ) const
			{
				gradPsi( r, z, qR, qZ );
				qR /= r;
				qZ /= r;
			}

			/// Delta*( psiExact ), in closed form. Derived in the file header and
			/// checked against deltaStarFD() by the first test case, because this
			/// is the one expression a mistake in would silently change which
			/// equation every rate below is measuring.
			double deltaStar( double r, double z ) const
			{
				return -( psiKr*psiKr + psiKz*psiKz )*psi( r, z )
				       - ( psiKr/r )*std::cos( psiKr*( r + psiR0 ) )*std::cos( psiKz*z );
			}

			/// Delta*( psiExact ) by central differences, arranged as
			/// r d_r( ( 1/r ) d_r psi ) + d_zz psi -- the same form
			/// Soloviev.hpp and ManufacturedNonlinear.hpp use, so that it is
			/// independent of the closed form above rather than a rearrangement
			/// of it.
			double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
			{
				auto innerR = [ & ]( double rr )
				{
					return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
				};

				double const dRInner = ( innerR( r + h ) - innerR( r - h ) )/( 2.0*h );
				double const dZZ = ( psi( r, z + h ) - 2.0*psi( r, z ) + psi( r, z - h ) )
				                   /( h*h );

				return r*dRInner + dZZ;
			}

			/// The remainder that makes psiExact solve the equation. A function of
			/// ( r, z ) ALONE, which is the entire reason dFdPsi() below can be
			/// the rotating source's own derivative with nothing added to it.
			double remainder( double r, double z ) const
			{
				return -deltaStar( r, z ) - rot.f( r, z, psi( r, z ) );
			}

			double f( double r, double z, double psiValue ) const override
			{
				return rot.f( r, z, psiValue ) + remainder( r, z );
			}

			/// EXACTLY meq::RotatingSource::dFdPsi. Nothing is added, because the
			/// remainder does not depend on psi.
			double dFdPsi( double r, double z, double psiValue ) const override
			{
				return rot.dFdPsi( r, z, psiValue );
			}

			meq::RotatingSource const & rotating() const
			{
				return rot;
			}

		private:
			meq::RotatingSource rot;
	};

	/// The range of psiExact over the box, which is what every profile has to be
	/// well behaved on. Sampled rather than reasoned about: the solve visits
	/// whatever is there.
	void psiRange( double & lowest, double & highest )
	{
		meq::tests::Rectangle const box = standardBox();
		ManufacturedRotatingSource const eq;

		lowest = 1.0e300;
		highest = -1.0e300;
		for ( int i = 0; i <= 80; ++i )
		{
			for ( int j = 0; j <= 80; ++j )
			{
				double const r = box.rMin + box.width()*i/80.0;
				double const z = box.zMin + box.height()*j/80.0;
				lowest = std::min( lowest, eq.psi( r, z ) );
				highest = std::max( highest, eq.psi( r, z ) );
			}
		}
	}

}

/*
 * THE GUARD, AND IT COMES FIRST BECAUSE EVERYTHING ELSE RESTS ON IT.
 *
 * Two independent statements, in the order a failure should be read:
 *
 *   1. the hand-derived Delta* agrees with a central difference of psiExact, so
 *      the closed form in the header is transcribed correctly;
 *   2. Delta*( psiExact ) = -F_total( r, z, psiExact ), so psiExact really is
 *      the solution of the equation the solver is about to be given.
 *
 * The second is the construction of the whole file -- the remainder cancelling
 * the rotating source at the exact solution -- and it is checked rather than
 * argued. A wrong sign convention converges, at the right rate, to the wrong
 * function; only this can see that.
 */
BOOST_AUTO_TEST_CASE( theManufacturedSolutionSatisfiesTheEquation )
{
	ManufacturedRotatingSource const eq;
	meq::tests::Rectangle const box = standardBox();

	double worstClosedForm = 0.0;
	double worstEquation = 0.0;

	for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.05 )
	{
		for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.075 )
		{
			double const fd = eq.deltaStarFD( r, z );
			double const closed = eq.deltaStar( r, z );
			double const minusF = -eq.f( r, z, eq.psi( r, z ) );

			worstClosedForm = std::max( worstClosedForm, std::fabs( closed - fd ) );
			worstEquation = std::max( worstEquation, std::fabs( fd - minusF ) );

			BOOST_TEST( std::fabs( closed - fd ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): the closed-form Delta* gives "
			            << closed << " where a central difference gives " << fd
			            << ", so the derivation in this file's header is mistyped and "
			            "every rate below is measuring the wrong equation" );
			BOOST_TEST( std::fabs( fd - minusF ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): Delta*( psi ) = " << fd
			            << " but -F_total = " << minusF << ", so the manufactured "
			            "remainder does not close the equation at the exact solution" );
		}
	}

	// Measured: 4.786e-07 and 4.786e-07, both the O( h^2 ) truncation floor of
	// the difference at h = 1e-4 rather than a disagreement. That the two agree to
	// every digit is the construction working rather than a coincidence: at the
	// exact solution -F_total IS the closed-form Delta*, by the definition of the
	// remainder, so both columns are measuring the same truncation error.
	std::printf( "\n  the manufactured rotating solution\n"
	             "    closed-form Delta* vs central difference: worst %.3e\n"
	             "    Delta*( psi ) vs -F_total:                worst %.3e\n",
	             worstClosedForm, worstEquation );
	std::fflush( stdout );
}

/*
 * THE CONTROL. If dF/dpsi came out near zero this whole file would pass and
 * test nothing at all -- which is exactly the position FL-4 is in, and the
 * reason FL-5 exists.
 *
 * Two things are asserted, and they are different:
 *
 *   * dF/dpsi is comfortably NON-ZERO, so there is a Jacobian mass term to get
 *     wrong at all;
 *   * dF/dpsi VARIES with psi, because a constant one would make this
 *     McCarthy's rung of the ladder -- affine source, one exact Newton step --
 *     rather than a new one.
 *
 * It also prints and asserts on the ranges every profile has to survive. The
 * temperatures are the ones that matter: the two-species closure divides by
 * Z_1 T_2 - Z_2 T_1, which for opposite charges is T_1 + T_2, and
 * meq::RotatingSource throws when that vanishes -- from inside a quadrature
 * loop, where it is a poor thing to debug. Asserting positivity here makes a
 * later change to the profiles fail with a sentence instead.
 */
BOOST_AUTO_TEST_CASE( theSourceIsGenuinelyNonlinearInPsi )
{
	ManufacturedRotatingSource const eq;
	meq::tests::Rectangle const box = standardBox();

	double lowest = 0.0, highest = 0.0;
	psiRange( lowest, highest );

	// A margin either side, because a Newton iterate is under no obligation to
	// stay inside the solution's own range.
	double const margin = 0.2;
	double const psiLow = lowest - margin;
	double const psiHigh = highest + margin;

	std::vector<meq::Species> const & species = eq.rotating().species();

	double worstTemperature = 1.0e300;
	double coldestSum = 1.0e300;
	double tMin[ 2 ] = { 1.0e300, 1.0e300 };
	double tMax[ 2 ] = { -1.0e300, -1.0e300 };

	for ( int i = 0; i <= 200; ++i )
	{
		double const psi = psiLow + ( psiHigh - psiLow )*i/200.0;

		double sum = 0.0;
		for ( std::size_t s = 0; s < species.size(); ++s )
		{
			double const t = ( *species[ s ].temperature )( psi );
			tMin[ s ] = std::min( tMin[ s ], t );
			tMax[ s ] = std::max( tMax[ s ], t );
			worstTemperature = std::min( worstTemperature, t );
			sum += t;
		}
		coldestSum = std::min( coldestSum, sum );
	}

	std::printf( "\n  the visited ranges\n"
	             "    psi over the box            [ %+.6f, %+.6f ]\n"
	             "    psi sampled with margin     [ %+.6f, %+.6f ]\n"
	             "    T_i over that range         [ %+.6f, %+.6f ]\n"
	             "    T_e over that range         [ %+.6f, %+.6f ]\n"
	             "    T_i + T_e, the closure's denominator, at least %+.6f\n",
	             lowest, highest, psiLow, psiHigh,
	             tMin[ 0 ], tMax[ 0 ], tMin[ 1 ], tMax[ 1 ], coldestSum );
	std::fflush( stdout );

	BOOST_TEST( worstTemperature > 0.0,
	            "a species temperature falls to " << worstTemperature << " over the "
	            "range this solve visits. A non-positive temperature has no "
	            "Maxwellian, and the closure's denominator Z_1 T_2 - Z_2 T_1 goes "
	            "singular before it gets there -- meq::RotatingSource will throw "
	            "from inside a quadrature loop" );
	BOOST_TEST( coldestSum > 0.1,
	            "Z_1 T_2 - Z_2 T_1 falls to " << coldestSum << ", which is close "
	            "enough to zero that the closed-form phi_0 is being evaluated near "
	            "its pole" );

	// dF/dpsi over the ( r, psi ) box, and separately over the solution's own
	// range, which is the one the reaction ratio should be read from -- sampling
	// a range the solve never visits is a mistake CLAUDE.md records making once
	// already.
	double smallest = 1.0e300;
	double largest = -1.0e300;
	double onSolution = 0.0;

	for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.02 )
	{
		for ( int i = 0; i <= 100; ++i )
		{
			double const psi = psiLow + ( psiHigh - psiLow )*i/100.0;
			double const d = eq.dFdPsi( r, 0.0, psi );
			smallest = std::min( smallest, d );
			largest = std::max( largest, d );
		}

		for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.05 )
			onSolution = std::max( onSolution, std::fabs( eq.dFdPsi( r, z, eq.psi( r, z ) ) ) );
	}

	// lambda_1 = pi^2( 1/w^2 + 1/h^2 ) for the Dirichlet Laplacian on this box:
	// the scale CLAUDE.md's reaction-ratio diagnostic is read against.
	double const lambdaOne = pi*pi*( 1.0/( box.width()*box.width() )
	                                 + 1.0/( box.height()*box.height() ) );

	std::printf( "    dF/dpsi over the ( r, psi ) box   [ %+.6f, %+.6f ]\n"
	             "    max | dF/dpsi | on the solution    %.6f\n"
	             "    lambda_1 = %.4f, so the reaction ratio is %.4f on the solution "
	             "and %.4f over the box\n",
	             smallest, largest, onSolution, lambdaOne,
	             onSolution/lambdaOne,
	             std::max( std::fabs( smallest ), std::fabs( largest ) )/lambdaOne );
	std::fflush( stdout );

	// Measured: dF/dpsi runs over [ +0.272, +16.985 ] on the ( r, psi ) box and
	// reaches 7.913 on the solution itself, against lambda_1 = 22.275 -- ratios
	// of 0.76 and 0.36. Real, and mild enough that the measurement below is of
	// Newton's order rather than of the edge of its basin.
	BOOST_TEST( onSolution > 1.0,
	            "max | dF/dpsi | on the solution is only " << onSolution
	            << ". The Jacobian's mass term is then negligible beside the "
	            "operator, this file's Newton order would be 2 whatever dFdPsi "
	            "returned, and FL-5 would be measuring nothing" );

	// And it must VARY, or this is McCarthy's rung and not a new one.
	BOOST_TEST( largest - smallest > 1.0,
	            "dF/dpsi varies by only " << largest - smallest << " over the whole "
	            "( r, psi ) box, so the source is effectively affine in psi and one "
	            "exact Newton step would finish it -- which is McCarthy.hpp's "
	            "statement, already made in NewtonConvergence.cpp" );

	/*
	 * AND THE VARIATION MUST BE THERE IN psi AT FIXED r, not merely inherited
	 * from the r^2 prefactor every source in meq carries.
	 *
	 * WHERE THAT IS MEASURED MATTERS, AND r = rRef IS EXACTLY THE WRONG PLACE.
	 * The gauge pins phi_0( rRef ) = 0, so the exponent ( r^2 - rRef^2 )/2 and
	 * both its psi-derivatives vanish there and d2p/dpsi2 collapses to P0''( psi )
	 * -- the answer a NON-rotating source would give. Measured on the reference
	 * curve dF/dpsi moves 0.333 across the solution's range, against 7.300 at the
	 * outboard edge, so a control placed on rRef would be blind to the entire
	 * rotation chain rule while looking like a reasonable check. It was placed
	 * there first, and read as a failure of the profiles rather than of the
	 * probe.
	 */
	double const atRefLow = eq.dFdPsi( referenceRadius, 0.0, lowest );
	double const atRefHigh = eq.dFdPsi( referenceRadius, 0.0, highest );
	double const atEdgeLow = eq.dFdPsi( box.rMax, 0.0, lowest );
	double const atEdgeHigh = eq.dFdPsi( box.rMax, 0.0, highest );

	std::printf( "    across psi = %+.4f to %+.4f, dF/dpsi moves %.6f at r = rRef "
	             "(where the exponent vanishes) and %.6f at r = rMax\n",
	             lowest, highest, std::fabs( atRefHigh - atRefLow ),
	             std::fabs( atEdgeHigh - atEdgeLow ) );
	std::fflush( stdout );

	BOOST_TEST( std::fabs( atEdgeHigh - atEdgeLow ) > 1.0,
	            "at fixed r = rMax, dF/dpsi moves only "
	            << std::fabs( atEdgeHigh - atEdgeLow ) << " across the solution's "
	            "range, so the psi-dependence of n_s0, T_s and omega is not "
	            "reaching the Jacobian" );
}

/*
 * THE STRONGEST CHECK IN THE FILE: THE ASSEMBLED JACOBIAN AGAINST A CENTRAL
 * DIFFERENCE OF THE ASSEMBLED RESIDUAL.
 *
 * This is NewtonConvergence.cpp's jacobianMatchesAFiniteDifferenceOfTheResidual
 * on the rotating source, and it is here for the reason that file gives: it
 * sees a dFdPsi that is wrong in a way a pointwise finite difference of f()
 * would not catch. RotatingSourceTests.cpp already differences f() -- what it
 * cannot see is the term reaching the operator wrongly, the hybridized
 * elimination, or the essential trace condition, all of which are in this one.
 *
 * The perturbation is zero on the essential trace dofs deliberately.
 * DarcyHybridization carries the Dirichlet condition the way a NonlinearForm
 * does: the values ride in the iterate, the residual is masked to zero on those
 * rows and the Jacobian is given a unit row. A difference quotient of a masked
 * row is zero while the Jacobian row is the identity, so probing those
 * directions would report a discrepancy that is not one.
 *
 * WHAT IS CHECKED IS THE CONDENSED OPERATOR'S JACOBIAN, which is the backup
 * ordering's and not the one meq's default NPC path drives. That is not a gap
 * this file can close and NewtonConvergence.cpp explains why: NPC's Jacobian
 * handle is solve-only, since after ComputeH() the local blocks are factored in
 * place, so there is no action to difference against. What stands in for it is
 * the assertion on observed Newton order in the next test case.
 *
 * The linearisation is established at `state` BEFORE any residual evaluation
 * and held there for the whole difference: the first Mult() on a freshly formed
 * system is what sets the local state the gradient is taken at, and differencing
 * without pinning it first measures a gradient taken at a third place. On
 * Example 5 that reads 4e-8 relative against 4e-11 at a fixed linearisation,
 * and none of it is the Jacobian's fault.
 */
BOOST_AUTO_TEST_CASE( theAssembledJacobianIsTheDerivativeOfTheAssembledResidual )
{
	ManufacturedRotatingSource const source;

	mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), 3 );
	mfem::FunctionCoefficient psiCoeff( [ &source ]( mfem::Vector const & x )
	{
		return source.psi( x( 0 ), x( 1 ) );
	} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( psiCoeff );
	solver.prepare();

	mfem::Operator & residual = solver.reducedOperator();
	int const size = solver.reducedSolution().Size();

	std::vector<bool> isEssential( size, false );
	mfem::Array<int> const & essential = solver.essentialTraceDofs();
	for ( int i = 0; i < essential.Size(); ++i )
		isEssential[ essential[ i ] ] = true;

	BOOST_TEST( essential.Size() > 0,
	            "no essential trace dofs: the Dirichlet condition is not being "
	            "imposed, which would make this check vacuous and the solver wrong" );

	// Away from the solution and away from zero, so that the rotating source's
	// exponential is genuinely being exercised. The wiggle is a deterministic
	// function of the dof index rather than a random draw: a test that fails only
	// sometimes is worse than no test.
	mfem::Vector state( solver.reducedSolution() );
	for ( int i = 0; i < size; ++i )
		if ( !isEssential[ i ] )
			state( i ) += 0.3*std::sin( 1.7*i + 0.4 );

	residual.GetGradient( state );

	double const step = 1.0e-5;
	int const numDirections = 4;

	std::vector<mfem::Vector> directions( numDirections );
	std::vector<mfem::Vector> differences( numDirections );

	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector & direction = directions[ trial ];
		direction.SetSize( size );
		for ( int i = 0; i < size; ++i )
			direction( i ) = isEssential[ i ]
			                 ? 0.0
			                 : std::cos( 2.3*i + 1.1*trial ) + 0.5*std::sin( 0.7*i - 0.3*trial );
		direction /= direction.Norml2();

		mfem::Vector forward( size );
		mfem::Vector backward( size );
		mfem::Vector shifted( state );

		shifted.Add( step, direction );
		residual.Mult( shifted, forward );

		shifted = state;
		shifted.Add( -step, direction );
		residual.Mult( shifted, backward );

		// A central difference: its truncation error is O( step^2 ), which keeps
		// the comparison well clear of the discrepancy a wrong Jacobian gives.
		mfem::Vector & difference = differences[ trial ];
		difference.SetSize( size );
		difference = forward;
		difference -= backward;
		difference /= 2.0*step;
	}

	mfem::Operator & jacobian = residual.GetGradient( state );

	double worst = 0.0;
	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector applied( size );
		jacobian.Mult( directions[ trial ], applied );

		mfem::Vector error( differences[ trial ] );
		error -= applied;

		double const relative = error.Norml2()/applied.Norml2();
		std::printf( "  rotating Jacobian vs central difference, direction %d: "
		             "|| J d || = %.6e, relative error %.3e\n",
		             trial, applied.Norml2(), relative );
		worst = std::max( worst, relative );
	}
	std::fflush( stdout );

	BOOST_TEST( worst < 1.0e-8,
	            "the assembled Jacobian differs from a central difference of the "
	            "assembled residual by " << worst << " relative. Since the "
	            "manufactured remainder carries no psi-dependence, ALL of this "
	            "Jacobian is meq::RotatingSource::dFdPsi -- so this is that chain "
	            "rule, through n_s0, T_s and omega. Newton will still converge, but "
	            "linearly, and no rate in this file would say so" );
}

/*
 * NEWTON'S OBSERVED ORDER, WHICH IS THE INSTRUMENT THAT SEES A WRONG dFdPsi
 * UNDER THE DEFAULT NPC ORDERING.
 *
 * The order is log( r2/r1 )/log( r1/r0 ): 2 for a quadratically convergent
 * iteration in its asymptotic regime, 1 for a linearly convergent one anywhere.
 * CLAUDE.md's measurement is the calibration -- a dF/dpsi 5% too large leaves
 * every error and every rate in this file unchanged to six significant figures
 * and drops the observed order to exactly 1.000.
 *
 * Per the testing stance, the assertion is on the BEST triple above the
 * round-off floor rather than on every one: Newton is not monotone on a source
 * like this, and only a monotone tail supports an order claim. The whole history
 * is printed so that a reader can see which triple it came from.
 */
BOOST_AUTO_TEST_CASE( newtonConvergesAtOrderTwo )
{
	ManufacturedRotatingSource const eq;

	std::vector<double> history;
	Measurement const point = meq::tests::measure( eq, standardBox(), 3, 8, &history );

	BOOST_TEST_REQUIRE( history.size() >= 3u,
	                    "Newton produced only " << history.size() << " residuals, so "
	                    "there is nothing to estimate an order from" );

	std::printf( "\n  Newton residual history, rotating manufactured case, "
	             "k = 3, h = %.5f, %d trace dofs\n", point.h, point.traceDofs );
	std::printf( "  %5s %16s %16s %8s\n", "it", "||r||", "||r||/||r_0||", "order" );
	for ( std::size_t i = 0; i < history.size(); ++i )
	{
		if ( i >= 2 )
		{
			std::printf( "  %5zu %16.6e %16.6e %8.3f\n", i, history[ i ],
			             history[ i ]/history[ 0 ],
			             newtonOrder( history[ i - 2 ], history[ i - 1 ], history[ i ] ) );
		}
		else
		{
			std::printf( "  %5zu %16.6e %16.6e %8s\n", i, history[ i ],
			             history[ i ]/history[ 0 ], "-" );
		}
	}
	std::fflush( stdout );

	double const best = bestNewtonOrder( history );

	BOOST_TEST( best >= 1.8,
	            "the best observed Newton order is " << best << ", not the ~2 a "
	            "Jacobian consistent with its residual gives. Every psi-dependent "
	            "term in this problem comes from meq::RotatingSource, so a run that "
	            "grinds down linearly means its dFdPsi and its f disagree -- and "
	            "the rates below would not move by a digit" );

	// AND BOUNDED ABOVE, which matters because this is the BEST over the triples
	// rather than the tail. This history is short and strictly monotone, so the
	// two agree here -- but RotatingNormalisedConvergence's is not, and there the
	// best triple reads 3.81 out of an iterate still walking into the basin,
	// which would pass a one-sided test for entirely the wrong reason. Anything
	// much above 2 is a triple straddling the round-off floor, not a better
	// Jacobian.
	BOOST_TEST( best <= 2.5,
	            "the best observed Newton order is " << best << ", above 2. That is not "
	            "a better Jacobian: it means the history has stopped being a clean "
	            "quadratic tail and the number is an artefact" );

	// The problem must actually be non-linear here, or the order above is
	// measuring nothing. Two residuals would mean a single Newton step sufficed,
	// which is what a Solov'ev or a Li & Zhu rotating source gives and is exactly
	// why FL-4 cannot stand in for this.
	BOOST_TEST( point.newtonIterations >= 2,
	            "Newton converged in " << point.newtonIterations << " iterations, so "
	            "the source is behaving linearly in psi and dFdPsi is untested" );

	// A quadratically convergent iteration reaching 1e-12 relative takes a
	// handful of steps, not a dozen.
	BOOST_TEST( point.newtonIterations <= 8,
	            "Newton took " << point.newtonIterations << " iterations to reach "
	            "the tolerance; the reference shape is four or five" );
}

/*
 * AND THE RATES, k+1 IN psi AND IN q OVER FOUR DYADIC MESHES.
 *
 * These say the right equation is being solved, which is a different and
 * independent statement from the Newton order above -- CLAUDE.md's whole point
 * is that the two fail separately. The ceilings sit at roughly three times the
 * measured finest-mesh error, as SolovievConvergence.cpp does, because a rate is
 * blind to a solution wrong by a constant factor or a sign and the ceiling is
 * the only assertion here that can see that.
 */

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Measured at h = 0.025: psi 4.3142e-04, q 7.4841e-04, rates 2.007 and 1.995.
	meq::tests::checkOrder( ManufacturedRotatingSource(),
	                        "manufactured rotating, two species", 1, 1.3e-3, 2.3e-3 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured at h = 0.025: psi 3.5796e-06, q 6.6266e-06, rates 2.999 and 2.998.
	meq::tests::checkOrder( ManufacturedRotatingSource(),
	                        "manufactured rotating, two species", 2, 1.1e-5, 2.0e-5 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured at h = 0.025: psi 2.2440e-08, q 4.0629e-08, rates 4.002 and 3.999.
	meq::tests::checkOrder( ManufacturedRotatingSource(),
	                        "manufactured rotating, two species", 3, 6.7e-8, 1.2e-7 );
}
