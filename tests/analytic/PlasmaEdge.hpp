#ifndef MEQ_TESTS_PLASMAEDGE_HPP
#define MEQ_TESTS_PLASMAEDGE_HPP

/*
 * A MANUFACTURED EQUILIBRIUM WHOSE SOURCE STOPS AT A CIRCLE NO MESH LINE
 * FOLLOWS -- the plasma edge of FREE-BOUNDARY-PLAN.md section 5.3, with
 * everything else about free boundary removed.
 *
 * Free boundary makes the source F = [ mu0 R^2 p'(Psi) + (g g')(Psi) ] carry a
 * factor chi_{Omega_p}, and Omega_p is a level set of the solution. So F stops
 * dead on a curve that cuts through elements, and it moves while Newton runs.
 * The question this fixture exists to answer is what that costs the ORDER --
 * MEQ is a k+1 code in psi_h and a k+2 code in psi*, and section 5.3 names this
 * as the one place a published code (CEDRES++) says it hit a wall.
 *
 * THE CONSTRUCTION. Take a disc of radius a about ( R_0, z0 ), strictly inside
 * the benchmark box, as the plasma:
 *
 *     phi( R, z ) = a^2 - ( R - R_0 )^2 - ( z - z0 )^2 ,   plasma = { phi > 0 }
 *
 * and take the exact solution to be a Delta*-HARMONIC vacuum field plus a term
 * supported in the plasma:
 *
 *     psi = w( R, z )  +  c ( phi_+ )^m ,        Delta* w = 0 ,   m = j + 2
 *
 * Then F = -Delta* psi is SUPPORTED IN THE PLASMA AND NOWHERE ELSE, because w
 * contributes nothing, and near the edge
 *
 *     F  ~  -4 c m ( m - 1 ) a^2 phi^( m - 2 )  =  O( phi^j ) ,
 *
 * so j is exactly the order to which the profiles vanish at the plasma edge --
 * p' ~ Psi^j, which is the modelling choice a user actually makes. w is what
 * makes the exterior field non-trivial; without it psi would be identically
 * zero outside the disc and the study would measure nothing there.
 *
 * WHY j IS THE PARAMETER AND NOT SOMETHING ELSE. The three rungs are the three
 * conventions that occur in practice:
 *
 *     j = 0    p'( 0 ) != 0        F JUMPS at the edge      psi in H^2.5
 *     j = 1    p' ~ Psi            F kinks                  psi in H^3.5
 *     j = 2    p' ~ Psi^2          F is C^1                 psi in H^4.5
 *
 * FreeGS's own ( 1 - Psi_n^alpha )^beta gives j = beta, and beta = 1 is its
 * default -- so j = 1 is the case a reader is most likely to meet.
 *
 * The regularity is what caps the achievable rate at m + 1/2 = j + 2.5, and
 * that cap is a property of the EXACT SOLUTION rather than of any method: the
 * singular part is | d |^m across a curve, whose best polynomial approximation
 * on an element of size h is O( h^m ) at every degree, over O( 1/h ) elements.
 * See CLAUDE.md's *The plasma edge caps the order, and it is the profiles that
 * set the cap* for the measurement.
 *
 * NOTHING HERE MOVES. The cut is prescribed, so this isolates the
 * APPROXIMATION and the QUADRATURE from the Jacobian. dFdPsi() is identically
 * zero and a Newton solve is affine -- one step, exactly, like Soloviev.hpp.
 * The moving cut is MovingPlasmaEdge below, which is the same geometry with
 * the support read off psi itself.
 *
 * AND IT CHECKS ITS OWN TRANSCRIPTION. deltaStarFD() recomputes Delta* psi by
 * central differences and the suite asserts it against -f(), which is what
 * caught the fourth VacuumHarmonic candidate and the Solov'ev coefficients.
 * It is only valid AWAY FROM THE EDGE: a central difference straddling the
 * circle differences a function whose m-th derivative jumps, and reports the
 * jump rather than a defect.
 */

#include <cmath>
#include <stdexcept>

namespace meq
{
namespace analytic
{

/// The manufactured equilibrium above, with a FIXED plasma edge.
class PlasmaEdge
{
	public:
		/// @param vanishingOrderIn  j: the order to which the source vanishes at
		///                          the plasma edge. 0, 1 and 2 are the rungs.
		/// @param amplitudeIn       c, the size of the plasma term.
		/// @param harmonicA         weight on R^2 in the vacuum field.
		/// @param harmonicB         weight on R^2 z.
		/// @param harmonicC         weight on R^4 - 4 R^2 z^2.
		PlasmaEdge( int vanishingOrderIn,
		            double centreRIn = 1.0, double centreZIn = 0.0,
		            double radiusIn = 0.23456789, double amplitudeIn = 4.0,
		            double harmonicA = 0.35, double harmonicB = 0.5,
		            double harmonicC = 0.3, double harmonicD = 0.45 )
			: jOrder( vanishingOrderIn ), mOrder( vanishingOrderIn + 2 ),
			  centreR( centreRIn ), centreZ( centreZIn ), discRadius( radiusIn ),
			  amplitude( amplitudeIn ),
			  weightA( harmonicA ), weightB( harmonicB ), weightC( harmonicC ),
			  weightD( harmonicD )
		{
			if ( vanishingOrderIn < 0 )
				throw std::invalid_argument(
					"meq::analytic::PlasmaEdge: the vanishing order must be "
					"non-negative; j = 0 is a source that jumps at the edge" );
		}

		/// phi. Positive inside the plasma, zero on the edge.
		double levelSet( double radius, double z ) const
		{
			double const dr = radius - centreR, dz = z - centreZ;
			return discRadius*discRadius - dr*dr - dz*dz;
		}

		/// The Delta*-harmonic vacuum field, which is psi outside the plasma.
		///
		/// THE LAST TERM IS NOT DECORATION AND THE FIRST VERSION OMITTED IT.
		/// R^2 ln R - z^2 is Delta*-harmonic (VacuumHarmonic::axisNonZero) and
		/// is NOT a polynomial, where the other three are and the highest is
		/// quartic. Without it P_4 represents the whole vacuum field exactly,
		/// so the control in a best-approximation study sits at 1.3e-15 from
		/// the coarsest mesh and its "rate" reads 0.03 -- an instrument
		/// reporting itself, which is the failure this tree keeps meeting.
		/// It is bounded away from the axis on any box this fixture is used on;
		/// its gradient carries a ln R and would not be.
		double vacuum( double radius, double z ) const
		{
			double const r2 = radius*radius;
			return weightA*r2 + weightB*r2*z + weightC*( r2*r2 - 4.0*r2*z*z )
			     + weightD*( r2*std::log( radius ) - z*z );
		}

		double psi( double radius, double z ) const
		{
			double const p = levelSet( radius, z );
			return vacuum( radius, z )
			     + ( p > 0.0 ? amplitude*std::pow( p, mOrder ) : 0.0 );
		}

		/// q = ( 1/R ) grad-bar psi, which is what meq::GradShafranovSolver
		/// solves for and what fluxError() is measured against.
		void flux( double radius, double z, double &qR, double &qZ ) const
		{
			double const r2 = radius*radius;
			double dpsiR = 2.0*weightA*radius + 2.0*weightB*radius*z
			             + weightC*( 4.0*r2*radius - 8.0*radius*z*z )
			             + weightD*( 2.0*radius*std::log( radius ) + radius );
			double dpsiZ = weightB*r2 - 8.0*weightC*r2*z - 2.0*weightD*z;

			double const p = levelSet( radius, z );
			if ( p > 0.0 )
			{
				double const dr = radius - centreR, dz = z - centreZ;
				double const g = amplitude*mOrder*std::pow( p, mOrder - 1 );
				dpsiR += g*( -2.0*dr );
				dpsiZ += g*( -2.0*dz );
			}

			qR = dpsiR/radius;
			qZ = dpsiZ/radius;
		}

		/// F as eq (2) writes it: -Delta* psi. Supported in the plasma alone.
		double f( double radius, double z ) const
		{
			double const p = levelSet( radius, z );
			if ( p <= 0.0 )
				return 0.0;

			double const dr = radius - centreR;
			double const m = static_cast<double>( mOrder );
			// |grad phi|^2 = 4( a^2 - phi ) and Delta* phi = -4 + 2 dr/R.
			double const deltaStarPhi = -4.0 + 2.0*dr/radius;
			double term = m*( m - 1.0 )*std::pow( p, mOrder - 2 )
			              *4.0*( discRadius*discRadius - p );
			term += m*std::pow( p, mOrder - 1 )*deltaStarPhi;
			return -amplitude*term;
		}

		/// meq::Source's spelling. The cut is FIXED here, so psi is ignored and
		/// the problem is affine -- see MovingPlasmaEdge for the other one.
		double f( double radius, double z, double /*psi*/ ) const
		{
			return f( radius, z );
		}

		double dFdPsi( double /*R*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// Delta* psi by central differences, RICHARDSON EXTRAPOLATED, for the
		/// transcription check.
		///
		/// MEANINGLESS WITHIN h OF THE EDGE, where it differences across the
		/// jump in the m-th derivative rather than across a smooth function.
		///
		/// AND THE EXTRAPOLATION IS NOT DECORATION. A plain central difference
		/// carries its own O( h^2 ) truncation and an O( eps/h^2 ) round-off,
		/// which together floor this at about 2e-7 ABSOLUTE however exact f()
		/// is -- and the source's own scale falls by a factor of 24 per rung of
		/// j, so a relative check against that floor passes at j = 0 and fails
		/// at j = 2 while nothing about the fixture has changed. Same finding as
		/// Zernike's derivative check and SurfaceAverage's d/dpsi, for the third
		/// time in this tree.
		double deltaStarFD( double radius, double z, double h = 1.0e-3 ) const
		{
			auto raw = [ & ]( double step )
			{
				double const dRR = ( psi( radius + step, z ) - 2.0*psi( radius, z ) + psi( radius - step, z ) )/( step*step );
				double const dZZ = ( psi( radius, z + step ) - 2.0*psi( radius, z ) + psi( radius, z - step ) )/( step*step );
				double const dR  = ( psi( radius + step, z ) - psi( radius - step, z ) )/( 2.0*step );
				return dRR - dR/radius + dZZ;
			};
			return ( 4.0*raw( 0.5*h ) - raw( h ) )/3.0;
		}

		/// True where deltaStarFD() may be believed: at least @a margin from the
		/// edge, measured in the level set's own units of distance.
		bool awayFromEdge( double radius, double z, double margin ) const
		{
			double const dr = radius - centreR, dz = z - centreZ;
			return std::fabs( std::sqrt( dr*dr + dz*dz ) - discRadius ) > margin;
		}

		int vanishingOrder() const { return jOrder; }
		int kinkOrder()      const { return mOrder; }
		double centre( int component ) const { return component == 0 ? centreR : centreZ; }
		double edgeRadius()  const { return discRadius; }

		/// The cap on the L2 rate this equilibrium admits, in psi_h and in psi*
		/// alike: no polynomial space on a mesh that does not follow the edge
		/// can beat it. m + 1/2 = j + 2.5.
		double rateCap() const { return static_cast<double>( mOrder ) + 0.5; }

	private:
		int jOrder, mOrder;
		double centreR, centreZ, discRadius, amplitude;
		double weightA, weightB, weightC, weightD;
};


/**
 * THE SAME GEOMETRY WITH THE SUPPORT READ OFF psi ITSELF, so the plasma edge
 * MOVES while Newton runs -- which is what free boundary actually does and
 * what PlasmaEdge above deliberately leaves out.
 *
 * The construction is arranged so that the plasma is exactly a level set of
 * the solution.  Take
 *
 *     psi = phi + c ( phi_+ )^m ,
 *
 * with phi the same circular level-set function.  Then psi > 0 exactly where
 * phi > 0, for any c > 0 and m >= 2, so
 *
 *     Omega_p = { psi > 0 } = { phi > 0 }
 *
 * IDENTICALLY -- the discrete support is whatever { psi_h > 0 } happens to be,
 * and it converges to the true circle as psi_h converges.  That is the whole
 * point: nothing tells the solver where the edge is.
 *
 * F IS THEN A FUNCTION OF ( R, z, psi ) AND meq::Source's interface survives
 * untouched, which is FREE-BOUNDARY-PLAN.md section 5.3's first claim made
 * concrete.  Recovering phi from psi means inverting
 *
 *     g( t ) = t + c t^m ,      g' = 1 + c m t^( m - 1 ) > 0 ,
 *
 * which is monotone on t >= 0 and so has one root, bracketed by [ 0, psi ]
 * because g( 0 ) = 0 and g( psi ) >= psi.  A safeguarded Newton gets it to
 * round-off in a handful of steps.  There is no closed form for m > 2 and that
 * is not a difficulty -- the fixture is allowed to root-find, the SOLVER is
 * not.
 *
 * WHY THERE IS A BACKGROUND SOURCE OUTSIDE THE PLASMA, and it is not laziness.
 * A compactly contained plasma cannot sit in a source-free field: Delta* has no
 * zeroth-order term, so it obeys a maximum principle and { w > w0 } can never
 * be compactly contained for a Delta*-harmonic w.  A real equilibrium confines
 * the plasma with COILS -- sources outside it -- and a manufactured
 * fixed-boundary problem has to put something there instead.  Here it is
 * -Delta* phi, which is smooth and O( 1 ), and it is the reason PlasmaEdge
 * above (whose exterior IS source-free, at the price of a fixed cut) is kept
 * as the other half of the pair.
 *
 * AND THE JACOBIAN IS EXACT WITHOUT A CUT-RULE DERIVATIVE, which is worth
 * saying because section 5.3 names that derivative as FB-4's one real gap.
 * meq::SourceIntegrator evaluates dF/dpsi pointwise on a rule whose points do
 * NOT move, so the assembled Jacobian is the exact derivative of the assembled
 * residual whatever the edge is doing.  What the moving edge costs is
 * SMOOTHNESS of that residual in the unknowns, not consistency:
 *
 *     j = 0    F jumps in psi     the residual is Lipschitz, not C^1
 *     j >= 1   F is C^0 in psi    the residual is C^1 and Newton is ordinary
 *
 * so the gap only opens if a cut rule is adopted, and adopting one is what
 * would create the problem it is meant to solve.
 */
class MovingPlasmaEdge
{
	public:
		MovingPlasmaEdge( int vanishingOrderIn,
		                  double centreRIn = 1.0, double centreZIn = 0.0,
		                  double radiusIn = 0.23456789, double amplitudeIn = 4.0 )
			: jOrder( vanishingOrderIn ), mOrder( vanishingOrderIn + 2 ),
			  centreR( centreRIn ), centreZ( centreZIn ), discRadius( radiusIn ),
			  amplitude( amplitudeIn )
		{
			if ( vanishingOrderIn < 0 )
				throw std::invalid_argument(
					"meq::analytic::MovingPlasmaEdge: the vanishing order must "
					"be non-negative" );
			if ( amplitudeIn <= 0.0 )
				throw std::invalid_argument(
					"meq::analytic::MovingPlasmaEdge: the amplitude sets the "
					"sign of g' and must be positive for the inversion to be "
					"single valued" );
		}

		double levelSet( double radius, double z ) const
		{
			double const dr = radius - centreR, dz = z - centreZ;
			return discRadius*discRadius - dr*dr - dz*dz;
		}

		double psi( double radius, double z ) const
		{
			double const p = levelSet( radius, z );
			return p + ( p > 0.0 ? amplitude*std::pow( p, mOrder ) : 0.0 );
		}

		void flux( double radius, double z, double &qR, double &qZ ) const
		{
			double const p = levelSet( radius, z );
			double const dr = radius - centreR, dz = z - centreZ;
			double scale = 1.0;
			if ( p > 0.0 )
				scale += amplitude*mOrder*std::pow( p, mOrder - 1 );
			qR = scale*( -2.0*dr )/radius;
			qZ = scale*( -2.0*dz )/radius;
		}

		/// phi recovered from psi: the root of t + c t^m = psi on t >= 0.
		/// Bracketed by [ 0, psi ] and bisected where Newton would leave it.
		double supportVariable( double psiValue ) const
		{
			if ( psiValue <= 0.0 )
				return psiValue;

			double lo = 0.0, hi = psiValue, t = psiValue;
			for ( int i = 0; i < 80; ++i )
			{
				double const g = t + amplitude*std::pow( t, mOrder ) - psiValue;
				if ( g > 0.0 ) hi = t; else lo = t;
				double const gp = 1.0 + amplitude*mOrder*std::pow( t, mOrder - 1 );
				double next = t - g/gp;
				if ( !( next > lo && next < hi ) )
					next = 0.5*( lo + hi );
				if ( std::fabs( next - t ) <= 1.0e-16*std::fabs( next ) )
					return next;
				t = next;
			}
			return t;
		}

		/// F( R, z, psi ). The support is { psi > 0 } and nothing outside this
		/// class knows where the edge is.
		double f( double radius, double z, double psiValue ) const
		{
			double const dr = radius - centreR;
			double const deltaStarPhi = -4.0 + 2.0*dr/radius;
			double value = -deltaStarPhi;              // the background term
			if ( psiValue <= 0.0 )
				return value;

			double const t = supportVariable( psiValue );
			double const m = static_cast<double>( mOrder );
			double bracket = m*( m - 1.0 )*std::pow( t, mOrder - 2 )
			                 *4.0*( discRadius*discRadius - t );
			bracket += m*std::pow( t, mOrder - 1 )*deltaStarPhi;
			return value - amplitude*bracket;
		}

		double dFdPsi( double radius, double /*z*/, double psiValue ) const
		{
			if ( psiValue <= 0.0 )
				return 0.0;

			double const dr = radius - centreR;
			double const deltaStarPhi = -4.0 + 2.0*dr/radius;
			double const t = supportVariable( psiValue );
			double const m = static_cast<double>( mOrder );

			// d/dt of the bracket in f(). The ( m - 2 ) t^( m - 3 ) term has a
			// vanishing coefficient at m = 2 and a singular factor beside it,
			// so it is dropped rather than evaluated.
			double dBracket = -4.0*m*( m - 1.0 )*std::pow( t, mOrder - 2 );
			if ( mOrder > 2 )
				dBracket += 4.0*m*( m - 1.0 )*( m - 2.0 )
				            *std::pow( t, mOrder - 3 )*( discRadius*discRadius - t );
			dBracket += m*( m - 1.0 )*std::pow( t, mOrder - 2 )*deltaStarPhi;

			double const dtdPsi = 1.0/( 1.0 + amplitude*mOrder
			                            *std::pow( t, mOrder - 1 ) );
			return -amplitude*dBracket*dtdPsi;
		}

		/// Richardson extrapolated, for the reason PlasmaEdge::deltaStarFD gives.
		double deltaStarFD( double radius, double z, double h = 1.0e-3 ) const
		{
			auto raw = [ & ]( double step )
			{
				double const dRR = ( psi( radius + step, z ) - 2.0*psi( radius, z ) + psi( radius - step, z ) )/( step*step );
				double const dZZ = ( psi( radius, z + step ) - 2.0*psi( radius, z ) + psi( radius, z - step ) )/( step*step );
				double const dR  = ( psi( radius + step, z ) - psi( radius - step, z ) )/( 2.0*step );
				return dRR - dR/radius + dZZ;
			};
			return ( 4.0*raw( 0.5*h ) - raw( h ) )/3.0;
		}

		bool awayFromEdge( double radius, double z, double margin ) const
		{
			double const dr = radius - centreR, dz = z - centreZ;
			return std::fabs( std::sqrt( dr*dr + dz*dz ) - discRadius ) > margin;
		}

		int vanishingOrder() const { return jOrder; }
		int kinkOrder()      const { return mOrder; }
		double rateCap()     const { return static_cast<double>( mOrder ) + 0.5; }

	private:
		int jOrder, mOrder;
		double centreR, centreZ, discRadius, amplitude;
};

}
}

#endif
