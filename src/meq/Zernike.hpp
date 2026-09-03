#ifndef MEQ_ZERNIKE_HPP
#define MEQ_ZERNIKE_HPP

#include <cstddef>
#include <vector>

/*
 * The Zernike polynomials on the unit disc.
 *
 *     Z_l^m( rho, theta ) = R_l^m( rho ) * { cos( m theta )    m >= 0
 *                                          { sin( |m| theta )  m <  0
 *
 * with R_l^m a polynomial of degree l in rho containing ONLY the powers
 * rho^|m|, rho^(|m|+2), ..., rho^l, admissible exactly when
 *
 *     l >= |m|    and    l - |m| is EVEN.
 *
 * THE INDEX CONSTRAINT IS THE ENTIRE POINT OF CHOOSING THIS BASIS, so it is
 * stated before anything else in the file. A naive tensor product -- "any
 * polynomial in rho" times "any Fourier mode in theta" -- admits terms such as
 * rho^2 cos( theta ) or rho cos( 2 theta ). Those are perfectly good functions
 * of ( rho, theta ) away from the origin and they are NOT smooth at rho = 0:
 * written in Cartesian coordinates, rho^2 cos( theta ) is x*sqrt( x^2 + y^2 ),
 * whose second derivative has no limit at the origin -- it reads
 * 3 cos( theta ) - cos^3( theta ), a different number down every ray. Such a
 * basis makes the centre of the disc a coordinate singularity that every
 * consumer has to special-case, and the special case is where the bugs live.
 *
 * The parity-and-minimum-power constraint excludes exactly those terms. What
 * survives is the statement that makes the basis worth having:
 *
 *     EVERY admissible Z_l^m IS A BIVARIATE POLYNOMIAL IN ( x, y ) OF DEGREE l.
 *
 * So a Zernike expansion is a polynomial in Cartesian coordinates, the centre
 * of the disc is an ordinary interior point of it, and derivatives there are
 * whatever the polynomial says they are -- no limit to take, no ray to choose,
 * no branch. tests/unit/ZernikeTests.cpp measures that as a contrast against a
 * deliberately parity-violating mode rather than asserting it as folklore.
 *
 * FOR MEQ THE CENTRE OF THE DISC IS THE MAGNETIC AXIS. INVERSION-PLAN.md
 * section 4.1 wants the flux surfaces as R( Psi, l ), z( Psi, l ) -- a map from
 * a disc, whose centre is the axis and whose boundary is the outermost surface.
 * The axis is precisely where the surfaces shrink to a point, where 1/|grad psi|
 * diverges, and where the in-surface coordinate l is undefined; it is the one
 * place a representation is most tempted to special-case and least able to
 * afford it. That is the whole argument for Zernike here, and it is a better
 * argument than "DESC does it".
 *
 * ---------------------------------------------------------------------------
 * THE RADIAL COORDINATE IS rho = sqrt( Psi_N ), NOT Psi_N. GETTING THIS WRONG
 * IS SILENT.
 * ---------------------------------------------------------------------------
 *
 * psi has a QUADRATIC maximum at the magnetic axis -- it is a smooth function
 * with a non-degenerate critical point there -- so psi - psi_ax falls off like
 * (geometric distance)^2 and the normalised flux Psi_N behaves like
 * (distance)^2 near the axis. The geometry R( . ), z( . ) is smooth in the
 * DISTANCE, so it is smooth in sqrt( Psi_N ) and carries a square-root branch
 * point in Psi_N itself.
 *
 * Parametrise by Psi_N directly and every basis -- Zernike, Chebyshev, anything
 * -- converges ALGEBRAICALLY against that branch point, worst near the axis,
 * and there is nothing in a convergence table to say why. The rate simply comes
 * out lower than the design order and looks like a discretisation problem. This
 * is the same species of defect as the ones CLAUDE.md records under "a wrong
 * Jacobian is invisible to a convergence table": the answer converges, to the
 * wrong rate, for a reason no assertion in the suite is looking at.
 *
 * DESC parametrises by the square root of the normalised toroidal flux for
 * exactly this reason, and its documentation notes that the quantity is
 * "proportional to the minor radius". So: the Zernike degree l counts powers of
 * sqrt( Psi_N ), and radiusFromNormalisedFlux() / normalisedFluxFromRadius() /
 * fluxDerivativeFromRadial() below are the conversions. Use them rather than
 * writing sqrt and a chain rule by hand at each call site, because the chain
 * rule factor is 1/( 2 rho ) and it is the thing that gets dropped.
 *
 * ---------------------------------------------------------------------------
 *
 * NOTHING IN THIS FILE INCLUDES MFEM, and that is deliberate in the same way
 * and for the same two reasons as src/meq/Profiles.hpp and src/meq/Source.hpp:
 * the Zernike polynomials are pure mathematics on plain doubles, so they are
 * unit-testable without the finite element library, and continuous integration
 * -- which cannot obtain the MFEM branch meq needs, see INSTALL.md -- can build
 * and test them. Anything that needs a mesh, a GridFunction or a Coefficient
 * belongs in the consumer, not here.
 *
 * WHAT IS DELIBERATELY ABSENT: there is no fitting routine and no quadrature.
 * A projection integral would need a rule on the disc and a callable f, and it
 * is not what INVERSION-PLAN.md IN-3 actually does -- IN-3 fits traced surface
 * POINTS, which is a least-squares problem against a point cloud, not an
 * integral against a function. Supplying a quadrature projection here would put
 * the wrong tool in the library and invite it to be used. The unit tests build
 * their own projection because a projection is the sharpest way to measure
 * orthogonality and coefficient decay, which is a different job.
 *
 * HOW IT IS EVALUATED: the radial polynomial IS a Jacobi polynomial under a
 * change of variable, so Zernike.cpp calls boost::math::jacobi() and
 * boost::math::jacobi_prime() rather than carrying a recurrence of its own.
 * That dependency is header-only and is confined to the .cpp -- this header
 * includes nothing but <cstddef> and <vector>, so a consumer of the basis takes
 * on nothing. The identity, the two sign conventions in it that are easy to get
 * wrong, and why the explicit factorial sum in every reference must not be used,
 * are all set out at the top of Zernike.cpp beside the code.
 */

namespace meq
{

	/// The largest degree this file will build a basis for. A cap rather than an
	/// open-ended size, in the spirit of meq::maxSpecies: the mode count grows
	/// like l^2/2, so a typo of 100000 where 30 was meant would ask for a vector
	/// of five thousand million entries and be diagnosed as a memory problem
	/// rather than as a typo. 128 is 8385 modes, far past anything a flux-surface
	/// fit has any business wanting, and Boost's Jacobi evaluation is still
	/// accurate there.
	inline constexpr int maxZernikeDegree = 128;

	/**
	 * One mode of the Zernike basis: the pair of integers that indexes it.
	 *
	 * A plain aggregate, in the style of meq::Knot and meq::Species -- a mode is
	 * data, and a caller enumerating modes wants to read the indices off, not to
	 * ask an object for them.
	 *
	 * The sign of m selects the angular function, which is the usual convention
	 * and is worth spelling out because it is a sign and this project has been
	 * bitten by signs: m >= 0 means cos( m theta ), m < 0 means sin( |m| theta ).
	 * m = 0 is therefore the cos branch and there is no sin( 0 ) mode, which is
	 * what makes the mode count ( l + 1 )( l + 2 )/2 rather than double it.
	 */
	struct ZernikeMode
	{
		int l;   ///< radial degree, l >= 0
		int m;   ///< azimuthal order; |m| <= l and l - |m| even
	};

	/// Modes compare equal when both indices do. Provided because tests and
	/// callers assembling mode lists want it; there is no ordering operator,
	/// because the ordering that matters is the enumeration order below and it
	/// is not a comparison of one mode against another.
	bool operator==( ZernikeMode left, ZernikeMode right );
	bool operator!=( ZernikeMode left, ZernikeMode right );

	/// True iff ( l, m ) is an admissible Zernike index pair: l >= 0,
	/// |m| <= l, and l - |m| even. This is the predicate every entry point in
	/// this file enforces; it is exposed so a caller can ask rather than catch.
	bool isValidZernikeMode( int l, int m );

	/// The number of admissible modes of degree at most maxDegree, which is
	/// ( maxDegree + 1 )( maxDegree + 2 )/2.
	/// @throws std::invalid_argument if maxDegree < 0 or > maxZernikeDegree.
	std::size_t zernikeModeCount( int maxDegree );

	/**
	 * Every admissible mode of degree at most maxDegree, in the standard
	 * ANSI/OSA order: by degree l ascending, and within a degree by m ascending
	 * from -l to +l in steps of two.
	 *
	 * The order is fixed and, more usefully, it is PREFIX STABLE: the first
	 * zernikeModeCount( d ) entries of zernikeModes( maxDegree ) are exactly
	 * zernikeModes( d ) for any d <= maxDegree. So truncating an expansion to a
	 * lower degree is truncating its coefficient vector, with no re-indexing,
	 * which is what a convergence study in mode number needs.
	 *
	 * @throws std::invalid_argument if maxDegree < 0 or > maxZernikeDegree.
	 */
	std::vector<ZernikeMode> zernikeModes( int maxDegree );

	/// Where ( l, m ) sits in zernikeModes(). Closed form -- ( l( l + 2 ) + m )/2
	/// -- and independent of maxDegree, by the prefix-stability above.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	std::size_t zernikeModeIndex( int l, int m );

	/// The radial polynomial R_l^m( rho ). Normalised so that R_l^m( 1 ) = 1,
	/// which is the universal convention and is what makes the coefficients of
	/// an expansion comparable across degrees. Depends on m only through |m|.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeRadial( int l, int m, double rho );

	/// d R_l^m / d rho, exactly -- the chain rule applied to Boost.Math's own
	/// analytic Jacobi derivative, with no differencing anywhere. A Newton
	/// method downstream will differentiate this again, and CLAUDE.md's standing
	/// finding is that a derivative that disagrees with its own function is
	/// invisible to every convergence table in the suite.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeRadialPrime( int l, int m, double rho );

	/// Z_l^m( rho, theta ).
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernike( int l, int m, double rho, double theta );

	/// d Z_l^m / d rho.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeRadialDerivative( int l, int m, double rho, double theta );

	/// d Z_l^m / d theta. Exactly zero for m = 0, not merely small.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeAngularDerivative( int l, int m, double rho, double theta );

	/// The radial-only square norm, integral over [ 0, 1 ] of
	/// R_l^m( rho )^2 rho d rho, which is 1/( 2( l + 1 ) ) -- independent of m.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeRadialNormSquared( int l, int m );

	/// The square norm on the disc, the integral of Z_l^m^2 rho d rho d theta
	/// over rho in [ 0, 1 ] and theta in [ 0, 2 pi ). It is pi/( l + 1 ) for
	/// m = 0 and pi/( 2( l + 1 ) ) otherwise -- THE FACTOR OF TWO ON m = 0 IS
	/// REAL and comes from cos( 0 ) integrating to 2 pi where cos( m theta )^2
	/// integrates to pi. Dropping it scales every m = 0 coefficient of a
	/// projection by a half, which looks like a fit that is merely poor.
	/// @throws std::invalid_argument if ( l, m ) is not admissible.
	double zernikeNormSquared( int l, int m );

	// -----------------------------------------------------------------------
	// The radial coordinate, and the chain rule that goes with it.
	// -----------------------------------------------------------------------

	/// rho = sqrt( Psi_N ). See the file header for why the square root is not
	/// optional.
	/// @throws std::invalid_argument if normalisedFlux is negative or not
	///         finite. Negative is rejected rather than clamped because it means
	///         the caller's normalisation is wrong -- Psi_N is a fraction of the
	///         flux between the axis and the boundary and cannot be negative --
	///         and a silently clamped zero would put the surface on the axis.
	double radiusFromNormalisedFlux( double normalisedFlux );

	/// Psi_N = rho^2. The inverse of the above, for the same reason.
	/// @throws std::invalid_argument if rho is negative or not finite.
	double normalisedFluxFromRadius( double rho );

	/**
	 * Convert a d/drho into a d/dPsi_N.
	 *
	 * THE FACTOR IS 1/( 2 rho ). It is neither 1 nor 2 rho, and both of those
	 * mistakes produce a quantity with the right units, the right sign and the
	 * right qualitative shape:
	 *
	 *     rho = sqrt( Psi_N )   =>   d Psi_N = 2 rho d rho
	 *                           =>   df/dPsi_N = ( df/drho ) / ( 2 rho ).
	 *
	 * This exists as a named function rather than as three characters at each
	 * call site because INVERSION-PLAN.md IN-3 exists to deliver dGeometry/dPsi
	 * for MANTA-COUPLING.md, that derivative is the deliverable, and a factor
	 * dropped in it would show up as a coupled run that converges to the wrong
	 * equilibrium rather than as anything failing.
	 *
	 * IT IS UNBOUNDED AT rho = 0 AND THAT IS THE COORDINATE, NOT A DEFECT. Even
	 * for a representation that is perfectly smooth in Cartesian coordinates,
	 * d/dPsi_N at fixed theta diverges at the axis whenever the m = 1 content is
	 * non-zero -- the m = 1 modes are the rigid shift of a surface, R ~ rho cos
	 * theta, and moving the axis sideways at a rate 1/( 2 rho ) is exactly what
	 * a Shafranov shift does as the surfaces shrink onto a point. A caller
	 * wanting a quantity that stays finite there wants d/drho, which this file
	 * also supplies. No exception is thrown: rho = 0 returns an infinity or a
	 * NaN by the ordinary floating point rules, because this is called per
	 * surface point and a throw in that loop would be worse than the infinity.
	 */
	double fluxDerivativeFromRadial( double radialDerivative, double rho );

	/**
	 * A truncated Zernike expansion: a maximum degree and one coefficient per
	 * mode, in zernikeModes() order.
	 *
	 * A plain value type -- copyable, movable, no hand-written special members.
	 * It holds no scratch and no cached state, so a single expansion may be
	 * evaluated from several threads at once, which a flux-surface consumer
	 * sampling a grid will want.
	 *
	 * Evaluation walks each |m| once rather than each mode independently, so the
	 * two trigonometric functions are evaluated once per azimuthal order rather
	 * than once per mode, and the cosine and sine coefficients of a degree are
	 * taken together. It allocates nothing and holds no cached state.
	 */
	class ZernikeExpansion
	{
		public:
			/// An expansion of the given degree with every coefficient zero.
			/// @throws std::invalid_argument if maxDegree is out of range.
			explicit ZernikeExpansion( int maxDegree );

			/// An expansion of the given degree with the supplied coefficients,
			/// which must be exactly zernikeModeCount( maxDegree ) of them, in
			/// zernikeModes() order.
			/// @throws std::invalid_argument if the size does not match or the
			///         degree is out of range. The size check is the whole reason
			///         to pass the degree separately rather than inferring it:
			///         the mode counts are triangular numbers, so a coefficient
			///         vector of the wrong length is usually still a valid length
			///         for some other degree and would be accepted silently.
			ZernikeExpansion( int maxDegree, std::vector<double> coefficients );

			/// The value of the expansion at ( rho, theta ).
			double operator()( double rho, double theta ) const;

			/// Value and both first derivatives in one pass. Cheaper than three
			/// calls and the form the tracer and the geometry consumer want.
			void evaluate( double rho, double theta,
			               double & value, double & dRho, double & dTheta ) const;

			/// d/drho of the expansion.
			double radialDerivative( double rho, double theta ) const;

			/// d/dtheta of the expansion.
			double angularDerivative( double rho, double theta ) const;

			/// d/dPsi_N of the expansion, at normalised flux Psi_N and angle
			/// theta. Exactly fluxDerivativeFromRadial() applied to
			/// radialDerivative() at rho = sqrt( Psi_N ); see that function for
			/// the factor and for its behaviour on the axis.
			/// @throws std::invalid_argument if normalisedFlux is negative.
			double fluxDerivative( double normalisedFlux, double theta ) const;

			/// The degree this expansion was built at.
			int maxDegree() const;

			/// The modes, in the order the coefficients are stored.
			std::vector<ZernikeMode> const & modes() const;

			/// The coefficients, in modes() order.
			std::vector<double> const & coefficients() const;

			/// The coefficient of Z_l^m.
			/// @throws std::invalid_argument if ( l, m ) is not admissible or its
			///         degree exceeds maxDegree().
			double coefficient( int l, int m ) const;

			/// Set the coefficient of Z_l^m.
			/// @throws std::invalid_argument as coefficient() does.
			void setCoefficient( int l, int m, double value );

			/// A copy truncated to a lower degree. Free, by the prefix stability
			/// of the mode order -- it is a prefix of the coefficient vector --
			/// and it is what a convergence study in mode number is made of.
			/// @throws std::invalid_argument if newDegree is negative or exceeds
			///         maxDegree().
			ZernikeExpansion truncated( int newDegree ) const;

		private:
			int degree;
			std::vector<ZernikeMode> modeList;
			std::vector<double> coefficientList;
	};

}

#endif // MEQ_ZERNIKE_HPP
