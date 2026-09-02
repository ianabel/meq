#ifndef MEQ_PROFILES_HPP
#define MEQ_PROFILES_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

/*
 * Profile functions for meq.
 *
 * meq solves the fixed-boundary Grad-Shafranov equation
 *
 *     -div_bar( ( 1/r ) grad_bar( psi ) ) = F( r, z, psi ) / r    in Omega
 *                                    psi  = 0                    on Gamma
 *
 * with  F( r, z, psi ) := mu0 r^2 dp/dpsi + g dg/dpsi  (Sanchez-Vizuet & Solano,
 * CPC 235 (2019) 120-132, eqs (1)-(4)). Everything in F that is not geometry is
 * a user-supplied function of the flux alone: the pressure p( psi ) and the
 * toroidal field function g( psi ). Those are the Profiles defined here.
 *
 * Nothing in this file knows about MFEM: profiles are plain functions of a
 * double, so they can be unit tested on their own. The mfem::Coefficient
 * adapter lives elsewhere.
 */

namespace meq
{

	/**
	 * A scalar function of the poloidal flux psi, together with its derivative.
	 *
	 * The Grad-Shafranov solve is nonlinear in psi and meq closes it with Newton,
	 * so every profile has to be able to report its own derivative: prime() is not
	 * an optional extra, it is half the interface.
	 *
	 * Normalisation convention: profiles are functions of the *normalised* flux,
	 * which by convention runs over [ 0, 1 ]. This class does not enforce that --
	 * SplineProfile is defined on whatever knot range it was handed -- but the
	 * whole library assumes it, and a table given on [ 0, 1 ] is what the
	 * configuration file is expected to supply. Mapping the solver's psi onto the
	 * profile's normalised flux is the caller's job (see meq::Source), because the
	 * map changes from one Newton iterate to the next and a profile is a fixed
	 * function.
	 *
	 * Out-of-range convention: implementations must not throw when handed a psi
	 * outside their natural range. A Newton iterate routinely overshoots, and an
	 * exception thrown out of a quadrature loop would kill a solve that would
	 * otherwise have converged. Implementations clamp instead; see SplineProfile
	 * for the exact policy.
	 */
	class Profile
	{
		public:
			virtual ~Profile() = default;

			/// The value of the profile at flux psi. Units are the profile's own
			/// (e.g. Pa/(Wb/rad) for a dp/dpsi profile); psi is dimensionless if
			/// the caller normalises it as described above.
			virtual double operator()( double psi ) const = 0;

			/// d/dpsi of the profile at psi, in the profile's units per unit psi.
			/// Required by the Newton solve; must be the exact derivative of
			/// operator(), not a finite difference of it.
			virtual double prime( double psi ) const = 0;

			/// d2/dpsi2 of the profile at psi. Must be the exact derivative of
			/// prime(), for the same reason prime() must be the exact derivative
			/// of operator().
			///
			/// WHY THIS EXISTS, because it is not needed by every source and its
			/// absence was not felt for a long time. meq::MHDSource stores the
			/// *products* mu0-free p' and g g', so F is one profile evaluation and
			/// dF/dpsi is one prime() -- two levels, and Profile supplied both. A
			/// source whose F is itself a psi-derivative of something built from
			/// profiles needs three: meq::RotatingSource has
			/// F = mu0 r^2 dp/dpsi with p = P0( psi ) exp( C( psi ) ( r^2 - rRef^2 )/2 ),
			/// so F already spends one derivative of P0 and C and the Jacobian
			/// spends a second. There is no reparametrisation that avoids it --
			/// p is not a flux function, so there is no product to pre-store.
			///
			/// CAVEAT FOR TABULATED PROFILES. A piecewise Hermite cubic is C^1 and
			/// no more, so its second derivative is piecewise linear and jumps at
			/// every knot. That is a set of measure zero and does not move a
			/// converged answer, but a Newton step taken exactly at a knot sees a
			/// Jacobian that is one-sided. Nothing has measured what that costs.
			virtual double doublePrime( double psi ) const = 0;

		protected:
			// Copying through a Profile& would slice; derived classes are free to be
			// copyable, and are.
			Profile() = default;
			Profile( Profile const & ) = default;
			Profile( Profile && ) = default;
			Profile & operator=( Profile const & ) = default;
			Profile & operator=( Profile && ) = default;
	};

	/**
	 * A profile that does not depend on psi at all.
	 *
	 * This is the Solov'ev case of a general profile: p' and g g' constant. Its
	 * prime() is identically zero, so a Source built from ConstantProfiles
	 * contributes nothing to the Newton Jacobian -- which is correct, the equation
	 * is then linear in psi.
	 */
	class ConstantProfile : public Profile
	{
		public:
			/// A profile equal to `value` everywhere.
			explicit ConstantProfile( double value );

			/// Returns the constant, for any psi.
			double operator()( double psi ) const override;

			/// Returns 0 exactly, for any psi.
			double prime( double psi ) const override;

			/// Returns 0 exactly, for any psi.
			double doublePrime( double psi ) const override;

			/// The constant this profile was built with.
			double value() const;

		private:
			double constantValue;
	};

	/// One tabulated point of a profile: the flux, the value there, and the
	/// derivative there. A Hermite spline interpolates value *and* slope, so the
	/// slope is data, not something inferred from the neighbours.
	struct Knot
	{
		double psi;         ///< abscissa, normalised flux
		double value;       ///< f( psi )
		double derivative;  ///< f'( psi )
	};

	/**
	 * The Hermite cubic on a single interval [ lower, upper ].
	 *
	 * The unique cubic matching a value and a derivative at each end. It therefore
	 * reproduces any cubic exactly, and a chain of them is C^1 across the shared
	 * knots by construction.
	 *
	 * Out-of-range: an argument below lower or above upper is clamped -- the value
	 * is that of the nearer endpoint, the derivative is zero, which is the exact
	 * derivative of that constant extension. (The old code returned a function
	 * *value* from prime() here, which is a units error as well as a wrong answer.)
	 */
	class HermiteCubicSpline
	{
		public:
			/// Build from the two endpoints and the value/derivative data at each.
			/// Throws std::invalid_argument unless lower < upper.
			HermiteCubicSpline( double lower, double upper, double fLower, double fUpper, double fPrimeLower, double fPrimeUpper );

			/// Build from two knots; the first must lie strictly to the left of the
			/// second. Throws std::invalid_argument otherwise.
			HermiteCubicSpline( Knot const & lower, Knot const & upper );

			/// The interpolated value at x, clamped to the endpoint values outside
			/// [ lower, upper ].
			double operator()( double x ) const;

			/// The derivative of the interpolant at x; zero outside
			/// [ lower, upper ], consistent with the constant extension used by
			/// operator().
			double prime( double x ) const;

			/// The second derivative of the interpolant at x; zero outside
			/// [ lower, upper ], for the same reason prime() is. Piecewise linear
			/// in x, so adjacent intervals disagree at a shared knot -- see
			/// Profile::doublePrime.
			double doublePrime( double x ) const;

			/// ( lower, upper ) -- the interval this cubic is defined on.
			std::pair<double,double> interval() const;

			/// ( f( lower ), f( upper ) ).
			std::pair<double,double> values() const;

			/// ( f'( lower ), f'( upper ) ).
			std::pair<double,double> derivatives() const;

		private:
			double xLower, xUpper;
			double delta;
			double fLower, fUpper;
			double fPrimeLower, fPrimeUpper;
	};

	/**
	 * Another profile, multiplied by a constant.
	 *
	 * WHY THIS EXISTS: a tabulated profile arrives in whatever units its author
	 * wrote it in, and that is frequently not meq's. A temperature table in keV
	 * has to become Joules, a density in 10^19 m^-3 has to become m^-3, and a
	 * profile from another code may be normalised to its own reference values.
	 * Editing the file is the wrong answer -- it makes the file a function of
	 * which code reads it -- and folding the factor into the *source* is worse,
	 * because then every source has to carry a units argument per profile.
	 *
	 * The scale multiplies all three derivative levels, which is what makes this
	 * safe: a scaled profile is exactly the profile of the scaled quantity, and
	 * prime() and doublePrime() stay the exact derivatives of operator().
	 */
	class ScaledProfile : public Profile
	{
		public:
			/// @param inner  the profile to scale. Must not be null.
			/// @param scale  the factor. May be negative; a sign convention is a
			///               scale like any other.
			/// @throws std::invalid_argument if @a inner is null or @a scale is
			///         not finite.
			ScaledProfile( std::shared_ptr<Profile const> inner, double scale );

			double operator()( double psi ) const override;
			double prime( double psi ) const override;
			double doublePrime( double psi ) const override;

			/// The profile underneath, unscaled.
			Profile const & unscaled() const;
			double scale() const;

		private:
			std::shared_ptr<Profile const> innerProfile;
			double scaleFactor;
	};

	/**
	 * A profile given as a table of ( psi, f, f' ) and interpolated by a
	 * piecewise Hermite cubic.
	 *
	 * Because both the value and the slope are tabulated, the interpolant is C^1
	 * everywhere and matches the data exactly at every knot, and prime() is the
	 * analytic derivative of operator() rather than a difference of it. That is
	 * what makes the Newton Jacobian consistent with the residual.
	 *
	 * Out-of-range policy (applied identically by operator() and prime()): the
	 * table is extended by a *constant* beyond its end knots. Below the first knot
	 * the value is f( first ) and the derivative is 0; above the last knot the
	 * value is f( last ) and the derivative is 0. Nothing throws. The rationale is
	 * in Profile's documentation: Newton overshoots, and a bounded extension keeps
	 * the residual finite where a linear extrapolation of a steep edge profile
	 * would not. The extension is continuous but not C^1 -- there is a kink in the
	 * derivative at each end knot.
	 *
	 * The class is a plain value type: copyable and movable, with no hand-written
	 * special members to get wrong.
	 */
	class SplineProfile : public Profile
	{
		public:
			using RealFunction = std::function<double( double )>;

			/// Build from a table, which must have at least two knots with strictly
			/// increasing psi. Throws std::invalid_argument otherwise.
			explicit SplineProfile( std::vector<Knot> data );

			/// Sample f and its derivative fPrime on `intervals` equal intervals of
			/// [ 0, 1 ] and spline the result. Handy for analytic profiles and for
			/// tests. Throws std::invalid_argument if intervals == 0 or either
			/// function is empty.
			SplineProfile( RealFunction f, RealFunction fPrime, unsigned int intervals );

			/**
			 * Read a table from a stream.
			 *
			 * The format is one knot per line, `psi f(psi) f'(psi)`, whitespace
			 * separated; lines whose first non-blank character is '#' are comments;
			 * a blank line terminates the table, so several profiles can be
			 * concatenated in one stream and read one after another. This is the
			 * format write() emits.
			 *
			 * Throws std::runtime_error on a malformed line or on fewer than two
			 * knots, and std::invalid_argument if the knots are not increasing.
			 */
			static SplineProfile fromStream( std::istream & is );

			/// As fromStream, reading the named file. Throws std::runtime_error if
			/// the file cannot be opened or does not parse -- never returns an
			/// empty profile. (The version of this that this class replaces had a
			/// body of `return;`, which silently produced an empty spline table
			/// that then indexed out of bounds on first use.)
			static SplineProfile fromFile( std::string const & fileName );

			/// The interpolated value at psi; clamped outside the knot range, see
			/// the class documentation.
			double operator()( double psi ) const override;

			/// The derivative of the interpolant at psi; zero outside the knot
			/// range, see the class documentation.
			double prime( double psi ) const override;

			/// The second derivative of the interpolant at psi; zero outside the
			/// knot range. Piecewise linear, so it jumps at every interior knot --
			/// findInterval() decides which side a knot belongs to and the answer
			/// there is that interval's, not an average. See Profile::doublePrime.
			double doublePrime( double psi ) const override;

			/// Write the table in the format fromStream() reads, at enough
			/// precision to round-trip a double exactly, terminated by a blank
			/// line.
			void write( std::ostream & os ) const;

			/// The tabulated knots, in increasing psi.
			std::vector<Knot> const & knots() const;

			/// ( first knot psi, last knot psi ) -- the range over which this
			/// profile interpolates rather than clamps.
			std::pair<double,double> domain() const;

			/// The number of Hermite cubics, i.e. knots().size() - 1.
			std::size_t numIntervals() const;

			/// The Hermite cubic on [ knots()[ i ].psi, knots()[ i + 1 ].psi ].
			/// Throws std::out_of_range if i >= numIntervals().
			HermiteCubicSpline intervalAt( std::size_t i ) const;

		private:
			// Index of the interval containing psi, clamped into
			// [ 0, numIntervals() - 1 ] for a psi outside the table.
			std::size_t findInterval( double psi ) const;

			std::vector<Knot> knotData;
	};

	/// Write a profile in the format SplineProfile::fromStream() reads.
	std::ostream & operator<<( std::ostream & os, SplineProfile const & profile );

}

#endif // MEQ_PROFILES_HPP
