#ifndef MEQ_PROFILES_HPP
#define MEQ_PROFILES_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

/*
 * Profile functions for MEQ.
 *
 * MEQ solves the fixed-boundary Grad-Shafranov equation
 *
 *     -div_bar( ( 1/R ) grad_bar( psi ) ) = F( R, z, psi ) / R    in Omega
 *                                    psi  = 0                    on Gamma
 *
 * with  F( R, z, psi ) := mu0 R^2 dp/dpsi + g dg/dpsi  (Sanchez-Vizuet & Solano,
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
	 * The Grad-Shafranov solve is nonlinear in psi and MEQ closes it with Newton,
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
			/// F = mu0 R^2 dp/dpsi with p = P0( psi ) exp( C( psi ) ( R^2 - R_ref^2 )/2 ),
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
	 * wrote it in, and that is frequently not MEQ's. A temperature table in keV
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
			/**
			 * This profile resampled onto @a intervals uniform intervals over
			 * its own range.
			 *
			 * **WHY YOU MIGHT WANT IT.** Interval lookup is `O( log n )` with
			 * unpredictable branches on arbitrary knots and one multiply on
			 * uniform ones -- which on a host is a mispredict per quadrature
			 * point and on a device is warp divergence in the innermost loop of
			 * the source kernel, where neighbouring points sit in different
			 * intervals and the whole warp walks every branch. Every profile
			 * MEQ ships is already uniform; this is for one that is not.
			 *
			 * **IT IS NOT THE SAME FUNCTION, AND THAT IS WHY IT IS A METHOD
			 * RATHER THAN SOMETHING THE LOADER DOES.** The result interpolates
			 * this profile's value AND derivative at each new knot, so it
			 * agrees exactly there and is a *different cubic* between them.
			 * Where a new interval sits inside an old one the two agree to
			 * round-off; where a new interval STRADDLES an old knot, one cubic
			 * is standing in for two and the error is the Hermite
			 * interpolation error of a function whose third derivative jumps
			 * inside the interval. `theResamplingErrorIsBoundedAndMeasured`
			 * in ProfilesTests puts a number on it rather than leaving this
			 * paragraph to be believed.
			 *
			 * **SO IT MUST BE DONE AT LOAD TIME AND TO THE PROFILE ITSELF**,
			 * never to a device view of one. MEQ asserts its POD source view
			 * against this class to the BIT; resampling one side and not the
			 * other would break that agreement, and the agreement is the only
			 * thing holding the two expressions of `F` together.
			 *
			 * @param intervals how many uniform intervals to use. More is
			 *        closer and none is exact.
			 * @throws std::invalid_argument if @a intervals is zero.
			 */
			SplineProfile resampledUniform( unsigned int intervals ) const;

			/// Whether the knots are equally spaced, to a relative 1e-9 -- the
			/// question resampledUniform() exists to turn from no into yes, and
			/// the one a caller should ask before spending anything on it.
			bool uniformKnots() const;

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


	/**
	 * THE NOT-A-KNOT CUBIC SPLINE THROUGH VALUES AT UNIFORM KNOTS, AS A LINEAR
	 * MAP OF THOSE VALUES: its CARDINAL BASIS, phi_k( x ) being the spline
	 * through one at knot k and zero at every other.
	 *
	 * WHY IT EXISTS BESIDE SplineProfile. A Hermite table takes its slopes as
	 * data, so a function of its knot VALUES alone has to make the slopes up,
	 * and the interpolant is then C^1 with a second derivative that jumps at
	 * every knot. This one is C^2 and depends on the values and nothing else:
	 * f( x ) = sum_k f_k phi_k( x ) exactly, so d f( x )/d f_k = phi_k( x ) at
	 * every derivative level. That is what the bordered midplane gauge needs
	 * of its crossings when R_0 at the knots are unknowns of the Newton
	 * (GradShafranovSolver::setMidplaneGaugeKnots()): the source's column for
	 * knot k is the chain rule through phi_k, and phi_k is this.
	 *
	 * NOT-A-KNOT and not natural: a natural spline puts f'' = 0 at both ends,
	 * which is an O( h^2 ) error there; not-a-knot keeps the O( h^4 ) of the
	 * interior to the ends. Each phi_k spans the whole range, decaying by about
	 * 2 - sqrt( 3 ) = 0.27 a knot away from k.
	 *
	 * Outside [ lower, upper ] the value is clamped to the nearer end and both
	 * derivatives are zero -- SplineProfile's convention.
	 */
	class CardinalSpline
	{
		public:
			/// @throws std::invalid_argument for fewer than four knots, which
			///         not-a-knot needs, or a range that is not finite and
			///         increasing.
			CardinalSpline( double lower, double upper, std::size_t knots );

			std::size_t size() const;
			double lower() const;
			double upper() const;
			/// The k-th knot's abscissa, lower + k ( upper - lower )/( size - 1 ).
			double knot( std::size_t k ) const;

			/**
			 * phi_k( x ), phi_k'( x ) and phi_k''( x ) for every k, into arrays
			 * of size() entries; any of the three may be null.
			 */
			void basis( double x, double *value, double *slope, double *curvature ) const;

			/// The spline's knot second derivatives, sum_j A_kj f_j: what
			/// evaluate() reads, so a caller holding values takes them once.
			std::vector<double> curvatures( std::vector<double> const &values ) const;

			/// The spline through @a values, whose curvatures() are
			/// @a curvatures, and its two derivatives at x.
			void evaluate( std::vector<double> const &values,
			               std::vector<double> const &curvatures, double x,
			               double &value, double &slope, double &curvature ) const;

		private:
			// The interval holding x, clamped into [ 0, size - 2 ].
			std::size_t interval( double x ) const;

			double xLower, xUpper, spacing;
			std::size_t count;
			// A_kj, row k: the knot curvature M_k as a combination of the
			// values. Dense and K^2, since every phi_j reaches every knot.
			std::vector<double> curvatureMap;
	};

	/**
	 * A meq::CardinalSpline through given values, as a meq::Profile.
	 *
	 * The curvatures are taken once at construction, so an evaluation is one
	 * interval's cubic. Shares its basis by pointer: the bordered midplane
	 * gauge builds one of these per Newton trial on the same knots.
	 */
	class CardinalSplineProfile : public Profile
	{
		public:
			/// @throws std::invalid_argument for a null basis or a value count
			///         that is not its size().
			CardinalSplineProfile( std::shared_ptr<CardinalSpline const> basis,
			                       std::vector<double> values );

			double operator()( double psi ) const override;
			double prime( double psi ) const override;
			/// Continuous: the spline is C^2, so unlike SplineProfile's this does
			/// not jump at a knot.
			double doublePrime( double psi ) const override;

			CardinalSpline const & basis() const;
			std::vector<double> const & values() const;

		private:
			std::shared_ptr<CardinalSpline const> splineBasis;
			std::vector<double> valueData;
			std::vector<double> curvatureData;
	};

	/**
	 * THE VARIABLE A TABLE IS TABULATED IN, against the normalised flux every
	 * profile here is a function of. PYTHON-PLAN.md, PY-2.
	 *
	 * `Psi` is MEQ's own, ( psi - psi_bnd )/( psi_ax - psi_bnd ): 0 at the
	 * boundary, 1 at the axis or a mirror's inner label
	 * (docs/normalised_flux.rst), and what a normalised source's profile file
	 * is in. `PsiN` is the conventional normalised poloidal flux, 1 - Psi, 0
	 * at the axis -- geq's psi_n -- and `SqrtPsiN` its square root, the usual
	 * rho. `PoloidalFlux` is psi itself, Wb/rad, the total flux the .nc's
	 * `psi` is: what an unnormalised source's profiles are functions of, and
	 * for a normalised one a map that moves with the labels
	 * (meq::FluxTableProfile). `RMidplane` is the radius, m, where the
	 * surface crosses z = 0: a mirror's, mapped through a crossing R_0( Psi )
	 * HELD for the solve (meq::CompositeProfile, makeNormalisedSource()).
	 *
	 * **A NAME, NOT A DEFAULT, WHEREVER A CALLER CHOOSES ONE**: tables read in
	 * psi_n where Psi was meant converge, at full order, to an answer 8.56e-04
	 * away (MEASUREMENTS.md M-263's mutation), so the swap is silent.
	 */
	enum class ProfileAbscissa
	{
		Psi,
		PsiN,
		SqrtPsiN,
		PoloidalFlux,
		RMidplane
	};

	/// One profile handed in as a table in place of a file: knots in
	/// `abscissa`, with the derivative in that variable too.
	struct ProfileTable
	{
		ProfileAbscissa abscissa = ProfileAbscissa::Psi;
		std::vector<Knot> knots;
	};

	/// Tables by the slot they replace, named by the key that names the file
	/// less `File`: "Omega", "GGPrime", "species.<Name>.Density",
	/// "species.<Name>.Temperature".
	using ProfileTables = std::map<std::string, ProfileTable>;

	/**
	 * A table as a profile of Psi -- ALWAYS A SplineProfile, the one a file of
	 * the mapped knots would make, so the device path reads every table
	 * (meq::SourceView) and the solver sees a C^1 cubic in Psi whatever the
	 * table was written in.
	 *
	 * `Psi` is the table's own SplineProfile, to the bit, and so is
	 * `PoloidalFlux`, a profile of psi for a source whose argument is psi; a
	 * normalised source wraps it in a FluxTableProfile. So is `RMidplane`, a
	 * profile of R, which a normalised source composes with its held
	 * crossing (CompositeProfile). `PsiN` is
	 * RE-KNOTTED EXACTLY: x = 1 - Psi is affine, and a Hermite cubic composed
	 * with an affine map is the Hermite cubic on the mapped knots, so the knots
	 * reverse, each abscissa becomes 1 - x and each slope changes sign.
	 *
	 * `SqrtPsiN` IS RE-KNOTTED AT ITS OWN KNOTS, Psi_i = 1 - x_i^2, with the
	 * value carried over and the slope f'( x )/( -2x ), and at the axis, x = 0,
	 * that slope's limit -f''( 0 )/2, f'' read off the first interval as the
	 * even function the profile is required to be.
	 * Not exact, as the affine map is: the cubic in Psi agrees with the table
	 * at every knot, value and slope, and between them to the tabulation's
	 * own order. **COMPOSITION, g( Psi ) = f( sqrt( 1 - Psi ) ), IS THE
	 * OTHER ROUTE AND IT IS NOT TAKEN FOR TWO REASONS**: the
	 * device path cannot read a composed profile, and its second derivative
	 * grows like ( 1 - Psi )^( -1/2 ) at the axis, which a cubic in x with
	 * zero slope there still carries through its x^3 term; the re-knotted
	 * cubic has neither. A table that reaches the axis is continued past it,
	 * linearly, to sqrtPsiNBeyondAxis; profileFromTable() says why.
	 *
	 * @throws std::invalid_argument for a table SplineProfile refuses -- fewer
	 *         than two knots, a non-finite entry, abscissae not increasing --
	 *         and for a `SqrtPsiN` table below x = 0, or that reaches x = 0
	 *         with a slope there: a profile smooth in the flux is even in rho,
	 *         and one that is not has a slope in Psi unbounded at the axis.
	 */
	std::shared_ptr<Profile const> profileFromTable( ProfileTable const & table );

	/**
	 * THE TWO LABELS A NORMALISED SOURCE IS AT, shared with the tables that read
	 * them. The source writes them in setNormalisation(), once per residual
	 * evaluation and before any thread reads them, and nothing else writes.
	 */
	struct FluxLabels
	{
		double axis = 1.0;
		double boundary = 0.0;
	};

	/**
	 * A TABLE IN psi, AS A PROFILE OF Psi FOR A NORMALISED SOURCE:
	 * g( Psi ) = span^k f( psi_bnd + span Psi ), span = psi_ax - psi_bnd, with
	 * g' = span^( k + 1 ) f' and g'' = span^( k + 2 ) f'', each the exact
	 * derivative of the one before, at whatever labels the source was last set
	 * to. k = 0 FOR A VALUE -- a density, a temperature, omega -- and k = 1
	 * FOR A DERIVATIVE IN THE FLUX: a normalised source's GGPrime is
	 * d( g^2/2 )/dPsi, which it divides by the span, so a table of
	 * d( g^2/2 )/dpsi reaches it multiplied by the span once.
	 *
	 * WHAT IT IS FOR. A profile in psi itself does not move with the labels at
	 * fixed psi, and a normalised source whose every profile is one, or a
	 * constant, has F label-free at fixed psi: its derivatives in the labels are
	 * zero, exactly, which meq::NormalisedRotatingSource::setProfilesInFlux()
	 * is how it is told. One such table beside a profile in Psi would make
	 * those derivatives neither the Psi form nor zero, and makeNormalisedSource()
	 * refuses the pair.
	 *
	 * Not a SplineProfile, so the device path does not read it and a source
	 * holding one evaluates on the host (meq::SourceView).
	 */
	class FluxTableProfile : public Profile
	{
		public:
			/// @param tableIn    the table, a profile of psi.
			/// @param labelsIn   the source's labels; borrowed by sharing.
			/// @param spanPower  k above: 0 for a value, 1 for a derivative.
			/// @throws std::invalid_argument for a null table or labels, or a
			///         power that is neither.
			FluxTableProfile( std::shared_ptr<Profile const> tableIn,
			                  std::shared_ptr<FluxLabels const> labelsIn,
			                  int spanPower = 0 );

			double operator()( double psi ) const override;
			double prime( double psi ) const override;
			double doublePrime( double psi ) const override;

		private:
			std::shared_ptr<Profile const> inPsi;
			std::shared_ptr<FluxLabels const> labels;
			int power;
	};

	/**
	 * A PROFILE OF A PROFILE: g( x ) = outer( inner( x ) ), with
	 * g' = outer'( inner ) inner' and g'' = outer''( inner ) inner'^2
	 * + outer'( inner ) inner'', each the exact derivative of the one before.
	 * What an `RMidplane` table is to a normalised source: outer the table in
	 * R, inner the held crossing R_0( Psi ). Outside inner's range inner
	 * clamps, so g holds outer's value there with no slope.
	 *
	 * Not a SplineProfile, so the device path does not read it (meq::
	 * SourceView) and a source holding one evaluates on the host.
	 */
	class CompositeProfile : public Profile
	{
		public:
			/// @throws std::invalid_argument for either null.
			CompositeProfile( std::shared_ptr<Profile const> outerIn,
			                  std::shared_ptr<Profile const> innerIn );

			double operator()( double psi ) const override;
			double prime( double psi ) const override;
			double doublePrime( double psi ) const override;

		private:
			std::shared_ptr<Profile const> outer;
			std::shared_ptr<Profile const> inner;
	};

	/// How far past the axis a `SqrtPsiN` table that reaches it is continued,
	/// linearly, by profileFromTable(): rho cannot name Psi > 1, and an
	/// overshooting iterate needs a profile there that is C^1 through the axis.
	/// rotating-normalised's own table reaches as far.
	inline constexpr double sqrtPsiNBeyondAxis = 1.5;

	/// "Psi", "psi_n", "sqrt_psi_n", "psi", "R_midplane": the names a front
	/// end takes. @throws std::invalid_argument naming the five for any other.
	ProfileAbscissa profileAbscissaNamed( std::string const & name );
	char const * profileAbscissaName( ProfileAbscissa abscissa );

}

#endif // MEQ_PROFILES_HPP
