#ifndef MEQ_ROTATINGSOURCE_HPP
#define MEQ_ROTATINGSOURCE_HPP

#include "Profiles.hpp"
#include "Source.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

/*
 * The Grad-Shafranov source for a plasma in sonic toroidal rotation.
 *
 * The equation is refs/RotatingGK.pdf eq (136) -- Abel, Plunk, Wang, Barnes,
 * Cowley, Dorland & Schekochihin, Rep. Prog. Phys. 76 (2013) 116201 -- closed by
 * its (96) for the poloidal density variation and (97) for the electrostatic
 * potential phi_0 that holds quasineutrality against it. docs/rotation.rst is
 * the derivation -- the equation written out, the collapse below, and the unit
 * and sign conventions against the source paper.
 *
 * WHAT ROTATION CHANGES, AND IT IS ONE THING: the density is no longer a flux
 * function. Centrifugal force sweeps heavy species to the outboard side, and an
 * electrostatic potential arises to stop that separating the charges. So
 *
 *     n_s( R, psi ) = n_s0( psi ) exp[ m_s omega^2 ( R^2 - R_ref^2 )/2T_s
 *                                      - Z_s e phi_0/T_s ]          (96)
 *     Sum_s Z_s n_s( R, psi ) = 0                                   (97)
 *
 * and (97) is what determines phi_0 at each ( R, psi ). Everything else about
 * MEQ is untouched: the operator, the discretisation, tau, the hybridization,
 * the estimator, the adaptive loop and the curved boundary all stay as they are.
 * meq::Source's signature already carries R, which is the whole reason a
 * non-flux-function density needs no interface change anywhere in the solver.
 *
 * THE GAUGE IS FREE AND THIS CLASS TAKES THE LOCAL ONE. RoPP fixes phi_0 by
 * <phi_0>_psi = 0 and says at its (59) that this is a convention rather than a
 * closure -- "we can add any function of psi to it". That choice would make
 * phi_0 proportional to R^2 - <R^2>_psi, a flux-surface average of the unknown
 * on every surface, which MEQ has no machinery for and which buys nothing
 * physical. This class instead pins
 *
 *     phi_0( R_ref, psi ) = 0
 *
 * with R_ref a CONSTANT -- the geometric axis, given to the constructor, not the
 * magnetic axis and not an average. Two consequences, both good: F is a
 * pointwise function of ( R, z, psi ), and n_s0 is then the PHYSICAL density of
 * species s on the curve R = R_ref, which is a quantity a user can state and
 * another code can be compared against. Li & Zhu (Comput. Phys. Commun. 260
 * (2021) 107264) make the same choice independently, referencing their exponent
 * to the magnetic axis.
 *
 * A MIRROR TAKES THE OTHER LOCAL GAUGE, Gauge::Midplane:
 *
 *     phi_0( R_0( psi ), psi ) = 0
 *
 * with R_0( psi ) each surface's own crossing of z = 0, so n_s0 is the density
 * where the surface crosses the midplane. MIRROR-PLAN.md section 5 measured why:
 * the fixed-radius transfer carries exp( m_s omega^2 ( R_ref^2 - R_0^2 )/2T_s ),
 * an exponent LINEAR IN THE SPECIES MASS, and at M = 6 carbon's n_s0 spans 1e216
 * against a double's exp( 709.8 ); under the midplane gauge that exponent is
 * zero on every surface by construction. It is a mirror's gauge only, because
 * z = 0 crosses a tokamak surface twice and the condition is then not one per
 * surface; a mirror's psi is monotone along the midplane, which is the
 * PRECONDITION tabulateMidplaneCrossing() checks rather than assumes.
 *
 * WHAT IT COSTS IS THE PROPERTY THE PARAGRAPH ABOVE ADVERTISES. R_0( psi ) is a
 * functional of the solution, so F stops being a pointwise function of
 * ( R, z, psi ) alone. The solver tabulates R_0 once per Newton iteration,
 * through NormalisedSource::setMidplaneField(), and FREEZES it within the
 * iteration -- the support freeze's precedent -- so that within one iteration F
 * is again pointwise and every derivative below is exact for it. What the
 * Jacobian does not carry is R_0's dependence on the field, which the freeze
 * leaves to the next iteration. GradShafranovSolver::setMidplaneGaugeKnots()
 * carries it instead: R_0 at a set of knots become unknowns of the Newton,
 * crossingThroughKnots() is the table, and crossingDerivatives() the columns.
 *
 * A FIXED CURVE IS A GAUGE TOO, AND IT KEEPS F POINTWISE: Gauge::FixedCurve,
 *
 *     phi_0( R_g( psi ), psi ) = 0
 *
 * on a curve R_g given once and never moved. A gauge pin does not have to lie
 * on the surface it labels -- quasineutrality fixes phi_0( R, psi ) pointwise
 * for any ( R, psi ), which is why R = R_ref is a gauge on a mirror whose
 * surfaces never reach it -- so any KNOWN R_g( psi ) is one, F is again a
 * pointwise function of ( R, z, psi ), and the derivatives below are exact
 * for it with nothing frozen. The solver gives a curveless normalised source
 * the VACUUM field's own midplane crossing, R_0vac( Psi ), before its first
 * residual and never again: n_s0 is then the density where the vacuum's
 * surface crossed z = 0, and the transfer to the true midplane is
 * exp( m_s omega^2 ( R_0vac^2 - R_0^2 )/2T_s ) -- only the plasma's own
 * displacement of its surfaces, not the machine's span of radii.
 * MIRROR-PLAN.md section 5, [M-VACGAUGE](MEASUREMENTS.md#m-vacgauge).
 *
 * THE PRICE OF THE GAUGE IS THAT TWO SETS OF n_s0 DIFFERING BY IT DESCRIBE THE
 * SAME PLASMA. Anyone comparing profiles against RoPP's own notation, against a
 * gyrokinetic code or against a transport code has to know which convention each
 * is in. That is the one part of the gauge freedom that is not a simplification.
 *
 * WHAT F IS. Differentiating p = Sum_s n_s T_s at fixed R, the dphi_0/dpsi terms
 * collect into -e ( dphi_0/dpsi ) Sum_s Z_s n_s, which vanishes IDENTICALLY by
 * (97). What is left is exactly (136)'s brace plus its omega omega' term, so
 *
 *     F( R, z, psi ) = mu0 R^2 dp/dpsi|_r + g g'
 *
 * which is meq::MHDSource's shape with an R-dependent p. The residual therefore
 * needs phi_0 but never its psi-derivative; only the Jacobian does. Checked two
 * ways: it is RoPP's own force balance (128), and at omega -> 0 it gives RoPP
 * (243), which is the static equation MEQ already solves.
 *
 * UNITS AND SIGN, pinned deliberately because two sign errors have already been
 * found in the papers MEQ does follow. RoPP is Gaussian, with the 4 pi of (136)
 * and the c of (135); MEQ is SI, so 4 pi -> mu0 and RoPP's I is MEQ's g = R B_phi.
 * RoPP (135) is j.grad(phi) = -( c/4 pi R^2 ) Delta* psi, so (136) reads
 * Delta* psi = -4 pi R^2 { ... }, so F = -Delta* psi = +mu0 R^2 { ... } -- the
 * same positive convention meq::MHDSource uses. RoPP's psi is poloidal flux per
 * radian, from its (31) B = I grad(phi) + grad(psi) x grad(phi), which is what
 * MEQ's g g' in T^2 m^2 per Wb/rad already assumes and what EQDSK tabulates.
 */

namespace meq
{

	/// The most species meq::RotatingSource will carry. A cap rather than a
	/// vector so that the per-quadrature-point work allocates nothing and needs
	/// no mutable scratch, which a source evaluated from a threaded assembly
	/// must not have. Eight covers a main ion, a bulk impurity, a seeded
	/// impurity, an alpha population and electrons with room to spare.
	inline constexpr std::size_t maxSpecies = 8;

	/**
	 * One plasma species: its constants, its temperature and its density on the
	 * reference curve.
	 *
	 * A plain value type with no MFEM in it, for the same reason meq::Profile and
	 * meq::Source have none -- so that the physics can be unit tested without the
	 * library.
	 */
	struct Species
	{
		/// Particle mass in kg. Must be finite and positive.
		double mass = 0.0;

		/// Charge in units of e, signed and dimensionless: +1 for a proton, -1
		/// for an electron, +6 for fully stripped carbon. Must be non-zero.
		double charge = 0.0;

		/// T_s( psi ) in JOULES, not eV and not keV. Must not be null, and must
		/// be positive everywhere the solve reaches -- a species with a
		/// non-positive temperature has no Maxwellian.
		std::shared_ptr<Profile const> temperature;

		/// n_s0( psi ) in m^-3: the density of this species ON THE CURVE
		/// R = R_ref -- or, under RotatingSource::Gauge::Midplane, where each
		/// flux surface crosses z = 0, and under Gauge::FixedCurve on R = R_g
		/// ( psi ), the vacuum's crossing -- which in the local gauges of this file
		/// is what N_s of RoPP (96) reduces to. Must not be null.
		std::shared_ptr<Profile const> density;
	};

	/**
	 * THE MIDPLANE CROSSING R_0( x ) OF EVERY FLUX SURFACE x IN [ lower, upper ],
	 * as a Hermite table on uniform knots: each knot's R_0 by a ROOT FIND on
	 * x( R, 0 ) = x_i, and its slope by the inverse function theorem,
	 * dR_0/dx = 1/( dx/dR ) at R_0, with dx/dR read off the field rather than
	 * differenced. MIRROR-PLAN.md section 5.
	 *
	 * @param field  x( R, 0 ) and dx/dR in the flux variable the caller's
	 *               profiles take -- PHYSICAL psi for meq::RotatingSource, Psi
	 *               for the normalised one, which converts before calling this.
	 * @param rMin, rMax  the radii along z = 0 to search.
	 * @param lower, upper  the range of x to tabulate; [ 0, 1 ] for a support.
	 * @param intervals  the knot count less one. The Hermite error is
	 *               h^4 R_0''''/384, and R_0 behaves like sqrt( x ) near the
	 *               inner label; the default holds it far below the
	 *               discretisation's error, which RotatingMirror's midplane
	 *               suite measures by doubling it.
	 * @param samples  how many points along [ rMin, rMax ] bracket the roots
	 *               and judge monotonicity.
	 *
	 * **THE PRECONDITION IS CHECKED, NOT ASSUMED.** The samples carrying x in
	 * [ lower, upper ], with one beyond each end, must lie in the domain without
	 * a gap, x must be STRICTLY monotone across them, and dx/dR must carry that
	 * sense at every one inside: then each surface crosses z = 0 once and R_0 is
	 * a function. A tokamak fails it -- z = 0 crosses each surface inboard and
	 * outboard of the axis -- and so does a mirror whose midplane flux has
	 * folded.
	 *
	 * @throws std::invalid_argument for a non-finite or empty range, or fewer
	 *         than two intervals or samples.
	 * @throws std::runtime_error if the precondition fails, or if the midplane
	 *         inside [ rMin, rMax ] does not reach every x in [ lower, upper ].
	 */
	std::shared_ptr<SplineProfile const> tabulateMidplaneCrossing(
		MidplaneField const &field, double rMin, double rMax,
		double lower, double upper,
		unsigned int intervals = 1024, unsigned int samples = 4096 );

	/**
	 * HOW LARGE F MAY STAND AT A LABEL, against its own peak, before
	 * NormalisedRotatingSource::setVacuumMidplaneField() refuses the profiles:
	 * the bound on RotatingSource::LabelCheck::relative, and the amplitude MR-4a
	 * accepts geq's own state at. What it guards is the ORDER: F jumping across
	 * a label's contour, which cuts elements, caps the convergence rate as a
	 * plasma edge does (CLAUDE_FB.md, *The plasma edge caps the order*). The
	 * Jacobian's missing edge term is the lesser cost -- MEASURED, F at a label
	 * at 1e-3 of its peak leaves Newton quadratic in three steps;
	 * [M-VACGAUGE](MEASUREMENTS.md#m-vacgauge).
	 */
	inline constexpr double labelTolerance = 1.0e-6;

	/**
	 * R_0( x ) THROUGH ITS VALUES AT @a basis's KNOTS, as the bordered midplane
	 * gauge's table (GradShafranovSolver::setMidplaneGaugeKnots()): the C^2
	 * spline of W = R_0^2/2 through W_k = R_0k^2/2, read back as
	 * R_0 = sqrt( 2W ), R_0' = W'/R_0 and R_0'' = ( W'' - R_0'^2 )/R_0. So
	 * RotatingSource's half-delta, which reads R_0^2/2, R_0 R_0' and
	 * R_0'^2 + R_0 R_0'', sees W, W' and W'' to round-off, and
	 * dW/dR_0k = R_0k phi_k exactly.
	 *
	 * W AND NOT R_0 BECAUSE W IS NEARLY LINEAR IN THE FLUX: psi ~ B R^2/2 in a
	 * uniform field, so the interpolant's h^4 error is set by the field's
	 * curvature and the plasma's own, not by the square root.
	 * [M-R0BORDER](MEASUREMENTS.md#m-r0border) measures the two.
	 *
	 * @throws std::invalid_argument for a null basis, a count that is not its
	 *         size, or a crossing that is not finite and positive.
	 */
	std::shared_ptr<Profile const> crossingThroughKnots(
		std::shared_ptr<CardinalSpline const> basis, std::vector<double> const &crossings );

	/**
	 * Sum_s Z_s n_s0( psi ), which charge neutrality on R = R_ref requires to
	 * vanish.
	 *
	 * Exposed because it is the quantity meq::RotatingSource's constructor
	 * checks, so a caller can see why a set was refused, and because it is what a
	 * configuration layer should drive to zero rather than asking the user to.
	 */
	double chargeNeutralityResidual( std::vector<Species> const & species, double psi );

	/**
	 * The density profile that charge neutrality implies for one species, given
	 * all the others.
	 *
	 * WHY THIS EXISTS RATHER THAN A CHECK ALONE. Fixing the gauge removes exactly
	 * one function's worth of freedom from the set of densities -- the
	 * transformation phi_0 -> phi_0 + delta( psi ), n_s0 -> n_s0 exp( Z_s e
	 * delta/T_s ) leaves every physical quantity alone -- so for n species there
	 * are n - 1 independent density flux functions. For two species that is one,
	 * which is the familiar statement. Asking a user to supply n profiles that
	 * happen to satisfy a constraint invites them to supply n that do not.
	 *
	 * The returned profile is -( 1/Z_index ) Sum_{s != index} Z_s n_s0, and its
	 * prime() and doublePrime() are the matching combinations, so it is exact at
	 * every level rather than differenced.
	 *
	 * The nominated species' own density is never read, so it may be null when
	 * this is called -- which is what lets a caller fill the set in one pass.
	 *
	 * @throws std::out_of_range if @a index is not a species.
	 * @throws std::invalid_argument if the nominated species has zero charge, if
	 *         fewer than two species are given, or if any density OTHER than the
	 *         nominated one is null.
	 */
	std::shared_ptr<Profile const> neutralisingDensity( std::vector<Species> const & species, std::size_t index );

	class RotatingSource;

	/**
	 * n_s( R_c( x ), x ) ALONG A CURVE, as a Hermite table on the curve's own
	 * knots: one species' density where the curve meets each surface, its
	 * slope the chain rule dn_s/dpsi + ( dn_s/dR ) R_c' with both partials
	 * analytic. On the solved state's own midplane crossing it is the TRUE
	 * midplane density whatever gauge the source was posed in -- which is the
	 * relabelling the output carries, and the table that poses the same
	 * plasma in Gauge::Midplane. [M-NCDENSITY](MEASUREMENTS.md#m-ncdensity).
	 *
	 * @param curve  R_c at each knot, in @a source's flux variable.
	 * @throws std::out_of_range if @a index is not a species.
	 */
	std::shared_ptr<SplineProfile const> densityOnCurve( RotatingSource const &source,
		std::size_t index, SplineProfile const &curve );

	/**
	 * The rotating Grad-Shafranov source.
	 *
	 * WHY THIS IS NOT A GENERALISED meq::MHDSource. That class stores the
	 * pre-multiplied PRODUCTS p' and g g' precisely so that no chain rule exists
	 * to get wrong: F is one profile evaluation and dF/dpsi is one prime(). The
	 * trick is unavailable here, because p depends on the density, the
	 * temperatures and omega separately and on R as well, so there is no product
	 * to pre-store. Two classes with different invariants, not one class with a
	 * flag -- and meq::MHDSource remains the omega = 0 path.
	 *
	 * IT COSTS A THIRD DERIVATIVE LEVEL, which is why meq::Profile grew
	 * doublePrime(). F is already dp/dpsi, so the Jacobian is d2p/dpsi2 and every
	 * input profile is asked for two derivatives. There is no reparametrisation
	 * that avoids it; see Profile::doublePrime.
	 *
	 * TWO CLOSURES, AND THE CHOICE IS NOT NUMERICAL. With exactly two species
	 * (97) is linear in phi_0 after taking logarithms, so the whole closure is a
	 * closed form: no root find, no inner tolerance, no implicit
	 * differentiation. Three or more needs a safeguarded scalar Newton per
	 * evaluation point and phi_0's two psi-derivatives by implicit
	 * differentiation. Both are implemented and Closure::Automatic picks the
	 * cheap one when it applies; Closure::RootFind can be forced at two species
	 * so that a test can drive the same problem down both paths and require the
	 * same answer, which is what says the general path is right.
	 */
	class RotatingSource : public Source
	{
		public:
			/// How (97) is solved for phi_0.
			enum class Closure
			{
				Automatic,  ///< ClosedForm at two species, RootFind above.
				ClosedForm, ///< Two species only; refused for more.
				RootFind    ///< Any number, and the cross-check at two.
			};

			/**
			 * WHERE phi_0 IS PINNED ON EACH FLUX SURFACE, and so where each
			 * n_s0 is the physical density. See the file comment.
			 */
			enum class Gauge
			{
				/// phi_0( R_ref, psi ) = 0 at a CONSTANT R_ref. The default,
				/// and every tokamak's.
				FixedRadius,
				/// phi_0( R_0( psi ), psi ) = 0 at each surface's own midplane
				/// crossing, tabulated by setMidplaneCrossing(). A mirror's.
				Midplane,
				/// phi_0( R_g( psi ), psi ) = 0 on a FIXED curve, given once by
				/// setMidplaneCrossing() and never replaced; the solver gives a
				/// normalised source the vacuum's midplane crossing. A mirror's,
				/// with an exact Jacobian.
				FixedCurve
			};

			/**
			 * THE BOX'S PRECONDITION ON A CURVE, F -> 0 AT BOTH LABELS, as
			 * labelCheck() measures it. |F|/R is mu0 j_phi.
			 */
			struct LabelCheck
			{
				/// max |F|/R along the curve, at ( R_g( x ), x ) for x in [ 0, 1 ].
				double peak = 0.0;
				/// max |F|/R at x = 0 and at x = 1, over R in ( 0, R_g( x ) ]: each label
				/// at every radius its surface can reach.
				double worst = 0.0;
				/// worst/peak; zero where both vanish.
				double relative = 0.0;
				/// Where worst is: the label, 0 or 1, and the radius.
				double label = 0.0;
				double radius = 0.0;
			};

			/**
			 * @param species  exactly two, of opposite charge sign, satisfying
			 *                 charge neutrality on R = R_ref.
			 * @param omega    the rigid rotation frequency omega( psi ) in rad/s.
			 *                 May be null, which means omega = 0 and reduces this
			 *                 source to meq::MHDSource's equation.
			 * @param ggPrime  the single product g dg/dpsi, in T^2 m^2 per Wb/rad,
			 *                 exactly as meq::MHDSource takes it. Must not be null.
			 * @param referenceRadius  R_ref in metres, where phi_0 vanishes and
			 *                 where each n_s0 is the physical density. Must be
			 *                 finite and positive.
			 * @param mu0      so that a run in normalised units can set it to 1.
			 *
			 * @param closure  which route through (97); see Closure.
			 * @param gauge    where phi_0 is pinned; see Gauge. Under
			 *                 Gauge::Midplane and Gauge::FixedCurve
			 *                 @a referenceRadius MUST BE ZERO: neither has a
			 *                 constant radius, and a radius accepted and not
			 *                 read would be a promise it does not keep. The
			 *                 source then answers nothing until
			 *                 setMidplaneCrossing() has given it its curve.
			 *
			 * @throws std::invalid_argument for a null profile, a non-finite or
			 *         non-positive mass, radius or mu0, a zero charge, fewer than
			 *         two species or more than maxSpecies, a species set whose
			 *         charges are all one sign, Closure::ClosedForm asked for
			 *         with more than two species, or a non-zero radius under
			 *         Gauge::Midplane or Gauge::FixedCurve.
			 * @throws std::invalid_argument if charge neutrality on R = R_ref is
			 *         violated. It is checked by sampling n_s0 over [ 0, 1 ], the
			 *         range meq::Profile documents itself on, because the closed
			 *         form below is derived from Z_1 n_10 = -Z_2 n_20 and is
			 *         simply wrong without it. Use neutralisingDensity() to build
			 *         a conforming set rather than hand-balancing one.
			 */
			RotatingSource( std::vector<Species> species,
				std::shared_ptr<Profile const> omega,
				std::shared_ptr<Profile const> ggPrime,
				double referenceRadius,
				double mu0 = vacuumPermeability,
				Closure closure = Closure::Automatic,
				Gauge gauge = Gauge::FixedRadius );

			/**
			 * R_0( psi ), THE MIDPLANE GAUGE'S CROSSINGS, in the flux variable
			 * this source takes, with its two derivatives -- which is what
			 * the midplane gauge differentiates: the exponent's centrifugal
			 * half is m_s omega^2 ( R^2 - R_0^2 )/2T_s, so dA_s/dpsi gains
			 * -m_s omega^2 R_0 R_0'/T_s, and d2A_s/dpsi2 the next term down.
			 * tabulateMidplaneCrossing() builds the one the solver installs;
			 * any meq::Profile serves, which is what lets a test install a
			 * closed form.
			 *
			 * Replaced once per Newton iteration by the solver, on the master
			 * thread, between assemblies, so a threaded assembly never sees it
			 * move. Held by shared pointer so the replacement is one store.
			 *
			 * UNDER Gauge::FixedCurve IT IS THE CURVE R_g( psi ), given once:
			 * the solver installs the vacuum's crossing only where none has
			 * been given, and replaces nothing.
			 *
			 * @throws std::logic_error under Gauge::FixedRadius, which has no
			 *         crossing to be given.
			 * @throws std::invalid_argument for a null table.
			 */
			void setMidplaneCrossing( std::shared_ptr<Profile const> table );

			/**
			 * F -> 0 AT BOTH LABELS, x = 0 AND x = 1, ON A CURVE: the box's
			 * precondition, which setPlasmaSupport() documents -- where F does
			 * not vanish at an edge of the support it jumps across a contour
			 * that cuts elements, which caps the order, and the Jacobian
			 * misses the edge term F delta( Psi ).
			 *
			 * CHECKABLE BEFORE THE SOLVE ON A FIXED CURVE ONLY. In the gauge
			 * on R_g a density flat at a label still carries
			 * -n_s m_s omega^2 R_g R_g'/T_s in dn_s/dpsi wherever omega does
			 * not vanish there, and on a fixed curve R_g' is known in advance,
			 * so whether a profile set meets the precondition is too. Under
			 * Gauge::Midplane the R_0' at a label is the ANSWER's, and the
			 * vacuum's is not an estimate of it good enough to refuse on --
			 * [M-VACGAUGE](MEASUREMENTS.md#m-vacgauge).
			 *
			 * @param curve  R at each x, in this source's flux variable. The
			 *               source's own table is neither read nor replaced.
			 * @throws std::logic_error under Gauge::FixedRadius.
			 * @throws std::invalid_argument for a null curve.
			 */
			LabelCheck labelCheck( std::shared_ptr<Profile const> curve ) const;

			/// The table setMidplaneCrossing() last gave, or null.
			Profile const *midplaneCrossing() const;

			/**
			 * dF/dW AND dF/dW' AT FIXED ( R, psi ), W( psi ) = R_0( psi )^2/2
			 * the HALF SQUARED CROSSING -- the derivatives of f() with respect
			 * to the table's value and slope at this point's psi, the other
			 * profiles and R held. MIRROR-PLAN.md section 5, and the columns of
			 * the bordered midplane gauge, GradShafranovSolver::
			 * setMidplaneGaugeKnots(): with R_0 at a set of knots as border
			 * unknowns, F depends on each knot through W and W' alone, and the
			 * chain rule through the interpolant is the solver's.
			 *
			 * W AND NOT R_0, BECAUSE THE EXPONENT READS W. The half-delta is
			 * h = R^2/2 - W and h' = -W', so f() reads the table through W and W'
			 * and nothing else -- R_0'' enters dFdPsi() only. Analytic: the
			 * closed form's shared exponent is C h, and the root find's is
			 * differentiated implicitly through (97) as dPotentialDPsi() is.
			 *
			 * @throws std::logic_error under Gauge::FixedRadius, whose F reads no
			 *         crossing, and as halfDelta() does before a table is set.
			 */
			void crossingDerivatives( double radius, double psi, double &dFdW,
			                          double &dFdWPrime ) const;

			/// F at ( R, z ) for the flux value psi, in the units of
			/// meq::Source::f: mu0 R^2 dp/dpsi|_r + g g', no 1/R applied.
			double f( double radius, double z, double psi ) const override;

			/// dF/dpsi at fixed ( R, z ), which is mu0 R^2 d2p/dpsi2|_r plus the
			/// derivative of g g'. Exact, and checked against a central difference
			/// of f() over a Mach sweep in RotatingSourceTests.
			double dFdPsi( double radius, double z, double psi ) const override;

			// ---- the closure, exposed because it is what FL-1 asserts on ----

			/// e phi_0( R, psi ) in JOULES -- the potential of (97) times the
			/// elementary charge, which is the combination that appears in every
			/// exponent and the one that carries no factor of e. Zero at R = R_ref
			/// -- on R_0( psi ) or R_g( psi ) under the other gauges -- by construction,
			/// exactly.
			double potential( double radius, double psi ) const;

			/// n_s( R, psi ) in m^-3, from (96). At R = R_ref, or R_0( psi ), this
			/// is n_s0( psi ).
			/// @throws std::out_of_range if @a index is not a species.
			double density( std::size_t index, double radius, double psi ) const;

			/// n_s( R, psi ) with its two partial derivatives, d/dR at fixed
			/// psi and d/dpsi at fixed R -- the first through d( e phi_0 )/dR,
			/// by implicit differentiation of (97) as dPotentialDPsi() is. What
			/// densityOnCurve() gives a table's slopes from.
			/// @throws std::out_of_range if @a index is not a species.
			void densityDerivatives( std::size_t index, double radius, double psi,
			                         double &value, double &dRadius,
			                         double &dPsi ) const;

			/// d( e phi_0 )/dpsi at fixed R, by implicit differentiation of (97)
			/// rather than by differencing it. Differencing an inner root find
			/// from outside gives a derivative whose accuracy is the inner
			/// tolerance, which is how a Newton iteration silently degrades from
			/// quadratic to linear with no wrong answer anywhere to point at.
			double dPotentialDPsi( double radius, double psi ) const;

			/// The total pressure Sum_s n_s T_s at ( R, psi ), in Pa. This is the
			/// p whose psi-derivative F is built from.
			double pressure( double radius, double psi ) const;

			/// dp/dpsi at fixed R, in Pa per Wb/rad. F is mu0 R^2 times this plus
			/// g g'; exposed separately so a test can see which half is wrong.
			double dPressureDPsi( double radius, double psi ) const;

			/// The exponent of (96) for one species. AT TWO SPECIES ALL SPECIES
			/// SHARE IT, which is what makes Sum_s Z_s n_s vanish at every R once
			/// it vanishes at R_ref; above two they differ, and the cancellation
			/// is the root find's job instead. Dimensionless, and equal to M^2/2
			/// at the outboard edge in the usual Mach-number sense.
			/// @throws std::out_of_range if @a index is not a species.
			double densityExponent( std::size_t index, double radius, double psi ) const;

			// ---- accessors ----

			std::vector<Species> const & species() const;
			Closure closure() const;
			Gauge gauge() const;
			Profile const * omega() const;
			Profile const & ggPrime() const;
			/// R_ref; zero under Gauge::Midplane and FixedCurve, which have none.
			double referenceRadius() const;
			double mu0() const;

		private:
			// Everything the source needs at one ( R, psi ), in one pass, on the
			// stack. f() and dFdPsi() are called at every quadrature point of
			// every element of every residual evaluation, so this allocates
			// nothing and holds no state between calls.
			struct State
			{
				double potential;             // e phi_0
				double potentialPrime;        // and its two psi-derivatives
				double potentialDoublePrime;
				std::array<double, maxSpecies> exponent;          // A_s of (96)
				std::array<double, maxSpecies> exponentPrime;     // dA_s/dpsi
				std::array<double, maxSpecies> exponentDoublePrime;
			};

			State stateAt( double radius, double psi ) const;
			State closedFormState( double radius, double psi ) const;
			State rootFindState( double radius, double psi ) const;

			// ( R^2 - R_*^2 )/2 and its two psi-derivatives at fixed R, R_* the
			// gauge's radius: the constant R_ref, whose derivatives vanish, or
			// the tabulated R_0( psi ), whose do not. The one place the two
			// gauges differ; everything downstream is the chain rule through it.
			void halfDelta( double radius, double psi, double & value,
				double & prime, double & doublePrime ) const;

			// p, dp/dpsi and d2p/dpsi2 at fixed R, from a State.
			void pressureFrom( State const & state, double psi,
				double & p, double & pPrime, double & pDoublePrime ) const;

			std::vector<Species> speciesData;
			std::shared_ptr<Profile const> omegaProfile;
			std::shared_ptr<Profile const> ggPrimeProfile;
			double refRadius;
			double permeability;
			Closure closureChoice;
			Gauge gaugeChoice;
			std::shared_ptr<Profile const> midplaneTable;
	};

	/**
	 * The rotating source with its profiles given against NORMALISED flux.
	 *
	 * Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ), and MEQ is fixed boundary
	 * with psi = 0 on Gamma, so Psi = psi/psi_ax -- which makes psi_ax a
	 * functional of the solution and therefore an unknown of the non-linear
	 * system, not an input. meq::NormalisedSource is the interface the solver
	 * closes by a bordered Newton; see CLAUDE_HDGGS.md, "Newton, and the obligation it
	 * creates", before using this.
	 *
	 * IT IS A WRAPPER AND NOT A REIMPLEMENTATION, which is the whole design. A
	 * normalised source has the form F( R, z, psi ) = H( R, z, psi/psi_ax )/psi_ax,
	 * so with the profiles read at Psi the answer is the unnormalised source's,
	 * divided by psi_ax once for f() and TWICE for dFdPsi() -- once for the
	 * profile argument and once for the derivative itself. Every property
	 * RotatingSourceTests establishes for meq::RotatingSource therefore carries
	 * over unchanged, and there is no second copy of the closure to keep in step.
	 *
	 * The profiles mean what they always meant, read at Psi: n_s0( Psi ) is the
	 * density on R = R_ref, T_s( Psi ) is in Joules, omega( Psi ) is in rad/s, and
	 * ggPrime is d( g^2/2 )/dPsi -- note dPsi, not dpsi, which is what the second
	 * factor of 1/psi_ax accounts for.
	 */
	class NormalisedRotatingSource : public NormalisedSource
	{
		public:
			/// Arguments as meq::RotatingSource's, plus the starting psi_ax.
			/// @throws std::invalid_argument for anything RotatingSource refuses,
			///         or a psi_ax that is not finite or is zero.
			NormalisedRotatingSource( std::vector<Species> species,
				std::shared_ptr<Profile const> omega,
				std::shared_ptr<Profile const> ggPrime,
				double referenceRadius,
				double psiAxis,
				double mu0 = vacuumPermeability,
				RotatingSource::Closure closure = RotatingSource::Closure::Automatic,
				RotatingSource::Gauge gauge = RotatingSource::Gauge::FixedRadius );

			/// True under RotatingSource::Gauge::Midplane, which is what makes
			/// the solver call setMidplaneField().
			bool readsMidplane() const override;

			/**
			 * TABULATE R_0( Psi ) OVER [ 0, 1 ] FROM THE FIELD ALONG z = 0, at
			 * the normalisation in force, and hand it to the source underneath.
			 *
			 * IN Psi, NOT IN psi, AND THAT KEEPS THIS CLASS A PURE WRAPPER. The
			 * inner source's argument is Psi, so a table in Psi keeps
			 * F = S H( Psi )/span with H a function of Psi alone, and
			 * normalisationDerivatives() stays exact for the frozen table
			 * without a term in R_0. Within an iteration it is R_0( Psi ) that
			 * is frozen; at convergence the two readings are one table.
			 *
			 * @throws std::logic_error under Gauge::FixedRadius.
			 * @throws std::runtime_error as tabulateMidplaneCrossing() does, for
			 *         a midplane on which psi is not monotone across the support.
			 */
			void setMidplaneField( MidplaneField const &field, double rMin,
			                       double rMax ) override;

			/// R_0( Psi ) handed to the source underneath as it stands -- the
			/// bordered gauge's interpolant through its knots. In Psi, for the
			/// reason setMidplaneField() tabulates in Psi.
			/// @throws std::logic_error under Gauge::FixedRadius.
			void setMidplaneCrossing( std::shared_ptr<Profile const> crossing ) override;

			/// RotatingSource::crossingDerivatives() at Psi, times the scale and
			/// over the span as f() is; zero outside the plasma.
			/// @throws std::logic_error under Gauge::FixedRadius.
			bool crossingDerivatives( double radius, double z, double psi,
			                          double &dFdW, double &dFdWPrime ) const override;

			/// Both honour setPlasmaSupport(): `F` and `dF/dpsi` are zero
			/// wherever `Psi <= 0`, as meq::NormalisedMHDSource's are -- and,
			/// under setPlasmaSupportShape( PlasmaSupportShape::Box ), wherever
			/// `Psi > 1` or `|z|` is beyond the throat, which is the support a
			/// MIRROR needs and the reason this source is the one the mirror
			/// path uses. normalisationDerivatives() is gated identically, so
			/// the three never disagree about where the plasma is.
			///
			/// THE DIAGNOSTICS BELOW DELIBERATELY DO NOT. potential(), density()
			/// and pressure() report what the profiles say at a flux value, which
			/// is a question about the profiles and not about where the plasma
			/// is; confining them would make a caller unable to ask.
			///
			/// BOTH CARRY setCurrentScale()'s FACTOR, as
			/// meq::NormalisedMHDSource's do. The contract is "scale the plasma
			/// term", and `[source] Type = "rotating"` accepts `PlasmaCurrent`,
			/// so a source that ignored it would give a current border that
			/// cannot respond to its own unknown.
			double f( double radius, double z, double psi ) const override;
			double dFdPsi( double radius, double z, double psi ) const override;

			/// ANALYTIC, AND THE CLOSURE IS NOT DIFFERENTIATED AGAIN TO GET IT.
			///
			/// This class is a pure wrapper -- `F = S H( Psi )/span` with
			/// `H = inner.f` and `inner.dFdPsi` already equal to `dH/dPsi` -- so
			/// the normalisation enters only through `Psi` and the overall
			/// `1/span`, and meq::NormalisedMHDSource's algebra applies verbatim.
			/// Whatever the species root find cost to differentiate has been
			/// paid inside `inner.dFdPsi()`.
			///
			/// It matters because the alternative is a DIFFERENCE, which
			/// perturbs the normalisation, which moves the edge, so it straddles
			/// a kink rather than measuring a derivative -- and because
			/// GradShafranovSolver::assembleCurrentNormalisationCorner() reaches
			/// for it too, where until it existed the fallback was a silent zero
			/// rather than a difference.
			bool normalisationDerivatives( double radius, double z, double psi,
			                               double &dFdAxis,
			                               double &dFdBoundary ) const override;

			/// @throws std::invalid_argument if psi_ax is not finite or is zero:
			///         Psi = psi/psi_ax is undefined there, and a solver that has
			///         wandered onto psi_ax = 0 should say so rather than return
			///         infinities.
			void setNormalisation( double psiAxis, double psiBoundary ) override;
			using NormalisedSource::setNormalisation;
			double normalisation() const override;
			double boundaryNormalisation() const override;

			// Diagnostics, taking PHYSICAL psi and converting internally, so that
			// a caller never has to remember which of the two it is holding.

			double potential( double radius, double psi ) const;
			double density( std::size_t index, double radius, double psi ) const;
			double pressure( double radius, double psi ) const;

			/// True under RotatingSource::Gauge::FixedCurve while the source
			/// underneath has no curve, which is what makes the solver call
			/// setVacuumMidplaneField() -- once.
			bool readsVacuumMidplane() const override;

			/**
			 * TABULATE THE GAUGE CURVE R_g( Psi ) OVER [ 0, 1 ] FROM THE
			 * VACUUM ALONG z = 0, NORMALISED BY THE VACUUM'S OWN LABELS, check
			 * the box's precondition on it, and give it to the source
			 * underneath for good.
			 *
			 * In Psi, for setMidplaneField()'s reason: the class stays a pure
			 * wrapper, so normalisationDerivatives() is exact -- and with the
			 * curve fixed there is nothing frozen for it to leave out. Against
			 * @a psiAxis and @a psiBoundary, the vacuum's flux at the label
			 * points, so R_g( 1 ) and R_g( 0 ) ARE those points where they lie
			 * on z = 0, whatever normalisation the solve starts from.
			 *
			 * @throws std::logic_error under any gauge but FixedCurve, or once
			 *         the curve is in.
			 * @throws std::runtime_error as tabulateMidplaneCrossing() does.
			 * @throws std::invalid_argument if labelCheck() on the curve finds
			 *         F at a label above meq::labelTolerance of its peak, naming
			 *         the label; the curve is then NOT installed.
			 */
			void setVacuumMidplaneField( MidplaneField const &field, double rMin,
			                             double rMax, double psiAxis,
			                             double psiBoundary ) override;

			/// The source underneath, whose arguments are Psi rather than psi.
			RotatingSource const & unnormalised() const;

			/// The same, to install a midplane table by hand -- what a test of
			/// the gauge does where there is no solver to call
			/// setMidplaneField().
			RotatingSource & unnormalised();

		private:
			RotatingSource inner;
			double psiAxisValue;
			double psiBoundaryValue = 0.0;
	};

}

#endif // MEQ_ROTATINGSOURCE_HPP
