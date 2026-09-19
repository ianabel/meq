#ifndef MEQ_COILS_HPP
#define MEQ_COILS_HPP

#include <cstddef>
#include <memory>
#include <vector>

#include "Source.hpp"

/*
 * Poloidal field coils of rectangular cross-section: stage FB-2 of
 * FREE-BOUNDARY-PLAN.md, and the thing FB-1 needs before it can run at all.
 *
 * WHAT THIS IS. A coil is a region of the ( r, z ) plane carrying a prescribed
 * toroidal current. Its current is DATA -- it does not depend on psi, it is not
 * an unknown, and nothing here iterates. Two things come out: the
 * Grad-Shafranov source term the assembly needs, and the EXACT field the coils
 * produce, which is what FB-1's acceptance compares a vacuum solve against.
 *
 * MFEM-FREE, DELIBERATELY, like Profiles, Source, Zernike, SurfaceFit and
 * ExteriorDtN: plain doubles in and out, geometry and special functions and
 * nothing else. That is what lets CI -- which cannot obtain the MFEM branch MEQ
 * needs -- build and test it, and it is why FB-2's coil half can have a unit
 * test rather than waiting on a convergence study CI cannot build. A
 * mfem::Coefficient wrapping a CoilSet belongs with the assembly that needs it,
 * exactly as it does for Source.
 *
 *
 * THE SOURCE TERM, DERIVED RATHER THAN QUOTED
 * -------------------------------------------
 *
 * FREE-BOUNDARY-PLAN.md section 5.4 prints F_coil = mu0 r I_k / |Omega_ck| and
 * this file does not take that on trust, because CLAUDE.md records two separate
 * occasions in this tree where a transcribed formula converged beautifully to
 * the wrong function. Here is the derivation; it agrees with the plan, and the
 * agreement is reported rather than assumed.
 *
 * MEQ solves -div_bar( ( 1/r ) grad_bar psi ) = F( r, z, psi )/r, and
 *
 *     div_bar( ( 1/r ) grad_bar psi )
 *         = d_r( ( 1/r ) d_r psi ) + d_z( ( 1/r ) d_z psi )
 *         = ( 1/r )[ d_rr psi - ( 1/r ) d_r psi + d_zz psi ]
 *         = ( 1/r ) Delta* psi,
 *
 * so the equation is Delta* psi = -F. That is MEQ's convention and it is the
 * one meq::Source::f() answers in: F as eq (2) writes it, with NO 1/r applied.
 *
 * Now the physics. MEQ's field convention, from CLAUDE.md, is
 * B = ( -q_z, +q_r ) with q = ( 1/r ) grad_bar psi, i.e.
 *
 *     B_r = -( 1/r ) d_z psi,      B_z = +( 1/r ) d_r psi,
 *
 * which is psi = r A_phi, the poloidal flux PER RADIAN. Ampere's law in its
 * toroidal component is mu0 j_phi = ( curl B )_phi = d_z B_r - d_r B_z, and
 *
 *     d_z B_r = -( 1/r ) d_zz psi,
 *     d_r B_z = ( 1/r ) d_rr psi - ( 1/r^2 ) d_r psi,
 *
 * so ( curl B )_phi = -( 1/r )[ d_rr psi - ( 1/r ) d_r psi + d_zz psi ]
 *                   = -( 1/r ) Delta* psi. Hence
 *
 *     Delta* psi = -mu0 r j_phi,     and comparing with Delta* psi = -F,
 *
 *     ===================================================================
 *     F( r, z ) = mu0 r j_phi( r, z ),   so   F_coil = mu0 r I / |Omega_c|
 *     ===================================================================
 *
 * on a coil of area |Omega_c| carrying a total current I. THAT IS THE PLAN'S
 * FORMULA, arrived at independently; there is no disagreement to report.
 *
 * The independent check that it is the same convention meq::Source uses:
 * MHDSource has F = mu0 r^2 p' + g g', and the equilibrium current is
 * j_phi = r p' + g g'/( mu0 r ), whose mu0 r j_phi is mu0 r^2 p' + g g'
 * exactly. So the coil term and the plasma term are in the same units and add.
 *
 * AN EXTRA OR MISSING r HERE CONVERGES AT FULL ORDER TO THE WRONG FUNCTION and
 * no rate table can see it, which is the failure CLAUDE.md records happening
 * twice in this tree already. The check that CAN see it is Delta* of the exact
 * field below, by central differences, against -f(). Measured on a coil at
 * r = 2.0 carrying 1 MA over 0.30 x 0.40 m, so that F = mu0 r J = 2.0944e+01 at
 * its centre; worst residual over nine sample points as a fraction of that F:
 *
 *     h        interior plain   interior Richardson   exterior Richardson
 *     2e-2        3.66e-03           1.24e-06              6.87e-09
 *     1e-2        9.16e-04           2.01e-08              1.83e-10
 *     5e-3        2.29e-04           3.36e-10              9.11e-11
 *
 * The plain column falls at exactly rate 2, which says the floor there is the
 * DIFFERENCE and not the quadrature -- so Richardson extrapolation sees past
 * it, as it does everywhere else in this tree.
 *
 * AND THE INTERIOR COLUMN IS THE ONE THAT PINS THE CONSTANT, which is why the
 * quadrature had to be made accurate inside the coil at all. Outside a coil the
 * equation is Delta* psi = 0, and any multiple of psi satisfies that; only
 * where F is non-zero does the residual say anything about the factor in front
 * of it. The control: at ( 2.06, 0.36 ) the extrapolated Delta*_FD reads
 * -2.157227e+01 against -F = -2.157227e+01, while -F/r would be
 * -1.047198e+01 and -F r would be -4.443888e+01. A missing or extra r is a
 * factor of 2.06 here, not a small bias -- but it is a factor that only a
 * comparison against a closed form can see, since it leaves every convergence
 * rate untouched.
 *
 *
 * THE SOURCE IS PIECEWISE CONSTANT AND DISCONTINUOUS, AND THAT COSTS SOMETHING
 * ---------------------------------------------------------------------------
 *
 * j_phi is I/|Omega_c| inside a coil and zero outside it, so f() JUMPS at every
 * coil edge. It is not a smooth right hand side and it is not meant to be: a
 * real conductor has a real boundary.
 *
 * What that costs is convergence rate. A discontinuous right hand side puts the
 * solution in H^{2+something small} rather than in the H^{k+2} an optimal
 * order needs, so unless the mesh RESOLVES the coil boundaries -- element
 * edges lying on them -- the achievable rate is limited by the cut elements and
 * not by k. That is the same species of statement as CLAUDE.md's warning about
 * the re-entrant corner of a rectangle, and it has the same remedy: align the
 * mesh, or accept the rate.
 *
 * IT IS ALSO WHY FB-2'S ACCEPTANCE IS THE OUTWARD FLUX AGAINST THE TOTAL
 * CURRENT rather than a rate against a closed form. That test is exact whatever
 * the mesh does inside the coil -- it is Ampere's law over the whole domain,
 * one number, and totalCurrent() is the number it is checked against.
 *
 *
 * THE EXACT FIELD
 * ---------------
 *
 * Delta* is linear and the coil is a continuum of filament loops, so
 *
 *     psi_coil( r, z ) = ( I / |Omega_c| ) * integral over the cross-section of
 *                        psi_filament( r, z ; r', z' ) dr' dz'
 *
 * with psi_filament the flux of a UNIT-current loop of radius r' at height z'.
 * That is the free-space field: it decays at infinity and it is what a coil
 * outside the computational domain actually produces.
 *
 * THE FILAMENT KERNEL IS IMPLEMENTED HERE RATHER THAN TAKEN FROM
 * tests/analytic/CurrentLoop.hpp, and the reason is not taste: this is library
 * code and a library must not depend on a test fixture. C++17 has the complete
 * elliptic integrals in <cmath> and Boost.Math -- which meq_core already links
 * for Zernike and ExteriorDtN -- has Carlson's symmetric forms, so there is
 * nothing to depend on. The two implementations are independent and agree; the
 * agreement is a check rather than a duplication, in the same way
 * SolovievGeometryConvergence checks coefficients the solver also uses.
 *
 * IT IS CARLSON'S FORM AND NOT std::comp_ellint, AND THAT IS THE ONE CHOICE IN
 * THIS FILE THAT MATTERS. K( k ) and E( k ) are needed with k -> 1, which is
 * where the field point approaches a source filament -- which, inside a coil,
 * is exactly where the quadrature has to put its points. Written the textbook
 * way, k = 2 sqrt( a r )/d, k rounds to exactly 1.0 at about 1e-8 of a coil
 * radius and std::comp_ellint_1( 1.0 ) is NaN while std::comp_ellint_2 throws.
 * Carlson's forms take the COMPLEMENTARY modulus squared,
 *
 *     k'^2 = ( ( a - r )^2 + ( z - z0 )^2 ) / d^2,
 *
 * which is the squared distance to the source over d^2 -- formed with no
 * cancellation whatever, since nothing is subtracted from one -- and
 * K = R_F( 0, k'^2, 1 ), E = K - ( k^2/3 ) R_D( 0, k'^2, 1 ). Measured, the two
 * routes agree to 7.0e-13 or better over k in [ 1e-8, 1 - 1e-12 ], and Carlson
 * keeps working to k'^2 = 1e-300 where the textbook form has been NaN for
 * nearly three hundred orders of magnitude. See Coils.cpp.
 *
 * THE TWO IMPLEMENTATIONS WERE COMPARED, AND WHERE THEY DISAGREE THIS ONE IS
 * RIGHT. Against tests/analytic/CurrentLoop.hpp -- an independent transcription
 * of the textbook form, written for FB-0 and verified there -- filamentPsi()
 * agrees to a worst 1.4e-12 relative over 4813 points of the benchmark box with
 * a disc of radius 0.2 excluded about the loop, reaching k = 0.9959. Walking in
 * to the conductor on the equator, psi at a distance eps outboard:
 *
 *     eps        CurrentLoop       this file      -(1/2) ln eps
 *     1e-03      3.4955956e+00   3.4955956e+00     3.4538776e+00
 *     1e-05      5.7962147e+00   5.7962150e+00     5.7564627e+00
 *     1e-07      8.0828343e+00   8.0987690e+00     8.0590478e+00
 *     1e-09              NaN     1.0401354e+01     1.0361633e+01
 *     1e-13              NaN     1.5006924e+01     1.4966803e+01
 *
 * -- and the third column settles the disagreement at 1e-07: the field there is
 * a line-current logarithm plus a constant, and this file's answer sits a
 * constant 0.0397 above -(1/2) ln eps at every one of the last four rows while
 * CurrentLoop's has already drifted 0.2% low before going NaN. That fixture
 * says so itself; the point of repeating it here is that the coil quadrature
 * evaluates in exactly that region, so the textbook form was not an option.
 *
 *
 * WHAT THE QUADRATURE ACHIEVES, INSIDE AND OUTSIDE
 * -----------------------------------------------
 *
 * Tensor Gauss-Legendre, with two refinements that are not decoration.
 *
 *   * THE RECTANGLE IS SPLIT AT THE FIELD POINT'S OWN COORDINATES, clamped to
 *     the rectangle, into at most four panels. That puts the singularity at a
 *     panel CORNER instead of in a panel interior, and for a point outside it
 *     aligns a panel edge with the nearest approach.
 *   * EACH PANEL IS THEN GRADED CUBICALLY TOWARD THAT CORNER, t = L u^3.
 *
 * The integrand is weakly singular when the field point is INSIDE the coil:
 * psi_filament has a logarithmic singularity where the source meets the field
 * point, exactly as a two-dimensional line current does. Ungraded, that costs
 * the rule its spectral accuracy and leaves it ALGEBRAIC. Measured on a coil at
 * r = 2.0, half-width 0.15, half-height 0.20, relative error against a
 * converged reference:
 *
 *     grading      interior rate in n      interior error at n = 32
 *     none                 3.97                     1.1e-07
 *     quadratic            7.90                     6.5e-11
 *     CUBIC (shipped)   about 10                    3.0e-12
 *
 * -- an observed law of about 4p for exponent p, and cubic grading reaches
 * ROUND-OFF by n = 48 where the ungraded rule is still at 1e-8. Outside the
 * coil the integrand is analytic and the rule is spectral either way: 5 cm
 * clear of the edge it is at round-off by n = 16, and far away by n = 12.
 * Cubically graded, in relative error against a converged reference:
 *
 *     n        inside      1 mm outside   5 cm outside   far outside
 *      8      1.80e-06       1.37e-06       1.56e-07       4.47e-13
 *     16      1.08e-09       1.63e-09       1.39e-14       6.0e-15
 *     32      2.99e-12       9.24e-14       8.7e-15        8.8e-15
 *     48      4.62e-14       2.3e-15        6.1e-15        7.5e-15
 *
 * SO THE CONTRACT IS: MACHINE PRECISION OUTSIDE, AND ABOUT 3e-12 INSIDE AT THE
 * DEFAULT ORDER, degrading to 1e-9 at n = 16 and reaching round-off at n = 48.
 * That is good enough to use the interior field as a reference, which is more
 * than the plan asked for -- and it is good enough that the Delta* check above
 * can be run INSIDE the coil, which is the only place it pins the constant.
 * A point sitting almost exactly on a corner of the coil, ( 2.1499, 0.4999 )
 * against a corner at ( 2.15, 0.50 ), behaves no differently: 2.42e-06,
 * 2.01e-09, 5.08e-12, 4.47e-14 down the same sequence.
 *
 * WHY THE GRADING EXPONENT IS FIXED AT 3 AND NOT EXPOSED. Quartic grading is
 * better on paper and does not work: at n = 96 the innermost node is within
 * 1e-16 of the field point in ABSOLUTE terms, so adding the offset to a
 * coordinate of order 2 rounds to the coordinate itself, the node lands exactly
 * on the field point, and k'^2 is exactly zero. The limit is floating point,
 * not mathematics. Cubic grading has that margin at every order this class
 * accepts, which is why the order is capped as well; and defensively, a node
 * whose distance underflows is dropped rather than evaluated -- see Coils.cpp.
 *
 * COST, on this machine, per point per coil, dominated entirely by the two
 * Carlson evaluations per node:
 *
 *     n        outside    inside
 *     16       14.0 us    70.0 us
 *     32       56.5 us   275.5 us
 *     48      126.0 us   644.4 us
 *
 * Outside is a quarter of inside because the panel split produces one panel
 * rather than four. This is a reference field and not an inner loop, and it is
 * priced accordingly -- but a consumer sweeping a grid should set the order
 * from the table above rather than leaving it at the default out of habit.
 *
 *
 * A NAMING WARNING. psi() here is the POLOIDAL FLUX, weber per radian, the
 * quantity the solver calls psi. It is NOT the HDG flux q = grad_bar psi / r,
 * which the rest of this tree also calls "the flux" and which
 * GradShafranovSolver::flux() returns. tests/analytic/CurrentLoop.hpp draws the
 * same line with the same two names, and this file follows it.
 */

namespace meq
{

	/// The Gauss order coilPsi() uses unless told otherwise: points per
	/// direction per panel, so up to 4 n^2 kernel evaluations.
	///
	/// 32 is chosen from the measured table in the file comment -- machine
	/// precision outside the coil and about 1e-12 inside it, which makes the
	/// interior usable as a reference rather than merely defined. 16 is enough
	/// for an exterior-only consumer at a quarter of the cost, and FB-2
	/// evaluates the exterior.
	inline constexpr int defaultCoilQuadratureOrder = 32;

	/**
	 * The order an INITIAL GUESS wants, which is not the order a field wants.
	 *
	 * **THE DEFAULT IS FOR A REFERENCE AND IT IS RUINOUS FOR A GUESS.** A guess
	 * built from `[[coils]]` evaluates every conductor at every nodal point of
	 * the potential space AND the trace space, and at order 32 that is
	 * **1392 us a point** for MAST-U's 23 conductors -- measured, and it made
	 * the conductor field **93.6% of that machine's whole run**, dwarfing the
	 * solve it was helping.
	 *
	 * 6 is chosen from a measurement rather than from taste. Against order 32
	 * over 1176 points of MAST-U's disc, where `max |psi_coil|` is
	 * 1.515779e-01 Wb/rad:
	 *
	 *     order   max abs err    relative    us/point
	 *         2   4.317102e-03   2.848e-02        6.6
	 *         4   6.521601e-04   4.302e-03       23.2
	 *         6   9.055783e-05   5.974e-04       50.0
	 *         8   1.168790e-05   7.711e-05       88.1
	 *        16   4.775856e-08   3.151e-07      342.7
	 *        32              --          --     1392.4
	 *
	 * At 6 the error is 9.1e-05 Wb/rad against a `psi_ax` of about 9.2e-02 --
	 * a thousandth of the quantity being guessed, and far below the 1.6e-04 m
	 * the X-point seed is already allowed to be wrong by. **28x cheaper than
	 * the default.** 4 would also do and 2 probably would; 6 is taken because
	 * the cost of being wrong here is a different equilibrium, and three orders
	 * of magnitude of headroom is worth 27 us a point.
	 *
	 * **THE GUESS IS NOT THE PROBLEM STATEMENT**, which is what makes this a
	 * free choice at all: nothing about the answer depends on it, only which
	 * basin Newton starts in. M-124 measures the same indifference from the
	 * other side -- a plasma column tuned to MAST-U and a round default reach
	 * the same `psi_ax` to every printed digit.
	 */
	inline constexpr int guessCoilQuadratureOrder = 6;

	/// The largest order coilPsi() will accept.
	///
	/// A REFUSAL RATHER THAN A LIMIT NOBODY MENTIONS. Cubic grading places the
	/// innermost node at about L( 1.4/n^2 )^3 from the field point, and once
	/// that falls below an ulp of the coordinate the node lands ON the field
	/// point and the kernel is infinite there. 256 keeps a wide margin for any
	/// coil geometry, and it is five times past the order at which the
	/// quadrature reaches round-off, so nothing useful is being refused.
	inline constexpr int maximumCoilQuadratureOrder = 256;

	/**
	 * One coil of rectangular cross-section, carrying a prescribed total
	 * current uniformly distributed over it.
	 *
	 * Geometry is ( centre, half-width in r, half-height in z ) rather than
	 * ( rMin, rMax, zMin, zMax ) because a coil is specified by where it is and
	 * how big it is, and because the half-widths are what the refusals are
	 * about. Both forms are available through the accessors.
	 *
	 * The current is SIGNED and may be zero: a coil carrying no current is a
	 * perfectly ordinary thing for a machine description to contain, and
	 * refusing it would make a scan over currents impossible to write. What is
	 * refused is a geometry that cannot be solved on -- see the constructor.
	 *
	 * IT IS NOT A meq::Source, deliberately. Source::f() takes psi, and a coil
	 * current does not depend on psi; the honest signature is f( r, z ) and it
	 * is on CoilSet. A Source adapter would have dFdPsi() identically zero,
	 * which makes a Newton solve driven by coils alone affine -- one step,
	 * exactly as Soloviev.hpp and CurrentLoop.hpp are -- and that adapter
	 * belongs with the assembly, not here.
	 */
	class Coil
	{
		public:
			/// @param centreRIn     the coil centre's major radius, metres.
			/// @param centreZIn     its height, metres.
			/// @param halfWidthIn   half the extent in r. Strictly positive.
			/// @param halfHeightIn  half the extent in z. Strictly positive.
			/// @param currentIn     the TOTAL current through the
			///                      cross-section, amperes. Signed; zero is
			///                      allowed.
			///
			/// @throws std::invalid_argument if any argument is not finite, if
			///         either half-extent is not positive, or if the coil
			///         reaches the axis: centreR - halfWidth <= 0 is refused
			///         because the Grad-Shafranov operator's 1/r is not
			///         integrable through r = 0, which is the same refusal
			///         meq::BoundaryShape makes and for the same reason.
			Coil( double centreRIn, double centreZIn, double halfWidthIn,
			      double halfHeightIn, double currentIn );

			double centreR() const;
			double centreZ() const;
			double halfWidth() const;
			double halfHeight() const;

			/// The total current through the cross-section, in amperes.
			double current() const;

			/// The bounding box, which is the coil itself.
			double rMin() const;
			double rMax() const;
			double zMin() const;
			double zMax() const;

			/// The cross-sectional area, 4 * halfWidth * halfHeight, in m^2.
			/// Strictly positive by construction.
			double area() const;

			/// The toroidal current density j_phi = I/area, in A/m^2. Uniform
			/// over the coil, which is the modelling assumption this whole
			/// class rests on: a real winding has turns and this does not.
			double currentDensity() const;

			/// Whether ( r, z ) is in the coil.
			///
			/// THE SET IS CLOSED -- the edges belong to it. That matters not at
			/// all for anything integrated, the boundary having measure zero,
			/// and it matters for a point-sampled f(): two coils sharing an
			/// edge both claim it, and a quadrature point landing exactly on an
			/// edge gets the full interior value rather than half of it.
			/// Neither is a defect and both are consequences of the source
			/// being genuinely discontinuous there.
			bool contains( double r, double z ) const;

		private:
			double centreRValue;
			double centreZValue;
			double halfWidthValue;
			double halfHeightValue;
			double currentValue;
	};

	/**
	 * The poloidal flux psi = r A_phi, in weber per radian, of a single
	 * circular filament of radius @a loopRadius at height @a loopHeight
	 * carrying @a current, evaluated at ( r, z ).
	 *
	 * Public because it is the exact kernel the cross-section integral is built
	 * from, because the thin-coil limit is checked against it, and because a
	 * caller modelling a coil as a filament -- which section 5.4 offers as the
	 * alternative to a subdomain -- wants exactly this and nothing more.
	 *
	 * EXACT AND FREE OF CANCELLATION AT THE CONDUCTOR: see the file comment for
	 * why it is Carlson's form. Two consequences worth knowing.
	 *
	 *   * psi( 0, z ) is 0.0 BIT EXACTLY, because the expression carries k^2 as
	 *     an explicit factor and k^2 = 4 a r/d^2 is exactly zero on the axis.
	 *     That is the boundary condition the free-boundary problem imposes
	 *     there, so it is worth having exactly rather than to round-off.
	 *   * NEAR the axis the RELATIVE accuracy degrades, because the bracket
	 *     ( 1 - k^2/2 )K - E vanishes like k^4 and is formed by cancellation.
	 *     Measured on a 1 A loop of radius 2.0 at z = 0.3 with the field point
	 *     at z = 0.7, this form and the textbook one agree to 1.6e-12 at
	 *     r = 1e-2 and to only 8.1e-5 at r = 1e-6 -- but psi itself is
	 *     1.481e-19 at that second point, so the ABSOLUTE disagreement is
	 *     1e-23 and nothing that consumes psi can see it. This is a property of
	 *     the closed form and not of either implementation;
	 *     tests/analytic/CurrentLoop.hpp records the same thing, and it is the
	 *     opposite end of the range from the k -> 1 breakdown above, where the
	 *     two forms genuinely differ.
	 *
	 * @throws std::invalid_argument if any argument is not finite, if
	 *         @a loopRadius is not positive, if @a r is negative, or if the
	 *         field point is exactly ON the loop, where psi is genuinely
	 *         infinite. The last is a refusal rather than an infinity because a
	 *         caller who has hit it has a geometry error, not a large number.
	 */
	double filamentPsi( double r, double z, double loopRadius,
	                    double loopHeight, double current,
	                    double mu0 = vacuumPermeability );

	/**
	 * grad_bar( psi ) = ( d_r psi, d_z psi ) of a single filament, analytic.
	 *
	 * The same geometry, current and conventions as filamentPsi(), and the same
	 * refusals. Not a difference: the chain rule is carried through in closed
	 * form and Carlson's symmetric forms are used throughout, for the reason the
	 * file comment gives for psi.
	 *
	 * THE ALGEBRA, BECAUSE IT IS WHERE THE ACCURACY IS WON. With
	 * A( k ) = ( 1 - k^2/2 )K - E and psi = ( mu0 I/2 pi ) d A, the chain rule
	 * needs dA/dk, and the standard dK/dk and dE/dk collapse it to
	 * dA/dk = ( k/2 )[ E/( 1 - k^2 ) - K ]. Substituting the geometry and
	 * letting the 1/( 2r ) in dk/dr cancel against the 2ar in d k^2 gives
	 *
	 *     d_r psi = ( mu0 I/2 pi )[ ( ( a + r )/d ) A
	 *                               + ( a( a^2 - r^2 + dz^2 )/d^3 ) B ]
	 *     d_z psi = ( mu0 I/2 pi )[ ( dz/d ) A - ( 2 a r dz/d^3 ) B ]
	 *
	 * with B := E/k'^2 - K, dz := z - z0. **Both bracketed forms are written
	 * that way to avoid a cancellation**, and neither is the arrangement the
	 * chain rule hands you:
	 *
	 *   * `a/d - 2ar( a + r )/d^3` is the coefficient of B as it falls out, and
	 *     it VANISHES at the loop -- which is the cancellation that makes the
	 *     gradient diverge like 1/distance rather than like 1/distance^2, since
	 *     B itself goes as 1/k'^2. Factored to `a( a^2 - r^2 + dz^2 )/d^3` it is
	 *     one subtraction of two comparable numbers instead of a difference of
	 *     products, and it is EXACTLY zero at r = a, dz = 0.
	 *   * B in Carlson's forms is `k^2( R_F - R_D/3 )/k'^2`, which carries k^2
	 *     as an explicit factor and is therefore exactly 0.0 at k = 0. Written
	 *     as `E/k'^2 - K` it is a difference of two numbers both near pi/2 on
	 *     the axis and in the far field, and loses its figures there.
	 *
	 * **SO THIS IS FINITE ON THE AXIS WHERE tests/analytic/CurrentLoop.hpp IS
	 * NaN**, and that is a difference from the fixture this was promoted from
	 * rather than an oversight in it. That fixture writes
	 * dk/dr = k[ 1/( 2r ) - ( a + r )/d^2 ] literally and its own header records
	 * the consequence -- *"the 1/( 2r ) is why this is NaN at r = 0. It is a
	 * real 1/r and not an artefact: psi ~ r^2 there, so d psi/d r ~ r and the
	 * limit exists, but the expression as written does not reach it."* The
	 * limit is what the factored form reaches: both A and B carry k^2, k^2 = 0
	 * on the axis, and `d_r psi( 0, z )` is 0.0 bit exactly. `d_z psi( 0, z )`
	 * is 0.0 as well and was already.
	 *
	 * IT STILL GIVES OUT AT THE CONDUCTOR, and sooner than psi does. B goes as
	 * 1/k'^2 and k'^2 is the squared distance to the filament over d^2, so the
	 * gradient diverges -- correctly, this being a line current. Carlson keeps
	 * the geometry cancellation-free far past where the textbook form is NaN,
	 * but the physical divergence is real. CurrentLoop.hpp's practical rule
	 * holds: keep evaluation points at least 1e-5 of a loop radius away if the
	 * gradient is wanted, 1e-7 if only psi is.
	 *
	 * @throws std::invalid_argument on the same conditions as filamentPsi().
	 */
	void filamentGradPsi( double r, double z, double loopRadius,
	                      double loopHeight, double current,
	                      double &dPsiDr, double &dPsiDz,
	                      double mu0 = vacuumPermeability );

	/**
	 * The HDG flux q = ( 1/r ) grad_bar( psi ) of a single filament.
	 *
	 * THIS is the primitive a free-boundary coupling wants: the transmission
	 * condition of FREE-BOUNDARY-PLAN.md section 4.2 is written in q, and its
	 * Neumann half needs q . nu on Gamma.
	 *
	 * **NaN ON THE AXIS IN BOTH COMPONENTS, AND THAT IS NOT THE SAME STATEMENT
	 * AS grad_bar psi's.** grad_bar psi is 0.0 there; dividing 0 by 0 is not.
	 * The limit exists and is finite -- q_r -> mu0 I a^2/( 2 d^3 ) -- but this
	 * expression does not reach it, exactly as CurrentLoop::flux() does not.
	 * **A caller sweeping Gamma must know this**, because a semicircle centred
	 * on the axis MEETS the axis at both ends: the endpoints of such a contour
	 * are precisely where q is unavailable, and a quadrature that samples them
	 * will get NaN rather than a large number.
	 *
	 * @throws std::invalid_argument on the same conditions as filamentPsi().
	 */
	void filamentFlux( double r, double z, double loopRadius,
	                   double loopHeight, double current,
	                   double &qR, double &qZ,
	                   double mu0 = vacuumPermeability );

	/**
	 * The poloidal flux of one rectangular coil: the cross-section integral of
	 * filamentPsi() over the coil, times its current density.
	 *
	 * @param order  Gauss points per direction per panel. Must be at least 2
	 *               and at most maximumCoilQuadratureOrder.
	 *
	 * Valid EVERYWHERE, including inside the coil -- see the file comment for
	 * what it achieves in each region, and note that the interior is where the
	 * weakly singular integrand lives and where the accuracy is 1e-12 rather
	 * than 1e-16.
	 *
	 * @throws std::invalid_argument on a non-finite point, a negative r, or an
	 *         order outside the accepted range.
	 */
	double coilPsi( Coil const &coil, double r, double z,
	                int order = defaultCoilQuadratureOrder,
	                double mu0 = vacuumPermeability );

	/**
	 * grad_bar( psi ) of one rectangular coil: the cross-section integral of
	 * filamentGradPsi() over the coil, times its current density.
	 *
	 * DIFFERENTIATED UNDER THE INTEGRAL SIGN, not differenced. The coil is a
	 * fixed region and the field point is the variable, so d/dr commutes with
	 * the integral over ( r', z' ) and this is the same panelled, cubically
	 * graded rule with the kernel replaced by its derivative. That matters: a
	 * difference of coilPsi() would cost two evaluations, floor at the
	 * difference's own O( h^2 ), and -- as CLAUDE.md records for every other
	 * differenced derivative in this tree -- need Richardson extrapolation to
	 * see past its own truncation.
	 *
	 * **THE INTERIOR IS WORSE THAN psi's AND THE EXTERIOR IS NOT, AND THAT IS
	 * MEASURED RATHER THAN ARGUED.** The derivative kernel is more singular
	 * than the kernel -- psi's integrand carries a logarithm where the field
	 * point meets the source and the gradient's carries a 1/distance -- so
	 * inside the coil the graded rule has a harder integrand to resolve, and
	 * outside it both integrands are analytic and both rules are spectral.
	 * Relative error against an order-128 reference, on the 0.30 x 0.40 coil
	 * at r = 2 with 1 MA:
	 *
	 *     probe                       n = 8    n = 16    n = 32    n = 48
	 *     inside,       grad       3.51e-03  9.12e-05  2.29e-06  2.28e-07
	 *                   psi        7.85e-07  1.59e-10  8.43e-13  1.24e-14
	 *     1 mm outside, grad       2.30e-04  2.64e-06  6.63e-10  2.78e-13
	 *                   psi        1.37e-06  1.63e-09  9.60e-14  1.04e-15
	 *     5 cm outside, grad       8.86e-06  3.70e-10  4.63e-15  2.96e-15
	 *                   psi        1.56e-07  2.48e-14  2.13e-15  4.82e-15
	 *
	 * **At the shipped default of 32 the exterior gradient is at round-off and
	 * the interior is at 2e-06** -- seven orders above what psi reaches on the
	 * same rule at the same point, and still falling only algebraically at
	 * n = 48 where psi has been flat since n = 32. Outside, the gradient is a
	 * constant behind psi rather than an order behind: it takes one more rung
	 * to reach round-off (n = 32 against n = 16 at 5 cm) and then stops, which
	 * is a spectral rule meeting the floor and not a rate.
	 *
	 * **FB-7 evaluates outside**, the conductor being beyond Gamma, so the
	 * exterior column is the one that governs -- but a caller reading the
	 * interior gradient should raise the order and check.
	 *
	 * @throws std::invalid_argument on the same conditions as coilPsi().
	 */
	void coilGradPsi( Coil const &coil, double r, double z,
	                  double &dPsiDr, double &dPsiDz,
	                  int order = defaultCoilQuadratureOrder,
	                  double mu0 = vacuumPermeability );

	/// The HDG flux q = ( 1/r ) grad_bar( psi ) of one rectangular coil.
	/// NaN on the axis in both components, for the reason filamentFlux() gives.
	void coilFlux( Coil const &coil, double r, double z,
	               double &qR, double &qZ,
	               int order = defaultCoilQuadratureOrder,
	               double mu0 = vacuumPermeability );

	/// The Gauss order ellipsePsi() uses unless told otherwise: nodes along
	/// each chord, with 4 times that many angles around it.
	///
	/// 16 rather than coilPsi()'s 32 because the two rules are not doing the
	/// same job. coilPsi() integrates a rectangle to 1e-12 INSIDE the
	/// conductor and is a reference field; ellipsePsi() exists for an initial
	/// GUESS, where the requirement is boundedness and smoothness rather than
	/// figures, and where the cost is paid once per quadrature point of a
	/// projection.
	inline constexpr int defaultEllipseQuadratureOrder = 16;

	/**
	 * The poloidal flux of a uniform current density over an ELLIPTICAL
	 * cross-section, coaxial with the axis: the cross-section integral of
	 * filamentPsi() over the ellipse, times its current density.
	 *
	 * @param centreR,centreZ  the centre of the ellipse, in metres.
	 * @param semiR,semiZ      its semi-axes. Both must be positive, and
	 *                         @a semiR must be strictly less than @a centreR.
	 * @param current          the SIGNED total current, which may be zero.
	 * @param order            Gauss nodes per chord. At least 2, at most
	 *                         maximumCoilQuadratureOrder.
	 *
	 * **WHY AN ELLIPSE AND NOT A RECTANGLE, WHICH meq::Coil ALREADY IS.** This
	 * is the shape a PLASMA COLUMN is, and the difference is not cosmetic: a
	 * uniform current density over a rectangle has four corners, and at each
	 * of them the second derivatives of psi carry a logarithm. Those are
	 * artefacts of the shape and nothing in the equilibrium puts them there.
	 * An ellipse has no corners and psi is C^infinity across its boundary in
	 * every direction but the normal, where it is C^1 -- which is what a
	 * current density with a jump gives and is all any guess needs.
	 *
	 * **AND NOT A FILAMENT, WHICH IS THE FAILURE THIS REPLACES.** psi of a
	 * filament diverges logarithmically AT the filament. Carrying I_p on one
	 * at the guessed magnetic axis therefore puts an unbounded spike exactly
	 * where the axis search has to look: measured on freegsnke's MAST-U,
	 * filamentPsi() reads 6.667e-01 at 3 mm from the guessed axis against a
	 * reference psi_axis of 9.187e-02, and it is still climbing. The same
	 * current spread over a finite cross-section stays bounded -- that is the
	 * whole of what a finite cross-section is for.
	 *
	 * **THE QUADRATURE IS POLAR ABOUT THE FIELD POINT, AND THAT IS WHAT MAKES
	 * ONE RULE SERVE INSIDE AND OUTSIDE.** The ellipse is swept as chords
	 * rho in [ rho-, rho+ ] along rays from the field point itself, so the
	 * AREA element rho d rho d theta carries a factor rho that meets the
	 * kernel's log( 1/rho ) head on; cubic grading toward rho- then leaves an
	 * integrand vanishing like t^5 log t, which Gauss takes to round-off. No
	 * node can land on the field point, because every node is at rho > 0.
	 * A field point inside the ellipse has rho- = 0 on every ray and is the
	 * case the grading is for; one outside has rho- > 0 and is graded toward
	 * its own nearest approach, which is where the kernel is largest.
	 *
	 * **THE ONE PLACE IT IS ONLY THREE FIGURES IS OUTSIDE, AND IT IS THE
	 * ANGLES RATHER THAN THE CHORDS.** Seen from an exterior point the ellipse
	 * subtends a cone, and at the two tangent rays the chord length vanishes
	 * like a square root -- so the midpoint rule in theta, which is spectral
	 * for an interior point, converges at about n^-1.5 for an exterior one.
	 * At the default order that is a few parts in 1e4 of the value, and rays
	 * that miss the ellipse cost nothing, so a far field point is both cheap
	 * and accurate. This is priced for a guess and is not a reference field;
	 * coilPsi() is the reference field.
	 *
	 * psi( 0, z ) is 0.0 BIT EXACTLY, as it is for both other conductor kinds:
	 * every chord's kernel carries k^2 = 4 a r/d^2 as a factor.
	 *
	 * @throws std::invalid_argument on a non-finite argument, a negative field
	 *         radius, a non-positive semi-axis, an ellipse reaching the axis,
	 *         or an order outside the accepted range.
	 */
	double ellipsePsi( double r, double z, double centreR, double centreZ,
	                   double semiR, double semiZ, double current,
	                   int order = defaultEllipseQuadratureOrder,
	                   double mu0 = vacuumPermeability );

	/**
	 * An IDEAL CIRCULAR FILAMENT: a ring current of zero cross-section,
	 * coaxial with the axis.
	 *
	 * Same conventions as meq::Coil in every respect -- psi = r A_phi in weber
	 * per radian, a SIGNED current in amperes that may be zero, geometry in
	 * metres, and no permeability of its own because mu0 belongs to whatever
	 * owns the conductor. It is deliberately shaped like Coil so that the two
	 * read alike at a call site.
	 *
	 * WHY IT EXISTS, AND IT IS NOT THAT A FILAMENT IS BETTER. meq::Coil models
	 * a real conductor: finite extent, uniform current density, a field that is
	 * finite everywhere including inside itself. A filament is an
	 * IDEALISATION, and FREE-BOUNDARY-PLAN.md section 7.19 records why MEQ
	 * wants one anyway -- **`../freegs4e`'s default `Coil` IS an exact
	 * filament**, `controlPsi` returning `Greens( R, Z, . )*turns` with its
	 * `area` used only for a current-density limit and never entering the
	 * field. So does FreeGS, and so do most codes of that family. MEQ could
	 * not model what the reference code models, and section 7.16's published
	 * 1.3e-04 agreement on `psi_ax` was reached ACROSS that mismatch rather
	 * than with it removed. This is what lets the two be matched, and the
	 * difference between two MEQ runs -- one filament, one rectangle, nothing
	 * else changed -- is then the finite-size effect measured on one code with
	 * one mesh and one solver.
	 *
	 * **IT STRUCTURALLY CANNOT DO WHAT Coil DOES, WHICH IS THE POINT OF HAVING
	 * BOTH.** psi diverges logarithmically at the filament and grad psi like
	 * 1/distance, so THERE IS NO SELF-FIELD AND NO SELF-FORCE -- the quantity a
	 * finite cross-section exists to make finite. A force calculation needs
	 * Coil whatever else is done. There is correspondingly no area(),
	 * currentDensity(), contains() or f(): a filament has infinite current
	 * density on a set of measure zero, so the Grad-Shafranov source term is
	 * not a function and CoilSet::f() has nothing to add. **A filament belongs
	 * OUTSIDE the computational domain**, which is exactly FB-7's
	 * configuration; one inside it would be evaluated at mesh points that may
	 * land on it, and the refusal in filamentPsi() is what that meets.
	 */
	class CurrentFilament
	{
		public:
			/// @param radiusIn   the ring's major radius, metres. Strictly
			///                   positive -- a filament at r <= 0 is the same
			///                   refusal meq::Coil makes for a coil reaching
			///                   the axis, and for the same reason: the
			///                   operator's 1/r is not integrable through it.
			/// @param heightIn   its height, metres.
			/// @param currentIn  the current, amperes. Signed; zero allowed,
			///                   because a scan over currents must be
			///                   writable.
			///
			/// @throws std::invalid_argument if any argument is not finite or
			///         if the radius is not positive.
			CurrentFilament( double radiusIn, double heightIn,
			                 double currentIn );

			double radius() const;
			double height() const;
			double current() const;

		private:
			double radiusValue;
			double heightValue;
			double currentValue;
	};

	/// psi of a filament. Identical to the five-argument filamentPsi(); this
	/// overload exists so that a caller holding a CurrentFilament need not
	/// unpack it.
	double filamentPsi( CurrentFilament const &filament, double r, double z,
	                    double mu0 = vacuumPermeability );

	/// grad_bar( psi ) of a filament. See the free function for the algebra.
	void filamentGradPsi( CurrentFilament const &filament, double r, double z,
	                      double &dPsiDr, double &dPsiDz,
	                      double mu0 = vacuumPermeability );

	/// q = ( 1/r ) grad_bar( psi ) of a filament. NaN on the axis.
	void filamentFlux( CurrentFilament const &filament, double r, double z,
	                   double &qR, double &qZ,
	                   double mu0 = vacuumPermeability );

	/**
	 * A set of coils, and the two things a free-boundary solve wants from them:
	 * the Grad-Shafranov source term, and the exact field.
	 *
	 * Coils may overlap -- nothing here prevents it and nothing here needs to,
	 * since f() sums the current densities of every coil containing the point,
	 * which is what a physical overlap would mean. indexContaining() is the one
	 * place the ambiguity shows, and it says so.
	 *
	 * mu0 is a constructor argument for the reason meq::Source gives: a run in
	 * normalised units sets it to 1, and it should be visible in one place
	 * rather than compiled in. It is meq::vacuumPermeability by default, which
	 * is the single 4 pi x 10^-7 in this tree.
	 */
	class CoilSet
	{
		public:
			/// @param mu0In  the permeability. Must be finite and positive; a
			///               zero would make every coil silently inert, which
			///               is worse than an error.
			explicit CoilSet( double mu0In = vacuumPermeability );

			/// Append a coil. The coil is copied; a Coil is five doubles.
			void add( Coil const &coil );

			std::size_t size() const;
			bool empty() const;

			/// @throws std::out_of_range naming the index and the size.
			Coil const &coil( std::size_t index ) const;

			/// The coils, in the order they were added.
			std::vector<Coil> const &coils() const;

			/**
			 * The Grad-Shafranov source term, IN MEQ'S F CONVENTION:
			 *
			 *     f( r, z ) = sum over coils containing ( r, z ) of
			 *                 mu0 * r * ( I / area ),
			 *
			 * and exactly zero outside every coil. The derivation of that
			 * factor is at the top of this file; the 1/r that turns F into the
			 * right hand side of the weak form belongs to the weak form and is
			 * NOT applied here, which is the same contract meq::Source::f()
			 * signs up to.
			 *
			 * Two arguments and not three: a coil current does not depend on
			 * psi. dF/dpsi is identically zero, so a solve driven by coils
			 * alone is affine.
			 *
			 * DISCONTINUOUS AT EVERY COIL EDGE, by construction -- see the file
			 * comment for what that costs a convergence rate and why FB-2's
			 * acceptance is a flux balance instead.
			 */
			double f( double r, double z ) const;

			/// The signed sum of the coil currents, in amperes.
			///
			/// This is the number FB-2 checks the outward flux against.
			/// Integrating Delta* psi = -F and using
			/// Delta* psi = r div_bar( ( 1/r ) grad_bar psi ) turns the area
			/// integral of j_phi into a boundary integral,
			///
			///     oint ( 1/r ) dpsi/dn dl = -mu0 * totalCurrent(),
			///
			/// so the total current is a property of the trace alone and the
			/// check is exact whatever the mesh does inside the coils.
			///
			/// THE SIGN IS NEGATIVE, and it is written out because the same
			/// identity as a circulation of B counterclockwise in ( r, z )
			/// comes out POSITIVE: phi-hat = z-hat x r-hat, so counterclockwise
			/// in the ( r, z ) plane has normal -phi-hat. Measured on one coil
			/// over an enclosing rectangle, by central differences of psi() and
			/// a midpoint rule refined from 200 to 400 points and Richardson
			/// extrapolated: -1.2566370749e+00 against -mu0 I =
			/// -1.2566370614e+00, 1.07e-08 relative. FREE-BOUNDARY-PLAN.md
			/// section 7 predicts this sign will be got wrong at least once;
			/// this is the version of it that is measured rather than argued.
			double totalCurrent() const;

			/// The index of a coil containing ( r, z ), or -1 if none does.
			///
			/// THE FIRST such coil in insertion order, which is only ambiguous
			/// if coils overlap -- and if they do, f() is the method that
			/// answers correctly, since it sums them. The return type is signed
			/// so that -1 can be the sentinel; a non-negative answer may be
			/// cast to std::size_t and passed to coil().
			int indexContaining( double r, double z ) const;

			/// The exact poloidal flux of the whole set: the sum over coils of
			/// coilPsi() at this set's quadrature order and permeability.
			///
			/// Delta* is linear, so a sum of coil fields is the field of the
			/// sum, and this is the closed form FB-1 compares a vacuum solve
			/// against. An empty set gives exactly 0.0.
			double psi( double r, double z ) const;

			/// One coil's contribution to psi(), for a caller separating them.
			/// @throws std::out_of_range as coil() does.
			double psiOf( std::size_t index, double r, double z ) const;

			/**
			 * grad_bar( psi ) of the whole set, and the HDG flux
			 * q = ( 1/r ) grad_bar( psi ).
			 *
			 * **q IS WHAT FB-7 NEEDS**: a conductor outside Gamma enters the
			 * coupling through the transmission condition, whose Neumann half
			 * is q . nu on Gamma. grad_bar psi is offered beside it because
			 * dividing by r is the caller's business on the axis -- see below
			 * -- and because a consumer computing B wants the pair undivided.
			 *
			 * Delta* is linear, so these sum over the coils exactly as psi()
			 * does, at the same quadrature order and permeability. An empty set
			 * gives exactly 0.0 in both components.
			 *
			 * **flux() IS NaN ON THE AXIS AND gradPsi() IS 0.0 THERE**, which
			 * is not the same statement twice. Every coil's grad_bar psi
			 * vanishes at r = 0 -- exactly, k^2 being an explicit factor -- and
			 * q divides that zero by zero. The limit is finite and this does
			 * not reach it. It matters because **a semicircle centred on the
			 * axis meets the axis at both ends**, so a sweep of such a Gamma
			 * hits the one place q is unavailable.
			 */
			void gradPsi( double r, double z,
			              double &dPsiDr, double &dPsiDz ) const;
			void flux( double r, double z, double &qR, double &qZ ) const;

			/// One coil's contribution to gradPsi().
			/// @throws std::out_of_range as coil() does.
			void gradPsiOf( std::size_t index, double r, double z,
			                double &dPsiDr, double &dPsiDz ) const;

			/// Gauss points per direction per panel for psi() and psiOf().
			/// @throws std::invalid_argument outside
			///         [ 2, maximumCoilQuadratureOrder ].
			void setQuadratureOrder( int order );
			int quadratureOrder() const;

			double mu0() const;

		private:
			std::vector<Coil> coilList;
			double permeability;
			int quadratureOrderValue;
	};

	/**
	 * CONDUCTORS OUTSIDE `Gamma`: the set FB-7 couples through the boundary
	 * rather than through the mesh. `FREE-BOUNDARY-PLAN.md` section 7.19.
	 *
	 * It holds meq::Coil and meq::CurrentFilament alike, sums their fields, and
	 * **HAS NO f()**. That absence is the whole reason it is a type of its own
	 * and not a flag on meq::CoilSet, and it is a modelling statement rather
	 * than a gap:
	 *
	 *   * **AN EXTERIOR CONDUCTOR CONTRIBUTES NOTHING TO THE INTERIOR
	 *     EQUATION.** Write psi = psi_coil + psi_tilde. The conductor is outside
	 *     Gamma, which contains dOmega, so Delta* psi_coil = 0 throughout Omega
	 *     and the source term is untouched. There is no interior current density
	 *     to return, not merely one nobody has written down.
	 *   * **AND A FILAMENT HAS NONE AT ALL**, anywhere: infinite current density
	 *     on a set of measure zero is not a function, so `CoilSet::f()` could not
	 *     be honest about one however the set were arranged. That is the design
	 *     question section 7.19 left open, and this is the answer to it -- the
	 *     two kinds of member can live in one set precisely because the one
	 *     method neither of them can answer is not on the interface.
	 *
	 * So the type expresses its own precondition. Handing an interior coil to
	 * GradShafranovSolver::setExteriorConductors() -- which silently drops its
	 * current from the source -- is now a compile error rather than a converged
	 * run describing a machine nobody asked for.
	 *
	 * **THERE IS DELIBERATELY NO CONVERSION FROM CoilSet.** The two sets answer
	 * opposite questions about the same conductor and the one error this type
	 * exists to prevent is a conductor appearing in BOTH -- once in the source
	 * and once through the coupling, which double-counts it while every residual
	 * and every border converges. Making that conversion one line would make the
	 * mistake one line. A caller who genuinely wants the same geometry both ways
	 * writes the coils out twice and can be asked why.
	 *
	 * **WHY IT CARRIES FILAMENTS AT ALL, AND IT IS NOT THAT THEY ARE BETTER.**
	 * See meq::CurrentFilament: `../freegs4e`'s default `Coil` IS an exact
	 * filament, so is FreeGS's, and section 7.16's published 1.3e-04 agreement
	 * in `psi_ax` was reached ACROSS that modelling difference rather than with
	 * it removed. A real coil has finite extent and finite current density, and
	 * near the plasma that matters; the filament is here so the mismatch can be
	 * MEASURED on one code with one mesh and one solver, not because it is the
	 * more accurate model.
	 *
	 * mu0 is a constructor argument for the reason meq::CoilSet gives, and the
	 * quadrature order is one for the reason coilPsi() gives -- **it governs the
	 * RECTANGLES only**. A filament's field is closed form and costs two Carlson
	 * evaluations whatever the order says.
	 */
	class ExteriorCoilSet
	{
		public:
			/// @param mu0In  the permeability. Must be finite and positive.
			explicit ExteriorCoilSet( double mu0In = vacuumPermeability );

			/// Append a conductor. Both are copied; each is a handful of
			/// doubles. The two overloads keep separate index spaces -- see
			/// coil() and filament() -- because a single index over a
			/// heterogeneous set would have to be decoded before it could be
			/// used for anything.
			void add( Coil const &coil );
			void add( CurrentFilament const &filament );

			/// Conductors of both kinds. size() is their sum, and is what
			/// empty() is about.
			std::size_t size() const;
			bool empty() const;
			std::size_t coilCount() const;
			std::size_t filamentCount() const;

			/// @throws std::out_of_range naming the index and the count.
			Coil const &coil( std::size_t index ) const;
			CurrentFilament const &filament( std::size_t index ) const;

			std::vector<Coil> const &coils() const;
			std::vector<CurrentFilament> const &filaments() const;

			/// The signed sum of every conductor's current, in amperes.
			///
			/// **IT IS NOT WHAT FB-2's FLUX BALANCE CHECKS.** That identity,
			/// `oint ( 1/r ) dpsi/dn dl = -mu0 I`, counts the current ENCLOSED,
			/// and by construction none of this is: an exterior conductor
			/// contributes exactly zero to a contour integral over Gamma. So
			/// this number is for reporting what a machine description was
			/// given, and for the acceptance that a rectangle and a filament
			/// being compared carry the same current.
			double totalCurrent() const;

			/// The poloidal flux of the whole set at ( r, z ): coilPsi() over
			/// the rectangles at this set's order, plus filamentPsi() over the
			/// filaments. Delta* is linear, so a sum of fields is the field of
			/// the sum. An empty set gives exactly 0.0.
			///
			/// **psi( 0, z ) IS 0.0 BIT EXACTLY** for both kinds of member, k^2
			/// being an explicit factor of the kernel -- which is the boundary
			/// condition the free-boundary problem imposes on the axis, and is
			/// worth having exactly rather than to round-off, since Gamma is a
			/// semicircle whose ends sit there.
			double psi( double r, double z ) const;

			/// grad_bar( psi ) of the whole set, and the HDG flux
			/// q = ( 1/r ) grad_bar( psi ). Summed exactly as psi() is.
			///
			/// **q IS WHAT THE TRANSMISSION ROW NEEDS** -- its Neumann half is
			/// q . nu on Gamma -- and **flux() IS NaN ON THE AXIS WHERE
			/// gradPsi() IS 0.0**, which is not the same statement twice: every
			/// member's grad_bar psi vanishes there exactly and q divides that
			/// zero by zero. The limit is finite and this does not reach it.
			/// GradShafranovSolver::conductorNormalFlux() carries the axis rule
			/// that deals with it.
			void gradPsi( double r, double z,
			              double &dPsiDr, double &dPsiDz ) const;
			void flux( double r, double z, double &qR, double &qZ ) const;

			/**
			 * How far the nearest conductor clears a semicircular `Gamma` of
			 * radius @a rhoGamma centred at ( 0, @a centreZ ), in metres.
			 * Positive when EVERY member lies strictly outside it, which is
			 * this type's precondition; an empty set gives infinity.
			 *
			 * **THIS IS THE CHECK GradShafranovSolver COULD NOT MAKE.** Its
			 * setExteriorConductors() documents that being outside Gamma "IS NOT
			 * CHECKED HERE", because at that point the solver does not know
			 * where Gamma is -- the DtN arrives separately. The geometry is
			 * knowable the moment both are in hand, and this is the half of it
			 * that needs no MFEM.
			 *
			 * A conductor inside Gamma is not a small error: psi_coil is then
			 * not Delta*-harmonic where the expansion assumes it is, so the
			 * Gegenbauer series does not represent the field it is being asked
			 * to represent, and the run converges to a machine nobody described.
			 *
			 * The distance is to the nearest POINT of each member -- the closest
			 * corner or edge of a rectangle, the ring itself for a filament --
			 * so a coil straddling Gamma reports a negative clearance rather
			 * than being judged by its centre.
			 */
			double clearance( double centreZ, double rhoGamma ) const;

			/// Gauss points per direction per panel for the RECTANGLES.
			/// @throws std::invalid_argument outside
			///         [ 2, maximumCoilQuadratureOrder ].
			void setQuadratureOrder( int order );
			int quadratureOrder() const;

			double mu0() const;

		private:
			std::vector<Coil> coilList;
			std::vector<CurrentFilament> filamentList;
			double permeability;
			int quadratureOrderValue;
	};



	/**
	 * A meq::Source that is a plasma source PLUS a coil set.
	 *
	 * This is the adapter Coil's own comment says belongs with the assembly
	 * rather than with the geometry, and it is here rather than in Source.hpp
	 * because it is the coils that need adapting: `Source::f()` takes psi and a
	 * coil current does not depend on it.
	 *
	 * `F = F_plasma( r, z, psi ) + F_coil( r, z )`, and
	 * `dF/dpsi = dF_plasma/dpsi` exactly -- the coils contribute nothing to the
	 * Jacobian, so a solve driven by coils ALONE is affine and Newton finishes
	 * in one step. That is the property `theCoilsMakeAVacuumSolveAffine`
	 * asserts, and it is the cheapest available statement that the coil term
	 * really is data.
	 *
	 * **SHARED OWNERSHIP OF BOTH, BECAUSE THE DRIVER NEEDS IT.** `apps/meq.cpp`
	 * builds a fresh solver every adaptive cycle and hands it the same source;
	 * the source must outlive every one of them, and the coil set is separately
	 * useful to whatever computes the exterior field. Neither may be null.
	 */
	class CoilAugmentedSource : public Source
	{
		public:
			/// @throws std::invalid_argument if either argument is null.
			CoilAugmentedSource( std::shared_ptr<Source const> plasmaIn,
			                     std::shared_ptr<CoilSet const> coilsIn );

			/// F_plasma( r, z, psi ) + F_coil( r, z ).
			double f( double r, double z, double psi ) const override;

			/// dF_plasma/dpsi. The coils contribute exactly zero.
			double dFdPsi( double r, double z, double psi ) const override;

			Source const & plasma() const;
			CoilSet const & coils() const;

			/// The set above. See meq::Source::conductors for what asks.
			CoilSet const * conductors() const override;

		private:
			std::shared_ptr<Source const> plasmaSource;
			std::shared_ptr<CoilSet const> coilSet;
	};

	/**
	 * The same sum, when the plasma source's profiles are functions of
	 * NORMALISED flux and psi_ax is therefore an unknown of the non-linear
	 * system.
	 *
	 * It is a separate class rather than a flag for the reason
	 * meq::NormalisedSource is a separate class from meq::MHDSource: the solver
	 * takes a `NormalisedSource &` on one path and a `Source const &` on the
	 * other, and which one it is decides whether there is a border row at all.
	 *
	 * **EVERY NORMALISATION CALL IS FORWARDED TO THE PLASMA SOURCE AND NONE IS
	 * ANSWERED HERE.** The normalisation belongs to the profiles, the profiles
	 * belong to the wrapped source, and duplicating the pair of doubles here
	 * would create two answers to `normalisation()` that could disagree. The
	 * coils are not normalised by anything -- their current is amperes.
	 *
	 * **AND setPlasmaSupport() IS FORWARDED TOO, WHICH IS WHY THE BASE MADE IT
	 * VIRTUAL.** `insidePlasma()` is consulted by whichever object evaluates
	 * the profiles, and that is the wrapped one. See
	 * NormalisedSource::setPlasmaSupport.
	 *
	 * **THE COIL TERM IS OUTSIDE THE PLASMA SUPPORT, DELIBERATELY.** With the
	 * support on, `F_plasma` vanishes wherever `Psi <= 0` and `F_coil` does
	 * not: a coil sits in the vacuum region by construction, so confining it to
	 * the plasma would switch off every coil in the machine. The sum is
	 * therefore taken after the plasma term has been confined and never before.
	 */
	class CoilAugmentedNormalisedSource : public NormalisedSource
	{
		public:
			/// @throws std::invalid_argument if either argument is null.
			///
			/// The plasma source is NON-const, as meq::makeNormalisedSource
			/// returns it and for the same reason: the solver calls
			/// setNormalisation() on it before every residual evaluation.
			CoilAugmentedNormalisedSource(
				std::shared_ptr<NormalisedSource> plasmaIn,
				std::shared_ptr<CoilSet const> coilsIn );

			double f( double r, double z, double psi ) const override;
			double dFdPsi( double r, double z, double psi ) const override;

			void setNormalisation( double psiAxis,
			                       double psiBoundary ) override;
			using NormalisedSource::setNormalisation;

			double normalisation() const override;
			double boundaryNormalisation() const override;

			void setPlasmaSupport( bool confined ) override;

			/// FORWARDED, AND FOR THE REASON setPlasmaSupport() IS. The plasma
			/// term is evaluated through the wrapped source, so it is that
			/// source's frozen edge insidePlasma() consults; freezing only this
			/// one would leave the support moving and say it was not.
			/// @see meq::NormalisedSource::freezePlasmaEdge
			void freezePlasmaEdge( double axis, double boundary ) override;

			/// @see freezePlasmaEdge
			void thawPlasmaEdge() override;

			/// FORWARDED, AND THIS IS THE FOURTH TIME. XP-1's connectivity test
			/// switches the PLASMA term off on any element the flood fill from
			/// the axis does not reach, and meq::SourceIntegrator asks the
			/// source what is left there. The base returns zero, which for this
			/// class would switch off every conductor in the machine wherever
			/// the fill did not reach -- and a coil sits in the vacuum region by
			/// construction, so that is everywhere it matters. The answer is the
			/// coil term, exactly, with no cancellation: computed as
			/// f() - scaledF() it would be ( plasma + coil ) - plasma, which is
			/// the coil term only to round-off.
			double fOutsidePlasma( double r, double z ) const override;

			/// FORWARDED, AND FOR THE SAME REASON setPlasmaSupport() IS. The
			/// normalisations belong to the profiles, which belong to the
			/// wrapped source, so the derivatives do too -- and the COIL term
			/// contributes exactly nothing to them, a coil current being
			/// amperes and not a function of any flux. Without this override the
			/// base's default returns false and every run carrying a [[coils]]
			/// block would silently fall back to a DIFFERENCED border column,
			/// which is the trap the base class documents one method up.
			bool normalisationDerivatives( double r, double z, double psi,
			                               double &dFdAxis,
			                               double &dFdBoundary ) const override;

			/// The set this wraps. See meq::Source::conductors for what asks.
			/// NOT forwarded: the conductors are THIS class's, not the plasma
			/// source's, which is the one thing here that does not delegate.
			CoilSet const * conductors() const override;

			/// FORWARDED, for the third time and for the same reason: the scale
			/// belongs to the plasma term and a COIL IS NOT SCALED. Its current
			/// is amperes and is prescribed input, so scaling it would make the
			/// prescribed plasma current move the conductors too.
			void setCurrentScale( double scale ) override;
			double currentScale() const override;

			/// The plasma term alone, which is what the current constraint is
			/// about -- this class's own f() is the SUM and would put the coil
			/// current inside the prescribed plasma current.
			double scaledF( double r, double z, double psi ) const override;
			double scaledDFdPsi( double r, double z, double psi ) const override;

			NormalisedSource & plasma() const;
			CoilSet const & coils() const;

		private:
			std::shared_ptr<NormalisedSource> plasmaSource;
			std::shared_ptr<CoilSet const> coilSet;
	};

}

#endif // MEQ_COILS_HPP
