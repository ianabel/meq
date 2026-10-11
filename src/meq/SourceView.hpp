#ifndef MEQ_SOURCEVIEW_HPP
#define MEQ_SOURCEVIEW_HPP

#include "Source.hpp"

#include "mfem.hpp"

#include <vector>

/*
 * `F( R, z, psi )` AND `dF/dpsi` AS PLAIN DATA, SO A DEVICE KERNEL CAN EVALUATE
 * THEM.
 *
 * WHAT THIS IS FOR. `meq::SourceIntegrator` is the one integrator MEQ owns and
 * the only thing between MEQ and upstream's batched local residual: everything
 * else in GradShafranovSolver's configuration already passes --
 * `CanBatchFluxMass`, `CanBatchDiv`, `CanBatchPotFaceAssembly`,
 * `CanBatchPotBdrFaceAssembly`, `CanBatchLocalFactor`, `CanBatchLocalSolve` and
 * `CanBatchTraceAssembly` are all YES on MEQ's own forms, and
 * `CanBatchLinearResidual()` is the single NO. It refuses on
 *
 *     if (m_nlfi_p && !HDGIntegratorIsLinear(m_nlfi_p)) { return false; }
 *
 * which is a correct refusal and not an oversight: meq::SourceIntegrator IS
 * nonlinear. What is missing is not a kernel but an INTERFACE --
 * `AssembleElementVector` and `AssembleElementGrad` are host virtuals over
 * dense element data, so there is nothing for a device to dispatch on. This
 * header is MEQ's half of that interface: the POD a caller can hand across, and
 * the two functions that read it.
 *
 * THIS IS A SECOND EXPRESSION OF THE PHYSICS AND A TEST IS WHAT MAKES THAT
 * SAFE. `CLAUDE.md` warns against exactly this -- "a second expression of
 * F( R, z, psi )" as "two implementations that must agree" -- and the warning
 * is right: the day one is edited and the other is not, MEQ solves a different
 * equilibrium and converges to it at the full rate, which no convergence table
 * can see. meq::PointLocator answers the same hazard by having ONE expression
 * that two compilers read, and that was available there because the arithmetic
 * had no state.
 *
 * IT IS NOT AVAILABLE HERE, AND THE REASON IS WHERE THE STATE LIVES.
 * meq::NormalisedMHDSource::f() reads psi_ax, psi_bnd, the current scale and
 * the frozen support, and three of those four are unknowns of the bordered
 * Newton that move between residual evaluations. Routing the virtual through a
 * SourceView would mean rebuilding that struct per call -- on MEQ's hottest
 * path, the residual leg at 48.2% serial and 23.1% threaded (M-80) -- or
 * caching it and inheriting a staleness bug whose symptom is a Jacobian taken
 * at last iteration's normalisation. **A cache that is wrong for one Newton
 * step is worse than a duplicate that a test pins**, because the test fails
 * loudly and the cache does not.
 *
 * So the two expressions stand and
 * `tests/unit/SourceViewTests.cpp` asserts they agree to the bit over a sweep
 * that includes the knots themselves, both sides of the plasma edge, and
 * outside the profile support -- the three places a spline and its POD form can
 * disagree. **That test is load-bearing: if it is ever deleted or relaxed, this
 * header becomes a second definition of the physics with nothing holding it to
 * the first.**
 *
 * WHAT IT DELIBERATELY DOES NOT COVER. `makeSourceView()` REFUSES rather than
 * approximates: a source or profile shape with no POD form returns false and
 * the caller keeps the host path. A refusal costs a device kernel; a POD that
 * is nearly right costs an equilibrium, and nothing downstream would say so.
 * Not covered today, and each would be a new `kind` rather than a change of
 * shape: the rotating source (`meq::RotatingSource`, whose closure is a root
 * find), the coil-augmented wrappers, and any profile read through
 * `meq::ScaledProfile` around something that is not itself representable.
 *
 * THE PLASMA SUPPORT IS TWO TESTS AND ONLY ONE OF THEM IS HERE. The pointwise
 * `insidePlasma( R, z, psi )` is arithmetic on the point and is in `sourceF()`.
 * XP-1's CONNECTIVITY test is over the element adjacency graph, is decided by
 * the solver, and reaches the integrator as a per-element flag -- so it stays
 * outside this header, exactly as it stays outside meq::Source.
 */

namespace meq
{

	/// One Hermite knot, as plain data. Trivial on purpose: mfem::Array<T>
	/// static_asserts `type T must be trivial` (general/array.hpp:61), so a
	/// default member initialiser here is a compile error in a device
	/// translation unit and nowhere else.
	struct KnotView
	{
		mfem::real_t psi;
		mfem::real_t value;
		mfem::real_t derivative;
	};

	/// What kind of profile a ProfileView describes.
	enum class ProfileKind : int
	{
		Constant = 0,   ///< meq::ConstantProfile
		Spline   = 1    ///< meq::SplineProfile, Hermite cubic over KnotView[]
	};

	/**
	 * A meq::Profile as plain data.
	 *
	 * `scale` folds meq::ScaledProfile in rather than giving it a kind of its
	 * own -- a scaled profile is its inner one times a constant, and carrying
	 * the factor here keeps the evaluator branch-free where it matters.
	 * `offset` and `count` index a flat KnotView array the caller owns, so one
	 * allocation carries every profile in a source.
	 */
	struct ProfileView
	{
		int kind;                  ///< a ProfileKind
		mfem::real_t scale;        ///< 1 unless a ScaledProfile wrapped it
		mfem::real_t constantValue;///< Constant only
		int offset;                ///< Spline only: first knot
		int count;                 ///< Spline only: number of knots

		/*
		 * THE UNIFORM FAST PATH, AND IT IS THE DIFFERENCE BETWEEN A BRANCHY
		 * BINARY SEARCH AND ONE MULTIPLY IN THE INNERMOST LOOP MEQ HAS.
		 *
		 * Every profile MEQ ships is uniform -- `diverted-tokamak-pprime.dat`
		 * is 40 knots at exactly 0.01 -- so the `log2 n` dependent loads and
		 * unpredictable branches of profileInterval()'s search buy nothing on
		 * the production path. On a host that is a mispredict per quadrature
		 * point; **on a device it is warp divergence in the innermost loop of
		 * the source kernel**, where neighbouring points in a warp sit in
		 * different intervals and the whole warp walks every branch.
		 *
		 * `uniform` is decided when the view is built, by checking the stored
		 * abscissae rather than by trusting the file, and a non-uniform profile
		 * keeps the search and is merely slower. Nothing about WHICH interval
		 * is chosen changes: see profileInterval().
		 */
		int uniform;               ///< Spline only: whether the fast index applies
		mfem::real_t firstKnot;    ///< Spline + uniform: knots[ offset ].psi
		mfem::real_t inverseSpacing;///< Spline + uniform: 1/h
	};

	/// What kind of source a SourceView describes.
	enum class SourceKind : int
	{
		Soloviev          = 0,  ///< F = -( ( 1 - a ) R^2 + a ); dF/dpsi == 0
		NormalisedMHD     = 1   ///< meq::NormalisedMHDSource
	};

	/**
	 * A meq::Source as plain data, sufficient to evaluate F and dF/dpsi.
	 *
	 * Everything is by value except the knots, which are indexed. There are no
	 * pointers in here at all, deliberately: a struct holding a device pointer
	 * is a struct that has to be rebuilt when the allocation moves, and M-79
	 * and M-130 are both MEQ holding a view onto something whose device copy
	 * went out from under it.
	 */
	struct SourceView
	{
		int kind;                      ///< a SourceKind
		mfem::real_t permeability;     ///< mu_0
		mfem::real_t psiAxis;          ///< psi_ax, an unknown of the bordered Newton
		mfem::real_t psiBoundary;      ///< psi_bnd, likewise
		mfem::real_t currentScale;     ///< the profile scale the I_p border sets
		mfem::real_t solovievA;        ///< Soloviev only: F = -( ( 1 - a ) R^2 + a )
		/*
		 * THE SUPPORT BOUNDS ARE SEPARATE FROM psiAxis/psiBoundary AND MUST BE.
		 * meq::Source::insidePlasma() tests against supportAxis() and
		 * supportBoundary(), which XP-2 FREEZES for the duration of a solve so
		 * that the support stops moving while psi_ax and psi_bnd remain
		 * unknowns of the same Newton. Reusing psiAxis here would make the
		 * support move with every border update -- which is precisely the
		 * oscillation the freeze exists to stop, and it would do it silently.
		 */
		mfem::real_t supportAxis;
		mfem::real_t supportBoundary;
		int confineToPlasma;           ///< whether insidePlasma() gates the term
		/*
		 * THE SUPPORT'S SHAPE, AND IT IS CARRIED RATHER THAN REFUSED.
		 * meq::PlasmaSupportShape::Box is a mirror's -- `0 < Psi <= 1` and
		 * `|z| <= axialEdge` -- and every caller of insidePlasmaValue() already
		 * holds `z`, so honouring it is two scalars and two comparisons. A
		 * refusal would have been one line, and it would also have left
		 * DeviceNPCDevice.cpp's refresh(), which keeps the last view it COULD
		 * build, integrating a one-sided support after the source had been
		 * boxed. Zero-initialised by makeSourceView(), which is OneSided.
		 */
		int supportShape;              ///< a meq::PlasmaSupportShape
		mfem::real_t axialEdge;        ///< Box only: the |z| bound
		ProfileView pPrime;            ///< p'( Psi )
		ProfileView ggPrime;           ///< gg'( Psi )
	};

	/**
	 * One rectangular conductor as plain data.
	 *
	 * **THE COIL TERM IS POD-ABLE AND THAT IS NOT OBVIOUS FROM THE CLASS.**
	 * meq::CoilSet also carries `psi()`, `flux()` and `psiOf()`, which are
	 * elliptic integrals and are emphatically not device arithmetic -- but the
	 * SOURCE term is not any of those. `CoilSet::f()` is
	 * `mu0 * R * sum over coils containing ( R, z ) of currentDensity`, a
	 * closed rectangle test and a constant, and that is all the residual
	 * integrand needs.
	 *
	 * The bounds are the closed ones `meq::Coil::contains()` tests, edges
	 * included, which its own comment says is stated rather than merely chosen.
	 */
	struct CoilView
	{
		mfem::real_t minRadius;
		mfem::real_t maxRadius;
		mfem::real_t minZ;
		mfem::real_t maxZ;
		mfem::real_t currentDensity;
	};

	/**
	 * `F_coil( R, z )`, the conductor half of a coil-augmented source.
	 *
	 * **SUMMED RATHER THAN SHORT-CIRCUITED ON THE FIRST HIT**, because
	 * overlapping coils really do add their current densities -- which is
	 * `CoilSet::f()`'s own comment, and the one place a "find the coil"
	 * rewrite would silently change an answer.
	 *
	 * Exactly zero with @a count zero, which is what keeps a source with no
	 * conductors bit-identical through this path.
	 */
	MFEM_HOST_DEVICE inline mfem::real_t coilF( CoilView const *coils, int count,
	                                            mfem::real_t permeability,
	                                            mfem::real_t radius,
	                                            mfem::real_t z )
	{
		mfem::real_t density = mfem::real_t( 0 );
		for ( int i = 0; i < count; ++i )
			if ( radius >= coils[ i ].minRadius && radius <= coils[ i ].maxRadius
			     && z >= coils[ i ].minZ && z <= coils[ i ].maxZ )
				density += coils[ i ].currentDensity;

		return permeability*radius*density;
	}

	/// The four Hermite basis functions on [ 0, 1 ], and their derivatives.
	/// @{
	MFEM_HOST_DEVICE inline mfem::real_t hermite00( mfem::real_t t )
	{ return ( mfem::real_t( 1 ) + mfem::real_t( 2 )*t )*( mfem::real_t( 1 ) - t )*( mfem::real_t( 1 ) - t ); }
	MFEM_HOST_DEVICE inline mfem::real_t hermite10( mfem::real_t t )
	{ return t*( mfem::real_t( 1 ) - t )*( mfem::real_t( 1 ) - t ); }
	MFEM_HOST_DEVICE inline mfem::real_t hermite01( mfem::real_t t )
	{ return t*t*( mfem::real_t( 3 ) - mfem::real_t( 2 )*t ); }
	MFEM_HOST_DEVICE inline mfem::real_t hermite11( mfem::real_t t )
	{ return t*t*( t - mfem::real_t( 1 ) ); }

	MFEM_HOST_DEVICE inline mfem::real_t hermite00Prime( mfem::real_t t )
	{ return mfem::real_t( 6 )*t*( t - mfem::real_t( 1 ) ); }
	/// EXPANDED, AS `Profiles.cpp`'s `h10Prime` WRITES IT. The factored form
	/// `( 3t - 1 )( t - 1 )` is the same polynomial and a different double:
	/// measured as a 1-ulp disagreement in every derivative-carrying device
	/// route against its host loop, while every value-only route agreed to the
	/// bit. ResidentIterate.hpp.
	MFEM_HOST_DEVICE inline mfem::real_t hermite10Prime( mfem::real_t t )
	{ return mfem::real_t( 3 )*t*t - mfem::real_t( 4 )*t + mfem::real_t( 1 ); }
	MFEM_HOST_DEVICE inline mfem::real_t hermite01Prime( mfem::real_t t )
	{ return mfem::real_t( 6 )*t*( mfem::real_t( 1 ) - t ); }
	MFEM_HOST_DEVICE inline mfem::real_t hermite11Prime( mfem::real_t t )
	{ return t*( mfem::real_t( 3 )*t - mfem::real_t( 2 ) ); }
	/// @}

	/**
	 * The interval of @a view containing @a x, by binary search.
	 *
	 * THE SAME ANSWER std::upper_bound GIVES, which is what the host spline
	 * uses -- and it has to be the same at a KNOT, not merely nearby. A profile
	 * with a flat segment meeting a rising one has a genuine kink there, so
	 * choosing the interval on the wrong side of an abscissa that lands exactly
	 * on a knot changes `prime` discontinuously. meq::Profile::doublePrime's
	 * own documentation makes the same point about adjacent intervals
	 * disagreeing at a shared knot.
	 */
	MFEM_HOST_DEVICE inline int profileInterval( ProfileView const &view,
	                                             KnotView const *knots,
	                                             int count, mfem::real_t x )
	{
		int interval;

		if ( view.uniform )
		{
			/*
			 * THE ARITHMETIC IS NOT TRUSTED, IT IS CORRECTED, and that is what
			 * makes this bit-identical rather than nearly so. `firstKnot` and
			 * the spacing are stored doubles read from a file, so they are
			 * uniform to within ulps and not exactly -- and
			 * `( x - x0 )*( 1/h )` can therefore land one interval either side
			 * of the answer `std::upper_bound` would give near a knot. Two
			 * comparisons put it back, which is still O(1) and branch-light,
			 * where being merely close would be a silently different cubic on
			 * a set of points nobody enumerates.
			 *
			 * At an EXACT knot either neighbour would answer identically --
			 * a Hermite interpolant matches value and derivative there, and
			 * h01(1) = h00(0) = 1 with the rest vanishing, so both sides
			 * return the knot datum to the bit. The correction is for the
			 * points NEAR one, where they would not.
			 */
			interval = static_cast<int>( ( x - view.firstKnot )*view.inverseSpacing );
			if ( interval < 0 )
				interval = 0;
			if ( interval > count - 2 )
				interval = count - 2;
			if ( interval > 0 && x < knots[ interval ].psi )
				--interval;
			else if ( interval < count - 2 && x >= knots[ interval + 1 ].psi )
				++interval;
			return interval;
		}

		int lo = 0;
		int hi = count;                       // first index with knots[ i ].psi > x
		while ( lo < hi )
		{
			int const mid = lo + ( hi - lo )/2;
			if ( knots[ mid ].psi > x )
				hi = mid;
			else
				lo = mid + 1;
		}
		interval = lo - 1;
		if ( interval < 0 )
			interval = 0;
		if ( interval > count - 2 )
			interval = count - 2;
		return interval;
	}

	/// @a view evaluated at @a x. @a knots is the flat array @a view indexes.
	MFEM_HOST_DEVICE inline mfem::real_t evaluateProfile( ProfileView const &view,
	                                                      KnotView const *knots,
	                                                      mfem::real_t x )
	{
		if ( view.kind == static_cast<int>( ProfileKind::Constant ) )
			return view.scale*view.constantValue;

		KnotView const *k = knots + view.offset;
		// Constant extension outside, matching SplineProfile::operator().
		if ( x <= k[ 0 ].psi )
			return view.scale*k[ 0 ].value;
		if ( x >= k[ view.count - 1 ].psi )
			return view.scale*k[ view.count - 1 ].value;

		int const i = profileInterval( view, k, view.count, x );
		mfem::real_t const lo = k[ i ].psi;
		mfem::real_t const hi = k[ i + 1 ].psi;
		mfem::real_t const delta = hi - lo;
		mfem::real_t const t = ( x - lo )/delta;

		return view.scale*( k[ i ].value*hermite00( t )
		                    + k[ i ].derivative*delta*hermite10( t )
		                    + k[ i + 1 ].value*hermite01( t )
		                    + k[ i + 1 ].derivative*delta*hermite11( t ) );
	}

	/// The derivative of @a view at @a x, with respect to @a x.
	MFEM_HOST_DEVICE inline mfem::real_t evaluateProfilePrime( ProfileView const &view,
	                                                           KnotView const *knots,
	                                                           mfem::real_t x )
	{
		if ( view.kind == static_cast<int>( ProfileKind::Constant ) )
			return mfem::real_t( 0 );

		KnotView const *k = knots + view.offset;
		// Zero outside, consistent with the constant extension above --
		// HermiteCubicSpline::prime() makes the same choice for the same reason.
		if ( x < k[ 0 ].psi || x > k[ view.count - 1 ].psi )
			return mfem::real_t( 0 );

		int const i = profileInterval( view, k, view.count, x );
		mfem::real_t const lo = k[ i ].psi;
		mfem::real_t const hi = k[ i + 1 ].psi;
		mfem::real_t const delta = hi - lo;
		mfem::real_t const t = ( x - lo )/delta;

		// The scale AFTER the division: meq::ScaledProfile::prime() is
		// `scale*inner.prime()`, and the inner derivative already divided.
		return view.scale*( ( k[ i ].value*hermite00Prime( t )
		                      + k[ i ].derivative*delta*hermite10Prime( t )
		                      + k[ i + 1 ].value*hermite01Prime( t )
		                      + k[ i + 1 ].derivative*delta*hermite11Prime( t ) )/delta );
	}

	/// Whether ( @a z, @a psi ) is inside the plasma, pointwise -- the value
	/// test only. meq::NormalisedSource::insidePlasma() written once more, in
	/// its order: normalised flux is positive inside by construction, whichever
	/// sign the span has, which is why the one-sided test is a product and not
	/// a comparison on psi, and the box's upper bound is on Psi for the same
	/// reason. `SourceViewTests` holds the two to the bit across both shapes.
	MFEM_HOST_DEVICE inline bool insidePlasmaValue( SourceView const &view,
	                                                mfem::real_t z,
	                                                mfem::real_t psi )
	{
		if ( !view.confineToPlasma )
			return true;
		mfem::real_t const span = view.supportAxis - view.supportBoundary;
		if ( view.supportShape == static_cast<int>( PlasmaSupportShape::OneSided ) )
			return ( psi - view.supportBoundary )*span > mfem::real_t( 0 );

		// PlasmaSupportShape::Box.
		mfem::real_t const height = z < mfem::real_t( 0 ) ? -z : z;
		if ( !( height <= view.axialEdge ) )
			return false;
		mfem::real_t const normalised = ( psi - view.supportBoundary )/span;
		return normalised > mfem::real_t( 0 ) && normalised <= mfem::real_t( 1 );
	}

	/**
	 * `F( R, z, psi )`.
	 *
	 * THE ONE EXPRESSION OF THE PHYSICS. meq::NormalisedMHDSource::f() calls
	 * this; it does not reimplement it.
	 */
	MFEM_HOST_DEVICE inline mfem::real_t sourceF( SourceView const &view,
	                                              KnotView const *knots,
	                                              mfem::real_t radius,
	                                              mfem::real_t z,
	                                              mfem::real_t psi )
	{
		if ( view.kind == static_cast<int>( SourceKind::Soloviev ) )
			return -( ( mfem::real_t( 1 ) - view.solovievA )*radius*radius
			          + view.solovievA );

		if ( !insidePlasmaValue( view, z, psi ) )
			return mfem::real_t( 0 );

		mfem::real_t const span = view.psiAxis - view.psiBoundary;
		mfem::real_t const psiN = ( psi - view.psiBoundary )/span;

		return view.currentScale
		       *( view.permeability*radius*radius
		            *evaluateProfile( view.pPrime, knots, psiN )
		          + evaluateProfile( view.ggPrime, knots, psiN ) )
		       /span;
	}

	/**
	 * `dF/dpsi` at the same point.
	 *
	 * ONE FACTOR OF 1/span FOR THE PROFILE ARGUMENT AND ONE ALREADY IN F, so
	 * this carries `1/span^2` where F carries `1/span`. The profiles are
	 * differentiated with respect to NORMALISED flux and psi is not, which is
	 * exactly the chain rule the host expression writes out; getting it wrong
	 * is a Jacobian that is merely scaled, and *A wrong Jacobian is invisible
	 * to a convergence table*.
	 */
	MFEM_HOST_DEVICE inline mfem::real_t sourceDFdPsi( SourceView const &view,
	                                                   KnotView const *knots,
	                                                   mfem::real_t radius,
	                                                   mfem::real_t z,
	                                                   mfem::real_t psi )
	{
		if ( view.kind == static_cast<int>( SourceKind::Soloviev ) )
			return mfem::real_t( 0 );

		if ( !insidePlasmaValue( view, z, psi ) )
			return mfem::real_t( 0 );

		mfem::real_t const span = view.psiAxis - view.psiBoundary;
		mfem::real_t const psiN = ( psi - view.psiBoundary )/span;

		return view.currentScale
		       *( view.permeability*radius*radius
		            *evaluateProfilePrime( view.pPrime, knots, psiN )
		          + evaluateProfilePrime( view.ggPrime, knots, psiN ) )
		       /( span*span );
	}

	/**
	 * `dF/dpsi_ax` and `dF/dpsi_bnd` at the same point, which is
	 * `meq::NormalisedMHDSource::normalisationDerivatives()` written once more
	 * as arithmetic a kernel can run.
	 *
	 * `F = scale g( Psi )/span` with `Psi = ( psi - psi_bnd )/span`, so
	 * `dPsi/dpsi_ax = -Psi/span` and `dPsi/dpsi_bnd = ( Psi - 1 )/span`, and the
	 * span moves by `+1` and `-1`. **Every operation is in the host's order**,
	 * because the plasma-current border's corner and the two normalisation
	 * columns are what the device route computes with it and the acceptance
	 * there is the host's digits rather than its neighbourhood.
	 *
	 * Solov'ev has no normalisation and is answered with zeros; the host has no
	 * NormalisedSource to ask in that case, so a caller never reaches it.
	 */
	MFEM_HOST_DEVICE inline void sourceNormalisationDerivatives(
		SourceView const &view, KnotView const *knots,
		mfem::real_t radius, mfem::real_t z, mfem::real_t psi,
		mfem::real_t &dFdAxis, mfem::real_t &dFdBoundary )
	{
		dFdAxis = mfem::real_t( 0 );
		dFdBoundary = mfem::real_t( 0 );
		if ( view.kind == static_cast<int>( SourceKind::Soloviev ) )
			return;
		if ( !insidePlasmaValue( view, z, psi ) )
			return;

		mfem::real_t const span = view.psiAxis - view.psiBoundary;
		mfem::real_t const psiN = ( psi - view.psiBoundary )/span;

		mfem::real_t const g = view.permeability*radius*radius
		                         *evaluateProfile( view.pPrime, knots, psiN )
		                       + evaluateProfile( view.ggPrime, knots, psiN );
		mfem::real_t const gPrime = view.permeability*radius*radius
		                              *evaluateProfilePrime( view.pPrime, knots, psiN )
		                            + evaluateProfilePrime( view.ggPrime, knots, psiN );

		dFdAxis = -view.currentScale*( gPrime*psiN + g )/( span*span );
		dFdBoundary = view.currentScale*( gPrime*( psiN - mfem::real_t( 1 ) ) + g )
		              /( span*span );
	}

	/**
	 * The element-local residual contribution of the source term, at one
	 * quadrature point.
	 *
	 * **THIS IS THE SIGNATURE MEQ IS ASKING UPSTREAM TO DISPATCH ON**, and it
	 * is deliberately the smallest one that works: everything it needs is a
	 * scalar or a stride into arrays the caller already owns, there is no MFEM
	 * object in it, and it is `MFEM_HOST_DEVICE`. A batched local residual can
	 * call it per point per element with no host round trip.
	 *
	 * It computes what `meq::SourceIntegrator::AssembleElementVector` computes
	 * inside its quadrature loop:
	 *
	 *     elvect -= weight * F( R, z, psi + psi_c ) / R * shape
	 *
	 * with `psi = shape . elfun` supplied by the caller rather than formed
	 * here, because the caller already has the element's dof values laid out
	 * the way its kernel wants them and MEQ has no opinion about that.
	 *
	 * **THE 1/R IS NOT COSMETIC AND IT IS WHY `radius` IS A PARAMETER RATHER
	 * THAN FOLDED INTO F.** `F/R` is `mu_0 j_phi`, so a source that does not
	 * vanish on the symmetry axis is an infinite current density there -- which
	 * MEQ refuses at startup rather than integrating. A kernel that receives
	 * `radius` can be given the same refusal; one that receives `F/R`
	 * pre-divided cannot tell a large value from a pole.
	 *
	 * @param psi     the SOLVED potential at the point, already summed against
	 *                the shape functions.
	 * @param psiShift `psi_c` under COIL-SUBTRACTION-PLAN.md's split, or zero.
	 *                **The source must see the TOTAL flux and the solver only
	 *                has the remainder**, so this is added before F is
	 *                evaluated and never afterwards. Zero makes the split
	 *                exactly inert.
	 * @param inPlasma XP-1's per-ELEMENT connectivity flag, which is the
	 *                solver's and cannot be computed from a point value. False
	 *                turns the plasma term off for the whole element; the
	 *                pointwise test inside F still runs independently.
	 */
	MFEM_HOST_DEVICE inline mfem::real_t sourceResidualAtPoint(
		SourceView const &view, KnotView const *knots,
		mfem::real_t radius, mfem::real_t z,
		mfem::real_t psi, mfem::real_t psiShift,
		mfem::real_t weight, bool inPlasma )
	{
		if ( !inPlasma )
			return mfem::real_t( 0 );
		return -weight*sourceF( view, knots, radius, z, psi + psiShift )/radius;
	}

	/// The matching Jacobian contribution, `-weight * dF/dpsi / R`, for
	/// `AssembleElementGrad`. Same parameters and the same shift, which is the
	/// half that is easy to talk yourself out of: no chain-rule FACTOR is owed
	/// since `d( psi_c + psi_p )/d psi_p` is 1, but `dF/dpsi` is EVALUATED AT a
	/// flux, and a Jacobian taken at the wrong point is invisible to a
	/// convergence table.
	MFEM_HOST_DEVICE inline mfem::real_t sourceGradientAtPoint(
		SourceView const &view, KnotView const *knots,
		mfem::real_t radius, mfem::real_t z,
		mfem::real_t psi, mfem::real_t psiShift,
		mfem::real_t weight, bool inPlasma )
	{
		if ( !inPlasma )
			return mfem::real_t( 0 );
		return -weight*sourceDFdPsi( view, knots, radius, z, psi + psiShift )/radius;
	}

	/**
	 * Builds a SourceView for @a source, or refuses.
	 *
	 * @param source  the live source. Read, not held.
	 * @param view    written only on success.
	 * @param knots   appended to only on success; @a view indexes into it, so
	 *                a caller building several views must keep one array and
	 *                pass the same one each time.
	 * @return false when @a source has no POD form. **A refusal is the safe
	 *         answer and the caller keeps the host path**; there is no
	 *         approximate mode, because a SourceView that is nearly right is a
	 *         different equilibrium reached at the full convergence rate.
	 *
	 * Declared here and defined in SourceView.cpp, which is ordinary host code
	 * -- it reads a polymorphic hierarchy and allocates, neither of which
	 * belongs in a device translation unit.
	 */
	bool makeSourceView( Source const &source, SourceView &view,
	                     std::vector<KnotView> &knots );

	/**
	 * A view of the PLASMA TERM ALONE, which is what meq::Source::scaledF()
	 * returns and what the prescribed-current constraint is about.
	 *
	 * **THIS IS A SECOND ENTRY POINT AND NOT A WIDENED CAST, AND THE
	 * DISTINCTION IS THE WHOLE SAFETY ARGUMENT.** makeSourceView()'s
	 * `dynamic_cast` to the exact production class is deliberate: its own
	 * comment records that accepting a derived class *"would build a view
	 * describing the plasma term alone and silently drop the rest"*. That
	 * reasoning is exactly right for `f()`, where a coil-augmented source's
	 * value is the SUM -- and exactly inverted for `scaledF()`, where the
	 * plasma term alone is the answer and the coil term is what must not be
	 * there. `CoilAugmentedNormalisedSource::scaledF()`'s own comment says so:
	 * *"this class's own f() is the SUM and would put the coil current inside
	 * the prescribed plasma current"*.
	 *
	 * So nothing is dropped silently here, because the caller asked for the
	 * plasma term BY NAME. A caller that wants `f()` must keep using
	 * makeSourceView() and keep taking its refusal.
	 *
	 * **WHY IT EXISTS**: `assemblePlasmaCurrent()` integrates
	 * `scaledF( R, z, psi )/R` over the plasma support, once per constraint
	 * evaluation, and it is the third item of DEVICE-PLAN.md section 3.4a's
	 * reduction list. Every production machine case in this tree is
	 * coil-augmented, so without this the POD covers the fixture ladder and
	 * refuses every case MEQ exists to solve.
	 *
	 * @return false on anything neither wrapper nor exact production class,
	 *         exactly as makeSourceView() does, and for the same reason.
	 */
	bool makePlasmaSourceView( Source const &source, SourceView &view,
	                           std::vector<KnotView> &knots );

	/**
	 * The conductors of @a source as plain data, or an empty list.
	 *
	 * @param permeability written with the set's own `mu0()`, so the caller
	 *        need not decide which of the several mu0s in this tree it is.
	 * @return false only on a source this cannot describe. **A source with NO
	 *         conductors is a success with an empty list**, not a refusal --
	 *         the coil term is then exactly zero and every non-coil path stays
	 *         bit-identical.
	 */
	bool makeCoilViews( Source const &source, std::vector<CoilView> &coils,
	                    mfem::real_t &permeability );

	/**
	 * `F` and `dF/dpsi` at @a count points, on the device.
	 *
	 * Returns false when there is no device path -- no mfem::Device, or MEQ
	 * built without the device translation unit -- and the caller then loops
	 * `sourceF()`/`sourceDFdPsi()` itself, which is the same arithmetic.
	 * Defined in SourceViewDevice.cpp, which is the one file carrying the
	 * device compile of this header.
	 */
	bool evaluateSourceOnDevice( SourceView const &view,
	                             KnotView const *knots, int knotCount,
	                             mfem::real_t const *radius,
	                             mfem::real_t const *z,
	                             mfem::real_t const *psi,
	                             mfem::real_t *value,
	                             mfem::real_t *derivative,
	                             int count );

}

#endif
