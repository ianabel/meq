#ifndef MEQ_AXISSEARCH_HPP
#define MEQ_AXISSEARCH_HPP

#include "DeviceBasis.hpp"

#include "mfem.hpp"

#include <cmath>

/**
 * @file
 *
 * THE MAGNETIC AXIS, FOUND WITHOUT BRINGING THE ITERATE TO THE HOST.
 *
 * **WHAT THIS IS FOR, AND IT IS ONE LINE OF ONE TABLE.** With the bordered
 * step's solved columns made device-capable, a `--device cuda` run of
 * `machine-f-diiid` moves **20 transfers and 15.1 MB** where the host-resident
 * arrangement moves 3 and 0.888 -- and the per-site table puts **15 of those
 * transfers and 12.5 of those megabytes on one `countedHostRead( state )`**,
 * the funnel at the top of `GradShafranovSolver::solveWithNormalisation()`'s
 * `locateAxisPoint()`. Every other site in the run put together is the
 * remaining five. [M-229](../MEASUREMENTS.md) is the measurement and its
 * conclusion is why this file exists: *residency is a property of a cycle and
 * not of a leg*, so the one host leg converts every other leg's kernel into an
 * upload, and this is the one host leg.
 *
 * `DEVICE-PLAN.md` section 3.4 calls it the hardest item on the page, and
 * section 3.3's table has it as the last one. Both are right about which it
 * is and the reason is below.
 *
 * **WHY IT IS THE HARD ONE: IT IS A SEARCH AND EVERY OTHER ITEM IS A
 * REDUCTION.** A reduction visits a fixed set of points -- an
 * `mfem::IntegrationRule`'s -- so MFEM's own
 * `FiniteElement::GetDofToQuad()` hands it shape functions on the device and
 * `meq::plasmaCurrentPerElementOnDevice()` needs nothing further. A search
 * evaluates wherever its iterate has got to, and
 * `meq::CriticalPointFinder::rootInElement()` is documented to let that
 * iterate leave the reference element by up to two element diameters. *A rule
 * is a set of points; a search needs a function.* meq::DeviceBasis is that
 * function and this file is what it was written for.
 *
 * WHAT MOVES AND WHAT DOES NOT.
 *
 *   **THE ROOT FINDING MOVES; THE SELECTION DOES NOT.** The device roots every
 *   element -- one thread each -- and returns the ACCEPTED candidates,
 *   compacted. Which candidate becomes the axis is then decided by the host,
 *   by the code that decides it today: the conductor filter, its two passes,
 *   the normalised-flux score, and the continuation from the previous axis.
 *   That is the half carrying the physics (*"a plasma has no magnetic axis
 *   inside a conductor"*, and *"following the extremum nearest the previous
 *   one is a continuation in the axis rather than a fresh competition between
 *   every O-point in the field at every step"*), it is a few dozen operations
 *   on a list of a handful, and moving it would buy nothing and risk the one
 *   thing this file must not change.
 *
 *   **AND THE DEDUPLICATION DOES NOT MOVE EITHER, BECAUSE THE SELECTION DOES
 *   NOT NEED IT.** `sweep()` merges candidates within `separation*diameter`
 *   of each other, serially and in element order, because doing it inside its
 *   OpenMP region would make the answer thread-order dependent. Both selection
 *   rules -- largest score, and nearest the previous axis -- are **invariant
 *   under duplicates**: two candidates at the same point score the same and
 *   are the same distance away, so whichever survives a merge is whichever
 *   wins without one. The merge is kept out rather than reproduced, and that
 *   is an argument rather than an omission.
 *
 * **THE ONE SEMANTIC DIFFERENCE, STATED RATHER THAN BURIED.** The host's warm
 * path is `tryFindCriticalPointFrom()`, a ring-by-ring walk out from the
 * element nearest the previous axis that stops at the **first ring** holding a
 * clean root. The device roots every element at once -- on a card the ring
 * limit buys nothing, since 4,848 Newtons and 19 Newtons are one launch -- so
 * the host picks the nearest clean root **globally** rather than within the
 * first productive ring. Those agree whenever the nearest clean root is in
 * that ring, which is the case the ring walk was written for, and can differ
 * where two O-points sit at comparable distances in different rings. *This is
 * a different algorithm with the same stated intent*, and the thing that
 * establishes it is not this comment: it is
 * `tests/convergence/DeviceParity.cpp`, where a different axis is a different
 * `psi_ax` and a different Newton count.
 *
 * WHAT IT REFUSES, AND A REFUSAL IS ALWAYS THE HOST PATH RATHER THAN AN
 * APPROXIMATION -- meq::makeSourceView()'s contract, for its reason.
 *
 *   * **A CONDUCTOR SPLIT.** `totalFlux()` adds `q_c` from a per-node table
 *     under `[conductors] Model`, and `fluxScale()` takes a second pass for
 *     it. Both are portable and neither is ported, because no shipped example
 *     that reaches this constraint uses the split -- so porting it would be
 *     writing a kernel nothing runs, and the tree already has five of those.
 *   * **A CURVED OR MIXED MESH**, which meq::AffineTriangle refuses by
 *     construction, and a degree above meq::triangleBasisMaxDegree.
 *   * **MORE ACCEPTED CANDIDATES THAN THE COMPACTION HOLDS.** The cap is a
 *     property of the download and not of the field, so exceeding it is
 *     reported rather than truncated -- a truncated candidate list is a
 *     competition missing an entrant, which is exactly how a coil's O-point
 *     would win one.
 */

namespace meq
{

	/**
	 * One root of `q_h` in one element, as a POD.
	 *
	 * The fields mirror meq::CriticalPoint's, because the host rebuilds that
	 * structure from this one and a field here that did not exist there would
	 * be a field the selection could not see.
	 *
	 * No default member initialisers: `mfem::Array<T>` static_asserts
	 * triviality, exactly as meq::ElementMap and meq::AffineTriangle record.
	 */
	struct AxisCandidate
	{
		mfem::real_t referenceX;
		mfem::real_t referenceY;
		mfem::real_t radius;
		mfem::real_t z;
		mfem::real_t psi;
		mfem::real_t fluxResidual;
		mfem::real_t determinant;
		mfem::real_t trace;
		mfem::real_t overshoot;
		int element;
		/// Matches meq::CriticalPointType's underlying value.
		int type;
		int accepted;
	};

	/**
	 * Everything one thread needs, as plain arithmetic types and pointers.
	 *
	 * `state` is the bordered Newton's iterate, device-side: the flux block
	 * first and the potential block at `potentialOffset`. **The flux is
	 * NEGATED on read** -- `mfem::DarcyForm` holds `-q`, and
	 * `CriticalPoints.hpp` records that a finder handed the raw block puts
	 * every Maximum where a Minimum is. Doing it in the kernel rather than
	 * through a scratch grid function is the whole saving: the host copy that
	 * negation used to need is the 12.5 MB.
	 */
	struct AxisSearchInputs
	{
		TriangleBasisView flux;
		TriangleBasisView potential;
		AffineTriangle const *maps;
		mfem::real_t const *state;
		/// The reference element's centre, the first of the two seeds
		/// `elementSeeds()` tries. Passed rather than written as `1/3` so that
		/// the kernel and `mfem::Geometries.GetCenter()` cannot disagree.
		mfem::real_t centreXi;
		mfem::real_t centreEta;
		/// `max | q |` over the flux dofs, already scaled by the tolerance.
		mfem::real_t const *target;
		int elements;
		/// `GetNDofs()` of the SCALAR flux space: component 1 sits that far on.
		int fluxScalarDofs;
		int potentialOffset;
		mfem::real_t jacobianStep;
		mfem::real_t containment;
		int maxIterations;
	};

	/**
	 * `mfem::Geometry::CheckPoint( TRIANGLE, ip, eps )`, as device arithmetic.
	 *
	 * Read out of `fem/geom.cpp:493` and its `internal::FuzzyLT`/`FuzzyGT` at
	 * `:420` rather than remembered, because the whole value of this routine
	 * is that it is the same predicate: the host accepts a root on exactly
	 * these three inequalities and a kernel that was off by a sign on the
	 * tolerance would accept roots of neighbouring elements.
	 */
	MFEM_HOST_DEVICE inline bool insideReferenceTriangle(
		mfem::real_t xi, mfem::real_t eta, mfem::real_t tolerance )
	{
		return !( xi < -tolerance || eta < -tolerance
		          || xi + eta > mfem::real_t( 1.0 ) + tolerance );
	}

	/**
	 * The element-local Newton on `q = 0`, ported whole.
	 *
	 * A line-for-line port of `meq::CriticalPointFinder::rootInElement()` and
	 * the `referenceJacobian()` it calls -- same seed, same central-difference
	 * Jacobian at the same step with the same stencil sliding, same Cramer
	 * solve, same unit step truncation, same two-diameter abandon, same
	 * `containment + 1e-9` acceptance, same overshoot, same degeneracy
	 * threshold, same classification. **The iteration is in REFERENCE
	 * coordinates throughout**, which is what makes the port exact where it
	 * matters: the physical map enters only at the end, and for a straight
	 * triangle that map is affine and bitwise constant over the element --
	 * measured, not assumed.
	 *
	 * @param coefficients `2*flux.dofs` monomial coefficients of this
	 *        element's flux, component 0 then component 1, as
	 *        triangleBasisCoefficients() produces them.
	 * @param potentialCoefficients `potential.dofs` of the potential.
	 * @return whether a root was found and accepted for this element.
	 */
	MFEM_HOST_DEVICE inline bool rootInTriangle(
		AxisSearchInputs const &in, int element,
		mfem::real_t const *coefficients,
		mfem::real_t const *potentialCoefficients,
		mfem::real_t seedXi, mfem::real_t seedEta,
		AxisCandidate &found )
	{
		AffineTriangle const &map = in.maps[ element ];
		int const degree = in.flux.degree;
		mfem::real_t const *second = coefficients + in.flux.dofs;
		mfem::real_t const target = *in.target;

		mfem::real_t xi = seedXi;
		mfem::real_t eta = seedEta;
		mfem::real_t first = 0.0, other = 0.0;

		bool converged = false;
		for ( int iteration = 0; iteration < in.maxIterations; ++iteration )
		{
			evaluateTriangleBasis( degree, coefficients, xi, eta,
			                       first, nullptr, nullptr );
			evaluateTriangleBasis( degree, second, xi, eta,
			                       other, nullptr, nullptr );

			mfem::real_t const residual =
				sqrt( first*first + other*other );
			if ( residual <= target )
			{
				converged = true;
				break;
			}

			// referenceJacobian(), inlined: a central difference per
			// direction, with the stencil slid until both ends are inside.
			mfem::real_t jacobian[ 2 ][ 2 ];
			mfem::real_t const base[ 2 ] = { xi, eta };
			for ( int d = 0; d < 2; ++d )
			{
				mfem::real_t shift = 0.0;
				for ( int attempt = 0; attempt < 3; ++attempt )
				{
					mfem::real_t hi[ 2 ] = { base[ 0 ], base[ 1 ] };
					mfem::real_t lo[ 2 ] = { base[ 0 ], base[ 1 ] };
					hi[ d ] += shift + in.jacobianStep;
					lo[ d ] += shift - in.jacobianStep;

					bool const hiIn = insideReferenceTriangle(
						hi[ 0 ], hi[ 1 ], mfem::real_t( 0.0 ) );
					bool const loIn = insideReferenceTriangle(
						lo[ 0 ], lo[ 1 ], mfem::real_t( 0.0 ) );

					if ( hiIn && loIn )
						break;
					if ( !hiIn )
						shift -= in.jacobianStep;
					else
						shift += in.jacobianStep;
				}

				mfem::real_t hi[ 2 ] = { base[ 0 ], base[ 1 ] };
				mfem::real_t lo[ 2 ] = { base[ 0 ], base[ 1 ] };
				hi[ d ] += shift + in.jacobianStep;
				lo[ d ] += shift - in.jacobianStep;

				mfem::real_t highFirst = 0.0, highOther = 0.0;
				mfem::real_t lowFirst = 0.0, lowOther = 0.0;
				evaluateTriangleBasis( degree, coefficients, hi[ 0 ], hi[ 1 ],
				                       highFirst, nullptr, nullptr );
				evaluateTriangleBasis( degree, second, hi[ 0 ], hi[ 1 ],
				                       highOther, nullptr, nullptr );
				evaluateTriangleBasis( degree, coefficients, lo[ 0 ], lo[ 1 ],
				                       lowFirst, nullptr, nullptr );
				evaluateTriangleBasis( degree, second, lo[ 0 ], lo[ 1 ],
				                       lowOther, nullptr, nullptr );

				jacobian[ 0 ][ d ] =
					( highFirst - lowFirst )/( mfem::real_t( 2.0 )*in.jacobianStep );
				jacobian[ 1 ][ d ] =
					( highOther - lowOther )/( mfem::real_t( 2.0 )*in.jacobianStep );
			}

			// solveTwoByTwo(), by Cramer. A singular Jacobian is a place to
			// stop rather than an error.
			mfem::real_t const determinant =
				jacobian[ 0 ][ 0 ]*jacobian[ 1 ][ 1 ]
				- jacobian[ 0 ][ 1 ]*jacobian[ 1 ][ 0 ];
			if ( !( fabs( determinant ) > mfem::real_t( 0.0 ) ) )
				return false;

			mfem::real_t const rhs[ 2 ] = { -first, -other };
			mfem::real_t step[ 2 ];
			step[ 0 ] = (  jacobian[ 1 ][ 1 ]*rhs[ 0 ]
			             - jacobian[ 0 ][ 1 ]*rhs[ 1 ] )/determinant;
			step[ 1 ] = ( -jacobian[ 1 ][ 0 ]*rhs[ 0 ]
			             + jacobian[ 0 ][ 0 ]*rhs[ 1 ] )/determinant;

			mfem::real_t const length =
				sqrt( step[ 0 ]*step[ 0 ] + step[ 1 ]*step[ 1 ] );
			/*
			 * isfinite(), WHICH IS WHAT rootInElement() CALLS AND IS THE ONLY
			 * SPELLING THAT SURVIVES A NARROWER real_t.
			 *
			 * The first draft wrote this as `length == length` for the NaN and
			 * `length > real_t( 1.0e300 )` for the infinity. The second half
			 * is CLAUDE.md's recorded `1e-300` trap at the other end of the
			 * range: **1e300 is above float's maximum**, so under
			 * `MFEM_PRECISION = single` the literal overflows to infinity and
			 * the test becomes `length > inf`, which is never true. *The guard
			 * does not loosen, it disappears* -- and a Newton that does not
			 * abandon on an infinite step walks out of its element with every
			 * printed number still looking plausible.
			 */
			if ( !isfinite( length ) )
				return false;
			if ( length > mfem::real_t( 1.0 ) )
			{
				step[ 0 ] /= length;
				step[ 1 ] /= length;
			}

			xi += step[ 0 ];
			eta += step[ 1 ];

			if ( !insideReferenceTriangle( xi, eta, mfem::real_t( 2.0 ) ) )
				return false;

			if ( length < mfem::real_t( 1.0e-15 ) )
			{
				evaluateTriangleBasis( degree, coefficients, xi, eta,
				                       first, nullptr, nullptr );
				evaluateTriangleBasis( degree, second, xi, eta,
				                       other, nullptr, nullptr );
				converged = true;
				break;
			}
		}

		if ( !converged )
			return false;
		if ( !insideReferenceTriangle( xi, eta,
		                               in.containment + mfem::real_t( 1.0e-9 ) ) )
			return false;

		mfem::real_t overshoot = 0.0;
		if ( !insideReferenceTriangle( xi, eta, mfem::real_t( 0.0 ) ) )
		{
			mfem::real_t const worst = -xi > -eta ? -xi : -eta;
			mfem::real_t const edge = xi + eta - mfem::real_t( 1.0 );
			overshoot = worst > edge ? worst : edge;
			if ( overshoot < mfem::real_t( 0.0 ) )
				overshoot = 0.0;
		}

		// The Jacobian AT THE ROOT, differenced again exactly as the host
		// does -- the loop above left `jacobian` at the previous iterate.
		mfem::real_t reference[ 2 ][ 2 ];
		{
			mfem::real_t const base[ 2 ] = { xi, eta };
			for ( int d = 0; d < 2; ++d )
			{
				mfem::real_t shift = 0.0;
				for ( int attempt = 0; attempt < 3; ++attempt )
				{
					mfem::real_t hi[ 2 ] = { base[ 0 ], base[ 1 ] };
					mfem::real_t lo[ 2 ] = { base[ 0 ], base[ 1 ] };
					hi[ d ] += shift + in.jacobianStep;
					lo[ d ] += shift - in.jacobianStep;
					bool const hiIn = insideReferenceTriangle(
						hi[ 0 ], hi[ 1 ], mfem::real_t( 0.0 ) );
					bool const loIn = insideReferenceTriangle(
						lo[ 0 ], lo[ 1 ], mfem::real_t( 0.0 ) );
					if ( hiIn && loIn )
						break;
					if ( !hiIn )
						shift -= in.jacobianStep;
					else
						shift += in.jacobianStep;
				}

				mfem::real_t hi[ 2 ] = { base[ 0 ], base[ 1 ] };
				mfem::real_t lo[ 2 ] = { base[ 0 ], base[ 1 ] };
				hi[ d ] += shift + in.jacobianStep;
				lo[ d ] += shift - in.jacobianStep;

				mfem::real_t highFirst = 0.0, highOther = 0.0;
				mfem::real_t lowFirst = 0.0, lowOther = 0.0;
				evaluateTriangleBasis( degree, coefficients, hi[ 0 ], hi[ 1 ],
				                       highFirst, nullptr, nullptr );
				evaluateTriangleBasis( degree, second, hi[ 0 ], hi[ 1 ],
				                       highOther, nullptr, nullptr );
				evaluateTriangleBasis( degree, coefficients, lo[ 0 ], lo[ 1 ],
				                       lowFirst, nullptr, nullptr );
				evaluateTriangleBasis( degree, second, lo[ 0 ], lo[ 1 ],
				                       lowOther, nullptr, nullptr );

				reference[ 0 ][ d ] =
					( highFirst - lowFirst )/( mfem::real_t( 2.0 )*in.jacobianStep );
				reference[ 1 ][ d ] =
					( highOther - lowOther )/( mfem::real_t( 2.0 )*in.jacobianStep );
			}
		}

		// dq/dx = ( dq/dxi )( dxi/dx ). The affine map's inverse Jacobian is
		// what MFEM's InverseJacobian() returns at any point of the element.
		mfem::real_t physicalJacobian[ 2 ][ 2 ] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
		for ( int i = 0; i < 2; ++i )
			for ( int j = 0; j < 2; ++j )
				for ( int m = 0; m < 2; ++m )
					physicalJacobian[ i ][ j ] +=
						reference[ i ][ m ]*map.inverseJacobian[ 2*m + j ];

		mfem::real_t const det =
			physicalJacobian[ 0 ][ 0 ]*physicalJacobian[ 1 ][ 1 ]
			- physicalJacobian[ 0 ][ 1 ]*physicalJacobian[ 1 ][ 0 ];
		mfem::real_t const tr =
			physicalJacobian[ 0 ][ 0 ] + physicalJacobian[ 1 ][ 1 ];

		mfem::real_t norm = 0.0;
		for ( int i = 0; i < 2; ++i )
			for ( int j = 0; j < 2; ++j )
			{
				mfem::real_t const a = fabs( physicalJacobian[ i ][ j ] );
				norm = a > norm ? a : norm;
			}
		mfem::real_t const degenerate = mfem::real_t( 1.0e-10 )*norm*norm;

		mfem::real_t radius = 0.0, z = 0.0;
		affineTransform( map, xi, eta, radius, z );

		mfem::real_t psi = 0.0;
		evaluateTriangleBasis( in.potential.degree, potentialCoefficients,
		                       xi, eta, psi, nullptr, nullptr );

		found.referenceX = xi;
		found.referenceY = eta;
		found.radius = radius;
		found.z = z;
		found.psi = psi;
		found.fluxResidual = sqrt( first*first + other*other );
		found.determinant = det;
		found.trace = tr;
		found.overshoot = overshoot;
		found.element = element;
		found.accepted = 1;

		// The same ladder as rootInElement()'s, and the same order: degenerate
		// first, then the sign of the determinant, then the trace. 0, 1, 2, 3
		// are meq::CriticalPointType's Maximum, Minimum, Saddle, Degenerate.
		if ( fabs( det ) <= degenerate )
			found.type = 3;
		else if ( det < mfem::real_t( 0.0 ) )
			found.type = 2;
		else if ( tr > mfem::real_t( 0.0 ) )
			found.type = 1;
		else
			found.type = 0;

		return true;
	}

	/**
	 * The largest accepted-candidate list this routine will return.
	 *
	 * A cap on the DOWNLOAD and not on the field, so exceeding it is a refusal
	 * rather than a truncation -- see the file comment. 256 against a handful
	 * observed; the compacted array is 14 KB, where the state it replaces is
	 * 830 KB.
	 */
	int constexpr axisCandidateCap = 256;

	/**
	 * `tolerance * max( | q |, 1 )` over the flux dofs, computed and LEFT on
	 * the device.
	 *
	 * @param out resized to 1 and left device-valid, so
	 *        meq::AxisSearchInputs::target points at it and the scale is never
	 *        downloaded. That is the point: a scalar brought to the host to be
	 *        handed straight back is a transfer, and this one would be taken
	 *        once per residual evaluation.
	 *
	 * A max is exactly associative, so unlike meq::DeviceRows' row dots there
	 * is nothing here a reassociation could cost -- and it is still summed in
	 * ascending dof order, because two runs of one kernel agreeing is worth
	 * having for free.
	 */
	/**
	 * The X-point border's field quantities, at ONE already-located point.
	 *
	 * **THE EASY HALF OF `DEVICE-PLAN.md` SECTION 3.4, AND IT IS THE LARGER
	 * LINE OF THE TABLE.** `refreshXPoint()` downloads the whole iterate --
	 * **15 transfers and 12.5 MB** of a `--device cuda` run of
	 * `machine-f-diiid`, one 830 KB download per call over twelve Newton
	 * steps and three support sweeps -- to come away with **eight doubles**.
	 * It is bigger than everything else in the run put together and it is not
	 * a search: `locateFieldPoint()` has already inverted the element map, and
	 * it reads no solution vector at all. One element, one point, one thread.
	 *
	 * **AND IT SHARES ITS PHASE B WITH THE AXIS SEARCH**, which is the whole
	 * reason the two live in one header: a change of basis per space and the
	 * affine map of every element are constants of the DISCRETISATION, so a
	 * second copy would be M-226's shape -- one upload counted twice and
	 * attributable to neither.
	 */
	struct XPointFieldInputs
	{
		TriangleBasisView flux;
		TriangleBasisView potential;
		AffineTriangle const *maps;
		mfem::real_t const *state;
		int element;
		/**
		 * The reference point, which MAY LIE OUTSIDE its element.
		 *
		 * `locateFieldPoint()` keeps a near miss within a fifth of an element,
		 * so the host evaluates the element's polynomial a little outside it
		 * on purpose -- and the device must evaluate THE SAME CONTINUATION.
		 * meq::evaluateTriangleBasis() takes plain powers of barycentrics that
		 * go negative there for exactly this reason; a basis assembled from
		 * clamped barycentrics would be a different function.
		 */
		mfem::real_t xi;
		mfem::real_t eta;
		/// `GetNDofs()` of the SCALAR flux space: component 1 sits that far on.
		int fluxScalarDofs;
		int potentialOffset;
	};

	/**
	 * `q`, `grad q` and `grad psi_h` at that point, all PHYSICAL.
	 *
	 * `fluxJacobian` is row major -- `d q_i / d x_j` is
	 * `fluxJacobian[ 2*i + j ]` -- which is
	 * `mfem::GridFunction::GetVectorGradient()`'s own layout, so the caller
	 * copies it out without transposing anything.
	 */
	struct XPointFieldValues
	{
		mfem::real_t flux[ 2 ];
		mfem::real_t fluxJacobian[ 4 ];
		mfem::real_t potentialGradient[ 2 ];
		/// 0 when the element was not affine; the caller keeps the host path.
		int ok;
	};

	/// @copydoc XPointFieldInputs
	bool xPointFieldOnDevice( XPointFieldInputs const &inputs,
	                          XPointFieldValues &out );

	bool axisFluxScaleOnDevice( mfem::Vector const &state, int scalarDofs,
	                            mfem::real_t tolerance, mfem::Vector &out );

	/**
	 * Whether `MEQ_AXIS_DEVICE=force` is set in the environment.
	 *
	 * **THE ROUTES ARE GATED ON `DeviceIsValid() && !HostIsValid()`
	 * ([M-232](../MEASUREMENTS.md)), AND THAT GATE WAS DECIDED ON CROSSINGS
	 * ALONE.** It is the right call for bus traffic and it makes the kernels
	 * unreachable on every shipped case -- which means the thing they were
	 * written for cannot be TIMED, and a leg nobody can time is a leg whose
	 * cost is an opinion. `constraint location` is **51% of a `--device cuda`
	 * solve** of the largest case in the tree, so that opinion is load
	 * bearing.
	 *
	 * An environment variable rather than a key, for `MEQ_STAGE_CHECKSUM`'s
	 * reason: it has to be reachable on a shipped `examples/*.toml` without
	 * editing the file, because the point is to measure the case somebody
	 * actually runs. Read once and cached, so a disabled check costs one
	 * predicted branch.
	 *
	 * **IT CHANGES WHICH ROUTE RUNS AND NOT WHAT EITHER COMPUTES.** Both
	 * arms are pinned against `CriticalPointFinder::sweep()` by
	 * `tests/unit/AxisSearchTests.cpp` and against each other by
	 * `DeviceParity`; what it buys is a clock reading, at a measured cost in
	 * transfers that M-232 records.
	 */
	bool axisDeviceRouteForced();

	/**
	 * The kernel halves, which live in the device translation unit.
	 *
	 * **THE SPLIT IS THE ONE meq::DeviceQuadrature ALREADY USES AND IT IS NOT
	 * TIDINESS.** `mfem::forall` does not compile under g++, so a file
	 * containing one is handed to nvcc -- and nvcc's command line does not
	 * parse under clang, which is why CMake excludes every such file from the
	 * `naming` check (CLAUDE.md records that exclusion costing a session when
	 * it was read as a rule being relaxed). Keeping the refusals, the
	 * marshalling and every name a reader looks up in an ORDINARY translation
	 * unit means all of it stays under the convention check and only the
	 * launches sit outside it.
	 *
	 * Both are no-ops returning false without a device language, so a
	 * host-only build is correct rather than short of a route.
	 */
	bool axisFluxScaleKernel( mfem::Vector const &state, int scalarDofs,
	                          mfem::real_t tolerance, mfem::Vector &out );
	/// @copydoc axisFluxScaleKernel
	bool axisCandidatesKernel( AxisSearchInputs const &inputs,
	                           mfem::Array<AxisCandidate> &perElement,
	                           mfem::Array<AxisCandidate> &out,
	                           mfem::Array<int> &total, int cap );
	/// @copydoc axisFluxScaleKernel
	bool xPointFieldKernel( XPointFieldInputs const &inputs,
	                        mfem::Array<XPointFieldValues> &out );

	/**
	 * Roots every element and returns the accepted candidates, compacted.
	 *
	 * @param out resized to meq::axisCandidateCap and **left host-valid**: it
	 *        is the answer and the host is what reads it. This is the one
	 *        download, and it is the whole of what the axis constraint costs
	 *        the bus.
	 * @param count how many of @a out are filled.
	 * @return false, having touched nothing, with no `mfem::Device`, with an
	 *         invalid basis, or when more candidates were accepted than the
	 *         cap holds -- in which case @a count reports how many there were,
	 *         so a caller can say so rather than guess.
	 */
	bool axisCandidatesOnDevice( AxisSearchInputs const &inputs,
	                             mfem::Array<AxisCandidate> &out, int &count );

}

#endif
