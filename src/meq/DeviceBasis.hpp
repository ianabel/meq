#ifndef MEQ_DEVICEBASIS_HPP
#define MEQ_DEVICEBASIS_HPP

#include "mfem.hpp"

/**
 * @file
 *
 * EVALUATING A GRID FUNCTION AT AN ARBITRARY REFERENCE POINT, ON A DEVICE.
 *
 * **THIS IS THE ONE THING STANDING BETWEEN MEQ AND A RESIDENT SOLVE, AND TWO
 * PLACES IN THE TREE ALREADY SAID SO.** `src/meq/PointLocatorDevice.cpp`
 * opens by naming what it cannot do -- *"IT LOCATES, IT DOES NOT INTERPOLATE.
 * Evaluating the field needs `FiniteElement::CalcShape`, which is a host
 * virtual taking an `mfem::IntegrationPoint`. A device evaluation needs a
 * `MFEM_HOST_DEVICE` basis"* -- and `DEVICE-PLAN.md` section 3.4 calls the
 * axis constraint's mesh search the hardest item on the page for the same
 * reason. They are the same blocker seen from two sides, and this file is it.
 *
 * **WHY A FIXED QUADRATURE RULE IS NOT THE ANSWER, WHICH IS WHY THIS EXISTS AT
 * ALL.** MFEM already hands shape functions to a device through
 * `FiniteElement::GetDofToQuad()`, and `meq::plasmaCurrentPerElementOnDevice()`
 * uses exactly that. It supplies `phi_i` at the points of ONE
 * `mfem::IntegrationRule`, decided before the launch. That is everything an
 * integral needs and nothing a SEARCH needs: `CriticalPointFinder`'s
 * element-local Newton evaluates `q` wherever its iterate has got to, and is
 * documented to let that iterate leave the reference element by up to two
 * element diameters. *A rule is a set of points; a search needs a function.*
 *
 * HOW IT WORKS, AND IT BORROWS NOTHING FROM MFEM'S INTERNALS.
 *
 * Every polynomial on a triangle is a polynomial, whatever basis a
 * `FiniteElement` happens to express it in. So the host asks the element for
 * its own shape functions at its own nodes -- `CalcShape`, once, at build time
 * -- and solves for the constant matrix `T` that re-expresses them in a basis
 * a kernel can evaluate with four lines of arithmetic:
 *
 *     u( xi ) = sum_i u_i phi_i( xi ) = sum_b c_b lambda^alpha_b,
 *     c_b = sum_i T[ b ][ i ] u_i.
 *
 * `T` depends on the geometry, the degree and the basis type and on **nothing
 * about the solution**, so it is built once and uploaded once. That makes it
 * `DEVICE-PLAN.md` section 3.3's phase B in its purest form -- a constant of
 * the discretisation, not a state -- and M-226's warning applies in full: a
 * fill that is safe in a constructor is not safe in a refresh, so
 * meq::DeviceTriangleBasis::build() is a construction and there is no
 * in-place refresh of it to get wrong.
 *
 * **THE KERNEL'S BASIS IS THE BARYCENTRIC MONOMIALS AND THE BERNSTEIN
 * NORMALISATION IS FOLDED INTO `T`.** `M` is formed against the genuine
 * Bernstein polynomials -- `C( p; i, j, k ) lambda0^i lambda1^j lambda2^k` --
 * because those are a well-conditioned basis to invert against, and the
 * multinomial factor is then multiplied into `T`'s rows so the kernel
 * evaluates plain powers. The conditioning is paid on the host once; the
 * device does `nd` multiply-adds per term and no factorials.
 *
 * **AND PLAIN POWERS ARE THE RIGHT FORM FOR THE CALLER THIS WAS WRITTEN FOR,
 * NOT MERELY THE CHEAPEST.** `lambda0 = 1 - xi - eta` goes NEGATIVE outside
 * the reference triangle, which is where `rootInElement()` deliberately
 * spends part of its Newton. `lambda0^i` by repeated multiplication is exactly
 * the polynomial's own continuation there; `std::pow` of a negative base is
 * not something to reach for and a basis assembled from clamped barycentrics
 * would be a different function outside the element. *A term added to an
 * extrapolating evaluation must be defined wherever the extrapolation goes* --
 * CLAUDE.md records that rule being learned the other way round, from a coil
 * term that was not.
 *
 * WHAT IT DOES NOT CLAIM.
 *
 *   **STRAIGHT-SIDED TRIANGLES.** The change of basis is about the reference
 *   element and is correct for any shape, but every caller here also wants the
 *   physical map, and meq::AffineTriangle is affine by construction.
 *   meq::DeviceTriangleBasis::build() refuses a geometry that is not a
 *   triangle and the callers refuse a curved mesh, each falling back to the
 *   host rather than approximating.
 *
 *   **IT IS NOT BIT-IDENTICAL TO `CalcShape`, AND NOTHING HERE PRETENDS
 *   OTHERWISE.** A change of basis is a dense matrix multiply, so the answer
 *   agrees with MFEM's to round-off and not to the bit -- unlike
 *   meq::DeviceRows, whose whole design is bent around reproducing a host
 *   summation order. The difference is which contract is available: a row dot
 *   has an obvious host ordering to match, and a polynomial evaluated through
 *   a different basis does not. `tests/unit/DeviceBasisTests.cpp` therefore
 *   asserts a TOLERANCE and says what it is measured at, which is the same
 *   stance `tests/unit/DeviceQuadratureTests.cpp` takes and for the same
 *   reason.
 *
 *   **AND THE CONSEQUENCE FOR THE ANSWER IS BOUNDED BY AN ENVELOPE ARGUMENT
 *   RATHER THAN HOPED FOR.** What this evaluator feeds is a search for a
 *   STATIONARY point, and M-220 measured that pairing: a stationary point's
 *   POSITION costs half the available digits and its VALUE costs none, because
 *   the field is flat there. So a basis perturbation at `eps` moves the
 *   located axis at `sqrt( eps )` -- 1.5e-08, against the six decimals the
 *   driver prints for a position -- and moves `psi_ax`, which is the number
 *   the constraint actually carries, at `eps`.
 */

namespace meq
{

	/**
	 * The largest polynomial degree this evaluator will take.
	 *
	 * **FOUR, AND IT IS A CONDITIONING BOUND RATHER THAN AN ARRAY SIZE.** The
	 * first draft said six on the grounds that the ceiling is the per-thread
	 * coefficient array and a kernel that spills is slow rather than wrong.
	 * That is true inside the element and false outside it, which is where
	 * this evaluator's caller spends part of its Newton: measured against
	 * `CalcShape` at two reference units out -- `rootInElement()`'s own
	 * abandon bound -- the worst ABSOLUTE disagreement grows by about twenty
	 * per degree, reading **4.4e-16, 1.6e-14, 6.3e-13 and 1.27e-11** at
	 * degrees 1 to 4 where inside the containment band it is 1.9e-16 to
	 * 5.3e-15. Taken with the ceiling temporarily at six, degrees 5 and 6 read
	 * **1.6e-09 and 4.8e-07**. `lambda0` is negative out there and
	 * `lambda0^p` carries the growth, while MFEM's own recursion degrades far
	 * more slowly. `tests/unit/DeviceBasisTests.cpp` prints the first four
	 * every run, so they are a measurement this file keeps rather than one it
	 * remembers.
	 *
	 * At degree 6 that figure is within sight of a stopping target of
	 * `1e-13 * max | q |`, and a search whose residual cannot fall below its
	 * target is a search that runs to its iteration cap and reports nothing.
	 * **So the ceiling is set where the evaluator was MEASURED rather than
	 * where its arrays end**, every degree it accepts is swept by
	 * `tests/unit/DeviceBasisTests.cpp`, and above it
	 * DeviceTriangleBasis::build() refuses and the caller keeps the host path
	 * -- meq::makeSourceView()'s contract, for its reason.
	 *
	 * MEQ's production configuration is `k = 2` and the convergence suite
	 * sweeps `k = 1` to `4`, so nothing in the tree meets the refusal today.
	 * **`psi*` WOULD**: it is `k + 2`, so a basis built on the post-processed
	 * space reaches 6 at `k = 4`. Raising the ceiling for that means a better
	 * conditioned device basis -- a Bernstein evaluation rather than a
	 * monomial one, which costs the kernel the multinomial fold this file
	 * currently does on the host -- and not a larger array.
	 */
	int constexpr triangleBasisMaxDegree = 4;

	/// Dofs of `P_degree` on a triangle: the length of a coefficient array.
	int constexpr triangleBasisMaxDofs =
		( triangleBasisMaxDegree + 1 )*( triangleBasisMaxDegree + 2 )/2;

	/// `( degree + 1 )( degree + 2 )/2`, as a host or device expression.
	MFEM_HOST_DEVICE inline int triangleBasisDofs( int degree )
	{
		return ( degree + 1 )*( degree + 2 )/2;
	}

	/**
	 * The per-element affine map, as plain arithmetic types.
	 *
	 * `x = origin + jacobian * ( xi, eta )`, and `inverseJacobian` is that
	 * 2x2's inverse so a kernel can push a reference gradient to a physical
	 * one without a solve. Both are held rather than one being derived,
	 * because the derivation is three flops and a branch on a determinant that
	 * the host has already checked.
	 *
	 * NO DEFAULT MEMBER INITIALISERS, for meq::ElementMap's reason:
	 * `mfem::Array<T>` static_asserts that `T` is trivial, and an initialiser
	 * here is a compile error in the device translation unit and nowhere else.
	 * Construct with `AffineTriangle{}`, which zeroes an aggregate.
	 */
	struct AffineTriangle
	{
		mfem::real_t origin[ 2 ];
		mfem::real_t jacobian[ 4 ];
		mfem::real_t inverseJacobian[ 4 ];
		mfem::real_t determinant;
		/// False for a curved or degenerate element; such an element is skipped.
		bool affine;
	};

	/**
	 * `x = origin + J ( xi, eta )`. The forward map, exactly.
	 */
	MFEM_HOST_DEVICE inline void affineTransform( AffineTriangle const &map,
	                                              mfem::real_t xi, mfem::real_t eta,
	                                              mfem::real_t &radius, mfem::real_t &z )
	{
		radius = map.origin[ 0 ] + map.jacobian[ 0 ]*xi + map.jacobian[ 1 ]*eta;
		z      = map.origin[ 1 ] + map.jacobian[ 2 ]*xi + map.jacobian[ 3 ]*eta;
	}

	/**
	 * The change of basis, as a POD a kernel can capture.
	 *
	 * `change` is `dofs x dofs`, row major, row `b` holding the coefficients
	 * that take a nodal dof vector to the `b`-th barycentric monomial's
	 * coefficient. It points into device memory owned by a
	 * meq::DeviceTriangleBasis and is valid for that object's lifetime.
	 */
	struct TriangleBasisView
	{
		int degree;
		int dofs;
		mfem::real_t const *change;
		/**
		 * The element's own nodes, `2*dofs` reference coordinates
		 * interleaved `( xi, eta )`.
		 *
		 * Carried because a SEARCH has to start somewhere and
		 * `meq::CriticalPointFinder::elementSeeds()` starts at the node
		 * carrying the smallest `| q |` -- which is `mfem::FiniteElement::
		 * GetNodes()`, a host `mfem::IntegrationRule`. It is a constant of the
		 * element exactly as `change` is, so it rides along rather than
		 * becoming a second upload with a second lifetime to get wrong.
		 */
		mfem::real_t const *nodes;
	};

	/**
	 * The coefficients of one component, from its nodal dofs.
	 *
	 * @param basis the change of basis.
	 * @param dofs  `basis.dofs` nodal values of this element and component.
	 * @param out   `basis.dofs` monomial coefficients.
	 *
	 * Summed in ASCENDING dof order, which is not a bit-exactness claim
	 * against anything -- there is no host loop with this shape to match -- but
	 * is what makes two runs of the same kernel agree with each other.
	 */
	/**
	 * The same, with the dof count a COMPILE-TIME constant.
	 *
	 * **THE TEMPLATE IS NOT A TIDINESS EXERCISE: IT IS WHERE THE
	 * COEFFICIENTS LIVE.** With a runtime bound nvcc cannot unroll these
	 * loops, so `out[ b ]` is a dynamic index and the array it belongs to is
	 * placed in LOCAL MEMORY rather than in registers -- measured with
	 * `-Xptxas -v`, the root-finding kernel carried a **608 byte stack
	 * frame** with zero register spills, which is what that combination
	 * means. At 1024 threads an SM then has 622 KB of local working set
	 * against 64 KB of L1. Constant bounds unroll, the indices become
	 * constants, and the arrays become registers.
	 */
	template <int DOFS>
	MFEM_HOST_DEVICE inline void triangleBasisCoefficients(
		mfem::real_t const *change, mfem::real_t const *dofs,
		mfem::real_t *out )
	{
		MFEM_UNROLL( DOFS )
		for ( int b = 0; b < DOFS; ++b )
		{
			mfem::real_t total = 0.0;
			mfem::real_t const *row = change + b*DOFS;
			MFEM_UNROLL( DOFS )
			for ( int i = 0; i < DOFS; ++i )
				total += row[ i ]*dofs[ i ];
			out[ b ] = total;
		}
	}

	/// The runtime form, which dispatches to the template above so that the
	/// arithmetic is written once. The host path and the tests take this one.
	MFEM_HOST_DEVICE inline void triangleBasisCoefficients(
		TriangleBasisView const &basis, mfem::real_t const *dofs,
		mfem::real_t *out )
	{
		switch ( basis.dofs )
		{
			case 1:  triangleBasisCoefficients<  1 >( basis.change, dofs, out ); return;
			case 3:  triangleBasisCoefficients<  3 >( basis.change, dofs, out ); return;
			case 6:  triangleBasisCoefficients<  6 >( basis.change, dofs, out ); return;
			case 10: triangleBasisCoefficients< 10 >( basis.change, dofs, out ); return;
			case 15: triangleBasisCoefficients< 15 >( basis.change, dofs, out ); return;
			default: break;
		}
		for ( int b = 0; b < basis.dofs; ++b )
			out[ b ] = 0.0;
	}

	/**
	 * The value at a reference point, and its reference gradient.
	 *
	 * @param degree the polynomial degree the coefficients belong to.
	 * @param coefficients `triangleBasisDofs( degree )` monomial coefficients,
	 *        as triangleBasisCoefficients() produces them.
	 * @param dValueDxi, dValueDeta may be null together, in which case only the
	 *        value is computed -- which is the case the axis search wants,
	 *        since it differences the value rather than reading the gradient.
	 *
	 * Defined for ANY `( xi, eta )`, inside the reference triangle or not.
	 * See the file comment: the caller this was written for spends part of its
	 * Newton outside, deliberately, and a clamped evaluation would silently be
	 * a different function there.
	 */
	/**
	 * The value at a reference point and its reference gradient, with the
	 * degree a COMPILE-TIME constant.
	 *
	 * See triangleBasisCoefficients() for why the template exists: a runtime
	 * degree leaves `coefficients[ at ]` a dynamic index, which puts the
	 * caller's array in local memory rather than in registers.
	 */
	template <int DEGREE>
	MFEM_HOST_DEVICE inline void evaluateTriangleBasis(
		mfem::real_t const *coefficients,
		mfem::real_t xi, mfem::real_t eta,
		mfem::real_t &value, mfem::real_t *dValueDxi,
		mfem::real_t *dValueDeta )
	{
		mfem::real_t const lambda0 = mfem::real_t( 1.0 ) - xi - eta;
		mfem::real_t const lambda1 = xi;
		mfem::real_t const lambda2 = eta;

		// Powers by repeated multiplication rather than pow(): the bases go
		// negative outside the element and the exponents are small integers.
		mfem::real_t p0[ DEGREE + 1 ];
		mfem::real_t p1[ DEGREE + 1 ];
		mfem::real_t p2[ DEGREE + 1 ];
		p0[ 0 ] = 1.0;
		p1[ 0 ] = 1.0;
		p2[ 0 ] = 1.0;
		MFEM_UNROLL( DEGREE )
		for ( int e = 1; e <= DEGREE; ++e )
		{
			p0[ e ] = p0[ e - 1 ]*lambda0;
			p1[ e ] = p1[ e - 1 ]*lambda1;
			p2[ e ] = p2[ e - 1 ]*lambda2;
		}

		mfem::real_t total = 0.0;
		mfem::real_t totalXi = 0.0;
		mfem::real_t totalEta = 0.0;
		bool const wantGradient = dValueDxi != nullptr && dValueDeta != nullptr;

		int at = 0;
		MFEM_UNROLL( DEGREE + 1 )
		for ( int j = 0; j <= DEGREE; ++j )
			for ( int k = 0; k <= DEGREE - j; ++k, ++at )
			{
				int const i = DEGREE - j - k;
				mfem::real_t const c = coefficients[ at ];
				total += c*p0[ i ]*p1[ j ]*p2[ k ];

				if ( !wantGradient )
					continue;

				// d/dxi : lambda0 contributes -1, lambda1 +1, lambda2 0.
				// d/deta: lambda0 contributes -1, lambda2 +1, lambda1 0.
				mfem::real_t const dropped0 =
					i > 0 ? c*mfem::real_t( i )*p0[ i - 1 ]*p1[ j ]*p2[ k ]
					      : mfem::real_t( 0.0 );
				mfem::real_t const raised1 =
					j > 0 ? c*mfem::real_t( j )*p0[ i ]*p1[ j - 1 ]*p2[ k ]
					      : mfem::real_t( 0.0 );
				mfem::real_t const raised2 =
					k > 0 ? c*mfem::real_t( k )*p0[ i ]*p1[ j ]*p2[ k - 1 ]
					      : mfem::real_t( 0.0 );
				totalXi  += raised1 - dropped0;
				totalEta += raised2 - dropped0;
			}

		value = total;
		if ( wantGradient )
		{
			*dValueDxi = totalXi;
			*dValueDeta = totalEta;
		}
	}

	/**
	 * The value at a reference point, and its reference gradient.
	 *
	 * @param degree the polynomial degree the coefficients belong to.
	 * @param coefficients `triangleBasisDofs( degree )` monomial coefficients,
	 *        as triangleBasisCoefficients() produces them.
	 * @param dValueDxi, dValueDeta may be null together, in which case only the
	 *        value is computed -- which is the case the axis search wants,
	 *        since it differences the value rather than reading the gradient.
	 *
	 * Defined for ANY `( xi, eta )`, inside the reference triangle or not.
	 * See the file comment: the caller this was written for spends part of its
	 * Newton outside, deliberately, and a clamped evaluation would silently be
	 * a different function there.
	 *
	 * A dispatcher onto the template above, so the arithmetic is written once.
	 * The bound is checked HERE and not only in build(), because the three
	 * stack arrays are sized by it, this function is public and
	 * MFEM_HOST_DEVICE, and a MFEM_HOST_DEVICE function cannot throw or
	 * assert -- so the answer out of range is a defined zero.
	 */
	MFEM_HOST_DEVICE inline void evaluateTriangleBasis(
		int degree, mfem::real_t const *coefficients,
		mfem::real_t xi, mfem::real_t eta,
		mfem::real_t &value, mfem::real_t *dValueDxi,
		mfem::real_t *dValueDeta )
	{
		switch ( degree )
		{
			case 0: evaluateTriangleBasis< 0 >( coefficients, xi, eta, value,
			                                    dValueDxi, dValueDeta ); return;
			case 1: evaluateTriangleBasis< 1 >( coefficients, xi, eta, value,
			                                    dValueDxi, dValueDeta ); return;
			case 2: evaluateTriangleBasis< 2 >( coefficients, xi, eta, value,
			                                    dValueDxi, dValueDeta ); return;
			case 3: evaluateTriangleBasis< 3 >( coefficients, xi, eta, value,
			                                    dValueDxi, dValueDeta ); return;
			case 4: evaluateTriangleBasis< 4 >( coefficients, xi, eta, value,
			                                    dValueDxi, dValueDeta ); return;
			default: break;
		}
		value = 0.0;
		if ( dValueDxi != nullptr ) { *dValueDxi = 0.0; }
		if ( dValueDeta != nullptr ) { *dValueDeta = 0.0; }
	}

	/**
	 * The host half: builds the change of basis, and owns it on the device.
	 *
	 * One object per `( geometry, degree, basis type )`, which in MEQ means one
	 * per finite element space. It holds an `mfem::Vector` left device-valid
	 * after build() and hands out a TriangleBasisView onto it.
	 *
	 * **THERE IS NO REFRESH, AND THAT IS THE M-226 LESSON APPLIED RATHER THAN
	 * RE-LEARNED.** The matrix is a constant of the discretisation, so the only
	 * event that can invalidate it is a new space -- at which point the caller
	 * builds a new object. A `refresh()` that filled this in place would be
	 * exactly the shape that shipped a neighbouring answer when phase B was
	 * hoisted out of its loop.
	 */
	class DeviceTriangleBasis
	{
		public:
			/**
			 * Builds the change of basis for @a space's element 0.
			 *
			 * @return false, having touched nothing, when the space is not
			 *         scalar-shaped triangles of a degree this file takes, or
			 *         when the inversion is too ill-conditioned to trust. **A
			 *         refusal is the safe answer and the caller keeps the host
			 *         path.**
			 *
			 * Element 0 because an `L2_FECollection` gives every element of one
			 * geometry the same element object; build() checks that by asking
			 * the mesh whether it is of a single geometry rather than assuming
			 * it, since a mixed mesh would otherwise be read off its first
			 * triangle and applied to its quadrilaterals.
			 */
			bool build( mfem::FiniteElementSpace const &space );

			/// Whether build() succeeded.
			bool valid() const;

			/**
			 * The view **for a kernel**: its pointers are `Read()`'s, which is
			 * a DEVICE address whenever an `mfem::Device` is configured.
			 *
			 * **DEREFERENCING THEM ON THE HOST IS A SEGFAULT UNDER `cuda` AND
			 * IS NOT CAUGHT BY `debug`**, which is why the two views are
			 * separate functions rather than one with a flag. Measured: a host
			 * loop reading `view().change` under `--device cuda` dies at the
			 * first access, while under `--device debug` it survives, because
			 * debug's device allocation is ordinary host memory and only the
			 * host copy is protected. That is M-206's recorded limit -- the
			 * debug backend catches a raw access to the WRONG LOCATION and is
			 * blind to an honest flag over a stale copy -- met here as an
			 * honest flag over a pointer into the wrong address space.
			 *
			 * `meq::DeviceRows` states the same policy the other way round and
			 * it is worth following rather than restating: it hands its
			 * kernel the CONTAINERS and lets `Read()` decide where they live.
			 * This class cannot, the view being a POD a kernel captures, so
			 * the next best thing is for the type to say which side it is for.
			 *
			 * Zeroed when `!valid()`, with `degree` set to -1 so that a caller
			 * testing it cannot mistake a refusal for degree 0.
			 */
			TriangleBasisView view() const;

			/**
			 * The same view with HOST pointers, for a host caller.
			 *
			 * `meq::rootInTriangle()` and `meq::evaluateTriangleBasis()` are
			 * MFEM_HOST_DEVICE and compile into both halves; this is what lets
			 * the host half be exercised -- by a test, or by a fallback -- with
			 * a Device configured. It funnels through `meq::countedHostRead()`
			 * so a host read of a phase B upload is counted as the download it
			 * is rather than being invisible.
			 */
			TriangleBasisView hostView() const;

			DeviceTriangleBasis() = default;
			/// Neither copied nor assigned: both views hand out raw pointers
			/// into this object's own storage, so a copy would leave two
			/// objects whose views silently disagree about which is live.
			DeviceTriangleBasis( DeviceTriangleBasis const & ) = delete;
			DeviceTriangleBasis &operator=( DeviceTriangleBasis const & ) = delete;

			/// Dofs per element per component, or 0.
			int dofsPerElement() const;

			/// The polynomial degree, or -1.
			int degree() const;

			/**
			 * The worst relative error with which the change of basis
			 * reproduces `CalcShape`, over the reference triangle **and the
			 * containment band around it**.
			 *
			 * NOT over the element's own nodes: the change of basis is solved
			 * to be exact there, so measuring it there measures the linear
			 * solve's residual rather than whether the two polynomials agree.
			 * The lattice is interior points plus the band out to
			 * `meq::CriticalPointFinder`'s containment allowance, which is the
			 * furthest from its element a root may be ACCEPTED -- so this is
			 * the error over the region a caller's answer can come from, and
			 * it is what `build()` refuses on.
			 */
			double reconstructionError() const;

			/**
			 * The worst ABSOLUTE disagreement out at the abandon bound -- two
			 * reference units, where `rootInElement()` gives up.
			 *
			 * Absolute where reconstructionError() is relative, because the
			 * question here is whether the Newton's own stopping test can
			 * still fire: `residual <= tolerance * max | q |` is an absolute
			 * threshold against a global scale, and normalising by shape
			 * functions that are themselves growing out there would hide the
			 * thing being asked about.
			 *
			 * **REPORTED AND NOT BUDGETED, BECAUSE IT IS ABOUT THE PATH AND
			 * NOT ABOUT THE ANSWER.** A Newton iterate is allowed to wander
			 * that far and the barycentric monomials are badly conditioned
			 * there -- measured at 5.5e-11 at degree 4 against 8.5e-15 inside,
			 * and 4.8e-07 at degree 6, because `lambda0^p` grows like `|lambda0|^p`
			 * and MFEM's own recursion degrades far more slowly. What that
			 * costs is a few iterates taking a different route to the same
			 * root, which is why it is not a refusal; what it would cost at a
			 * degree where it approached the stopping target is a search that
			 * cannot converge, which is why the ceiling is where it is.
			 */
			double extrapolationError() const;

		private:
			mfem::Vector change;
			mfem::Vector nodes;
			int dofs = 0;
			int order = -1;
			double worstError = 0.0;
			double worstOutside = 0.0;
			bool built = false;
	};

	/**
	 * The affine map of every element of @a mesh, for a device kernel.
	 *
	 * @param out resized to `mesh.GetNE()` and left **device-valid**.
	 * @return false when the mesh has any element that is not a straight-sided
	 *         triangle, having filled @a out anyway with `affine` false on
	 *         those -- so a caller that can skip them may, and a caller that
	 *         cannot has one test to make. Under a Device, on a mesh of
	 *         straight triangles, the maps are built on the device and this is
	 *         true; an element with a zero determinant is still marked
	 *         `affine` false, and the kernels skip it as they skip a curved one.
	 *
	 * Geometry only: nothing here depends on the solution, so this is a phase B
	 * upload and is rebuilt only when the mesh is.
	 */
	bool affineTrianglesOf( mfem::Mesh &mesh, mfem::Array<AffineTriangle> &out );

	/**
	 * affineTrianglesOf()'s device half, out of MFEM's geometric factors at
	 * the reference corner: `X` is `( 1 x 2 x NE )` and `J` `( 1 x 2 x 2 x
	 * NE )`. In AxisSearchDevice.cpp; false without a device language.
	 */
	bool affineTrianglesKernel( mfem::Vector const &coordinates,
	                            mfem::Vector const &jacobians, int elements,
	                            mfem::Array<AffineTriangle> &out );

}

#endif
