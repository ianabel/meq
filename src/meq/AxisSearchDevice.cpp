#include "AxisSearch.hpp"

#include "Crossings.hpp"

/*
 * THE KERNELS BEHIND meq::AxisSearch, AND A SEPARATE TRANSLATION UNIT FOR
 * meq::PointLocatorDevice's REASON RATHER THAN FOR STYLE.
 *
 * `mfem::forall` does not compile under g++ -- general/backends.hpp:144 is a
 * static_assert saying so -- and MEQ compiles with g++ throughout, so CMake
 * marks this file alone as CUDA when the linked MFEM has it. Without a device
 * language the two entry points below are stubs returning false, and
 * `locateAxisPoint()` keeps the host route it has always had. That is what
 * makes a host-only build correct rather than short of one.
 *
 * NONE OF THE ARITHMETIC IS HERE. `rootInTriangle()` is in the header,
 * MFEM_HOST_DEVICE, and is compiled twice -- once by g++ for the unit test's
 * reference arm and once by nvcc into the kernel below. meq::PointLocator's
 * `locateAffine()` is in exactly that arrangement for exactly that reason: it
 * is the only one in which the two halves cannot drift apart, and a search
 * that disagrees with itself between builds is an equilibrium that depends on
 * whether the machine has a card.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	/**
	 * How many threads the two-stage reductions here spread over.
	 *
	 * **A COMPROMISE BETWEEN OCCUPANCY AND A SERIAL TAIL, AND IT IS NEITHER
	 * TUNED NOR ARBITRARY.** Stage one gives each thread a contiguous chunk
	 * and stage two walks the per-chunk results on ONE thread, so the serial
	 * tail is exactly this number while each thread's work is `N` over it.
	 * The balance point is `sqrt( N )` -- about 1,100 for the 1.23 M flux
	 * dofs of the largest case in the tree -- and this card has 40 SMs that
	 * want a good deal more than that in flight. 8,192 costs a 8 K serial
	 * tail, which is microseconds, and keeps stage one wide.
	 */
	int constexpr reductionChunks = 8192;

	/**
	 * The root find, with the degree a compile-time constant.
	 *
	 * See the dispatcher in axisCandidatesKernel() for why this is a template
	 * and not a loop bound.
	 */
	template <int DEGREE>
	static void axisRootLaunch( AxisSearchInputs const &in,
	                            AxisCandidate *d_perElement )
	{
		// lowerCamelCase: the naming ctest cannot read this file --
		// nvcc's flags do not parse under clang -- but CLAUDE.md records
		// that exclusion as an instrument failing rather than a rule
		// being relaxed, so the convention still binds here.
		constexpr int basisDofs = ( DEGREE + 1 )*( DEGREE + 2 )/2;

		mfem::forall( in.selectedCount, [ = ] MFEM_HOST_DEVICE ( int s )
		{
			// The OUTPUT is indexed by the slot and the INPUT by the element,
			// which is what lets the compaction scan the neighbourhood rather
			// than the mesh. The list is ascending, so slot order is element
			// order and the merge's tie-break is the host's.
			int const e = in.selected[ s ];
			d_perElement[ 2*s ].accepted = 0;
			d_perElement[ 2*s + 1 ].accepted = 0;

			if ( e < 0 || e >= in.elements || !in.maps[ e ].affine )
				return;

			mfem::real_t local[ basisDofs ];
			mfem::real_t coefficients[ 2*basisDofs ];
			mfem::real_t potentialCoefficients[ basisDofs ];

			/*
			 * AN L2 SPACE LAYS ITS ELEMENT DOFS OUT BLOCKED -- `e*nd + i` --
			 * and a vdim-2 space is `Ordering::byNODES`, so component 1 is
			 * one `GetNDofs()` further on. Both were checked by running MFEM
			 * rather than read off its documentation, at every degree this
			 * file takes. That is what lets the kernel index the state
			 * directly and skip the dof list, which is the other host object
			 * it would otherwise need.
			 *
			 * NEGATED, because mfem::DarcyForm holds `-q`.
			 */
			for ( int c = 0; c < 2; ++c )
			{
				MFEM_UNROLL( basisDofs )
				for ( int i = 0; i < basisDofs; ++i )
					local[ i ] = -in.state[ c*in.fluxScalarDofs
					                        + e*basisDofs + i ];
				triangleBasisCoefficients< basisDofs >(
					in.flux.change, local, coefficients + c*basisDofs );
			}
			MFEM_UNROLL( basisDofs )
			for ( int i = 0; i < basisDofs; ++i )
				local[ i ] = in.state[ in.potentialOffset + e*basisDofs + i ];
			triangleBasisCoefficients< basisDofs >(
				in.potential.change, local, potentialCoefficients );

			// Seed one: the element centre.
			if ( rootInTriangle< DEGREE >(
				     in, e, coefficients, potentialCoefficients,
				     in.centreXi, in.centreEta, d_perElement[ 2*s ] ) )
				d_perElement[ 2*s ].accepted = 1;

			// Seed two: the flux node carrying the smallest | q |, which is
			// elementSeeds()' second seed and is why the basis carries its
			// element's nodes.
			int best = 0;
			mfem::real_t smallest = -1.0;
			MFEM_UNROLL( basisDofs )
			for ( int i = 0; i < basisDofs; ++i )
			{
				mfem::real_t const xi = in.flux.nodes[ 2*i ];
				mfem::real_t const eta = in.flux.nodes[ 2*i + 1 ];
				mfem::real_t first = 0.0, other = 0.0;
				evaluateTriangleBasis< DEGREE >( coefficients, xi, eta,
				                                 first, nullptr, nullptr );
				evaluateTriangleBasis< DEGREE >( coefficients + basisDofs, xi, eta,
				                                 other, nullptr, nullptr );
				mfem::real_t const size = sqrt( first*first + other*other );
				if ( smallest < mfem::real_t( 0.0 ) || size < smallest )
				{
					smallest = size;
					best = i;
				}
			}
			if ( rootInTriangle< DEGREE >(
				     in, e, coefficients, potentialCoefficients,
				     in.flux.nodes[ 2*best ], in.flux.nodes[ 2*best + 1 ],
				     d_perElement[ 2*s + 1 ] ) )
				d_perElement[ 2*s + 1 ].accepted = 1;
		} );
	}

	bool axisFluxScaleKernel( mfem::Vector const &state, int scalarDofs,
	                          mfem::real_t tolerance, mfem::Vector &out )
	{
		out.SetSize( 1 );
		out.UseDevice( true );

		int const chunks = scalarDofs < reductionChunks ? scalarDofs
		                                                : reductionChunks;
		mfem::Vector partial( chunks );
		partial.UseDevice( true );

		mfem::real_t const *const d_state = countedDeviceRead( state, __LINE__ );
		mfem::real_t *const d_partial = partial.Write();
		mfem::real_t *const d_out = out.Write();

		/*
		 * TWO STAGES, AND THE FIRST DRAFT OF THIS WAS ONE THREAD.
		 *
		 * It walked every flux dof on a single CUDA core, justified in the
		 * comment as *"thirty microseconds of serial walk"* -- an estimate
		 * made against the 29 K dofs of `machine-f-diiid` and wrong by forty
		 * times at `fixed-f-diiid`, where it is **1.23 million iterations on
		 * one core**. The same comment then argued that `mfem::reduce` with
		 * `MaxReducer` was *"the textbook case"* and that a max is exactly
		 * associative so there would be nothing to lose by it, and used one
		 * thread anyway. *Having the argument is not the same as taking it.*
		 *
		 * **IT IS NOT `mfem::reduce`, AND THE REASON IS THE BUS RATHER THAN
		 * THE ARITHMETIC.** That routine takes `T &res` -- a HOST reference --
		 * so the scalar would come to the host and have to be pushed back,
		 * which is a download and an upload per axis evaluation of the one
		 * number this kernel exists to keep on the device.
		 * `meq::AxisSearchInputs::target` points at `out`.
		 *
		 * **AND NOTHING HERE DEPENDS ON THE ORDER.** `fmax` is associative AND
		 * commutative and exact, so unlike `meq::DeviceRows`' row dots there
		 * is no summation order to reproduce: any grouping gives the same bits.
		 * That is what makes a two-stage reduction free here where it would
		 * cost that class its whole contract.
		 */
		mfem::forall( chunks, [ = ] MFEM_HOST_DEVICE ( int c )
		{
			mfem::real_t worst = 0.0;
			for ( int i = c; i < scalarDofs; i += chunks )
			{
				mfem::real_t const a = d_state[ i ];
				mfem::real_t const b = d_state[ scalarDofs + i ];
				mfem::real_t const size = sqrt( a*a + b*b );
				worst = size > worst ? size : worst;
			}
			d_partial[ c ] = worst;
		} );

		/*
		 * THE SIGN DOES NOT MATTER HERE AND IT DOES EVERYWHERE ELSE.
		 * `locateAxisPoint()` negates the flux block because mfem::DarcyForm
		 * holds `-q`; this is a magnitude, and `( -a )*( -a )` is `a*a` to the
		 * bit, so the raw block gives the same scale. Said rather than left to
		 * be noticed, because every other reader of this block must negate.
		 *
		 * **AND THE EXACTNESS BELONGS TO THE REDUCTION AND NOT TO THE TERM**,
		 * which a double build cannot tell you. `fluxScale()` reads each pair
		 * into `double` locals and rounds the root back to `real_t`; this
		 * does the whole term in `real_t`. Identical in double, one rounding
		 * apart per term under `MFEM_PRECISION = single` -- measured at
		 * 0.1 eps on the agreement, and `tests/unit/AxisSearchTests.cpp`
		 * asserts the two claims separately for that reason.
		 */
		mfem::forall( 1, [ = ] MFEM_HOST_DEVICE ( int )
		{
			mfem::real_t worst = 0.0;
			for ( int c = 0; c < chunks; ++c )
				worst = d_partial[ c ] > worst ? d_partial[ c ] : worst;
			d_out[ 0 ] = tolerance
				*( worst > mfem::real_t( 0.0 ) ? worst : mfem::real_t( 1.0 ) );
		} );

		return true;
	}

	/// How the slots divide between threads. Written once because the count
	/// half and the scatter half must agree about it exactly -- a chunking
	/// that differed between them would scatter into the wrong offsets.
	static int axisChunkCount( int slots )
	{
		return slots < reductionChunks ? slots : reductionChunks;
	}

	bool axisCandidatesCountKernel( AxisSearchInputs const &inputs,
	                                mfem::Array<AxisCandidate> &perElement,
	                                mfem::Array<int> &offsets,
	                                mfem::Array<int> &total )
	{
		perElement.UseDevice( true );
		total.UseDevice( true );

		AxisSearchInputs const in = inputs;
		AxisCandidate *const d_perElement = perElement.Write();
		int *const d_total = total.Write();

		/*
		 * ONE THREAD PER ELEMENT. This is the parallelism and it is the whole
		 * reason the ring walk goes away: a host sweep is 40x a seeded search
		 * ([M-100](../MEASUREMENTS.md)) and on a card 4,848 Newtons and 19
		 * Newtons are one launch, so the device can afford the answer the
		 * host has to approximate.
		 *
		 * **TEMPLATED ON THE DEGREE, AND THAT IS WHERE THE COEFFICIENTS
		 * LIVE.** The first draft sized the three per-thread arrays by
		 * `triangleBasisMaxDofs` and bounded every loop by `basis.dofs`, a
		 * runtime value. nvcc cannot unroll a runtime bound, so every
		 * `coefficients[ at ]` stayed a dynamic index and the arrays were
		 * placed in LOCAL MEMORY -- measured with `-Xptxas -v` as a **608
		 * byte stack frame with zero register spills**, which is exactly what
		 * that combination means. At 1024 threads an SM then carries 622 KB
		 * of local working set against 64 KB of L1. Halving
		 * `triangleBasisMaxDegree` as a diagnostic took the frame to 272
		 * bytes, accounting for the difference to the byte, which is what
		 * said the arrays were the whole of it.
		 *
		 * With the degree a compile-time constant the loops unroll, the
		 * indices become constants and the arrays become registers. The
		 * launch is dispatched on `flux.degree` below; MEQ builds both spaces
		 * from one `orderValue`, and a caller whose two degrees differ is
		 * refused rather than served by the wrong instantiation.
		 */
		switch ( in.flux.degree )
		{
			case 1: axisRootLaunch< 1 >( in, d_perElement ); break;
			case 2: axisRootLaunch< 2 >( in, d_perElement ); break;
			case 3: axisRootLaunch< 3 >( in, d_perElement ); break;
			case 4: axisRootLaunch< 4 >( in, d_perElement ); break;
			default: return false;
		}

		/*
		 * AND THE COMPACTION IS AN ORDER-PRESERVING STREAM COMPACTION, WHICH
		 * THE FIRST DRAFT DID ON ONE THREAD.
		 *
		 * The order is load bearing and that part of the first draft's
		 * argument was right: `sweep()` deduplicates OUTSIDE its OpenMP
		 * region, with a comment recording that doing it inside would make
		 * the answer thread-order dependent at the 1e-3 level, and
		 * `mergeCandidates()` keeps the FIRST of two equals -- so the order a
		 * candidate list arrives in is the order a tie is broken in. An
		 * atomic-append scan would hand the host a list ordered by which
		 * blocks finished first.
		 *
		 * **WHAT DOES NOT FOLLOW IS ONE THREAD.** Chunks are CONTIGUOUS and
		 * numbered in slot order, each thread walks its own chunk ascending,
		 * and the offsets come from an exclusive prefix sum over the chunk
		 * counts -- so the output is exactly the order a serial walk produces,
		 * every run, with the parallelism in the two passes that do the work.
		 * The serial part is the prefix sum alone: `reductionChunks` entries,
		 * against the 205,680 slots the first draft walked on one core at the
		 * largest case in the tree.
		 */
		int const slots = 2*in.selectedCount;
		int const chunks = axisChunkCount( slots );
		offsets.SetSize( chunks + 1 );
		offsets.UseDevice( true );
		int *const d_offsets = offsets.Write();

		// long, because `c*slots` overflows an int on a mesh this routine is
		// otherwise happy to take: 8192 chunks times two million slots.
		auto chunkLow = [ = ] MFEM_HOST_DEVICE ( int c )
		{
			return static_cast<int>(
				static_cast<long>( c )*slots/chunks );
		};

		mfem::forall( chunks, [ = ] MFEM_HOST_DEVICE ( int c )
		{
			int const lo = chunkLow( c );
			int const hi = chunkLow( c + 1 );
			int count = 0;
			for ( int slot = lo; slot < hi; ++slot )
				if ( d_perElement[ slot ].accepted != 0 )
					++count;
			d_offsets[ c ] = count;
		} );

		mfem::forall( 1, [ = ] MFEM_HOST_DEVICE ( int )
		{
			int running = 0;
			for ( int c = 0; c < chunks; ++c )
			{
				int const count = d_offsets[ c ];
				d_offsets[ c ] = running;
				running += count;
			}
			d_offsets[ chunks ] = running;
			d_total[ 0 ] = running;
		} );

		return true;
	}

	void axisCandidatesScatterKernel(
		AxisSearchInputs const &inputs,
		mfem::Array<AxisCandidate> const &perElement,
		mfem::Array<int> const &offsets,
		mfem::Array<AxisCandidate> &out )
	{
		out.UseDevice( true );

		int const slots = 2*inputs.selectedCount;
		int const chunks = axisChunkCount( slots );
		int const room = out.Size();
		AxisCandidate const *const d_perElement = perElement.Read();
		int const *const d_offsets = offsets.Read();
		AxisCandidate *const d_out = out.Write();

		auto chunkLow = [ = ] MFEM_HOST_DEVICE ( int c )
		{
			return static_cast<int>( static_cast<long>( c )*slots/chunks );
		};

		mfem::forall( chunks, [ = ] MFEM_HOST_DEVICE ( int c )
		{
			int const lo = chunkLow( c );
			int const hi = chunkLow( c + 1 );
			int at = d_offsets[ c ];
			for ( int slot = lo; slot < hi; ++slot )
			{
				if ( d_perElement[ slot ].accepted == 0 )
					continue;
				if ( at < room )
					d_out[ at ] = d_perElement[ slot ];
				++at;
			}
		} );
	}


	bool xPointFieldKernel( XPointFieldInputs const &inputs,
	                        mfem::Array<XPointFieldValues> &out )
	{
		out.SetSize( 1 );
		out.UseDevice( true );

		XPointFieldInputs const in = inputs;
		XPointFieldValues *const d_out = out.Write();

		/*
		 * ONE THREAD, BECAUSE THERE IS ONE POINT. There is no parallelism to
		 * find here and none is wanted: this kernel exists to avoid an 830 KB
		 * download, not to go fast, and the eight numbers it produces are the
		 * whole of what the X-point border needs from the iterate.
		 */
		mfem::forall( 1, [ = ] MFEM_HOST_DEVICE ( int )
		{
			XPointFieldValues got;
			got.flux[ 0 ] = 0.0;
			got.flux[ 1 ] = 0.0;
			for ( int i = 0; i < 4; ++i )
				got.fluxJacobian[ i ] = 0.0;
			got.potentialGradient[ 0 ] = 0.0;
			got.potentialGradient[ 1 ] = 0.0;
			got.ok = 0;

			AffineTriangle const &map = in.maps[ in.element ];
			if ( !map.affine )
			{
				d_out[ 0 ] = got;
				return;
			}

			int const fluxDofs = in.flux.dofs;
			int const potentialDofs = in.potential.dofs;
			mfem::real_t local[ triangleBasisMaxDofs ];
			mfem::real_t coefficients[ 2*triangleBasisMaxDofs ];
			mfem::real_t potentialCoefficients[ triangleBasisMaxDofs ];

			// The same gather the axis kernel makes, with the element fixed:
			// blocked L2 dofs, byNODES components, and the flux NEGATED
			// because mfem::DarcyForm holds -q.
			for ( int c = 0; c < 2; ++c )
			{
				for ( int i = 0; i < fluxDofs; ++i )
					local[ i ] = -in.state[ c*in.fluxScalarDofs
					                        + in.element*fluxDofs + i ];
				triangleBasisCoefficients( in.flux, local,
				                           coefficients + c*fluxDofs );
			}
			for ( int i = 0; i < potentialDofs; ++i )
				local[ i ] = in.state[ in.potentialOffset
				                       + in.element*potentialDofs + i ];
			triangleBasisCoefficients( in.potential, local,
			                           potentialCoefficients );

			/*
			 * REFERENCE GRADIENTS, THEN ONE CHAIN RULE EACH.
			 * mfem::GridFunction::GetVectorGradient() is
			 * `grad_hat * T.InverseJacobian()` and GetGradient() is
			 * `InverseJacobian().MultTranspose( gh )`; meq::AffineTriangle
			 * stores that same inverse row major, so both contractions are
			 * the one meq::rootInTriangle() already performs and no new index
			 * convention is introduced.
			 */
			mfem::real_t value[ 2 ];
			mfem::real_t dXi[ 2 ];
			mfem::real_t dEta[ 2 ];
			for ( int c = 0; c < 2; ++c )
				evaluateTriangleBasis( in.flux.degree,
				                       coefficients + c*fluxDofs,
				                       in.xi, in.eta,
				                       value[ c ], &dXi[ c ], &dEta[ c ] );

			got.flux[ 0 ] = value[ 0 ];
			got.flux[ 1 ] = value[ 1 ];
			for ( int i = 0; i < 2; ++i )
				for ( int j = 0; j < 2; ++j )
					got.fluxJacobian[ 2*i + j ] =
						dXi[ i ]*map.inverseJacobian[ j ]
						+ dEta[ i ]*map.inverseJacobian[ 2 + j ];

			mfem::real_t psi = 0.0;
			mfem::real_t psiXi = 0.0;
			mfem::real_t psiEta = 0.0;
			evaluateTriangleBasis( in.potential.degree, potentialCoefficients,
			                       in.xi, in.eta, psi, &psiXi, &psiEta );
			for ( int j = 0; j < 2; ++j )
				got.potentialGradient[ j ] =
					psiXi*map.inverseJacobian[ j ]
					+ psiEta*map.inverseJacobian[ 2 + j ];

			got.ok = 1;
			d_out[ 0 ] = got;
		} );

		return true;
	}

	bool affineTrianglesKernel( mfem::Vector const &coordinates,
	                            mfem::Vector const &jacobians, int elements,
	                            mfem::Array<AffineTriangle> &out )
	{
		auto const *d_x = coordinates.Read();
		auto const *d_j = jacobians.Read();
		auto *d_out = out.Write();
		mfem::forall( elements, [ = ] MFEM_HOST_DEVICE ( int e )
		{
			// NQ = 1: X[ i + 2 e ], J[ i + 2 ( j + 2 e ) ] = d x_i / d xi_j.
			AffineTriangle map{};
			map.origin[ 0 ] = d_x[ 0 + 2*e ];
			map.origin[ 1 ] = d_x[ 1 + 2*e ];
			map.jacobian[ 0 ] = d_j[ 0 + 2*( 0 + 2*e ) ];
			map.jacobian[ 1 ] = d_j[ 0 + 2*( 1 + 2*e ) ];
			map.jacobian[ 2 ] = d_j[ 1 + 2*( 0 + 2*e ) ];
			map.jacobian[ 3 ] = d_j[ 1 + 2*( 1 + 2*e ) ];
			mfem::real_t const determinant =
				map.jacobian[ 0 ]*map.jacobian[ 3 ] - map.jacobian[ 1 ]*map.jacobian[ 2 ];
			map.affine = determinant != mfem::real_t( 0 );
			if ( map.affine )
			{
				mfem::real_t const inverse = mfem::real_t( 1 )/determinant;
				map.determinant = determinant;
				map.inverseJacobian[ 0 ] =  map.jacobian[ 3 ]*inverse;
				map.inverseJacobian[ 1 ] = -map.jacobian[ 1 ]*inverse;
				map.inverseJacobian[ 2 ] = -map.jacobian[ 2 ]*inverse;
				map.inverseJacobian[ 3 ] =  map.jacobian[ 0 ]*inverse;
			}
			d_out[ e ] = map;
		} );
		return true;
	}

#else

	bool axisFluxScaleKernel( mfem::Vector const &, int, mfem::real_t,
	                          mfem::Vector & )
	{
		// No device translation unit in this build; locateAxisPoint() keeps
		// meq::CriticalPointFinder, which is the whole implementation.
		return false;
	}

	bool axisCandidatesCountKernel( AxisSearchInputs const &,
	                                mfem::Array<AxisCandidate> &,
	                                mfem::Array<int> &,
	                                mfem::Array<int> & )
	{
		return false;
	}

	void axisCandidatesScatterKernel( AxisSearchInputs const &,
	                                  mfem::Array<AxisCandidate> const &,
	                                  mfem::Array<int> const &,
	                                  mfem::Array<AxisCandidate> & )
	{
	}

	bool xPointFieldKernel( XPointFieldInputs const &,
	                        mfem::Array<XPointFieldValues> & )
	{
		return false;
	}

	bool affineTrianglesKernel( mfem::Vector const &, mfem::Vector const &,
	                            int, mfem::Array<AffineTriangle> & )
	{
		return false;
	}

#endif

}
