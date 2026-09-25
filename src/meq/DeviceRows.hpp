#ifndef MEQ_DEVICEROWS_HPP
#define MEQ_DEVICEROWS_HPP

#include "mfem.hpp"

#include <vector>

/**
 * @file
 *
 * meq::CompressedRows' DOT PRODUCTS, ON THE DEVICE, IN ONE LAUNCH.
 *
 * WHAT THIS IS FOR, AND IT IS A TRANSFER RATHER THAN A FLOP COUNT.
 * `CompressedRows::gather()` opens with `v.HostRead()`. That is the whole
 * point of this file: to read a few hundred entries of a support it brings the
 * ENTIRE vector to the host. On `machine-f-diiid` at `k = 2` the pair of
 * solution blocks is 87,264 doubles, 698 KB, and it is `O( elements )`. The
 * bordered step does this **21.1 times per Newton step** -- 15.8 of them the
 * exterior transmission contraction and the rest the border's own row dots --
 * measured over a twelve-step run. DEVICE-PLAN.md section 3.4a is the
 * inventory and this is the first item off it.
 *
 * **THE KERNEL IS ONE THREAD PER ROW AND THAT IS DELIBERATE, NOT A FIRST
 * DRAFT.** meq::CompressedRows documents a bit-exactness contract in terms --
 * *"summing the surviving terms in ASCENDING INDEX ORDER reproduces the dense
 * loop's every partial sum"* -- and a tree reduction inside a row would
 * reassociate it and give that up for a number that differs in its last bits.
 * MEQ cannot afford that: `--device cuda` reproducing `--device cpu` to every
 * printed digit is the instrument that root-caused M-79, M-130 and M-187, and
 * DEVICE-PLAN.md section 3.5 records that cuDSS will take it away later on.
 * Giving it up HERE, voluntarily, for a reduction over at most a few thousand
 * terms, would be trading the only cheap diagnostic MEQ has for occupancy on a
 * launch that exists to avoid a 698 KB copy. **The rows are the parallelism**;
 * there are `N + 4` of them and they are independent.
 *
 * SO THE ANSWER IS BIT-IDENTICAL TO THE HOST'S, AND
 * `tests/unit/DeviceRowsTests.cpp` ASSERTS EXACTLY THAT -- `==` on the double,
 * not a tolerance. If that test is ever relaxed to a tolerance, this file has
 * stopped being what it claims to be.
 *
 * IT IS NOT WIRED INTO THE SOLVE YET, AND THAT IS ALSO DELIBERATE.
 * A device reduction over the state is CORRECT as soon as it is written and
 * only PAYS once the state is already device-resident. Today the trace solve
 * and `NPCRecover` leave it host-valid, so calling this from the bordered step
 * would push the state TO the device and add a transfer rather than remove
 * one. It is built, compiled by the device compiler and pinned now because it
 * can be -- against a host reference, with the digit-for-digit instrument
 * still available -- and switched on when DEVICE-PLAN.md's item 1 lands. That
 * is the same arrangement meq::SourceView is in, for the same reason.
 */

namespace meq
{

	/**
	 * A device-resident mirror of a meq::CompressedRows, and its reduction.
	 *
	 * The mirror is uploaded when the rows are built -- once per support sweep
	 * on the free-boundary path, not once per Newton step -- so the per-step
	 * cost is the launch alone.
	 */
	class DeviceRows
	{
		public:
			/**
			 * Mirror a compressed row set into device memory.
			 *
			 * @param support   the shared support, ascending and duplicate
			 *                  free, exactly as meq::CompressedRows holds it.
			 * @param values    `rowCount * support.size()` values, row major,
			 *                  which is meq::CompressedRows' own layout.
			 * @param rowCount  how many rows.
			 *
			 * @throws std::invalid_argument if @a values is not
			 *         `rowCount * support.size()` long, which is the one error
			 *         that would otherwise read off the end of a device buffer
			 *         with nothing to say about it.
			 */
			void upload( std::vector<int> const &support,
			             std::vector<double> const &values, int rowCount );

			/// Forget everything and release the device memory.
			void clear();

			/// Whether upload() has been called with a non-empty row set.
			bool valid() const;

			int rowCount() const;

			/// The support's length, which is the row width.
			int width() const;

			/**
			 * Every row's dot product against @a v, in one launch.
			 *
			 * `out( r ) = sum_k values[ r*width + k ] * v( support[ k ] )`,
			 * summed ascending in `k`, which is **bit-identical** to
			 * `meq::CompressedRows::dot( r, v, scratch )`.
			 *
			 * @param out resized to rowCount() and left **DEVICE-VALID**. The
			 *        caller funnels if it wants these on the host, which is
			 *        MEQ's convention everywhere else and is deliberate here:
			 *        the bordered step's arithmetic around them is host
			 *        arithmetic *today*, and a kernel that forced the read
			 *        would make it impossible for a resident border ever to
			 *        keep the result. Whether the host needs them is the
			 *        caller's question.
			 *
			 * @return false, having touched nothing, when no mfem::Device is
			 *         configured or the mirror is empty. **A refusal is the
			 *         safe answer and the caller keeps the host path** -- the
			 *         same contract meq::makeSourceView() has, and for the same
			 *         reason: there is no approximate mode worth having.
			 */
			bool dotAll( mfem::Vector const &v, mfem::Vector &out ) const;

			/**
			 * Every row against every column, in ONE launch.
			 *
			 * `out( r, c ) = sum_k values[ r*width + k ]
			 *                * columns[ c*length + support[ k ] ]`,
			 * summed ascending in `k`, so each entry is **bit-identical** to
			 * `meq::CompressedRows::dot()` exactly as dotAll()'s are.
			 *
			 * **THIS IS THE BORDERED STEP'S `( N + 4 )²` FILL.** The host does
			 * `dense( i, j ) = cornerEntry( i, j ) - rowDot( i, columnZ( j ) )`
			 * with one gather per column and then a dot per pair; this is the
			 * dot per pair, with the gather implicit because the columns are
			 * already where the support indexes them.
			 *
			 * @param columns `columnCount` vectors laid out CONTIGUOUSLY, column
			 *        `c` at offset `c*length`. That is not an imposition on the
			 *        caller: the bordered step already produces its columns as
			 *        one batch through `DarcyNPCSolver::ArrayMult()` -- which is
			 *        what M-98's 1.34× is -- so they arrive packed and a
			 *        per-column `mfem::Vector` would be the copy, not this.
			 * @param length the length of one column, which must be at least
			 *        one past the support's largest index.
			 * @param out resized to `rowCount() x columnCount` and left
			 *        **DEVICE-VALID**, for the reason dotAll() gives. Bringing
			 *        it over is 196 doubles on the shipped machine case against
			 *        the state's 698 KB -- a large win and still a per-step
			 *        SYNCHRONISATION, so it is the caller's to choose rather
			 *        than the kernel's to impose.
			 *
			 * @return false, having touched nothing, on the same refusals
			 *         dotAll() makes, plus a `columns` that is not
			 *         `columnCount*length` long.
			 */
			bool dotColumns( mfem::Vector const &columns, int columnCount,
			                 int length, mfem::DenseMatrix &out ) const;

		private:
			mfem::Array<int> indices;
			mfem::Vector values;
			int rows = 0;
			int columns = 0;

			/*
			 * THE LARGEST SUPPORT INDEX, KEPT AS A PLAIN int RATHER THAN READ
			 * BACK OFF `indices`.
			 *
			 * upload() leaves the mirror DEVICE-VALID, which is the whole
			 * point of it, and mfem::Array<int>::operator[] is as raw an
			 * accessor as mfem::Vector::operator() -- meq::CompressedRows' own
			 * comment about rowDot() says exactly that. So the bounds check in
			 * dotAll() may not look at the array it is checking. It read
			 * `indices[ columns - 1 ]` in the first draft of this class and
			 * mfem::Device( "debug" ) faulted on it at the first run, with
			 * "invalid permissions" on the protected host page: the same defect
			 * class as M-79 and M-187, committed while writing the file whose
			 * header argues about them.
			 */
			int largestIndex = -1;
	};


	/**
	 * The launch itself, split out so that exactly one translation unit in MEQ
	 * is compiled by the device compiler per feature.
	 *
	 * `mfem::forall` `static_assert`s under `g++`, so a kernel cannot live in
	 * an ordinary MEQ source; see PointLocatorDevice.cpp. Compiles to a stub
	 * returning false where no device language was enabled.
	 *
	 * @a indices and @a values are the mirror's own containers, taken by
	 * reference rather than as raw pointers so that `Read()` decides where they
	 * live -- a raw pointer here would be exactly the class of defect M-187
	 * records, a host address believed about a buffer whose live copy is
	 * elsewhere.
	 */
	bool dotCompressedRowsOnDevice( mfem::Vector const &values,
	                                mfem::Array<int> const &indices,
	                                int rowCount, int width,
	                                mfem::Vector const &v,
	                                mfem::Vector &out );

	/// dotColumns()'s launch; see dotCompressedRowsOnDevice() for why the
	/// containers are taken by reference rather than as raw pointers.
	bool dotCompressedColumnsOnDevice( mfem::Vector const &values,
	                                   mfem::Array<int> const &indices,
	                                   int rowCount, int width,
	                                   mfem::Vector const &columns,
	                                   int columnCount, int length,
	                                   mfem::DenseMatrix &out );

}

#endif // MEQ_DEVICEROWS_HPP
