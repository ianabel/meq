#ifndef MEQ_CONTINUOUS_FIELD_HPP
#define MEQ_CONTINUOUS_FIELD_HPP

/*
 * THE CONTINUOUS PAIR THE FLUX-SURFACE EXTRACTION ROOTS AND WEIGHTS BY --
 * CLAUDE_INVERSION.md, IN-0; MEASUREMENTS.md M-364 to M-367.
 *
 * psi* and q* are element-wise polynomials, broken across every face, so a
 * flux-surface node whose ray root slides across a face lands on the
 * neighbour's polynomial and the geometry JUMPS -- by the DG jump, which is
 * inside the discretisation's accuracy but which a caller differencing the
 * geometry at a relative step of 1e-8 divides by that step (M-364). The
 * fix is a pair of fields that are continuous, built from the solved state
 * with no side chosen anywhere:
 *
 *   psi~ in H1 of degree k+2, the MOMENT interpolant: vertex values the mean
 *   of psi* over the elements meeting there; on each interior edge the P_k
 *   moments of psihat, the HDG solve's single-valued face unknown; on each
 *   boundary edge psi*'s own trace; in each cell psi*'s moments against
 *   P_(k-1) (triangles) or Q_k (quadrilaterals). Degree k+2 because that is
 *   the degree whose edge -- two vertex values and k+1 more -- can carry ALL
 *   of psihat's k+1 moments: psi~'s trace then projects onto psihat on every
 *   interior edge. On a boundary edge there is one side, so nothing to choose,
 *   and on Gamma_h the trace unknown is a pinned zero rather than psihat
 *   (GradShafranovSolver::setExtension()).
 *
 *   q~ in [H1]^2 of q*'s degree k+1, the same interpolant: vertex means of
 *   q*, edge moments of the two sides' mean trace, cell moments of q*. The
 *   weight 1/( R | q | ) needs a continuous flux as much as the root needs a
 *   continuous potential, and grad psi~ is not one: its normal component jumps.
 *
 * The interpolant is unisolvent, reproduces polynomials of its degree, and is
 * linear in ( psi*, psihat, q* ) -- so in the HDG state, Reconstruct() being
 * linear, which is what the extraction's derivative needs. Measured
 * (M-365): a central difference at 1e-8 straddling a crossing is off by
 * the noise rather than by up to 700 times the derivative, and the geometry
 * moves by less than its own discretisation error.
 *
 * EVERY DOF IS WRITTEN BY ONE ENTITY -- a vertex, an edge or a cell -- in that
 * order, each reading only what the order before it wrote; the vertex means
 * are summed serially in element order. So the threaded build is the serial
 * one bit for bit, which FluxSurfaceConvergence asserts.
 *
 * AND THE MAP IS KEPT. Every coefficient of it belongs to the spaces and the
 * mesh, so ContinuousInterpolant does the quadratures and the moment solves
 * once and applies the result as matrices, on the device under one: the
 * solver keeps one, and buildContinuousFields() builds one per call
 * (MEASUREMENTS.md M-385).
 */

#include <functional>
#include <memory>

#include "mfem.hpp"

namespace meq
{
	/// A field's value in an element at a reference point.
	using CellValue = std::function<double( int, mfem::IntegrationPoint const & )>;

	/// A field's value on a face at the face's own reference point; the face
	/// transformation has had SetAllIntPoints() called at it, so both sides'
	/// element points are set.
	using FaceValue = std::function<double( int, mfem::FaceElementTransformations &,
	                                        mfem::IntegrationPoint const & )>;

	/**
	 * Fill @a out, a scalar H1 GridFunction of uniform degree p on a mesh of
	 * triangles or quadrilaterals, with the moment interpolant: vertex values
	 * the mean of @a cell at that vertex over the elements meeting there; edge
	 * dofs from @a face's moments against P_(p-2); cell dofs from @a cell's
	 * moments against P_(p-3) (triangles) or Q_(p-2) (quadrilaterals). Both
	 * callables are called from several threads at once and must be reentrant.
	 * @throws std::invalid_argument on any other element geometry.
	 */
	void momentInterpolant( mfem::GridFunction &out, CellValue const &cell, FaceValue const &face );

	/// The pair, and what keeps their spaces alive: the spaces are shared
	/// with the interpolant that made them, so that applying it again builds
	/// no space (`spaces` is declared first and so outlives the fields).
	struct ContinuousFields
	{
		std::shared_ptr<void const> spaces;
		std::unique_ptr<mfem::GridFunction> potential;
		std::unique_ptr<mfem::GridFunction> flux;
	};

	/**
	 * psi~ and q~ from psi* ( @a potential, an L2 field of degree k+1 ), q*
	 * ( @a flux, L2 of vdim 2 ) and psihat ( @a trace, the face space ), on
	 * their mesh; each GridFunction owns its space. The inputs must have VALUE
	 * map type, which is what lets them be evaluated from several threads.
	 * One ContinuousInterpolant, built and applied once: a caller applying the
	 * map again on the same spaces keeps the interpolant instead.
	 * @throws std::invalid_argument otherwise.
	 */
	ContinuousFields buildContinuousFields( mfem::GridFunction const &potential,
	                                        mfem::GridFunction const &flux,
	                                        mfem::GridFunction const &trace );
	/**
	 * A LINEAR MAP APPLIED BY BLOCKS, which is what the moment interpolant is
	 * once its quadratures are done: each block writes its own outputs as a
	 * small dense matrix times a gathered column -- inputs first, then output
	 * entries an earlier application wrote -- and no two blocks write the same
	 * entry. Index `i` of an input reads `first[ i ]` below the caller's split
	 * and `second[ i - split ]` above it. `matrix` is row-major per block.
	 */
	struct LinearBlocks
	{
		int count = 0;
		mfem::Array<int> outOffset, inOffset, knownOffset, matrixOffset;
		mfem::Array<int> outIndex, inIndex, knownIndex;
		mfem::Vector matrix;
	};

	/// The host route: every block, in order, threaded where the host may
	/// thread (no two blocks write one entry, so the order is immaterial and
	/// the result the serial one to the bit).
	void applyLinearBlocks( LinearBlocks const &blocks, mfem::Vector const &first,
	                        mfem::Vector const *second, int split, mfem::Vector &out, int outShift );

	/// The same loop as a kernel (ContinuousFieldDevice.cpp), false where it
	/// cannot run: no Device, or a build whose device files are host stubs.
	/// It reads its inputs and writes @a out on the device, and leaves both
	/// there.
	bool applyLinearBlocksOnDevice( LinearBlocks const &blocks, mfem::Vector const &first,
	                                mfem::Vector const *second, int split, mfem::Vector &out,
	                                int outShift );

	/**
	 * THE MOMENT INTERPOLANT FOR ONE SET OF SPACES, ITS QUADRATURES DONE ONCE.
	 *
	 * psi~ and q~ are linear in ( psi*, psihat, q* ), and every coefficient of
	 * that map is a property of the spaces and the mesh: so the constructor
	 * does each entity's quadrature once, solves its moment system once, and
	 * keeps the result as three stages of LinearBlocks -- vertex means, edge
	 * moments, cell moments -- for psi~ and for each component of q~. apply()
	 * is then three gathered matrix products per field, on the device under
	 * one (applyLinearBlocksOnDevice()) and on the host otherwise, the same
	 * loop either way, so `--device debug` reproduces the host to the bit.
	 * MEASUREMENTS.md M-385.
	 *
	 * The vertex value is the mean of psi*'s corner values over the elements
	 * meeting there, folded into one row per vertex whose terms run in element
	 * order; the edge and cell systems are buildContinuousFields()'s, solved in
	 * double. The spaces must be the ones buildContinuousFields() accepts.
	 */
	class ContinuousInterpolant
	{
		public:
			/// Where apply() ran. Automatic takes the device under any Device
			/// that allows one, and the host otherwise.
			enum class Route
			{
				Automatic,
				Host,
				Device
			};

			/// @throws std::invalid_argument as buildContinuousFields() does.
			ContinuousInterpolant( mfem::FiniteElementSpace const &potential,
			                       mfem::FiniteElementSpace const &flux,
			                       mfem::FiniteElementSpace const &trace );

			/// psi~ and q~ of fields on the spaces this was built for, on output
			/// spaces this interpolant shares with every result. Under a Device
			/// the result is handed back to the host for its raw readers (the
			/// tracer).
			ContinuousFields apply( mfem::GridFunction const &potential,
			                        mfem::GridFunction const &flux,
			                        mfem::GridFunction const &trace,
			                        Route route = Route::Automatic ) const;

			/// Whether this was built for exactly these spaces, unchanged since.
			bool builtFor( mfem::FiniteElementSpace const &potential,
			               mfem::FiniteElementSpace const &flux,
			               mfem::FiniteElementSpace const &trace ) const;

			/// Whether the last apply() ran its stages on the device.
			bool lastRanOnDevice() const { return ranOnDevice; }

		private:
			mfem::FiniteElementSpace const *potentialSpace;
			mfem::FiniteElementSpace const *fluxSpace;
			mfem::FiniteElementSpace const *traceSpace;
			long potentialSequence, fluxSequence, traceSequence, meshSequence;
			int potentialOrder, fluxOrder;
			int potentialSize, fluxSize, traceSize, elementCount;
			/// [ stage ] for psi~; [ component ][ stage ] for q~.
			LinearBlocks potentialBlocks[ 3 ];
			LinearBlocks fluxBlocks[ 2 ][ 3 ];
			mutable bool ranOnDevice = false;
			/// The output collections and spaces, shared with every result.
			struct OutputSpaces;
			std::shared_ptr<OutputSpaces> output;
	};

}

#endif // MEQ_CONTINUOUS_FIELD_HPP
