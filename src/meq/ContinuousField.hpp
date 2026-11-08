#ifndef MEQ_CONTINUOUS_FIELD_HPP
#define MEQ_CONTINUOUS_FIELD_HPP

/*
 * THE CONTINUOUS PAIR THE FLUX-SURFACE EXTRACTION ROOTS AND WEIGHTS BY --
 * MANTA-EQUILIBRIUM-PLAN.md section 4.5, M-SHAPE-1 and M-SHAPE-2.
 *
 * psi* and q* are element-wise polynomials, broken across every face, so a
 * flux-surface node whose ray root slides across a face lands on the
 * neighbour's polynomial and the geometry JUMPS -- by the DG jump, which is
 * inside the discretisation's accuracy but which a caller differencing the
 * geometry at a relative step of 1e-8 divides by that step (M-SHAPE-1). The
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
 * (M-SHAPE-2): a central difference at 1e-8 straddling a crossing is off by
 * the noise rather than by up to 700 times the derivative, and the geometry
 * moves by less than its own discretisation error.
 *
 * EVERY DOF IS WRITTEN BY ONE ENTITY -- a vertex, an edge or a cell -- in that
 * order, each reading only what the order before it wrote; the vertex means
 * are summed serially in element order. So the threaded build is the serial
 * one bit for bit, which FluxSurfaceConvergence asserts.
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

	struct ContinuousFields
	{
		std::unique_ptr<mfem::GridFunction> potential;
		std::unique_ptr<mfem::GridFunction> flux;
	};

	/**
	 * psi~ and q~ from psi* ( @a potential, an L2 field of degree k+1 ), q*
	 * ( @a flux, L2 of vdim 2 ) and psihat ( @a trace, the face space ), on
	 * their mesh; each GridFunction owns its space. The inputs must have VALUE
	 * map type, which is what lets them be evaluated from several threads.
	 * @throws std::invalid_argument otherwise.
	 */
	ContinuousFields buildContinuousFields( mfem::GridFunction const &potential,
	                                        mfem::GridFunction const &flux,
	                                        mfem::GridFunction const &trace );
}

#endif // MEQ_CONTINUOUS_FIELD_HPP
