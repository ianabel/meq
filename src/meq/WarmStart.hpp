#ifndef MEQ_WARMSTART_HPP
#define MEQ_WARMSTART_HPP

#include "mfem.hpp"

/*
 * Carrying a solved potential onto a different mesh, at the order it was
 * computed with.
 *
 * DRIVER-PLAN.md section 4 names three restart routes and they are not the same
 * problem:
 *
 *   EXACT          same mesh, same degree. Read the stored GridFunction and hand
 *                  it straight to setInitialGuess(). Bitwise resumable, and it
 *                  needs nothing in this file.
 *   FULL ORDER     a different mesh or degree, same code. What this file is.
 *   INTERCHANGE    another code entirely: psi( R, Z ) on a structured grid, read
 *                  back by bilinear interpolation. That is Output.hpp's NetCDF
 *                  file, and it is genuinely the right format for foreign input
 *                  -- no MFEM, no mesh, no agreement about element types.
 *
 * THE THIRD IS NOT A SUBSTITUTE FOR THE SECOND, WHICH IS THE WHOLE ARGUMENT FOR
 * TAKING THE GSLIB DEPENDENCY. Bilinear interpolation on a structured grid is
 * SECOND ORDER however fine the grid is, so restarting a k = 3 solve through it
 * throws away almost everything the previous solve computed. It still converges
 * -- a guess is only a guess -- but "warm" is then doing less work than it
 * looks. Staying inside the finite element representation carries the full
 * k+1. WarmStartConvergence.cpp measures the difference rather than asserting
 * it, which is what justifies the dependency.
 *
 * WHY THIS IS A BATCHED TRANSFER AND NOT AN mfem::Coefficient. The obvious
 * design is a Coefficient that calls FindPoints per evaluation, and it would be
 * unusable: gslib's search amortises a hash and a crystal-router exchange over
 * the whole query, so one point at a time pays all of that per point. So the
 * transfer runs once, over every node of the target space, and the result is an
 * ordinary GridFunction on the target's own mesh -- which
 * GradShafranovSolver::setInitialGuess( GridFunction const & ) already accepts,
 * because a GridFunctionCoefficient over it evaluates on the mesh the solver is
 * actually using.
 *
 * POINTS THE SOURCE MESH DOES NOT COVER ARE COUNTED AND REPORTED. A refined
 * domain can reach beyond the stored one, and on the extension path D_h moves
 * outwards every adaptive cycle by construction. The guess is only a guess, so
 * the fallback need not be clever -- but a restart that silently found no data
 * for most of the domain is a restart that quietly became a cold start, and the
 * iteration count will not obviously say so.
 */

namespace meq
{

	/**
	 * Evaluates a stored field at the nodes of another finite element space.
	 *
	 * Setup is done once per source mesh and is the expensive part, so a caller
	 * transferring several fields from the same mesh -- psi and the two flux
	 * components, say -- should keep one of these rather than build one each
	 * time.
	 */
	class FieldTransfer
	{
		public:
			/// @param sourceMeshIn the mesh the stored field lives on. Borrowed,
			///                     and must outlive this object. Not const because
			///                     FindPointsGSLIB::Setup() is not.
			explicit FieldTransfer( mfem::Mesh &sourceMeshIn );

			~FieldTransfer();

			FieldTransfer( FieldTransfer const & ) = delete;
			FieldTransfer &operator=( FieldTransfer const & ) = delete;

			/**
			 * Fill @a target with @a source evaluated at @a target's nodes.
			 *
			 * @param source    a field on the mesh this was constructed with.
			 * @param fallback  evaluated at any node the source mesh does not
			 *                  cover. The Dirichlet datum is the sensible choice:
			 *                  it is what a cold start would have had there.
			 * @param target    filled. Its space decides where the source is
			 *                  sampled, so it must already be sized.
			 * @return          how many nodes fell outside the source mesh.
			 *
			 * The nodal reading is exact for the Gauss-Legendre L2 spaces MEQ
			 * uses, where a dof IS a point value. It is an interpolation rather
			 * than a projection for anything else, which for a starting guess is
			 * the right trade -- one pass, no mass matrix.
			 */
			int transfer( mfem::GridFunction const &source,
			              mfem::Coefficient &fallback,
			              mfem::GridFunction &target );

			/// Nodes outside the source mesh on the last transfer().
			int missed() const;

			/// Nodes on the last transfer(), so a caller can report a fraction
			/// rather than a bare count.
			int queried() const;

			/// The largest distance gslib had to move a point to land it in an
			/// element, on the last transfer(). Large values mean the two meshes
			/// disagree about where the domain is, which a miss count of zero
			/// will not tell you on its own.
			double worstDistance() const;

		private:
			mfem::Mesh &sourceMesh;
			mfem::FindPointsGSLIB finder;
			int misses;
			int points;
			double worstDist;
	};

}

#endif // MEQ_WARMSTART_HPP
