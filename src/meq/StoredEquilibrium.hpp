#ifndef MEQ_STOREDEQUILIBRIUM_HPP
#define MEQ_STOREDEQUILIBRIUM_HPP

/*
 * AN EQUILIBRIUM READ BACK OUT OF THE `.nc` MEQ WROTE, AS FINITE ELEMENT
 * FIELDS ON THE MESH THEY WERE SOLVED ON.
 *
 * The `.nc` is MEQ's default and only routine output, and it carries everything
 * the exact formats do: the solve's mesh as MFEM's own text (the `mesh`
 * variable, the bytes `<stem>.mesh` would hold), every coefficient of psi_h, of
 * the flux q and of the post-processed psi*, with the space each lives in, and
 * the boundary curve when the run had one. So it is the RESTART --
 * `[initialguess] Type = "gridfunction"`, `File = "<stem>.nc"` -- and it is what
 * `meq-convert` makes the GLVis and ParaView files from.
 *
 * THE MESH IS READ WITH NEITHER REFINEMENT MARKING NOR ORIENTATION FIXING.
 * Both are allowed to reorder a triangle's vertices, and an L2 coefficient is
 * indexed by its element's LOCAL basis, which follows the vertex order -- so a
 * mesh reordered on reading would carry the coefficients in the wrong places.
 * The text is exactly what the solve held, and reading it as written
 * reproduces it whatever made the mesh. Measured on the shipped curved and
 * conductor cases, MFEM's default reading reorders nothing there (with it,
 * OutputConversion's byte comparison still passes), so this is a guarantee for
 * a mesh not already in MFEM's marked order rather than a defect seen.
 *
 * Needs MFEM and netcdf-cxx4; without the second the constructor throws.
 */

#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "BoundaryShape.hpp"

namespace meq
{
	class StoredEquilibrium
	{
		public:
			/// @throws std::runtime_error if @a path cannot be read, carries no
			///         embedded mesh -- a file from before the mesh was embedded
			///         -- or holds coefficients that do not fit their space.
			explicit StoredEquilibrium( std::string const &path );

			StoredEquilibrium( StoredEquilibrium const & ) = delete;
			StoredEquilibrium &operator=( StoredEquilibrium const & ) = delete;

			/// The solve's mesh, as read; meq-convert bends it for the VTK.
			mfem::Mesh &mesh() { return *meshPointer; }
			/// The text it was read from, byte for byte what `<stem>.mesh` holds.
			std::string const &meshText() const { return text; }

			/// psi_h in P_k: the restart field, and `<stem>_psi.gf`.
			mfem::GridFunction const &potential() const { return *potentialField; }
			/// q, vdim 2: `<stem>_grad_psi.gf`.
			mfem::GridFunction const &flux() const { return *fluxField; }
			/// Its space, for a field built beside it -- B, which the VTK draws.
			mfem::FiniteElementSpace &fluxSpaceOf() { return *fluxSpace; }
			/// psi* in P_(k+1): `<stem>_psistar.gf`, and what the VTK draws.
			mfem::GridFunction const &postProcessed() const { return *postField; }

			/// The boundary curve, or null when the run had no [boundary.shape].
			BoundaryShape const *shape() const { return curve.get(); }

			/// The `content` attribute: whether the coefficients are psi or the
			/// remainder psi - psi_c of a subtracted run.
			bool holdsTotal() const { return total; }
			std::string const &content() const { return contentText; }

			/// The converged axis flux, psi_ax, which a normalised run records
			/// as its `psi_axis` attribute and a restart starts its bordered
			/// unknown from; hasPsiAxis() is false for a run without one.
			bool hasPsiAxis() const { return axisGiven; }
			double psiAxis() const { return axis; }

			/// [discretisation] PolynomialDegree, k.
			int polynomialDegree() const { return degree; }

			/// The element record's four drawn columns, one per element, or
			/// empty when the file has no record.
			std::vector<double> const &regionColumn() const { return region; }
			std::vector<double> const &plasmaFractionColumn() const { return plasmaFraction; }
			std::vector<double> const &conductorMaterialColumn() const { return conductorMaterial; }
			std::vector<double> const &conductorFractionColumn() const { return conductorFraction; }
			bool hasElementColumns() const { return !region.empty(); }

		private:
			std::string text;
			std::unique_ptr<mfem::Mesh> meshPointer;

			std::unique_ptr<mfem::FiniteElementCollection> potentialCollection;
			std::unique_ptr<mfem::FiniteElementCollection> fluxCollection;
			std::unique_ptr<mfem::FiniteElementCollection> postCollection;
			std::unique_ptr<mfem::FiniteElementSpace> potentialSpace;
			std::unique_ptr<mfem::FiniteElementSpace> fluxSpace;
			std::unique_ptr<mfem::FiniteElementSpace> postSpace;
			std::unique_ptr<mfem::GridFunction> potentialField;
			std::unique_ptr<mfem::GridFunction> fluxField;
			std::unique_ptr<mfem::GridFunction> postField;

			std::unique_ptr<BoundaryShape> curve;
			bool total = true;
			bool axisGiven = false;
			double axis = 0.0;
			std::string contentText;
			int degree = 0;

			std::vector<double> region, plasmaFraction, conductorMaterial, conductorFraction;
	};
}

#endif // MEQ_STOREDEQUILIBRIUM_HPP
