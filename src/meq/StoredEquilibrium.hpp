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

			/**
			 * THE SOLVER'S WHOLE ITERATE, ( q, psi, psihat ) as
			 * GradShafranovSolver::state() lays it out -- the flux block being
			 * DarcyForm's -q, which flux() negates -- assembled from the file's
			 * flux, potential and trace coefficients, bit for bit what the run
			 * held. Empty when the file carries no trace, from before it did.
			 *
			 * VALID ONLY ON THE MESH IT WAS SOLVED ON, NUMBERED AS THAT SOLVE
			 * NUMBERED IT: the trace lives on the faces, and a mesh read back
			 * from text may number its faces differently. So a caller hands it
			 * to setInitialState() only for a solve whose own mesh prints as
			 * meshText() does, and that solver then numbers the faces as the
			 * one that wrote this did.
			 */
			bool hasTrace() const { return !trace.empty(); }
			void wholeState( mfem::Vector &state ) const;

			/// The exterior coupling's coefficients, `exterior_a2` ... in order,
			/// empty without one: setInitialState()'s second argument.
			std::vector<double> const &exteriorCoefficients() const { return exterior; }

			/**
			 * THE BORDER'S OTHER UNKNOWNS, AS THE RUN SOLVED THEM: what the
			 * adaptive loop hands one cycle's solver from the last, read from
			 * the file so a restart starts there rather than at the TOML's
			 * seeds. Each is present only where the run had it.
			 *
			 *   xPoint      the null, when it was an unknown and still located
			 *   profileScale  the I_p border's profile scale
			 *   inverseCurrents  the free circuits' currents, in the run's
			 *                order of `inverse_circuits`
			 *   gSquaredCoefficients  a q-driven run's g^2 polynomial
			 */
			bool hasXPoint() const { return xPointGiven; }
			double xPointR() const { return xPointRadius; }
			double xPointZ() const { return xPointHeight; }
			/// The located magnetic axis the run reported, `axis_r`/`axis_z`,
			/// when it found one: what seeds a restart's first plasma fill.
			bool hasAxisPoint() const { return axisPointGiven; }
			double axisPointR() const { return axisPointRadius; }
			double axisPointZ() const { return axisPointHeight; }

			/// psi_bnd, the boundary label of a limited or diverted run.
			bool hasPsiBoundary() const { return boundaryGiven; }
			double psiBoundary() const { return boundary; }
			bool hasProfileScale() const { return scaleGiven; }
			double profileScale() const { return scale; }
			std::vector<double> const &inverseCurrents() const { return currents; }
			/// A q-driven run's g^2( Psi ) coefficients, its outer unknowns.
			std::vector<double> const &gSquaredCoefficients() const { return gSquared; }

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
			std::vector<double> trace, exterior;
			bool xPointGiven = false, scaleGiven = false, boundaryGiven = false;
			bool axisPointGiven = false;
			double axisPointRadius = 0.0, axisPointHeight = 0.0;
			double boundary = 0.0;
			double xPointRadius = 0.0, xPointHeight = 0.0, scale = 1.0;
			std::vector<double> currents, gSquared;
			bool axisGiven = false;
			double axis = 0.0;
			std::string contentText;
			int degree = 0;

			std::vector<double> region, plasmaFraction, conductorMaterial, conductorFraction;
	};
}

#endif // MEQ_STOREDEQUILIBRIUM_HPP
