#ifndef MEQ_OUTPUT_HPP
#define MEQ_OUTPUT_HPP

/*
 * Writing a solved equilibrium out.
 *
 * Two audiences and two formats, and they are not interchangeable.
 *
 * MFEM's own .mesh and .gf carry the discrete solution exactly: the same spaces,
 * the same degree, every coefficient. They are what GLVis reads, what an exact
 * restart reads, and what nothing outside MFEM reads.
 *
 * The NetCDF file carries psi and both components of B on a uniform ( R, Z )
 * grid. It is lossy by construction -- a k+1 field sampled onto a rectangle --
 * and it is the format every downstream tool actually wants. It is also the
 * INTERCHANGE format of DRIVER-PLAN section 4: a structured grid interpolates
 * back in O( 1 ) per point with no mesh search, so a foreign code needs to
 * produce nothing but psi on a rectangle to warm-start meq.
 *
 * WHAT IS NOT WRITTEN, and why the file says so rather than omitting it
 * silently: psi* when the solve did not produce one. The post-processed
 * potential is only defined after postProcess(), which is refused on paths where
 * MFEM's reconstruction cannot be trusted -- see GradShafranovSolver. A file
 * that quietly lacks a field is worse than one that says the field was not
 * computed, so the writer records its absence as an attribute.
 */

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "Sampler.hpp"

namespace meq
{
	/// True when meq was built with netcdf-cxx4. Everything in
	/// NetCDFWriter throws without it.
	bool hasNetCDF();

	/// The mesh and the grid functions, in MFEM's own formats. Exact, and
	/// readable only by MFEM.
	///
	/// @param stem  the path stem: "<stem>.mesh", "<stem>_psi.gf" and
	///              "<stem>_grad_psi.gf" are written.
	/// @param flux  q in meq's sign convention -- what
	///              GradShafranovSolver::flux() returns, not the raw block.
	void writeMfem( std::string const &stem, mfem::Mesh &mesh,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux );

	/**
	 * The same discrete solution as VTK, for ParaView and VisIt.
	 *
	 * Written through mfem::ParaViewDataCollection, which produces a
	 * directory rather than a file:
	 *
	 *     <stem>/<name>.pvd              <-- THE FILE TO OPEN
	 *     <stem>/Cycle000000/data.pvtu
	 *     <stem>/Cycle000000/proc000000.vtu
	 *
	 * where <name> is the last path component of @a stem. The .pvd is the
	 * index; the Cycle directory is an implementation detail of the format and
	 * is not meant to be opened piece by piece. **The .pvd is INSIDE the
	 * directory, not beside it** -- which is worth stating because every other
	 * writer here produces a file at the stem, and a caller printing
	 * "<stem>.pvd" would be naming something that does not exist.
	 *
	 * HIGH-ORDER OUTPUT IS ON, AND THAT IS THE WHOLE POINT OF WRITING IT.
	 * VTK's native cells are linear, so the default path samples a P_k field
	 * at the element vertices and throws away exactly the accuracy this
	 * discretisation exists to buy -- a k = 3 solution would be drawn as if it
	 * were k = 1, and the picture would be wrong in a way that looks like a
	 * coarse mesh rather than like a bug. SetHighOrderOutput() writes VTK
	 * Lagrange cells instead, and @a levelsOfDetail is the subdivision they
	 * carry. Pass the polynomial degree.
	 *
	 * @param stem            path stem, as writeMfem().
	 * @param potential       psi.
	 * @param field           the POLOIDAL FIELD B, not the HDG flux q: this
	 *                        file is for looking at, and B is the physical
	 *                        quantity. meq::poloidalField() converts. The
	 *                        exact q is in "<stem>_grad_psi.gf".
	 * @param levelsOfDetail  subdivision per element; the polynomial degree is
	 *                        the right value. Clamped to at least 1.
	 */
	void writeVtu( std::string const &stem, mfem::Mesh &mesh,
	               mfem::GridFunction const &potential,
	               mfem::GridFunction const &field,
	               int levelsOfDetail );

	/**
	 * The domain boundary of @a mesh, as an ordered closed polyline.
	 *
	 * For a FITTED run this is Gamma itself: the plasma boundary is the mesh
	 * boundary, and psi = 0 is imposed on it. For a run on an extracted
	 * subdomain it is Gamma_h, the polygonal approximation the datum is
	 * actually imposed on -- which is a different curve from the smooth Gamma
	 * the user specified, and the difference is the whole subject of GS-2.
	 *
	 * ORDERED, because the point of it is to be drawn. The boundary elements
	 * come out of MFEM in no particular order, so a caller plotting them as
	 * given gets a star of chords rather than an outline; this walks the
	 * vertex adjacency instead. @a r and @a z are cleared first and the loop is
	 * NOT repeated at the end -- a reader closing the curve appends the first
	 * point, which is what tools/plot_equilibrium.py does.
	 *
	 * ONLY THE LOOP CONTAINING THE FIRST BOUNDARY ELEMENT IS RETURNED. meq
	 * solves a simply connected domain, so there is one; a mesh with a hole
	 * would have two and this would quietly describe the wrong one, so it
	 * reports how many boundary vertices it did not reach.
	 *
	 * @param unreached  set to the number of boundary vertices not on the
	 *                   returned loop. Non-zero means the domain is not
	 *                   simply connected and the answer is partial.
	 */
	void boundaryPolyline( mfem::Mesh &mesh,
	                       std::vector<double> &r, std::vector<double> &z,
	                       int &unreached );

	/**
	 * Bend the mesh's boundary out onto the true Gamma, for OUTPUT ONLY.
	 *
	 * On the curved path the solve happens on Omega_h, whose boundary Gamma_h
	 * is a polygon inscribed in the smooth Gamma. Every picture drawn from it
	 * therefore has a faceted edge that is not the boundary anybody asked for.
	 * This installs a curvature of @a order on the mesh and moves the nodes of
	 * each boundary face out onto Gamma, so the outer layer of elements gains a
	 * curved edge and the drawn domain is Omega rather than Omega_h.
	 *
	 * IT NEEDS NO ADVANCED VTK, which is the pleasant surprise here: writeVtu()
	 * already emits VTK Lagrange cells, and a curvilinear MFEM mesh is exactly
	 * what those represent. The two features were independent and turn out to
	 * compose.
	 *
	 * **THIS CHANGES THE GEOMETRY AND MUST BE THE LAST THING DONE.** The
	 * GridFunction coefficients are untouched, but the map from reference to
	 * physical space is not, so anything that samples, integrates or writes the
	 * mesh afterwards sees a different domain. Call it after writeMfem() and
	 * after the grid sampling, immediately before writeVtu().
	 *
	 * MOVING A BOUNDARY BY O( h ) CAN TURN AN ELEMENT INSIDE OUT, and a tangled
	 * element renders as a black spike rather than as nothing. So the whole
	 * displacement is applied, every element's Jacobian determinant is checked,
	 * and on failure the displacement is halved and the check repeated. If no
	 * fraction works the mesh is left exactly as it was found.
	 *
	 * @param order    curvature to install; the solve's polynomial degree is
	 *                 the right value, since that is what writeVtu() subdivides
	 *                 to anyway.
	 * @param project  ( r, z ) on Gamma_h -> the corresponding point on Gamma.
	 *                 A radial projection is what the shape supports and what
	 *                 the driver passes.
	 * @param applied  set to the fraction of the displacement that survived the
	 *                 tangling check: 1 normally, 0 if the mesh was left alone.
	 * @return the number of boundary nodes moved.
	 */
	int curveBoundaryOnto( mfem::Mesh &mesh, int order,
	                       std::function<void( double, double,
	                                           double &, double & )> const &project,
	                       double &applied );

	/*
	 * The ( R, Z ) grid file.
	 *
	 * Layout, which is fixed and is what a reader may rely on:
	 *
	 *     dimensions   R = nodesR, Z = nodesZ, boundary = however many
	 *     double  R( R )              major radius, metres
	 *     double  Z( Z )              height, metres
	 *     double  psi( Z, R )         poloidal flux per radian, Wb/rad
	 *     double  B_R( Z, R )         radial field, T
	 *     double  B_Z( Z, R )         vertical field, T
	 *     byte    inside( Z, R )      1 where the node lies in the domain
	 *     double  boundary_R( boundary ), boundary_Z( boundary )   optional
	 *
	 * ( Z, R ) with R fastest is C row-major, and matches MaNTA's ( t, x ).
	 *
	 * A NODE OUTSIDE THE DOMAIN gets both a NaN _FillValue and a zero in
	 * `inside`. Both, deliberately: some tools honour the fill attribute and
	 * some do not, and `inside` is the one a reader can always rely on.
	 */
	class NetCDFWriter
	{
		public:
			/// @throws std::runtime_error if meq was built without netcdf-cxx4,
			///         or if @a path cannot be created.
			NetCDFWriter( std::string const &path, GridSampler const &sampler );
			~NetCDFWriter();

			NetCDFWriter( NetCDFWriter const & ) = delete;
			NetCDFWriter &operator=( NetCDFWriter const & ) = delete;

			/// A global attribute. Provenance belongs in the file: which
			/// configuration produced it, which source, what degree, how many
			/// Newton steps it took.
			void attribute( std::string const &name, std::string const &value );
			void attribute( std::string const &name, double value );
			void attribute( std::string const &name, int value );

			/// One ( Z, R ) field. @a values is the sampler's layout, R fastest,
			/// with NaN where the node was not located.
			void field( std::string const &name, std::vector<double> const &values,
			            std::string const &longName, std::string const &units );

			/// The prescribed boundary, sampled, so a plot can draw Gamma without
			/// re-deriving it from the configuration.
			void boundary( std::vector<double> const &r,
			               std::vector<double> const &z );

			/// Flush and close. Called by the destructor; call it explicitly to
			/// see an error rather than have it thrown from a destructor.
			void close();

		private:
			struct State;
			std::unique_ptr<State> state;
	};
}

#endif // MEQ_OUTPUT_HPP
