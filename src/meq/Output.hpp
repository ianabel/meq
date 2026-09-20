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
 * INTERCHANGE format -- see docs/running.rst -- because a structured grid
 * interpolates back in O( 1 ) per point with no mesh search, so a foreign code
 * needs to produce nothing but psi on a rectangle. What that costs is second
 * order in the GRID spacing at every k, which docs/output.rst states as a limit
 * rather than leaving to be discovered.
 *
 * WHICH POTENTIAL IS IN WHICH FILE, because there are now two and they are not
 * interchangeable. psi_h is the solved potential in P_k; psi* is the
 * element-local post-processing of it in P_(k+1), which converges one order
 * faster. The driver reports psi*: the VTK and the NetCDF carry it, and the
 * NetCDF says so in a `potential` attribute so a reader never has to guess.
 * "<stem>_psi.gf" keeps psi_h, because that trio is the exact restart format and
 * is read back into a degree-k space; psi* goes beside it in
 * "<stem>_psistar.gf". See writePostProcessed().
 */

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "FluxFamily.hpp"
#include "Sampler.hpp"

namespace meq
{
	/// True when MEQ was built with netcdf-cxx4. Everything in
	/// NetCDFWriter throws without it.
	bool hasNetCDF();

	/// The mesh and the grid functions, in MFEM's own formats. Exact, and
	/// readable only by MFEM.
	///
	/// @param stem       the path stem: "<stem>.mesh", "<stem>_psi.gf" and
	///                   "<stem>_grad_psi.gf" are written.
	/// @param potential  psi_h, THE SOLVED POTENTIAL IN P_k -- not psi*. This
	///                   trio is the exact restart format, read back into a
	///                   degree-k potential space, and psi* is degree k+1. See
	///                   writePostProcessed().
	/// @param flux       q in MEQ's sign convention -- what
	///                   GradShafranovSolver::flux() returns, not the raw block.
	void writeMfem( std::string const &stem, mfem::Mesh &mesh,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux );

	/**
	 * psi*, the post-processed potential, as "<stem>_psistar.gf".
	 *
	 * A SEPARATE FUNCTION RATHER THAN A FOURTH ARTEFACT OF writeMfem(), and
	 * deliberately. psi* exists only after GradShafranovSolver::postProcess(),
	 * which is a step of its own; folding it into writeMfem() would make the
	 * mesh and the exact restart pair depend on a post-processing they have
	 * nothing to do with, and would oblige every caller -- the tests included --
	 * to have one. Keeping the two apart also keeps writeMfem()'s three files
	 * byte for byte what they were.
	 *
	 * IT IS NOT A RESTART FILE. The exact restart reads "<stem>_psi.gf" back
	 * into the degree-k potential space of a solver built from the same
	 * configuration; this is degree k+1 and does not fit it. What it is for is
	 * GLVis, and for anyone who wants the reported field at full precision
	 * rather than sampled onto the (R, Z) grid.
	 *
	 * Precision 16, as writeMfem(), for the same reason.
	 *
	 * @param postProcessed  GradShafranovSolver::postProcessedPotential(),
	 *                       valid only after postProcess().
	 */
	void writePostProcessed( std::string const &stem,
	                         mfem::GridFunction const &postProcessed );

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
	 * carry. Pass THE DEGREE OF THE POTENTIAL BEING DRAWN -- which for the
	 * driver is k+1, not k, because what it draws is psi*.
	 *
	 * @param stem            path stem, as writeMfem().
	 * @param potential       psi. The driver passes psi*, the post-processed
	 *                        potential, which is a degree richer and an order
	 *                        more accurate than psi_h; hence the note on
	 *                        @a levelsOfDetail.
	 * @param field           the POLOIDAL FIELD B, not the HDG flux q: this
	 *                        file is for looking at, and B is the physical
	 *                        quantity. meq::poloidalField() converts. The
	 *                        exact q is in "<stem>_grad_psi.gf".
	 * @param levelsOfDetail  subdivision per element; the degree of @a potential
	 *                        is the right value. Clamped to at least 1.
	 */
	void writeVtu( std::string const &stem, mfem::Mesh &mesh,
	               mfem::GridFunction const &potential,
	               mfem::GridFunction const &field,
	               int levelsOfDetail );

	/**
	 * A VTK TIME SERIES, one frame per adaptive cycle.
	 *
	 * The adaptive loop produces a mesh and a solution per cycle and until now
	 * threw all but the last away. ParaView reads a collection with several
	 * cycles as a time series and will scrub through it, mesh and all, so the
	 * refinement can be watched rather than inferred from a table of element
	 * counts.
	 *
	 * WRITTEN TO "<stem>_cycles", A SEPARATE COLLECTION FROM THE ANSWER, and
	 * deliberately: "<stem>" is the converged state and has its boundary bent
	 * onto Gamma by curveBoundaryOnto(), which mutates the mesh. Doing that
	 * mid-loop would hand the next refinement a geometry the estimator never
	 * saw. **The frames here are therefore uncurved** -- Gamma_h, faceted, as
	 * solved -- which is also the honest thing to animate, since it is the
	 * domain each cycle actually used.
	 *
	 * Each append() takes its fields by reference and writes immediately; it
	 * retains nothing, so the caller may destroy the solver on the next line.
	 */
	class VtuSeries
	{
		public:
			/// @param stem            "<stem>_cycles" is the collection.
			/// @param levelsOfDetail  as writeVtu(): the polynomial degree.
			VtuSeries( std::string const &stem, int levelsOfDetail );
			~VtuSeries();

			VtuSeries( VtuSeries const & ) = delete;
			VtuSeries &operator=( VtuSeries const & ) = delete;

			/// One frame. @a field is the poloidal field B, as writeVtu().
			/// @a time is what ParaView's slider shows; the cycle index is the
			/// natural choice and there is no physical time here.
			void append( mfem::Mesh &mesh,
			             mfem::GridFunction const &potential,
			             mfem::GridFunction const &field,
			             int cycle, double time );

			/// How many frames have been written.
			int frames() const;

		private:
			struct State;
			std::unique_ptr<State> state;
	};

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
	 * ONLY THE LOOP CONTAINING THE FIRST BOUNDARY ELEMENT IS RETURNED. MEQ
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
	 *                 the right value -- it is the geometry that is being bent,
	 *                 and Gamma_h is a facet of the degree-k solve. (writeVtu()
	 *                 now subdivides one degree finer than this, because it
	 *                 draws psi*; the two numbers used to coincide and no longer
	 *                 do.)
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
	 *     byte    extrapolated( Z, R ) 1 where it was continued, not solved on
	 *     double  boundary_R( boundary ), boundary_Z( boundary )   optional
	 *
	 * ( Z, R ) with R fastest is C row-major, and matches MaNTA's ( t, x ).
	 *
	 * A NODE OUTSIDE THE DOMAIN gets both a NaN _FillValue and a zero in
	 * `inside`. Both, deliberately: some tools honour the fill attribute and
	 * some do not, and `inside` is the one a reader can always rely on.
	 *
	 * `inside` AND `extrapolated` ANSWER DIFFERENT QUESTIONS and a reader that
	 * conflates them will be misled. `inside` is "is there data here"; a band
	 * node scores 1 and should, since it is inside the plasma and carries a real
	 * value. `extrapolated` is "was this solved for", and the same node scores 1
	 * there too -- its value was continued outward from Gamma_h and is an order
	 * less accurate than its neighbours. Drop those nodes before computing an
	 * error norm or differencing two runs at different resolutions.
	 */
	class NetCDFWriter
	{
		public:
			/// @throws std::runtime_error if MEQ was built without netcdf-cxx4,
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

			/**
			 * A ONE-DIMENSIONAL ARRAY ON A DIMENSION OF ITS OWN, WHICH IS WHAT
			 * CS-6 NEEDS AND field() CANNOT GIVE.
			 *
			 * `field()` writes on the ( Z, R ) lattice, which is the whole of
			 * what this file held: a RASTERIZATION. `COIL-SUBTRACTION-PLAN.md`
			 * §9 makes it carry a second representation beside it -- every
			 * `P_k` coefficient of the solved fields, and the conductor set
			 * that says what a REMAINDER is a remainder from -- and neither is
			 * grid shaped.
			 *
			 * The dimension is created on first use and reused afterwards, so
			 * several arrays share one; a second array of a different length
			 * on the same name is refused rather than silently truncated.
			 *
			 * @throws std::runtime_error if @a dimension already exists at a
			 *         different length, or if @a values is empty -- NetCDF's
			 *         zero-length dimension is UNLIMITED, which is not what a
			 *         caller writing no conductors means.
			 */
			void vector( std::string const &dimension, std::string const &name,
			             std::vector<double> const &values,
			             std::string const &longName,
			             std::string const &units );

			/// The same, for an integer column -- a conductor's KIND, which is
			/// an enumeration and not a measurement.
			void vector( std::string const &dimension, std::string const &name,
			             std::vector<int> const &values,
			             std::string const &longName );

			/// Flush and close. Called by the destructor; call it explicitly to
			/// see an error rather than have it thrown from a destructor.
			void close();

		private:
			struct State;
			std::unique_ptr<State> state;
	};

	/*
	 * The ( Psi, theta ) flux-surface file: INVERSION-PLAN.md stage IN-6, and
	 * the fourth of meq's output formats.
	 *
	 * WHAT IT IS FOR, AND WHY THE OTHER THREE DO NOT COVER IT. The .mesh/.gf
	 * trio is exact and readable only by MFEM; the .vtu is a picture; the .nc
	 * grid file is psi and B on a uniform ( R, Z ) lattice, which is what a
	 * plotting tool wants and is the wrong shape entirely for a 1-D transport
	 * code. A transport code reads scalar functions of a FLUX LABEL -- V'( rho ),
	 * < R^-2 >( rho ), the geometry of each surface -- and reconstructing those
	 * from a rasterised psi means re-doing the whole of the inversion item at the
	 * far end. This file is that reduction, done once, by the code that has q in
	 * hand.
	 *
	 * Layout, which is fixed and is what a reader may rely on:
	 *
	 *     dimensions  flux = surfaces, theta = angles
	 *     double  rho( flux )                the flux label, sqrt( Psi_N )
	 *     double  normalised_flux( flux )    Psi_N, 0 on the axis, 1 at the edge
	 *     double  psi( flux )                the level itself, Wb/rad
	 *     double  theta( theta )             poloidal angle about the axis, rad
	 *     double  R( flux, theta )           major radius of the surface, m
	 *     double  Z( flux, theta )           height, m
	 *     byte    extrapolated( flux, theta ) 1 where the node is band data
	 *     double  V_prime( flux )            dV/dpsi, m^3 / ( Wb/rad )
	 *     double  volume( flux )             enclosed volume, m^3
	 *     double  cross_section_area( flux ) enclosed poloidal area, m^2
	 *     double  arc_length( flux )         poloidal circumference, m
	 *     double  surface_area( flux )       area of the toroidal surface, m^2
	 *     double  inverse_R_squared( flux )          < R^-2 >, m^-2
	 *     double  grad_psi_squared_over_R_squared( flux )
	 *     double  abs_grad_psi( flux )               < | grad psi | >
	 *     double  grad_psi_squared( flux )           < | grad psi |^2 >
	 *     double  safety_factor( flux )      RoPP (142). PRESENT ONLY when the
	 *                                        family carries a g( psi )
	 *     byte    band( flux )               1 where any node of the surface is
	 *                                        band data
	 *     double  worst_residual( flux )     worst | psi_h - level | on it
	 *     double  transversality( flux )     min | u x t | over the fit
	 *
	 * ( flux, theta ) with theta fastest is C row-major, the same convention the
	 * ( Z, R ) grid file uses.
	 *
	 * THE LABEL IS rho = sqrt( Psi_N ) AND THE FILE SAYS SO, in the global
	 * attribute `flux_label`. Psi_N is carried beside it because a consumer
	 * whose own grid is in normalised flux should not have to square anything,
	 * but rho is what the geometry is smooth in and what meq interpolates
	 * against. FluxFamily.hpp has the measurement.
	 *
	 * `extrapolated` IS A MASK AND `extrapolated_nodes` IS A COUNT, and the
	 * distinction is the one the ( R, Z ) writer had to be repaired for:
	 * CLAUDE.md records that carrying only the count meant nothing downstream
	 * could tell WHICH of 1667 nodes had been continued into the band between
	 * Gamma_h and Gamma. `band( flux )` is the per-surface summary of the same
	 * thing, so a consumer can drop a whole surface rather than a node.
	 *
	 * THERE IS NO FILL VALUE AND NO `inside`, WHICH IS THE DIFFERENCE FROM THE
	 * GRID FILE. A ( R, Z ) lattice has nodes outside the domain and has to say
	 * so; a family has none. Every surface in the file was traced, fitted and
	 * integrated, and a level that could not be is not in the file at all --
	 * meq::extractFluxSurfaces() abandons the family rather than leaving a hole
	 * an interpolation would bridge. So a reader does not have to check for
	 * missing data; it has to check `extrapolated` and `worst_residual`, which
	 * are about how good the data is rather than whether it is there.
	 */
	class FluxGridWriter
	{
		public:
			/// @param path   the file to create.
			/// @param family what to write. Everything is written here, at
			///               construction; attributes are added afterwards, as
			///               NetCDFWriter does.
			/// @throws std::runtime_error if MEQ was built without netcdf-cxx4,
			///         or if @a path cannot be created.
			/// @throws std::invalid_argument if the family is empty or its
			///         surfaces disagree about how many nodes they carry, since
			///         a rectangular array cannot be written from a ragged one.
			FluxGridWriter( std::string const &path,
			                FluxSurfaceFamily const &family );
			~FluxGridWriter();

			FluxGridWriter( FluxGridWriter const & ) = delete;
			FluxGridWriter &operator=( FluxGridWriter const & ) = delete;

			/// A global attribute, as NetCDFWriter's. The family's own
			/// provenance -- the axis, the two flux values, the cut, the label
			/// -- is written by the constructor and needs no help.
			void attribute( std::string const &name, std::string const &value );
			void attribute( std::string const &name, double value );
			void attribute( std::string const &name, int value );

			/// Flush and close. Called by the destructor; call it explicitly to
			/// see an error rather than have it thrown from a destructor.
			void close();

		private:
			struct State;
			std::unique_ptr<State> state;
	};
}

#endif // MEQ_OUTPUT_HPP
