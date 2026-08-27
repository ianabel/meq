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
