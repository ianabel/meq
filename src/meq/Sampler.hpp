#ifndef MEQ_SAMPLER_HPP
#define MEQ_SAMPLER_HPP

/*
 * Evaluating a finite element field on a uniform ( R, Z ) grid.
 *
 * The driver's output is a rectangular grid of psi and B, because that is what
 * every downstream tool expects and because a structured grid can be
 * interpolated back in O( 1 ) per point -- which is what makes a warm start from
 * a file cheap, and is why DRIVER-PLAN section 4 keeps the NetCDF file as the
 * interchange format.
 *
 * THE OBVIOUS ROUTE IS TOO SLOW, and this is recorded in CLAUDE.md as a trap:
 * mfem::Mesh::FindPoints is O( elements x points ), a brute-force scan over
 * element centres. A 129x129 grid against a 20k-element mesh is 3.3e8 element
 * tests for a job that should be instant, and this MFEM is built with
 * MFEM_USE_GSLIB = NO, so FindPointsGSLIB is not available as an escape.
 *
 * SO THE LOOP IS INVERTED. Rather than asking, for each point, which element
 * holds it, this asks for each element which points it might hold: the element's
 * bounding box converts to a grid index range by arithmetic, because the grid is
 * uniform, and only those few candidates are tested. The cost is
 * O( elements x points per element ) -- linear in both -- and the test asserts
 * that rather than the comment claiming it.
 *
 * A POINT ON AN INTER-ELEMENT FACE is claimed by whichever element reaches it
 * first. psi_h and q_h are discontinuous, so the two answers differ by the jump,
 * which is O( h^(k+1) ) and converges away. That is a real choice and not a
 * principled one; it is recorded here rather than hidden.
 *
 * BEING OUTSIDE IS NOT AN ERROR. A grid node in the bounding box but outside the
 * mesh is simply not found, and that is exactly the mask the output wants: on
 * the extension path the mesh IS Omega^h, so "found in an element" and "inside
 * Gamma" are the same question and the locator answers both at once.
 */

#include <functional>
#include <vector>

#include "mfem.hpp"

namespace meq
{
	/// A uniform ( R, Z ) grid, and where each of its nodes sits in a mesh.
	///
	/// Locating is done once at construction; sampling a field afterwards is a
	/// loop over located nodes with no searching.
	class GridSampler
	{
		public:
			/// @param mesh  the mesh to locate in. Borrowed, and must outlive
			///              this; the located element indices refer to it.
			/// @param rMinIn, rMaxIn, zMinIn, zMaxIn  the grid's extent, metres.
			/// @param nRIn, nZIn  node counts, >= 2 in each direction. These are
			///                    NODES and not cells, so the spacing is
			///                    ( rMax - rMin )/( nR - 1 ).
			GridSampler( mfem::Mesh &mesh,
			             double rMinIn, double rMaxIn, int nRIn,
			             double zMinIn, double zMaxIn, int nZIn );

			int nodesR() const { return nR; }
			int nodesZ() const { return nZ; }

			/// Node coordinates. Index j*nR + i, R fastest -- the ordering the
			/// NetCDF file uses, which is C row-major with R as the last
			/// dimension.
			double rAt( int i ) const;
			double zAt( int j ) const;

			/// True where the node was located in an element.
			bool located( int i, int j ) const;

			/// How many nodes were located, out of nodesR() * nodesZ().
			int locatedCount() const { return found; }

			/// Sample a scalar field. @a values is resized to nR*nZ, R fastest,
			/// and set to @a fill wherever the node was not located.
			void sample( mfem::GridFunction const &field,
			             std::vector<double> &values,
			             double fill ) const;

			/// Sample one component of a vector field, same layout.
			void sampleComponent( mfem::GridFunction const &field, int component,
			                      std::vector<double> &values,
			                      double fill ) const;

			/// Sample a Coefficient, same layout. This is how B is written
			/// without building a GridFunction for it.
			void sampleCoefficient( mfem::Coefficient &coefficient,
			                        std::vector<double> &values,
			                        double fill ) const;

			/**
			 * Fill nodes just OUTSIDE the mesh by extrapolating the element
			 * that owns the nearest boundary face.
			 *
			 * THIS EXISTS FOR THE CURVED PATH AND THE SLIVER IT LEAVES.
			 * `Omega_h` is the union of background elements lying inside
			 * `Gamma`, so `Gamma_h` is INSCRIBED in `Gamma` and there is a band
			 * `O(h)` wide that is inside the plasma and outside the mesh.
			 * Rasterised without this, that band is NaN, and a picture of a
			 * smooth boundary comes out with a ragged polygonal edge -- which
			 * looks like a meshing bug and is not one.
			 *
			 * IT IS EXTRAPOLATION AND THE FILLED VALUES ARE NOT THE SOLUTION.
			 * A polynomial evaluated outside its own element is only as good as
			 * its Taylor remainder, so @a reach is a hard limit rather than a
			 * tuning knob: at `O(h)` the error is the discretisation error, and
			 * far outside it is nothing at all. That `O(h)` band is exactly the
			 * regime Cockburn & Solano's transfer analysis covers, which is
			 * what makes this the same operation the extension technique
			 * already performs -- not a new approximation invented for a plot.
			 *
			 * @param reach   how far to extrapolate, as a multiple of the
			 *                boundary face's own length. Values above about 1
			 *                are outside anything anybody has analysed.
			 * @param accept  optional predicate on ( R, Z ): a node is filled
			 *                only where this returns true. The curved path
			 *                passes "inside Gamma", without which this would
			 *                also paint a band OUTSIDE the plasma, where the
			 *                solve makes no claim at all.
			 * @return how many nodes were newly filled.
			 */
			int extendOutward( double reach,
			                   std::function<bool( double, double )> const &accept
			                       = std::function<bool( double, double )>() );

			/// How many of locatedCount() were extrapolated by extendOutward()
			/// rather than found inside an element. Written to the output file,
			/// because a reader is entitled to know which nodes are the
			/// solution and which are a continuation of it.
			int extendedCount() const { return extended; }

		private:
			int index( int i, int j ) const { return j*nR + i; }

			mfem::Mesh &mesh;
			double rMin, rMax, zMin, zMax;
			int nR, nZ;
			int found;
			int extended = 0;

			/// Per node: the element holding it, or -1, and where in it.
			std::vector<int> element;
			std::vector<mfem::IntegrationPoint> point;
	};
}

#endif // MEQ_SAMPLER_HPP
