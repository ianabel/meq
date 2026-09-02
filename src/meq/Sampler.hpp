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
 * tests for a job that should be instant.
 *
 * GSLIB IS NOT THE ANSWER HERE, AND THIS COMMENT USED TO SAY IT WAS
 * UNAVAILABLE, WHICH IS FALSE -- the install has MFEM_USE_GSLIB = YES and
 * meq::FieldTransfer already uses mfem::FindPointsGSLIB. It is not used in this
 * class because the inversion below is already linear and needs no search
 * structure at all, not because it could not be.
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
			///
			/// THIS READS THE FOOT ON Gamma_h FOR A BAND NODE, so it is `O( h )`
			/// there. For `B`, or anything else written to the interchange file,
			/// use sampleComponentWithGradient() instead; this one is for a
			/// caller that wants the value in the element and nothing else.
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
			 * @param gapToBoundary  optional distance from ( R, Z ) to the true
			 *                Gamma. Supplying it makes blendWeight() meaningful,
			 *                which is what lets a caller pull the extrapolation
			 *                back onto a boundary value it knows -- see there.
			 * @return how many nodes were newly filled.
			 */
			int extendOutward( double reach,
			                   std::function<bool( double, double )> const &accept
			                       = std::function<bool( double, double )>(),
			                   std::function<double( double, double )> const
			                       &gapToBoundary
			                       = std::function<double( double, double )>() );

			/**
			 * Continue the potential across the band using the FLUX, which the
			 * mixed method computes at the SAME order as the potential itself.
			 *
			 * THIS IS WHY THE BAND IS WORTH FILLING AT ALL, and it replaces
			 * extrapolating psi_h's own polynomial outside its element -- which
			 * is bounded by nothing and was measured crossing a value known
			 * exactly (see blendWeight). Here nothing is evaluated outside an
			 * element: extendOutward() records the FOOT of each band node on
			 * Gamma_h, which lies on its element's own boundary, and both fields
			 * are read there. The step outward is then a Taylor extension
			 *
			 *     psi( p ) = psi( x0 ) + grad psi( x0 ) . ( p - x0 )
			 *
			 * with `grad psi = r q` -- meq's flux convention, Field.hpp. Its
			 * error is O( |p - x0|^2 ) against an extrapolation with no bound,
			 * and `q` carries no extra error to spend: getting the derivative at
			 * the potential's own order rather than one down is exactly what
			 * this discretisation is for.
			 *
			 * @param flux  q as GradShafranovSolver::flux() returns it, NOT the
			 *              raw block and NOT the poloidal field.
			 */
			void samplePotentialWithFlux( mfem::GridFunction const &potential,
			                              mfem::GridFunction const &flux,
			                              std::vector<double> &values,
			                              double fill ) const;

			/**
			 * Continue one component of a VECTOR field across the band using
			 * that field's OWN gradient. This is what `B` needs, and it is what
			 * sampleComponent() does not do.
			 *
			 * The step and the guarantee are samplePotentialWithFlux()'s: the
			 * foot recorded by extendOutward() lies on its element's own
			 * closure, so `GetVectorGradient` is evaluated inside the element
			 * and nothing is read outside one.
			 *
			 *     u_c( p ) = u_c( x0 ) + grad u_c( x0 ) . ( p - x0 )
			 *
			 * IT DOES NOT REACH THE FLUX'S OWN ORDER, AND THAT IS STRUCTURAL
			 * RATHER THAN A SHORTCUT -- the two band continuations are not
			 * equally good and it would be wrong to imply they are.
			 * samplePotentialWithFlux() is as accurate as it is because `q` is a
			 * SOLVED variable carrying the potential's own order, which is the
			 * whole point of the mixed method. There is no solved variable for
			 * `grad q`: differentiating an L2 field of degree `k` leaves degree
			 * `k-1`, converging one order down, so this step is `O( h^2 )`
			 * whatever `k` is. That is still a full order better than the
			 * `O( h )` of reading the foot, which is what it replaces, and
			 * `theBandVectorContinuesAtItsGradientsOrder` measures both.
			 *
			 * IF SOMEONE WANTS BETTER, the route is known and is not this one.
			 * The divergence and the curl of `q` are both available in closed
			 * form -- `div q = -F/r` is the equation being solved, and
			 * `d_r q_z - d_z q_r = -q_z/r` follows from `r q = grad psi`. Those
			 * pin two of the four entries of `grad q` exactly and leave the
			 * symmetric traceless part still to be differentiated, so they buy
			 * STRUCTURE rather than an order, at the cost of plumbing the source
			 * into this class. Not done, and not obviously worth it.
			 */
			void sampleComponentWithGradient( mfem::GridFunction const &field,
			                                  int component,
			                                  std::vector<double> &values,
			                                  double fill ) const;

			/// Whether this node was filled by extendOutward() rather than found
			/// inside an element -- the per-node form of extendedCount().
			///
			/// A READER OF THE OUTPUT FILE IS ENTITLED TO THIS, which is why it
			/// is public: the `inside` mask says 1 in the band, because the node
			/// really was located and really does carry data, so nothing else in
			/// the file distinguishes a solved node from a continued one.
			bool wasExtended( int i, int j ) const;

			/**
			 * How far through the Gamma_h-to-Gamma band a node sits: 0 on
			 * Gamma_h, 1 on Gamma, and 0 for every node that was found inside an
			 * element rather than extrapolated.
			 *
			 * THIS EXISTS BECAUSE THE EXTRAPOLATION IS NOT TRUSTWORTHY ON ITS
			 * OWN, and that is measured rather than suspected. On
			 * examples/miller-curved.toml, where the datum makes psi exactly
			 * zero on Gamma and strictly negative inside, the raw extrapolation
			 * puts **17 nodes at positive psi**, worst +1.06e-02 against a peak
			 * of 2.5e-01 -- it crosses a value that is known exactly. A
			 * polynomial fitted inside an element is not constrained outside it,
			 * and one element-length out it is already wrong by 4% of the
			 * solution.
			 *
			 * A caller that knows the boundary value g can therefore do far
			 * better than the extrapolation alone:
			 *
			 *     v = ( 1 - t )*extrapolated + t*g
			 *
			 * which is exact at both ends -- the solve's own value where the
			 * band meets Gamma_h, the known datum on Gamma -- and monotone
			 * between them, so it cannot overshoot. Its error is O( h^2 ) in the
			 * band's width against an extrapolation bounded by nothing.
			 *
			 * Zero unless extendOutward() was given a @a gapToBoundary.
			 */
			double blendWeight( int i, int j ) const;

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
			/// How far through the band, for blendWeight().
			std::vector<double> blend;
			/// For a band node, p - x0: the step from its foot on Gamma_h.
			/// Zero for every node found inside an element.
			std::vector<double> offsetR, offsetZ;
	};
}

#endif // MEQ_SAMPLER_HPP
