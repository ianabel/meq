#ifndef MEQ_SAMPLER_HPP
#define MEQ_SAMPLER_HPP

/*
 * Evaluating a finite element field on a uniform ( R, Z ) grid.
 *
 * The driver's output is a rectangular grid of psi and B, because that is what
 * every downstream tool expects and because a structured grid can be
 * interpolated back in O( 1 ) per point -- which is what makes a warm start from
 * a file cheap, and is why the NetCDF file is the INTERCHANGE format rather than
 * a restart format. docs/running.rst names the three restart routes and
 * docs/output.rst says what this one costs: bilinear interpolation on a grid is
 * SECOND ORDER however fine the grid is, so its accuracy is a property of the
 * file where the .gf route's is a property of the solve.
 *
 * THE OBVIOUS ROUTE IS TOO SLOW, and this is recorded in CLAUDE.md as a trap:
 * mfem::Mesh::FindPoints is O( elements x points ), a brute-force scan over
 * element centres. A 129x129 grid against a 20k-element mesh is 3.3e8 element
 * tests for a job that should be instant.
 *
 * A SEARCH STRUCTURE IS NOT THE ANSWER HERE, AND IT IS WORTH SAYING WHY RATHER
 * THAN LEAVING THE ABSENCE TO LOOK LIKE ONE. MEQ has a point locator of its own
 * -- meq::PointLocator, which is what meq::FieldTransfer is -- and MFEM offers
 * FindPointsGSLIB wherever MFEM_USE_GSLIB is set. Neither is reached from this
 * class, because the inversion below is already linear and needs no search
 * structure at all; not because one could not be had.
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

#include <algorithm>
#include <functional>
#include <limits>
#include <vector>

#include "mfem.hpp"

namespace meq
{
	/**
	 * How far outside the reference element a located point may sit.
	 *
	 * `mfem::InverseElementTransformation`'s own default -- `ip_tol`, applied
	 * to the same `Geometry::CheckPoint()` test -- so the affine path accepts
	 * exactly the candidates `TransformBack()` accepts rather than nearly them.
	 *
	 * **AND UPSTREAM ALREADY SCALES IT WITH THE FORMAT, WHICH IS WHY THIS IS A
	 * COPIED PAIR AND NOT AN `eps` EXPRESSION.** `eltrans.hpp`'s constructor
	 * carries `ref_tol( 1e-15 ), phys_rtol( 4e-15 ), ip_tol( 1e-8 )` in double
	 * and `ref_tol( 4e-7 ), phys_rtol( 1e-6 ), ip_tol( 1e-4 )` in single -- the
	 * opposite of the `ExtensionBoundaryQuadrature` finding of M-183, and worth
	 * recording as such. A floor like `128 eps` would be `1.5e-05`, SIX TIMES
	 * TIGHTER than what MFEM's own inverse accepts, so the affine path would
	 * refuse points the Newton path takes -- which is exactly the disagreement
	 * `SamplerConvergence` reported at node ( 24, 0 ), on the domain's bottom
	 * edge.
	 *
	 * **IT IS HERE RATHER THAN IN `Sampler.cpp` BECAUSE TWO PARTIES NEED IT.**
	 * That test compares MEQ's affine inverse against MFEM's Newton one and has
	 * to ask both the SAME question, so it hands this number to
	 * `InverseElementTransformation::SetElementTol()`. A constant kept in step
	 * by hand across two files IS the defect -- M-165 and M-189 -- so it is one
	 * shared expression. M-184, M-220.
	 */
	inline double samplerInsideTolerance()
	{
#ifdef MFEM_USE_SINGLE
		return 1.0e-4;
#else
		return 1.0e-8;
#endif
	}

	/**
	 * The physical relative tolerance to hand `TransformBack()`, because its
	 * own default is BELOW ONE `float` EPSILON and therefore unreachable.
	 *
	 * **`ElementTransformation::TransformBack( v, ip, phys_rel_tol = tol_0 )`
	 * DEFAULTS TO A TIGHTER NUMBER THAN THE CLASS IT CONSTRUCTS.**
	 * `eltrans.hpp:62` carries `tol_0 = 1e-15` in double and `1e-7` in single,
	 * while `InverseElementTransformation`'s own constructor carries
	 * `phys_rtol( 4e-15 )` and `phys_rtol( 1e-6 )`. The two disagree in both
	 * formats and only in one does the disagreement bite: `1e-15` is **4.5
	 * eps** of a double and `1e-7` is **0.84 eps** of a float. `tol_0` was
	 * scaled by the ratio of the two formats' decimal digits and not by the
	 * ratio of their epsilons, and it landed under the floor.
	 *
	 * **WHAT THAT COSTS IS THE WHOLE VERDICT, NOT A DIGIT.** `NewtonSolve()`
	 * exits on `err_phys < phys_rtol*|pt|_inf` or on `dx_norm < ref_tol`, and
	 * in single an affine triangle reaches NEITHER. The physical residual of an
	 * exactly-inverted affine map is round-off, `eps*|pt|`, which is not
	 * strictly below `0.84 eps*|pt|`; the reference increment that round-off
	 * then produces is `eps*|pt|/h`, which exceeds `ref_tol = 4e-7` for every
	 * `h` below `0.31`. So the solve runs to `max_iter` and returns
	 * **`Unknown`** -- the third answer, which a caller testing `== Inside`
	 * reads as *not this element* for the element the point is actually in.
	 * Measured on `ContourTracer`'s element walk: **108, 240 and 95** falls
	 * back to `Mesh::FindPoints` over three surfaces where the double build
	 * falls back none, which is the trace going quadratic in the mesh.
	 *
	 * So this is `samplerInsideTolerance()`'s sibling with the sign of the
	 * finding reversed. There, MFEM already scaled `ip_tol` with the format and
	 * MEQ's own floor was six times TIGHTER than upstream's; here upstream's
	 * own default is the thing that does not survive the format, and copying it
	 * is what must not be done. The class default is taken instead -- 8.4 eps,
	 * comfortably reachable -- and the double arm is `tol_0` itself, so every
	 * double call is the call it always was, to the bit. M-222.
	 */
	inline double inverseTransformTolerance()
	{
#ifdef MFEM_USE_SINGLE
		return 1.0e-6;
#else
		return 1.0e-15;
#endif
	}



	/**
	 * The offset to hand `mfem::MarkLevelSetSubdomain()`, so that no vertex of
	 * `Gamma_h` can land ON `Gamma` to within the format's own ability to place
	 * it.
	 *
	 * @param gradient    the level set's magnitude of change per unit length.
	 *                    **1 for a distance or a radial gap**, which is what
	 *                    `meq::BoundaryShape::levelSet()` returns; `| grad psi |`
	 *                    for a flux level set, because the offset is in the
	 *                    level set's OWN units and those differ per caller.
	 * @param coordinate  a representative coordinate magnitude, whose ulp is
	 *                    what a vertex can be placed to.
	 *
	 * **THE MARKING IS A DISCRETE DECISION ON A CONTINUOUS PREDICATE, AND AT
	 * `offset = 0` THE SET IT DECIDES ON IS EXACTLY THE SET ROUND-OFF CANNOT
	 * DECIDE.** `MarkLevelSetSubdomain()` keeps an element when `phi <= -offset`
	 * at every vertex and at every refined point, so with `offset` zero an
	 * element whose closest approach to `Gamma` is under one ulp is admitted or
	 * refused by the last bits of `phi`. Two things then go wrong, and both are
	 * measured rather than feared.
	 *
	 * **A ZERO-LENGTH PATH.** `mfem::VertexConePath` builds the direction at a
	 * vertex by shooting a fan, averaging the nearest endpoints and DIFFERENCING
	 * that mean from the vertex. A vertex on `Gamma` has every endpoint at
	 * distance zero, so in single the endpoints round to the vertex itself, the
	 * difference is the zero vector, `len <= 0` sends all three widening passes
	 * to `continue`, and the constructor aborts with *"no ray from the vertex
	 * ... reaches Gamma"*. Measured on `FreeBoundaryCoupling`'s circle of radius
	 * 0.40 about `R = 1.10`: the Cartesian grid puts a vertex at
	 * `( 1.1, -0.4 )`, the circle's own lowest point, at `phi = -2.4e-08`. In
	 * double the same vertex is on `Gamma` too and the same arithmetic works,
	 * because a bisected crossing at `6e-14` is still representable as a
	 * displacement from a coordinate of 1.1 -- a double's ulp there is
	 * `2.2e-16` and a float's is `1.2e-07`.
	 *
	 * **AND TWO ROUTES TO ONE MARKING DISAGREE.** `DriverAcceptance` marks the
	 * same shape twice -- once through the driver and once by hand -- and reads
	 * **451 elements against 450**. The two build their background through
	 * different arithmetic, so their vertices differ in the last bits, which is
	 * nothing until an element straddles.
	 *
	 * **EXACTLY ZERO IN DOUBLE**, so `D_h` is the set it has always been and
	 * every published table stands. In single it shrinks `Omega_h` by at most
	 * the sliver that was undecidable, which is a direction the analysis
	 * already assumes: `Gamma_h` is inscribed. M-222.
	 */
	inline double markingOffset( double gradient, double coordinate )
	{
#ifdef MFEM_USE_SINGLE
		return 64.0*std::abs( gradient )*std::abs( coordinate )
		       *static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() );
#else
		( void ) gradient;
		( void ) coordinate;
		return 0.0;
#endif
	}


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
			///                    ( R_max - R_min )/( nR - 1 ).
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
			 * with `grad psi = R q` -- MEQ's flux convention, Field.hpp. Its
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
			 * form -- `div q = -F/R` is the equation being solved, and
			 * `d_r q_z - d_z q_r = -q_z/R` follows from `R q = grad psi`. Those
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

			/**
			 * The located nodes, grouped by the element that holds them.
			 *
			 * EVERY SAMPLING PASS IS A LOOP OVER ELEMENTS AND NOT OVER NODES,
			 * and this is what makes that possible. A field's value at a node
			 * needs three things from its element -- the dof list, the field's
			 * coefficients on those dofs, and the shape function -- and only
			 * the last of them depends on the node. Asking
			 * GridFunction::GetValue() per node, which is the obvious loop,
			 * therefore re-fetches the first two every time: an Array<int>, a
			 * gathered Vector and a Vector for the shape, which is THREE HEAP
			 * ALLOCATIONS and a gather for a dot product of ten numbers.
			 * Measured on the DIII-D case, that is 0.35 us per node against
			 * about 0.35 us of arithmetic per ELEMENT.
			 *
			 * AND IT IS ALSO WHAT WOULD MAKE THE PASSES PARALLEL. Grouped this
			 * way each element writes a disjoint set of node indices, so the
			 * outer loop is a forall with no reduction and no contention; a
			 * node-major loop over a scattered element map is neither.
			 *
			 * CSR over element index: groupStart has GetNE() + 1 entries and
			 * groupNode holds node indices, increasing within each group, so
			 * the order a pass visits nodes in does not depend on how the
			 * groups were built.
			 */
			void buildGroups() const;

			mutable std::vector<int> groupStart;
			mutable std::vector<int> groupNode;
			mutable bool groupsValid = false;

			mfem::Mesh &mesh;
			double minRadius, maxRadius, zMin, zMax;
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
