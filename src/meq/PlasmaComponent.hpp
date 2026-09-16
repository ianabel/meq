#ifndef MEQ_PLASMACOMPONENT_HPP
#define MEQ_PLASMACOMPONENT_HPP

#include <vector>

/*
 * THE PLASMA IS A CONNECTED SET AND `{ Psi > 0 }` IS NOT: stage XP-1 of
 * FREE-BOUNDARY-PLAN.md section 10.6.
 *
 * meq::NormalisedSource::insidePlasma() is a POINTWISE test on the value,
 * ( psi - psi_bnd )*span > 0, with no connectivity in it at all. That is
 * correct for a plasma whose level set happens to have one component and wrong
 * whenever it does not, and section 10.3 records the direction it is wrong in:
 * across an X-point the level psi = psi_X cuts a neighbourhood into four
 * sectors and TWO OPPOSITE ONES carry Psi > 0 -- the plasma, and the private
 * flux region under the divertor. A run with `ConfineToPlasma` on then
 * converges and describes a machine with a second current channel nobody asked
 * for.
 *
 * **AND IT IS LIVE RATHER THAN LATENT, WHICH SECTION 10.3 SAID IT WAS NOT.**
 * That paragraph reads "it cannot fire today: no shipped example sets
 * ConfineToPlasma and MEQ has no diverted case". The first half is still true
 * and the second is beside the point. Measured on the LIMITER fixture of
 * tests/convergence/FreeBoundaryCoupling.cpp -- a half-disc, an exterior
 * coupling, two borders, and NO X-POINT ANYWHERE -- the converged
 * { psi > psi_bnd } is 705 elements in more than one piece, and a ray across
 * the midplane crosses two lobes of it. A pointwise support test is switching
 * the source on in several places at once on a case that ships and converges.
 *
 * The reason is more ordinary than a saddle, and it is why this is not a
 * diverted-plasma problem at all: psi_bnd is a LEVEL, and a level of any field
 * that is not monotone in radius cuts the domain into as many pieces as it
 * likes. Coil pockets, the vacuum region between the plasma and Gamma, and the
 * iterate's own transients all produce lobes.
 *
 *
 * WHAT THIS FILE IS
 * -----------------
 *
 * The flood fill, and only the flood fill: a face-neighbour connected-component
 * labelling of a graph, seeded at one node. Adjacency arrives as a CSR pair of
 * plain integer vectors, so this is a graph algorithm on plain data and NOTHING
 * HERE KNOWS WHAT A MESH IS.
 *
 * **MFEM-FREE, DELIBERATELY**, like Profiles, Source, Coils, ExteriorDtN,
 * Zernike and SurfaceFit. CI cannot obtain the MFEM branch MEQ builds against,
 * so a component that needs the library gets no test there at all; a component
 * that takes a CSR graph gets a unit test that runs on every push, and the
 * configurations that matter -- two lobes touching at a vertex, two lobes
 * bridged by a band, a band with only one component to belong to -- are far
 * easier to construct by hand as graphs than as meshes.
 *
 * The MFEM half is meq::GradShafranovSolver::refreshPlasmaComponent(): build
 * the CSR from mfem::Mesh::ElementToElementTable(), mark the elements whose
 * potential dofs reach Psi > 0, and seed at the magnetic axis.
 *
 *
 * WHY FACES AND NOT VERTICES, AND THE PREDICTION THAT DID NOT SURVIVE IT
 * ---------------------------------------------------------------------
 *
 * `../freegs4e` does this fill on its uniform ( R, z ) grid
 * ( critical.core_mask ) and has to EXPLICITLY BLOCK a neighbourhood of each
 * X-point first, because a grid fill leaks diagonally through the saddle: the
 * two lobes meet at a point, and on a grid a point is an eight-neighbour hop.
 *
 * Section 10.3 predicts MEQ needs no such blocking, because two lobes meeting
 * at a shared VERTEX are not face neighbours, and that "the expected leak is
 * ONE ELEMENT WIDE" -- the saddle's own element, which is cut by both branches.
 * MEASURED on Soloviev::iterExample2(), whose saddle is known in closed form,
 * BOTH HALVES ARE WRONG IN THE SAME DIRECTION.
 *
 *   * A vertex-touching pair really is separated -- the graph carries the
 *     topology a grid destroys, and two_lobes_meeting_at_a_vertex_are_not
 *     _connected pins that. So the premise holds.
 *   * But the two lobes are not joined at the saddle. They are joined through
 *     the BAND of elements STRADDLING the separatrix, every one of which
 *     carries Psi > 0 at some vertex and is therefore a candidate, and which
 *     near a saddle is several elements wide. The plain fill takes the WHOLE
 *     private flux region: 2275 elements of 16688 and 10.4% of int |F|.
 *   * And blocking the saddle's own element -- section 10.3's cure, which needs
 *     an X-point finder that stage XP-0 has not built -- works on a fine mesh
 *     and NOT on a coarse one: 161 of 2304 elements still reached below the
 *     saddle at the coarsest of three resolutions.
 *
 * So the answer is not a better blocking rule but a different TRAVERSAL rule.
 * The two-argument fill() below walks only elements that are unambiguously
 * inside, where the lobes are genuinely not adjacent, and shares the straddling
 * band out between the interior components by a watershed. It locates nothing,
 * it has no parameter, and it reaches the same mask the blocked fill reaches
 * where the blocked fill works -- int |F| agreeing to every printed digit.
 *
 * blockNode is kept regardless: it is one element rather than a neighbourhood,
 * which IS what the face-neighbour graph buys over a grid, and a caller who has
 * located a saddle may still want it.
 */

namespace meq
{
	/**
	 * A connected-component labelling of a graph, with one component named.
	 *
	 * The graph is the mesh's element adjacency; a node is "a candidate" when
	 * the caller says the element carries plasma; the named component is the one
	 * containing the seed. holds() is then the connectivity test the pointwise
	 * one is missing.
	 *
	 * **componentCount() IS THE DIAGNOSTIC AND NOT A BY-PRODUCT.** The whole
	 * labelling costs the same as the single fill -- one pass, O( nodes +
	 * edges ) -- and the number of components is the quantity that says whether
	 * the connectivity test is doing anything at all on a given configuration.
	 * A fill that reports one component has changed nothing, and a measurement
	 * that cannot distinguish "the fix works" from "there was nothing to fix"
	 * is not a measurement.
	 */
	class PlasmaComponent
	{
		public:
			/// An empty graph: holds() is true everywhere, which is the
			/// "no connectivity test" state and is what a default-constructed
			/// solver carries.
			PlasmaComponent() = default;

			/**
			 * The element adjacency, as CSR: the neighbours of node @a e are
			 * @a neighbours[ @a offsets[e] .. @a offsets[e+1] ).
			 *
			 * @throws std::invalid_argument if @a offsets is empty, is not
			 *         non-decreasing, does not end at @a neighbours.size(), or
			 *         if any neighbour is out of range. A malformed graph would
			 *         otherwise produce a plausible mask, and a plausible mask
			 *         is exactly what nothing downstream can check.
			 */
			void setAdjacency( std::vector< int > offsets,
			                   std::vector< int > neighbours );

			/// Nodes in the graph. Zero before setAdjacency().
			int nodeCount() const;

			/**
			 * Label every component of { e : carriesPlasma[ e ] } and name the
			 * one containing @a seed.
			 *
			 * @param carriesPlasma one entry per node; nonzero means the element
			 *        has plasma in it somewhere.
			 * @param seed the node to start from -- the element holding the
			 *        magnetic axis. If it is not itself a candidate the fill is
			 *        EMPTY and holds() is false everywhere, which is a loud
			 *        failure rather than a quiet one: a source confined to
			 *        nothing gives a vacuum field, and a vacuum field where a
			 *        plasma was asked for is visible in one glance.
			 * @param blockNode a node the fill may not pass through, or -1.
			 *        This is the X-point's own element; see the header comment.
			 *
			 * @throws std::invalid_argument if @a carriesPlasma is not one entry
			 *         per node, or @a seed is out of range.
			 */
			void fill( std::vector< char > const &carriesPlasma, int seed,
			           int blockNode = -1 );

			/**
			 * THE TWO-RULE FILL, and it is what makes an X-point separate
			 * WITHOUT anybody having to find it.
			 *
			 * MEASURED on Soloviev::iterExample2(), whose saddle is known in
			 * closed form, over three candidate rules:
			 *
			 *     any vertex Psi > 0      4286 candidates, ONE component
			 *     element centre Psi > 0  4058 candidates, TWO
			 *     every vertex Psi > 0    3821 candidates, TWO
			 *
			 * So the leak is a property of the INCLUSIVE candidate rule and not
			 * of the graph. Two lobes meeting at a point are joined by the band
			 * of elements STRADDLING the separatrix, which near the saddle is
			 * several elements wide -- section 10.3 predicted a leak "one
			 * element wide" and the plain fill in fact takes the whole private
			 * flux region, 2275 elements and 10.4% of int |F|.
			 *
			 * The exclusive rule separates them and costs an O( h ) band of the
			 * plasma edge, where the source is switched off although Psi > 0
			 * there. That is not affordable: FB-4's result is that psi* keeps
			 * k+2 exactly when k <= j, and an element-aligned support would put
			 * an O( h^(1+j) ) perturbation under it.
			 *
			 * So the two rules are used for two different things. The fill
			 * TRAVERSES the strictly interior elements, where the lobes really
			 * are separate. The straddling band is then shared out between the
			 * components it touches by a WATERSHED -- a breadth-first wave from
			 * every interior component at once, each straddling element going to
			 * whichever wave reaches it first -- and the mask is the seed's
			 * component together with the band assigned to it.
			 *
			 * **A WATERSHED AND NOT A FIXED NUMBER OF RINGS, WHICH WAS TRIED AND
			 * MEASURED FIRST.** Rings of a chosen depth work on the diverted
			 * fixture -- two is the only depth that both keeps every plasma
			 * element and reaches none below the saddle, at three resolutions --
			 * and they fail on an ORDINARY confined rectangle, where the band
			 * along psi = 0 is three elements thick at the corners: depth 2
			 * dropped 4 of 512 elements, moved psi_ax by 1.5e-04 and cost seven
			 * Newton steps. A depth that has to be right on every geometry is a
			 * tuning parameter, and this tree does not ship one where a
			 * parameter-free rule is available.
			 *
			 * The watershed has neither failure. Where there is ONE interior
			 * component the whole band goes to it and nothing is dropped; where
			 * there are two, the band splits at the pinch, which near a saddle is
			 * the saddle itself. What it cannot keep is a straddling element
			 * connected to no interior element at all -- a plasma thinner than
			 * one element everywhere, which is a mesh too coarse to resolve the
			 * plasma rather than a case to accommodate.
			 *
			 * @param interior elements whose plasma content is unambiguous, one
			 *        entry per node. The fill walks these and no others.
			 * @param carriesPlasma elements with any plasma in them at all. Must
			 *        be a superset of @a interior; the band is taken from here.
			 * @param seed as fill().
			 * @param blockNode as fill(). Not needed for a saddle -- that is the
			 *        point of this overload -- and kept because a caller who has
			 *        located one may still want it.
			 *
			 * @throws std::invalid_argument as fill(), and if @a carriesPlasma
			 *         is not a superset of @a interior: two rules that disagree
			 *         about which is stricter would give a mask with holes in
			 *         it, which is a plausible answer and so is refused here.
			 */
			void fill( std::vector< char > const &interior,
			           std::vector< char > const &carriesPlasma, int seed,
			           int blockNode = -1 );

			/**
			 * ELEMENTS THAT CAN NEVER BE PLASMA, WHATEVER THE FLUX SAYS THERE.
			 *
			 * A geometric statement about the DEVICE rather than about the
			 * iterate -- the far side of a vessel wall, a port, a region the
			 * mesh carries for the coils' sake -- so it is set once and read on
			 * every sweep, which is the whole economy: deciding inside/outside
			 * from a polygon costs a point-in-polygon test per element per
			 * evaluation, and an attribute costs a lookup.
			 *
			 * **IT IS NOT A SUBSTITUTE FOR THE FILL AND IT IS NOT A LEVEL SET.**
			 * `psi_bnd` already confines the plasma and the fill already
			 * separates the lobes of `{ Psi > 0 }`; what neither can do is know
			 * that a lobe is on the far side of a wall. Where the fill's
			 * judgement is right this changes nothing and only saves it work; it
			 * earns its place where several O-points sit across a saddle and
			 * connectivity alone cannot say which is the plasma.
			 *
			 * @param excluded  one entry per node, non-zero to exclude. Empty
			 *                  clears it. Sized against the adjacency when one
			 *                  is set.
			 * @throws std::invalid_argument on a size that is neither empty nor
			 *         nodeCount().
			 */
			void setExcluded( std::vector< char > excluded );

			/// Whether setExcluded() named this node. False when nothing was
			/// excluded, and for a node out of range.
			bool isExcluded( int element ) const;

			/// How many nodes setExcluded() named. Zero when it was not called.
			int excludedNodes() const;

			/// True where the plasma is. Always true before fill() EXCEPT where
			/// setExcluded() says otherwise, so a caller may test
			/// unconditionally and an unconfigured solver is unchanged.
			bool holds( int element ) const;

			/// Whether fill() has been called and holds() is therefore a real
			/// test rather than the constant true.
			bool active() const;

			/// Nodes holds() is true at: the seed's component, plus the
			/// straddling elements the watershed gave it.
			int componentNodes() const;

			/// Straddling nodes alone. Zero for the one-rule fill.
			int ringNodes() const;

			/// Nodes with carriesPlasma set, over every component. The
			/// difference from componentNodes() is what the connectivity test
			/// removed.
			int candidateNodes() const;

			/// Components among the candidates. One means the pointwise test was
			/// already right on this configuration.
			int componentCount() const;

			/// Which component a node is in, or -1 for a non-candidate. Numbered
			/// in ascending node order of their first member, so it is stable.
			int label( int element ) const;

			/// The seed's component number, or -1 if the seed was not a
			/// candidate.
			int seedLabel() const;

			/// Back to the no-connectivity state.
			void clear();

		private:
			std::vector< int > rowOffsets;
			std::vector< int > neighbourList;
			std::vector< int > labels;

			/// The two-rule fill's share of the straddling band: elements the
			/// watershed gave to the seed's component. Empty for the one-rule
			/// fill.
			std::vector< char > ring;
			std::vector< char > excludedNode;
			int excludedSize = 0;
			int labelCount = 0;
			int seedLabelValue = -1;
			int componentSize = 0;
			int ringSize = 0;
			int candidateSize = 0;
			bool filled = false;
	};
}

#endif // MEQ_PLASMACOMPONENT_HPP
