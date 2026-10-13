#ifndef MEQ_MATERIALS_HPP
#define MEQ_MATERIALS_HPP

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "Coils.hpp"

/*
 * WHAT EACH ELEMENT OF THE MESH IS MADE OF: `ELEMENT-REGIONS-PLAN.md` §2.1.
 *
 * Every element stores the fraction of it occupied by each MATERIAL -- vacuum,
 * plasma, each coil, and, reserved for later, passive and ferromagnetic
 * structures -- and a REGION is derived from those fractions for plotting. The
 * fractions are the record and the region is a convenience; §2.1's table says
 * which fractions name which region.
 *
 * **THIS HALF IS MFEM-FREE, DELIBERATELY.** The table, the region rule and the
 * bookkeeping that turns per-element weights into fractions are plain doubles
 * in and out, so they are unit-tested where CI can run them. What needs a mesh
 * -- asking the solver's own support test at the solver's own quadrature
 * points -- is meq::GradShafranovSolver::elementMaterials(), which feeds
 * ElementMaterials::addElement() below and decides nothing this file decides.
 */

namespace meq
{
	struct CoilConfig;

	/// The kind of a material-table row. The numbers are the `.nc`'s
	/// `material_kind` and do not move.
	enum class MaterialKind
	{
		Vacuum = 0,
		Plasma = 1,
		/// One row per `[[coils]]` block, whatever the conductor model.
		Coil = 2,
		/// RESERVED: no input produces one yet. Readable before it is produced,
		/// which is why it is a kind now rather than when the input arrives.
		Passive = 3,
		/// RESERVED, as Passive is.
		Ferromagnetic = 4
	};

	/// Whether a kind is a conductor -- a coil, a passive structure or iron --
	/// for §2.1b's rules: a conductor wins a point over vacuum, and a point
	/// that is plasma and conductor at once is a fault.
	bool isConductor( MaterialKind kind );

	/// "vacuum", "plasma", "coil", "passive", "ferromagnetic".
	char const *materialKindName( MaterialKind kind );

	/**
	 * THE REGION DERIVED FROM AN ELEMENT'S FRACTIONS, for plotting. The
	 * numbers are the `.nc`'s `element_region` flag values and do not move.
	 *
	 * The first five are an element holding one material whole. The other
	 * three are §2.1b's: a plasma boundary is unavoidable, a conductor boundary
	 * happens only when the mesh is not aligned to the copper, and Mixed is
	 * either a fault or a cell `[mesh.generate] AlignConductors` would split.
	 */
	enum class ElementRegion
	{
		Vacuum = 0,
		Plasma = 1,
		Coil = 2,
		Passive = 3,
		Ferromagnetic = 4,
		/// Plasma and vacuum share the element, and nothing else is in it.
		PlasmaBoundary = 5,
		/// Conductors and vacuum share the element, and nothing else is in it.
		ConductorBoundary = 6,
		/// Anything else.
		Mixed = 7
	};

	/// How many ElementRegion values there are.
	inline constexpr int elementRegionCount = 8;

	/// "vacuum", ..., "plasma_boundary", "conductor_boundary", "mixed" --
	/// the `.nc`'s `flag_meanings` spelling.
	char const *elementRegionName( ElementRegion region );

	/**
	 * ONE ROW OF THE MATERIAL TABLE, with every column of §2.1. A column that
	 * does not apply to a kind holds its fill value, so a new kind adds rows
	 * and no new columns.
	 */
	struct Material
	{
		MaterialKind kind = MaterialKind::Vacuum;
		std::string name;
		/// The circuit this row belongs to, or -1 for none. Passive
		/// structures have none; every coil has one, its `[[circuits]]`
		/// block or `default`, when the table is buildMaterialTable()'s.
		int circuit = -1;
		/// Amperes, the resolved current of a coil. NaN elsewhere.
		double current = std::numeric_limits< double >::quiet_NaN();
		/// Ohm metres, for a passive structure. NaN elsewhere.
		double resistivity = std::numeric_limits< double >::quiet_NaN();
		/// 1 except iron: a linear iron's constant, and NaN -- the fill -- for
		/// an iron whose permeability is a curve.
		double relativePermeability = 1.0;
		/// An index into the permeability-curve table for non-linear iron, or
		/// -1 for none.
		int permeabilityCurve = -1;
	};

	/**
	 * THE NAME OF THE CIRCUIT EVERY UNWIRED COIL BELONGS TO, in the output.
	 *
	 * A GROUPING, NOT A SUPPLY: its members carry independent currents, so its
	 * `circuit_current` is the fill value and `circuit_supplied` is 0. The
	 * configuration refuses `[[circuits]] Name = "default"`, since a supplied
	 * circuit of that name would merge with the unwired coils here. For the
	 * SOLVE an unwired coil is still its own one-member circuit
	 * (meq::CircuitParameters); this group exists only in the file.
	 */
	inline constexpr char const *defaultCircuitName = "default";

	/// One row of the circuit table: `circuit_name`, `circuit_current`,
	/// `circuit_supplied`.
	struct Circuit
	{
		std::string name;
		/// Amperes: the one current a supplied circuit drives, each member
		/// carrying its multiplier times it. NaN for the `default` grouping.
		double current = std::numeric_limits< double >::quiet_NaN();
		/// A `[[circuits]]` block, as against the `default` grouping.
		bool supplied = false;
	};

	/**
	 * A NON-LINEAR IRON'S B( H ), stored in the file so that a `.nc` describes
	 * its own equilibrium without a library of steels beside it. RESERVED, as
	 * the ferromagnetic kind is: nothing produces one yet, and a reader must be
	 * able to hold one.
	 *
	 * **IT CARRIES ITS MATERIAL'S NAME**, which repeats the pointing row's
	 * `material_name` on purpose: a reader checks that the index and the name
	 * agree. A curve belongs to one material, and two regions of one steel
	 * hold two copies.
	 */
	struct PermeabilityCurve
	{
		std::string material;
		/// H, A/m, from 0 and strictly increasing.
		std::vector< double > fieldStrength;
		/// B, T, from 0, one per H.
		std::vector< double > fluxDensity;
	};

	/// The rule between a curve's points, which is the `.nc`'s
	/// `permeability_curve_b:interpolation`. Piecewise linear in H; the stage
	/// that solves with iron owns the rule, and changing it changes this.
	inline constexpr char const *permeabilityInterpolation = "piecewise_linear";

	/**
	 * THE MATERIAL TABLE: row 0 is vacuum, row 1 plasma, and the conductors
	 * follow in the order they were added -- for coils, `[[coils]]` order.
	 *
	 * **A COIL ROW CARRIES ITS RECTANGLE**, because conductor membership is
	 * geometric (§2.2): meq::Coil::contains() at the solver's points, whatever
	 * the conductor model. The cross-section is the copper whether its current
	 * is meshed, subtracted or carried as filaments -- so this takes the
	 * `[[coils]]` block's rectangle and never the filaments it was split into.
	 */
	class MaterialTable
	{
		public:
			static constexpr int vacuum = 0;
			static constexpr int plasma = 1;

			/// The vacuum and plasma rows, and nothing else.
			MaterialTable();

			/// Append one `[[coils]]` block's row. Its current is the
			/// rectangle's own. @return the row's index.
			/// @param circuit  a row of circuits(), or -1 for none.
			/// @throws std::invalid_argument for a circuit that is not a row.
			int addCoil( Coil const &shape, std::string const &name,
			             int circuit = -1 );

			/// Append a row with no shape -- a RESERVED kind, which nothing
			/// produces yet but a reader must be able to hold. Such a row is
			/// never matched by conductorContaining().
			///
			/// A PASSIVE row has no circuit and no curve, and a relative
			/// permeability of 1. A FERROMAGNETIC row holds EITHER a constant
			/// relative permeability (finite and positive, no curve) OR a
			/// curve, whose own column is then NaN, the fill -- never both --
			/// and the curve must name this row.
			/// @throws std::invalid_argument for Vacuum, Plasma or Coil, which
			///         have their own constructors above, and for a row
			///         breaking the rules just stated.
			int addMaterial( Material const &row );

			/// Append a circuit row. @return its index.
			/// @throws std::invalid_argument for an empty name or one already
			///         in the table.
			int addCircuit( Circuit const &row );
			std::vector< Circuit > const &circuits() const;

			/// Append a permeability curve, before the row that points to it.
			/// @return its index.
			/// @throws std::invalid_argument unless it has a material name and
			///         at least two points, starts at ( 0, 0 ), and is
			///         strictly increasing in H with B finite and
			///         non-decreasing.
			int addPermeabilityCurve( PermeabilityCurve const &curve );
			std::vector< PermeabilityCurve > const &permeabilityCurves() const;

			int size() const;

			/// @throws std::out_of_range naming the index and the size.
			Material const &material( int index ) const;
			std::vector< Material > const &materials() const;

			/// The coil rows, in the order they were added.
			std::size_t coilCount() const;
			Coil const &coil( std::size_t which ) const;
			int coilMaterial( std::size_t which ) const;

			/// The row of the FIRST conductor containing ( R, z ), or -1. A
			/// coil's set is closed (meq::Coil::contains), so a point on an
			/// edge two coils share goes to the earlier one.
			int conductorContaining( double radius, double z ) const;

		private:
			std::vector< Material > rows;
			std::vector< Coil > coilShapes;
			std::vector< int > coilRows;
			std::vector< Circuit > circuitRows;
			std::vector< PermeabilityCurve > curves;
	};

	/**
	 * THE TABLE A RUN'S `[[coils]]` AND `[[circuits]]` MAKE: one coil row per
	 * block, in block order, and one circuit row per `[[circuits]]` block, in
	 * block order, followed by `default` when some coil is unwired -- so a
	 * file whose coils are all wired has no empty group.
	 *
	 * @param coils     the configuration's coil set.
	 * @param currents  the SOLVED current of each coil, in block order, when it
	 *                  differs from the file's -- an inverse solve moves them.
	 *                  Empty takes the file's. A supplied circuit's current is
	 *                  its first member's current over that member's
	 *                  multiplier, which is how the driver seeds an inverse.
	 * @throws std::invalid_argument if @a currents is neither empty nor one
	 *         per coil.
	 */
	MaterialTable buildMaterialTable( CoilConfig const &coils,
	                                  std::vector< double > const &currents = {} );

	/// One non-zero entry of an element's record.
	struct MaterialFraction
	{
		int material = 0;
		double fraction = 0.0;
	};

	/// The region an element's fractions derive, by §2.1's table.
	/// @param entries the element's non-zero fractions.
	ElementRegion deriveRegion( std::vector< MaterialFraction > const &entries,
	                            MaterialTable const &table );

	/**
	 * EVERY ELEMENT'S MATERIALS: the ragged record of §2.3 and what is derived
	 * from it. Elements are in mesh order and the offsets are a prefix sum, so
	 * element `e`'s entries are `entries[ offsets[ e ] .. offsets[ e + 1 ] )`,
	 * sorted by material.
	 */
	struct ElementMaterials
	{
		std::vector< int > offsets{ 0 };
		std::vector< MaterialFraction > entries;

		/// Dense, one per element.
		std::vector< ElementRegion > region;
		std::vector< double > plasmaFraction;
		/// `sum of ip.weight * |J|` over the element's points: its area.
		std::vector< double > area;
		/// Whether a point of the element was plasma and conductor at once.
		std::vector< char > faulted;
		/// `[source] ExcludeAttributes` named the element's attribute.
		std::vector< char > excluded;
		std::vector< int > attribute;

		/// Per material-table row: `sum of ip.weight * |J|` over the points
		/// assigned to it.
		std::vector< double > materialArea;

		/// Elements with a faulted point, which is `region_faults`.
		int faults = 0;

		/// One weight per material present, the weights being the
		/// quadrature weights of the element's points assigned to it.
		struct Weight
		{
			int material = 0;
			double weight = 0.0;
		};

		/**
		 * APPEND THE NEXT ELEMENT. The fractions are each weight over their
		 * sum, so an element holding one material reads exactly 1 and the
		 * fractions sum to 1 to round-off; entries are merged by material and
		 * sorted.
		 *
		 * @throws std::invalid_argument for an empty list, a material outside
		 *         the table, or a weight that is negative or not finite, or a
		 *         sum that is not positive.
		 */
		void addElement( MaterialTable const &table, std::vector< Weight > weights,
		                 bool faultedIn, bool excludedIn, int attributeIn );

		int elementCount() const;

		/// The fraction of @a material in @a element, zero if absent.
		double fraction( int element, int material ) const;

		/**
		 * The conductor holding the largest share of @a element -- a coil
		 * row today, and a passive or ferromagnetic row when those exist --
		 * with its fraction; `{ -1, 0 }` where no conductor is present. A tie
		 * goes to the lower row, so the answer does not depend on the order
		 * the entries were stored in.
		 *
		 * The VTK's `conductor_material` and `conductor_fraction`: one
		 * conductor per cell is what a cell array can carry, and the full
		 * ragged record stays in the `.nc` (`ELEMENT-REGIONS-PLAN.md` §2.3).
		 */
		MaterialFraction dominantConductor( int element, MaterialTable const &table ) const;

		/// Elements whose derived region is @a which.
		int regionCount( ElementRegion which ) const;
	};
}

#endif // MEQ_MATERIALS_HPP
