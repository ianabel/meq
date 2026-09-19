#ifndef MEQ_CONFIG_HPP
#define MEQ_CONFIG_HPP

/*
 * Configuration for the fixed-boundary Grad-Shafranov solver.
 *
 * MEQ solves
 *
 *     -div_bar( (1/r) grad_bar( psi ) ) = F( r, z, psi ) / r    in Omega,
 *                                   psi = psi_D                 on Gamma,
 *
 * by HDG, with the nonlinearity in F handled by Newton. A run is described by a
 * TOML file with eight tables. THIS COMMENT SAID SIX FOR A LONG TIME, and the
 * two it omitted are the two that arrived last:
 *
 *     [mesh]            the background box and its subdivision
 *     [discretisation]  polynomial degree and the HDG stabilisation tau
 *     [source]          which F, and its parameters
 *     [boundary]        the Dirichlet data psi_D
 *     [solver]          Newton controls
 *     [output]          where the mesh and grid functions are written
 *     [initialguess]    where Newton starts
 *     [adaptivity]      the refinement loop, when there is one
 *
 * and three nested forms beneath those:
 *
 *     [mesh.generate]   how to MAKE the mesh named above, when it is this
 *                       file's own build product rather than an input
 *     [boundary.shape]  the analytic curve Gamma, on the extension path
 *     [[source.species]] an ARRAY of tables, one per species, on the rotating
 *                       path -- the only array of tables in the schema
 *
 * [mesh], [discretisation] and [source] are required; the rest are optional and
 * every key in them has a documented default. See examples/soloviev-nstx.toml
 * and examples/manufactured.toml for annotated, working files.
 *
 * References below are to the two papers in refs/ (see refs/Refs.md):
 *   refs/HDG-GradShafranov.pdf           Sanchez-Vizuet & Solano, CPC 235 (2019)
 *   refs/HDG-GradShafranov-Adaptive.pdf  ... & Cerfon, CPC 255 (2020)
 *
 * NAMING. Two conventions meet in this file and they deliberately disagree:
 * TOML key names are UpperCamelCase (RMin, PolynomialDegree, GGPrimeFile), as
 * they are across the sibling projects, while C++ identifiers -- including the
 * members that hold those keys' values -- are lowerCamelCase (rMin,
 * polynomialDegree, ggPrimeFile). So the string literals in Config.cpp are
 * capitalised and the members they are read into are not. That is intended;
 * please do not "fix" either side into the other.
 *
 * SCHEMA EVOLUTION. Unknown keys are rejected, with a did-you-mean naming the
 * nearest accepted key, because a silently ignored key is how a configuration
 * format grows two readers with different defaults and never says so. The
 * corollary, when a key is renamed later: keep the old spelling as a deprecated
 * alias that warns and still works, rather than deleting it. There are no
 * aliases yet -- this schema is new, and shares no key with the free-boundary
 * configuration it replaces -- so the mechanism is a TODO rather than dead code
 * nothing exercises. It belongs in the Table reader in Config.cpp, beside
 * rejectUnknownKeys().
 *
 * Layering note: this class parses and validates, and nothing more. It does not
 * build a meq::Source, an mfem::Mesh or a solver -- it reports what the file
 * said, as a small tree of plain structs, and leaves construction to a factory
 * that knows about the numerics. That keeps Config free of MFEM and of the
 * source hierarchy, so the configuration can be unit tested on its own (and so
 * that a Source interface still in flux does not force a change here).
 *
 * Every quantity carrying units is documented with them below. Lengths are
 * metres unless a particular benchmark is posed in normalised units, in which
 * case the example file says so.
 */

#include <stdexcept>
#include <string>
#include <variant>

#include <toml.hpp>

namespace meq
{

	// The single exception type thrown for every configuration failure: an
	// unreadable file, a TOML syntax error, a missing or misspelt key, a value
	// of the wrong type, or a value that parses but cannot describe a run
	// (RMax below RMin, a negative polynomial degree, ...).
	//
	// what() always names the file, and names the offending key when there is
	// one. getKey() is the fully qualified key ("mesh.RMin"), or empty for
	// errors that are not about one particular key.
	class ConfigError : public std::runtime_error
	{
		public:
			ConfigError( std::string const & file, std::string const & key, std::string const & message );

			std::string const & getFile() const noexcept { return fileName; };
			std::string const & getKey() const noexcept { return keyName; };

		private:
			std::string fileName;
			std::string keyName;
	};

	/**
	 * `[mesh.generate]` -- THE MESH AS A BUILD PRODUCT OF THIS FILE.
	 *
	 * Everything else in `[mesh]` reads a mesh somebody already made. A free
	 * boundary run cannot be made by `mfem::Mesh::MakeCartesian2D` at all --
	 * it needs a semicircle reaching `r = 0` exactly, with the conductors
	 * fragmented in -- so it comes from `tools/mesh/halfdisc.py`. Without this
	 * block the geometry is stated TWICE: once on that script's command line
	 * and once in the TOML the solve reads, with nothing checking that the two
	 * agree. `examples/diverted-tokamak.toml`'s header is what that looks like.
	 *
	 * **THE COILS ARE NOT REPEATED HERE, AND THAT IS THE POINT.** The
	 * generator's `--coil` rectangles are derived from the `[[coils]]` blocks,
	 * in file order, which is the order `tools/mesh/halfdisc.py` documents its
	 * `10 + i` element attributes in. So a machine's conductors are written
	 * once and the mesh is aligned to the conductors the solve will actually
	 * integrate over. FB-2 measured what that alignment is worth: rates of
	 * **1.99 / 2.88 / 3.01 aligned against 1.33 / 1.27 / 1.09 with the
	 * conductor cut by elements**, so a coil the mesh does not know about
	 * costs a full order, silently.
	 *
	 * **MEQ DOES NOT RUN THE GENERATOR.** `meq` is linked against MFEM and not
	 * against gmsh, deliberately -- `tools/mesh/README.md` records why -- so
	 * what this block buys inside the driver is a REFUSAL and a QUERY:
	 * `meq --mesh-command` prints the generator's argument list from these
	 * values, and `meq` will not solve a configuration carrying this block
	 * unless the caller passes `--mesh-ready` to say the mesh has been made
	 * from it. `meq-run` is the caller that does both.
	 *
	 * **THE UNITS AND THE CONVENTIONS ARE MEQ'S, NOT THE SCRIPT'S.** The box
	 * is `RMin`/`RMax`/`ZMin`/`ZMax` as `[mesh]`'s own box is, where
	 * `halfdisc.py --plasma` takes a corner and two extents; the driver
	 * converts. A schema that exposed the tool's convention would put two
	 * meanings of four numbers in one file.
	 */
	struct MeshGeneratorConfig
	{
		/// Whether the file asked for one. False is every configuration that
		/// names a mesh made somewhere else, which is all of them but one.
		bool given = false;

		/// `Tool` -- which generator. "halfdisc" is the only one there is;
		/// naming it is what lets a second arrive without a second table.
		std::string tool;

		/// `Radius` -- the DISC's radius, metres, strictly positive. The
		/// semicircle is centred on the origin and reaches `r = 0` exactly,
		/// because that is the geometry the exterior expansion of
		/// `meq::ExteriorDtN` is a statement about.
		///
		/// **IT IS NOT GAMMA, AND READING IT AS GAMMA IS THE EASY MISTAKE.**
		/// `D_h` is cut FROM this mesh as the elements lying inside
		/// `[boundary.exterior] Radius`, so the arc gmsh draws is the
		/// BACKGROUND's outer edge and Gamma is the smaller semicircle within
		/// it, with the band between the resulting staircase and Gamma bridged
		/// by the Cockburn-Solano transfer. Gamma must therefore fit STRICTLY
		/// inside this, which the `[boundary.exterior]` parse checks -- before
		/// gmsh has run at all.
		double radius = 0.0;

		/// `Size` -- the background element size, metres, strictly positive.
		double size = 0.0;

		/// `Order` -- the geometric order, >= 1. Above 1 the arc is curved:
		/// the mid-edge nodes are placed on the true circle rather than on the
		/// chord, and the axis stays exact.
		int order = 1;

		/// `CoilSize` -- the element size inside the conductors, metres.
		/// Zero means "the background `Size`", which is the script's own
		/// default and not a special case here.
		double coilSize = 0.0;

		/// `PlasmaRMin` and friends -- a box to refine inside, where the
		/// plasma is expected, with `PlasmaSize` the size in it. Gamma has to
		/// stand far enough out for the exterior expansion to converge and the
		/// plasma occupies a small part of what that encloses, so meshing the
		/// whole disc at the plasma's resolution spends most of the elements on
		/// vacuum. Measured, that is 2.6x at `Radius = 1.5` and 9.0x at 3.0.
		bool plasmaGiven = false;
		double plasmaRMin = 0.0;
		double plasmaRMax = 0.0;
		double plasmaZMin = 0.0;
		double plasmaZMax = 0.0;
		double plasmaSize = 0.0;

		/// `LimiterR`, `LimiterZ`, `LimiterRadius` -- a circular limiter
		/// FRAGMENTED IN rather than cut, so its edges are mesh faces, and
		/// written as element attribute 20. That is what
		/// `[boundary.limiter] SurfaceAttribute` reads, and the two are checked
		/// against each other: a `SurfaceAttribute` naming a region this mesh
		/// will not carry is refused rather than discovered at the solve.
		bool limiterGiven = false;
		double limiterR = 0.0;
		double limiterZ = 0.0;
		double limiterRadius = 0.0;

		/**
		 * `Vessel` -- A CLOSED POLYGON, ALTERNATING `R` AND `Z` IN METRES,
		 * FRAGMENTED IN. Everything inside `Gamma`, outside this, and not a
		 * conductor takes element attribute **30**.
		 *
		 * It exists so that `[source] ExcludeAttributes = [ 30 ]` can say *the
		 * plasma can never be here* ONCE, at mesh time. The support is
		 * re-decided on every residual evaluation and a vessel does not move, so
		 * a per-element geometric test would be paying repeatedly for an answer
		 * that cannot change.
		 *
		 * **IT IS NOT CONFINEMENT.** `psi_bnd` confines the plasma -- `F` is
		 * zero wherever `Psi <= 0` and setting `psi_bnd` IS the confinement --
		 * and meq::PlasmaComponent separates the lobes of `{ Psi > 0 }` that a
		 * level set leaves joined. This covers the one case neither can:
		 * several O-points across a saddle, where connectivity cannot say which
		 * is the plasma.
		 *
		 * Empty is the default. At least three points, so at least six numbers,
		 * every `R >= 0`, and a non-zero area.
		 */
		std::vector<double> vessel;

		/// `Transition` -- the width of the graded transition out of a refined
		/// region, metres. Zero means the script's default of four background
		/// sizes.
		double transition = 0.0;

		/// `Check` -- re-read the written file and assert MEQ's preconditions
		/// on it: that `r` reaches 0 EXACTLY, that Gamma and the axis are the
		/// outer boundary and nothing else, and that each coil attribute covers
		/// its rectangle. **It defaults to TRUE here and to false on the
		/// script's own command line**, which is a deliberate disagreement: a
		/// human meshing interactively reads the printed report, and a mesh
		/// generated inside a run has nobody looking at it.
		bool check = true;

		/**
		 * `Symmetric` -- mesh `z >= 0` and REFLECT it, so the mesh is exactly
		 * mirror-symmetric about `z = 0`.
		 *
		 * **THIS IS WHAT MAKES `[solver] UpDownSymmetry` REACHABLE AT ALL.**
		 * That projection averages every dof with the dof at its own
		 * reflection, so it needs a mesh whose dofs are mirror-paired, and
		 * gmsh's triangulation of a symmetric geometry is NOT symmetric --
		 * it picks a diagonal and it picks freely. Measured on MAST-U's
		 * committed mesh, 3588 of 4735 vertices have no partner.
		 *
		 * The generator REFUSES this on a geometry that is not itself
		 * mirror-symmetric -- an unpaired conductor, a limiter off the
		 * midplane, a vessel outline that is not its own image -- rather than
		 * reflecting a machine into a different machine. `Plasma*` is the
		 * exception and is allowed to be asymmetric: it is a size field and
		 * not geometry, so the half that is meshed is refined as asked and the
		 * other half gets the reflection.
		 *
		 * Defaults false, which is the mesh every existing example has.
		 */
		bool symmetric = false;
	};

	/// The element attribute `tools/mesh/halfdisc.py` writes the interior of a
	/// `--limiter` circle as. Named here because `[boundary.limiter]
	/// SurfaceAttribute` is checked against it, and a bare 20 in a comparison
	/// says nothing about where it came from.
	inline constexpr int generatedLimiterAttribute = 20;

	// [mesh]
	//
	// Following HDG-GradShafranov-Adaptive.pdf section 2.1, the computational
	// domain is cut out of a uniform, shape-regular background triangulation of
	// a box containing the plasma region. The box is deliberately NOT fitted to
	// the plasma boundary, so these bounds should enclose it with a margin.
	struct MeshConfig
	{
		// Background box [RMin,RMax] x [ZMin,ZMax], in metres. RMin >= 0: the
		// Grad-Shafranov operator carries a 1/r, so a box reaching r = 0
		// contains the coordinate singularity. That is allowed, but is rarely
		// what is wanted.
		double rMin = 0.0;
		double rMax = 0.0;
		double zMin = 0.0;
		double zMax = 0.0;

		// NR, NZ: cells across the box in each direction before refinement
		// (dimensionless counts, >= 1). Each cell is split into triangles, so
		// the initial mesh diameter is the cell diagonal,
		// h = sqrt( ((RMax-RMin)/NR)^2 + ((ZMax-ZMin)/NZ)^2 ).
		int nR = 1;
		int nZ = 1;

		// RefinementLevels: levels of uniform refinement applied to the
		// background mesh (dimensionless count, >= 0). Each level halves h.
		int refinementLevels = 0;

		// File: an optional alternative to the box -- a mesh in any format MFEM
		// reads. Empty (the usual case) means "generate the box above". When it
		// is set the box bounds are not required and NR/NZ are ignored;
		// RefinementLevels still applies.
		std::string file;

		bool fromFile() const { return !file.empty(); };

		/// `[mesh.generate]` -- how to MAKE the file named above, when the file
		/// is this configuration's own build product. See MeshGeneratorConfig.
		MeshGeneratorConfig generate;
	};

	// [discretisation]
	struct DiscretisationConfig
	{
		// PolynomialDegree: the degree k of the HDG spaces (dimensionless,
		// >= 0). Both papers report results for k = 1 ... 5.
		int polynomialDegree = 1;

		// Tau: the HDG stabilisation. Dimensionless, > 0. Both papers take
		// tau = 1 and note that optimal convergence needs only tau = O(1), so
		// 1.0 is the default and there is rarely a reason to change it. (The
		// code this replaces used 5.0, with no justification recorded.)
		double tau = 1.0;
	};

	// [source] -- which right-hand side F( r, z, psi ), and its parameters.
	enum class SourceType
	{
		Soloviev,      // "soloviev"
		MHD,           // "mhd"
		Manufactured,  // "manufactured"
		Rotating       // "rotating"
	};

	// F( r, z, psi ) = -( (1 - A) r^2 + A ), independent of psi and hence
	// linear. See HDG-GradShafranov.pdf eq (10), the NSTX case of
	// HDG-GradShafranov-Adaptive.pdf section 4.1, and Cerfon & Freidberg,
	// Phys. Plasmas 17, 032502 (2010) for the geometry.
	struct SolovievParameters
	{
		// A: dimensionless, and with the flux normalised so that A + C = 1 it
		// is the only free parameter; it fixes the plasma pressure against the
		// magnetic pressure. The NSTX benchmark uses A = -0.52.
		double a = 0.0;
	};

	// Static MHD equilibrium: F is built from the two tabulated flux functions
	//
	//     F( r, z, psi ) = mu0 r^2 p'(psi) + (g g')(psi)
	//
	// (meq::MHDSource, which documents the conventions the tables must follow).
	// Both files are read by meq::SplineProfile, which documents the format:
	// one knot per line, `psi f(psi) f'(psi)`, '#' comments. Relative paths are
	// taken as given, i.e. resolved against the working directory of the run,
	// not against the directory holding the configuration file.
	struct MHDParameters
	{
		// PPrimeFile: path to the tabulated dp/dpsi  [Pa / (Wb per radian)]
		std::string pPrimeFile;
		// SafetyFactorFile: path to a tabulated TARGET q( Psi ), which makes
		// gg' an OUTPUT of the run rather than an input and puts an outer
		// Newton around the whole solve. ROADMAP.md item 10.
		//
		// IT IS AN ALTERNATIVE TO GGPrimeFile AND NAMING BOTH IS REFUSED.
		// The two say opposite things about the same quantity -- one prescribes
		// the toroidal field, the other asks for whatever field delivers a
		// given q -- so a precedence rule would decide which physics a run did
		// on the strength of key order, and both spellings converge.
		//
		// THE TABLE IS IN THE SOURCE'S Psi, one on the axis, like every other
		// profile table in examples/. meq::ToroidalFieldMap owns the
		// reflection to the family's Psi_N; see SafetyFactor.hpp section 2 on
		// why getting it backwards converges to a reversed shear rather than
		// failing.
		std::string safetyFactorFile;

		// SafetyFactorDegree: the degree of the g^2 polynomial the loop
		// iterates on. The fit is what makes the outer Newton affordable -- the
		// unknown is d + 1 coefficients rather than nFieldDOF -- and what stops
		// a noisy inversion near the axis becoming SHAPE in gg'.
		unsigned int safetyFactorDegree = 2;

		// ToroidalFieldGuess: g = R B_phi to start from, as a constant, so the
		// loop opens at gg' = 0. Required with SafetyFactorFile and has no
		// default: q determines g through the geometry, but only once there IS
		// a geometry, and a machine's vacuum R0 B0 is the number a user has.
		double toroidalFieldGuess = 0.0;

		// GGPrimeFile: path to the tabulated g dg/dpsi
		//                                     [T^2 m^2 / (Wb per radian)]
		std::string ggPrimeFile;
		// PPrimeScale, GGPrimeScale: constants multiplying the tables as read.
		// A table arrives in whatever units its author wrote it in, and that is
		// frequently not MEQ's; editing the file would make the file a function
		// of which code reads it. See meq::ScaledProfile.
		double pPrimeScale = 1.0;
		double ggPrimeScale = 1.0;
		// Mu0: the vacuum permeability multiplying the r^2 p' term [H/m]. The
		// SI value by default; set it to 1 for a problem posed in normalised
		// units.
		double mu0 = 4.0e-7*3.14159265358979323846;
		// Normalised: whether the two tables are functions of NORMALISED flux,
		// Psi = psi/psi_ax, rather than of psi itself. That makes psi_ax a
		// functional of the solution and therefore an UNKNOWN of the non-linear
		// system, which the solver closes by a bordered Newton -- see
		// meq::NormalisedSource and CLAUDE_HDGGS.md, "Newton, and the obligation it
		// creates". False by default, because it changes what the tables mean.
		bool normalised = false;
		// PsiAxis: the starting value of psi_ax [Wb per radian]. REQUIRED when
		// Normalised is true and refused otherwise. It is a guess in the Newton
		// sense and not a scale factor: the iteration has to start inside a
		// basin, and at a fixed psi_ax the equation generally has a second,
		// non-physical solution the iteration can reach instead.
		double psiAxis = 0.0;
		// ConfineToPlasma: whether F and dF/dpsi are ZERO wherever the
		// normalised flux is non-positive, so that the plasma's SUPPORT moves
		// with the solution instead of being the whole domain. That is what a
		// free-boundary run needs and what a fixed-boundary run must not have:
		// there the domain IS the plasma, Psi > 0 throughout by construction,
		// and switching this on would only put a branch in the inner loop.
		//
		// REQUIRES Normalised = true, because the test is on Psi and psi_bnd
		// and psi_ax are what define it. Refused otherwise rather than ignored.
		//
		// AND IT REQUIRES A PROFILE THAT VANISHES AT THE EDGE, p'( 0 ) = 0.
		// With p'( 0 ) != 0 the source JUMPS across the plasma boundary, the
		// assembled residual is discontinuous in the unknowns, and there is no
		// Jacobian to iterate with -- measured, Newton fails at every degree
		// and every mesh, from the exact solution, and under PicardThenNewton.
		// Nothing here can check that, the profile being a table; see
		// meq::NormalisedSource::setPlasmaSupport.
		bool confineToPlasma = false;

		// [source] PlasmaCurrent, in AMPERES and signed. Non-zero makes the
		// profile SCALE an unknown of the bordered Newton and prescribes the
		// current instead: the profiles then give the current's SHAPE and this
		// gives its SIZE. Zero, the default, fixes the amplitude as before.
		//
		// AMPERES HERE AND mu0 I_p AT THE LIBRARY, WHICH IS A DELIBERATE
		// DIFFERENCE. meq::GradShafranovSolver::setPlasmaCurrent takes mu0 I_p
		// because everything in the solver already speaks in it -- Ampere's law
		// reads the flux integral as -mu0 I_p and the constraint is assembled as
		// int F/r, which IS mu0 I_p -- and taking amperes THERE would mean the
		// solver knowing a mu0, which could disagree with the source's own and
		// scale two terms of one equation differently. The CONFIGURATION layer
		// has no such problem: the file names exactly one mu0, under [source],
		// and the driver multiplies by it. That is the same rule [[coils]]
		// already follows for its Current, and it is why there is no Mu0 key
		// anywhere but [source].
		double plasmaCurrent = 0.0;
	};

	// The nonlinear manufactured solution of HDG-GradShafranov.pdf Example 5,
	//
	//     psi = sin( Kr ( r + R0 ) ) cos( Kz z ),
	//
	// with F chosen so that psi solves the equation. Its Dirichlet data is not
	// zero, so a run using this source normally sets [boundary] Type = "exact".
	struct ManufacturedParameters
	{
		// R0: the radial offset r0 in the expression above, in metres. NOT the
		// major radius. Example 5 uses r0 = -0.5.
		double r0 = 0.0;
		// Kr: radial wavenumber, in radians per metre. Example 5: 1.15 pi.
		double kr = 0.0;
		// Kz: vertical wavenumber, in radians per metre. Example 5: 1.15.
		double kz = 0.0;
	};

	// One species of a rotating plasma. Each profile is given EITHER as a
	// constant, in the scalar key, OR as a path to a table in the *File key --
	// exactly one of the two. Two keys rather than one that changes meaning by
	// node type, for the reason recorded at the top of Config.cpp: TOML
	// distinguishes 1 from 1.0 and a type-dispatched key would inherit that trap.
	struct SpeciesParameters
	{
		// Name: what this species is called in the output file's variable names.
		// Defaults to "species<i>" if omitted.
		std::string name;
		// Mass: particle mass [kg]. Must be positive.
		double mass = 0.0;
		// Charge: Z_s, signed and DIMENSIONLESS -- +1 for a proton, -1 for an
		// electron, +6 for fully stripped carbon. Not a charge in coulombs.
		double charge = 0.0;
		// Temperature / TemperatureFile: T_s [JOULES, not eV and not keV].
		// TemperatureScale multiplies whichever was given -- so a table in keV
		// becomes Joules with TemperatureScale = 1.602176634e-16 rather than by
		// rewriting the file.
		double temperature = 0.0;
		std::string temperatureFile;
		double temperatureScale = 1.0;
		// Density / DensityFile: n_s0 [m^-3], the density of this species ON
		// THE CURVE r = ReferenceRadius. Both are absent when Neutralising, and
		// so is DensityScale.
		double density = 0.0;
		std::string densityFile;
		double densityScale = 1.0;
		// Neutralising: whether this species' density is DERIVED from the
		// others by charge neutrality rather than given. Exactly one species
		// must set it. Fixing the gauge removes exactly one function's worth of
		// freedom from the densities, so for n species there are n - 1
		// independent ones; this is where that missing one comes from, and
		// asking for n profiles that happen to balance invites n that do not.
		bool neutralising = false;
	};

	// A plasma in sonic toroidal rotation: refs/RotatingGK.pdf eq (136), closed
	// by its (96) for the poloidal density variation and (97) for the
	// electrostatic potential phi_0 that holds quasineutrality against it. The
	// density is NOT a flux function -- centrifugal force sweeps heavy species
	// outboard -- which is the whole content of this source. See
	// meq::RotatingSource, and docs/rotation.rst for the derivation.
	struct RotatingParameters
	{
		// The species, [[source.species]]. Between two and meq::maxSpecies, of
		// both charge signs, exactly one of them Neutralising.
		std::vector<SpeciesParameters> species;
		// Omega / OmegaFile: the rigid rotation frequency omega(psi) [rad/s].
		// BOTH ABSENT MEANS NO ROTATION, and the source then reduces to the
		// static equation meq::MHDSource solves.
		double omega = 0.0;
		std::string omegaFile;
		double omegaScale = 1.0;
		bool omegaGiven = false;
		// GGPrime / GGPrimeFile: the single product g dg/dpsi
		//                                      [T^2 m^2 / (Wb per radian)],
		// exactly as the "mhd" source takes it.
		double ggPrime = 0.0;
		std::string ggPrimeFile;
		double ggPrimeScale = 1.0;
		// ReferenceRadius: rRef [m]. THE GAUGE. phi_0 vanishes on this curve,
		// which is what makes each Density the physical density there, and it
		// is a CONSTANT -- the geometric axis -- not the magnetic axis and not
		// a flux-surface average. Two sets of densities differing by the gauge
		// describe the same plasma, so a comparison against another code has to
		// agree on this first.
		double referenceRadius = 1.0;
		// Mu0, Normalised, PsiAxis, ConfineToPlasma: as MHDParameters, and
		// meaning the same.
		double mu0 = 4.0e-7*3.14159265358979323846;
		bool normalised = false;
		double psiAxis = 0.0;
		bool confineToPlasma = false;
		double plasmaCurrent = 0.0;
	};

	using SourceParameters = std::variant< SolovievParameters, MHDParameters,
	                                       ManufacturedParameters, RotatingParameters >;

	struct SourceConfig
	{
		SourceType type = SourceType::Soloviev;
		SourceParameters parameters = SolovievParameters{};

		/**
		 * `[source] ExcludeAttributes` -- MESH ELEMENT ATTRIBUTES THAT CAN NEVER
		 * BE PLASMA, WHATEVER THE FLUX SAYS THERE.
		 *
		 * **A STATEMENT ABOUT THE DEVICE AND NOT ABOUT THE SOLUTION**, which is
		 * why it lives beside `parameters` rather than inside one of them: it
		 * means the same thing for every source type that has a plasma support
		 * at all. `psi_bnd` confines the plasma -- `F` is zero wherever
		 * `Psi <= 0`, and setting `psi_bnd` IS the confinement -- and
		 * meq::PlasmaComponent separates the lobes of `{ Psi > 0 }` that a level
		 * set leaves connected. **Neither can know that a lobe is behind a
		 * wall.** That is the gap, and it is a real one only where several
		 * O-points sit across a saddle, since connectivity handles every case
		 * where the lobes are genuinely disjoint.
		 *
		 * AN ATTRIBUTE AND NOT A POLYGON: the support is re-decided on every
		 * residual evaluation and this answer cannot move, so the mesher says it
		 * once and the solve reads a lookup.
		 *
		 * Empty is the default and changes nothing. Refused unless
		 * `Normalised = true`, for ConfineToPlasma's reason.
		 */
		std::vector<int> excludeAttributes;

		// Typed access. Each throws ConfigError if the configured type is not
		// the matching one, so a factory that has already switched on type()
		// can use them without a second check.
		SolovievParameters const & getSoloviev() const;
		MHDParameters const & getMHD() const;
		ManufacturedParameters const & getManufactured() const;
		RotatingParameters const & getRotating() const;

		/// Whether the profiles are functions of normalised flux, so that
		/// psi_ax is an unknown and the solver must be handed the source
		/// through setSource( NormalisedSource &, double ) rather than through
		/// setSource( Source const & ). Exposed so that a driver can branch
		/// without duplicating the switch over type.
		bool isNormalised() const;

		/// The starting psi_ax, meaningful only when isNormalised().
		double psiAxisGuess() const;

		/// Whether the source is confined to `{ Psi > 0 }`, so that the
		/// plasma's support moves with the solution. False for every source
		/// that is not normalised, since the test is on Psi.
		bool confinesToPlasma() const;

		/// `[source] PlasmaCurrent` in AMPERES, signed, or zero if the file did
		/// not prescribe one. Zero means the profile amplitude is fixed and the
		/// current is whatever it comes out as; non-zero makes the amplitude an
		/// unknown of the bordered Newton. See MHDParameters::plasmaCurrent for
		/// why this is amperes where the library's setPlasmaCurrent is `mu0 I_p`.
		double plasmaCurrent() const;

		/// The permeability this source multiplies its pressure term by.
		///
		/// **EXPOSED SO THE COILS CAN SHARE IT, AND THAT IS THE WHOLE REASON.**
		/// A coil's contribution is `mu0 r I/|Omega_c|` and the plasma's is
		/// `mu0 r^2 p' + g g'`; they are ADDED, so a run in normalised units
		/// that sets `[source] Mu0 = 1` and leaves the coils at the SI value
		/// would be summing two terms scaled a million-fold apart -- and it
		/// would converge, at full order, to a machine nobody described. There
		/// is deliberately no `Mu0` key on `[[coils]]` for that reason: two
		/// keys that must agree are a way of writing down a disagreement.
		///
		/// The benchmark sources (soloviev, manufactured) carry no mu0 of
		/// their own -- it is folded into their coefficients -- and answer with
		/// meq::vacuumPermeability, which is what SI coil currents want.
		double permeability() const;
	};

	// [boundary] -- the Dirichlet data psi_D on Gamma.
	enum class BoundaryDataType
	{
		Zero,   // "zero":  psi = 0, the fixed-boundary problem proper
		Exact   // "exact": psi = the source's known exact solution, for
		        //          convergence studies. Only a source that has one
		        //          (soloviev, manufactured) may ask for this.
	};

	// [boundary.shape] -- the curve Gamma, when it is not the mesh boundary.
	enum class ShapeType
	{
		None,     // "none":   Gamma is the mesh boundary; the fitted path
		Miller,   // "miller": the D-shape, in physical parameters
		Mxh       // "mxh":    the general Miller extended harmonic surface
	};

	// Parameters as written; meq::BoundaryShape is what validates them, since
	// star-shapedness is a property of the curve rather than of any one number.
	// Miller is converted to MXH on construction -- see BoundaryShape.hpp.
	/**
	 * ONE RECTANGULAR-CROSS-SECTION COIL CARRYING A UNIFORM CURRENT DENSITY.
	 *
	 * FB-6, and section 5.4 of FREE-BOUNDARY-PLAN.md calls the coils "ordinary"
	 * -- coil currents are data, and the source adds
	 * `F_coil = mu0 r I_k / |Omega_ck|` on each coil subdomain. This is the
	 * configuration half of that; `meq::Coil` is the library half and has been
	 * measured since FB-2.
	 *
	 * **THE CURRENT MAY BE GIVEN EITHER WAY AND EXACTLY ONE MUST BE.** `Current`
	 * is the TOTAL through the cross-section in amperes; `CurrentDensity` is the
	 * uniform `j_phi` in A/m^2, and the two are related by the area
	 * `4 * HalfWidth * HalfHeight`. Naming both is refused rather than resolved
	 * by precedence: an author who writes both has two numbers in mind and
	 * silently honouring one of them is how a coil set ends up carrying a
	 * current nobody chose. Naming neither is refused for the same reason.
	 *
	 * The parse resolves whichever was given to a total current, because that is
	 * what meq::Coil takes; `densityGiven` records which the author wrote, so a
	 * diagnostic can quote it back in the units it arrived in.
	 */
	struct CoilParameters
	{
		/// For diagnostics. Defaults to "coil<i>" if the file does not name it.
		std::string name;

		/// The centre, in metres.
		double centreR = 0.0;
		double centreZ = 0.0;

		/// Half-extents, in metres. Strictly positive, and
		/// `centreR - halfWidth` must be strictly positive too: a coil reaching
		/// the axis is refused, because the operator's 1/r is not integrable
		/// through r = 0. Same refusal meq::Coil and meq::BoundaryShape make.
		double halfWidth = 0.0;
		double halfHeight = 0.0;

		/// The TOTAL current through the cross-section, in amperes. Signed.
		/// Resolved at parse time from whichever of Current / CurrentDensity
		/// the file gave.
		double current = 0.0;

		/// True if the file wrote CurrentDensity rather than Current.
		bool densityGiven = false;
	};

	/// The coil set, as a sequence of `[[coils]]` blocks. Empty is legal and is
	/// what every fixed-boundary configuration has.
	struct CoilConfig
	{
		std::vector< CoilParameters > coils;
	};

	struct ShapeConfig
	{
		ShapeType type = ShapeType::None;

		// R0, Z0, r, kappa: MXH's bounding-box parameters, metres and
		// dimensionless. Required unless Type is "none".
		double majorRadius = 0.0;
		double centreHeight = 0.0;
		double minorRadius = 0.0;
		double elongation = 1.0;

		// Type = "miller" only. Triangularity delta, |delta| < 1, entering as
		// s_1 = arcsin( delta ); squareness zeta, entering as s_2 = -zeta.
		double triangularity = 0.0;
		double squareness = 0.0;

		// Type = "mxh" only. CosCoefficients is c_0, c_1, ... c_N and STARTS AT
		// c_0, the tilt; SinCoefficients is s_1, s_2, ... s_N and starts at s_1,
		// there being no s_0. That asymmetry is real and is checked on load.
		std::vector< double > cosCoefficients;
		std::vector< double > sinCoefficients;
	};

	/**
	 * `[boundary.limiter]` -- the contact point that pins `psi_bnd`, FB-3.
	 *
	 * THE PLASMA EDGE IS WHERE THE PROFILES STOP, AND IN A FREE BOUNDARY IT IS
	 * NOT KNOWN IN ADVANCE. The profiles are functions of
	 * `Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd )`, so `psi_bnd` is a
	 * functional of the solution exactly as `psi_ax` is, and it gets a border
	 * row of its own. This block is the point it is pinned at:
	 * meq::GradShafranovSolver::setBoundaryFluxPoint.
	 *
	 * **IT IS MEANINGLESS WITHOUT A NORMALISATION**, and is refused without one.
	 * `psi_bnd` only enters through `Psi`, so a file naming a limiter on a
	 * source that is not `Normalised = true` has asked for an unknown nothing
	 * reads -- which would converge, at full order, to the equilibrium the file
	 * did not describe.
	 */
	struct LimiterConfig
	{
		/// Whether the file named one. False leaves `psi_bnd` fixed at zero and
		/// the bordered solve is the 1x1 it always was.
		bool given = false;

		/// The contact, in metres. `R` must be strictly positive: the axis is
		/// not a limiter, and the nearest potential dof to a point on `r = 0`
		/// is in an element whose flux mass degenerates.
		double r = 0.0;
		double z = 0.0;

		/**
		 * `SurfaceAttribute` -- THE LIMITER AS A CURVE RATHER THAN A POINT.
		 *
		 * The ELEMENT attribute of the region the limiter encloses, so the
		 * limiter is that region's boundary and `psi_bnd = max psi_h` over it,
		 * with the contact FOUND rather than prescribed:
		 * meq::GradShafranovSolver::setLimiterSurface. Zero when the file did
		 * not name one.
		 *
		 * **IT IS AN ALTERNATIVE TO `R`/`Z` AND NAMING BOTH IS REFUSED.** A
		 * prescribed contact and a found one are different constraints on the
		 * same unknown, and both converge -- to equilibria that differ by the
		 * `O( h )` the point version costs. A file that says the answer twice
		 * has not asked a question with one answer.
		 *
		 * **THE MESH MUST BE FITTED TO THE LIMITER**, which is what makes the
		 * polygon exact rather than approximate: `tools/mesh/halfdisc.py
		 * --limiter` fragments the limiter circle into the geometry, so its
		 * edges are mesh faces and its vertices sit on the true circle -- to
		 * 3.3e-16, measured -- and it writes that region as attribute 20.
		 */
		int surfaceAttribute = 0;
	};

	/**
	 * `[boundary.xpoint]` -- THE X-POINT AS TWO UNKNOWNS OF THE SAME NEWTON,
	 * XP-3, and what `[boundary.limiter]` cannot be for a divertor.
	 *
	 * A LIMITER CONTACT IS HARDWARE AND AN X-POINT IS NOT. `[boundary.limiter]`
	 * pins `psi_bnd` at a PRESCRIBED point, which is exactly right for a
	 * material limiter -- the tile is where the drawings say it is -- and wrong
	 * for a diverted plasma, whose null is a functional of the solution and
	 * moves as Newton moves. This block makes `( r_X, z_X )` unknowns beside
	 * `psi_bnd`, closing
	 *
	 *     q_r( r_X, z_X )             = 0
	 *     q_z( r_X, z_X )             = 0
	 *     psi_bnd - psi_h( r_X, z_X ) = 0
	 *
	 * on the same factorisation as everything else:
	 * meq::GradShafranovSolver::setXPointBoundary, whose documentation carries
	 * the derivation and the two structural facts that make it cheap.
	 *
	 * **`R` AND `Z` ARE A SEED AND NOT AN ANSWER**, which is the whole
	 * difference from the block above. They are the initial value of an unknown,
	 * so a file whose numbers are a few centimetres out describes the same
	 * equilibrium as one whose numbers are exact -- and `meq` reports where the
	 * solve actually put it. Getting them badly wrong selects a different
	 * saddle, since this follows ONE null exactly as the axis constraint follows
	 * one O-point.
	 *
	 * **IT IS AN ALTERNATIVE TO `[boundary.limiter]` AND NAMING BOTH IS
	 * REFUSED**, for that block's own reason: all three routes pin the one
	 * unknown `psi_bnd`, and a precedence rule here would decide which
	 * equilibrium a run reports on the strength of key order.
	 *
	 * **AND IT IS MEANINGLESS WITHOUT A NORMALISATION**, again as the limiter
	 * is: `psi_bnd` enters only through `Psi`.
	 */
	struct XPointConfig
	{
		/// Whether the file named one. False leaves `psi_bnd` to
		/// `[boundary.limiter]` or to zero.
		bool given = false;

		/// Where the X-point is BELIEVED to be, in metres. `R` must be strictly
		/// positive: the symmetry axis carries near-zeros of `q` that are not
		/// X-points, and the flux mass `( r q, v )` degenerates there.
		double r = 0.0;
		double z = 0.0;
	};

	/**
	 * `[boundary.exterior]` -- the exact exterior Dirichlet-to-Neumann map on
	 * `Gamma`, FB-5, and the block that makes a run FREE boundary.
	 *
	 * **THIS BLOCK DEFINES `Gamma` ITSELF, WHICH `[boundary.shape]` CANNOT.**
	 * meq::BoundaryShape refuses a surface reaching `r <= 0` -- rightly, since a
	 * closed plasma surface through the axis has a non-integrable `1/r` on it --
	 * and the artificial boundary this block describes is a SEMICIRCLE CENTRED
	 * ON THE AXIS, whose flat side IS the axis. That is not a degenerate MXH
	 * surface; it is a different object, and the axis half of it is ordinary
	 * fitted boundary needing no transfer at all. So the two blocks are
	 * alternatives rather than layers, and naming both is refused.
	 *
	 * **THE SEMICIRCLE IS A REQUIREMENT AND NOT A CONVENIENCE.** meq::ExteriorDtN
	 * is diagonal because the Gegenbauer separation holds on a semicircle about
	 * the axis and nowhere else; on any other curve the exterior map is a dense
	 * boundary-integral operator and this class does not represent it. The
	 * `[mesh]` box must therefore reach `r = 0` exactly, which the driver checks
	 * rather than assumes.
	 */
	struct ExteriorConfig
	{
		/// Whether the file named one.
		bool given = false;

		/// `rho_Gamma`, the radius of the semicircle, in metres. > 0, and it
		/// must fit strictly inside the `[mesh]` box.
		double radius = 0.0;

		/// The axial position of its centre, in metres. Zero is the common case.
		double centreZ = 0.0;

		/// How many Gegenbauer modes, degrees 2 .. modes + 1. The truncation is
		/// the ONLY approximation in the exterior -- the map is exact mode by
		/// mode -- and its error converges spectrally in the smoothness of the
		/// trace, so this is small in practice. >= 1.
		int modes = 0;
	};

	struct BoundaryConfig
	{
		BoundaryDataType type = BoundaryDataType::Zero;
		ShapeConfig shape;
		LimiterConfig limiter;
		XPointConfig xpoint;
		ExteriorConfig exterior;
	};

	// [solver] AssemblyMode -- who computes the element-local work.
	//
	// A PERFORMANCE choice and never a numerical one: MFEM guarantees the two
	// modes agree BIT FOR BIT, and MEQ asserts it on both a linear and a
	// nonlinear source. So this key can be changed between two runs of the same
	// configuration and the answers must not move.
	//
	// It is spelled as MEQ's own enum rather than as
	// GradShafranovSolver::AssemblyMode because this header is deliberately
	// MFEM-free -- Config, Profiles, Source and SourceFactory are what CI can
	// build without the MFEM branch it cannot obtain. The driver maps it and
	// asks GradShafranovSolver::assemblyModeAvailable() whether the build can
	// honour it, which is a question only the linked library can answer.
	// [solver] LocalFactorMode. See the member in SolverConfig.
	enum class LocalFactorModeType
	{
		Serial,
		Batched
	};

	// [solver] TraceAssemblyMode. See the member in SolverConfig.
	enum class TraceAssemblyModeType
	{
		Serial,
		Batched
	};

	enum class AssemblyModeType
	{
		// One thread.
		Serial,
		// Thread every element-local loop in the hybridization -- assembly, and
		// since MFEM threaded MultNL(), the residual and the Jacobian of every
		// NPC step too. **THE DEFAULT since 2026-09-04**, matching the library's
		// own: worth 2.4x on HighBetaConvergence, which is the case that used to
		// argue against it, and a wash at one thread.
		//
		// Needs an MFEM built with MFEM_USE_OPENMP and MFEM_THREAD_SAFE. A build
		// without them defaults to Serial, and a file that ASKS for "threaded"
		// there is refused by the driver with a message about the build rather
		// than relaying an exception.
		Threaded,
		// The interior-face potential term through one batched kernel instead
		// of a host call per face. **A DEVICE mode**: MFEM's own note is that
		// its D accumulation goes through AtomicAdd, which costs on a host
		// where the per-face loop's plain += does not, so on the CPU path this
		// is a prerequisite for the offload work rather than a speedup. It
		// falls back silently, and the run reports whether it was taken.
		Batched
	};

	// [solver] TraceSolver -- which direct solver factorises the trace system.
	//
	// Also a performance choice, and a LICENCE choice. All three reach the same
	// equilibrium; `theTraceSolversAgree` pins them to 1e-10 and they measure
	// 1e-14 or better.
	enum class TraceSolverType
	{
		// SuiteSparse UMFPACK, METIS ordering. The default, and the only one
		// present in every build -- every rate in the suite was measured with it.
		UMFPack,
		// oneMKL PARDISO. Faster than UMFPack even single-threaded, and it takes
		// MKL threads where UMFPack cannot. Needs MFEM_USE_MKL_PARDISO.
		Pardiso,
		// NVIDIA cuDSS. Correct, and not recommended on the strength of any
		// timing taken on this machine. Needs MFEM_USE_CUDSS and an mfem::Device.
		cuDSS // NOLINT(readability-identifier-naming)
	};

	// [solver] -- Newton on the outside, a linear solve on the inside.
	/// `[solver] LineSearchMerit`. Mirrors
	/// meq::GradShafranovSolver::LineSearchMerit, which this layer cannot name:
	/// Config is deliberately MFEM-free.
	enum class LineSearchMeritChoice
	{
		Augmented,
		Field
	};

	struct SolverConfig
	{
		// Who computes the element-local work, and which direct solver
		// factorises the trace system. Neither changes the answer.
		//
		// THIS DEFAULT IS A PLAIN Threaded RATHER THAN A BUILD-CONDITIONAL ONE,
		// deliberately, and it is the one place the two defaults differ. This
		// header is MFEM-free -- it is one of the four translation units CI
		// compiles without the library -- so it cannot ask whether this build
		// has OpenMP. apps/meq.cpp asks assemblyModeAvailable() and falls back
		// to Serial when the answer is no, which is what keeps a file that says
		// nothing working on every build.
		AssemblyModeType assemblyMode = AssemblyModeType::Threaded;

	// [solver] LocalFactorMode -- how the element-local blocks are factored.
	//
	// A FOURTH AXIS, independent of AssemblyMode, and the one of the batched
	// modes with a measured HOST win: upstream's in-situ figure is 10-12%
	// faster at order 2 and 24% slower at order 6. MEQ's machine cases run at
	// k = 2. Bit exact either way.
	//
	// Serial is the default, which is MFEM's, so a file that says nothing gets
	// what it always got.
	LocalFactorModeType localFactorMode = LocalFactorModeType::Serial;

	// [solver] TraceAssemblyMode -- how the element blocks reach the global
	// trace matrix.
	//
	// **THE ONE KEY HERE THAT IS NOT BIT EXACT**, and it is opt-in for that
	// reason. The two modes agree on the pattern and on every value to the bit,
	// but a row's columns come out in a different ORDER, so SparseMatrix::Mult
	// reassociates and the trace solve differs in its last bits -- and Serial's
	// order cannot be reproduced, being a function of the VALUES. What it buys
	// is the sparse third of ComputeH(), which on a host rebuilds a linked-list
	// matrix every linearisation.
	TraceAssemblyModeType traceAssemblyMode = TraceAssemblyModeType::Serial;
		// PARDISO, matching the library's own defaultTraceSolver(). It is faster
		// than UMFPack single-threaded and it is the only one of the two whose
		// MKL threads are spendable under the threaded assembly above, the two
		// responding to MKL_NUM_THREADS in OPPOSITE directions.
		TraceSolverType traceSolver = TraceSolverType::Pardiso;

		// Whether the file SAID "threaded" or "pardiso" or merely inherited
		// them, and the two must behave differently on a build that lacks the
		// backing package. A file that asks is refused, because a caller naming
		// a mode or a solver has a reason; a file that says nothing falls back
		// -- to Serial, and to UMFPack -- and runs. Without these flags the
		// defaults would make MEQ unusable on any build without OpenMP or
		// without oneMKL, which is most of them and is exactly what CI builds.
		bool assemblyModeWasGiven = false;
		bool traceSolverWasGiven = false;

		// Newton stops when either ||R|| <= NewtonAbsoluteTolerance or
		// ||R|| <= NewtonRelativeTolerance * ||R_0||, and fails after
		// NewtonMaxIterations. Both tolerances are in the units of the residual
		// -- dimensionless, once scaled by the initial residual.
		int newtonMaxIterations = 20;
		double newtonRelativeTolerance = 1.0e-8;
		double newtonAbsoluteTolerance = 1.0e-12;

		/*
		 * PlasmaSupportSweeps -- THE SUPPORT'S OWN OUTER LOOP, and the one
		 * discrete state the bordered Newton structurally cannot carry.
		 *
		 * With `ConfineToPlasma = true` the set of elements carrying `F` is a
		 * functional of the iterate: the pointwise `Psi > 0` test moves with
		 * `psi`, and so does the connected component the flood fill reaches.
		 * Its derivative is a surface term on a moving edge and the Jacobian
		 * does not carry it -- which is not a small error. **MEASURED, one key
		 * changed and nothing else: on the diverted machine the bootstrap
		 * stalls at the 200-iteration cap with the support moving and converges
		 * in 21 steps with the confinement off.** MEASUREMENTS.md M-82.
		 *
		 * The answer is FREE-BOUNDARY-PLAN.md section 10.5's, applied to the
		 * support rather than to the bounding point: fix the topology within a
		 * solve and re-decide it between solves. This is how many times the
		 * driver may re-decide it -- freeze at the answer just reached, solve
		 * again, and stop when the support stops moving.
		 *
		 * **ZERO IS THE DEFAULT AND MEANS TODAY'S BEHAVIOUR**, the support
		 * moving inside Newton, so every existing file is bit-unchanged. It is
		 * opt-in rather than implied by `ConfineToPlasma` because the LIMITED
		 * machine converges perfectly well without it -- examples/
		 * limited-tokamak.toml takes 11 Newton steps -- and a key that silently
		 * changed which equilibrium those runs report would be exactly the kind
		 * of thing CLAUDE.md refuses to expose.
		 *
		 * One sweep is a freeze and a solve, so `1` is "decide the support from
		 * the initial guess and hold it", and the loop below that is the whole
		 * of the outer iteration XP-3 leaves standing.
		 */
		int plasmaSupportSweeps = 0;

		// [solver] XPointMeritWeight -- a multiplier on the length that puts
		// XP-3's two rows into the LINE SEARCH's merit, and into nothing else.
		// The border solves q_r = q_z = 0 whatever this is, so it changes how
		// many iterations a solve costs and must not change the answer.
		// One is the natural scale ( r h, which turns q into a flux across the
		// X-point's own element ) and is the default. See
		// meq::GradShafranovSolver::setXPointMeritWeight for the plateau it
		// exists to attack.
		//
		// AND THERE IS NO GOOD UNIVERSAL VALUE: the sensitivity INVERTS between
		// cases. MAST goes 56 iterations to 35 at weight 20 and FAILS at 30,
		// while examples/diverted-tokamak.toml goes 14 to 82 at weight 4. The
		// natural scale is right on one and wrong on the other, so this is a
		// per-case knob and 1.0 stays the default. MEASUREMENTS.md M-113.
		double xPointMeritWeight = 1.0;

		/**
		 * `[solver] BorderRegularisation` and
		 * `[solver] BorderCollinearityRegularisation` -- Levenberg damping on
		 * the DENSE border solve, so a near-singular Schur complement gives a
		 * damped step rather than throwing. Both default to zero, which is off
		 * and bit-identical.
		 *
		 * EXPOSED ON THE SAME GROUNDS AS `AssemblyMode`: they change the
		 * Jacobian and so the WORK, and they cannot change the fixed point,
		 * because the fixed points of `x - alpha M^-1 G( x )` are the zeros of
		 * `G` for any non-singular `M`. That is the line `Globalisation` falls
		 * the wrong side of and these do not.
		 *
		 * MEASUREMENTS.md M-119 is why they exist: `machine-c-mast-shaped`
		 * throws *"the bordered Jacobian is singular in ( psi_ax, psi_bnd, a )"*
		 * once the field block is made non-singular, so what is left degenerate
		 * is the border itself.
		 */
		double borderRegularisation = 0.0;
		double borderCollinearityRegularisation = 0.0;

		/**
		 * `[solver] TopologyRetry` -- how many 0.75 reductions a trial that
		 * BROKE THE TOPOLOGY gets before the halving ladder takes over. Zero is
		 * off and bit-identical.
		 *
		 * A constraint that cannot be evaluated at all and a merit that came
		 * back worse are different failures wanting different responses, and
		 * MEQ's line search reported them identically -- which is why five of
		 * M-119's six cold failures print one message that cannot say which
		 * they are.
		 */
		int topologyRetry = 0;

		/**
		 * `[solver] UpDownSymmetry` -- find the up-down symmetric solution by
		 * projecting every iterate onto the symmetric subspace. Off by default.
		 *
		 * **A STATEMENT ABOUT WHICH EQUILIBRIUM IS WANTED, exposed for the
		 * reason `[boundary.xpoint]`'s seed is** -- a free boundary has more
		 * than one solution and something must say which. What it replaces is
		 * worse: seeding ONE null of a double null asks the user to break a
		 * symmetry the machine does not break, and MEQ then follows whichever
		 * saddle its search wanders to. On freegsnke's MAST-U, symmetric to
		 * 2.4e-10 in its currents, that lands 1.574 m away on the other side.
		 *
		 * It is also the vertical instability's cure -- see
		 * GradShafranovSolver::setUpDownSymmetry() -- and it REFUSES an
		 * asymmetric mesh rather than deleting the asymmetry.
		 */
		bool upDownSymmetry = false;

		// [solver] LineSearchMerit -- what the Armijo backtracking compares,
		// and NOT what is solved or when it stops. See
		// meq::GradShafranovSolver::LineSearchMerit.
		LineSearchMeritChoice lineSearchMerit = LineSearchMeritChoice::Augmented;


		/**
		 * `PicardSweeps` -- FIND THE BASIN BEFORE DRIVING THE NEWTON.
		 *
		 * A free-boundary Grad-Shafranov problem has SEVERAL solutions and
		 * which one is reported is decided by where the iteration starts.
		 * Measured on freegs4e's TestTokamak from a cold start -- the
		 * conductors at their given currents, a current blob shaped by the
		 * design profiles, and nothing from a converged equilibrium -- MEQ's
		 * bordered Newton converges cleanly to `psi_ax = 1.30e-01` where the
		 * reference has `8.27e-02`, at the same `I_p`. Seeding `PsiAxis` with
		 * the reference's own converged value changes it by not one digit, so
		 * it is the basin and not the normalisation.
		 *
		 * **AND THERE IS NO GLOBALISATION TO REACH FOR.**
		 * `GradShafranovSolver::solve` refuses every
		 * `Globalisation` but `None` once `psi_ax` is an unknown: the KINSOL
		 * paths drive a residual of their own and the Picard ones build no
		 * Jacobian to border. The bordered loop's own Armijo backtracking is
		 * what there is, and on these machines the step it is damping is
		 * already leaving the branch.
		 *
		 * **SO THE PICARD IS OUTSIDE THE BORDER RATHER THAN INSIDE IT.** With
		 * `( psi_ax, psi_bnd )` held FIXED, a normalised source is an ordinary
		 * `meq::Source` -- `F( r, z, psi )` with no unknowns in it -- and the
		 * field solve is the unbordered problem MEQ has always been able to
		 * solve. This key counts sweeps of
		 *
		 *     freeze the normalisation at the current estimate
		 *     solve the unbordered problem
		 *     re-read `psi_ax` at the located O-point and `psi_bnd` at the
		 *         bounding point, and rescale the profiles to the target `I_p`
		 *
		 * which is freegs4e's own algorithm, and which converges globally where
		 * a Newton converges locally. The state it reaches is handed to the
		 * bordered solve as its initial guess.
		 *
		 * **IT IS AN INITIALISER AND NOT A SOLVER.** It is not asked to meet
		 * any tolerance and its answer is not the run's; what it has to do is
		 * choose the branch, which is a topological question and survives a
		 * loose sweep. `Globalisation::PicardThenNewton`'s own documentation
		 * makes the same distinction for the unbordered problem.
		 *
		 * Zero is the default, so every existing configuration is unchanged,
		 * and a run that converges without it does not want it.
		 */
		int picardSweeps = 0;

		/**
		 * `PicardBlend` -- how much of each Picard sweep's answer is taken.
		 *
		 * `1.0` is the plain fixed point and `0.5` is what freegs4e's own
		 * adaptive blending averages out at on these cases. Under-relaxation
		 * is what stops the normalisation and the support chasing each other:
		 * the core is `{ psi > psi_bnd }`, so a `psi_bnd` that overshoots
		 * shrinks the plasma, which concentrates the current, which moves
		 * `psi_bnd` further -- measured here on a standalone reimplementation,
		 * that runaway takes the core from 310 cells to 1 in 150 sweeps.
		 */
		double picardBlend = 0.5;

		// THERE ARE NO INNER-LINEAR-SOLVE CONTROLS HERE, AND THAT IS THE POINT.
		// LinearMaxIterations and LinearTolerance used to sit in this struct,
		// parsed and validated and read by nothing. MEQ's trace solve is
		// DIRECT, so an iteration count and a tolerance have nothing to
		// control; the keys are now refused at parse time with a message that
		// says so. See refuseIterativeSolverKeys in Config.cpp.
	};

	// [initialguess] -- where Newton starts. See docs/running.rst, which names
	// the three restart routes and what each one carries.
	enum class InitialGuessType
	{
		// The Dirichlet datum extended inward, which is what prepare() does
		// anyway. The default, and a cold start.
		None,
		// A ramp putting psi = 0 in the interior rather than on the boundary.
		// Not a nicety: EVERY GS-2 section 4.2-4.5 source vanishes at psi = 0,
		// so with homogeneous data psi = 0 SOLVES the problem and Newton stops
		// on it in zero iterations. See CLAUDE.md under Traps.
		Ramp,
		// THE CONDUCTORS' OWN FIELD, computed from this file's [[coils]] and
		// nothing else. `Type = "conductors"`.
		//
		// **IT EXISTS BECAUSE BUILDING A GOOD GUESS IS MEQ'S WORK.** Every
		// machine example in this tree hands the solver a .gf that
		// tools/freegs4e-benchmark/mkexactguess.py reconstructed by summing
		// Green's functions over the source -- which is MEQ asking the user to
		// do its convergence work, and which no user of a free-boundary code
		// should have to do. freegsnke solves its MAST-U from its own default
		// initialisation with nothing supplied.
		//
		// MEQ ALREADY HAS EVERYTHING THAT GUESS IS MADE OF. The conductors and
		// their currents are in this very file; meq::coilPsi() integrates one
		// analytically by Carlson's elliptic integrals, agreeing with the
		// filament limit to 7.0e-13; and the sum is Delta*-HARMONIC off the
		// conductors at measured rate 2.00, which
		// `the_conductor_field_is_delta_star_harmonic_off_the_conductors`
		// asserts. So the vacuum half of the guess needs no external file and
		// no helper script.
		//
		// THE PLASMA IS OPTIONAL AND IS AN ELLIPTICAL COLUMN. With no
		// `CentreR` this is the vacuum field alone, which has the right scale
		// and topology outside the plasma and says nothing about the core; on
		// MAST-U that converges in 9 iterations to the WRONG branch, which is
		// what a guess describing no plasma buys. `CentreR`/`CentreZ` then
		// place a guessed magnetic axis and `[source] PlasmaCurrent` is spread
		// over an ellipse of semi-axes `RadiusR`/`RadiusZ` about it, defaulting
		// to half the major radius and circular.
		//
		// **NOT A FILAMENT, AND THAT IS MEASURED RATHER THAN PREFERRED.**
		// meq::filamentPsi() diverges logarithmically at the filament, so I_p
		// carried on one at the guessed axis reads 6.667e-01 at 3 mm from it
		// against MAST-U's reference psi_axis of 9.187e-02 -- an unbounded
		// spike at exactly the point the axis search exists to find, and the
		// solve fails. **And an ellipse rather than a rectangle**, which
		// meq::Coil already is and which would be bounded: a uniform current
		// density over a rectangle carries a logarithm in its second
		// derivatives at each of four corners, and those are artefacts of the
		// shape rather than anything the equilibrium puts there.
		Conductors,
		// A paraboloid BUMP: a core, positive inside an ellipse about a
		// prescribed centre and zero outside it.
		//
		// THE RAMP IS THE WRONG SHAPE FOR A FREE BOUNDARY AND THIS IS WHY THERE
		// ARE TWO. A ramp is antisymmetric in z -- it exists to put psi = 0 in
		// the interior so the trivial branch is not a fixed point -- and it
		// describes no plasma at all. A free-boundary solve has to be told
		// WHICH EQUILIBRIUM to find: the fixed-boundary rehearsal against
		// freegs4e measured three converged solutions of one discrete problem,
		// with the physical one lying BETWEEN the two a ramp sweep reaches. So
		// WHICH BRANCH is genuinely the user's to state, and a core is what
		// selects the one that is physical.
		//
		// **THAT IS NOT A LICENCE TO REQUIRE A FIELD FILE.** Stating which
		// equilibrium is wanted -- an X-point near here, this much current --
		// is a physical input and MEQ already takes it as [boundary.xpoint],
		// [source] PsiAxis and PlasmaCurrent. Handing MEQ a .gf reconstructed
		// by a helper script is a different thing: it is the solver asking the
		// user to do its convergence work. Every machine example in this tree
		// does it, and that is a DEFECT those files now name as one rather
		// than a precondition they document.
		Bump,
		// An MFEM GridFunction and its mesh, from a previous MEQ run.
		GridFunction
	};

	// [initialguess]
	struct InitialGuessConfig
	{
		InitialGuessType type = InitialGuessType::None;

		// GridFunction: the stored potential and the mesh it lives on. The mesh
		// must match the one being solved on -- this is the EXACT restart of
		// docs/running.rst, not the interpolating one.
		//
		// THE INTERPOLATING ONE IS NOT "NOT WRITTEN", WHICH IS WHAT THIS SAID,
		// AND IT IS NO LONGER UNWIRED EITHER. meq::FieldTransfer in
		// WarmStart.hpp is built on mfem::FindPointsGSLIB, the install has
		// MFEM_USE_GSLIB = YES, and WarmStartConvergence is a registered ctest.
		// Since 2026-09-02 apps/meq.cpp uses it: a mesh-count mismatch
		// INTERPOLATES rather than being refused, so restarting from a run at
		// another resolution -- the ordinary way to use a stored answer -- works.
		// The exact restart is still taken when the meshes match, since it is
		// every coefficient rather than an interpolation.
		std::string file;
		std::string meshFile;

		// Ramp: psi runs from -Amplitude to +Amplitude across z, so that the
		// interior crosses zero and the trivial branch is not a fixed point of
		// the iteration.
		//
		// Bump: the PEAK value, at the centre. Positive.
		double amplitude = 0.3;

		// Bump and Conductors. The centre, in metres, and the extent --
		// psi = Amplitude*( 1 - ( dr/RadiusR )^2 - ( dz/RadiusZ )^2 ) where that
		// is positive and zero elsewhere. RadiusZ defaults to RadiusR, which is
		// the circular case.
		//
		// On Conductors the same four numbers place the ELLIPTICAL PLASMA
		// COLUMN instead: CentreR/CentreZ are the guessed magnetic axis,
		// RadiusR/RadiusZ its semi-axes, and RadiusR defaults to half of
		// CentreR -- a default that cannot reach the axis whatever CentreR is.
		// A zero centreR there means no column at all.
		double centreR = 0.0;
		double centreZ = 0.0;
		double radiusR = 0.0;
		double radiusZ = 0.0;
	};

	// [adaptivity] -- the stage-6 loop, exposed rather than rebuilt.
	enum class MarkingStrategy
	{
		Doerfler,
		Maximum
	};

	// [adaptivity]
	struct AdaptivityConfig
	{
		bool enabled = false;
		int maxIterations = 10;
		MarkingStrategy strategy = MarkingStrategy::Doerfler;
		// The Doerfler fraction: refine the smallest set of elements carrying
		// this share of the total estimated error. Ignored for Maximum.
		double theta = 0.6;
		// Stop once the estimate falls below this. Absolute, in the estimator's
		// own norm.
		double targetError = 1.0e-6;
	};

	// [output] -- one directory and one prefix; the file names follow from
	// them. (The scheme this replaces named all five output files separately,
	// and drifted out of step with itself.)
	struct OutputConfig
	{
		// Directory: where the files are written. Created by the caller if it
		// does not exist.
		std::string directory = ".";
		// Prefix: the stem of every output file name.
		std::string prefix = "meq";

		// GridNR, GridNZ: the ( R, Z ) sampling grid for the NetCDF file, in
		// NODES. Nothing to do with [mesh] NR/NZ, which are pre-refinement
		// CELLS of the solve -- deriving the output grid from those would tie
		// the resolution of the picture to the coarsest description of the
		// mesh, and give a 4x5 grid for a 1536-element solve.
		int gridNR = 129;
		int gridNZ = 129;

		// <Directory>/<Prefix>.mesh, the mesh actually solved on.
		std::string getMeshFile() const;
		// <Directory>/<Prefix>_psi.gf, the flux function psi [Wb per radian].
		std::string getPsiFile() const;
		// <Directory>/<Prefix>_grad_psi.gf, the HDG flux variable q.
		std::string getGradPsiFile() const;
		// <Directory>/<Prefix>_psistar.gf, the POST-PROCESSED potential psi*.
		//
		// A SEPARATE FILE FROM _psi.gf AND NOT A REPLACEMENT FOR IT. psi* lives
		// in P_(k+1) and converges at k+2, so it is the better field to look at
		// and the better field to sample -- and it does not fit the degree-k
		// potential space an exact restart reads back into. _psi.gf therefore
		// keeps psi_h and this carries psi*; see apps/meq.cpp's write block.
		std::string getPsiStarFile() const;

		// ------------------------------------------------------------------
		// [output] FluxSurfaces and its four sizes: the ( Psi, theta ) file,
		// INVERSION-PLAN.md stage IN-6.
		//
		// OFF BY DEFAULT, AND NOT BECAUSE IT IS EXPERIMENTAL. It costs a
		// contour trace and an angle fit per surface -- comparable with the
		// solve on a coarse mesh -- and it is the one output that can FAIL on a
		// run that solved perfectly well: a level whose surface is not closed,
		// not star-shaped about the axis, or off the mesh has no flux-surface
		// average, and meq refuses rather than inventing one. A run asking for
		// the answer should not be made to pay for a reduction it did not ask
		// for, and should not be made to fail for one either.
		//
		// FluxSurfaces: write <Directory>/<Prefix>_surfaces.nc.
		bool fluxSurfaces = false;

		// FluxSurfaceCount, FluxAngleCount: surfaces in the family, and nodes
		// on each. Both are resolutions of the OUTPUT and have nothing to do
		// with the mesh, in the same way GridNR and GridNZ do not.
		int fluxSurfaceCount = 24;
		int fluxAngleCount = 128;

		// FluxInnerCut, FluxOuterCut: the range of normalised flux the family
		// covers, 0 on the magnetic axis and 1 on the plasma boundary.
		//
		// BOTH ENDS ARE CUT AND FOR DIFFERENT REASONS -- see
		// src/meq/FluxFamily.hpp, which carries the decision, and CLAUDE.md's
		// IN-6 section, which carries the measurement it was made from. In
		// short: at the inner end a surface shrinks to a point and drho/dpsi is
		// unbounded; at the outer end the surface stops being made of solved
		// data, because on a curved boundary Omega_h is inscribed in Gamma and
		// the outermost surfaces cross the band. Nothing FAILS at either end,
		// which is exactly why the cut has to be a decision.
		double fluxInnerCut = 0.05;
		double fluxOuterCut = 0.95;

		// <Directory>/<Prefix>_surfaces.nc, the ( Psi, theta ) grid.
		std::string getFluxSurfaceFile() const;
	};

	// A parsed, validated configuration. Construction either succeeds and
	// leaves every accessor meaningful, or throws ConfigError.
	class Configuration
	{
		public:
			// Parse a TOML file. Throws ConfigError if it cannot be read, does
			// not parse, or does not describe a runnable problem.
			explicit Configuration( std::string const & fileName );

			// Parse an already-loaded document. `source` is used only in error
			// messages, to say where the document came from.
			Configuration( toml::value const & document, std::string const & source );

			// Parse TOML held in memory -- for tests, and for embedding.
			static Configuration fromString( std::string const & text, std::string const & source = "<string>" );

			MeshConfig const & getMesh() const noexcept { return meshOptions; };

			/// The `[[coils]]` blocks, in file order. Empty unless the file has
			/// any, which every fixed-boundary configuration does not.
			CoilConfig const & getCoils() const noexcept { return coilOptions; };
			DiscretisationConfig const & getDiscretisation() const noexcept { return discretisationOptions; };
			SourceConfig const & getSource() const noexcept { return sourceOptions; };
			BoundaryConfig const & getBoundary() const noexcept { return boundaryOptions; };
			SolverConfig const & getSolver() const noexcept { return solverOptions; };
			OutputConfig const & getOutput() const noexcept { return outputOptions; };
			InitialGuessConfig const & getInitialGuess() const noexcept
			{ return initialGuessOptions; };
			AdaptivityConfig const & getAdaptivity() const noexcept
			{ return adaptivityOptions; };

			// Where this configuration came from: the file name, or the
			// `source` label given for an in-memory document.
			std::string const & getFileName() const noexcept { return sourceName; };

		private:
			void parse( toml::value const & document );

			std::string sourceName;

			MeshConfig meshOptions;
			CoilConfig coilOptions;
			DiscretisationConfig discretisationOptions;
			SourceConfig sourceOptions;
			BoundaryConfig boundaryOptions;
			SolverConfig solverOptions;
			OutputConfig outputOptions;
			InitialGuessConfig initialGuessOptions;
			AdaptivityConfig adaptivityOptions;
	};

}

#endif // MEQ_CONFIG_HPP
