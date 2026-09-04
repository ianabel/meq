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
 * and two nested forms beneath those:
 *
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
		// meq::NormalisedSource and CLAUDE.md, "Newton, and the obligation it
		// creates". False by default, because it changes what the tables mean.
		bool normalised = false;
		// PsiAxis: the starting value of psi_ax [Wb per radian]. REQUIRED when
		// Normalised is true and refused otherwise. It is a guess in the Newton
		// sense and not a scale factor: the iteration has to start inside a
		// basin, and at a fixed psi_ax the equation generally has a second,
		// non-physical solution the iteration can reach instead.
		double psiAxis = 0.0;
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
	// meq::RotatingSource and FLOW-PLAN.md.
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
		// Mu0, Normalised, PsiAxis: as MHDParameters, and meaning the same.
		double mu0 = 4.0e-7*3.14159265358979323846;
		bool normalised = false;
		double psiAxis = 0.0;
	};

	using SourceParameters = std::variant< SolovievParameters, MHDParameters,
	                                       ManufacturedParameters, RotatingParameters >;

	struct SourceConfig
	{
		SourceType type = SourceType::Soloviev;
		SourceParameters parameters = SolovievParameters{};

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

	struct BoundaryConfig
	{
		BoundaryDataType type = BoundaryDataType::Zero;
		ShapeConfig shape;
	};

	// [solver] -- Newton on the outside, a linear solve on the inside.
	struct SolverConfig
	{
		// Newton stops when either ||R|| <= NewtonAbsoluteTolerance or
		// ||R|| <= NewtonRelativeTolerance * ||R_0||, and fails after
		// NewtonMaxIterations. Both tolerances are in the units of the residual
		// -- dimensionless, once scaled by the initial residual.
		int newtonMaxIterations = 20;
		double newtonRelativeTolerance = 1.0e-8;
		double newtonAbsoluteTolerance = 1.0e-12;

		// THERE ARE NO INNER-LINEAR-SOLVE CONTROLS HERE, AND THAT IS THE POINT.
		// LinearMaxIterations and LinearTolerance used to sit in this struct,
		// parsed and validated and read by nothing. MEQ's trace solve is
		// DIRECT, so an iteration count and a tolerance have nothing to
		// control; the keys are now refused at parse time with a message that
		// says so. See refuseIterativeSolverKeys in Config.cpp.
	};

	// [initialguess] -- where Newton starts. See DRIVER-PLAN.md section 4.
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
		// An MFEM GridFunction and its mesh, from a previous MEQ run.
		GridFunction
	};

	// [initialguess]
	struct InitialGuessConfig
	{
		InitialGuessType type = InitialGuessType::None;

		// GridFunction: the stored potential and the mesh it lives on. The mesh
		// must match the one being solved on -- this is the EXACT restart of
		// DRIVER-PLAN.md section 4, not the interpolating one.
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
		double amplitude = 0.3;
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
