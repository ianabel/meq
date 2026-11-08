#ifndef MEQ_RUN_HPP
#define MEQ_RUN_HPP

/*
 * ONE RUN OF MEQ: THE CONFIGURATION READ, AND THE RUN SET UP, SOLVED AND
 * WRITTEN.
 *
 * A LIBRARY CLASS SO THAT EVERY FRONT END RUNS THE SAME RUN. `meq` is one
 * (apps/meq.cpp, which keeps argv and the exit code and nothing else), and the
 * Python module PYTHON-PLAN.md builds is another. A front end assembled out of
 * the library's pieces -- Configuration, makeSource(), GradShafranovSolver --
 * would be a second driver: it would make none of the refusals this one makes
 * on the user's behalf, and two drivers can each converge, to different
 * equilibria, which is what M-264 asserts this one and the library do not.
 *
 * Its report goes to stdout and its diagnostics to stderr, as `meq`'s always
 * have.
 */

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Dataset.hpp"
#include "EquilibriumSnapshot.hpp"
#include "Profiles.hpp"

namespace meq
{
	/// The exit codes docs/running.rst gives the user. apps/meq.cpp's header
	/// says what each means.
	enum ExitCode
	{
		Solved = 0,
		ConfigurationError = 1,
		SolveFailed = 2,
		OutputFailed = 3,
		// An inverse solve converged and was written, to the BEST FIT its
		// free currents allow: the targets were not met. Not a plain success,
		// by INVERSE-AND-RECONSTRUCTION-PLAN.md section 4.5 -- a reader
		// excuses a shape residual they would refuse if it were named.
		TargetsNotMet = 4
	};

	/**
	 * WHAT A RUN TAKES BESIDE ITS CONFIGURATION FILE: the choices that are
	 * about the PROCESS rather than the equilibrium. `meq`'s command line
	 * sets them, and apps/meq.cpp's main() says why each is not a key.
	 */
	struct RunOptions
	{
		/// The TOML file -- or, with `configText`, only the name messages and
		/// the output's `config_file` give the configuration.
		std::string configPath;

		/// The configuration AS TEXT, in place of reading `configPath`: a
		/// front end that holds the document and has no file, such as the
		/// Python module's from_dict(). The text is what the .nc records as
		/// `input_toml`, as it records a file's.
		std::string configText;

		/// With `configText`, the directory its relative input paths resolve
		/// against, as a file's own directory does. Empty leaves them as
		/// given, relative to the working directory.
		std::string baseDirectory;

		/// Profiles as tables in place of the files the configuration names,
		/// by slot: meq::makeSource()'s `tables`. The .nc lists them as
		/// `profile_tables`.
		ProfileTables profiles;

		/// A .nc a previous run of this configuration wrote, to start from in
		/// place of the configuration's [initialguess]: on the same mesh that
		/// is the exact restart, the whole iterate ( q, psi, psihat ), labels
		/// and exterior coefficients, which is how a front end re-solves warm
		/// (PYTHON-PLAN.md, PY-3). A configuration with [adaptivity] is
		/// refused beside it: its mesh would be rebuilt from the file and
		/// re-adapted, not the one the state was solved on.
		std::string restartFrom;

		/// R_0( Psi ), the midplane crossing an `RMidplane` table is mapped
		/// through, HELD for the whole solve: a front end hands in the last
		/// answer's (Run::midplaneCrossing()). Empty refuses such a table.
		std::vector<Knot> heldCrossing;


		/// MEQ_VERSION of the front end's build, which the files written carry.
		/// The front end's because the commit is stamped into the binary at
		/// link time (cmake/GitCommit.cmake), and meq_core is not re-built
		/// when the commit moves.
		std::string version;

		/// `--device`: an mfem::Device backend string, and whether one was
		/// given.
		std::string deviceName;
		bool deviceGiven = false;

		/// `--mesh-ready`: the caller says the [mesh.generate] mesh is built.
		bool meshReady = false;

		/// `--profile` and `--trace-solver-stats`.
		bool wantProfile = false;
		bool wantTraceSolverStats = false;

		/// `--trace-refinement-cap`; unset leaves the solver's own default.
		std::optional<int> traceRefinementCap;

		/// `--q-iteration-cap`; unset leaves meq::OuterNewtonOptions' own.
		std::optional<int> qIterationCap;

		/// `--gf` and `--vtk`: formats beside the .nc, either asking enough.
		bool wantGridFunctions = false;
		bool wantVtk = false;
	};

	class Run
	{
		public:
			explicit Run( RunOptions const &options );
			~Run();

			Run( Run const & ) = delete;
			Run &operator=( Run const & ) = delete;

			/// Configure, set up, solve and write, in that order, and the exit
			/// code: meq::ExitCode's. Once per object, in place of the three
			/// below.
			int execute();

			/**
			 * THE SAME RUN A PHASE AT A TIME, for a front end that wants to
			 * act between them. Each returns meq::ExitCode's code: Solved
			 * when the phase succeeded, the run's own failure code when it did
			 * not, and write() TargetsNotMet for an inverse best fit, exactly
			 * as execute() would have.
			 *
			 * IN ORDER AND EACH ONCE, AND OUT OF ORDER IS REFUSED: a phase
			 * called before the one it follows has succeeded, again, or after
			 * any phase has failed throws std::logic_error. Calling them in
			 * order is execute(), to the bit.
			 */
			int prepare();          ///< read the configuration and set the run up
			int solve();            ///< solve, and refine if the file asks
			int write();            ///< write what [output] and the options ask for

			/**
			 * THE .nc OF THE SOLVED STATE AT `stem`, AND NOTHING ELSE: no
			 * [output] directory or prefix, no grid functions, no VTK, no flux
			 * surfaces. What a re-solve restarts from (RunOptions::
			 * restartFrom). Any time after a successful solve, and it does not
			 * count as write(): the phase order is untouched.
			 * @throws std::logic_error before a successful solve.
			 */
			int writeRestart( std::string const &stem );

			/**
			 * THE SOLVED STATE'S MIDPLANE CROSSING, R_0( Psi ) on Psi in
			 * [ 0, 1 ], as knots: where each surface crosses z = 0, read off
			 * the answer's own field along the midplane and its labels, as the
			 * .nc's true-midplane densities are (meq::
			 * tabulateMidplaneCrossing()). What an `RMidplane` table of the
			 * next solve is held at.
			 * @throws std::logic_error before a successful solve, or for a
			 *         support that is not a mirror's box.
			 * @throws std::runtime_error as tabulateMidplaneCrossing() does,
			 *         for a midplane flux that is not monotone.
			 */
			std::vector<Knot> midplaneCrossing();

			/**
			 * THE .nc THE LAST write() OR writeRestart() WROTE, IN MEMORY: its
			 * root group, every dimension, variable and attribute, recorded by
			 * the calls that wrote the file (meq::NetCDFWriter::dataset()).
			 * Shared, so it outlives the run that wrote it.
			 * @throws std::logic_error before anything was written.
			 */
			std::shared_ptr<Dataset const> dataset() const;

			/**
			 * psi, B_R AND B_Z AT POINTS OF ONE'S CHOOSING, AS THE .nc's GRID
			 * HAS THEM: the same location, band continuation and sampling the
			 * write phase applies to its lattice (meq::GridSampler's point
			 * set), so at a lattice node's coordinates the values are that
			 * node's, to the bit. `inside` and `extrapolated` are the .nc's
			 * masks; outside, the values are NaN.
			 * @throws std::logic_error before a successful solve, or for two
			 *         lengths or none.
			 */
			struct Samples
			{
				std::vector<double> psi, bR, bZ;
				std::vector<signed char> inside, extrapolated;
			};
			Samples sampleAt( std::vector<double> const &radius,
			                  std::vector<double> const &z );

			/**
			 * THE _surfaces.nc OF THE SOLVED STATE, IN MEMORY: the family the
			 * write phase extracts at [output]'s FluxSurfaceCount,
			 * FluxAngleCount and cuts, as its file holds it -- write()'s, where
			 * [output] FluxSurfaces wrote one, or written now to @a file.
			 * @throws std::logic_error before a successful solve, or beside a
			 *         box support, whose field lines are open.
			 * @throws std::runtime_error where the extraction fails, as the
			 *         write phase's warning does.
			 */
			std::shared_ptr<Dataset const> fluxSurfaces( std::string const &file );

			/**
			 * THE FAMILY'S 1-D QUANTITIES AT LABELS rho = sqrt( Psi_N ) OF
			 * ONE'S CHOOSING, under the _surfaces.nc names: V_prime, volume,
			 * cross_section_area, arc_length, surface_area, inverse_R_squared,
			 * grad_psi_squared_over_R_squared, abs_grad_psi, grad_psi_squared,
			 * and safety_factor where the file has it, by the family's monotone
			 * cubic (meq::FluxSurfaceFamily::at()); rho, normalised_flux and psi
			 * exactly. Of the family fluxSurfaces() holds, which is called
			 * first. @throws std::domain_error for a label outside the cut.
			 */
			std::vector<std::pair<std::string, std::vector<double>>>
			fluxSurfacesAt( std::vector<double> const &labels, std::string const &file );

			/**
			 * RE-SOLVE IN MEMORY WITH NEW TABLES: A SESSION
			 * (MANTA-EQUILIBRIUM-PLAN.md section 4.1).
			 *
			 * After a successful solve(), any number of times: @a tables
			 * replace the profile tables of the slots they name -- the rest
			 * stand as they were -- through the same factories and refusals
			 * prepare() applies, and the run solves again on its own mesh,
			 * spaces, forms and psi_c, warm from @a warm, or from the state
			 * it holds when @a warm is null. That is the .nc restart of
			 * PYTHON-PLAN.md's PY-3, field for field, held in memory: no file,
			 * no second setup.
			 *
			 * A FRESH SOLVER EACH TIME, ON PURPOSE. The answer is then a
			 * function of ( tables, warm ) and nothing else -- not of which
			 * state this run held when it was called -- which is the
			 * determinism a transport code's memoised Newton is pinned to.
			 * Its cost is the solver's construction; MEASUREMENTS.md
			 * M-SESSION-3 prices it against PY-3's whole setup.
			 *
			 * TRANSACTIONAL. On any failure the run is left exactly as it
			 * was -- its solver, sources, tables, snapshot and what write()
			 * would write -- and the code says why: SolveFailed, or
			 * ConfigurationError for tables prepare() would have refused.
			 * On success the run holds the new state, write() may be called
			 * again, and dataset() keeps what the last write wrote until then.
			 *
			 * REFUSED, as ConfigurationError with the run untouched, where a
			 * re-solve cannot restart on the state it holds: [adaptivity],
			 * whose mesh is the loop's; a mirror's box, whose R_midplane
			 * tables are iterated across solves; [inverse], whose coil set is
			 * shared state; and an mfem::Device, since a session is MaNTA's and
			 * MaNTA is double and host (the plan's D5).
			 *
			 * @throws std::logic_error before a successful solve.
			 * @throws std::invalid_argument for a snapshot of another run.
			 */
			int resolve( ProfileTables const &tables,
			             std::shared_ptr<EquilibriumSnapshot const> const &warm = nullptr );

			/**
			 * THE SOLVED STATE AS A VALUE: meq::EquilibriumSnapshot, with its
			 * flux-surface family extracted at [output]'s FluxSurfaceCount,
			 * FluxAngleCount and cuts (none on a box support). The same object
			 * until the next successful resolve(), and never changed by one.
			 * @throws std::logic_error before a successful solve.
			 * @throws std::runtime_error where the extraction fails.
			 */
			std::shared_ptr<EquilibriumSnapshot const> snapshot();

			/// Whether the configuration's support is a mirror's box, the one
			/// support midplaneCrossing() answers for. After prepare().
			/// @throws std::logic_error before it.
			bool supportIsBox() const;

		private:
			struct Impl;
			std::unique_ptr<Impl> impl;

			/// Which phase is next, and whether one has failed.
			enum class Stage
			{
				Fresh,
				Prepared,
				Solved,
				Written,
				Failed
			};
			Stage stage = Stage::Fresh;

			/// Refuse a phase whose predecessor is not `from`.
			void advance( Stage from, char const *what );
	};
}

#endif // MEQ_RUN_HPP
