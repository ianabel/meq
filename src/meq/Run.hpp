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
