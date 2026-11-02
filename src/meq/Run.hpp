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
		/// The TOML file.
		std::string configPath;

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
			/// code: meq::ExitCode's. Once per object.
			int execute();

		private:
			struct Impl;
			std::unique_ptr<Impl> impl;
	};
}

#endif // MEQ_RUN_HPP
