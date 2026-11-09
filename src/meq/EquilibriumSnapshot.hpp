#ifndef MEQ_EQUILIBRIUM_SNAPSHOT_HPP
#define MEQ_EQUILIBRIUM_SNAPSHOT_HPP

/*
 * ONE CONVERGED EQUILIBRIUM OF A SESSION, BY VALUE -- MANTA-EQUILIBRIUM-PLAN.md
 * section 4.1.
 *
 * WHAT A TRANSPORT CODE HOLDS AS ITS OPAQUE `Equilibrium`. MaNTA keeps several
 * at once -- the last, its Jacobian build's, perturbed ones while it
 * differences -- and asks geometry and responses of one that is not the
 * latest, so a snapshot is a VALUE: nothing in it points into the solver that
 * produced it, and a later meq::Run::resolve() cannot change it. A handle onto
 * one live solver would answer for the newest state, silently.
 *
 * WHAT IT HOLDS IS WHAT PUTS A SOLVER BACK AT THE ANSWER: the whole iterate
 * ( q, psi, psihat ) and every border value the bordered Newton carries, so
 * that meq::Run::resolve() warm from it starts AT that equilibrium -- the
 * in-memory form of the .nc restart, field for field -- and so that
 * `manta/linearised-solve`'s linearisation solver (section 4.3) can load it
 * and assemble its Jacobian there. Beside them, the tables it was solved at,
 * since a linearisation needs its source, and the family its geometry is read
 * from. ~1-2 MB at 46,000 trace dofs; the family is ~55 KB.
 *
 * MFEM-FREE, like the family, so a front end can hold one without the library.
 */

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "FluxFamily.hpp"
#include "Profiles.hpp"

namespace meq
{
	struct EquilibriumSnapshot
	{
		/// The meq::Run it belongs to. Its state is numbered by that run's
		/// mesh and spaces, so it means nothing to another, and
		/// meq::Run::resolve() refuses it there.
		std::uint64_t session = 0;

		/// Which of that run's solves it is, from one; a linearisation is
		/// keyed on ( session, sequence ).
		std::uint64_t sequence = 0;

		/// The profile tables it was solved at, beside the configuration's
		/// files: meq::RunOptions::profiles and every resolve() since.
		ProfileTables profiles;

		/// GradShafranovSolver::state(): ( q, psi, psihat ), exactly.
		std::vector<double> state;

		/// Where state's blocks end: q is [ 0, fluxEnd ), psi
		/// [ fluxEnd, potentialEnd ) and the trace psihat
		/// [ potentialEnd, state.size() ).
		std::size_t fluxEnd = 0;
		std::size_t potentialEnd = 0;

		/// GradShafranovSolver::exteriorCoefficients(); empty on a fixed
		/// boundary.
		std::vector<double> exteriorCoefficients;

		/// The normalisation, psi_ax and psi_bnd, as the solve closed them.
		double psiAxis = 0.0;
		double psiBoundary = 0.0;

		/// Whether psi_bnd is an unknown of the solve -- a limiter or an
		/// X-point -- and so a seed a restart carries, as the .nc's
		/// `psi_boundary` is written only then.
		bool psiBoundaryIsUnknown = false;

		/// The located magnetic axis, where there was one.
		bool axisLocated = false;
		double axisR = 0.0;
		double axisZ = 0.0;

		/// The other border unknowns, where the run has them: the X-point,
		/// the I_p profile scale, the free coil currents, the q-driven
		/// loop's g^2 coefficients.
		bool xPointLocated = false;
		double xPointR = 0.0;
		double xPointZ = 0.0;
		bool hasProfileScale = false;
		double profileScale = 1.0;
		std::vector<double> inverseCurrents;
		std::vector<double> gSquaredCoefficients;

		/// The Newton steps the solve that produced it took: zero for a
		/// resolve that started at its answer.
		int newtonIterations = 0;

		/// The flux-surface family at the run's [output] FluxSurfaceCount,
		/// FluxAngleCount and cuts, which is where its geometry is read: null
		/// until meq::Run::snapshot() has extracted it, and on a box support,
		/// which has no closed surfaces.
		std::shared_ptr<FluxSurfaceFamily const> family;
	};
}

#endif // MEQ_EQUILIBRIUM_SNAPSHOT_HPP
