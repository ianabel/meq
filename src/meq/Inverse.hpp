#ifndef MEQ_INVERSE_HPP
#define MEQ_INVERSE_HPP

#include <memory>

#include "Coils.hpp"
#include "Config.hpp"
#include "GradShafranov.hpp"

/**
 * @file Inverse.hpp
 * A CONFIGURATION'S INVERSE PROBLEM, POSED ON A SOLVER, BY THE ONE ROUTE.
 * INVERSE-AND-RECONSTRUCTION-PLAN.md stage IV-5.
 *
 * The solver knows nothing of TOML and the configuration knows nothing of a
 * solver, so something has to translate `Free`, `[[circuits]]` and `[inverse]`
 * into setFreeCurrents() and the target calls. **IT IS THIS FUNCTION AND
 * NOTHING ELSE**: the driver calls it, and so does every test that pins the
 * driver against the library, so the two cannot wire a circuit or order the
 * targets differently and both converge -- the reason meq::Configuration is the
 * only reader of the schema.
 */
namespace meq
{
	/**
	 * @brief Make `config`'s `Free` currents unknowns of `solver` and add
	 *        `inverse`'s targets.
	 *
	 * The currents' seeds are what `coils` holds now, which is the
	 * configuration's own unless a caller has carried an answer forward.
	 * `coils` must be the set the solver's source carries; solve() refuses
	 * any other. Isoflux targets are added first, then X-point targets, in file
	 * order -- the order meq::GradShafranovSolver::targetResiduals() reports
	 * their rows in. `Gamma` applies to every free circuit.
	 *
	 * @throws whatever the solver's setters throw. A configuration that parsed
	 *         has already been checked for every pairing the setters refuse.
	 */
	void poseInverse( GradShafranovSolver &solver, std::shared_ptr< CoilSet > coils,
	                  CoilConfig const &config, InverseConfig const &inverse );
}

#endif
