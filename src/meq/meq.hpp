#ifndef MEQ_HPP
#define MEQ_HPP

/*
 * Umbrella header for the meq library.
 *
 * Include this to get everything that is currently ported. A header appears here
 * once it compiles.
 *
 * Not yet included, because it still targets the Waterloo HDGBilinearForm API
 * from MFEM 4.5.1 and does not compile:
 *
 *     meq/Estimator.hpp       residual estimator   -- stage 6
 *
 * meq/Solution.hpp is not missing but gone: DarcyForm owns the spaces and the
 * block structure it used to wrap, and GradShafranovSolver owns the rest.
 *
 * See CLAUDE.md for the stage plan.
 */

#include "meq/Config.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#endif // MEQ_HPP
