#ifndef MEQ_HPP
#define MEQ_HPP

/*
 * Umbrella header for the meq library.
 *
 * Include this to get everything that is currently ported. It is deliberately
 * short: the solver core is still being moved onto MFEM 4.9.1's DarcyForm, and
 * a header only appears here once it compiles.
 *
 * Not yet included, because they still target the Waterloo HDGBilinearForm API
 * from MFEM 4.5.1 and do not compile:
 *
 *     meq/GradShafranov.hpp   the HDG assembly     -- stage 2
 *     meq/Solution.hpp        solution container   -- stage 2
 *     meq/Estimator.hpp       residual estimator   -- stage 6
 *
 * See CLAUDE.md for the stage plan.
 */

#include "meq/Config.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#endif // MEQ_HPP
