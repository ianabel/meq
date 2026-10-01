#ifndef MEQ_BUILD_INFO_HPP
#define MEQ_BUILD_INFO_HPP

#include <string>

/**
 * @file
 *
 * WHICH BUILD PRODUCED A RUN: THE PRECISION AND THE DEVICE BACKENDS, AS
 * STRINGS A FILE CAN CARRY.
 *
 * **MEQ SHIPS AS MORE THAN ONE BINARY AND ITS OUTPUT COULD NOT SAY WHICH.**
 * Precision is not a MEQ option and cannot be: `MFEM_PRECISION` decides
 * `real_t` throughout the library, so it is a property of the MFEM an install
 * links and MEQ inherits it through `MFEM_DIR`. A double install and a single
 * one therefore produce two binaries with the same name, the same `--version`
 * and the same output attributes, differing in the only thing that matters for
 * reading their answers — `DEVICE-PLAN.md` §7.5 measures the same equilibrium
 * at `2.2 × eps` on one machine and `1,260 × eps` on another, so *which build
 * wrote this file* is part of what the number means.
 *
 * **AND THE WHOLE DEVICE PROGRAMME IS A COMPARISON BETWEEN TWO BUILDS.**
 * `--device cuda` reproducing `--device cpu` is the instrument M-79, M-130,
 * M-187 and M-207 were root-caused with. An instrument whose two arms cannot
 * be told apart afterwards from what they wrote is not one.
 *
 * **THE IMPLEMENTATION IS IN THE `.cpp` AND THAT IS LOAD BEARING, NOT
 * TIDINESS.** These are answers about MFEM's own preprocessor, and
 * `CLAUDE.md` records the trap exactly: `tests/convergence/MachineFixedBoundary.cpp`
 * carried a `#ifdef MFEM_USE_SINGLE` in a translation unit that included no
 * MFEM header, so the guard compiled, ran, and **silently took the other
 * branch**. `mfem/config/config.hpp` is what defines those macros. Keeping the
 * header MFEM-free and the answers behind a function means no caller can
 * repeat that mistake by including this and forgetting the other — the same
 * reason `meq::Config`, `Profiles` and `Source` are MFEM-free in *Layout*,
 * applied to a question that is entirely MFEM's.
 */
namespace meq
{

	/// What `mfem::real_t` is in the linked MFEM: `"double"` or `"single"`.
	///
	/// Not a MEQ choice. `MFEM_PRECISION` is a cache variable of MFEM's own
	/// build and MEQ has no key for it, so this reports what was linked rather
	/// than what was asked for.
	char const *realPrecision();

	/// The device backends this build could configure, as a comma-separated
	/// list, or `"none"`.
	///
	/// **AVAILABLE IS NOT CONFIGURED.** This says what `--device` could be
	/// given, not what a particular run was given; `mfem::Device` is
	/// process-wide runtime state and the driver reports that separately.
	char const *deviceBackends();

	/// The two above, joined, for a one-line banner or a file attribute:
	/// `"double, cuda"`.
	std::string buildDescription();

}

#endif // MEQ_BUILD_INFO_HPP
