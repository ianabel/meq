#ifndef MEQ_THREADING_HPP
#define MEQ_THREADING_HPP

#include "mfem.hpp"

/**
 * @file
 *
 * MEQ'S OWN OPENMP DIRECTIVES, GATED IN ONE PLACE.
 *
 * MEQ has element loops of its own -- the bordered Newton's five plasma
 * assemblers in GradShafranov.cpp and CriticalPointFinder::sweep() -- which are
 * MEQ's parallelism rather than MFEM's. They are written through MEQ_OMP() so
 * that the condition they are gated on is stated once.
 *
 * **`MFEM_USE_OPENMP` ALONE IS NOT THE CONDITION, AND THAT IS THE POINT OF THE
 * PAIR.** `MFEM_THREAD_SAFE` is what removes the mutable scratch from
 * mfem::FiniteElement::CalcShape and its kind, and what puts the OpenMP lock
 * around mfem::IntegrationRules::Get's lazy generation. Without it an element
 * loop reading shape functions from several threads is a silent wrong answer
 * rather than a build failure. mfem::DarcyHybridization::SetAssemblyMode(
 * Threaded ) refuses the same pair for the same reason, and
 * tests/unit/SolverContract.cpp guards on it.
 *
 * With either absent the directives vanish and every loop is the serial loop it
 * was -- which is also what makes the bit-exactness claim checkable, since the
 * serial build IS the reference.
 *
 * WHAT THE DIRECTIVES DO NOT EXCUSE. Three hazards are the caller's and no
 * macro can see them:
 *
 *   * `mfem::Mesh::GetElementTransformation( int )` hands out the mesh's own
 *     shared scratch. So does `Mesh::GetElementSize( int, int )`, which calls
 *     it, and `Mesh::GetBdrFaceTransformations( int )`. Take the
 *     caller-allocated overload.
 *   * `mfem::FiniteElementSpace::GetElementDofs( int, Array<int> & )` writes
 *     into `FiniteElementSpace::DoFTrans`, a mutable member of the SPACE. Take
 *     the three-argument overload with a local `mfem::DofTransformation`.
 *   * Scratch declared above a loop is shared by every thread in it, and
 *     `mfem::Vector::SetSize` is a reallocation -- so the failure is memory
 *     corruption, not a stale read. Declare it inside the body.
 */

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
	#define MEQ_PRAGMA( x ) _Pragma( #x )
	/// One OpenMP directive, or nothing at all on a build without the pair.
	#define MEQ_OMP( x ) MEQ_PRAGMA( omp x )
#else
	#define MEQ_OMP( x )
#endif

#endif // MEQ_THREADING_HPP
