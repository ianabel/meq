#ifndef MEQ_THREADING_HPP
#define MEQ_THREADING_HPP

#include "mfem.hpp"

#include <initializer_list>

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
 *
 * **AND THE PAIR IS NOT THE WHOLE CONDITION: A CONFIGURED `mfem::Device` TAKES
 * THESE LOOPS BACK TO ONE THREAD, AT RUN TIME.** MFEM registers buffers lazily
 * in its memory manager once a Device exists -- `Memory<T>::MakeAlias()`,
 * `SetDeviceMemoryType()` and the device arms of `Read`/`Write`/`ReadWrite` --
 * into `general/mem_manager.cpp`'s `maps->memories` and `maps->aliases`, which
 * are bare `std::unordered_map`s carrying no mutex, no critical section and no
 * atomic. Every loop below reaches one: `GridFunction::GetVectorValue()` alone
 * goes through `GetElementVDofs` to `Ordering::DofsToVDofs`, which allocates an
 * `Array<int>` per call, which registers.
 *
 * **MEASURED RATHER THAN ARGUED, AND THE INSTRUMENT WAS A THREAD COUNT.**
 * `machine-f-diiid` under `--device debug` with `AssemblyMode = "serial"` --
 * so none of MFEM's own assembly threading in the process at all -- runs clean
 * at `OMP_NUM_THREADS=1` and **segfaults at 8**, in
 * `MemoryManager::Insert()` beneath `CriticalPointFinder::sweep()`'s own
 * `GOMP_parallel` frame. `--device cuda` happens not to fire on the same
 * configuration, because `GetHostMemoryType()` is `HOST` there and leaves
 * owners unregistered where `debug`'s `HOST_DEBUG` registers them; that is
 * luck about how often the map is touched, not safety.
 *
 * **THIS IS MFEM'S OWN ANSWER FOR ITS OWN LOOPS, COPIED.**
 * `DarcyHybridization::ThreadHostLoops()` opens
 * `if ( Device::Allows( Backend::DEVICE_MASK ) ) { return false; }` for exactly
 * this reason. The library turns its host threading off under a Device; MEQ's
 * loops are MEQ's to turn off, and nothing else would have.
 *
 * **IT COSTS NO ANSWER AND IT IS NOT A WORKAROUND FOR A DEVICE RUN BEING
 * SLOW.** These loops are bit exact against the serial build by construction --
 * that is the claim MEQ_OMP's compile-time gate exists to keep checkable -- so
 * a Device run computes what it always did and does it on one thread. Under
 * residency it is what is wanted anyway: host threads are work the device is
 * supposed to be doing. DEVICE-PLAN.md section 3.3, item 7.
 */

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
	#define MEQ_PRAGMA( x ) _Pragma( #x )
	/// One OpenMP directive, or nothing at all on a build without the pair.
	#define MEQ_OMP( x ) MEQ_PRAGMA( omp x )
#else
	#define MEQ_OMP( x )
#endif

namespace meq
{
	/**
	 * @brief May MEQ's own element loops take more than one thread?
	 *
	 * False whenever an mfem::Device is configured, whatever the backend --
	 * `cuda`, `debug` and any other all reach the same unguarded registration
	 * maps. `--device cpu` is NOT a device by this test and keeps its threads,
	 * which is what leaves every host measurement in MEASUREMENTS.md intact.
	 *
	 * Safe to call before any Device is constructed: mfem::Device::Allows reads
	 * the singleton's backend mask, which is the CPU backend until one is.
	 */
	inline bool hostLoopsAreSafe()
	{
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
		return !mfem::Device::Allows( mfem::Backend::DEVICE_MASK );
#else
		return false;
#endif
	}

	/**
	 * @brief May a SETUP loop take more than one thread, under a Device too?
	 *
	 * hostLoopsAreSafe() says no to every loop under a Device, and for a loop
	 * that belongs on the card that is the right answer: it should be a kernel.
	 * A setup loop is a different case. The warm start's search of a curved
	 * guess mesh is host work on every arm, as meshing is, and holding it to
	 * one core cost an A100 run 27% of its wall (M-319).
	 *
	 * **IT THREADS UNDER A Device ONLY WHERE THE LOOP CANNOT REACH MFEM'S
	 * MEMORY BOOKKEEPING AT ALL.** That is narrower than the race upstream
	 * documents above mfem::MemoryManager, where `cuda` merely touches the maps
	 * less often. Here they are not touched:
	 *
	 *   * the host memory type is plain HOST, so every Vector, Array and
	 *     DenseMatrix a thread allocates is new[] and delete[] and never enters
	 *     the maps. `debug`'s HOST_DEBUG registers every one, and is refused;
	 *   * none of @a shared -- the host pointers of every object the loop
	 *     reads or writes through an MFEM accessor -- is registered, so
	 *     Memory::Read() returns the host pointer before it reaches the maps or
	 *     writes its own flags. A registered one keeps the loop serial.
	 *
	 * Beyond those two it relies on nothing inside the loop asking for DEVICE
	 * memory, which would register; the loops it gates are element
	 * transformations, shape functions and dense algebra. A null pointer in
	 * @a shared is an object that does not exist, and is skipped.
	 */
	inline bool setupLoopMayThread( std::initializer_list<void const *> shared )
	{
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return true;
		if ( mfem::MemoryManager::GetHostMemoryType() != mfem::MemoryType::HOST )
			return false;
		for ( void const *pointer : shared )
			if ( pointer != nullptr
			     && ( mfem::mm.IsKnown( pointer ) || mfem::mm.IsAlias( pointer ) ) )
				return false;
		return true;
#else
		static_cast<void>( shared );
		return false;
#endif
	}
}

/**
 * One parallel element loop, gated on BOTH conditions -- the build's at compile
 * time and the Device's at run time, through OpenMP's own `if` clause, which
 * runs the region on a team of one rather than skipping it.
 *
 * The `atomic` and `critical` directives inside these loops stay plain
 * MEQ_OMP(): both are correct on a one-thread team, and neither takes an `if`.
 */
#define MEQ_OMP_ELEMENT_LOOP \
	MEQ_OMP( parallel for schedule( dynamic ) if ( ::meq::hostLoopsAreSafe() ) )

#endif // MEQ_THREADING_HPP
