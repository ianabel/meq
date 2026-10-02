#ifndef MEQ_POINTLOCATOR_HPP
#define MEQ_POINTLOCATOR_HPP

#include "mfem.hpp"

#include <vector>

/*
 * Locating arbitrary points in a 2D mesh, and evaluating a field at them.
 *
 * WHY THIS EXISTS WHEN mfem::FindPointsGSLIB DOES THE SAME JOB. Three reasons,
 * and the third is the one that decides it -- MEASUREMENTS.md M-175:
 *
 *   PRECISION   gslib has no floating-point abstraction at all: no `real`
 *               typedef, no build switch, 1,153 `double`s across 28 files of C.
 *               MFEM's own reimplementation under fem/gslib/ is 6,019 lines
 *               with 542 `double`s and zero real_t. So MFEM_USE_GSLIB and
 *               MFEM_PRECISION = single are mutually exclusive, and the warm
 *               start is what MEQ loses -- which is load-bearing rather than a
 *               convenience, since every bordered free-boundary case in the
 *               tree fails cold even in double.
 *   COST        gslib's simplex path splits every triangle into THREE quads
 *               (fem/gslib.cpp, NEsplit = 3) and resamples the field onto
 *               them, so machine-f-diiid's 4848 elements become 14,544 inside
 *               it. MEQ is triangles on every mesh it ships.
 *   THE DEVICE  MFEM's device findpts is gated on `tensor_product_only`,
 *               which requires QUADRILATERAL or HEXAHEDRON
 *               (fem/gslib.cpp:1388). MEQ is triangles, so that gate is FALSE
 *               and those 6,019 lines of device kernels are unreachable here
 *               whatever the precision. A device warm start for MEQ has to be
 *               MEQ's own or it does not exist.
 *
 * WHAT THIS IS NOT. It is not a gslib replacement -- it is a replacement for
 * MEQ's USE of gslib: serial, 2D, one field at a time, no MPI, no communicator,
 * no crystal router. Of the external symbols MFEM calls, sarray_transfer (18
 * sites), comm_init (7) and crystal_init/free (7) are exactly that router, and
 * a build with MFEM_USE_MPI = NO never enters it. What is left is a hash, a
 * bounding box and an evaluation, which is what this file is.
 *
 * THE SEARCH IS A UNIFORM BIN GRID over element bounding boxes, in CSR form.
 * meq::GridSampler solves the same location problem the other way round --
 * element bounding box to GRID INDEX RANGE -- because its targets are a regular
 * grid and the loop inverts. These targets are a target mesh's Gauss points and
 * do not, so the index has to exist. Everything downstream of finding the
 * candidate element is GridSampler's: the affine inverse for a straight-sided
 * triangle is a 2x2 solve measured at 69x and 81x against TransformBack (M-93),
 * and the curved path keeps Newton.
 *
 * WHY THE KERNEL IS IN THE HEADER AND TAKES RAW POINTERS. mfem::forall does not
 * compile under g++ -- general/backends.hpp:144 is a static_assert reading "The
 * calling function needs to be compiled with CUDA/HIP language!" -- so a device
 * kernel from MEQ needs ITS OWN translation unit compiled by nvcc. locateAffine
 * below is therefore MFEM_HOST_DEVICE, POD in and POD out, and is compiled
 * TWICE: once by g++ into PointLocator.cpp's host loop and once by nvcc into
 * PointLocatorDevice.cpp's mfem::forall. One expression of the arithmetic, two
 * compilers, which is the only arrangement that cannot drift.
 */

namespace meq
{

	/**
	 * The per-element data the search kernel needs, as plain arithmetic types
	 * so the same struct can be read on a host or a device.
	 *
	 * `inverse` is the 2x2 inverse of the affine map's Jacobian and `originR`,
	 * `originZ` its first vertex, so a physical point maps to reference
	 * coordinates by one subtraction and one 2x2 multiply. `affine` is false
	 * for a curved element or a degenerate one, and such an element is skipped
	 * by the kernel and left to the host's Newton route.
	 */
	/*
	 * NO DEFAULT MEMBER INITIALISERS, AND THAT IS A REQUIREMENT RATHER THAN A
	 * STYLE. mfem::Array<T> static_asserts "type T must be trivial"
	 * (general/array.hpp:61), and a default member initialiser makes a struct
	 * non-trivially-default-constructible -- so writing `= 0.0` here is a
	 * compile error in the device translation unit and nowhere else, which is
	 * a confusing place to meet it. Every construction below is value
	 * initialised (`ElementMap{}`), which zeroes an aggregate.
	 */
	struct ElementMap
	{
		mfem::real_t originR;
		mfem::real_t originZ;
		mfem::real_t inverse[ 4 ];
		bool affine;
	};

	/**
	 * Reference coordinates of @a r, @a z in the affine element @a map, and
	 * how far outside the reference triangle they fall.
	 *
	 * Returns the overshoot: zero or less when the point is inside, and
	 * otherwise the largest violation of the three barycentric constraints, so
	 * a caller comparing overshoots across candidates picks the element the
	 * point is least outside. That is what makes a point on a shared face land
	 * in one of its two elements rather than in neither.
	 *
	 * MFEM_HOST_DEVICE, POD only, no MFEM object: this is the one expression of
	 * the arithmetic and both the host loop and the device kernel call it.
	 */
	MFEM_HOST_DEVICE inline mfem::real_t locateAffine( ElementMap const &map,
	                                                   mfem::real_t r, mfem::real_t z,
	                                                   mfem::real_t &xi, mfem::real_t &eta )
	{
		mfem::real_t const dr = r - map.originR;
		mfem::real_t const dz = z - map.originZ;
		xi  = map.inverse[ 0 ]*dr + map.inverse[ 1 ]*dz;
		eta = map.inverse[ 2 ]*dr + map.inverse[ 3 ]*dz;

		// The three constraints of the reference triangle, as violations.
		mfem::real_t const a = -xi;
		mfem::real_t const b = -eta;
		mfem::real_t const c = xi + eta - mfem::real_t( 1.0 );
		mfem::real_t worst = a > b ? a : b;
		return worst > c ? worst : c;
	}

	/**
	 * Finds points in a mesh and evaluates fields at them.
	 *
	 * Setup is done once per mesh and is the expensive part, so a caller
	 * locating several point clouds in one mesh, or evaluating several fields
	 * at one cloud, should keep one of these.
	 *
	 * THE CODES ARE gslib'S, DELIBERATELY, because meq::FieldTransfer was
	 * written against them and a silent change of convention in a warm start is
	 * a solve that begins somewhere else: 0 found strictly inside, 1 found by
	 * snapping to the nearest element within tolerance, 2 not found.
	 */
	class PointLocator
	{
		public:

			/// @param mesh  the mesh to locate in. Borrowed, and must outlive
			///              this; the element indices reported refer to it.
			explicit PointLocator( mfem::Mesh &mesh );

			/**
			 * Locates every point in @a positions, which is `byNODES` -- every
			 * R, then every z -- matching the ordering meq::FieldTransfer hands
			 * gslib and for the same reason: the two orderings differ silently
			 * and a transposed cloud gives a plausible wrong answer.
			 */
			void locate( mfem::Vector const &positions );

			/// Evaluates @a field at the located points, in the order they were
			/// given. Entries whose code is 2 are left at zero; the caller
			/// decides what a missing point means, as gslib's callers must.
			void interpolate( mfem::GridFunction const &field,
			                  mfem::Vector &values ) const;

			/// 0 interior, 1 snapped, 2 not found -- one per point.
			mfem::Array<unsigned int> const &codes() const { return code; }

			/// How far each point moved to land in its element; 0 when interior.
			mfem::Vector const &distances() const { return distance; }

			/// The element each point landed in, or -1.
			mfem::Array<int> const &elements() const { return element; }

			/// How many points were not found.
			int missed() const { return misses; }

			/// Whether the device path was taken by the last locate().
			bool usedDevice() const { return deviceUsed; }

			/**
			 * The bin grid itself, for a caller that wants a spatial
			 * NEIGHBOURHOOD rather than the element containing a point.
			 *
			 * **EXPOSED SO THAT meq::AxisSearch DOES NOT BUILD A SECOND
			 * ONE.** It wants "the elements near where the axis was last
			 * time", which is a block of bins and is O( 1 ) to name; what it
			 * does NOT want is `locate()`, which answers a different
			 * question. The grid is a constant of the MESH -- about one
			 * element per bin by construction -- so it is built once and is
			 * `DEVICE-PLAN.md` section 3.3's phase B, and a second copy of it
			 * would be two uploads of one table, attributable to neither.
			 *
			 * The pointers are HOST pointers and are valid for this object's
			 * lifetime. A caller that wants them on a device uploads them; a
			 * caller that wants to walk a bin block on the host, as the axis
			 * search does when choosing its neighbourhood, reads them here.
			 */
			struct BinGrid
			{
				int const *start;    ///< `binsR*binsZ + 1` prefix offsets
				int const *element;  ///< element ids, grouped by bin
				mfem::real_t rMin, zMin, dR, dZ;
				int binsR, binsZ;
				/// False when any element is not a straight-sided triangle,
				/// which is the same refusal locate() makes.
				bool affine;
			};

			/// @copydoc BinGrid
			BinGrid grid() const;

		private:

			void buildBins();

			mfem::Mesh &meshRef;

			// The bin grid, in CSR form so it is one allocation per array and
			// can be read by a device kernel without an indirection per bin.
			std::vector<ElementMap> maps;
			mfem::Array<int> binStart;      ///< size bins + 1
			mfem::Array<int> binElement;    ///< element ids, grouped by bin
			mfem::real_t rMin = 0.0, zMin = 0.0, dR = 1.0, dZ = 1.0;
			int binsR = 1, binsZ = 1;
			bool anyCurved = false;

			// Results of the last locate().
			mfem::Array<int> element;
			mfem::Array<unsigned int> code;
			mfem::Vector reference;         ///< 2 per point
			mfem::Vector distance;
			int points = 0;
			int misses = 0;
			bool deviceUsed = false;
	};

	/**
	 * The device half of locate(), in its own translation unit because
	 * mfem::forall will not compile under g++.
	 *
	 * Returns false when there is no device path to take -- no mfem::Device
	 * configured, or MEQ built without the CUDA translation unit -- and the
	 * caller then runs the host loop. It handles AFFINE elements only; a mesh
	 * with any curved element is refused here and located on the host, because
	 * the curved route is an mfem::ElementTransformation Newton and that is not
	 * a device object.
	 */
	bool locateOnDevice( ElementMap const *maps, int elements,
	                     int const *binStart, int const *binElement,
	                     mfem::real_t rMin, mfem::real_t zMin,
	                     mfem::real_t dR, mfem::real_t dZ,
	                     int binsR, int binsZ,
	                     mfem::real_t const *positions, int points,
	                     int *element, unsigned int *code,
	                     mfem::real_t *reference, mfem::real_t *distance );

}

#endif
