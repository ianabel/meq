#include "PointLocator.hpp"

/*
 * THE DEVICE HALF OF meq::PointLocator, AND IT IS A SEPARATE TRANSLATION UNIT
 * FOR A REASON THAT IS NOT STYLE.
 *
 * mfem::forall DOES NOT COMPILE UNDER g++. general/backends.hpp:144 is a
 * static_assert reading "The calling function needs to be compiled with
 * CUDA/HIP language!", so a device kernel written anywhere in MEQ's ordinary
 * sources is a compile error rather than a slow path. MEQ compiles with g++
 * throughout -- that is recorded in CLAUDE.md as a property of the build and it
 * is what this file works around: CMake marks THIS FILE ONLY as CUDA when the
 * linked MFEM has it, and every other file in src/meq is untouched.
 *
 * WHAT IT DOES NOT DO, so nobody reads the file name as a capability:
 *
 *   IT LOCATES, IT DOES NOT INTERPOLATE. Evaluating the field needs
 *   FiniteElement::CalcShape, which is a host virtual taking an
 *   mfem::IntegrationPoint. A device evaluation needs a MFEM_HOST_DEVICE basis
 *   and that is DEVICE-PLAN.md item 3's blocker, not this file's -- the two
 *   share it rather than competing, and until it is answered the located points
 *   come back to the host to be read.
 *
 *   AFFINE ELEMENTS ONLY. The curved fallback is an ElementTransformation
 *   Newton, which is not a device object. PointLocator::locate() refuses the
 *   device for a mesh with any curved or degenerate element rather than
 *   splitting one query across two locators.
 *
 * THE ARITHMETIC IS NOT REPEATED HERE. locateAffine() is in the header,
 * MFEM_HOST_DEVICE, POD in and POD out, and is compiled twice -- once by g++
 * into the host loop and once by nvcc into the kernel below. That is the only
 * arrangement in which the two halves cannot drift apart, which matters more
 * than it sounds: a locator that disagrees with itself between builds is a warm
 * start that puts the solve somewhere else on a machine with a GPU.
 */

namespace meq
{

/*
 * THE GUARD IS MFEM'S OWN AND MEQ INVENTS NOTHING. general/backends.hpp:54
 * defines MFEM_USE_CUDA_OR_HIP_LANG exactly when the translation unit is being
 * compiled by the device-language compiler, and :134 turns its absence into the
 * static_assert this file exists to avoid. So the test is "was I compiled as a
 * device TU", not "is CUDA enabled" -- which is the difference between a build
 * that works under HIP for free and one that names a vendor. cuDSS pins MEQ to
 * CUDA today, so nothing here is exercised under HIP; the macro costs nothing
 * over naming CUDA and is the reason this file never has to change if that
 * stops being true.
 */
#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	bool locateOnDevice( ElementMap const *maps, int elements,
	                     int const *binStart, int const *binElement,
	                     mfem::real_t rMin, mfem::real_t zMin,
	                     mfem::real_t dR, mfem::real_t dZ,
	                     int binsR, int binsZ,
	                     mfem::real_t const *positions, int points,
	                     int *element, unsigned int *code,
	                     mfem::real_t *reference, mfem::real_t *distance )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) || points == 0
		     || elements == 0 )
			return false;

		/*
		 * The five arrays the kernel reads are POD and are registered with the
		 * memory manager here rather than being members, so nothing in this
		 * locator holds a device allocation between calls. That is deliberate:
		 * M-79 and M-130 are both MEQ holding a long-lived alias onto something
		 * whose device copy moved underneath it, and the cheapest way not to
		 * write that defect a third time is not to keep the buffer.
		 */
		int const bins = binsR*binsZ;

		mfem::Array<ElementMap> mapBuf( const_cast<ElementMap *>( maps ), elements );
		mfem::Array<int> startBuf( const_cast<int *>( binStart ), bins + 1 );
		mfem::Array<int> listBuf( const_cast<int *>( binElement ),
		                          binStart[ bins ] );
		mfem::Vector posBuf( const_cast<mfem::real_t *>( positions ), 2*points );
		mfem::Array<int> elementBuf( element, points );
		mfem::Array<unsigned int> codeBuf( code, points );
		mfem::Vector referenceBuf( reference, 2*points );
		mfem::Vector distanceBuf( distance, points );

		auto const *d_maps  = mapBuf.Read();
		auto const *d_start = startBuf.Read();
		auto const *d_list  = listBuf.Read();
		auto const *d_pos   = posBuf.Read();
		auto *d_element     = elementBuf.Write();
		auto *d_code        = codeBuf.Write();
		auto *d_reference   = referenceBuf.Write();
		auto *d_distance    = distanceBuf.Write();

		constexpr mfem::real_t tolerance = 1.0e-8;
		int const maxRing = binsR > binsZ ? binsR : binsZ;

		/*
		 * ONE POINT PER THREAD, and the ring walk is inside the thread rather
		 * than being a second kernel. The candidate lists are short -- about
		 * one element per bin by construction -- so the divergence between a
		 * thread that answers on ring 0 and one that walks outwards is bounded
		 * by how far the query cloud reaches outside the mesh, which is the
		 * band an adaptive cycle adds and not the mesh's diameter.
		 */
		mfem::forall( points, [ = ] MFEM_HOST_DEVICE ( int p )
		{
			mfem::real_t const r = d_pos[ p ];
			mfem::real_t const z = d_pos[ points + p ];

			int i = static_cast<int>( ( r - rMin )/dR );
			int j = static_cast<int>( ( z - zMin )/dZ );
			i = i < 0 ? 0 : ( i > binsR - 1 ? binsR - 1 : i );
			j = j < 0 ? 0 : ( j > binsZ - 1 ? binsZ - 1 : j );

			int found = -1;
			mfem::real_t foundXi = 0.0, foundEta = 0.0;

			for ( int ring = 0; ring <= maxRing && found < 0; ++ring )
			{
				int const i0 = ( i - ring ) < 0 ? 0 : i - ring;
				int const i1 = ( i + ring ) > binsR - 1 ? binsR - 1 : i + ring;
				int const j0 = ( j - ring ) < 0 ? 0 : j - ring;
				int const j1 = ( j + ring ) > binsZ - 1 ? binsZ - 1 : j + ring;

				for ( int jj = j0; jj <= j1 && found < 0; ++jj )
					for ( int ii = i0; ii <= i1 && found < 0; ++ii )
					{
						if ( ring > 0 && ii != i0 && ii != i1 && jj != j0 && jj != j1 )
							continue;

						int const b = jj*binsR + ii;
						for ( int k = d_start[ b ]; k < d_start[ b + 1 ]; ++k )
						{
							int const e = d_list[ k ];
							mfem::real_t xi = 0.0, eta = 0.0;
							mfem::real_t const overshoot =
								locateAffine( d_maps[ e ], r, z, xi, eta );
							if ( overshoot <= tolerance )
							{
								found = e;
								foundXi = xi;
								foundEta = eta;
								break;
							}
						}
					}
			}

			d_element[ p ] = found;
			d_code[ p ] = found < 0 ? 2u : 0u;
			d_reference[ 2*p ]     = foundXi;
			d_reference[ 2*p + 1 ] = foundEta;
			d_distance[ p ] = 0.0;
		} );

		/*
		 * THE RESULTS ARE PULLED BACK EXPLICITLY, because everything that reads
		 * them is host arithmetic -- the L2 projection in meq::FieldTransfer
		 * reads through operator() and operator[], which neither sync nor
		 * invalidate. This is the same species as the four sites behind M-79
		 * and the two behind M-130, and the comment in WarmStart.cpp beside the
		 * three HostRead() calls is the longer form. A device write that is
		 * never pulled back is a warm start built from whatever the host buffer
		 * last held, and under CUDA that is silent.
		 */
		elementBuf.HostRead();
		codeBuf.HostRead();
		referenceBuf.HostRead();
		distanceBuf.HostRead();

		return true;
	}

#else

	bool locateOnDevice( ElementMap const *, int,
	                     int const *, int const *,
	                     mfem::real_t, mfem::real_t,
	                     mfem::real_t, mfem::real_t,
	                     int, int,
	                     mfem::real_t const *, int,
	                     int *, unsigned int *,
	                     mfem::real_t *, mfem::real_t * )
	{
		// No CUDA translation unit in this build. The host loop is the whole
		// implementation and PointLocator::locate() runs it unconditionally.
		return false;
	}

#endif

}
