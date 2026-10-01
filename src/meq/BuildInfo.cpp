#include "BuildInfo.hpp"

// THE ONE INCLUDE THIS FILE EXISTS FOR. Every answer below is a question about
// MFEM's own preprocessor, and a precision or backend macro tested in a
// translation unit that has not included MFEM is not false -- it is UNDEFINED,
// which reads as false and takes the other branch without a diagnostic. See
// BuildInfo.hpp.
#include "mfem.hpp"

namespace meq
{

	char const *realPrecision()
	{
		/*
		 * ASKED OF THE MACRO AND NOT OF sizeof( mfem::real_t ).
		 *
		 * The two agree today and they answer different questions: `sizeof`
		 * measures a type this translation unit was compiled against, and the
		 * macro states what the LIBRARY was configured as. Where those can
		 * disagree -- a stale installed header against a rebuilt archive, which
		 * `CLAUDE.md`'s mtime trap produces silently -- the configured value is
		 * the one a reader of the output wants, and the disagreement is a
		 * defect to see rather than to paper over.
		 */
#ifdef MFEM_USE_SINGLE
		return "single";
#else
		return "double";
#endif
	}


	char const *deviceBackends()
	{
		/*
		 * SPELLED OUT RATHER THAN BUILT UP, because there are four reachable
		 * combinations and a std::string assembled from fragments would need a
		 * static to return a `char const *` from.
		 */
#if defined( MFEM_USE_CUDA ) && defined( MFEM_USE_HIP )
		return "cuda, hip";
#elif defined( MFEM_USE_CUDA )
		return "cuda";
#elif defined( MFEM_USE_HIP )
		return "hip";
#else
		return "none";
#endif
	}


	std::string buildDescription()
	{
		return std::string( realPrecision() ) + ", " + deviceBackends();
	}

}
