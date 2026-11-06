#ifndef MEQ_RULEPOINT_HPP
#define MEQ_RULEPOINT_HPP

#include "mfem.hpp"

#include <cstdint>

namespace meq
{
	/**
	 * @brief Which of @a rule's points @a point is, by address; -1 for none.
	 *
	 * For a table tabulated at one rule's points and looked up from inside an
	 * integrator, which sees only `ElementTransformation::GetIntPoint()`: the
	 * point is the table's when it is one of the rule's own IntegrationPoint
	 * objects, and anything else -- another rule's point, a copy, a point on a
	 * stack -- is not, and the caller evaluates it instead.
	 *
	 * **IN BYTES, AND NOT `&point - &rule.IntPoint( 0 )`, WHICH WAS A CRASH.**
	 * Pointer subtraction is defined only within one array, so g++ divides the
	 * byte difference by sizeof( IntegrationPoint ) = 40 as an EXACT division,
	 * a multiply by 0xcccc...cd. A point at a higher address and not a multiple
	 * of 40 away came out as a garbage index; about half of those are negative
	 * as signed, pass an `i < GetNPoints()` check, and index far outside the
	 * table. Whether MFEM's rule for some other order sits above or below the
	 * table's is allocation order, and under OpenMP that is whichever thread's
	 * arena asked first -- so the read was out of bounds intermittently, and
	 * more often under load. RotatingMirrorAnisotropic's SIGSEGV in
	 * GradShafranovSolver::conductorAtLawPoint() inside Reconstruct(), found at
	 * an index register of 0x6666bbbbce21b000.
	 */
	inline int indexInRule( mfem::IntegrationRule const &rule,
	                        mfem::IntegrationPoint const &point )
	{
		int const count = rule.GetNPoints();
		if ( count == 0 )
			return -1;
		std::uintptr_t const address = reinterpret_cast<std::uintptr_t>( &point );
		std::uintptr_t const first = reinterpret_cast<std::uintptr_t>( &rule.IntPoint( 0 ) );
		std::uintptr_t const stride = sizeof( mfem::IntegrationPoint );
		if ( address < first )
			return -1;
		std::uintptr_t const bytes = address - first;
		if ( bytes % stride != 0 || bytes/stride >= static_cast<std::uintptr_t>( count ) )
			return -1;
		return static_cast<int>( bytes/stride );
	}
}

#endif // MEQ_RULEPOINT_HPP
