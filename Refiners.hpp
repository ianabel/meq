#ifndef REFINERS_HPP
#define REFINERS_HPP

#include "mfem.hpp"

namespace mfem {
/** @brief Mesh refinement using Maximum Marking

    This class uses the given ErrorEstimator to estimate local element errors
    and then marks for refinement all elements i such that loc_err_i > threshold.
    The threshold is computed as
    \code
       threshold = gamma * Max_{i in Mesh} { local_error_i }
    \endcode
	 where gamma is a configurable parameter
*/
class MaximumMarkingRefiner : public MeshOperator
{
protected:
   ErrorEstimator &estimator;

   long   max_elements;

   double threshold;
   long num_marked_elements;

   Array<Refinement> marked_elements;
   long current_sequence;

   int non_conforming;
   int nc_limit;

   /** @brief Apply the operator to the mesh.
       @return STOP if a stopping criterion is satisfied or no elements were
       marked for refinement; REFINED + CONTINUE otherwise. */
   virtual int ApplyImpl(Mesh &mesh);

public:
   /// Construct a ThresholdRefiner using the given ErrorEstimator.
   MaximumMarkingRefiner(ErrorEstimator &est);

   // default destructor (virtual)

   void SetMaxElements(long max_elem) { max_elements = max_elem; }

   /// Use nonconforming refinement, if possible (triangles, quads, hexes).
   void PreferNonconformingRefinement() { non_conforming = 1; }

   /** @brief Use conforming refinement, if possible (triangles, tetrahedra)
       -- this is the default. */
   void PreferConformingRefinement() { non_conforming = -1; }

   /** @brief Set the maximum ratio of refinement levels of adjacent elements
       (0 = unlimited). */
   void SetNCLimit(int nc_limit)
   {
      MFEM_ASSERT(nc_limit >= 0, "Invalid NC limit");
      this->nc_limit = nc_limit;
   }

   /// Get the number of marked elements in the last Apply() call.
   long GetNumMarkedElements() const { return num_marked_elements; }

   /// Get the threshold used in the last Apply() call.
   double GetThreshold() const { return threshold; }

   /// Reset the associated estimator.
	virtual void Reset() {
		estimator.Reset();
		current_sequence = -1;
		num_marked_elements = 0;
	};
};

/** @brief Mesh refinement using Dörfler Marking

    This class uses the given ErrorEstimator to estimate local element errors
    and then marks for refinement a set M such that
	 \code
       gamma Sum_{i in M} local_error_i^2 <= Sum_{i in Mesh} local_error_i^2
    \endcode
	 where gamma is a configurable parameter.
	 We use the QuickMark algorithm of Pfeiler and Praetorius ( arXiv:1907.13078 )
	 to compute the set M.
*/
class DoerflerMarkingRefiner : public MeshOperator
{
protected:
   ErrorEstimator &estimator;

   long   max_elements;

   double threshold;
   long num_marked_elements;

   Array<Refinement> marked_elements;
   long current_sequence;

   int non_conforming;
   int nc_limit;

   /** @brief Apply the operator to the mesh.
       @return STOP if a stopping criterion is satisfied or no elements were
       marked for refinement; REFINED + CONTINUE otherwise. */
   virtual int ApplyImpl(Mesh &mesh);

public:
   /// Construct a ThresholdRefiner using the given ErrorEstimator.
   DoerflerMarkingRefiner(ErrorEstimator &est);

   // default destructor (virtual)

   void SetMaxElements(long max_elem) { max_elements = max_elem; }

   /// Use nonconforming refinement, if possible (triangles, quads, hexes).
   void PreferNonconformingRefinement() { non_conforming = 1; }

   /** @brief Use conforming refinement, if possible (triangles, tetrahedra)
       -- this is the default. */
   void PreferConformingRefinement() { non_conforming = -1; }

   /** @brief Set the maximum ratio of refinement levels of adjacent elements
       (0 = unlimited). */
   void SetNCLimit(int nc_limit)
   {
      MFEM_ASSERT(nc_limit >= 0, "Invalid NC limit");
      this->nc_limit = nc_limit;
   }

   /// Get the number of marked elements in the last Apply() call.
   long GetNumMarkedElements() const { return num_marked_elements; }

   /// Get the threshold used in the last Apply() call.
   double GetThreshold() const { return threshold; }

   /// Reset the associated estimator.
	virtual void Reset() {
		estimator.Reset();
		current_sequence = -1;
		num_marked_elements = 0;
	};
}

}
#endif // REFINERS_HPP

