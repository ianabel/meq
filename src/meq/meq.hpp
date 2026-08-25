#ifndef MEQ_HPP
#define MEQ_HPP
/*
 * Main include file for MEQ
 */

#include <vector>
#include <memory>
#include <tuple>

#include "mfem.hpp"
#include "SplineInterpolant.hpp"

namespace meq {

	using RealScalarField = std::function<double( const mfem::Vector & )>;
	using RealVectorField = std::function<void( const mfem::Vector &, mfem::Vector & )>;

};

#include "utility.hpp"
#include "model.hpp"
#include "config.hpp"
#include "solver.hpp"

#endif // MEQ_HPP
