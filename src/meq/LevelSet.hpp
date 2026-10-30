#ifndef MEQ_LEVELSET_HPP
#define MEQ_LEVELSET_HPP

#include "BoundaryShape.hpp"

#include "mfem.hpp"

#include <cmath>

/**
 * @file
 *
 * GAMMA AS A LEVEL SET, AND ITS GRADIENT, HANDED OUT TOGETHER.
 *
 * mfem's path families find each path's endpoint on Gamma by bisecting the
 * level set, and ExtensionBoundaryQuadrature() and ExtensionRegionQuadrature()
 * need that endpoint's derivative along the face. Given the level set's
 * gradient (`SetLevelSetGradient()`), the families differentiate the endpoint
 * in closed form and the sweeps take no finite-difference step at all.
 * Without it they central-difference the endpoint, which reads the
 * bisection's round-off divided by the step: on the rotating mirror, a 1/step
 * dependence of psi_h near 1e-9 (TODO, the extension's endpoint staircase).
 *
 * **THE GRADIENT MUST BE OF THE SAME FUNCTION.** Nothing in mfem can check
 * that, and a gradient of anything else gives a wrong derivative without an
 * error. So both halves are built here, in one place, and the driver and every
 * test that compares against it take them from here: a comparison between two
 * routes is exact only if both differentiate the same endpoints the same way.
 */

namespace meq
{
	/// A level set of Gamma, negative inside, and its exact gradient.
	struct LevelSet
	{
		mfem::PositionFunction value;
		mfem::VectorPositionFunction gradient;
	};

	/**
	 * The circle of @a radius about ( @a centreR, @a centreZ ), as a signed
	 * distance. With @a centreR = 0 it is `[boundary.exterior]`'s semicircle,
	 * whose flat side is the axis.
	 */
	inline LevelSet circleLevelSet( double centreR, double centreZ, double radius )
	{
		LevelSet out;
		out.value = [ centreR, centreZ, radius ]( mfem::Vector const &x )
		{
			return std::hypot( x( 0 ) - centreR, x( 1 ) - centreZ ) - radius;
		};
		out.gradient = [ centreR, centreZ ]( mfem::Vector const &x, mfem::Vector &gradient )
		{
			double const dr = x( 0 ) - centreR;
			double const dz = x( 1 ) - centreZ;
			double const distance = std::hypot( dr, dz );
			gradient.SetSize( 2 );
			gradient( 0 ) = distance > 0.0 ? dr/distance : 0.0;
			gradient( 1 ) = distance > 0.0 ? dz/distance : 0.0;
		};
		return out;
	}

	/**
	 * meq::BoundaryShape's radial gap and its gradient. The shape is held by
	 * pointer and must outlive the result, as it must outlive anything built
	 * from its level set.
	 */
	inline LevelSet shapeLevelSet( BoundaryShape const &shape )
	{
		BoundaryShape const *shapePointer = &shape;
		LevelSet out;
		out.value = [ shapePointer ]( mfem::Vector const &x )
		{
			return shapePointer->levelSet( x( 0 ), x( 1 ) );
		};
		out.gradient = [ shapePointer ]( mfem::Vector const &x, mfem::Vector &gradient )
		{
			double dRadius = 0.0, dZ = 0.0;
			shapePointer->levelSetGradient( x( 0 ), x( 1 ), dRadius, dZ );
			gradient.SetSize( 2 );
			gradient( 0 ) = dRadius;
			gradient( 1 ) = dZ;
		};
		return out;
	}
}

#endif // MEQ_LEVELSET_HPP
