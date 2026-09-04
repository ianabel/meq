#ifndef MEQ_HPP
#define MEQ_HPP

/*
 * Umbrella header for the MEQ library: every public header, in one include.
 *
 * IT USED TO CARRY FIVE OF SEVENTEEN AND SAY NOTHING ABOUT IT, which is the
 * worst of the three states this file could be in -- a consumer who included it
 * got less than a third of the library and had no way to tell. It also claimed
 * that meq/Estimator.hpp "still targets the Waterloo HDGBilinearForm API from
 * MFEM 4.5.1 and does not compile". That stopped being true at stage 6, and the
 * header contains no reference to that API at all.
 *
 * So the rule now is the one the name implies: A HEADER APPEARS HERE UNLESS
 * THERE IS A REASON IN THIS COMMENT WHY NOT, and there are currently none.
 *
 * TWO OF THESE DO NOT NEED MFEM AND CAN BE INCLUDED ON THEIR OWN.
 * meq/Zernike.hpp and meq/SurfaceFit.hpp are plain doubles in and coefficients
 * out, deliberately, so that they are unit-testable without the finite element
 * library and so that continuous integration -- which cannot obtain the MFEM
 * branch MEQ needs, see INSTALL.md -- can build and test them. Including them
 * through this header is fine and costs nothing beyond the MFEM everything else
 * here already pulls in; a consumer who wants only the disc basis should
 * include it directly instead.
 *
 * meq/Solution.hpp is not missing but gone: DarcyForm owns the spaces and the
 * block structure it used to wrap, and GradShafranovSolver owns the rest.
 * v0-legacy has the original.
 */

#include "meq/BoundaryShape.hpp"
#include "meq/Config.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/Estimator.hpp"
#include "meq/Field.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Sampler.hpp"
#include "meq/Source.hpp"
#include "meq/SourceFactory.hpp"
#include "meq/SurfaceAverage.hpp"
#include "meq/SurfaceFit.hpp"
#include "meq/WarmStart.hpp"
#include "meq/Zernike.hpp"

#endif // MEQ_HPP
