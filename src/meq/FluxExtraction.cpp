#include "FluxExtraction.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "GradShafranov.hpp"
#include "SurfaceAverage.hpp"

namespace meq
{

	namespace
	{
		double const pi = 3.141592653589793238462643383279;
		double const twoPi = 2.0*pi;

		std::string number( double value )
		{
			std::ostringstream out;
			out.precision( 6 );
			out << value;
			return out.str();
		}

		/// The enclosed volume and cross-section area of one fitted surface, by
		/// Green's theorem on its own nodes, with dz/dtheta POINTWISE from rho'.
		///
		/// theta_j = 2 pi j / N increases, and the ray direction is
		/// ( cos theta, sin theta ), so the traversal is counter-clockwise in
		/// ( R, z ) about the axis and both integrals come out positive. A
		/// negative answer is not a sign convention to absorb: it says the fit's
		/// nodes are not going round the axis once, which is the star-shapedness
		/// hypothesis having failed without the transversality test noticing.
		void enclosed( AngleParametrisation const &fit, bool differenced,
		               double &volume, double &area )
		{
			std::size_t const n = fit.count();
			double const step = twoPi/static_cast<double>( n );

			volume = 0.0;
			area = 0.0;
			for ( std::size_t j = 0; j < n; ++j )
			{
				double const theta = step*static_cast<double>( j );
				double const cosine = std::cos( theta );
				double const sine = std::sin( theta );

				double slope = fit.radiusPrime[ j ];
				if ( differenced )
				{
					// THE CONTROL, and it is the metric trap on a third
					// integrand: rho' by central difference of the neighbouring
					// rho_j instead of pointwise from the solved flux. Second
					// order, and nothing about the periodic trapezoid it feeds
					// says so.
					std::size_t const next = ( j + 1 )%n;
					std::size_t const previous = ( j + n - 1 )%n;
					slope = ( fit.radius[ next ] - fit.radius[ previous ] )
					        /( 2.0*step );
				}

				double const dzdtheta = slope*sine + fit.radius[ j ]*cosine;
				double const rNode = fit.pointR[ j ];

				volume += pi*rNode*rNode*dzdtheta*step;
				area += rNode*dzdtheta*step;
			}
		}
	}

	FluxSurfaceFamily extractFluxSurfaces( ContourTracer const &tracer,
	                                       CriticalPoint const &axis,
	                                       double psiBoundary,
	                                       FluxFamilyOptions const &options )
	{
		if ( options.surfaces < 2 )
			throw std::invalid_argument(
				"meq::extractFluxSurfaces: a family needs at least two surfaces "
				"to be interpolated between, and " + std::to_string( options.surfaces )
				+ " were asked for" );

		if ( options.angles < 3 )
			throw std::invalid_argument(
				"meq::extractFluxSurfaces: a surface needs at least three angles "
				"to enclose anything, and " + std::to_string( options.angles )
				+ " were asked for" );

		if ( !( options.innerCut > 0.0 ) || !( options.outerCut < 1.0 )
		     || !( options.innerCut < options.outerCut ) )
			throw std::invalid_argument(
				"meq::extractFluxSurfaces: the cut [ " + number( options.innerCut )
				+ ", " + number( options.outerCut ) + " ] is not a range of "
				"normalised flux strictly inside ( 0, 1 ). Psi_N = 0 is the "
				"magnetic axis, where a surface is a point and drho/dpsi is "
				"unbounded, and Psi_N = 1 is the plasma boundary, where "
				"1/| grad psi | diverges if it is a separatrix. See "
				"FluxFamily.hpp on why both ends are cut and why the two numbers "
				"are different" );

		double const span = axis.psi - psiBoundary;
		if ( span == 0.0 )
			throw std::runtime_error(
				"meq::extractFluxSurfaces: psi at the magnetic axis equals psi on "
				"the boundary, so there is no flux span to normalise by and no "
				"flux label to express a geometry against" );

		FluxSurfaceFamily family;
		family.axisR = axis.r;
		family.axisZ = axis.z;
		family.psiAxis = axis.psi;
		family.psiBoundary = psiBoundary;
		family.innerCut = options.innerCut;
		family.outerCut = options.outerCut;
		family.angles = options.angles;
		family.safetyFactorAvailable = static_cast<bool>( options.toroidalField );
		family.surfaces.reserve( options.surfaces );

		double const innerLabel = std::sqrt( options.innerCut );
		double const outerLabel = std::sqrt( options.outerCut );
		double const last = static_cast<double>( options.surfaces - 1 );

		for ( std::size_t i = 0; i < options.surfaces; ++i )
		{
			double const t = static_cast<double>( i )/last;

			double label = 0.0;
			if ( options.spacing == FluxLevelSpacing::Radial )
				label = innerLabel + t*( outerLabel - innerLabel );
			else
				label = std::sqrt( options.innerCut
				                   + t*( options.outerCut - options.innerCut ) );

			double const normalised = label*label;
			double const level = fluxAtNormalised( normalised, axis.psi,
			                                       psiBoundary );

			Contour contour;
			AngleParametrisation fit;
			try
			{
				contour = tracer.traceFromAxis( level, axis );
				if ( !contour.closed() )
					throw std::runtime_error(
						"the contour did not close, so there is no closed surface "
						"to average over" );
				fit = tracer.fitByAngle( contour, axis, options.angles );
			}
			catch ( std::exception const &error )
			{
				// ONE BAD LEVEL TAKES THE FAMILY WITH IT, and the message says
				// which so that the cut can be moved rather than guessed at. A
				// family with a hole would be bridged silently by
				// FluxSurfaceFamily::at()'s interpolation -- see the header.
				throw std::runtime_error(
					std::string( "meq::extractFluxSurfaces: the surface at "
					             "Psi_N = " ) + number( normalised )
					+ " ( rho = " + number( label ) + ", psi = " + number( level )
					+ " ) could not be extracted: " + error.what()
					+ ". The family is abandoned rather than returned with a hole "
					"in it, since an interpolation would bridge the hole "
					"silently. Move the cut, or look at the mesh there" );
			}

			SurfaceAverages const averages = surfaceAverages( tracer, fit );

			FluxSurface surface;
			surface.level = level;
			surface.normalisedFlux = normalised;
			surface.radial = label;
			surface.r = fit.pointR;
			surface.z = fit.pointZ;
			surface.extended = fit.extended;

			surface.vPrime = averages.vPrime;
			surface.arcLength = averages.arcLength();
			surface.inverseRSquared = averages.inverseRSquared();
			surface.gradPsiSquaredOverRSquared =
				averages.gradPsiSquaredOverRSquared();

			// closed-integral 2 pi R dl, the toroidal surface area. Every weight
			// is 2 pi R | dx/ds | ds / | grad psi |, so an integrand of
			// | grad psi | recovers 2 pi R dl exactly -- the same device
			// SurfaceAverages::arcLength() uses, and the same caveat applies:
			// it is the weights being read back, not an independent measurement
			// of them.
			surface.surfaceArea = averages.integrate(
				[]( SurfaceNode const &node ) { return node.gradient; } );

			surface.absGradPsi = averages.average(
				[]( SurfaceNode const &node ) { return node.gradient; } );
			surface.gradPsiSquared = averages.average(
				[]( SurfaceNode const &node )
				{ return node.gradient*node.gradient; } );

			if ( options.toroidalField )
				surface.safetyFactor =
					averages.safetyFactor( options.toroidalField( level ) );

			double volume = 0.0;
			double area = 0.0;
			enclosed( fit, false, volume, area );
			if ( !( volume > 0.0 ) )
				throw std::runtime_error(
					"meq::extractFluxSurfaces: the surface at Psi_N = "
					+ number( normalised ) + " encloses a volume of "
					+ number( volume ) + ", which is not positive. Green's "
					"theorem on a fit whose nodes go once round the axis "
					"counter-clockwise cannot do that, so the fit is not going "
					"round the axis -- star-shapedness has failed in a way the "
					"transversality test ( " + number( fit.transversality )
					+ " ) did not catch" );
			surface.volume = volume;
			surface.crossSectionArea = area;

			surface.extendedNodes = fit.extendedNodes;
			surface.crossesBand = fit.crossesBand();
			surface.deepestBandNode = fit.deepestBandNode;
			surface.worstResidual = fit.worstResidual;
			surface.transversality = fit.transversality;
			surface.transverse = fit.transverse;
			surface.stalledRays = fit.stalledRays;
			surface.fallbackLocations = fit.fallbackLocations;

			family.surfaces.push_back( std::move( surface ) );
		}

		return family;
	}

	FluxSurfaceFamily extractFluxSurfaces( GradShafranovSolver const &solver,
	                                       FluxFamilyOptions const &options,
	                                       Potential which )
	{
		ContourTracer const tracer( solver, which );
		CriticalPointFinder const finder( solver );
		CriticalPoint const axis = finder.findAxis();
		return extractFluxSurfaces( tracer, axis, solver.psiBoundary(), options );
	}

	double enclosedVolume( AngleParametrisation const &fit )
	{
		if ( fit.count() < 3 )
			throw std::invalid_argument(
				"meq::enclosedVolume: the fit carries "
				+ std::to_string( fit.count() )
				+ " nodes, which encloses nothing" );

		double volume = 0.0;
		double area = 0.0;
		enclosed( fit, false, volume, area );
		return volume;
	}

	double enclosedVolumeByDifferencedSlope( AngleParametrisation const &fit )
	{
		if ( fit.count() < 3 )
			throw std::invalid_argument(
				"meq::enclosedVolumeByDifferencedSlope: the fit carries "
				+ std::to_string( fit.count() )
				+ " nodes, which encloses nothing" );

		double volume = 0.0;
		double area = 0.0;
		enclosed( fit, true, volume, area );
		return volume;
	}

}
