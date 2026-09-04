#include "SurfaceAverage.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "mfem.hpp"

/*
 * INVERSION-PLAN.md stage IN-2. The header carries the conventions, the design
 * argument and every trap; this file is the arithmetic, and the comments in it
 * are only for the places where the code makes a choice the header does not
 * already explain -- where the band flag is read from, why the differenced
 * control needs a second loop, and why the contour builder marks a whole
 * segment.
 */

namespace
{

	double const twoPi = 6.283185307179586476925286766559;

}

namespace meq
{

	double SurfaceAverages::integrate( SurfaceIntegrand const &f ) const
	{
		if ( nodes.empty() )
			throw std::runtime_error(
				"SurfaceAverages::integrate: the surface has no nodes" );

		double sum = 0.0;
		for ( SurfaceNode const &node : nodes )
			sum += node.weight*f( node );
		return sum;
	}

	double SurfaceAverages::average( SurfaceIntegrand const &f ) const
	{
		if ( !( vPrime > 0.0 ) )
			throw std::runtime_error(
				"SurfaceAverages::average: V' is not positive, so there is "
				"nothing to divide by. Either the surface has no nodes or every "
				"weight came out zero, which on a real surface it cannot" );

		return integrate( f )/vPrime;
	}

	double SurfaceAverages::integrateDifferenced( SurfaceIntegrand const &f ) const
	{
		if ( !differencedAvailable() )
			throw std::runtime_error(
				"SurfaceAverages::integrateDifferenced: this surface carries no "
				"differenced control. It was built from a traced contour, where "
				"the nodes are curvature-controlled rather than equispaced and "
				"there is no neighbour to difference against. Returning the "
				"pointwise answer instead would make the control agree with the "
				"answer for the worst possible reason" );

		double sum = 0.0;
		for ( SurfaceNode const &node : nodes )
			sum += node.differencedWeight*f( node );
		return sum;
	}

	double SurfaceAverages::averageDifferenced( SurfaceIntegrand const &f ) const
	{
		// integrateDifferenced() is called FIRST so that a surface carrying no
		// differenced control gets that message rather than a complaint about a
		// zero denominator, which is a symptom of it and reads like a different
		// fault.
		double const weighted = integrateDifferenced( f );

		if ( !( vPrimeDifferenced > 0.0 ) )
			throw std::runtime_error(
				"SurfaceAverages::averageDifferenced: the differenced V' is not "
				"positive" );

		return weighted/vPrimeDifferenced;
	}

	double SurfaceAverages::inverseRSquared() const
	{
		return average( []( SurfaceNode const &node )
		{
			return 1.0/( node.r*node.r );
		} );
	}

	double SurfaceAverages::gradPsiSquaredOverRSquared() const
	{
		return average( []( SurfaceNode const &node )
		{
			return node.gradient*node.gradient/( node.r*node.r );
		} );
	}

	double SurfaceAverages::arcLength() const
	{
		// dl = | grad psi | / ( 2 pi R ) times the weight, by construction of
		// the weight. Written this way rather than as a second loop over the
		// metric so that it goes through the same facility every other quantity
		// does -- which is the whole point of there being a facility.
		return integrate( []( SurfaceNode const &node )
		{
			return node.gradient/( twoPi*node.r );
		} );
	}

	double SurfaceAverages::safetyFactor( double g ) const
	{
		return vPrime*g*inverseRSquared()/( twoPi*twoPi );
	}

	SurfaceAverages surfaceAverages( ContourTracer const &tracer,
	                                 AngleParametrisation const &fit )
	{
		std::size_t const count = fit.count();
		if ( count < 4 )
			throw std::invalid_argument(
				"meq::surfaceAverages: the angle fit has fewer than four nodes" );

		SurfaceAverages surface;
		surface.level = fit.level;
		surface.transversality = fit.transversality;
		surface.transverse = fit.transverse;
		surface.differencedControl = true;
		surface.nodes.resize( count );

		double const dTheta = twoPi/static_cast<double>( count );

		for ( std::size_t j = 0; j < count; ++j )
		{
			SurfaceNode &node = surface.nodes[ j ];
			node.parameter = dTheta*static_cast<double>( j );
			node.r = fit.pointR[ j ];
			node.z = fit.pointZ[ j ];
			node.metric = fit.speed[ j ];

			// EVERYTHING THIS LOOP NEEDS IS ALREADY IN THE FIT, SO IT SAMPLES
			// NOTHING. The ray Newton that placed each node evaluated the field
			// there and AngleParametrisation keeps what it found -- the
			// potential, the flux and the band flag alike. Re-deriving them here
			// would cost one location per node in the one place CLAUDE.md
			// records Mesh::FindPoints as O( elements x points ), and on a
			// discontinuous field it could land on the other side of a face from
			// the fit and return a different number. Reading them from the fit
			// is cheaper AND keeps the two in agreement by construction.
			node.psi = fit.potential[ j ];
			node.qR = fit.fluxR[ j ];
			node.qZ = fit.fluxZ[ j ];
			node.extended = j < fit.extended.size() && fit.extended[ j ] != 0;

			node.gradient = node.r*std::sqrt( node.qR*node.qR + node.qZ*node.qZ );
			if ( !( node.gradient > 0.0 ) )
			{
				std::ostringstream message;
				message << "meq::surfaceAverages: | grad psi | vanishes at ( "
				        << node.r << ", " << node.z << " ), so the weight "
				        << "dl / | grad psi | is not defined. The level is at or "
				        << "beyond a critical point";
				throw std::runtime_error( message.str() );
			}

			node.residual = std::abs( node.psi - surface.level );
			node.weight = twoPi*node.r*node.metric*dTheta/node.gradient;

			surface.vPrime += node.weight;
			surface.worstResidual = std::max( surface.worstResidual,
			                                  node.residual );
			if ( node.extended )
				++surface.extendedNodes;
		}

		// The differenced control, in its own loop because it needs both
		// neighbours and therefore every position first. It differences the
		// POSITIONS rather than the radii, which is what
		// meq::analytic::SurfaceQuadrature::differencedMetric() does, so that
		// the control on the exact field and the control on the discrete one are
		// the same control.
		for ( std::size_t j = 0; j < count; ++j )
		{
			std::size_t const next = ( j + 1 )%count;
			std::size_t const previous = ( j + count - 1 )%count;

			double const dr = ( surface.nodes[ next ].r
			                    - surface.nodes[ previous ].r )/( 2.0*dTheta );
			double const dz = ( surface.nodes[ next ].z
			                    - surface.nodes[ previous ].z )/( 2.0*dTheta );

			SurfaceNode &node = surface.nodes[ j ];
			node.differencedMetric = std::sqrt( dr*dr + dz*dz );
			node.differencedWeight = twoPi*node.r*node.differencedMetric*dTheta
			                         /node.gradient;

			surface.vPrimeDifferenced += node.differencedWeight;
		}

		surface.extended = surface.extendedNodes > 0;
		surface.deepestBandNode = fit.deepestBandNode;
		surface.bandExtension = tracer.bandExtension();

		// THE FIT'S COUNT, NOT THIS QUADRATURE'S, and the difference is worth
		// stating rather than glossing. ContourTracer::sampleAt() does not
		// report whether its element walk fell back on Mesh::FindPoints -- the
		// private seam counts them and the public one does not expose the
		// counter -- so what is carried here is the count from the ray Newton
		// that placed these nodes. It is the right order of magnitude and the
		// right diagnosis; it is not this loop's own tally.
		surface.fallbackLocations = fit.fallbackLocations;
		return surface;
	}

	SurfaceAverages surfaceAverages( ContourTracer const &tracer,
	                                 Contour const &contour, int gaussPoints )
	{
		if ( !contour.closed() )
			throw std::invalid_argument(
				"meq::surfaceAverages: the contour is not closed, and every rule "
				"here is a closed one. An open surface wants IN-5" );
		if ( contour.segments() == 0 )
			throw std::invalid_argument(
				"meq::surfaceAverages: the contour has no segments" );
		if ( gaussPoints < 1 )
			throw std::invalid_argument(
				"meq::surfaceAverages: need a quadrature" );

		// The same rule Contour::hermiteLength() uses, for the same reason: the
		// integrand is a smooth function of the segment parameter within a
		// segment, and a Gauss rule of n points is exact to degree 2n-1 on it.
		// Across segments the Hermite is only C^0 in general, so this is a
		// COMPOSITE rule and not one global one -- which is right, and is why
		// the accuracy comes from the segment count rather than from n.
		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*gaussPoints - 1 );

		SurfaceAverages surface;
		surface.level = contour.level;
		surface.differencedControl = false;
		surface.transverse = true;

		int hint = contour.points.front().element;

		for ( std::size_t i = 0; i < contour.segments(); ++i )
		{
			for ( int j = 0; j < rule.GetNPoints(); ++j )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( j );

				SurfaceNode node;
				node.parameter = static_cast<double>( i ) + ip.x;

				contour.pointOnSegment( i, ip.x, node.r, node.z );

				double dr = 0.0;
				double dz = 0.0;
				contour.tangentOnSegment( i, ip.x, dr, dz );
				node.metric = std::sqrt( dr*dr + dz*dz );

				// PER POINT, NOT PER SEGMENT. sampleAt()'s seven-argument
				// overload exists for exactly this caller: a Gauss node lies
				// strictly between two accepted points, so the band does not
				// respect the segment it sits in and either endpoint can be
				// inside Omega_h while the node is not. Marking the whole
				// segment from its endpoints -- which this loop used to do --
				// over-reports at the ends of a band excursion and, worse,
				// UNDER-reports for a node whose own segment straddles Gamma_h
				// with both endpoints inside. An under-report is a band
				// quantity presented as a solved one.
				bool nodeExtended = false;

				if ( !tracer.sampleAt( node.r, node.z, node.psi, node.qR,
				                       node.qZ, hint, nodeExtended ) )
				{
					std::ostringstream message;
					message << "meq::surfaceAverages: the Hermite quadrature "
					        << "point at ( " << node.r << ", " << node.z
					        << " ) on segment " << i << " is not in the mesh. "
					        << "The interpolant has left the domain between two "
					        << "accepted points, which on the fitted path means "
					        << "the contour is not interior";
					throw std::runtime_error( message.str() );
				}

				node.extended = nodeExtended;

				node.gradient = node.r*std::sqrt( node.qR*node.qR
				                                  + node.qZ*node.qZ );
				if ( !( node.gradient > 0.0 ) )
					throw std::runtime_error(
						"meq::surfaceAverages: | grad psi | vanishes on the "
						"Hermite interpolant, so dl / | grad psi | is not "
						"defined there" );

				node.residual = std::abs( node.psi - surface.level );
				node.weight = twoPi*node.r*node.metric*ip.weight/node.gradient;

				surface.vPrime += node.weight;
				surface.worstResidual = std::max( surface.worstResidual,
				                                  node.residual );
				if ( node.extended )
					++surface.extendedNodes;

				surface.nodes.push_back( node );
			}
		}

		surface.extended = surface.extendedNodes > 0;
		surface.deepestBandNode = contour.deepestBandPoint;
		surface.bandExtension = contour.bandExtension;
		surface.fallbackLocations = contour.fallbackLocations;
		return surface;
	}

	AveragedEquationResidual averagedGradShafranovResidual(
		ContourTracer const &tracer, CriticalPoint const &axis, double level,
		std::size_t angles, double step,
		std::function<double( double r, double z, double psi )> const &f,
		FluxDerivative derivative )
	{
		if ( !( step > 0.0 ) )
			throw std::invalid_argument(
				"meq::averagedGradShafranovResidual: the flux step must be "
				"positive" );

		AveragedEquationResidual out;

		// V' < | grad psi |^2 / R^2 >, which is the integral rather than the
		// average times V' -- the same number, one fewer division, and it is the
		// quantity the identity differentiates.
		auto weighted = [ & ]( double c, double &worst, bool &band )
		{
			Contour const contour = tracer.traceFromAxis( c, axis );
			SurfaceAverages const surface =
				surfaceAverages( tracer, tracer.fitByAngle( contour, axis, angles ) );

			worst = std::max( worst, surface.worstResidual );
			band = band || surface.extended;

			return surface.integrate( []( SurfaceNode const &node )
			{
				return node.gradient*node.gradient/( node.r*node.r );
			} );
		};

		// Sequenced by hand. The lambda accumulates into out, and the order in
		// which the operands of a subtraction are evaluated is unspecified; the
		// accumulation is a max and an or so it would not matter, but a reader
		// should not have to establish that.
		double const above = weighted( level + step, out.worstResidual,
		                               out.extended );
		double const below = weighted( level - step, out.worstResidual,
		                               out.extended );
		double const coarse = ( above - below )/( 2.0*step );

		if ( derivative == FluxDerivative::Richardson )
		{
			double const halfAbove = weighted( level + 0.5*step,
			                                   out.worstResidual, out.extended );
			double const halfBelow = weighted( level - 0.5*step,
			                                   out.worstResidual, out.extended );
			out.derivative = ( 4.0*( halfAbove - halfBelow )/step - coarse )/3.0;
		}
		else
		{
			out.derivative = coarse;
		}

		Contour const here = tracer.traceFromAxis( level, axis );
		SurfaceAverages const surface =
			surfaceAverages( tracer, tracer.fitByAngle( here, axis, angles ) );
		out.worstResidual = std::max( out.worstResidual, surface.worstResidual );
		out.extended = out.extended || surface.extended;

		out.vPrime = surface.vPrime;
		out.leftHandSide = out.derivative/surface.vPrime;

		// Delta* psi = -F, so < Delta* psi / R^2 > = -< F / R^2 >. F is
		// evaluated at the level, which is what psi is on this surface -- the
		// caller's F is free to depend on psi and several of meq's do.
		out.rightHandSide = surface.average( [ & ]( SurfaceNode const &node )
		{
			return -f( node.r, node.z, level )/( node.r*node.r );
		} );

		out.residual = out.leftHandSide - out.rightHandSide;
		return out;
	}

}
