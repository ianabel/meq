#include "FluxFamily.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "Profiles.hpp"

namespace meq
{

	namespace
	{
		/// Fritsch-Carlson limited slopes, the standard shape-preserving choice.
		///
		/// The published rule, and the reason it is here rather than an ordinary
		/// cubic, is in FluxFamily.hpp's section 3: an overshoot in V' near the
		/// axis is a NEGATIVE volume derivative, which is a nonsense value
		/// arriving through the one route MANTA-COUPLING.md section 8 forbids.
		/// The interior slope is the weighted harmonic mean of the neighbouring
		/// secants, zeroed at an extremum; the ends take the one-sided
		/// three-point slope, limited so they cannot introduce one.
		std::vector<double> monotoneSlopes( std::vector<double> const &x,
		                                    std::vector<double> const &y )
		{
			std::size_t const n = x.size();
			std::vector<double> slope( n, 0.0 );
			if ( n < 2 )
				return slope;

			std::vector<double> h( n - 1, 0.0 );
			std::vector<double> secant( n - 1, 0.0 );
			for ( std::size_t i = 0; i + 1 < n; ++i )
			{
				h[ i ] = x[ i + 1 ] - x[ i ];
				secant[ i ] = ( y[ i + 1 ] - y[ i ] )/h[ i ];
			}

			if ( n == 2 )
			{
				slope[ 0 ] = secant[ 0 ];
				slope[ 1 ] = secant[ 0 ];
				return slope;
			}

			for ( std::size_t i = 1; i + 1 < n; ++i )
			{
				double const left = secant[ i - 1 ];
				double const right = secant[ i ];
				if ( left*right <= 0.0 )
				{
					slope[ i ] = 0.0;
					continue;
				}

				double const w1 = 2.0*h[ i ] + h[ i - 1 ];
				double const w2 = h[ i ] + 2.0*h[ i - 1 ];
				slope[ i ] = ( w1 + w2 )/( w1/left + w2/right );
			}

			// The ends: the three-point one-sided slope, then limited so that a
			// monotone data set stays monotone across the first and last
			// intervals. Without the limiter the end slope can turn the first
			// cubic round on itself, which is exactly the overshoot the rule
			// exists to prevent and is worst where the data is steepest -- which
			// here is the axis end.
			auto endSlope = []( double hFirst, double hSecond, double sFirst,
			                    double sSecond )
			{
				double d = ( ( 2.0*hFirst + hSecond )*sFirst - hFirst*sSecond )
				           /( hFirst + hSecond );
				if ( d*sFirst <= 0.0 )
					return 0.0;
				if ( sFirst*sSecond <= 0.0 && std::abs( d ) > std::abs( 3.0*sFirst ) )
					return 3.0*sFirst;
				return d;
			};

			slope[ 0 ] = endSlope( h[ 0 ], h[ 1 ], secant[ 0 ], secant[ 1 ] );
			slope[ n - 1 ] = endSlope( h[ n - 2 ], h[ n - 3 ], secant[ n - 2 ],
			                           secant[ n - 3 ] );
			return slope;
		}

		/// Evaluate the monotone cubic through ( x, y ) at @a at. @a x must be
		/// strictly increasing and @a at inside its range; both are the caller's
		/// to have checked.
		double monotoneCubic( std::vector<double> const &x,
		                      std::vector<double> const &y, double at )
		{
			std::size_t const n = x.size();
			std::vector<double> const slope = monotoneSlopes( x, y );

			std::size_t interval =
				static_cast<std::size_t>(
					std::upper_bound( x.begin(), x.end(), at ) - x.begin() );
			if ( interval == 0 )
				interval = 1;
			if ( interval > n - 1 )
				interval = n - 1;
			std::size_t const lower = interval - 1;

			// The evaluation itself is meq::HermiteCubicSpline rather than five
			// lines of Horner beside it: the standing preference is to take the
			// maintained implementation, and what this file adds is the slope
			// rule above and nothing else.
			HermiteCubicSpline const cubic( x[ lower ], x[ lower + 1 ],
			                                y[ lower ], y[ lower + 1 ],
			                                slope[ lower ], slope[ lower + 1 ] );
			return cubic( at );
		}

		std::string number( double value )
		{
			std::ostringstream out;
			out.precision( 6 );
			out << value;
			return out.str();
		}
	}

	double normalisedFlux( double psi, double psiAxis, double psiBoundary )
	{
		double const span = psiAxis - psiBoundary;
		if ( span == 0.0 )
			throw std::invalid_argument(
				"meq::normalisedFlux: psi_ax equals psi_bnd, so there is no span "
				"to normalise by and no normalised flux to return" );
		return ( psiAxis - psi )/span;
	}

	double fluxAtNormalised( double normalisedFluxIn, double psiAxis,
	                         double psiBoundary )
	{
		double const span = psiAxis - psiBoundary;
		if ( span == 0.0 )
			throw std::invalid_argument(
				"meq::fluxAtNormalised: psi_ax equals psi_bnd, so there is no "
				"span to normalise by" );
		return psiAxis - normalisedFluxIn*span;
	}

	double FluxSurfaceFamily::innerLabel() const
	{
		if ( surfaces.empty() )
			throw std::runtime_error(
				"meq::FluxSurfaceFamily::innerLabel: the family is empty, so it "
				"covers no range of the flux label" );
		return surfaces.front().radial;
	}

	double FluxSurfaceFamily::outerLabel() const
	{
		if ( surfaces.empty() )
			throw std::runtime_error(
				"meq::FluxSurfaceFamily::outerLabel: the family is empty, so it "
				"covers no range of the flux label" );
		return surfaces.back().radial;
	}

	double FluxSurfaceFamily::normalisedFluxDerivative() const
	{
		double const span = psiAxis - psiBoundary;
		if ( span == 0.0 )
			throw std::runtime_error(
				"meq::FluxSurfaceFamily::normalisedFluxDerivative: psi_ax equals "
				"psi_bnd, so the family carries no normalisation" );
		return -1.0/span;
	}

	double FluxSurfaceFamily::radialDerivative( double label ) const
	{
		if ( label <= 0.0 )
			throw std::runtime_error(
				"meq::FluxSurfaceFamily::radialDerivative: drho/dpsi is unbounded "
				"at rho = 0, which is the magnetic axis. That is the coordinate "
				"and not a defect -- see FluxFamily.hpp -- and it is one of the "
				"reasons the family is cut away from the axis" );
		return normalisedFluxDerivative()/( 2.0*label );
	}

	bool FluxSurfaceFamily::covers( double label ) const
	{
		if ( surfaces.empty() )
			return false;
		return label >= innerLabel() && label <= outerLabel();
	}

	double FluxSurfaceFamily::at( double label, SurfaceQuantity const &of ) const
	{
		if ( !of )
			throw std::invalid_argument(
				"meq::FluxSurfaceFamily::at: the quantity is empty, so there is "
				"nothing to interpolate" );

		if ( surfaces.size() < 2 )
			throw std::runtime_error(
				"meq::FluxSurfaceFamily::at: the family has "
				+ std::to_string( surfaces.size() )
				+ " surface(s), so there is nothing to interpolate between" );

		if ( !covers( label ) )
			throw std::domain_error(
				"meq::FluxSurfaceFamily::at: the flux label " + number( label )
				+ " is outside the family's range [ " + number( innerLabel() )
				+ ", " + number( outerLabel() ) + " ], which is Psi_N in [ "
				+ number( innerCut ) + ", " + number( outerCut )
				+ " ]. THIS IS A REFUSAL AND NOT A FAILURE TO EXTRAPOLATE: "
				"everything degrades approaching the separatrix and the surfaces "
				"shrink to a point at the axis, so meq cuts both ends and says "
				"so. Move the query inside the cut, or widen it and pay for the "
				"surfaces" );

		std::vector<double> x( surfaces.size(), 0.0 );
		std::vector<double> y( surfaces.size(), 0.0 );
		for ( std::size_t i = 0; i < surfaces.size(); ++i )
		{
			x[ i ] = surfaces[ i ].radial;
			y[ i ] = of( surfaces[ i ] );

			// STRICTLY INCREASING, CHECKED. meq::extractFluxSurfaces() builds
			// them so, and a family assembled by hand out of order would
			// otherwise divide by a zero interval width and return a NaN --
			// which is the plausible-looking nonsense MANTA-COUPLING.md
			// section 8 forbids, arriving through a caller's mistake rather
			// than through the physics.
			if ( i > 0 && !( x[ i ] > x[ i - 1 ] ) )
				throw std::runtime_error(
					"meq::FluxSurfaceFamily::at: surface " + std::to_string( i )
					+ " is at flux label " + number( x[ i ] )
					+ " and surface " + std::to_string( i - 1 ) + " at "
					+ number( x[ i - 1 ] )
					+ ", so the family is not ordered by strictly increasing "
					"label and there is nothing to interpolate along" );
		}

		return monotoneCubic( x, y, label );
	}

	double FluxSurfaceFamily::atNormalisedFlux( double normalisedFluxIn,
	                                            SurfaceQuantity const &of ) const
	{
		if ( normalisedFluxIn < 0.0 )
			throw std::domain_error(
				"meq::FluxSurfaceFamily::atNormalisedFlux: the normalised flux "
				+ number( normalisedFluxIn )
				+ " is negative, so it has no square root and names no surface" );
		return at( std::sqrt( normalisedFluxIn ), of );
	}

	double FluxSurfaceFamily::vPrimeAt( double label ) const
	{
		return at( label, []( FluxSurface const &s ) { return s.vPrime; } );
	}

	double FluxSurfaceFamily::volumeAt( double label ) const
	{
		return at( label, []( FluxSurface const &s ) { return s.volume; } );
	}

	double FluxSurfaceFamily::inverseRSquaredAt( double label ) const
	{
		return at( label, []( FluxSurface const &s )
			{ return s.inverseRSquared; } );
	}

	double FluxSurfaceFamily::gradPsiSquaredOverRSquaredAt( double label ) const
	{
		return at( label, []( FluxSurface const &s )
			{ return s.gradPsiSquaredOverRSquared; } );
	}

	double FluxSurfaceFamily::arcLengthAt( double label ) const
	{
		return at( label, []( FluxSurface const &s ) { return s.arcLength; } );
	}

	double FluxSurfaceFamily::surfaceAreaAt( double label ) const
	{
		return at( label, []( FluxSurface const &s ) { return s.surfaceArea; } );
	}

	double FluxSurfaceFamily::absGradLabelAt( double label ) const
	{
		// | grad rho | = | drho/dpsi | | grad psi |, with drho/dpsi ANALYTIC and
		// constant on the surface, so it comes out of the average exactly. No
		// difference anywhere -- see FluxFamily.hpp on the metric trap.
		double const chain = std::abs( radialDerivative( label ) );
		return chain*at( label, []( FluxSurface const &s )
			{ return s.absGradPsi; } );
	}

	double FluxSurfaceFamily::gradLabelSquaredAt( double label ) const
	{
		double const chain = radialDerivative( label );
		return chain*chain*at( label, []( FluxSurface const &s )
			{ return s.gradPsiSquared; } );
	}

	bool FluxSurfaceFamily::crossesBand() const
	{
		for ( FluxSurface const &s : surfaces )
			if ( s.crossesBand )
				return true;
		return false;
	}

	int FluxSurfaceFamily::extendedNodes() const
	{
		int total = 0;
		for ( FluxSurface const &s : surfaces )
			total += s.extendedNodes;
		return total;
	}

	double FluxSurfaceFamily::worstResidual() const
	{
		double worst = 0.0;
		for ( FluxSurface const &s : surfaces )
			worst = std::max( worst, s.worstResidual );
		return worst;
	}

	GeometryCache::GeometryCache( Extractor extractorIn )
		: extractor( std::move( extractorIn ) )
	{
		if ( !extractor )
			throw std::invalid_argument(
				"meq::GeometryCache: the extractor is empty, so the cache could "
				"never fill itself -- every query would miss and then fail with "
				"something unrelated to the cause" );
	}

	FluxSurfaceFamily const &
	GeometryCache::family( std::vector<double> const &psi )
	{
		++queryCount;

		if ( psi.empty() )
			throw std::invalid_argument(
				"meq::GeometryCache::family: the psi vector is empty, so there is "
				"no state to extract a geometry from" );

		// THE NON-FINITE GUARD. MANTA-COUPLING.md section 1 guarantees meq is
		// called at states far from equilibrium as a normal part of the
		// consumer's Newton iteration, and a NaN or an infinity in psi must come
		// back as a refusal the integrator can retry from rather than as a
		// geometry. Cheap: one pass over a vector that is about to be compared
		// entry by entry anyway.
		for ( std::size_t i = 0; i < psi.size(); ++i )
		{
			if ( !std::isfinite( psi[ i ] ) )
				throw std::invalid_argument(
					"meq::GeometryCache::family: psi[ " + std::to_string( i )
					+ " ] is not finite, so there is no equilibrium here to "
					"extract surfaces from. Refusing rather than returning a "
					"geometry: see MANTA-COUPLING.md section 8" );
		}

		if ( valid && key.size() == psi.size()
		     && std::memcmp( key.data(), psi.data(),
		                     key.size()*sizeof( double ) ) == 0 )
		{
			++hitCount;
			return heldFamily;
		}

		// INVALIDATE BEFORE EXTRACTING, so that a throw commits nothing. Leaving
		// the previous family under the previous key would serve a different
		// state's geometry the moment the caller retried at the old psi; leaving
		// it under the new key would serve it immediately. See FluxFamily.hpp's
		// section 6.
		valid = false;
		key.clear();
		heldFamily = FluxSurfaceFamily{};

		++extractionCount;
		FluxSurfaceFamily built = extractor( psi );

		heldFamily = std::move( built );
		key = psi;
		valid = true;
		return heldFamily;
	}

	double GeometryCache::at( std::vector<double> const &psi, double label,
	                          SurfaceQuantity const &of )
	{
		return family( psi ).at( label, of );
	}

	void GeometryCache::resetForRun()
	{
		valid = false;
		key.clear();
		heldFamily = FluxSurfaceFamily{};
	}

	bool GeometryCache::held() const
	{
		return valid;
	}

	std::size_t GeometryCache::extractions() const
	{
		return extractionCount;
	}

	std::size_t GeometryCache::queries() const
	{
		return queryCount;
	}

	std::size_t GeometryCache::hits() const
	{
		return hitCount;
	}

}
