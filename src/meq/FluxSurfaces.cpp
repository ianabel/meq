#include "FluxSurfaces.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

/*
 * The implementation of INVERSION-PLAN.md stages IN-0 and IN-1.
 * FluxSurfaces.hpp carries the argument; this file carries the arithmetic, and
 * comments here are confined to the places where the arithmetic is not the
 * obvious transcription of it.
 */

namespace
{

	double const twoPi = 6.283185307179586476925286766559;

	/// The signed angular increment from @a from to @a to, folded into
	/// ( -pi, pi ]. Identical in intent to CriticalPoints.cpp's angleStep(): it
	/// is the branch choice that makes an accumulated turning a turning number,
	/// and it assumes the step resolves the rotation. The tracer's own
	/// safeguard -- halving a step whose turning exceeds a third of pi -- is
	/// what keeps that assumption true.
	double angleStep( double from, double to )
	{
		double delta = to - from;
		while ( delta > M_PI )
			delta -= twoPi;
		while ( delta <= -M_PI )
			delta += twoPi;
		return delta;
	}

	/// The cubic Hermite basis and its derivative, on [ 0, 1 ].
	void hermiteBasis( double t, double basis[ 4 ] )
	{
		double const t2 = t*t;
		double const t3 = t2*t;
		basis[ 0 ] =  2.0*t3 - 3.0*t2 + 1.0;
		basis[ 1 ] =        t3 - 2.0*t2 + t;
		basis[ 2 ] = -2.0*t3 + 3.0*t2;
		basis[ 3 ] =        t3 -     t2;
	}

	void hermiteBasisPrime( double t, double basis[ 4 ] )
	{
		double const t2 = t*t;
		basis[ 0 ] =  6.0*t2 - 6.0*t;
		basis[ 1 ] =  3.0*t2 - 4.0*t + 1.0;
		basis[ 2 ] = -6.0*t2 + 6.0*t;
		basis[ 3 ] =  3.0*t2 - 2.0*t;
	}

}

namespace meq
{

	char const *contourStatusName( ContourStatus status )
	{
		switch ( status )
		{
			case ContourStatus::Closed:
				return "closed";
			case ContourStatus::LeftMesh:
				return "left mesh";
			case ContourStatus::Stalled:
				return "stalled";
			case ContourStatus::TooLong:
				break;
		}
		return "too long";
	}

	char const *bandExtensionName( BandExtension which )
	{
		switch ( which )
		{
			case BandExtension::None:
				return "none";
			case BandExtension::FluxTaylor:
				return "flux Taylor";
			case BandExtension::TransferLift:
				break;
		}
		return "transfer lift";
	}

	std::size_t Contour::segments() const
	{
		return points.size() > 1 ? points.size() - 1 : 0;
	}

	double Contour::length() const
	{
		return points.empty() ? 0.0 : points.back().arcLength;
	}

	void Contour::pointOnSegment( std::size_t i, double t, double &r,
	                              double &z ) const
	{
		if ( i + 1 >= points.size() )
			throw std::out_of_range( "Contour::pointOnSegment: no such segment" );

		ContourPoint const &a = points[ i ];
		ContourPoint const &b = points[ i + 1 ];

		// The tangent magnitudes are the CHORD LENGTH. See the header: this is
		// not the optimal scaling and it is still fourth order.
		double const chord = std::sqrt( ( b.r - a.r )*( b.r - a.r )
		                                + ( b.z - a.z )*( b.z - a.z ) );

		double basis[ 4 ];
		hermiteBasis( t, basis );

		r = basis[ 0 ]*a.r + basis[ 1 ]*chord*a.tangentR
		    + basis[ 2 ]*b.r + basis[ 3 ]*chord*b.tangentR;
		z = basis[ 0 ]*a.z + basis[ 1 ]*chord*a.tangentZ
		    + basis[ 2 ]*b.z + basis[ 3 ]*chord*b.tangentZ;
	}

	void Contour::tangentOnSegment( std::size_t i, double t, double &r,
	                                double &z ) const
	{
		if ( i + 1 >= points.size() )
			throw std::out_of_range( "Contour::tangentOnSegment: no such segment" );

		ContourPoint const &a = points[ i ];
		ContourPoint const &b = points[ i + 1 ];

		double const chord = std::sqrt( ( b.r - a.r )*( b.r - a.r )
		                                + ( b.z - a.z )*( b.z - a.z ) );

		double basis[ 4 ];
		hermiteBasisPrime( t, basis );

		r = basis[ 0 ]*a.r + basis[ 1 ]*chord*a.tangentR
		    + basis[ 2 ]*b.r + basis[ 3 ]*chord*b.tangentR;
		z = basis[ 0 ]*a.z + basis[ 1 ]*chord*a.tangentZ
		    + basis[ 2 ]*b.z + basis[ 3 ]*chord*b.tangentZ;
	}

	void Contour::chordOnSegment( std::size_t i, double t, double &r,
	                              double &z ) const
	{
		if ( i + 1 >= points.size() )
			throw std::out_of_range( "Contour::chordOnSegment: no such segment" );

		ContourPoint const &a = points[ i ];
		ContourPoint const &b = points[ i + 1 ];

		r = ( 1.0 - t )*a.r + t*b.r;
		z = ( 1.0 - t )*a.z + t*b.z;
	}

	double Contour::hermiteLength( int gaussPoints ) const
	{
		if ( gaussPoints < 1 )
			throw std::invalid_argument( "Contour::hermiteLength: need a quadrature" );

		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*gaussPoints - 1 );

		double total = 0.0;
		for ( std::size_t i = 0; i + 1 < points.size(); ++i )
		{
			for ( int j = 0; j < rule.GetNPoints(); ++j )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( j );
				double dr = 0.0;
				double dz = 0.0;
				tangentOnSegment( i, ip.x, dr, dz );
				total += ip.weight*std::sqrt( dr*dr + dz*dz );
			}
		}
		return total;
	}

	void Contour::pointAtArcLength( double s, double &r, double &z ) const
	{
		if ( points.size() < 2 )
			throw std::out_of_range( "Contour::pointAtArcLength: no segments" );

		double const total = length();
		if ( total > 0.0 )
		{
			s = std::fmod( s, total );
			if ( s < 0.0 )
				s += total;
		}

		std::size_t lo = 0;
		std::size_t hi = points.size() - 1;
		while ( hi - lo > 1 )
		{
			std::size_t const mid = ( lo + hi )/2;
			if ( points[ mid ].arcLength <= s )
				lo = mid;
			else
				hi = mid;
		}

		double const span = points[ lo + 1 ].arcLength - points[ lo ].arcLength;
		double const t = span > 0.0 ? ( s - points[ lo ].arcLength )/span : 0.0;
		pointOnSegment( lo, t, r, z );
	}

	namespace
	{

		/// Which pair of fields Potential names. Free rather than a member so
		/// that the two constructors can delegate.
		mfem::GridFunction const &potentialFor( GradShafranovSolver const &solver,
		                                        Potential which )
		{
			return which == Potential::Raw ? solver.potential()
			                               : solver.postProcessedPotential();
		}

		mfem::GridFunction const &fluxFor( GradShafranovSolver const &solver,
		                                   Potential which )
		{
			return which == Potential::Raw ? solver.flux()
			                               : solver.postProcessedFlux();
		}

	}

	ContourTracer::ContourTracer( GradShafranovSolver const &solverIn,
	                              Potential which )
		: ContourTracer( potentialFor( solverIn, which ), fluxFor( solverIn, which ) )
	{
		// THE DEFAULT NEEDS A GUARD BECAUSE ITS FAILURE IS SILENT.
		// postProcessedPotential() returns the GridFunction whether or not
		// postProcess() has filled it, so a tracer built on it too early would
		// root a field that is identically zero -- and the corrector would
		// simply fail to find the level, which reads as "this level is not
		// attained" rather than as "the solver was not asked for psi*".
		if ( which == Potential::PostProcessed && !solverIn.isPostProcessed() )
			throw std::invalid_argument(
				"ContourTracer: Potential::PostProcessed needs "
				"GradShafranovSolver::postProcess() to have been called since the "
				"last solve. It is the default because rooting psi* traces a "
				"contour one order closer to the true one; Potential::Raw is the "
				"other door and needs nothing" );
	}

	ContourTracer::ContourTracer( mfem::GridFunction const &potentialIn,
	                              mfem::GridFunction const &fluxIn )
		: potentialField( potentialIn ),
		  fluxField( fluxIn ),
		  meshRef( *fluxIn.FESpace()->GetMesh() ),
		  neighbours( ( *fluxIn.FESpace()->GetMesh() ).ElementToElementTable() )
	{
		if ( fluxIn.FESpace()->GetVDim() != 2 )
			throw std::invalid_argument( "ContourTracer: the flux must have vdim 2" );

		if ( meshRef.Dimension() != 2 )
			throw std::invalid_argument( "ContourTracer: meq is two dimensional" );

		if ( potentialIn.FESpace()->GetMesh() != &meshRef )
			throw std::invalid_argument(
				"ContourTracer: the potential and the flux are on different meshes" );
	}

	void ContourTracer::setStep( double stepIn )
	{
		if ( stepIn < 0.0 || !std::isfinite( stepIn ) )
			throw std::invalid_argument( "ContourTracer: the step must be finite and >= 0" );
		step = stepIn;
	}

	void ContourTracer::setTargetTurn( double turnIn )
	{
		if ( turnIn < 0.0 || turnIn > M_PI )
			throw std::invalid_argument( "ContourTracer: the target turn must be in [ 0, pi ]" );
		targetTurn = turnIn;
	}

	void ContourTracer::setMinStepFraction( double fractionIn )
	{
		if ( !( fractionIn > 0.0 ) || fractionIn > 1.0 )
			throw std::invalid_argument(
				"ContourTracer: the step floor must be a fraction in ( 0, 1 ]" );
		minStepFraction = fractionIn;
	}

	void ContourTracer::setLocalStepCeiling( double fractionIn )
	{
		if ( !std::isfinite( fractionIn ) )
			throw std::invalid_argument( "ContourTracer: the step ceiling must be finite" );
		localStepCeiling = fractionIn;
	}

	void ContourTracer::setTolerance( double toleranceIn )
	{
		if ( !( toleranceIn > 0.0 ) )
			throw std::invalid_argument( "ContourTracer: the tolerance must be positive" );
		tolerance = toleranceIn;
	}

	void ContourTracer::setMaxCorrectorIterations( int maxIterationsIn )
	{
		if ( maxIterationsIn < 1 )
			throw std::invalid_argument( "ContourTracer: need at least one iteration" );
		maxCorrectorIterations = maxIterationsIn;
	}

	void ContourTracer::setMaxPoints( int maxPointsIn )
	{
		if ( maxPointsIn < 4 )
			throw std::invalid_argument( "ContourTracer: need at least four points" );
		maxPoints = maxPointsIn;
	}

	void ContourTracer::setCircuits( int circuitsIn )
	{
		if ( circuitsIn < 1 )
			throw std::invalid_argument( "ContourTracer: need at least one circuit" );
		circuits = circuitsIn;
	}

	void ContourTracer::setMeasureFaceJumps( bool measureIn )
	{
		measureFaceJumps = measureIn;
	}

	void ContourTracer::setWalkDepth( int depthIn )
	{
		if ( depthIn < 1 )
			throw std::invalid_argument( "ContourTracer: need at least one ring" );
		walkDepth = depthIn;
	}

	void ContourTracer::setBandExtension( BandExtension which,
	                                      mfem::Array<int> const &gammaHMarkerIn,
	                                      mfem::TransferPath const *pathIn,
	                                      mfem::PositionFunction g,
	                                      double reach, int lineOrderIn )
	{
		if ( which == BandExtension::None )
		{
			clearBandExtension();
			return;
		}

		if ( which == BandExtension::TransferLift && pathIn == nullptr )
			throw std::invalid_argument(
				"ContourTracer::setBandExtension: BandExtension::TransferLift is "
				"the extension technique's own construction and needs the same "
				"mfem::TransferPath the solve was given -- a different family is a "
				"different lifting" );

		if ( !( reach > 0.0 ) || !std::isfinite( reach ) )
			throw std::invalid_argument(
				"ContourTracer::setBandExtension: the reach must be finite and positive" );

		if ( gammaHMarkerIn.Size() < meshRef.bdr_attributes.Max() )
			throw std::invalid_argument(
				"ContourTracer::setBandExtension: the Gamma_h marker must be sized "
				"by the largest boundary attribute of the mesh" );

		bool any = false;
		for ( int i = 0; i < gammaHMarkerIn.Size(); ++i )
			any = any || gammaHMarkerIn[ i ] != 0;
		if ( !any )
			throw std::invalid_argument(
				"ContourTracer::setBandExtension: the Gamma_h marker selects no "
				"boundary attribute, so there is no band to extend into" );

		// The faces of Gamma_h, with everything the band search needs taken now
		// rather than per call: this is O( boundary faces ) once, against a
		// nearest-face scan per band point.
		std::vector<BoundaryFace> faces;
		faces.reserve( meshRef.GetNBE() );

		mfem::IntegrationPoint at;
		mfem::Vector x( 2 );
		mfem::Vector centroid( 2 );

		for ( int be = 0; be < meshRef.GetNBE(); ++be )
		{
			int const attribute = meshRef.GetBdrAttribute( be );
			if ( attribute < 1 || attribute > gammaHMarkerIn.Size() )
				continue;
			if ( !gammaHMarkerIn[ attribute - 1 ] )
				continue;

			mfem::FaceElementTransformations *ftr =
				meshRef.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;

			BoundaryFace face;
			face.boundaryElement = be;
			face.element = ftr->Elem1No;

			// THE FACE TRANSFORMATION'S OWN PARAMETRISATION, read off rather
			// than assumed. mfem::VertexConePath interpolates the directions of
			// the two vertices along xi, and it decides which vertex sits at
			// xi = 0 by transforming and comparing -- precisely because the
			// mesh's vertex order need not agree. Taking the endpoints the same
			// way is what makes the face parameter this class computes the one
			// Endpoint() expects.
			at.Set1w( 0.0, 1.0 );
			ftr->Transform( at, x );
			face.r0 = x( 0 );
			face.z0 = x( 1 );
			at.Set1w( 1.0, 1.0 );
			ftr->Transform( at, x );
			face.r1 = x( 0 );
			face.z1 = x( 1 );

			double const dr = face.r1 - face.r0;
			double const dz = face.z1 - face.z0;
			face.length = std::sqrt( dr*dr + dz*dz );
			if ( !( face.length > 0.0 ) )
				continue;

			// AFTER the face transformations are finished with. Mesh keeps ONE
			// ElementTransformation and hands out a pointer to it, so asking for
			// an element's transformation invalidates the face's Elem1.
			mfem::ElementTransformation *trans =
				meshRef.GetElementTransformation( face.element );
			mfem::IntegrationPoint const &centre =
				mfem::Geometries.GetCenter( trans->GetGeometryType() );
			trans->Transform( centre, centroid );

			double normalR = dz/face.length;
			double normalZ = -dr/face.length;
			double const outR = 0.5*( face.r0 + face.r1 ) - centroid( 0 );
			double const outZ = 0.5*( face.z0 + face.z1 ) - centroid( 1 );
			if ( normalR*outR + normalZ*outZ < 0.0 )
			{
				normalR = -normalR;
				normalZ = -normalZ;
			}
			face.normalR = normalR;
			face.normalZ = normalZ;

			faces.push_back( face );
		}

		if ( faces.empty() )
			throw std::runtime_error(
				"ContourTracer::setBandExtension: the marked attributes carry no "
				"boundary face, so the marker names an attribute this mesh does "
				"not have. Extending nothing would read as the band being absent, "
				"which is the failure this refuses rather than survives" );

		bandFaces.swap( faces );
		bandMethod = which;
		bandPath = pathIn;
		bandReach = reach;
		bandLineOrder = lineOrderIn;
		bandDatum = g ? std::move( g )
		              : mfem::PositionFunction(
		                    []( mfem::Vector const & ) { return 0.0; } );
	}

	void ContourTracer::clearBandExtension()
	{
		bandMethod = BandExtension::None;
		bandFaces.clear();
		bandPath = nullptr;
		bandDatum = mfem::PositionFunction();
		bandReach = 2.0;
		bandLineOrder = -1;
	}

	BandExtension ContourTracer::bandExtension() const
	{
		return bandMethod;
	}

	std::size_t ContourTracer::bandFaceCount() const
	{
		return bandFaces.size();
	}

	double ContourTracer::elementSize( int element ) const
	{
		mfem::ElementTransformation *trans = meshRef.GetElementTransformation( element );
		mfem::IntegrationPoint const &centre =
			mfem::Geometries.GetCenter( trans->GetGeometryType() );
		trans->SetIntPoint( &centre );

		// The square root of the transformation's own volume measure. For a
		// straight triangle with legs a and b that is sqrt( a b ), which is the
		// leg length on the meshes meq runs; the point is that it is a LENGTH,
		// local, and costs one Jacobian.
		double const weight = trans->Weight();
		return weight > 0.0 ? std::sqrt( weight ) : 0.0;
	}

	double ContourTracer::potentialScale() const
	{
		double worst = 0.0;
		for ( int i = 0; i < potentialField.Size(); ++i )
			worst = std::max( worst, std::abs( potentialField( i ) ) );
		return worst > 0.0 ? worst : 1.0;
	}

	bool ContourTracer::locate( double r, double z, int hint, int &element,
	                            mfem::IntegrationPoint &ip, int &fallbacks,
	                            bool allowFallback ) const
	{
		mfem::Vector point( 2 );
		point( 0 ) = r;
		point( 1 ) = z;

		int const elements = meshRef.GetNE();

		// The walk. Every candidate is tried by inverting that element's own
		// transformation, which for a straight-sided element is exact and costs
		// a 2x2 solve.
		auto tryElement = [ & ]( int candidate ) -> bool
		{
			if ( candidate < 0 || candidate >= elements )
				return false;
			mfem::ElementTransformation *trans =
				meshRef.GetElementTransformation( candidate );
			return trans->TransformBack( point, ip )
			       == mfem::InverseElementTransformation::Inside;
		};

		if ( hint >= 0 && tryElement( hint ) )
		{
			element = hint;
			return true;
		}

		if ( hint >= 0 && hint < elements )
		{
			// A widening search over face neighbours, breadth first, to
			// walkDepth rings. One ring is not enough for the same reason
			// CriticalPointFinder::axisSeeds() takes two: a step that clips a
			// vertex leaves the element diagonally, and the element on the far
			// side of a vertex is not a face neighbour of the one left behind.
			// The visited list is a plain vector searched linearly rather than a
			// mark array over the mesh, because the frontier is a few dozen
			// elements and an O( elements ) allocation per location would
			// reintroduce exactly the cost this walk exists to avoid.
			std::vector<int> visited;
			visited.reserve( 64 );
			visited.push_back( hint );

			std::size_t frontier = 0;
			for ( int depth = 0; depth < walkDepth; ++depth )
			{
				std::size_t const end = visited.size();
				if ( frontier >= end )
					break;

				for ( std::size_t i = frontier; i < end; ++i )
				{
					int const *row = neighbours.GetRow( visited[ i ] );
					int const rowSize = neighbours.RowSize( visited[ i ] );
					for ( int j = 0; j < rowSize; ++j )
					{
						int const candidate = row[ j ];
						if ( candidate < 0 )
							continue;
						if ( std::find( visited.begin(), visited.end(), candidate )
						     != visited.end() )
							continue;
						visited.push_back( candidate );
						if ( tryElement( candidate ) )
						{
							element = candidate;
							return true;
						}
					}
				}
				frontier = end;
			}
		}

		// The last resort. Counted, because it is O( elements ) and a trace that
		// needs it repeatedly is quadratic in the mesh. See the header.
		//
		// With a band extension configured the caller asks for the walk WITHOUT
		// it first, because most failed walks are then band points and paying an
		// O( elements ) scan for each would reintroduce exactly that cost. It is
		// still owed to a point the band turns down; see sampleField().
		if ( !allowFallback )
			return false;

		++fallbacks;

		mfem::DenseMatrix matrix( 2, 1 );
		matrix( 0, 0 ) = r;
		matrix( 1, 0 ) = z;

		mfem::Array<int> found;
		mfem::Array<mfem::IntegrationPoint> ips;
		meshRef.FindPoints( matrix, found, ips, false );

		if ( found.Size() < 1 || found[ 0 ] < 0 )
			return false;

		element = found[ 0 ];
		ip = ips[ 0 ];
		return true;
	}

	int ContourTracer::nearestBandFace( double r, double z, double &footR,
	                                    double &footZ, double &parameter,
	                                    double &depth ) const
	{
		// THE NEAREST FACE THE POINT IS OUTSIDE OF, WHICH IS NOT THE NEAREST
		// FACE, AND THE DIFFERENCE COST A TRACE.
		//
		// The obvious version takes the nearest face and then checks the point
		// is on its outward side, refusing if not -- the check being there to
		// stop a point the element walk merely LOST from being answered by an
		// extension when it is sitting inside the mesh with a perfectly good
		// element of its own. The walk failing is not evidence that a point is
		// outside; the geometry is.
		//
		// That version truncates real contours, and Gamma_h being non-convex is
		// why. A staircase cut from a diagonally split Cartesian mesh PINCHES:
		// two triangles of D_h meeting at a single vertex are face-disconnected
		// there, so the exterior has a pinch point and the two lobes' faces are
		// EQUIDISTANT from a point just outside it. Which of them comes out
		// "nearest" is then a tie broken by loop order, and half the time it is
		// the lobe whose outward normal points the other way -- so a point that
		// really is in the band is refused, and the trace ends as LeftMesh.
		// Measured on the benchmark at k = 2, n = 64: the trace stopped after 85
		// points of about 320, at a point 0.15 h outside Gamma_h whose nearest
		// face reported a normal component of -2.2e-03.
		//
		// So the outward test comes FIRST and selects the candidates, and the
		// nearest is taken among those. Every point outside a polygon is outside
		// across some face, so this never refuses a genuine band point; and it
		// still refuses a point inside a convex D_h, where no face has the point
		// on its outward side at all.
		int best = -1;
		double bestDistance = 0.0;
		double bestParameter = 0.0;
		double bestFootR = r;
		double bestFootZ = z;

		for ( std::size_t f = 0; f < bandFaces.size(); ++f )
		{
			BoundaryFace const &face = bandFaces[ f ];
			double const dr = face.r1 - face.r0;
			double const dz = face.z1 - face.z0;
			double const square = dr*dr + dz*dz;

			double t = square > 0.0
				? ( ( r - face.r0 )*dr + ( z - face.z0 )*dz )/square : 0.0;
			t = std::min( 1.0, std::max( 0.0, t ) );

			double const onR = face.r0 + t*dr;
			double const onZ = face.z0 + t*dz;

			if ( ( r - onR )*face.normalR + ( z - onZ )*face.normalZ <= 0.0 )
				continue;

			double const distance = std::sqrt( ( r - onR )*( r - onR )
			                                   + ( z - onZ )*( z - onZ ) );
			if ( best < 0 || distance < bestDistance )
			{
				best = static_cast<int>( f );
				bestDistance = distance;
				bestParameter = t;
				bestFootR = onR;
				bestFootZ = onZ;
			}
		}

		if ( best < 0 )
			return -1;

		// Per face rather than mesh wide, so a graded Gamma_h reaches further
		// where its faces are longer -- which is where the band is wider anyway.
		if ( bestDistance > bandReach*bandFaces[ best ].length )
			return -1;

		footR = bestFootR;
		footZ = bestFootZ;
		parameter = bestParameter;
		depth = bestDistance;
		return best;
	}

	bool ContourTracer::extendField( double r, double z, FieldSample &sample ) const
	{
		if ( bandMethod == BandExtension::None || bandFaces.empty() )
			return false;

		double footR = 0.0;
		double footZ = 0.0;
		double parameter = 0.0;
		double depth = 0.0;
		int const which = nearestBandFace( r, z, footR, footZ, parameter, depth );
		if ( which < 0 )
			return false;

		BoundaryFace const &face = bandFaces[ which ];

		sample.element = face.element;
		sample.extended = true;
		sample.depth = depth;
		sample.footR = footR;
		sample.footZ = footZ;

		if ( bandMethod == BandExtension::FluxTaylor )
		{
			// THE FOOT, NOT THE POINT, which is the whole reason this version is
			// safe: x0 lies on the element's own closure, so both fields are read
			// INSIDE the element that owns them and nothing is ever evaluated
			// outside one. What it costs is that grad psi is then frozen at the
			// foot, so the extended field is affine and a contour in the band is
			// a straight line whatever k is.
			mfem::ElementTransformation *trans =
				meshRef.GetElementTransformation( face.element );
			mfem::ElementExtension extender;
			extender.SetElement( *trans );

			mfem::Vector foot( 2 );
			foot( 0 ) = footR;
			foot( 1 ) = footZ;

			mfem::IntegrationPoint ip;
			if ( !extender.TransformBack( foot, ip ) )
				return false;

			mfem::Vector q( 2 );
			fluxField.GetVectorValue( face.element, ip, q );

			// grad_bar( psi ) = r q, with r AT THE FOOT -- the point the flux was
			// actually read at. Using the band point's own r instead is the trap
			// CLAUDE.md records against sampleCoefficient(), and it is a factor
			// of 1.7e5 there on a quantity carrying r^2.
			sample.ip = ip;
			sample.psi = potentialField.GetValue( face.element, ip )
			             + footR*( q( 0 )*( r - footR ) + q( 1 )*( z - footZ ) );
			sample.qR = q( 0 );
			sample.qZ = q( 1 );
			return true;
		}

		// THE TRANSFER-PATH LIFT. C u = -grad_bar( psi ) is a gradient, so its
		// line integral from the band point to ANY point of Gamma is
		// psi( p ) - g there, whatever path is taken. That is the identity the
		// whole curved path rests on, and it is stated for arbitrary endpoints:
		// mfem::PathIntegral takes x and xbar and nothing else.
		//
		// mfem::PathLiftCoefficient -- which INVERSION-PLAN.md section 4.3 names
		// as the primitive -- is NOT usable here: it evaluates on a
		// FaceElementTransformations at that face's own integration point and
		// answers "what is phi_h on Gamma_h", which is eta_5's question. The
		// usable primitive is one level down and public.
		mfem::Vector xbar( 2 );
		{
			mfem::FaceElementTransformations *ftr =
				meshRef.GetBdrFaceTransformations( face.boundaryElement );
			if ( !ftr )
				return false;

			mfem::IntegrationPoint at;
			at.Set1w( parameter, 1.0 );
			bandPath->Endpoint( *ftr, at, xbar );
		}

		// AFTER Endpoint() AND NEVER BEFORE. It resets the mesh's shared
		// transformations, which HDGExtensionIntegrator::ComputeLift() says in as
		// many words -- "the element is handed to the extension only below".
		mfem::ElementTransformation *trans =
			meshRef.GetElementTransformation( face.element );
		mfem::ElementExtension extender;
		extender.SetElement( *trans );

		int const element = face.element;
		int const order = bandLineOrder >= 0
			? bandLineOrder
			: 2*fluxField.FESpace()->GetOrder( element ) + 2;
		mfem::IntegrationRule const &lineRule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, order );

		bool reached = true;
		auto flux = [ & ]( mfem::Vector const &y, mfem::Vector &value )
		{
			// E_h( q_h ): the element's OWN polynomial, evaluated outside it.
			// mfem::ElementExtension is what does not clamp the reference point,
			// which the ordinary TransformBack does -- and a clamped point turns
			// the extension into a constant, silently.
			mfem::IntegrationPoint eip;
			if ( !extender.TransformBack( y, eip ) )
			{
				reached = false;
				value = 0.0;
				return;
			}
			fluxField.GetVectorValue( element, eip, value );

			// C u = -grad_bar( psi ) = -r q, in MEQ's sign convention. NOT
			// DarcyForm's: transferredDatum() hands mfem::PathLiftCoefficient the
			// raw block because that class re-runs the integrator the assembly
			// used, and CLAUDE.md records that feeding it flux() instead returns
			// -psi rather than psi. Here the sign is written out, so the
			// convention is chosen rather than inherited.
			value *= -y( 0 );
		};

		mfem::Vector here( 2 );
		here( 0 ) = r;
		here( 1 ) = z;

		double const lift = mfem::PathIntegral( flux, here, xbar, lineRule );
		if ( !reached || !std::isfinite( lift ) )
			return false;

		mfem::IntegrationPoint eip;
		if ( !extender.TransformBack( here, eip ) )
			return false;

		mfem::Vector q( 2 );
		fluxField.GetVectorValue( element, eip, q );

		sample.ip = eip;
		sample.psi = bandDatum( xbar ) + lift;
		sample.qR = q( 0 );
		sample.qZ = q( 1 );
		return true;
	}

	bool ContourTracer::sampleField( double r, double z, int hint,
	                                 FieldSample &sample, int &fallbacks ) const
	{
		// THE SEAM. Everything in this file that wants psi or q at a physical
		// point comes through here, and the band of IN-0's second half is
		// answered here and nowhere else.
		auto inElement = [ & ]( int element, mfem::IntegrationPoint const &ip )
		{
			mfem::Vector value( 2 );
			fluxField.GetVectorValue( element, ip, value );

			sample.element = element;
			sample.ip = ip;
			sample.psi = potentialField.GetValue( element, ip );
			sample.qR = value( 0 );
			sample.qZ = value( 1 );
			sample.extended = false;
			sample.depth = 0.0;
			sample.footR = 0.0;
			sample.footZ = 0.0;
		};

		bool const band = bandMethod != BandExtension::None;

		int element = -1;
		mfem::IntegrationPoint ip;
		if ( locate( r, z, hint, element, ip, fallbacks, !band ) )
		{
			inElement( element, ip );
			return true;
		}

		if ( band )
		{
			FieldSample extendedSample;
			if ( extendField( r, z, extendedSample ) )
			{
				sample = extendedSample;
				return true;
			}

			// The band turned it down, so the failed walk above was a failed
			// walk and not a point outside the mesh. The last resort is still
			// owed, and it is counted the way it always was.
			if ( locate( r, z, hint, element, ip, fallbacks, true ) )
			{
				inElement( element, ip );
				return true;
			}
		}

		return false;
	}

	bool ContourTracer::sampleAt( double r, double z, double &psi, double &qR,
	                              double &qZ, int &hint ) const
	{
		FieldSample sample;
		int fallbacks = 0;
		if ( !sampleField( r, z, hint, sample, fallbacks ) )
			return false;

		hint = sample.element;
		psi = sample.psi;
		qR = sample.qR;
		qZ = sample.qZ;
		return true;
	}

	bool ContourTracer::sampleAt( double r, double z, double &psi, double &qR,
	                              double &qZ ) const
	{
		int hint = -1;
		return sampleAt( r, z, psi, qR, qZ, hint );
	}

	bool ContourTracer::sampleAt( double r, double z, double &psi, double &qR,
	                              double &qZ, int &hint, bool &extended ) const
	{
		FieldSample sample;
		int fallbacks = 0;
		extended = false;
		if ( !sampleField( r, z, hint, sample, fallbacks ) )
			return false;

		hint = sample.element;
		psi = sample.psi;
		qR = sample.qR;
		qZ = sample.qZ;
		extended = sample.extended;
		return true;
	}

	bool ContourTracer::correct( double level, double target, double maxMove,
	                             double &r, double &z, int hint,
	                             FieldSample &sample, int &iterations,
	                             bool &stalled, int &fallbacks ) const
	{
		// @a target is the ABSOLUTE residual to stop at, computed once by the
		// caller from potentialScale(). It is a parameter rather than a lookup
		// for the reason CriticalPointFinder::rootInElement() gives for the
		// same choice: the scale costs a pass over every potential dof, and
		// taking it once per corrector iteration would make a trace quadratic
		// in the mesh for no gain whatever.
		//
		// THE CORRECTOR CANNOT ALWAYS REACH ITS TOLERANCE, AND THE REASON IS THE
		// DISCONTINUITY RATHER THAN THE ITERATION. { psi_h = c } is a union of
		// per-element arcs offset from each other by the DG jump, so a point
		// that lands within jump/|grad psi| of a face may be on NEITHER arc: the
		// Newton step computed in element A pushes it across the face into B,
		// the step computed in B pushes it back, and the two residuals alternate
		// without ever falling below a tolerance tighter than the jump.
		// Measured, it is rare -- it needs the point to land inside a band about
		// 2e-5 wide on a contour of length 1.7 -- and it becomes commoner as
		// Delta_s falls, simply because more points are placed.
		//
		// So the iteration keeps the BEST iterate it has seen and accepts it
		// when four steps in a row fail to improve on it, reporting @a stalled
		// so that Contour::stalledCorrections can count them. That is honest
		// rather than convenient: on a discontinuous field a point cannot be
		// closer to "the" level set than the field's own ambiguity at a face,
		// and ContourPoint::residual says how close it got. What it must NOT do
		// is wander, so the iterate is refused if it leaves @a maxMove of where
		// the predictor put it -- a corrector that has travelled a whole step is
		// not correcting.
		double const originR = r;
		double const originZ = z;

		double best = std::numeric_limits<double>::infinity();
		double bestR = r;
		double bestZ = z;
		FieldSample bestSample;
		bool haveBest = false;
		int sinceImprovement = 0;

		stalled = false;
		iterations = 0;
		for ( int iteration = 0; iteration <= maxCorrectorIterations; ++iteration )
		{
			if ( !sampleField( r, z, hint, sample, fallbacks ) )
				break;
			hint = sample.element;

			double const residual = level - sample.psi;
			if ( std::abs( residual ) <= target )
				return true;

			if ( std::abs( residual ) < best )
			{
				best = std::abs( residual );
				bestR = r;
				bestZ = z;
				bestSample = sample;
				haveBest = true;
				sinceImprovement = 0;
			}
			else if ( ++sinceImprovement >= 4 )
			{
				break;
			}

			// grad_bar( psi ) = r q, and the minimum-norm Newton step is
			// x <- x + grad( c - psi ) / | grad |^2. Written in terms of q the r
			// appears once in the numerator and twice in the denominator.
			double const magnitude = std::sqrt( sample.qR*sample.qR
			                                    + sample.qZ*sample.qZ );
			if ( !( magnitude > 0.0 ) || !( r > 0.0 ) )
				break;

			double const factor = residual/( r*magnitude*magnitude );
			r += factor*sample.qR;
			z += factor*sample.qZ;
			++iterations;

			if ( !std::isfinite( r ) || !std::isfinite( z ) )
				break;

			double const goneR = r - originR;
			double const goneZ = z - originZ;
			if ( goneR*goneR + goneZ*goneZ > maxMove*maxMove )
				break;
		}

		if ( !haveBest )
			return false;

		r = bestR;
		z = bestZ;
		sample = bestSample;
		stalled = true;
		return true;
	}

	double ContourTracer::faceJump( double r0, double z0, int element0,
	                                double r1, double z1, int element1,
	                                int &fallbacks ) const
	{
		if ( element0 == element1 )
			return 0.0;

		// Bisect the segment for the parameter at which the located element
		// stops being element0. Sixty halvings puts the two ends of the bracket
		// about 1e-18 of a segment apart, so the two evaluations below are at
		// the same point for every purpose -- which is what lets each be taken
		// INSIDE its own element rather than extrapolated across the face.
		double lo = 0.0;
		double hi = 1.0;
		int loElement = element0;
		int hiElement = element1;
		mfem::IntegrationPoint loIp;
		mfem::IntegrationPoint hiIp;

		{
			int element = -1;
			if ( !locate( r0, z0, element0, element, loIp, fallbacks ) )
				return 0.0;
			if ( !locate( r1, z1, element1, element, hiIp, fallbacks ) )
				return 0.0;
		}

		for ( int i = 0; i < 60; ++i )
		{
			double const mid = 0.5*( lo + hi );
			double const rm = ( 1.0 - mid )*r0 + mid*r1;
			double const zm = ( 1.0 - mid )*z0 + mid*z1;

			int element = -1;
			mfem::IntegrationPoint ip;
			if ( !locate( rm, zm, loElement, element, ip, fallbacks ) )
				return 0.0;

			if ( element == loElement )
			{
				lo = mid;
				loIp = ip;
			}
			else
			{
				hi = mid;
				hiElement = element;
				hiIp = ip;
			}
		}

		double const a = potentialField.GetValue( loElement, loIp );
		double const b = potentialField.GetValue( hiElement, hiIp );
		return std::abs( a - b );
	}

	Contour ContourTracer::trace( double level, double startR,
	                              double startZ ) const
	{
		// The ONE legitimate Mesh::FindPoints call: locating the seed, where
		// there is by definition no previous element to walk from. It is
		// deliberately not counted in Contour::fallbackLocations, which is a
		// measurement of whether the WALK is working and would otherwise read
		// one on every trace and mean nothing.
		int seed = -1;
		{
			mfem::DenseMatrix matrix( 2, 1 );
			matrix( 0, 0 ) = startR;
			matrix( 1, 0 ) = startZ;

			mfem::Array<int> found;
			mfem::Array<mfem::IntegrationPoint> ips;
			meshRef.FindPoints( matrix, found, ips, false );
			if ( found.Size() > 0 )
				seed = found[ 0 ];
		}

		if ( seed < 0 )
		{
			std::ostringstream message;
			message << "ContourTracer::trace: ( " << startR << ", " << startZ
			        << " ) is not in the mesh";
			throw std::runtime_error( message.str() );
		}

		return traceFrom( level, startR, startZ, seed );
	}

	Contour ContourTracer::traceFrom( double level, double startR, double startZ,
	                                  int seed ) const
	{
		Contour contour;
		contour.level = level;
		contour.bandExtension = bandMethod;

		int fallbacks = 0;
		int iterations = 0;
		FieldSample sample;

		double const target = tolerance*potentialScale();
		contour.correctorTarget = target;

		double r = startR;
		double z = startZ;
		bool stalled = false;
		if ( !correct( level, target, elementSize( seed ), r, z, seed, sample,
		               iterations, stalled, fallbacks ) )
		{
			std::ostringstream message;
			message << "ContourTracer::trace: the corrector could not reach psi = "
			        << level << " from ( " << startR << ", " << startZ << " )";
			throw std::runtime_error( message.str() );
		}

		double const nominal = step > 0.0 ? step : 0.5*elementSize( sample.element );
		if ( !( nominal > 0.0 ) )
			throw std::runtime_error( "ContourTracer::trace: the step came out non-positive" );

		contour.nominalStep = nominal;
		contour.shortestStep = std::numeric_limits<double>::infinity();
		contour.longestStep = 0.0;

		auto appendPoint = [ & ]( FieldSample const &at, double rIn, double zIn,
		                          double arcLength, int correctorIterations )
		{
			double const magnitude = std::sqrt( at.qR*at.qR + at.qZ*at.qZ );
			ContourPoint p;
			p.r = rIn;
			p.z = zIn;
			p.tangentR = -at.qZ/magnitude;
			p.tangentZ = at.qR/magnitude;
			p.fluxMagnitude = magnitude;
			p.arcLength = arcLength;
			p.element = at.element;
			p.correctorIterations = correctorIterations;
			p.residual = std::abs( at.psi - level );

			// THE BAND FLAG, PER POINT. Contour::extendedPoints is a summary of
			// it and not a substitute for it: a consumer dropping band points
			// before an error norm needs to know WHICH, and CLAUDE.md records a
			// count where a mask was wanted as the other half of a real defect
			// in the .nc.
			p.extended = at.extended;
			p.bandDepth = at.depth;
			if ( at.extended )
			{
				++contour.extendedPoints;
				contour.deepestBandPoint =
					std::max( contour.deepestBandPoint, at.depth );
			}

			contour.points.push_back( p );

			contour.correctorIterationsTotal += correctorIterations;
			contour.worstCorrectorIterations =
				std::max( contour.worstCorrectorIterations, correctorIterations );
			contour.worstResidual = std::max( contour.worstResidual, p.residual );
		};

		if ( stalled )
			++contour.stalledCorrections;
		appendPoint( sample, r, z, 0.0, iterations );

		double const firstR = contour.points.front().r;
		double const firstZ = contour.points.front().z;
		double const firstTangentR = contour.points.front().tangentR;
		double const firstTangentZ = contour.points.front().tangentZ;

		double curvature = 0.0;
		double arcLength = 0.0;
		double jumpTotal = 0.0;
		contour.status = ContourStatus::TooLong;

		while ( static_cast<int>( contour.points.size() ) < maxPoints )
		{
			ContourPoint const &here = contour.points.back();
			double const tangentR = here.tangentR;
			double const tangentZ = here.tangentZ;

			// The step. Curvature control equidistributes TURNING, which is the
			// geometric variant of the header; with targetTurn at zero the step
			// is exactly the nominal one, which is what a convergence study in
			// Delta_s needs.
			double stepLength = nominal;
			if ( targetTurn > 0.0 )
			{
				if ( curvature > 0.0 )
					stepLength = std::min( nominal, targetTurn/curvature );
				stepLength = std::max( stepLength, minStepFraction*nominal );
				if ( localStepCeiling > 0.0 )
					stepLength = std::min( stepLength,
					                       localStepCeiling*elementSize( here.element ) );
			}

			// Closure. The turning gate is what stops a multi-circuit trace
			// closing on its first pass; the proximity gate is what stops it
			// closing at the far side of the contour when the level set happens
			// to come back near the start.
			double const gapR = firstR - here.r;
			double const gapZ = firstZ - here.z;
			double const along = gapR*tangentR + gapZ*tangentZ;
			double const gap = std::sqrt( gapR*gapR + gapZ*gapZ );

			bool const turned = std::abs( contour.turning )
			                    >= circuits*twoPi - 0.5*M_PI;

			if ( turned && along > 0.0 && gap <= 1.5*stepLength )
			{
				// THE FINAL STEP IS SHORTENED TO AIM AT THE START POINT, not
				// taken at full length and left as a stub. The predictor is
				// given the projection of the gap onto the tangent, so what the
				// corrector leaves is O( kappa^2 Delta_s^3 ) along the curve
				// rather than O( Delta_s ).
				double closeR = here.r + along*tangentR;
				double closeZ = here.z + along*tangentZ;

				FieldSample closing;
				int closingIterations = 0;
				bool closingStalled = false;
				if ( correct( level, target, std::max( along, 1.0e-300 ), closeR,
				              closeZ, here.element, closing, closingIterations,
				              closingStalled, fallbacks ) )
				{
					// Measured BEFORE the curve is closed. The normal component
					// is the one that says the corrector is doing its job; see
					// the header on why "machine precision" needed refining.
					double const errorR = closeR - firstR;
					double const errorZ = closeZ - firstZ;
					contour.closureTangential =
						std::abs( errorR*firstTangentR + errorZ*firstTangentZ );
					contour.closureNormal =
						std::abs( -errorR*firstTangentZ + errorZ*firstTangentR );
				}

				double const finalR = firstR - here.r;
				double const finalZ = firstZ - here.z;
				double const finalStep = std::sqrt( finalR*finalR + finalZ*finalZ );

				// A BAND POINT HAS NO SECOND SIDE TO DISAGREE WITH, so the jump
				// is not measured across it. FieldSample::element for a band
				// point is the element whose polynomial was extended rather than
				// one containing the point, so two such points would look like a
				// face crossing when they are not, and faceJump()'s bisection
				// would be locating points that are in no element at all.
				if ( measureFaceJumps && !here.extended
				     && !contour.points.front().extended
				     && contour.points.front().element != here.element )
				{
					double const jump = faceJump( here.r, here.z, here.element,
					                              firstR, firstZ,
					                              contour.points.front().element,
					                              fallbacks );
					if ( jump > 0.0 )
					{
						++contour.faceCrossings;
						jumpTotal += jump;
						contour.worstFaceJump = std::max( contour.worstFaceJump, jump );
					}
				}

				// The closing step turns too, and leaving it out would make the
				// reported turning short of a circuit by whatever the last step
				// was worth -- which reads as a trace that did not quite get
				// round.
				contour.turning += angleStep( std::atan2( tangentZ, tangentR ),
				                              std::atan2( firstTangentZ, firstTangentR ) );

				arcLength += finalStep;
				contour.shortestStep = std::min( contour.shortestStep, finalStep );
				contour.longestStep = std::max( contour.longestStep, finalStep );

				// The closing point IS the first point, exactly.
				ContourPoint closingPoint = contour.points.front();
				closingPoint.arcLength = arcLength;
				contour.points.push_back( closingPoint );

				contour.status = ContourStatus::Closed;
				break;
			}

			// The predictor, and the safeguard on it. Halving a step whose
			// turning came out larger than a third of pi is what keeps
			// angleStep()'s branch choice honest; it is skipped where curvature
			// control is off, because there the caller has asked for a fixed
			// step and silently changing it would corrupt the study that asked.
			bool accepted = false;
			for ( int attempt = 0; attempt < 8; ++attempt )
			{
				double nextR = here.r + stepLength*tangentR;
				double nextZ = here.z + stepLength*tangentZ;

				FieldSample next;
				int nextIterations = 0;
				bool nextStalled = false;
				if ( !correct( level, target, stepLength, nextR, nextZ,
				               here.element, next, nextIterations, nextStalled,
				               fallbacks ) )
				{
					// Distinguish "left the mesh" from "the corrector gave up
					// where it still had a field". The first is the outcome
					// v0-legacy printed a message about and then hid.
					int probe = -1;
					mfem::IntegrationPoint ip;
					contour.status =
						locate( nextR, nextZ, here.element, probe, ip, fallbacks )
						  ? ContourStatus::Stalled : ContourStatus::LeftMesh;
					break;
				}

				double const magnitude = std::sqrt( next.qR*next.qR + next.qZ*next.qZ );
				if ( !( magnitude > 0.0 ) )
				{
					contour.status = ContourStatus::Stalled;
					break;
				}

				double const nextTangentR = -next.qZ/magnitude;
				double const nextTangentZ = next.qR/magnitude;
				double const turn = angleStep( std::atan2( tangentZ, tangentR ),
				                               std::atan2( nextTangentZ, nextTangentR ) );

				if ( targetTurn > 0.0 && std::abs( turn ) > M_PI/3.0
				     && stepLength > minStepFraction*nominal )
				{
					stepLength = std::max( 0.5*stepLength, minStepFraction*nominal );
					continue;
				}

				double const moveR = nextR - here.r;
				double const moveZ = nextZ - here.z;
				double const moved = std::sqrt( moveR*moveR + moveZ*moveZ );

				if ( measureFaceJumps && !here.extended && !next.extended
				     && next.element != here.element )
				{
					double const jump = faceJump( here.r, here.z, here.element,
					                              nextR, nextZ, next.element, fallbacks );
					if ( jump > 0.0 )
					{
						++contour.faceCrossings;
						jumpTotal += jump;
						contour.worstFaceJump = std::max( contour.worstFaceJump, jump );
					}
				}

				contour.turning += turn;
				curvature = moved > 0.0 ? std::abs( turn )/moved : 0.0;
				arcLength += moved;
				contour.shortestStep = std::min( contour.shortestStep, moved );
				contour.longestStep = std::max( contour.longestStep, moved );

				if ( nextStalled )
					++contour.stalledCorrections;
				appendPoint( next, nextR, nextZ, arcLength, nextIterations );
				accepted = true;
				break;
			}

			if ( !accepted )
			{
				if ( contour.status == ContourStatus::TooLong )
					contour.status = ContourStatus::Stalled;
				break;
			}
		}

		if ( contour.faceCrossings > 0 )
			contour.meanFaceJump = jumpTotal/contour.faceCrossings;
		if ( !std::isfinite( contour.shortestStep ) )
			contour.shortestStep = 0.0;

		contour.circuits = static_cast<int>(
			std::lround( std::abs( contour.turning )/twoPi ) );
		contour.fallbackLocations = fallbacks;
		return contour;
	}

	Contour ContourTracer::traceFromAxis( double level,
	                                      CriticalPoint const &axis ) const
	{
		// Walk out from the axis until psi_h brackets the level, then bisect.
		// A ray is the cheapest bracket there is and the axis is where psi is
		// extremal, so a level between the axis value and the boundary value is
		// crossed along a ray that has room for it.
		//
		// SEVERAL RAYS, AND THAT IS A MEASUREMENT RATHER THAN CAUTION. An
		// earlier version took +r alone, which is the obvious choice and fails
		// on an ordinary benchmark: the Solov'ev NSTX axis sits at r = 1.318 on
		// the standard box [ 0.6, 1.4 ] x [ -0.6, 0.6 ], so there is 0.08 of
		// room to the right of it and an outer surface reaches the edge of the
		// mesh along +r while being comfortably inside it in every other
		// direction. The failure is loud -- the level is simply not attained --
		// but it is a property of where the box was drawn and not of the
		// surface, and refusing to trace a surface that exists would be wrong.
		//
		// This is NOT an assumption of star-shapedness: only ONE crossing on ONE
		// ray is needed, and it is a seed for the tracer rather than a
		// parametrisation. fitByAngle() is where star-shapedness becomes a
		// hypothesis, and it reports it.
		int hint = -1;
		double psi = 0.0;
		double qR = 0.0;
		double qZ = 0.0;

		if ( !sampleAt( axis.r, axis.z, psi, qR, qZ, hint ) )
			throw std::runtime_error( "ContourTracer::traceFromAxis: the axis is not in the mesh" );

		double const axisPsi = psi;
		int const home = hint;
		double const probe = 0.25*elementSize( hint );
		if ( !( probe > 0.0 ) )
			throw std::runtime_error( "ContourTracer::traceFromAxis: the mesh has no scale" );

		int const rays = 8;
		for ( int ray = 0; ray < rays; ++ray )
		{
			double const theta = twoPi*ray/rays;
			double const cosine = std::cos( theta );
			double const sine = std::sin( theta );

			hint = home;
			double lo = 0.0;
			double hi = 0.0;
			bool bracketed = false;

			for ( int i = 1; i < 1000000; ++i )
			{
				double const radius = i*probe;
				if ( !sampleAt( axis.r + radius*cosine, axis.z + radius*sine,
				                psi, qR, qZ, hint ) )
					break;

				if ( ( axisPsi - level )*( psi - level ) <= 0.0 )
				{
					lo = ( i - 1 )*probe;
					hi = radius;
					bracketed = true;
					break;
				}
			}

			if ( !bracketed )
				continue;

			for ( int i = 0; i < 60; ++i )
			{
				double const mid = 0.5*( lo + hi );
				if ( !sampleAt( axis.r + mid*cosine, axis.z + mid*sine,
				                psi, qR, qZ, hint ) )
					break;
				if ( ( axisPsi - level )*( psi - level ) <= 0.0 )
					hi = mid;
				else
					lo = mid;
			}

			// The bracket search has already located the seed, so the walk
			// starts with a hint and the trace makes no FindPoints call at all.
			double const found = 0.5*( lo + hi );
			return traceFrom( level, axis.r + found*cosine, axis.z + found*sine,
			                  hint );
		}

		std::ostringstream message;
		message << "ContourTracer::traceFromAxis: psi = " << level
		        << " is not attained on any of " << rays << " rays from the axis "
		        << "at ( " << axis.r << ", " << axis.z << " ), where psi = "
		        << axisPsi << ". Either the level is outside the range psi_h "
		        << "takes on this mesh, or every ray leaves the mesh before "
		        << "reaching it";
		throw std::runtime_error( message.str() );
	}

	double AngleParametrisation::length() const
	{
		std::size_t const n = speed.size();
		if ( n == 0 )
			return 0.0;

		double total = 0.0;
		for ( std::size_t j = 0; j < n; ++j )
			total += speed[ j ];
		return ( twoPi/n )*total;
	}

	double AngleParametrisation::differencedLength() const
	{
		std::size_t const n = radius.size();
		if ( n < 3 )
			return 0.0;

		// THE TRAP. The same periodic trapezoid rule, fed a Jacobian obtained by
		// differencing neighbouring node positions. Nothing about the rule says
		// it has been reduced to second order, which is the whole point of
		// keeping it.
		double const dTheta = twoPi/n;
		double total = 0.0;
		for ( std::size_t j = 0; j < n; ++j )
		{
			std::size_t const next = ( j + 1 )%n;
			std::size_t const previous = ( j + n - 1 )%n;
			double const derivative = ( radius[ next ] - radius[ previous ] )
			                          /( 2.0*dTheta );
			total += std::sqrt( derivative*derivative + radius[ j ]*radius[ j ] );
		}
		return dTheta*total;
	}

	double AngleParametrisation::chordLength() const
	{
		std::size_t const n = pointR.size();
		if ( n < 2 )
			return 0.0;

		double total = 0.0;
		for ( std::size_t j = 0; j < n; ++j )
		{
			std::size_t const next = ( j + 1 )%n;
			double const dr = pointR[ next ] - pointR[ j ];
			double const dz = pointZ[ next ] - pointZ[ j ];
			total += std::sqrt( dr*dr + dz*dz );
		}
		return total;
	}

	std::vector<double> AngleParametrisation::spectralRadiusPrime() const
	{
		// The Fourier differentiation matrix on an equispaced periodic grid,
		// applied directly. O( N^2 ) and needs no transform; at the sizes IN-1
		// runs it is free. Trefethen, Spectral Methods in MATLAB, ch. 3:
		// even N uses the cotangent form and odd N the cosecant one, and using
		// the wrong one is a silently wrong derivative rather than a failure.
		std::size_t const n = radius.size();
		std::vector<double> derivative( n, 0.0 );
		if ( n < 2 )
			return derivative;

		bool const even = ( n%2 == 0 );
		for ( std::size_t j = 0; j < n; ++j )
		{
			double total = 0.0;
			for ( std::size_t k = 0; k < n; ++k )
			{
				if ( k == j )
					continue;
				double const half = M_PI*( static_cast<double>( j )
				                           - static_cast<double>( k ) )/n;
				double const sign = ( ( j + k )%2 == 0 ) ? 1.0 : -1.0;
				double const weight = even ? 0.5*sign/std::tan( half )
				                           : 0.5*sign/std::sin( half );
				total += weight*radius[ k ];
			}
			derivative[ j ] = total;
		}
		return derivative;
	}

	AngleParametrisation ContourTracer::fitByAngle( Contour const &contour,
	                                               CriticalPoint const &axis,
	                                               std::size_t count,
	                                               double floor ) const
	{
		if ( !contour.closed() )
			throw std::invalid_argument(
				"ContourTracer::fitByAngle: the contour is not closed" );
		if ( count < 4 )
			throw std::invalid_argument(
				"ContourTracer::fitByAngle: need at least four angles" );

		AngleParametrisation fit;
		fit.axisR = axis.r;
		fit.axisZ = axis.z;
		fit.level = contour.level;

		// The traced curve, in polar coordinates about the axis, unwrapped. The
		// unwrapping is what turns "is this star-shaped" into a measurement: a
		// non-monotone unwrapped angle is a ray that meets the curve more than
		// once, and it is reported rather than survived.
		std::size_t const traced = contour.points.size();
		std::vector<double> angle( traced, 0.0 );
		std::vector<double> radius( traced, 0.0 );

		double running = 0.0;
		for ( std::size_t i = 0; i < traced; ++i )
		{
			double const dr = contour.points[ i ].r - axis.r;
			double const dz = contour.points[ i ].z - axis.z;
			radius[ i ] = std::sqrt( dr*dr + dz*dz );

			double const raw = std::atan2( dz, dr );
			if ( i == 0 )
				running = raw;
			else
				running += angleStep( running, raw );
			angle[ i ] = running;
		}

		// Orient increasing, so that the bracket search below is a plain binary
		// search whichever way round the trace went.
		bool const descending = angle.back() < angle.front();
		if ( descending )
		{
			std::reverse( angle.begin(), angle.end() );
			std::reverse( radius.begin(), radius.end() );
		}

		for ( std::size_t i = 1; i < traced; ++i )
		{
			double const increment = angle[ i ] - angle[ i - 1 ];
			if ( increment < 0.0 )
			{
				fit.angleMonotone = false;
				fit.worstBacktrack = std::max( fit.worstBacktrack, -increment );
			}
		}

		if ( !fit.angleMonotone )
		{
			std::ostringstream message;
			message << "ContourTracer::fitByAngle: the traced contour is not "
			        << "star-shaped about ( " << axis.r << ", " << axis.z
			        << " ) -- the polar angle turns back by " << fit.worstBacktrack
			        << " rad. INVERSION-PLAN.md section 3.4 records that ray "
			        << "methods fail on indented cross-sections; this is that "
			        << "hypothesis failing rather than a defect";
			throw std::runtime_error( message.str() );
		}

		double const start = angle.front();
		double const span = angle.back() - angle.front();
		if ( !( span > 0.0 ) )
			throw std::runtime_error(
				"ContourTracer::fitByAngle: the traced contour subtends no angle" );

		fit.radius.assign( count, 0.0 );
		fit.pointR.assign( count, 0.0 );
		fit.pointZ.assign( count, 0.0 );
		fit.fluxR.assign( count, 0.0 );
		fit.fluxZ.assign( count, 0.0 );
		fit.radiusPrime.assign( count, 0.0 );
		fit.speed.assign( count, 0.0 );
		fit.crossing.assign( count, 0.0 );
		fit.extended.assign( count, 0 );
		fit.transversality = std::numeric_limits<double>::infinity();

		double const target = tolerance*potentialScale();
		int fallbacks = 0;
		int hint = contour.points.front().element;

		for ( std::size_t j = 0; j < count; ++j )
		{
			double const theta = twoPi*j/count;

			// Fold the target angle into the span the trace covers, then find
			// the two traced points that straddle it. The bracket comes from the
			// TRACE, which is why this does not assume star-shapedness globally.
			double folded = theta;
			while ( folded < start )
				folded += twoPi;
			while ( folded > start + span )
				folded -= twoPi;
			if ( folded < start )
				folded = start;

			std::size_t lo = 0;
			std::size_t hi = traced - 1;
			while ( hi - lo > 1 )
			{
				std::size_t const mid = ( lo + hi )/2;
				if ( angle[ mid ] <= folded )
					lo = mid;
				else
					hi = mid;
			}

			double const width = angle[ lo + 1 ] - angle[ lo ];
			double const weight = width > 0.0 ? ( folded - angle[ lo ] )/width : 0.0;
			double guess = ( 1.0 - weight )*radius[ lo ] + weight*radius[ lo + 1 ];

			double const margin = 0.05;
			double bracketLo = std::min( radius[ lo ], radius[ lo + 1 ] )*( 1.0 - margin );
			double bracketHi = std::max( radius[ lo ], radius[ lo + 1 ] )*( 1.0 + margin );

			double const cosine = std::cos( theta );
			double const sine = std::sin( theta );

			// The band flag of the LAST evaluation, which after the search is
			// the accepted node's. It is carried out of the lambda rather than
			// returned, because both the Newton and the bisection branch below
			// leave the answer in whichever evaluation happened to converge.
			bool nodeExtended = false;
			double nodeDepth = 0.0;

			// THE BEST ITERATE, AND WHY THIS RAY NEWTON KEEPS ONE.
			//
			// The corrector above already does this and the reason is the same:
			// { psi_h = c } is a union of per-element arcs offset by the DG jump,
			// so a RAY crossing a face where c falls inside the jump has NO point
			// on it with psi_h = c at all, and no tolerance tighter than the jump
			// is attainable there. Demanding one and throwing is refusing to
			// answer a question that has an answer to within the field's own
			// ambiguity -- measured, at k = 1 on the raw pairing it threw on
			// EVERY mesh from n = 12 to 32, and the probability rises with the
			// angle count simply because more rays are more chances.
			//
			// So the search keeps the closest iterate it has seen anywhere along
			// the ray and accepts it if neither the Newton nor the bisection
			// reaches the tolerance, counting it in AngleParametrisation::
			// stalledRays and leaving worstResidual to say how close it got.
			//
			// WHAT IS NOT SOFTENED IS THE BRACKET. A ray that never brackets the
			// level at all is a different thing entirely -- the surface is not
			// star-shaped about the axis, or the level is not attained on that
			// ray -- and it still throws, because there is no best iterate to
			// speak of and a number returned from it would be fiction.
			double bestRadius = 0.0;
			double bestResidual = std::numeric_limits<double>::infinity();
			double bestPsi = 0.0;
			double bestQR = 0.0;
			double bestQZ = 0.0;
			bool bestExtended = false;
			double bestDepth = 0.0;
			bool haveBest = false;

			auto evaluate = [ & ]( double rho, double &psi, double &qR,
			                       double &qZ ) -> bool
			{
				FieldSample sample;
				if ( !sampleField( axis.r + rho*cosine, axis.z + rho*sine, hint,
				                   sample, fallbacks ) )
					return false;
				hint = sample.element;
				psi = sample.psi;
				qR = sample.qR;
				qZ = sample.qZ;
				nodeExtended = sample.extended;
				nodeDepth = sample.depth;

				double const distance = std::abs( sample.psi - contour.level );
				if ( distance < bestResidual )
				{
					bestResidual = distance;
					bestRadius = rho;
					bestPsi = sample.psi;
					bestQR = sample.qR;
					bestQZ = sample.qZ;
					bestExtended = sample.extended;
					bestDepth = sample.depth;
					haveBest = true;
				}

				return true;
			};

			double psi = 0.0;
			double qR = 0.0;
			double qZ = 0.0;
			bool converged = false;
			int iterations = 0;

			for ( ; iterations < maxCorrectorIterations; ++iterations )
			{
				if ( !evaluate( guess, psi, qR, qZ ) )
					break;

				double const residual = psi - contour.level;
				if ( std::abs( residual ) <= target )
				{
					converged = true;
					break;
				}

				// The derivative along the ray, POINTWISE FROM q and not
				// differenced: d psi / d rho = grad_bar( psi ) . u = r ( q . u ).
				double const rCoord = axis.r + guess*cosine;
				double const slope = rCoord*( qR*cosine + qZ*sine );
				if ( !( std::abs( slope ) > 0.0 ) )
					break;

				double next = guess - residual/slope;
				if ( !std::isfinite( next ) || next < bracketLo || next > bracketHi )
					break;
				guess = next;
			}

			if ( !converged )
			{
				// Newton left its bracket or ran out. Bisection on the same
				// bracket, which is why the bracket is kept rather than only the
				// seed.
				++fit.bisections;

				double aLo = bracketLo;
				double aHi = bracketHi;
				double psiLo = 0.0;
				double psiHi = 0.0;
				if ( !evaluate( aLo, psiLo, qR, qZ ) || !evaluate( aHi, psiHi, qR, qZ )
				     || ( psiLo - contour.level )*( psiHi - contour.level ) > 0.0 )
				{
					std::ostringstream message;
					message << "ContourTracer::fitByAngle: no crossing of psi = "
					        << contour.level << " on the ray at theta = " << theta
					        << " within [ " << bracketLo << ", " << bracketHi
					        << " ] taken from the traced curve";
					throw std::runtime_error( message.str() );
				}

				for ( int i = 0; i < 200; ++i )
				{
					double const mid = 0.5*( aLo + aHi );
					if ( !evaluate( mid, psi, qR, qZ ) )
						break;
					if ( std::abs( psi - contour.level ) <= target )
					{
						guess = mid;
						converged = true;
						break;
					}
					if ( ( psiLo - contour.level )*( psi - contour.level ) <= 0.0 )
					{
						aHi = mid;
						psiHi = psi;
					}
					else
					{
						aLo = mid;
						psiLo = psi;
					}
					guess = mid;
				}

				if ( !converged )
				{
					// Neither route met the tolerance. Take the closest point the
					// ray ever reached and say so; see the note beside the best
					// iterate above for why that is honest rather than
					// convenient. haveBest is false only if every evaluation on
					// this ray left the mesh, which is not a tolerance problem.
					if ( !haveBest )
					{
						std::ostringstream message;
						message << "ContourTracer::fitByAngle: every evaluation on"
						        << " the ray at theta = " << theta << " fell"
						        << " outside the field, so there is no point on it"
						        << " at all";
						throw std::runtime_error( message.str() );
					}

					guess = bestRadius;
					psi = bestPsi;
					qR = bestQR;
					qZ = bestQZ;
					nodeExtended = bestExtended;
					nodeDepth = bestDepth;
					++fit.stalledRays;
				}
			}

			double const magnitude = std::sqrt( qR*qR + qZ*qZ );
			if ( !( magnitude > 0.0 ) )
				throw std::runtime_error(
					"ContourTracer::fitByAngle: the flux vanishes on the surface" );

			double const tangentR = -qZ/magnitude;
			double const tangentZ = qR/magnitude;

			// rho' = rho ( u . t ) / ( u' . t ), with u' . t = u x t. See the
			// header for the derivation and for the two closed forms it was
			// checked against.
			double const uDotT = cosine*tangentR + sine*tangentZ;
			double const cross = cosine*tangentZ - sine*tangentR;

			fit.radius[ j ] = guess;
			fit.pointR[ j ] = axis.r + guess*cosine;
			fit.pointZ[ j ] = axis.z + guess*sine;
			fit.fluxR[ j ] = qR;
			fit.fluxZ[ j ] = qZ;
			fit.crossing[ j ] = std::abs( cross );
			fit.transversality = std::min( fit.transversality, std::abs( cross ) );

			fit.extended[ j ] = nodeExtended ? 1 : 0;
			if ( nodeExtended )
			{
				++fit.extendedNodes;
				fit.deepestBandNode = std::max( fit.deepestBandNode, nodeDepth );
			}

			fit.radiusPrime[ j ] = std::abs( cross ) > 0.0 ? guess*uDotT/cross : 0.0;
			fit.speed[ j ] = std::sqrt( fit.radiusPrime[ j ]*fit.radiusPrime[ j ]
			                            + guess*guess );

			fit.worstResidual = std::max( fit.worstResidual,
			                              std::abs( psi - contour.level ) );
			fit.worstIterations = std::max( fit.worstIterations, iterations );
			fit.totalIterations += iterations;
		}

		fit.fallbackLocations = fallbacks;
		fit.transverse = fit.transversality >= floor;
		if ( !fit.transverse )
		{
			std::ostringstream message;
			message << "ContourTracer::fitByAngle: min | u x t | = "
			        << fit.transversality << " is below the floor " << floor
			        << ", so a ray is tangent to the surface and rho' is not "
			        << "defined there. This is the star-shapedness hypothesis "
			        << "degenerating, and a number returned from it would be "
			        << "worse than no number";
			throw std::runtime_error( message.str() );
		}

		return fit;
	}

}
