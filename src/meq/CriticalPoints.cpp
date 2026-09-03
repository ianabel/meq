#include "CriticalPoints.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

/*
 * The implementation of INVERSION-PLAN.md stage IN-A. CriticalPoints.hpp
 * carries the argument; this file carries the arithmetic, and comments here are
 * confined to the places where the arithmetic is not the obvious transcription
 * of it.
 */

namespace
{

	double const twoPi = 6.283185307179586476925286766559;

	/// A 2x2 solve, by Cramer. Returns false on a singular matrix, which for
	/// Newton on q = 0 means the point is degenerate and there is no step to
	/// take -- not an error, just a place to stop.
	bool solveTwoByTwo( double const matrix[ 2 ][ 2 ], double const rhs[ 2 ],
	                    double solution[ 2 ] )
	{
		double const det = matrix[ 0 ][ 0 ]*matrix[ 1 ][ 1 ]
		                   - matrix[ 0 ][ 1 ]*matrix[ 1 ][ 0 ];
		if ( !( std::abs( det ) > 0.0 ) )
			return false;

		solution[ 0 ] = (  matrix[ 1 ][ 1 ]*rhs[ 0 ] - matrix[ 0 ][ 1 ]*rhs[ 1 ] )/det;
		solution[ 1 ] = ( -matrix[ 1 ][ 0 ]*rhs[ 0 ] + matrix[ 0 ][ 0 ]*rhs[ 1 ] )/det;
		return true;
	}

	/// The signed angular increment from @a from to @a to, folded into
	/// ( -pi, pi ]. This is the branch choice that makes a turning number a
	/// turning number: it assumes the walk resolves the rotation, and
	/// IndexAudit::worstTurn is what says whether that assumption held.
	double angleStep( double from, double to )
	{
		double delta = to - from;
		while ( delta > M_PI )
			delta -= twoPi;
		while ( delta <= -M_PI )
			delta += twoPi;
		return delta;
	}

	mfem::IntegrationPoint referencePoint( double x, double y )
	{
		mfem::IntegrationPoint ip;
		ip.Init( 0 );
		ip.Set2( x, y );
		return ip;
	}

}

namespace meq
{

	char const *criticalPointName( CriticalPointType type )
	{
		switch ( type )
		{
			case CriticalPointType::Maximum:
				return "maximum";
			case CriticalPointType::Minimum:
				return "minimum";
			case CriticalPointType::Saddle:
				return "saddle";
			case CriticalPointType::Degenerate:
				break;
		}
		return "degenerate";
	}

	int criticalPointIndex( CriticalPointType type )
	{
		switch ( type )
		{
			case CriticalPointType::Maximum:
			case CriticalPointType::Minimum:
				return 1;
			case CriticalPointType::Saddle:
				return -1;
			case CriticalPointType::Degenerate:
				break;
		}
		return 0;
	}

	int eulerCharacteristic( mfem::Mesh &mesh )
	{
		// V - E + F. GetNEdges() is only populated once the edge table has been
		// built, which GetNEdges() itself does not do -- Mesh::GetNEdges returns
		// the cached count, and on a freshly constructed 2D mesh that count is
		// already set by Mesh::FinalizeTopology. Asking for the edge-to-vertex
		// table is the cheap way to be sure rather than to assume.
		mesh.GetEdgeVertexTable();
		return mesh.GetNV() - mesh.GetNEdges() + mesh.GetNE();
	}

	CriticalPointFinder::CriticalPointFinder( GradShafranovSolver const &solverIn )
		: CriticalPointFinder( solverIn.flux(), solverIn.potential() )
	{
	}

	CriticalPointFinder::CriticalPointFinder( mfem::GridFunction const &fluxIn,
	                                          mfem::GridFunction const &potentialIn )
		: fluxField( fluxIn ),
		  potentialField( potentialIn ),
		  meshRef( *fluxIn.FESpace()->GetMesh() )
	{
		if ( fluxIn.FESpace()->GetVDim() != 2 )
			throw std::invalid_argument( "CriticalPointFinder: the flux must have vdim 2" );

		if ( meshRef.Dimension() != 2 )
			throw std::invalid_argument( "CriticalPointFinder: meq is two dimensional" );

		if ( potentialIn.FESpace()->GetMesh() != &meshRef )
			throw std::invalid_argument(
				"CriticalPointFinder: the flux and the potential are on different meshes" );
	}

	void CriticalPointFinder::setTolerance( double toleranceIn )
	{
		if ( !( toleranceIn > 0.0 ) )
			throw std::invalid_argument( "CriticalPointFinder: tolerance must be positive" );
		tolerance = toleranceIn;
	}

	void CriticalPointFinder::setMaxIterations( int maxIterationsIn )
	{
		if ( maxIterationsIn < 1 )
			throw std::invalid_argument( "CriticalPointFinder: need at least one iteration" );
		maxIterations = maxIterationsIn;
	}

	void CriticalPointFinder::setJacobianStep( double stepIn )
	{
		if ( !( stepIn > 0.0 ) || stepIn > 0.25 )
			throw std::invalid_argument(
				"CriticalPointFinder: the Jacobian step must be in ( 0, 0.25 ]" );
		jacobianStep = stepIn;
	}

	void CriticalPointFinder::setBoundarySamples( int samplesIn )
	{
		if ( samplesIn < 1 )
			throw std::invalid_argument( "CriticalPointFinder: need at least one sample" );
		boundarySamples = samplesIn;
	}

	void CriticalPointFinder::setSeparation( double separationIn )
	{
		if ( !( separationIn > 0.0 ) )
			throw std::invalid_argument( "CriticalPointFinder: separation must be positive" );
		separation = separationIn;
	}

	void CriticalPointFinder::setContainment( double containmentIn )
	{
		if ( containmentIn < 0.0 || containmentIn > 0.5 )
			throw std::invalid_argument(
				"CriticalPointFinder: containment must be in [ 0, 0.5 ]" );
		containment = containmentIn;
	}

	double CriticalPointFinder::fluxScale() const
	{
		// The largest | q | over the flux dofs. For a nodal basis those are
		// nodal values, so this is a genuine sup over the nodes; for any other
		// basis it is still an O( 1 ) scale, which is all a relative stopping
		// rule needs.
		mfem::FiniteElementSpace const &space = *fluxField.FESpace();
		int const scalarDofs = space.GetNDofs();

		double worst = 0.0;
		for ( int i = 0; i < scalarDofs; ++i )
		{
			double const a = fluxField( space.DofToVDof( i, 0 ) );
			double const b = fluxField( space.DofToVDof( i, 1 ) );
			worst = std::max( worst, std::sqrt( a*a + b*b ) );
		}
		return worst > 0.0 ? worst : 1.0;
	}

	void CriticalPointFinder::referenceJacobian( int element,
	                                             mfem::IntegrationPoint const &ip,
	                                             double jacobian[ 2 ][ 2 ] ) const
	{
		mfem::Geometry::Type const geom = meshRef.GetElementBaseGeometry( element );

		double const base[ 2 ] = { ip.x, ip.y };
		mfem::Vector high( 2 );
		mfem::Vector low( 2 );

		for ( int d = 0; d < 2; ++d )
		{
			// A central difference about the point wanted, unless that would put
			// a sample outside the reference element -- in which case the whole
			// stencil slides until both ends are inside. Sliding changes the
			// point the derivative belongs to by at most one step, which at the
			// default 1e-4 is far below anything that reaches the answer; see the
			// header on why the Jacobian's accuracy does not reach it at all.
			double shift = 0.0;
			for ( int attempt = 0; attempt < 3; ++attempt )
			{
				double hi[ 2 ] = { base[ 0 ], base[ 1 ] };
				double lo[ 2 ] = { base[ 0 ], base[ 1 ] };
				hi[ d ] += shift + jacobianStep;
				lo[ d ] += shift - jacobianStep;

				bool const hiIn = mfem::Geometry::CheckPoint(
					geom, referencePoint( hi[ 0 ], hi[ 1 ] ), 0.0 );
				bool const loIn = mfem::Geometry::CheckPoint(
					geom, referencePoint( lo[ 0 ], lo[ 1 ] ), 0.0 );

				if ( hiIn && loIn )
					break;
				if ( !hiIn )
					shift -= jacobianStep;
				else
					shift += jacobianStep;
			}

			double hi[ 2 ] = { base[ 0 ], base[ 1 ] };
			double lo[ 2 ] = { base[ 0 ], base[ 1 ] };
			hi[ d ] += shift + jacobianStep;
			lo[ d ] += shift - jacobianStep;

			fluxField.GetVectorValue( element, referencePoint( hi[ 0 ], hi[ 1 ] ), high );
			fluxField.GetVectorValue( element, referencePoint( lo[ 0 ], lo[ 1 ] ), low );

			jacobian[ 0 ][ d ] = ( high( 0 ) - low( 0 ) )/( 2.0*jacobianStep );
			jacobian[ 1 ][ d ] = ( high( 1 ) - low( 1 ) )/( 2.0*jacobianStep );
		}
	}

	bool CriticalPointFinder::rootInElement( int element,
	                                         mfem::IntegrationPoint const &seed,
	                                         double target,
	                                         CriticalPoint &found ) const
	{
		mfem::Geometry::Type const geom = meshRef.GetElementBaseGeometry( element );

		mfem::IntegrationPoint ip = seed;
		mfem::Vector value( 2 );

		bool converged = false;
		for ( int iteration = 0; iteration < maxIterations; ++iteration )
		{
			fluxField.GetVectorValue( element, ip, value );
			double const residual = std::sqrt( value( 0 )*value( 0 )
			                                   + value( 1 )*value( 1 ) );
			if ( residual <= target )
			{
				converged = true;
				break;
			}

			double jacobian[ 2 ][ 2 ];
			referenceJacobian( element, ip, jacobian );

			double const rhs[ 2 ] = { -value( 0 ), -value( 1 ) };
			double step[ 2 ];
			if ( !solveTwoByTwo( jacobian, rhs, step ) )
				return false;

			// The reference element has diameter one, so a step longer than that
			// is Newton leaving rather than converging. Truncating it keeps the
			// iterate in a region where the polynomial is the one being rooted;
			// the acceptance test below throws away anything that ends up outside
			// regardless.
			double const length = std::sqrt( step[ 0 ]*step[ 0 ] + step[ 1 ]*step[ 1 ] );
			if ( !std::isfinite( length ) )
				return false;
			if ( length > 1.0 )
			{
				step[ 0 ] /= length;
				step[ 1 ] /= length;
			}

			ip.x += step[ 0 ];
			ip.y += step[ 1 ];

			if ( !mfem::Geometry::CheckPoint( geom, ip, 2.0 ) )
				return false;

			if ( length < 1.0e-15 )
			{
				fluxField.GetVectorValue( element, ip, value );
				converged = true;
				break;
			}
		}

		if ( !converged )
			return false;

		// The root belongs to this element only if it lies in it, or close enough
		// to a face that q_h's own jump there makes the question meaningless. A
		// root of this element's polynomial well outside this element is a root of
		// nothing: psi_h is a different polynomial there. See setContainment() and
		// CriticalPoint::overshoot.
		if ( !mfem::Geometry::CheckPoint( geom, ip, containment + 1.0e-9 ) )
			return false;

		double overshoot = 0.0;
		if ( !mfem::Geometry::CheckPoint( geom, ip, 0.0 ) )
		{
			// How far outside, measured the way CheckPoint measures inside: the
			// worst violated face of the reference element.
			if ( geom == mfem::Geometry::TRIANGLE )
			{
				overshoot = std::max( std::max( -ip.x, -ip.y ), ip.x + ip.y - 1.0 );
			}
			else
			{
				overshoot = std::max( std::max( -ip.x, ip.x - 1.0 ),
				                      std::max( -ip.y, ip.y - 1.0 ) );
			}
			overshoot = std::max( overshoot, 0.0 );
		}

		mfem::ElementTransformation *trans = meshRef.GetElementTransformation( element );
		trans->SetIntPoint( &ip );

		mfem::Vector physical( 2 );
		trans->Transform( ip, physical );

		double reference[ 2 ][ 2 ];
		referenceJacobian( element, ip, reference );

		// dq/dx = ( dq/dxi )( dxi/dx ). MFEM's InverseJacobian() is dxi/dx.
		mfem::DenseMatrix const &inverse = trans->InverseJacobian();
		double physicalJacobian[ 2 ][ 2 ] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
		for ( int i = 0; i < 2; ++i )
			for ( int j = 0; j < 2; ++j )
				for ( int m = 0; m < 2; ++m )
					physicalJacobian[ i ][ j ] += reference[ i ][ m ]*inverse( m, j );

		double const det = physicalJacobian[ 0 ][ 0 ]*physicalJacobian[ 1 ][ 1 ]
		                   - physicalJacobian[ 0 ][ 1 ]*physicalJacobian[ 1 ][ 0 ];
		double const tr = physicalJacobian[ 0 ][ 0 ] + physicalJacobian[ 1 ][ 1 ];

		double norm = 0.0;
		for ( int i = 0; i < 2; ++i )
			for ( int j = 0; j < 2; ++j )
				norm = std::max( norm, std::abs( physicalJacobian[ i ][ j ] ) );

		// Degeneracy is judged against the square of the matrix's own scale,
		// which is the only scale a determinant can be compared with. The
		// threshold is loose on purpose: a point this close to degenerate has no
		// classification worth reporting, and saying so is better than picking
		// one of the three answers.
		double const degenerate = 1.0e-10*norm*norm;

		found.r = physical( 0 );
		found.z = physical( 1 );
		found.psi = potentialField.GetValue( element, ip );
		found.element = element;
		found.fluxResidual = std::sqrt( value( 0 )*value( 0 ) + value( 1 )*value( 1 ) );
		found.determinant = det;
		found.trace = tr;
		found.overshoot = overshoot;

		if ( std::abs( det ) <= degenerate )
			found.type = CriticalPointType::Degenerate;
		else if ( det < 0.0 )
			found.type = CriticalPointType::Saddle;
		else if ( tr > 0.0 )
			found.type = CriticalPointType::Minimum;
		else
			found.type = CriticalPointType::Maximum;

		found.index = criticalPointIndex( found.type );
		return true;
	}

	void CriticalPointFinder::elementSeeds(
		int element, std::vector<mfem::IntegrationPoint> &seeds ) const
	{
		seeds.clear();

		mfem::Geometry::Type const geom = meshRef.GetElementBaseGeometry( element );
		seeds.push_back( mfem::Geometries.GetCenter( geom ) );

		// The node where | q_h | is smallest is the best starting point the
		// element offers for free, and on a coarse mesh it is much better than
		// the centre. Skipped for a basis whose dofs are not nodal values, where
		// GetNodes() is empty and the centre is all there is.
		mfem::FiniteElement const *fe = fluxField.FESpace()->GetFE( element );
		mfem::IntegrationRule const &nodes = fe->GetNodes();
		if ( nodes.Size() != fe->GetDof() )
			return;

		int best = -1;
		double bestValue = std::numeric_limits<double>::infinity();
		mfem::Vector value( 2 );
		for ( int i = 0; i < nodes.Size(); ++i )
		{
			fluxField.GetVectorValue( element, nodes[ i ], value );
			double const magnitude = std::sqrt( value( 0 )*value( 0 )
			                                    + value( 1 )*value( 1 ) );
			if ( magnitude < bestValue )
			{
				bestValue = magnitude;
				best = i;
			}
		}
		if ( best >= 0 )
			seeds.push_back( nodes[ best ] );
	}

	void CriticalPointFinder::axisSeeds( std::vector<int> &elements ) const
	{
		elements.clear();

		mfem::FiniteElementSpace const &space = *potentialField.FESpace();

		// The largest and the smallest NODAL value of psi_h, and which element
		// each sits in. The larger one is what GradShafranovSolver::psiAxis()
		// reports when psi_ax is an unknown -- deliberately, because a nodal value
		// is one entry of the discrete unknown and the bordered Newton needs a
		// differentiable constraint. It is used here only as a place to start
		// looking; the answer this file returns is the critical point, which is a
		// different quantity. See the header.
		//
		// BOTH ENDS, because meq's psi is not sign-normalised: with F single
		// signed negative -- which is what every Solov'ev fixture has -- the axis
		// is an interior MINIMUM and the largest nodal value is a corner of the
		// mesh.
		int maximumElement = -1;
		int minimumElement = -1;
		double largest = -std::numeric_limits<double>::infinity();
		double smallest = std::numeric_limits<double>::infinity();

		mfem::Array<int> dofs;
		for ( int element = 0; element < meshRef.GetNE(); ++element )
		{
			space.GetElementDofs( element, dofs );
			for ( int i = 0; i < dofs.Size(); ++i )
			{
				double const value = potentialField( dofs[ i ] );
				if ( value > largest )
				{
					largest = value;
					maximumElement = element;
				}
				if ( value < smallest )
				{
					smallest = value;
					minimumElement = element;
				}
			}
		}

		if ( maximumElement < 0 || minimumElement < 0 )
			return;

		std::vector<bool> chosen( meshRef.GetNE(), false );
		chosen[ maximumElement ] = true;
		chosen[ minimumElement ] = true;

		// Two rings of face neighbours. See the header for the measurement that
		// says one is not enough.
		mfem::Table const &neighbours = meshRef.ElementToElementTable();
		for ( int ring = 0; ring < 2; ++ring )
		{
			std::vector<bool> grown = chosen;
			for ( int element = 0; element < meshRef.GetNE(); ++element )
			{
				if ( !chosen[ element ] )
					continue;
				int const *row = neighbours.GetRow( element );
				for ( int i = 0; i < neighbours.RowSize( element ); ++i )
					if ( row[ i ] >= 0 )
						grown[ row[ i ] ] = true;
			}
			chosen.swap( grown );
		}

		for ( int element = 0; element < meshRef.GetNE(); ++element )
			if ( chosen[ element ] )
				elements.push_back( element );
	}

	std::vector<CriticalPoint>
	CriticalPointFinder::extremaFrom( std::vector<int> const &elements,
	                                  AxisSense sense ) const
	{
		std::vector<CriticalPoint> extrema;
		std::vector<mfem::IntegrationPoint> seeds;
		double const target = tolerance*fluxScale();

		for ( std::size_t e = 0; e < elements.size(); ++e )
		{
			int const element = elements[ e ];
			elementSeeds( element, seeds );

			for ( std::size_t i = 0; i < seeds.size(); ++i )
			{
				CriticalPoint point;
				if ( !rootInElement( element, seeds[ i ], target, point ) )
					continue;
				if ( point.type != CriticalPointType::Maximum
				     && point.type != CriticalPointType::Minimum )
					continue;
				if ( sense == AxisSense::Maximum
				     && point.type != CriticalPointType::Maximum )
					continue;
				if ( sense == AxisSense::Minimum
				     && point.type != CriticalPointType::Minimum )
					continue;

				// One physical extremum found from two neighbouring elements gives
				// two answers O( h^(k+1) ) apart, because that is how far q_h
				// disagrees with itself across a face. They are the same object, so
				// they are merged whenever they are of the same type and within an
				// element of each other -- a scale on which a seeded search cannot
				// tell two extrema apart in any case. A maximum and a saddle are
				// never merged, so a spurious PAIR survives this and is reported,
				// which is what sweep() is for.
				double const reach = meshRef.GetElementSize( element );
				bool duplicate = false;
				for ( std::size_t j = 0; j < extrema.size(); ++j )
				{
					double const dr = extrema[ j ].r - point.r;
					double const dz = extrema[ j ].z - point.z;
					if ( extrema[ j ].type == point.type
					     && std::sqrt( dr*dr + dz*dz ) < reach )
					{
						duplicate = true;
						if ( point.overshoot < extrema[ j ].overshoot )
							extrema[ j ] = point;
						break;
					}
				}
				if ( !duplicate )
					extrema.push_back( point );
			}
		}

		return extrema;
	}

	bool CriticalPointFinder::tryFindAxis( CriticalPoint &found,
	                                       AxisSense sense ) const
	{
		std::vector<int> elements;
		axisSeeds( elements );

		std::vector<CriticalPoint> extrema = extremaFrom( elements, sense );

		// THE SEEDED SEARCH IS THE FAST PATH AND NOT THE ONLY ONE, because the
		// seed can be wrong and it is cheap to find out. A zero of q_h sitting on
		// a mesh line is credited to whichever element the L2 jump at that node
		// happens to favour, and neither ring reaches the element that actually
		// holds it; measured at k = 1, n = 4 on the Solov'ev benchmark, the seeded
		// search returns nothing at all. Falling back to a full sweep costs one
		// Newton per element -- microseconds beside the solve that produced the
		// field -- and makes the answer independent of a heuristic. The seeded
		// path is kept because it is what runs every other time.
		//
		// AND AN OUT-OF-ELEMENT ROOT IS NOT GOOD ENOUGH TO STOP AT, which is the
		// less obvious half. A root accepted on the containment allowance is an
		// answer from a polynomial evaluated outside its own element, and some
		// OTHER element may hold the same root properly -- a seed set of two rings
		// cannot know. Measured at k = 1, n = 8: the seeded search returns a root
		// 9.2e-2 outside its element and 2.8e-3 from the true axis, while a full
		// sweep finds one strictly inside an element and 6.4e-4 away, four times
		// better. So a non-zero overshoot buys the sweep too.
		if ( extrema.size() != 1 || extrema.front().overshoot > 0.0 )
		{
			CriticalPoint const seeded
				= extrema.size() == 1 ? extrema.front() : CriticalPoint();
			bool const haveSeeded = ( extrema.size() == 1 );

			std::vector<CriticalPoint> const all = sweep();
			extrema.clear();
			for ( std::size_t i = 0; i < all.size(); ++i )
			{
				if ( all[ i ].type != CriticalPointType::Maximum
				     && all[ i ].type != CriticalPointType::Minimum )
					continue;
				if ( sense == AxisSense::Maximum
				     && all[ i ].type != CriticalPointType::Maximum )
					continue;
				if ( sense == AxisSense::Minimum
				     && all[ i ].type != CriticalPointType::Minimum )
					continue;
				extrema.push_back( all[ i ] );
			}

			// The sweep is a superset of the seeded search, so it should never do
			// worse. If it somehow finds nothing where the seeded search found
			// something, keep what was found rather than throwing it away.
			if ( extrema.empty() && haveSeeded )
				extrema.push_back( seeded );
		}

		if ( extrema.size() != 1 )
			return false;

		found = extrema.front();
		return true;
	}

	CriticalPoint CriticalPointFinder::findAxis( AxisSense sense ) const
	{
		CriticalPoint found;
		if ( tryFindAxis( found, sense ) )
			return found;

		// Say which of the two failures it was, because they call for opposite
		// responses: nothing found means the mesh carries no interior extremum of
		// q_h at all, and more than one found means the caller has to say which
		// one they mean.
		int maxima = 0;
		int minima = 0;
		int saddles = 0;
		std::vector<CriticalPoint> const all = sweep();
		for ( std::size_t i = 0; i < all.size(); ++i )
		{
			if ( all[ i ].type == CriticalPointType::Maximum )
				++maxima;
			if ( all[ i ].type == CriticalPointType::Minimum )
				++minima;
			if ( all[ i ].type == CriticalPointType::Saddle )
				++saddles;
		}

		std::ostringstream message;
		message << "CriticalPointFinder::findAxis: no unique interior extremum. "
		        << "A sweep of the mesh found " << maxima << " maxima, "
		        << minima << " minima and " << saddles << " saddles";
		if ( sense == AxisSense::Either && maxima > 0 && minima > 0 )
			message << "; pass AxisSense::Maximum or AxisSense::Minimum to choose";
		if ( maxima + minima == 0 )
			message << ". A sweep is a seeded search and not an exhaustive one, so "
			           "this is not proof that there is none";
		message << ".";
		throw std::runtime_error( message.str() );
	}

	std::vector<CriticalPoint> CriticalPointFinder::sweep() const
	{
		// The bounding box, purely as the length scale the deduplication
		// tolerance is relative to.
		double diameter = 0.0;
		if ( meshRef.GetNV() > 0 )
		{
			double lo[ 2 ] = { meshRef.GetVertex( 0 )[ 0 ], meshRef.GetVertex( 0 )[ 1 ] };
			double hi[ 2 ] = { lo[ 0 ], lo[ 1 ] };
			for ( int v = 1; v < meshRef.GetNV(); ++v )
			{
				double const *vertex = meshRef.GetVertex( v );
				for ( int d = 0; d < 2; ++d )
				{
					lo[ d ] = std::min( lo[ d ], vertex[ d ] );
					hi[ d ] = std::max( hi[ d ], vertex[ d ] );
				}
			}
			diameter = std::sqrt( ( hi[ 0 ] - lo[ 0 ] )*( hi[ 0 ] - lo[ 0 ] )
			                      + ( hi[ 1 ] - lo[ 1 ] )*( hi[ 1 ] - lo[ 1 ] ) );
		}

		std::vector<CriticalPoint> points;
		std::vector<mfem::IntegrationPoint> seeds;
		double const target = tolerance*fluxScale();

		for ( int element = 0; element < meshRef.GetNE(); ++element )
		{
			elementSeeds( element, seeds );
			double const reach = 0.5*meshRef.GetElementSize( element );

			for ( std::size_t i = 0; i < seeds.size(); ++i )
			{
				CriticalPoint point;
				if ( !rootInElement( element, seeds[ i ], target, point ) )
					continue;

				bool duplicate = false;
				for ( std::size_t j = 0; j < points.size(); ++j )
				{
					double const dr = points[ j ].r - point.r;
					double const dz = points[ j ].z - point.z;
					double const distance = std::sqrt( dr*dr + dz*dz );
					if ( distance < separation*diameter )
					{
						duplicate = true;
						break;
					}
					// Same type and within half an element: one object seen from
					// both sides of a face. Different types are never merged, so a
					// spurious maximum-and-saddle pair -- which is the only shape
					// numerical noise produces, and the one thing audit() is blind
					// to -- is reported rather than tidied away.
					if ( points[ j ].type == point.type && distance < reach )
					{
						duplicate = true;
						// KEEP THE ONE LEAST OUTSIDE ITS OWN ELEMENT, which without
						// this rule is decided by element numbering. Measured on the
						// Solov'ev benchmark at k = 1, n = 6: the candidate strictly
						// inside its element is 2.7e-3 from the true axis and the one
						// sitting 8.5e-2 outside a neighbour is 6.1e-3, and the
						// neighbour has the lower element index. Element order is not
						// a tie break.
						if ( point.overshoot < points[ j ].overshoot )
							points[ j ] = point;
						break;
					}
				}
				if ( !duplicate )
					points.push_back( point );
			}
		}

		return points;
	}

	IndexAudit CriticalPointFinder::audit() const
	{
		IndexAudit result;
		result.eulerCharacteristic = eulerCharacteristic( meshRef );

		int const boundaryElements = meshRef.GetNBE();
		if ( boundaryElements == 0 )
			throw std::runtime_error( "CriticalPointFinder::audit: the mesh has no boundary" );

		// --- assemble the boundary into closed loops -------------------------
		//
		// The winding number is an integral along the boundary, so the boundary
		// has to be walked in order and with a consistent orientation. MFEM
		// stores boundary elements in no particular order, so they are threaded
		// here by shared vertex.

		std::vector<std::array<int, 2> > endpoints( boundaryElements );
		std::vector<std::vector<int> > incident( meshRef.GetNV() );

		mfem::Array<int> vertices;
		for ( int b = 0; b < boundaryElements; ++b )
		{
			meshRef.GetBdrElementVertices( b, vertices );
			if ( vertices.Size() != 2 )
				throw std::runtime_error(
					"CriticalPointFinder::audit: a boundary element is not a segment" );
			endpoints[ b ][ 0 ] = vertices[ 0 ];
			endpoints[ b ][ 1 ] = vertices[ 1 ];
			incident[ vertices[ 0 ] ].push_back( b );
			incident[ vertices[ 1 ] ].push_back( b );
		}

		struct Segment
		{
			int boundaryElement;
			bool forward;
		};

		std::vector<std::vector<Segment> > loops;
		std::vector<bool> used( boundaryElements, false );

		for ( int start = 0; start < boundaryElements; ++start )
		{
			if ( used[ start ] )
				continue;

			std::vector<Segment> loop;
			int current = start;
			bool forward = true;
			int const first = endpoints[ start ][ 0 ];

			while ( true )
			{
				used[ current ] = true;
				loop.push_back( Segment{ current, forward } );

				int const last = endpoints[ current ][ forward ? 1 : 0 ];
				if ( last == first )
					break;

				int next = -1;
				bool nextForward = true;
				for ( std::size_t i = 0; i < incident[ last ].size(); ++i )
				{
					int const candidate = incident[ last ][ i ];
					if ( used[ candidate ] )
						continue;
					next = candidate;
					nextForward = ( endpoints[ candidate ][ 0 ] == last );
					break;
				}

				if ( next < 0 )
					throw std::runtime_error(
						"CriticalPointFinder::audit: the boundary does not close into loops" );

				current = next;
				forward = nextForward;
			}

			loops.push_back( loop );
		}

		result.boundaryLoops = static_cast<int>( loops.size() );

		// Orient: the outer loop counter-clockwise, every hole clockwise, which
		// is the orientation of dOmega as the boundary of Omega -- the domain on
		// the left throughout, so the right-hand normal is the outward one.
		std::vector<double> areas( loops.size(), 0.0 );
		for ( std::size_t l = 0; l < loops.size(); ++l )
		{
			double area = 0.0;
			for ( std::size_t s = 0; s < loops[ l ].size(); ++s )
			{
				Segment const &segment = loops[ l ][ s ];
				double const *from = meshRef.GetVertex(
					endpoints[ segment.boundaryElement ][ segment.forward ? 0 : 1 ] );
				double const *to = meshRef.GetVertex(
					endpoints[ segment.boundaryElement ][ segment.forward ? 1 : 0 ] );
				area += from[ 0 ]*to[ 1 ] - to[ 0 ]*from[ 1 ];
			}
			areas[ l ] = 0.5*area;
		}

		std::size_t outer = 0;
		for ( std::size_t l = 1; l < loops.size(); ++l )
			if ( std::abs( areas[ l ] ) > std::abs( areas[ outer ] ) )
				outer = l;

		for ( std::size_t l = 0; l < loops.size(); ++l )
		{
			bool const wantPositive = ( l == outer );
			if ( ( areas[ l ] > 0.0 ) != wantPositive )
			{
				std::reverse( loops[ l ].begin(), loops[ l ].end() );
				for ( std::size_t s = 0; s < loops[ l ].size(); ++s )
					loops[ l ][ s ].forward = !loops[ l ][ s ].forward;
				areas[ l ] = -areas[ l ];
			}
		}

		// --- walk, accumulating the turning of q ------------------------------

		double turning = 0.0;
		double worstTurn = 0.0;
		double smallest = std::numeric_limits<double>::infinity();
		double transversality = std::numeric_limits<double>::infinity();
		int normalSign = 0;
		bool signConsistent = true;

		mfem::Vector value( 2 );
		mfem::Vector physical( 2 );

		for ( std::size_t l = 0; l < loops.size(); ++l )
		{
			bool haveAngle = false;
			double previous = 0.0;
			double firstAngle = 0.0;

			for ( std::size_t s = 0; s < loops[ l ].size(); ++s )
			{
				Segment const &segment = loops[ l ][ s ];
				int const b = segment.boundaryElement;

				mfem::FaceElementTransformations *face
					= meshRef.GetBdrFaceTransformations( b );
				if ( face == nullptr )
					throw std::runtime_error(
						"CriticalPointFinder::audit: no transformation for a boundary face" );

				// Which way the face's own reference coordinate runs relative to
				// the direction this loop traverses it. Asked of the geometry
				// rather than assumed of the vertex ordering, because the two are
				// related by MFEM's face orientation bookkeeping and a wrong guess
				// here reverses a segment silently.
				double const *fromVertex = meshRef.GetVertex(
					endpoints[ b ][ segment.forward ? 0 : 1 ] );

				face->Transform( referencePoint( 0.0, 0.0 ), physical );
				double const atZero = std::abs( physical( 0 ) - fromVertex[ 0 ] )
				                      + std::abs( physical( 1 ) - fromVertex[ 1 ] );
				face->Transform( referencePoint( 1.0, 0.0 ), physical );
				double const atOne = std::abs( physical( 0 ) - fromVertex[ 0 ] )
				                     + std::abs( physical( 1 ) - fromVertex[ 1 ] );
				bool const faceForward = ( atZero <= atOne );

				double const *toVertex = meshRef.GetVertex(
					endpoints[ b ][ segment.forward ? 1 : 0 ] );
				double tangent[ 2 ] = { toVertex[ 0 ] - fromVertex[ 0 ],
				                        toVertex[ 1 ] - fromVertex[ 1 ] };
				double const length = std::sqrt( tangent[ 0 ]*tangent[ 0 ]
				                                 + tangent[ 1 ]*tangent[ 1 ] );
				if ( length > 0.0 )
				{
					tangent[ 0 ] /= length;
					tangent[ 1 ] /= length;
				}
				// The domain is on the left, so the outward normal is the right
				// one.
				double const normal[ 2 ] = { tangent[ 1 ], -tangent[ 0 ] };

				int const firstSample = haveAngle ? 1 : 0;
				for ( int j = firstSample; j <= boundarySamples; ++j )
				{
					double const t = static_cast<double>( j )
					                 /static_cast<double>( boundarySamples );
					double const s2 = faceForward ? t : 1.0 - t;

					mfem::IntegrationPoint faceIp = referencePoint( s2, 0.0 );
					mfem::IntegrationPoint elementIp;
					face->Loc1.Transform( faceIp, elementIp );
					fluxField.GetVectorValue( face->Elem1No, elementIp, value );

					double const magnitude = std::sqrt( value( 0 )*value( 0 )
					                                    + value( 1 )*value( 1 ) );
					smallest = std::min( smallest, magnitude );

					if ( magnitude > 0.0 )
					{
						double const projection = ( value( 0 )*normal[ 0 ]
						                            + value( 1 )*normal[ 1 ] )/magnitude;
						transversality = std::min( transversality, std::abs( projection ) );
						int const sign = projection > 0.0 ? 1 : -1;
						if ( normalSign == 0 )
							normalSign = sign;
						else if ( normalSign != sign )
							signConsistent = false;
					}

					double const angle = std::atan2( value( 1 ), value( 0 ) );
					if ( !haveAngle )
					{
						firstAngle = angle;
						previous = angle;
						haveAngle = true;
					}
					else
					{
						double const delta = angleStep( previous, angle );
						worstTurn = std::max( worstTurn, std::abs( delta ) );
						turning += delta;
						previous = angle;
					}
				}
			}

			// Close the loop: the last sample is the last face's far end, which is
			// the first face's near end, so the two are the same point on the
			// boundary and there is no further increment to take -- except that
			// they are evaluated from different elements, and q_h disagrees with
			// itself by O( h^(k+1) ) across a face. That disagreement is the final
			// increment, and dropping it would leave the turning off by it.
			if ( haveAngle )
			{
				double const delta = angleStep( previous, firstAngle );
				worstTurn = std::max( worstTurn, std::abs( delta ) );
				turning += delta;
			}
		}

		result.turning = turning/twoPi;
		result.windingNumber = static_cast<int>( std::lround( result.turning ) );
		result.windingDefect = std::abs( result.turning
		                                 - static_cast<double>( result.windingNumber ) );
		result.worstTurn = worstTurn;
		result.smallestFlux = ( smallest == std::numeric_limits<double>::infinity() )
		                      ? 0.0 : smallest;
		result.transversality
			= ( transversality == std::numeric_limits<double>::infinity() )
			  ? 0.0 : transversality;
		result.transverse = signConsistent && result.transversality > 0.0;

		return result;
	}

}
