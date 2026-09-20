#include "CriticalPoints.hpp"

#include "ConductorField.hpp"
#include "Threading.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

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
		// The solver knows whether the split is in use, so a finder built from
		// one never has to be told separately -- which is what stops a caller
		// searching the remainder for a physical axis by forgetting a line.
		setConductorField( solverIn.conductorField() );
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
			throw std::invalid_argument( "CriticalPointFinder: MEQ is two dimensional" );

		if ( potentialIn.FESpace()->GetMesh() != &meshRef )
			throw std::invalid_argument(
				"CriticalPointFinder: the flux and the potential are on different meshes" );

		/*
		 * THIS CLASS IS HOST CODE THROUGHOUT, SO THE FIELDS COME TO THE HOST
		 * ONCE, HERE.
		 *
		 * fluxScale() and the potential sweeps read these through operator(),
		 * which is RAW -- it neither syncs nor invalidates -- and every other
		 * reader goes through GetValue()/GetVectorValue(), which evaluate on
		 * the host as well. With an mfem::Device configured the fields arrive
		 * from a solve that left them device-resident, so those raw reads take
		 * a stale host copy.
		 *
		 * THE CONSTRUCTOR IS THE RIGHT PLACE because the finder does not own
		 * the fields and never writes them: one sync at the point they are
		 * taken covers every method, where a sync per reader is a list that has
		 * to be kept correct as methods are added. Found by
		 * mfem::Device( "debug" ) faulting in fluxScale(), reached from
		 * checkAxis() -- the driver's post-solve axis check, which is to say on
		 * a path that had already produced the right answer and was about to
		 * report on it.
		 *
		 * No-ops with no Device configured.
		 */
		fluxIn.HostRead();
		potentialIn.HostRead();
	}

	void CriticalPointFinder::setExcluded(
		std::function< bool( double, double ) > excludedIn )
	{
		excluded = std::move( excludedIn );
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

	void CriticalPointFinder::setSeedRings( int ringsIn )
	{
		if ( ringsIn < 0 )
			throw std::invalid_argument(
				"CriticalPointFinder::setSeedRings: the number of rings cannot be "
				"negative" );
		seedRings = ringsIn;
	}

	bool CriticalPointFinder::senseAccepts( AxisSense sense,
	                                        CriticalPointType type )
	{
		// Degenerate is refused by every sense, including Saddle. A determinant
		// at round-off says the classification is not entitled, not that the
		// point is indefinite -- see CriticalPointType's own comment -- and a
		// caller following an X-point through a continuation would take one as
		// the X-point having moved rather than as a point that could not be
		// classified.
		switch ( sense )
		{
			case AxisSense::Either:
				return type == CriticalPointType::Maximum
				       || type == CriticalPointType::Minimum;
			case AxisSense::Maximum:
				return type == CriticalPointType::Maximum;
			case AxisSense::Minimum:
				return type == CriticalPointType::Minimum;
			case AxisSense::Saddle:
				break;
		}
		return type == CriticalPointType::Saddle;
	}

	long CriticalPointFinder::newtonSolves() const
	{
		return newtonSolveCount;
	}

	long CriticalPointFinder::elementsRooted() const
	{
		return elementCount;
	}

	void CriticalPointFinder::resetCounters()
	{
		newtonSolveCount = 0;
		elementCount = 0;
	}

	void CriticalPointFinder::setContainment( double containmentIn )
	{
		if ( containmentIn < 0.0 || containmentIn > 0.5 )
			throw std::invalid_argument(
				"CriticalPointFinder: containment must be in [ 0, 0.5 ]" );
		containment = containmentIn;
	}

	void CriticalPointFinder::setConductorField( ConductorField const *c )
	{
		conductors = c;
	}

	ConductorField const *CriticalPointFinder::conductorField() const
	{
		return conductors;
	}

	namespace
	{
		/*
		 * The physical point of a reference point, through the REENTRANT
		 * transformation overload. CLAUDE.md records Mesh::GetElementTransformation( int )
		 * handing out shared scratch, and this file already paid for that once:
		 * its own site was among the six fixed into function-local thread_locals.
		 */
		void pointOf( mfem::Mesh &mesh, int element,
		              mfem::IntegrationPoint const &ip, double &r, double &z )
		{
			thread_local mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( element, &transformation );
			transformation.SetIntPoint( &ip );

			double coordinates[ 3 ] = { 0.0, 0.0, 0.0 };
			mfem::Vector position( coordinates, 3 );
			transformation.Transform( ip, position );
			r = position( 0 );
			z = position( 1 );
		}
	}

	void CriticalPointFinder::totalFlux( int element,
	                                     mfem::IntegrationPoint const &ip,
	                                     mfem::Vector &out ) const
	{
		fluxField.GetVectorValue( element, ip, out );
		if ( !conductors )
			return;

		double r = 0.0;
		double z = 0.0;
		pointOf( meshRef, element, ip, r, z );

		double qR = 0.0;
		double qZ = 0.0;
		conductors->flux( r, z, qR, qZ );
		out( 0 ) += qR;
		out( 1 ) += qZ;
	}

	double CriticalPointFinder::totalPotential(
		int element, mfem::IntegrationPoint const &ip ) const
	{
		double const solved = potentialField.GetValue( element, ip );
		if ( !conductors )
			return solved;

		double r = 0.0;
		double z = 0.0;
		pointOf( meshRef, element, ip, r, z );
		return solved + conductors->psi( r, z );
	}

	double CriticalPointFinder::nodeShift( int element, int localDof ) const
	{
		if ( !conductors )
			return 0.0;

		// A coefficient IS the value at its node here: meq's spaces are
		// BasisType::GaussLobatto, so this shift is exact rather than a
		// convenient approximation. On a non-nodal basis it would still be the
		// right O( 1 ) correction for a SCREEN, which is all these callers are.
		mfem::FiniteElement const *element_fe
			= potentialField.FESpace()->GetFE( element );
		mfem::IntegrationPoint const &ip
			= element_fe->GetNodes().IntPoint( localDof );

		double r = 0.0;
		double z = 0.0;
		pointOf( meshRef, element, ip, r, z );
		return conductors->psi( r, z );
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

		/*
		 * AND THE SCALE MUST BE THE PHYSICAL FIELD'S UNDER THE SPLIT, WHICH IS
		 * NOT A TIDINESS POINT. The residual this scale sets a target for is
		 * | q_p + q_c |, so a target of tolerance * max | q_p | is measured
		 * against the wrong field -- and where a conductor dominates, q_p is
		 * small while the residual is not, giving a target the Newton cannot
		 * reach and a search that reports failure on a perfectly good axis. In
		 * the vacuum limit q_p is identically zero and the fallback below would
		 * make the rule absolute at `tolerance`, which is the same defect in
		 * its most extreme form.
		 *
		 * One extra pass over the nodes, only when the split is in use.
		 */
		if ( conductors )
		{
			mfem::FiniteElementSpace const &scalarSpace = *potentialField.FESpace();
			for ( int element = 0; element < meshRef.GetNE(); ++element )
			{
				mfem::FiniteElement const *element_fe
					= scalarSpace.GetFE( element );
				mfem::IntegrationRule const &nodes = element_fe->GetNodes();

				for ( int i = 0; i < nodes.GetNPoints(); ++i )
				{
					double r = 0.0;
					double z = 0.0;
					pointOf( meshRef, element, nodes.IntPoint( i ), r, z );
					if ( !( r > 0.0 ) )
						continue;

					double qR = 0.0;
					double qZ = 0.0;
					conductors->flux( r, z, qR, qZ );
					worst = std::max( worst,
					                  std::sqrt( qR*qR + qZ*qZ ) );
				}
			}
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

			totalFlux( element, referencePoint( hi[ 0 ], hi[ 1 ] ), high );
			totalFlux( element, referencePoint( lo[ 0 ], lo[ 1 ] ), low );

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

		// THE UNIT OF WORK, COUNTED AT THE TOP AND NOT AT THE BOTTOM. An attempt
		// that diverges or lands outside its element costs the same iterations
		// as one that succeeds, so counting only the accepted roots would make a
		// search that fails everywhere look free. newtonSolves() is what the
		// seeded entry points are measured against a sweep by.
		//
		// ATOMIC, because sweep() calls this from an OpenMP region. A count is
		// order independent, so the total is exactly the serial one -- which is
		// why this is an atomic increment rather than a per-thread partial: the
		// number is the same and there is nothing to reassociate.
		MEQ_OMP( atomic )
		++newtonSolveCount;

		mfem::IntegrationPoint ip = seed;
		mfem::Vector value( 2 );

		bool converged = false;
		for ( int iteration = 0; iteration < maxIterations; ++iteration )
		{
			totalFlux( element, ip, value );
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
				totalFlux( element, ip, value );
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

		// THREAD LOCAL, NOT THE MESH'S SHARED SCRATCH.
		// Mesh::GetElementTransformation( int ) returns a pointer to one
		// IsoparametricTransformation owned by the Mesh -- "calling this
		// function resets pointers obtained from previous calls" -- so two
		// threads refining two candidates would silently be reading each
		// other's element, and the failure is a Jacobian and a position taken
		// from the wrong map rather than a crash. The
		// ( i, IsoparametricTransformation * ) overload writes only into what
		// it is given. A function static rather than a local so that a sweep
		// over every element does not construct a DenseMatrix per candidate;
		// see the same pattern in ContourTracer::locate().
		thread_local mfem::IsoparametricTransformation scratch;
		meshRef.GetElementTransformation( element, &scratch );
		mfem::ElementTransformation *trans = &scratch;
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
		found.psi = totalPotential( element, ip );
		found.element = element;
		found.fluxResidual = std::sqrt( value( 0 )*value( 0 ) + value( 1 )*value( 1 ) );
		found.determinant = det;
		found.trace = tr;
		found.overshoot = overshoot;

		// The Newton iterate itself, in reference coordinates. See
		// CriticalPoint::referenceX for why this is reported rather than left to
		// be recovered by inverting the element map.
		found.referenceX = ip.x;
		found.referenceY = ip.y;

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
			totalFlux( element, nodes[ i ], value );
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
		// BOTH ENDS, because MEQ's psi is not sign-normalised: with F single
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
				// SCREENING ON THE TOTAL, not the remainder: a nearby
				// conductor's psi_c can dominate psi_p, and seeding in the
				// wrong basin finds a different equilibrium's axis.
				double const value = potentialField( dofs[ i ] )
				                     + nodeShift( element, i );
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
	CriticalPointFinder::pointsFrom( std::vector<int> const &elements,
	                                 AxisSense sense ) const
	{
		std::vector<CriticalPoint> points;
		std::vector<mfem::IntegrationPoint> seeds;
		double const target = tolerance*fluxScale();

		for ( std::size_t e = 0; e < elements.size(); ++e )
		{
			int const element = elements[ e ];
			++elementCount;
			elementSeeds( element, seeds );

			for ( std::size_t i = 0; i < seeds.size(); ++i )
			{
				CriticalPoint point;
				if ( !rootInElement( element, seeds[ i ], target, point ) )
					continue;
				if ( !senseAccepts( sense, point.type ) )
					continue;

				// One physical critical point found from two neighbouring elements
				// gives two answers O( h^(k+1) ) apart, because that is how far q_h
				// disagrees with itself across a face. They are the same object, so
				// they are merged whenever they are of the same type and within an
				// element of each other -- a scale on which a seeded search cannot
				// tell two of them apart in any case. A maximum and a saddle are
				// never merged, so a spurious PAIR survives this and is reported,
				// which is what sweep() is for. That the same-type rule is what
				// keeps the pair is why this filters AFTER classifying rather than
				// asking rootInElement() for one type: the merge needs to see both.
				double const reach = meshRef.GetElementSize( element );
				bool duplicate = false;
				for ( std::size_t j = 0; j < points.size(); ++j )
				{
					double const dr = points[ j ].r - point.r;
					double const dz = points[ j ].z - point.z;
					if ( points[ j ].type == point.type
					     && std::sqrt( dr*dr + dz*dz ) < reach )
					{
						duplicate = true;
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

	bool CriticalPointFinder::tryFindAxis( CriticalPoint &found,
	                                       AxisSense sense ) const
	{
		// REFUSED RATHER THAN ANSWERED false, and the header says why: this
		// entry point's seeds are the extreme NODAL values of psi_h, which an
		// X-point is nowhere near, so a false here would report the absence of a
		// saddle on the strength of never having looked for one.
		if ( sense == AxisSense::Saddle )
			throw std::invalid_argument(
				"CriticalPointFinder::tryFindAxis: AxisSense::Saddle. This entry "
				"point seeds from the extreme nodal values of psi_h, which is "
				"where an axis is and is not where an X-point is. Use "
				"tryFindCriticalPointFrom() with a prior, or sweep() without one" );

		std::vector<int> elements;
		axisSeeds( elements );

		std::vector<CriticalPoint> extrema = pointsFrom( elements, sense );

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
				if ( senseAccepts( sense, all[ i ].type ) )
					extrema.push_back( all[ i ] );

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

	int CriticalPointFinder::nearestElementCentre( double r, double z ) const
	{
		int best = -1;
		double bestDistance = std::numeric_limits<double>::infinity();

		mfem::Vector centre;
		for ( int element = 0; element < meshRef.GetNE(); ++element )
		{
			meshRef.GetElementCenter( element, centre );
			double const dr = centre( 0 ) - r;
			double const dz = centre( 1 ) - z;
			double const distance = dr*dr + dz*dz;
			if ( distance < bestDistance )
			{
				bestDistance = distance;
				best = element;
			}
		}
		return best;
	}

	bool CriticalPointFinder::tryFindAxisFrom( double r, double z,
	                                          AxisSense sense,
	                                          CriticalPoint &found ) const
	{
		// The whole of this function, and deliberately so: the search below
		// serves a saddle perfectly well, and what is withheld here is only the
		// NAME. See the header on why the axis keeps a name of its own over a
		// search that finds more than one kind of point.
		if ( sense == AxisSense::Saddle )
			throw std::invalid_argument(
				"CriticalPointFinder::tryFindAxisFrom: AxisSense::Saddle is not an "
				"axis. The search itself serves it -- call "
				"tryFindCriticalPointFrom(), which is this function without the "
				"refusal" );

		return tryFindCriticalPointFrom( r, z, sense, found );
	}

	bool CriticalPointFinder::tryFindCriticalPointFrom( double r, double z,
	                                                   AxisSense sense,
	                                                   CriticalPoint &found ) const
	{
		int const seed = nearestElementCentre( r, z );
		if ( seed < 0 )
			return false;

		/*
		 * RING BY RING, STOPPING AS SOON AS SOMETHING IS FOUND, WHICH IS WHAT
		 * MAKES A GENEROUS CAP AFFORDABLE.
		 *
		 * The cost of this entry point is then set by how far the axis actually
		 * moved rather than by seedRings: a warm start whose axis is still in
		 * the seed element roots ONE element and stops, and only a seed that has
		 * fallen behind pays for the outer rings. A fixed count would charge the
		 * worst case every time, and the worst case is the rare one.
		 *
		 * AND THE CAP HAS TO BE GENEROUS BECAUSE A RING IS NOT A DISTANCE.
		 * Face-neighbour hops on a diagonally split Cartesian mesh advance about
		 * a quarter of an element each: measured on the Solov'ev benchmark at
		 * k = 1, n = 16, a seed HALF an element away diagonally needs THREE
		 * rings, and two -- axisSeeds()' own choice -- returns nothing at all.
		 * That is not a defect in the seed; axisSeeds() grows its two rings from
		 * the element holding the extreme NODAL value, which is already within a
		 * node of the answer, and a general seed is not.
		 */
		mfem::Table const &neighbours = meshRef.ElementToElementTable();
		std::vector<bool> chosen( static_cast<std::size_t>( meshRef.GetNE() ),
		                          false );
		chosen[ static_cast<std::size_t>( seed ) ] = true;

		std::vector<int> frontier( 1, seed );
		std::vector<CriticalPoint> points;

		for ( int ring = 0; ring <= seedRings && !frontier.empty(); ++ring )
		{
			std::vector<CriticalPoint> const reached =
				pointsFrom( frontier, sense );
			points.insert( points.end(), reached.begin(), reached.end() );

			/*
			 * AN OUT-OF-ELEMENT ROOT IS NOT GOOD ENOUGH TO STOP AT, AND THAT IS
			 * MEASURED RATHER THAN CAUTIOUS. tryFindAxis() records the same rule
			 * for its own seeded path; this is the same hazard one entry point
			 * along, and the test for it caught this stopping on the first ring
			 * that yielded anything.
			 *
			 * q_h jumps across a face, so a root lying within that jump belongs
			 * to NEITHER neighbour strictly: each side's polynomial puts its own
			 * version of it a little way into the other's territory, and
			 * setContainment() lets both through so that an axis landing on a
			 * mesh line is found at all. The two versions are O( h^(k+1) )
			 * apart. Measured on the Solov'ev benchmark at n = 8, seeded half an
			 * element away: element 62 offers the root 3.1e-03 / 1.0e-04 /
			 * 2.6e-07 away at k = 1 / 2 / 3 with a non-zero overshoot, and
			 * element 79 -- one ring further out -- holds it properly.
			 *
			 * FOR A POSITION THAT GAP IS HARMLESS. For a caller building a
			 * BORDERED JACOBIAN ROW out of this element's shape functions it is
			 * not: the row would sit on the wrong element's dofs. So a root that
			 * is strictly inside its element stops the search and one that is
			 * merely admissible does not.
			 */
			bool clean = false;
			for ( std::size_t i = 0; i < points.size(); ++i )
				clean = clean || ( points[ i ].overshoot <= 0.0 );
			if ( clean )
				break;

			// GROWN FROM THE FRONTIER RATHER THAN BY RE-SCANNING THE MESH.
			// axisSeeds() re-walks every element per ring, which it can afford
			// because it is already paying a pass over every nodal value; here
			// that would make the growth the dominant cost and there would be no
			// point to the entry point at all.
			std::vector<int> next;
			for ( std::size_t f = 0; f < frontier.size(); ++f )
			{
				int const element = frontier[ f ];
				int const *row = neighbours.GetRow( element );
				for ( int i = 0; i < neighbours.RowSize( element ); ++i )
				{
					int const other = row[ i ];
					if ( other < 0 || chosen[ static_cast<std::size_t>( other ) ] )
						continue;
					chosen[ static_cast<std::size_t>( other ) ] = true;
					next.push_back( other );
				}
			}
			frontier.swap( next );
		}

		if ( points.empty() )
			return false;

		// STRICTLY INSIDE ITS ELEMENT FIRST, THEN NEAREST THE SEED.
		//
		// The first key is the one above: a root its own element actually holds
		// is worth more than a nearer one seen from across a face, because what
		// the caller does with `element` is evaluate its shape functions.
		//
		// The second is the tie-break, and it is where this entry point differs
		// from tryFindAxis(). A caller with a prior is FOLLOWING one critical
		// point, so among equally admissible roots the continuation of the one
		// being followed is the nearest to where it was last seen. tryFindAxis()
		// has no prior and refuses the ambiguity instead.
		std::size_t best = 0;
		double bestDistance = std::numeric_limits<double>::infinity();
		bool bestClean = false;
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			double const dr = points[ i ].r - r;
			double const dz = points[ i ].z - z;
			double const distance = dr*dr + dz*dz;
			bool const isClean = ( points[ i ].overshoot <= 0.0 );

			if ( i == 0 || ( isClean && !bestClean )
			     || ( isClean == bestClean && distance < bestDistance ) )
			{
				bestDistance = distance;
				bestClean = isClean;
				best = i;
			}
		}

		found = points[ best ];
		return true;
	}

	CriticalPoint CriticalPointFinder::findAxis( AxisSense sense ) const
	{
		// Named here as well as in tryFindAxis(), so that the message a caller
		// meets says findAxis() -- which is what they wrote -- rather than
		// naming the helper it delegates to.
		if ( sense == AxisSense::Saddle )
			throw std::invalid_argument(
				"CriticalPointFinder::findAxis: AxisSense::Saddle. An X-point is "
				"not an extreme nodal value of psi_h and is not near one, so this "
				"entry point's seeds cannot reach it. Use "
				"tryFindCriticalPointFrom() with a prior, or sweep() without one" );

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

		double const target = tolerance*fluxScale();
		int const elementTotal = meshRef.GetNE();

		/*
		 * THE ROOT FINDING IS THREADED AND THE DEDUPLICATION IS NOT, AND THAT
		 * SPLIT IS THE WHOLE DESIGN.
		 *
		 * The merge rule below keeps "the one least outside its own element",
		 * and the comment on it records why: on the Solov'ev benchmark at
		 * k = 1, n = 6 the candidate strictly inside its element is 2.7e-3 from
		 * the true axis and the neighbour sitting 8.5e-2 outside is 6.1e-3, so
		 * WHICH candidate survives is the answer rather than a detail. The
		 * serial loop compares each candidate against the list built SO FAR,
		 * and that list is in element order. Deduplicating inside the parallel
		 * region would make it thread-arrival order, and the reported axis
		 * would move at the 1e-3 level with OMP_NUM_THREADS -- quietly, at a
		 * magnitude that reads as mesh noise.
		 *
		 * So each element writes its own candidates, the buffers are walked in
		 * ELEMENT ORDER afterwards, and the existing merge runs serially over
		 * that concatenation. It sees exactly the sequence the serial loop saw,
		 * so the returned vector is identical entry for entry and field for
		 * field at any thread count.
		 *
		 * `reach` travels WITH the candidate: it is the size of the element the
		 * root was found in, so a candidate compared later must be compared
		 * against its own element's reach and not the one being merged into.
		 */
		struct Candidate
		{
			CriticalPoint point;
			double reach = 0.0;
		};
		std::vector<std::vector<Candidate>> perElement(
			static_cast<std::size_t>( elementTotal ) );
		std::exception_ptr failure;

		MEQ_OMP( parallel for schedule( dynamic ) )
		for ( int element = 0; element < elementTotal; ++element )
		{
			// Nothing on this path is documented to throw, but the region is
			// structured and an escaping exception would be undefined
			// behaviour rather than a diagnosable failure. Caught and rethrown
			// after the region, as the border assemblers do.
			try
			{
				// PER THREAD, not the one vector the serial loop reused across
				// elements: elementSeeds() clears and refills it, and a shared
				// std::vector resized from several threads is a reallocation
				// under another thread's iterator.
				std::vector<mfem::IntegrationPoint> seeds;
				elementSeeds( element, seeds );

				/*
				 * AND `Mesh::GetElementSize( int, int )` IS THE MESH'S SHARED
				 * SCRATCH BY ANOTHER NAME, WHICH NOTHING HAD NAMED.
				 *
				 * It is `GetElementSize( GetElementTransformation( i ), type )`
				 * -- the one-argument overload, the one CLAUDE.md records as
				 * "the returned object is owned by the class and is shared" --
				 * and it then calls SetIntPoint() and Jacobian() on it. Two
				 * threads asking two elements for their size would each read
				 * the other's geometry, with no crash and no error. The
				 * const overload taking a transformation is the reentrant
				 * route, and filling a local with the same
				 * GetElementTransformation( i, & ) the shared one uses makes
				 * the answer bit-identical.
				 */
				thread_local mfem::IsoparametricTransformation sizeScratch;
				meshRef.GetElementTransformation( element, &sizeScratch );
				double const reach = 0.5*meshRef.GetElementSize( &sizeScratch );

				for ( std::size_t i = 0; i < seeds.size(); ++i )
				{
					CriticalPoint point;
					if ( !rootInElement( element, seeds[ i ], target, point ) )
						continue;
					perElement[ static_cast<std::size_t>( element ) ].push_back(
						Candidate{ point, reach } );
				}
			}
			catch ( ... )
			{
				MEQ_OMP( critical( meqCriticalPointSweep ) )
				{
					if ( !failure )
						failure = std::current_exception();
				}
			}
		}

		if ( failure )
			std::rethrow_exception( failure );

		// Unconditional in the serial loop and unconditional here: every
		// element is visited whether or not it yields a candidate.
		elementCount += elementTotal;

		std::vector<CriticalPoint> points;

		for ( int element = 0; element < elementTotal; ++element )
		{
			for ( Candidate const &candidate :
			      perElement[ static_cast<std::size_t>( element ) ] )
			{
				CriticalPoint const &point = candidate.point;
				double const reach = candidate.reach;

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

	void CriticalPointFinder::nodalExtreme( bool wantMaximum, double &value,
	                                        int &element, double &r,
	                                        double &z ) const
	{
		// The same loop GradShafranovSolver runs to produce psi_ax: every element,
		// every dof, strict comparison so a tie goes to the first seen. Written
		// out rather than shared because the solver's version reads a block of the
		// recovery scratch mid-Newton and this one reads a GridFunction after it.
		mfem::FiniteElementSpace const *space = potentialField.FESpace();
		double best = wantMaximum ? -std::numeric_limits<double>::infinity()
		                          : std::numeric_limits<double>::infinity();
		int bestElement = -1;
		int bestLocal = -1;

		mfem::Array<int> dofs;
		for ( int e = 0; e < meshRef.GetNE(); ++e )
		{
			space->GetElementDofs( e, dofs );
			for ( int i = 0; i < dofs.Size(); ++i )
			{
				int const index = dofs[ i ] >= 0 ? dofs[ i ] : -1 - dofs[ i ];
				// The same screen and the same reason; see nodeShift().
				double const here = potentialField( index )
				                    + nodeShift( e, i );
				if ( wantMaximum ? here > best : here < best )
				{
					best = here;
					bestElement = e;
					bestLocal = i;
				}
			}
		}

		value = best;
		element = bestElement;
		r = 0.0;
		z = 0.0;
		if ( bestElement < 0 )
			return;

		// The node's own position, from the element's nodal IntegrationRule. Every
		// space MEQ builds for the potential is a nodal L2 collection, so the rule
		// has one point per dof; a basis where it does not gets the element centre
		// rather than a throw, because this is a diagnostic and a basis choice is
		// not the thing it is guarding.
		mfem::FiniteElement const *fe = space->GetFE( bestElement );
		mfem::IntegrationRule const &nodes = fe->GetNodes();

		// THREAD LOCAL, NOT THE MESH'S SHARED SCRATCH -- see rootInElement().
		thread_local mfem::IsoparametricTransformation scratch;
		meshRef.GetElementTransformation( bestElement, &scratch );

		mfem::IntegrationPoint node;
		if ( nodes.GetNPoints() == fe->GetDof() )
			node = nodes.IntPoint( bestLocal );
		else
			node = mfem::Geometries.GetCenter(
				meshRef.GetElementBaseGeometry( bestElement ) );

		mfem::Vector physical( 2 );
		scratch.Transform( node, physical );
		r = physical( 0 );
		z = physical( 1 );
	}

	AxisAgreement CriticalPointFinder::checkAxis( double psiAxisIn,
	                                             double psiBoundaryIn,
	                                             double toleranceIn ) const
	{
		double const span = psiAxisIn - psiBoundaryIn;
		if ( !( std::abs( span ) > 0.0 ) )
			throw std::invalid_argument(
				"CriticalPointFinder::checkAxis: psi_ax equals psi_bnd, so the "
				"normalised flux is undefined and there is nothing to compare" );

		AxisAgreement result;
		result.psiAxis = psiAxisIn;
		result.psiBoundary = psiBoundaryIn;

		// THE SENSE FOLLOWS THE SPAN AND IS NOT GUESSED. The plasma is where
		// ( psi - psi_bnd ) carries the span's sign -- NormalisedSource's own
		// insidePlasma() -- so a positive span puts the axis at a maximum. That is
		// what lets this dodge AxisSense::Either's refusal, which is the right
		// behaviour for a caller who does not know the sign of F and the wrong one
		// here, where psi_ax itself says which way round the plasma is.
		CriticalPointType const wanted = span > 0.0 ? CriticalPointType::Maximum
		                                            : CriticalPointType::Minimum;

		nodalExtreme( span > 0.0, result.nodalExtreme, result.nodeElement,
		              result.nodeR, result.nodeZ );

		// A FULL SWEEP, NOT tryFindAxis(). Two reasons, and both are about the
		// failure this exists to catch. tryFindAxis() seeds from the extreme nodal
		// values -- the very quantity under suspicion -- and it returns false
		// outright wherever more than one extremum is reachable, which on a real
		// machine case with a spurious ridge near the axis is the ordinary
		// outcome. A sweep costs two Newtons per element and is negligible beside
		// the solve that produced the field.
		std::vector<CriticalPoint> const all = sweep();
		for ( std::size_t i = 0; i < all.size(); ++i )
		{
			if ( all[ i ].type == CriticalPointType::Saddle )
				++result.saddles;
			if ( all[ i ].type != wanted )
				continue;

			// WHERE AN AXIS CANNOT BE. setExcluded() says at length what this
			// is for; the short version is that a conductor's O-point is an
			// extremum of psi and would otherwise win on flux. It is NOT
			// counted in extrema either: the count is of candidates, and a
			// point that cannot be an axis is not one.
			if ( excluded && excluded( all[ i ].r, all[ i ].z ) )
				continue;

			++result.extrema;
			double const flux = ( all[ i ].psi - psiBoundaryIn )/span;
			// THE LARGEST Psi WINS, deliberately: a spurious extremum then costs a
			// missed detection rather than a false alarm, which is the right way
			// round for a warning. AxisAgreement's comment carries the argument.
			if ( !result.located || flux > result.normalisedFlux )
			{
				result.located = true;
				result.axis = all[ i ];
				result.normalisedFlux = flux;
			}
		}

		if ( !result.located )
			return result;

		double const dr = result.axis.r - result.nodeR;
		double const dz = result.axis.z - result.nodeZ;
		result.separation = std::sqrt( dr*dr + dz*dz );

		double const size = result.axis.element >= 0
			? meshRef.GetElementSize( result.axis.element ) : 0.0;
		result.separationInElements = size > 0.0 ? result.separation/size : 0.0;

		result.agrees = result.normalisedFlux >= 1.0 - toleranceIn;
		return result;
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
					totalFlux( face->Elem1No, elementIp, value );

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
