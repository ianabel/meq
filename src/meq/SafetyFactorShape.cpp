#include "meq/SafetyFactorTangent.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "meq/ConductorField.hpp"
#include "meq/FluxSurfaces.hpp"

/*
 * THE q-DRIVEN JACOBIAN BY SHAPE DERIVATIVE (MEASUREMENTS.md M-QS-3).
 *
 * The extraction's quadrature, as fitByAngle() and surfaceAverages() compute
 * it, is
 *
 *     ( q/g )_i = ( 1/N ) sum_j a_j/( | u_j . Q( x_j ) | R_j^2 ),
 *
 * over the N rays at theta_j = 2 pi j/N about the origin O (the located
 * axis). u_j = ( cos theta_j, sin theta_j ), and x_j = O + a_j u_j is the ray's
 * root of psi( x ) = c_i. psi is psi~ plus the conductors' psi_c, and Q is
 * q~ plus q_c. The weight 2 pi R | dx/dtheta | dtheta/( R | Q | ) and
 * < R^-2 > fold into that one sum, because | dx/dtheta | = a | Q |/| u . Q |
 * when the tangent is Q's perpendicular, as the tracer takes it.
 *
 * Along a tangent the fields move by ( dpsi, dQ ), the origin by dO, and the
 * label by dc_i = ( 1 - Psi_N,i ) dpsi_ax + Psi_N,i dpsi_bnd. So, at each
 * node:
 *
 *     da  = ( dc - dpsi( x ) - grad psi . dO )/( grad psi . u ),
 *     dx  = dO + da u,
 *     d( u . Q ) = u . ( dQ( x ) + ( grad Q ) dx ),
 *     df  = f ( da/a - d( u . Q )/( u . Q ) - 2 dR/R ).
 *
 * That is the derivative of the map as implemented, not of its continuous
 * idealisation. The trapezoid rule's own dependence on the origin and its
 * error are both in it, so it matches the differenced extraction up to that
 * difference's truncation and noise. The label is read on psi~ and the metric
 * on q~, as the tracer reads them.
 *
 * WHAT IS GATHERED. Each interior node's values -- psi~, grad psi~, q~,
 * grad q~ and, per tangent, dpsi~ and dq~ -- are element polynomials at a
 * fixed point. Each is one row of a LinearBlocks, one block per node, applied
 * by applyLinearBlocks() or its kernel (applyLinearBlocksOnDevice()), the
 * same loop either way. The conductors' psi_c gradient, q_c and q_c's
 * gradient are added on the host in double, from ConductorField, q_c's
 * gradient centrally differenced as the axis root's already is.
 *
 * A BAND NODE (a curved boundary's, past Gamma_h) has no element
 * polynomial. Its values come from the tracer's own sampleAt(), and the
 * tangent's from a tracer on the tangent fields, since the transfer lift is
 * linear. Its gradients are central differences over cbrt( eps ) of the
 * family's size, so in single they carry the noise the gathered nodes do not.
 */

namespace meq
{

	namespace
	{
		int const baseRows = 9;
		int const tangentRows = 3;

		/// Where a value's rows sit in a base block.
		enum BaseRow
		{
			Psi = 0,
			PsiR,
			PsiZ,
			FluxR,
			FluxZ,
			FluxRR,
			FluxRZ,
			FluxZR,
			FluxZZ
		};

		/// The total field at one node, and what a tangent moves there.
		struct NodeValues
		{
			double psiR = 0.0, psiZ = 0.0;
			double qR = 0.0, qZ = 0.0;
			double qRR = 0.0, qRZ = 0.0, qZR = 0.0, qZZ = 0.0;
		};

		bool gather( LinearBlocks const &blocks, mfem::GridFunction const &potential,
		             mfem::GridFunction const &flux, ContinuousInterpolant::Route route, mfem::Vector &out )
		{
			int const split = potential.Size();
			bool const deviceAllowed = mfem::Device::Allows( mfem::Backend::DEVICE_MASK );
			if ( route != ContinuousInterpolant::Route::Host && deviceAllowed
			     && applyLinearBlocksOnDevice( blocks, potential, &flux, split, out, 0 ) )
				return true;
			if ( route == ContinuousInterpolant::Route::Device )
				throw std::runtime_error( "meq::SafetyFactorEquilibrium::shapeDerivative: the device "
				                          "gather was asked for and cannot run" );
			applyLinearBlocks( blocks, potential, &flux, split, out, 0 );
			return false;
		}
	}

	bool SafetyFactorEquilibrium::shapeDerivative( ContinuousFields const &base, CriticalPoint const &axis,
	                                               double boundary, std::vector<mfem::Vector> const &dx,
	                                               std::vector<double> const &dAxis,
	                                               std::vector<double> const &dBoundary,
	                                               std::vector<double> const &originR,
	                                               std::vector<double> const &originZ,
	                                               std::vector<double> &out )
	{
		std::size_t const m = dx.size();
		FluxSurfaceFamily const family =
			extractFamily( *base.potential, *base.flux, axis, boundary );
		extractionCount += 1;
		if ( family.surfaces.size() != m )
			return false;

		/*
		 * THE NODES. Every surface's, in order; an interior node's element
		 * point found by an unrestricted Newton inverse, as the axis root's is
		 * (a root on a face sits a little outside its element).
		 */
		mfem::FiniteElementSpace const &potentialSpace = *base.potential->FESpace();
		mfem::FiniteElementSpace const &fluxSpace = *base.flux->FESpace();
		mfem::Mesh &mesh = *potentialSpace.GetMesh();
		int const potentialSize = base.potential->Size();

		std::vector<std::size_t> surfaceOf;
		std::vector<double> nodeR, nodeZ, nodeA, cosine, sine;
		std::vector<int> gathered;
		std::vector<int> bandNodes;
		LinearBlocks baseBlocks, tangentBlocks;
		std::vector<int> outOffset{ 0 }, tangentOutOffset{ 0 }, inOffset{ 0 }, matrixOffset{ 0 },
			tangentMatrixOffset{ 0 };
		std::vector<int> outIndex, tangentOutIndex, inIndex;
		std::vector<double> matrix, tangentMatrix;
		double size = 0.0;

		mfem::IsoparametricTransformation transformation;
		mfem::Array<int> potentialDofs, fluxDofs;
		mfem::Vector shape, fluxShape;
		mfem::DenseMatrix dshape, fluxDShape;
		int interior = 0;
		for ( std::size_t i = 0; i < m; ++i )
		{
			FluxSurface const &surface = family.surfaces[ i ];
			std::size_t const count = surface.radius.size();
			if ( count < 4 )
				return false;
			for ( std::size_t j = 0; j < count; ++j )
			{
				double const theta = 2.0*M_PI*static_cast<double>( j )/static_cast<double>( count );
				double const c = std::cos( theta ), sn = std::sin( theta );
				double const radius = surface.radius[ j ], z = surface.z[ j ];
				surfaceOf.push_back( i );
				nodeR.push_back( radius );
				nodeZ.push_back( z );
				cosine.push_back( c );
				sine.push_back( sn );
				double const a = ( radius - axis.radius )*c + ( z - axis.z )*sn;
				nodeA.push_back( a );
				size = std::max( size, a );
				int const element = j < surface.elements.size() ? surface.elements[ j ] : -1;
				bool const band = ( j < surface.extended.size() && surface.extended[ j ] != 0 ) || element < 0;
				if ( band )
				{
					gathered.push_back( -1 );
					bandNodes.push_back( static_cast<int>( nodeR.size() - 1 ) );
					continue;
				}
				gathered.push_back( interior++ );

				mesh.GetElementTransformation( element, &transformation );
				mfem::InverseElementTransformation inverse( &transformation );
				inverse.SetSolverType( mfem::InverseElementTransformation::Newton );
				mfem::IntegrationPoint ip;
				mfem::Vector physical( 2 );
				physical( 0 ) = radius;
				physical( 1 ) = z;
				inverse.Transform( physical, ip );
				transformation.SetIntPoint( &ip );

				mfem::FiniteElement const &potentialElement = *potentialSpace.GetFE( element );
				mfem::FiniteElement const &fluxElement = *fluxSpace.GetFE( element );
				int const np = potentialElement.GetDof();
				int const nq = fluxElement.GetDof();
				shape.SetSize( np );
				dshape.SetSize( np, 2 );
				fluxShape.SetSize( nq );
				fluxDShape.SetSize( nq, 2 );
				potentialElement.CalcPhysShape( transformation, shape );
				potentialElement.CalcPhysDShape( transformation, dshape );
				fluxElement.CalcPhysShape( transformation, fluxShape );
				fluxElement.CalcPhysDShape( transformation, fluxDShape );
				potentialSpace.GetElementDofs( element, potentialDofs );
				fluxSpace.GetElementVDofs( element, fluxDofs );

				// The block's columns: psi~'s dofs, then q~'s two components'.
				int const columns = np + 2*nq;
				std::vector<double> sign( static_cast<std::size_t>( columns ), 1.0 );
				for ( int k = 0; k < np; ++k )
				{
					int dof = potentialDofs[ k ];
					if ( dof < 0 )
					{
						dof = -1 - dof;
						sign[ static_cast<std::size_t>( k ) ] = -1.0;
					}
					inIndex.push_back( dof );
				}
				for ( int k = 0; k < 2*nq; ++k )
				{
					int dof = fluxDofs[ k ];
					if ( dof < 0 )
					{
						dof = -1 - dof;
						sign[ static_cast<std::size_t>( np + k ) ] = -1.0;
					}
					inIndex.push_back( potentialSize + dof );
				}
				inOffset.push_back( static_cast<int>( inIndex.size() ) );

				std::vector<double> rows( static_cast<std::size_t>( baseRows*columns ), 0.0 );
				auto entry = [ &rows, columns ]( int row, int column ) -> double &
				{
					return rows[ static_cast<std::size_t>( row*columns + column ) ];
				};
				for ( int k = 0; k < np; ++k )
				{
					entry( Psi, k ) = shape( k );
					entry( PsiR, k ) = dshape( k, 0 );
					entry( PsiZ, k ) = dshape( k, 1 );
				}
				for ( int k = 0; k < nq; ++k )
				{
					entry( FluxR, np + k ) = fluxShape( k );
					entry( FluxZ, np + nq + k ) = fluxShape( k );
					entry( FluxRR, np + k ) = fluxDShape( k, 0 );
					entry( FluxRZ, np + k ) = fluxDShape( k, 1 );
					entry( FluxZR, np + nq + k ) = fluxDShape( k, 0 );
					entry( FluxZZ, np + nq + k ) = fluxDShape( k, 1 );
				}
				for ( int row = 0; row < baseRows; ++row )
					for ( int column = 0; column < columns; ++column )
						matrix.push_back( entry( row, column )*sign[ static_cast<std::size_t>( column ) ] );
				// The tangent's: the values alone -- psi, q_R, q_Z.
				for ( int row : { static_cast<int>( Psi ), static_cast<int>( FluxR ), static_cast<int>( FluxZ ) } )
					for ( int column = 0; column < columns; ++column )
						tangentMatrix.push_back( entry( row, column )*sign[ static_cast<std::size_t>( column ) ] );
				for ( int row = 0; row < baseRows; ++row )
					outIndex.push_back( ( interior - 1 )*baseRows + row );
				for ( int row = 0; row < tangentRows; ++row )
					tangentOutIndex.push_back( ( interior - 1 )*tangentRows + row );
				outOffset.push_back( static_cast<int>( outIndex.size() ) );
				tangentOutOffset.push_back( static_cast<int>( tangentOutIndex.size() ) );
				matrixOffset.push_back( static_cast<int>( matrix.size() ) );
				tangentMatrixOffset.push_back( static_cast<int>( tangentMatrix.size() ) );
			}
		}
		std::size_t const nodes = nodeR.size();
		lastShapeNodes = static_cast<int>( nodes );
		lastShapeBandNodes = static_cast<int>( bandNodes.size() );

		auto fill = [ interior ]( LinearBlocks &blocks, std::vector<int> const &outs, std::vector<int> const &outIdx,
		                          std::vector<int> const &ins, std::vector<int> const &inIdx,
		                          std::vector<int> const &matrixOffsets, std::vector<double> const &values )
		{
			blocks.count = interior;
			auto copy = []( mfem::Array<int> &to, std::vector<int> const &from )
			{
				to.SetSize( static_cast<int>( from.size() ) );
				for ( std::size_t k = 0; k < from.size(); ++k )
					to[ static_cast<int>( k ) ] = from[ k ];
			};
			copy( blocks.outOffset, outs );
			copy( blocks.outIndex, outIdx );
			copy( blocks.inOffset, ins );
			copy( blocks.inIndex, inIdx );
			copy( blocks.matrixOffset, matrixOffsets );
			blocks.knownOffset.SetSize( interior + 1 );
			blocks.knownOffset = 0;
			blocks.knownIndex.SetSize( 0 );
			blocks.matrix.SetSize( static_cast<int>( values.size() ) );
			for ( std::size_t k = 0; k < values.size(); ++k )
				blocks.matrix( static_cast<int>( k ) ) = static_cast<mfem::real_t>( values[ k ] );
		};
		fill( baseBlocks, outOffset, outIndex, inOffset, inIndex, matrixOffset, matrix );
		fill( tangentBlocks, tangentOutOffset, tangentOutIndex, inOffset, inIndex, tangentMatrixOffset,
		      tangentMatrix );

		/*
		 * THE BASE AT EVERY NODE: gathered where it can be, the conductors
		 * added in double; sampled and differenced in the band.
		 */
		std::vector<NodeValues> values( nodes );
		mfem::Vector baseGathered( interior*baseRows );
		lastShapeOnDevice = gather( baseBlocks, *base.potential, *base.flux, options.gatherRoute, baseGathered );
		{
			mfem::real_t const *g = baseGathered.HostRead();
			ConductorField const *conductors = solver.conductorField();
			for ( std::size_t n = 0; n < nodes; ++n )
			{
				if ( gathered[ n ] < 0 )
					continue;
				mfem::real_t const *row = g + gathered[ n ]*baseRows;
				NodeValues &v = values[ n ];
				v.psiR = row[ PsiR ];
				v.psiZ = row[ PsiZ ];
				v.qR = row[ FluxR ];
				v.qZ = row[ FluxZ ];
				v.qRR = row[ FluxRR ];
				v.qRZ = row[ FluxRZ ];
				v.qZR = row[ FluxZR ];
				v.qZZ = row[ FluxZZ ];
				if ( conductors )
				{
					double const radius = nodeR[ n ], z = nodeZ[ n ];
					double dR = 0.0, dZ = 0.0, cR = 0.0, cZ = 0.0;
					conductors->gradPsi( radius, z, dR, dZ );
					conductors->flux( radius, z, cR, cZ );
					v.psiR += dR;
					v.psiZ += dZ;
					v.qR += cR;
					v.qZ += cZ;
					double const step = 1.0e-6*std::max( 1.0, std::abs( radius ) );
					double aR = 0.0, aZ = 0.0, bR = 0.0, bZ = 0.0;
					conductors->flux( radius + step, z, aR, aZ );
					conductors->flux( radius - step, z, bR, bZ );
					v.qRR += ( aR - bR )/( 2.0*step );
					v.qZR += ( aZ - bZ )/( 2.0*step );
					conductors->flux( radius, z + step, aR, aZ );
					conductors->flux( radius, z - step, bR, bZ );
					v.qRZ += ( aR - bR )/( 2.0*step );
					v.qZZ += ( aZ - bZ )/( 2.0*step );
				}
			}
		}

		std::unique_ptr<ContourTracer> baseTracer;
		double const bandStep = std::cbrt( static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() ) )
		                        *std::max( size, 1.0e-12 );
		auto bandTracer = [ this ]( mfem::GridFunction const &potential, mfem::GridFunction const &flux,
		                            bool conductors )
		{
			auto tracer = std::make_unique<ContourTracer>( potential, flux );
			if ( conductors )
				tracer->setConductorField( solver.conductorField() );
			if ( options.gammaHMarker != nullptr )
				tracer->setBandExtension( BandExtension::TransferLift, *options.gammaHMarker, options.path );
			return tracer;
		};
		if ( !bandNodes.empty() )
		{
			baseTracer = bandTracer( *base.potential, *base.flux, true );
			for ( int n : bandNodes )
			{
				double const radius = nodeR[ static_cast<std::size_t>( n ) ];
				double const z = nodeZ[ static_cast<std::size_t>( n ) ];
				double psi[ 4 ], qR[ 5 ], qZ[ 5 ];
				double const dR[ 5 ] = { bandStep, -bandStep, 0.0, 0.0, 0.0 };
				double const dZ[ 5 ] = { 0.0, 0.0, bandStep, -bandStep, 0.0 };
				for ( int k = 0; k < 5; ++k )
				{
					double value = 0.0;
					if ( !baseTracer->sampleAt( radius + dR[ k ], z + dZ[ k ], value, qR[ k ], qZ[ k ] ) )
						return false;
					if ( k < 4 )
						psi[ k ] = value;
				}
				NodeValues &v = values[ static_cast<std::size_t>( n ) ];
				v.psiR = ( psi[ 0 ] - psi[ 1 ] )/( 2.0*bandStep );
				v.psiZ = ( psi[ 2 ] - psi[ 3 ] )/( 2.0*bandStep );
				v.qR = qR[ 4 ];
				v.qZ = qZ[ 4 ];
				v.qRR = ( qR[ 0 ] - qR[ 1 ] )/( 2.0*bandStep );
				v.qZR = ( qZ[ 0 ] - qZ[ 1 ] )/( 2.0*bandStep );
				v.qRZ = ( qR[ 2 ] - qR[ 3 ] )/( 2.0*bandStep );
				v.qZZ = ( qZ[ 2 ] - qZ[ 3 ] )/( 2.0*bandStep );
			}
		}

		// Each node's own terms: f, u . grad psi, u . Q.
		std::vector<double> f( nodes ), gradAlong( nodes ), fluxAlong( nodes );
		for ( std::size_t n = 0; n < nodes; ++n )
		{
			NodeValues const &v = values[ n ];
			gradAlong[ n ] = v.psiR*cosine[ n ] + v.psiZ*sine[ n ];
			fluxAlong[ n ] = v.qR*cosine[ n ] + v.qZ*sine[ n ];
			if ( !( std::abs( gradAlong[ n ] ) > 0.0 ) || !( std::abs( fluxAlong[ n ] ) > 0.0 )
			     || !( nodeA[ n ] > 0.0 ) || !( nodeR[ n ] > 0.0 ) )
				return false;
			f[ n ] = nodeA[ n ]/( std::abs( fluxAlong[ n ] )*nodeR[ n ]*nodeR[ n ] );
		}

		/*
		 * EACH TANGENT: its continuous pair, gathered at the interior nodes and
		 * sampled at the band's, then each node's df summed per surface.
		 */
		out.assign( m*m, 0.0 );
		mfem::Vector tangentGathered( interior*tangentRows );
		for ( std::size_t j = 0; j < m; ++j )
		{
			Fields const direction = along( dx[ j ] );
			ContinuousFields const tangent =
				solver.continuousInterpolant().apply( *direction.psiStar, *direction.qStar, *direction.trace );
			gather( tangentBlocks, *tangent.potential, *tangent.flux, options.gatherRoute, tangentGathered );
			mfem::real_t const *g = tangentGathered.HostRead();
			std::unique_ptr<ContourTracer> tangentTracer;
			if ( !bandNodes.empty() )
				tangentTracer = bandTracer( *tangent.potential, *tangent.flux, false );

			std::vector<double> sums( m, 0.0 );
			std::vector<std::size_t> counts( m, 0 );
			for ( std::size_t n = 0; n < nodes; ++n )
			{
				double dPsi = 0.0, dQR = 0.0, dQZ = 0.0;
				if ( gathered[ n ] >= 0 )
				{
					mfem::real_t const *row = g + gathered[ n ]*tangentRows;
					dPsi = row[ 0 ];
					dQR = row[ 1 ];
					dQZ = row[ 2 ];
				}
				else if ( !tangentTracer->sampleAt( nodeR[ n ], nodeZ[ n ], dPsi, dQR, dQZ ) )
					return false;
				std::size_t const i = surfaceOf[ n ];
				double const normalised = family.surfaces[ i ].normalisedFlux;
				double const dLevel = ( 1.0 - normalised )*dAxis[ j ] + normalised*dBoundary[ j ];
				NodeValues const &v = values[ n ];
				double const dA = ( dLevel - dPsi - v.psiR*originR[ j ] - v.psiZ*originZ[ j ] )/gradAlong[ n ];
				double const dR = originR[ j ] + dA*cosine[ n ];
				double const dZ = originZ[ j ] + dA*sine[ n ];
				double const dFluxAlong = cosine[ n ]*( dQR + v.qRR*dR + v.qRZ*dZ )
				                          + sine[ n ]*( dQZ + v.qZR*dR + v.qZZ*dZ );
				sums[ i ] += f[ n ]*( dA/nodeA[ n ] - dFluxAlong/fluxAlong[ n ] - 2.0*dR/nodeR[ n ] );
				counts[ i ] += 1;
			}
			for ( std::size_t i = 0; i < m; ++i )
				out[ i*m + j ] = sums[ i ]/static_cast<double>( counts[ i ] );
		}
		return true;
	}

}
