#include "PointLocator.hpp"
#include "Crossings.hpp"
#include "Sampler.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>

namespace meq
{

	namespace
	{

		/// Whether the mesh's geometry map is affine on simplices.
		///
		/// THE SAME TEST meq::GridSampler MAKES, AND FOR THE SAME REASON: a
		/// mesh acquires nodes without becoming curved, so testing
		/// Mesh::GetNodes() against null drops every later query onto the
		/// Newton route for the life of the process -- measured at 159x and 83x
		/// there, with the answer unchanged, so only a clock can see it. The
		/// question is the geometry's DEGREE, and the values are compared with
		/// the vertex array rather than assumed equal to it, since MFEM lets a
		/// caller deform a mesh through its nodes and leave the vertices behind.
		bool affineGeometry( mfem::Mesh &mesh )
		{
			mfem::GridFunction const *nodes = mesh.GetNodes();
			if ( nodes == nullptr )
				return true;
			// ON THE HOST, before this and every element transformation after
			// it reads them raw. MFEM's geometric factors read the nodes on the
			// device -- the batched reconstruction and assembly do, as MEQ's own
			// device routes do -- and under the debug backend that leaves the
			// host page protected until something asks for it back. Free where
			// both copies are valid, as under CUDA.
			nodes->HostRead();

			mfem::FiniteElementSpace const *space = nodes->FESpace();
			if ( space == nullptr || space->GetMaxElementOrder() > 1 )
				return false;

			// An order-1 nodal field agrees with the vertices, or somebody has
			// moved one without the other.
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				mfem::Array<int> vertices;
				mesh.GetElementVertices( e, vertices );
				mfem::Array<int> dofs;
				space->GetElementVDofs( e, dofs );
				if ( dofs.Size() != 2*vertices.Size() )
					return false;

				for ( int v = 0; v < vertices.Size(); ++v )
				{
					auto const *x = mesh.GetVertex( vertices[ v ] );
					for ( int d = 0; d < 2; ++d )
					{
						mfem::real_t const moved =
							( *nodes )( dofs[ d*vertices.Size() + v ] ) - x[ d ];
						// Same rule as insideTolerance below: a bare 1e-12 is
						// under FP32's epsilon, so in single precision an
						// undeformed mesh would read as deformed and every
						// element would take the Newton route.
						if ( std::abs( moved ) > std::max(
							mfem::real_t( 1.0e-12 ),
							mfem::real_t( 16 )*std::numeric_limits<mfem::real_t>::epsilon() ) )
							return false;
					}
				}
			}
			return true;
		}

		/// Geometry::CheckPoint()'s own slack, so this path accepts exactly the
		/// candidates TransformBack() accepts rather than nearly them -- OR the
		/// format's own noise floor, whichever is larger.
		///
		/// THE SECOND HALF IS NOT DEFENSIVE, IT IS A MEASURED DEFECT IN THIS
		/// FILE. Written as a bare `1.0e-8` it is below FP32's epsilon of
		/// 1.19e-07, so under MFEM_PRECISION = single the 2x2 solve's own
		/// round-off exceeds the tolerance and a point genuinely inside its
		/// element is rejected: measured on examples/limited-tokamak.toml,
		/// **1232 of 25616 warm-start nodes fell outside** where the double
		/// build missed none, and the solve then converged to the wrong
		/// topology -- caught only because checkAxis() refused the answer.
		///
		/// It is the same defect as MFEM's own `MFEM_VERIFY( err < 1e-10 )` in
		/// DarcyHybridization::TraceEmbedding(), which MEQ reported upstream
		/// the same afternoon this was written. **A tolerance that does not
		/// follow real_t is a tolerance that is wrong in one of the two
		/// builds**, and which one is decided by a constant nobody re-reads.
		/// 128 epsilon is the slack of a barycentric coordinate coming out of
		/// a 2x2 solve, and in double it is 2.8e-14 -- far below 1e-8, so the
		/// double path keeps exactly MFEM's own number and is unchanged.
		mfem::real_t const insideTolerance = std::max(
			mfem::real_t( 1.0e-8 ),
			mfem::real_t( 128 )*std::numeric_limits<mfem::real_t>::epsilon() );

	}

	PointLocator::PointLocator( mfem::Mesh &mesh )
		: meshRef( mesh )
	{
		buildBins();
	}

	PointLocator::BinGrid PointLocator::grid() const
	{
		BinGrid out{};
		out.start = binStart.HostRead();
		out.element = binElement.HostRead();
		out.rMin = rMin;
		out.zMin = zMin;
		out.dR = dR;
		out.dZ = dZ;
		out.binsR = binsR;
		out.binsZ = binsZ;
		out.affine = !anyCurved;
		return out;
	}

	void PointLocator::buildBins()
	{
		int const elements = meshRef.GetNE();
		maps.assign( static_cast<std::size_t>( elements ), ElementMap{} );

		bool const straightSided = affineGeometry( meshRef );
		anyCurved = false;

		// One pass for the element maps and the global bounding box.
		mfem::real_t rLo =  std::numeric_limits<mfem::real_t>::max();
		mfem::real_t rHi = -std::numeric_limits<mfem::real_t>::max();
		mfem::real_t zLo =  std::numeric_limits<mfem::real_t>::max();
		mfem::real_t zHi = -std::numeric_limits<mfem::real_t>::max();

		std::vector<mfem::real_t> boxes( 4*static_cast<std::size_t>( elements ) );

		for ( int e = 0; e < elements; ++e )
		{
			mfem::Array<int> vertices;
			meshRef.GetElementVertices( e, vertices );

			mfem::real_t eRLo =  std::numeric_limits<mfem::real_t>::max();
			mfem::real_t eRHi = -std::numeric_limits<mfem::real_t>::max();
			mfem::real_t eZLo =  std::numeric_limits<mfem::real_t>::max();
			mfem::real_t eZHi = -std::numeric_limits<mfem::real_t>::max();
			for ( int v = 0; v < vertices.Size(); ++v )
			{
				auto const *p = meshRef.GetVertex( vertices[ v ] );
				eRLo = std::min( eRLo, p[ 0 ] );
				eRHi = std::max( eRHi, p[ 0 ] );
				eZLo = std::min( eZLo, p[ 1 ] );
				eZHi = std::max( eZHi, p[ 1 ] );
			}

			/*
			 * A CURVED ELEMENT BULGES OUTSIDE THE HULL OF ITS VERTICES, so a
			 * box drawn on the vertices alone can exclude a point that is
			 * genuinely inside -- and the failure is a MISS rather than a wrong
			 * answer, which is the quiet kind here because a missed warm-start
			 * point silently falls back to the Dirichlet datum. The box is
			 * therefore padded by a tenth of its own size on a curved mesh.
			 */
			bool affine = straightSided && vertices.Size() == 3
			              && meshRef.GetElementBaseGeometry( e ) == mfem::Geometry::TRIANGLE;
			if ( !affine )
			{
				anyCurved = true;
				mfem::real_t const padR = mfem::real_t( 0.1 )*( eRHi - eRLo );
				mfem::real_t const padZ = mfem::real_t( 0.1 )*( eZHi - eZLo );
				eRLo -= padR; eRHi += padR;
				eZLo -= padZ; eZHi += padZ;
			}
			else
			{
				auto const *p0 = meshRef.GetVertex( vertices[ 0 ] );
				auto const *p1 = meshRef.GetVertex( vertices[ 1 ] );
				auto const *p2 = meshRef.GetVertex( vertices[ 2 ] );
				mfem::real_t const j00 = p1[ 0 ] - p0[ 0 ], j01 = p2[ 0 ] - p0[ 0 ];
				mfem::real_t const j10 = p1[ 1 ] - p0[ 1 ], j11 = p2[ 1 ] - p0[ 1 ];
				mfem::real_t const determinant = j00*j11 - j01*j10;
				if ( determinant == mfem::real_t( 0.0 ) )
				{
					// Degenerate. Newton will call it Outside, which is right.
					affine = false;
					anyCurved = true;
				}
				else
				{
					ElementMap &m = maps[ static_cast<std::size_t>( e ) ];
					m.affine = true;
					m.originR = p0[ 0 ];
					m.originZ = p0[ 1 ];
					m.inverse[ 0 ] =  j11/determinant;
					m.inverse[ 1 ] = -j01/determinant;
					m.inverse[ 2 ] = -j10/determinant;
					m.inverse[ 3 ] =  j00/determinant;
				}
			}

			boxes[ 4*static_cast<std::size_t>( e ) + 0 ] = eRLo;
			boxes[ 4*static_cast<std::size_t>( e ) + 1 ] = eRHi;
			boxes[ 4*static_cast<std::size_t>( e ) + 2 ] = eZLo;
			boxes[ 4*static_cast<std::size_t>( e ) + 3 ] = eZHi;

			rLo = std::min( rLo, eRLo ); rHi = std::max( rHi, eRHi );
			zLo = std::min( zLo, eZLo ); zHi = std::max( zHi, eZHi );
		}

		if ( elements == 0 )
		{
			binsR = binsZ = 1;
			binStart.SetSize( 2 );
			binStart[ 0 ] = binStart[ 1 ] = 0;
			binElement.SetSize( 0 );
			return;
		}

		/*
		 * ABOUT ONE ELEMENT PER BIN, which is the usual choice and is a trade
		 * rather than a tuning: fewer bins means longer candidate lists and
		 * more 2x2 solves per point, more bins means a bigger index and more
		 * empty-bin walking for a point near nothing. sqrt( NE ) per axis puts
		 * a mesh of N elements in N bins.
		 */
		int const perAxis = std::max( 1, static_cast<int>( std::sqrt(
			static_cast<double>( elements ) ) ) );
		binsR = binsZ = perAxis;

		rMin = rLo;
		zMin = zLo;
		dR = ( rHi - rLo )/static_cast<mfem::real_t>( binsR );
		dZ = ( zHi - zLo )/static_cast<mfem::real_t>( binsZ );
		if ( dR <= mfem::real_t( 0.0 ) ) dR = mfem::real_t( 1.0 );
		if ( dZ <= mfem::real_t( 0.0 ) ) dZ = mfem::real_t( 1.0 );

		auto binOf = [ & ]( mfem::real_t x, mfem::real_t origin, mfem::real_t step,
		                    int count )
		{
			int const b = static_cast<int>( std::floor( ( x - origin )/step ) );
			return std::min( count - 1, std::max( 0, b ) );
		};

		// Counting sort into CSR: one pass to count, one to fill.
		int const bins = binsR*binsZ;
		binStart.SetSize( bins + 1 );
		binStart = 0;

		auto forEachBin = [ & ]( int e, auto &&visit )
		{
			mfem::real_t const *b = &boxes[ 4*static_cast<std::size_t>( e ) ];
			int const i0 = binOf( b[ 0 ], rMin, dR, binsR );
			int const i1 = binOf( b[ 1 ], rMin, dR, binsR );
			int const j0 = binOf( b[ 2 ], zMin, dZ, binsZ );
			int const j1 = binOf( b[ 3 ], zMin, dZ, binsZ );
			for ( int j = j0; j <= j1; ++j )
				for ( int i = i0; i <= i1; ++i )
					visit( j*binsR + i );
		};

		for ( int e = 0; e < elements; ++e )
			forEachBin( e, [ & ]( int b ) { ++binStart[ b + 1 ]; } );

		for ( int b = 0; b < bins; ++b )
			binStart[ b + 1 ] += binStart[ b ];

		binElement.SetSize( binStart[ bins ] );
		mfem::Array<int> cursor( binStart );
		for ( int e = 0; e < elements; ++e )
			forEachBin( e, [ & ]( int b ) { binElement[ cursor[ b ]++ ] = e; } );
	}

	void PointLocator::locate( mfem::Vector const &positions )
	{
		points = positions.Size()/2;
		element.SetSize( points );
		code.SetSize( points );
		reference.SetSize( 2*points );
		distance.SetSize( points );
		misses = 0;
		deviceUsed = false;

		if ( points == 0 )
			return;


		element = -1;
		code = 2u;
		reference = 0.0;
		distance = 0.0;

		/*
		 * THE DEVICE PATH IS TRIED FIRST AND ONLY WHERE IT IS WHOLE. A mesh
		 * with any curved or degenerate element keeps the host, because the
		 * fallback for those is an mfem::ElementTransformation Newton and that
		 * is not a device object. Splitting one cloud across both would mean
		 * two locators answering one query, which is the arrangement this file
		 * exists to avoid.
		 */
		if ( !anyCurved )
		{
			countedHostRead( positions, __LINE__ );
			deviceUsed = locateOnDevice(
				maps.data(), meshRef.GetNE(),
				binStart.GetData(), binElement.GetData(),
				rMin, zMin, dR, dZ, binsR, binsZ,
				positions.GetData(), points,
				element.GetData(), code.GetData(),
				reference.GetData(), distance.GetData() );
		}

		if ( !deviceUsed )
		{
			/*
			 * ONE POINT AT A TIME, AND EVERY POINT IS ITS OWN QUESTION.
			 *
			 * A point reads the bins, the element maps and the mesh, none of
			 * which it writes, and writes only its own entries of `element`,
			 * `code`, `reference` and `distance`. So the loop below is threaded
			 * over points with nothing shared but read-only state, and the
			 * answer does not depend on the thread count: each point's search
			 * is the same arithmetic in the same order whichever thread runs
			 * it. The two pieces of scratch the search needs -- the element
			 * transformation and the physical point -- are per thread, and the
			 * transformation is filled through the CALLER-ALLOCATED
			 * GetElementTransformation( e, & ), which CLAUDE.md's traps require
			 * of anything threaded.
			 *
			 * WHAT IT IS FOR: the cold start's warm-start transfer locates every
			 * quadrature point of the new mesh in the stored guess -- 61,116 of
			 * them on TCV shaped at k = 2 -- serially, inside makeSolver(),
			 * which --profile reads at one core.
			 */
			auto locateOne = [ & ]( int p, mfem::IsoparametricTransformation &transformation,
			                        mfem::Vector &physical )
			{
				mfem::real_t const r = positions( p );
				mfem::real_t const z = positions( points + p );

				int const i = std::min( binsR - 1, std::max( 0,
					static_cast<int>( std::floor( ( r - rMin )/dR ) ) ) );
				int const j = std::min( binsZ - 1, std::max( 0,
					static_cast<int>( std::floor( ( z - zMin )/dZ ) ) ) );

				mfem::real_t best = std::numeric_limits<mfem::real_t>::max();
				int bestElement = -1;      ///< strictly inside, or -1
				int bestCandidate = -1;    ///< least-outside, for the snap
				mfem::real_t bestXi = 0.0, bestEta = 0.0;

				/*
				 * THE RING WIDENS, AND ONE RING IS NOT ENOUGH ON ITS OWN. A
				 * point outside the mesh entirely -- which a refined target
				 * domain really does produce, since D_h grows outwards every
				 * adaptive cycle -- has no element in its own bin, and the
				 * nearest one may be several bins away. Widening until
				 * something is found, capped, costs nothing for the common
				 * case where the first ring answers and keeps the snap
				 * meaningful for the rest.
				 */
				int const maxRing = std::max( binsR, binsZ );
				for ( int ring = 0; ring <= maxRing && bestElement < 0; ++ring )
				{
					int const i0 = std::max( 0, i - ring ), i1 = std::min( binsR - 1, i + ring );
					int const j0 = std::max( 0, j - ring ), j1 = std::min( binsZ - 1, j + ring );

					for ( int jj = j0; jj <= j1; ++jj )
						for ( int ii = i0; ii <= i1; ++ii )
						{
							// Only the new ring, not the filled square again.
							if ( ring > 0 && ii != i0 && ii != i1 && jj != j0 && jj != j1 )
								continue;

							int const b = jj*binsR + ii;
							for ( int k = binStart[ b ]; k < binStart[ b + 1 ]; ++k )
							{
								int const e = binElement[ k ];
								ElementMap const &m = maps[ static_cast<std::size_t>( e ) ];

								mfem::real_t xi = 0.0, eta = 0.0;
								mfem::real_t overshoot =
									std::numeric_limits<mfem::real_t>::max();
								bool isInside = false;

								if ( m.affine )
								{
									overshoot = locateAffine( m, r, z, xi, eta );
									isInside = overshoot <= insideTolerance;
								}
								else
								{
									/*
									 * THE REFERENCE POINT IS TAKEN WHATEVER THE
									 * STATUS, which is the half that makes the
									 * snap possible on a curved element.
									 * TransformBack() answers Inside or Outside
									 * and writes its iterate either way, so an
									 * Outside result still says WHERE the
									 * Newton landed -- and for a point a hair
									 * beyond a curved face that is the right
									 * neighbourhood. Reading only the status
									 * throws that away and is what the first
									 * version of this file did.
									 */
									physical( 0 ) = r;
									physical( 1 ) = z;
									meshRef.GetElementTransformation( e, &transformation );
									mfem::IntegrationPoint found;
									int status =
										transformation.TransformBack(
											physical, found,
											inverseTransformTolerance() );

									/*
									 * A SECOND, MORE ROBUST ATTEMPT WHERE THE
									 * CHEAP ONE GIVES UP -- AND ONLY THERE, SO
									 * THE COMMON PATH IS UNCHANGED TO THE BIT.
									 *
									 * TransformBack() uses
									 * InverseElementTransformation's defaults,
									 * which are `Center` for the initial guess.
									 * MFEM's own note ranks that second of four
									 * for reliability and says `ClosestPhysNode`
									 * is "much more reliable, especially in the
									 * case of highly distorted elements" -- and
									 * adds, for all four, that **none of them
									 * guarantees a point inside an element will
									 * be found**.
									 *
									 * That is not a theoretical caveat here.
									 * Under MFEM_PRECISION = single on
									 * limited-tokamak, whose stored guess is a
									 * CURVED mesh, the default lost **458 of
									 * 25616** warm-start nodes where the double
									 * build lost none -- the Newton landing
									 * just outside a reference element on a
									 * point that is genuinely inside it, which
									 * is a round-off decision at FP32's floor.
									 *
									 * Retrying costs nothing where the default
									 * succeeds, which is every point on the
									 * double path, so this cannot move an
									 * answer that was already being found.
									 */
									if ( status != mfem::InverseElementTransformation::Inside )
									{
										mfem::InverseElementTransformation inverse;
										inverse.SetTransformation( transformation );
										inverse.SetInitialGuessType(
											mfem::InverseElementTransformation::ClosestPhysNode );
										inverse.SetSolverType(
											mfem::InverseElementTransformation::NewtonElementProject );
										mfem::IntegrationPoint retry;
										int const second = inverse.Transform( physical, retry );
										if ( second == mfem::InverseElementTransformation::Inside )
										{
											status = second;
											found = retry;
										}
									}

									xi = found.x;
									eta = found.y;

									/*
									 * TransformBack()'S OWN VERDICT IS THE
									 * AUTHORITY ON "INSIDE", AND OVERRIDING IT
									 * WITH THE CLAMPED OVERSHOOT COST TWELVE
									 * NEWTON ITERATIONS. A version of this loop
									 * accepted any candidate whose reference
									 * coordinates were within tolerance,
									 * without asking the status -- so on a
									 * curved mesh the first element in bin
									 * order could claim a point belonging to
									 * its neighbour, the warm start came back
									 * subtly wrong, and limited-tokamak went
									 * from 12 iterations to 32 with the
									 * constraint falling from 8.5e-16 to
									 * 1.0e-11. The answer was still right,
									 * which is what made it a cost rather than
									 * a failure and is exactly why it needs
									 * saying: a locator can be wrong about
									 * WHICH element without being wrong about
									 * the value, and only the iteration count
									 * shows it.
									 *
									 * The overshoot below is therefore a
									 * CANDIDATE measure for the snap alone, and
									 * never promotes a point to interior.
									 */
									if ( status == mfem::InverseElementTransformation::Inside )
									{
										overshoot = mfem::real_t( 0.0 );
										isInside = true;
									}
									else
									{
										mfem::real_t const a = -xi;
										mfem::real_t const b = -eta;
										mfem::real_t const c =
											xi + eta - mfem::real_t( 1.0 );
										mfem::real_t worst = a > b ? a : b;
										overshoot = worst > c ? worst : c;
										// A Newton that ran away rather than
										// converging is not a candidate.
										if ( !( overshoot < mfem::real_t( 1.0 ) ) )
											overshoot =
												std::numeric_limits<mfem::real_t>::max();
									}
								}

								if ( overshoot < best )
								{
									best = overshoot;
									bestCandidate = e;
									bestXi = xi;
									bestEta = eta;
								}
								if ( isInside )
								{
									bestElement = e;
									bestXi = xi;
									bestEta = eta;
									break;
								}
							}
							if ( bestElement >= 0 )
								break;
						}
				}

				/*
				 * THE SNAP, AND LEAVING IT OUT WAS A MEASURED DEFECT RATHER
				 * THAN AN OMISSION. The loop above accepts only a point
				 * STRICTLY inside an element, and a first version of this file
				 * stopped there -- which is correct on a mesh whose elements
				 * are straight, and wrong the moment they are not. Measured on
				 * examples/limited-tokamak.toml, whose stored guess is curved:
				 * **1232 of 25616 nodes, exactly 77 whole target elements,
				 * found nothing** under MFEM_PRECISION = single while the
				 * double build found all of them, and the warm start then put
				 * the solve on the wrong topology.
				 *
				 * THE CAUSE IS NOT A TOLERANCE AND THAT IS THE POINT.
				 * InverseElementTransformation already scales its own by the
				 * format -- ref_tol is 1e-15 in double and 4e-7 in single
				 * (fem/eltrans.hpp:306-311) -- so TransformBack is doing what
				 * it promises. A curved element's Newton simply lands a hair
				 * outside the reference triangle for a point on its boundary,
				 * and whether it lands inside or outside is decided by
				 * round-off. A locator that has no answer for "outside every
				 * element, but only just" therefore reports a MISS at the
				 * domain boundary, and a warm start's misses are silent: they
				 * fall back to the Dirichlet datum and the run merely starts
				 * somewhere else.
				 *
				 * gslib's code 1 is exactly this case, which is why the codes
				 * were kept rather than simplified to found/not-found.
				 */
				if ( bestElement < 0 && bestCandidate >= 0 )
				{
					// Clamp to the reference triangle and measure how far the
					// point actually moved, in PHYSICAL space -- the distance a
					// caller can compare against a mesh size. A clamp in
					// reference coordinates is not a distance in metres.
					mfem::real_t xi = bestXi, eta = bestEta;
					xi  = xi  < mfem::real_t( 0.0 ) ? mfem::real_t( 0.0 ) : xi;
					eta = eta < mfem::real_t( 0.0 ) ? mfem::real_t( 0.0 ) : eta;
					mfem::real_t const sum = xi + eta;
					if ( sum > mfem::real_t( 1.0 ) )
					{
						xi  /= sum;
						eta /= sum;
					}

					meshRef.GetElementTransformation( bestCandidate, &transformation );
					mfem::IntegrationPoint snapped;
					snapped.Set2( xi, eta );
					transformation.SetIntPoint( &snapped );
					mfem::Vector landed( 2 );
					transformation.Transform( snapped, landed );

					mfem::real_t const dr = landed( 0 ) - r;
					mfem::real_t const dz = landed( 1 ) - z;
					mfem::real_t const moved = std::sqrt( dr*dr + dz*dz );

					/*
					 * A CAP, BECAUSE A SNAP WITHOUT ONE IS NOT A MISS REPORT.
					 * Every point would then find SOME element -- the nearest
					 * one, however far -- and `misses` would read zero on a
					 * target that reaches a metre outside the stored mesh. The
					 * count is what tells a caller their restart quietly became
					 * a cold start, so the snap is allowed a bin and no more.
					 */
					mfem::real_t const reach = dR > dZ ? dR : dZ;
					if ( moved <= reach )
					{
						bestElement = bestCandidate;
						reference( 2*p )     = xi;
						reference( 2*p + 1 ) = eta;
						distance( p ) = moved;
						element[ p ] = bestElement;
						code[ p ] = 1u;
					}
				}
				else if ( bestElement >= 0 )
				{
					element[ p ] = bestElement;
					code[ p ] = 0u;
					reference( 2*p )     = bestXi;
					reference( 2*p + 1 ) = bestEta;
					distance( p ) = 0.0;
				}

				if ( element[ p ] < 0 )
					code[ p ] = 2u;
			};

			// Serial under an mfem::Device, which reaches this loop only for a
			// curved or degenerate mesh: the element transformations read mesh
			// nodes through MFEM's memory manager, which nothing documents as
			// reentrant, and a device run is a correctness instrument anyway.
			bool const hostThreads = !mfem::Device::Allows( mfem::Backend::DEVICE_MASK );
			int notFound = 0;
			std::exception_ptr failure;
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
			#pragma omp parallel if ( hostThreads ) reduction( + : notFound )
#endif
			{
				mfem::IsoparametricTransformation transformation;
				mfem::Vector physical( 2 );
				// Dynamic, because a point near nothing walks rings of bins and
				// one in a full bin stops at the first element: the work per
				// point is uneven, and a static split would leave threads idle.
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
				#pragma omp for schedule( dynamic, 256 )
#endif
				for ( int p = 0; p < points; ++p )
				{
					// An exception may not leave a parallel region, so the first
					// one is kept and rethrown once the region has closed.
					try
					{
						locateOne( p, transformation, physical );
						if ( element[ p ] < 0 )
							++notFound;
					}
					catch ( ... )
					{
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
						#pragma omp critical( meqPointLocatorFailure )
#endif
						if ( !failure )
							failure = std::current_exception();
					}
				}
			}
			if ( failure )
				std::rethrow_exception( failure );
			misses += notFound;
		}
		else
		{
			for ( int p = 0; p < points; ++p )
				if ( code[ p ] == 2u )
					++misses;
		}
	}

	void PointLocator::interpolate( mfem::GridFunction const &field,
	                                mfem::Vector &values ) const
	{
		values.SetSize( points );
		values = 0.0;
		if ( points == 0 )
			return;

		mfem::FiniteElementSpace const *space = field.FESpace();
		MFEM_VERIFY( space != nullptr,
		             "meq::PointLocator::interpolate: the field has no finite element space" );

		/*
		 * GROUPED BY ELEMENT, which is meq::GridSampler's 1.6x and the same
		 * argument: the dof list and the gathered coefficients are three heap
		 * allocations per point through GridFunction::GetValue(), and they are
		 * per ELEMENT. What is left inside the point loop is CalcShape, which
		 * genuinely is per point.
		 */
		std::vector<int> order( static_cast<std::size_t>( points ) );
		for ( int p = 0; p < points; ++p )
			order[ static_cast<std::size_t>( p ) ] = p;
		std::sort( order.begin(), order.end(),
		           [ & ]( int a, int b ) { return element[ a ] < element[ b ]; } );

		mfem::Array<int> dofs;
		mfem::Vector shape, local;
		countedHostRead( field, __LINE__ );

		std::size_t at = 0;
		while ( at < order.size() )
		{
			int const e = element[ order[ at ] ];
			std::size_t end = at;
			while ( end < order.size() && element[ order[ end ] ] == e )
				++end;

			if ( e < 0 )
			{
				at = end;      // the not-found points sort to the front
				continue;
			}

			mfem::FiniteElement const *fe = space->GetFE( e );
			space->GetElementDofs( e, dofs );
			local.SetSize( dofs.Size() );
			field.GetSubVector( dofs, local );
			shape.SetSize( fe->GetDof() );

			for ( std::size_t k = at; k < end; ++k )
			{
				int const p = order[ k ];
				mfem::IntegrationPoint ip;
				ip.Set2( reference( 2*p ), reference( 2*p + 1 ) );
				fe->CalcShape( ip, shape );
				values( p ) = shape*local;
			}

			at = end;
		}
	}

}
