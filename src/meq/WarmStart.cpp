#include "WarmStart.hpp"
#include "Crossings.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace meq
{

	FieldTransfer::FieldTransfer( mfem::Mesh &sourceMeshIn )
		: sourceMesh( sourceMeshIn ), locator( sourceMeshIn ), misses( 0 ),
		  points( 0 ), worstDist( 0.0 )
	{
		if ( sourceMesh.Dimension() != 2 )
			throw std::invalid_argument( "meq::FieldTransfer: the source mesh must be two dimensional ( R, z )" );

		// A mesh from Mesh::MakeCartesian2D or read from a .mesh file carries no
		// explicit nodal grid function until it is asked for one. EnsureNodes()
		// builds it at the mesh's own order, which for MEQ's straight-sided
		// triangles is degree 1 and exact. This is why the constructor takes a
		// non-const Mesh &: setting up the search modifies the mesh, and
		// pretending otherwise would need a copy of it.
		/*
		 * EnsureNodes() IS KEPT THOUGH THE LOCATOR DOES NOT NEED IT, and the
		 * reason is a measured trap rather than caution. gslib refused a mesh
		 * without nodes, so this call was made for it; meq::PointLocator asks
		 * the geometry's DEGREE instead and is happy either way. But the call
		 * has a second effect nobody asked for -- it installs a nodal field on
		 * the caller's mesh, permanently -- and meq::GridSampler's affine fast
		 * path keys on that same question. Removing the call would therefore
		 * change GridSampler's behaviour on every mesh a warm start has
		 * touched, which is a 159x change in a component this file does not
		 * own (M-93). It stays until that is deliberate rather than incidental.
		 */
		sourceMesh.EnsureNodes();

		/*
		 * NONE, AND THAT IS THE DIFFERENCE BETWEEN AN EXACT TRANSFER AND A LOSSY
		 * ONE. Measured before it was chosen.
		 *
		 * FindPointsGSLIB::Interpolate() evaluates every point exactly through
		 * InterpolateGeneral() and then, for an L2 field with an averaging type
		 * set, goes back over the points gslib classified as being ON AN ELEMENT
		 * BORDER (code 1), projects the whole field to H1 and re-interpolates
		 * those. For a continuous field that is a sensible tie-break. For a
		 * DISCONTINUOUS one it replaces the value with a smoothed one, and the
		 * smoothing is not small: with ARITHMETIC, transferring a converged
		 * solution onto a mesh whose space CONTAINS it -- where the exact answer
		 * exists -- came back 28%, 10% and 12% wrong in L2 at k = 1, 2, 3.
		 *
		 * MEQ does not need the tie-break. The target nodes are Gauss-Legendre
		 * points, strictly interior to their own element, so on a nested
		 * refinement each is interior to exactly one source element and the value
		 * there is unambiguous. Where a node genuinely does land on a border --
		 * two meshes that are not nested -- picking one side of a jump is a
		 * better guess than averaging across it, because the average belongs to
		 * neither element's polynomial.
		 *
		 * MEQ'S OWN LOCATOR HAS NO SUCH MODE AND THAT IS THE POINT OF THE
		 * PARAGRAPH SURVIVING THE SWITCH. meq::PointLocator returns the value
		 * of the element it located the point in and never averages across a
		 * face, so the behaviour above is the default rather than a setting --
		 * and the measurement that justified choosing NONE is what says the
		 * default is the right one.
		 */
	}

	FieldTransfer::~FieldTransfer() = default;

	int FieldTransfer::transfer( mfem::GridFunction const &source,
	                             mfem::Coefficient &fallback,
	                             mfem::GridFunction &target )
	{
		mfem::FiniteElementSpace *targetSpace = target.FESpace();
		if ( targetSpace == nullptr )
			throw std::invalid_argument( "meq::FieldTransfer::transfer: the target has no finite element space" );
		if ( source.FESpace() == nullptr || source.FESpace()->GetMesh() != &sourceMesh )
			throw std::invalid_argument( "meq::FieldTransfer::transfer: the source field does not live on the mesh this was constructed with" );

		mfem::Mesh &targetMesh = *targetSpace->GetMesh();
		int const elements = targetMesh.GetNE();

		/*
		 * ELEMENT-LOCAL L2 PROJECTION, NOT NODAL INTERPOLATION, AND THE REASON IS
		 * MEASURED RATHER THAN STYLISTIC.
		 *
		 * The obvious transfer reads the source at each target dof point. It is
		 * wrong here for two compounding reasons. MEQ's spaces are L2 on a
		 * GAUSS-LOBATTO basis, so the dof points include the element boundary --
		 * and an L2 field is DISCONTINUOUS across a source element boundary, so a
		 * target node landing on one has two values and gslib picks a side.
		 * Measured: transferring a converged solution onto a mesh whose space
		 * CONTAINS it, where the exact answer exists, came back 9%, 6% and 11%
		 * wrong in L2 at k = 1, 2, 3. Setting FindPointsGSLIB's L2 averaging did
		 * not fix it either -- it replaces the jump with a smoothed value that
		 * belongs to neither element, and measured 28%, 10%, 12%.
		 *
		 * A projection avoids the discontinuities entirely: GAUSS quadrature
		 * points are strictly interior, so every evaluation is unambiguous, and
		 * the projection of a function the target space can represent is that
		 * function. It also drops the requirement that the target be nodal at
		 * all, which the previous version had to check for and refuse.
		 *
		 * The mass matrix is element-local because the space is discontinuous:
		 * ndof x ndof, dense, one per element, and no global solve.
		 *
		 * WHAT IS LEFT after that is 1e-6 to 1e-9 relative on a transfer that has
		 * an exact answer, and it is NOT gslib's point search: tightening
		 * Setup()'s Newton tolerance from its default 1e-12 to 1e-15 changes not
		 * one digit. It is not chased further, because a guess reproduced to
		 * around 1e-8 relative is not the limiting factor in anything a starting
		 * point is used for.
		 */
		mfem::IntegrationRule const **rules = nullptr;
		(void)rules;

		std::vector<double> radius, z;
		std::vector<int> firstPoint( elements + 1, 0 );

		mfem::Vector node( 2 );
		for ( int e = 0; e < elements; ++e )
		{
			mfem::FiniteElement const *fe = targetSpace->GetFE( e );
			mfem::ElementTransformation *tr = targetMesh.GetElementTransformation( e );
			// 2k + 2 integrates the mass matrix of a degree-k space exactly and
			// leaves room for the source, which is the same degree on a mesh that
			// does not line up with this one.
			mfem::IntegrationRule const &ir =
				mfem::IntRules.Get( fe->GetGeomType(), 2*fe->GetOrder() + 2 );

			firstPoint[ e ] = static_cast<int>( radius.size() );
			for ( int q = 0; q < ir.GetNPoints(); ++q )
			{
				tr->Transform( ir.IntPoint( q ), node );
				radius.push_back( node( 0 ) );
				z.push_back( node( 1 ) );
			}
		}
		firstPoint[ elements ] = static_cast<int>( radius.size() );

		points = static_cast<int>( radius.size() );
		misses = 0;
		worstDist = 0.0;
		if ( points == 0 )
			return 0;

		// byNODES: every R, then every z. gslib's own default, passed explicitly
		// because the two orderings differ silently and a transposed point cloud
		// would give a plausible wrong answer rather than an error.
		mfem::Vector positions( 2*points );
		for ( int i = 0; i < points; ++i )
		{
			positions( i ) = radius[ i ];
			positions( points + i ) = z[ i ];
		}

		locator.locate( positions );

		mfem::Vector values;
		locator.interpolate( source, values );

		mfem::Array<unsigned int> const &code = locator.codes();
		mfem::Vector const &distance = locator.distances();

		/*
		 * THE THREE HOST READS THAT FOLLOW, MADE HOST READS EXPLICITLY.
		 *
		 * Everything below this line is host arithmetic -- CalcShape, a dense
		 * element mass matrix, DenseMatrixInverse -- and it reads `values`,
		 * `code` and `distance` through operator() and operator[], which are
		 * RAW: they neither sync nor invalidate. With an mfem::Device
		 * configured, meq::PointLocator runs its search as a kernel --
		 * PointLocatorDevice.cpp is one of the five translation units nvcc owns
		 * -- and leaves its output on the device, so those raw reads take a
		 * STALE HOST COPY. The hazard belongs to the LOCATOR being device-aware
		 * and not to any one implementation of it: gslib had it too, through
		 * fem/gslib.hpp's InterpolateOnDevice, and replacing gslib did not
		 * remove it.
		 *
		 * WHAT THAT COSTS IS NOT A CRASH, WHICH IS WHY IT IS WORTH THE COMMENT.
		 * Under CUDA it is silent: the projection is built from whatever the
		 * host buffer last held, so the warm start comes back wrong and the
		 * solve begins from a different iterate. Measured on
		 * examples/limited-tokamak.toml, ||r_0|| read 1.048568e-01 against the
		 * host's 2.411595e-01 at the SAME stored guess, and the bordered Newton
		 * then reported a singular Jacobian at iteration zero -- a failure that
		 * names the border and has nothing to do with it.
		 *
		 * mfem::Device( "debug" ) is what turns it into a diagnosis: it
		 * mprotects the host page, so the first of these reads is a named fault
		 * with a backtrace instead of a wrong number.
		 *
		 * SAME SPECIES AS THE FOUR SITES BEHIND M-79, and dormant for the same
		 * reason: a contract MEQ never had to honour until a Device existed.
		 * These calls are no-ops with none configured.
		 */
		countedHostRead( values, __LINE__ );
		countedHostRead( code, __LINE__ );
		countedHostRead( distance, __LINE__ );

		mfem::Array<int> elementDofs;
		mfem::DenseMatrix mass;
		mfem::Vector shape, rhs, coefficients;

		for ( int e = 0; e < elements; ++e )
		{
			mfem::FiniteElement const *fe = targetSpace->GetFE( e );
			mfem::ElementTransformation *tr = targetMesh.GetElementTransformation( e );
			mfem::IntegrationRule const &ir =
				mfem::IntRules.Get( fe->GetGeomType(), 2*fe->GetOrder() + 2 );

			int const ndof = fe->GetDof();
			mass.SetSize( ndof );
			mass = 0.0;
			rhs.SetSize( ndof );
			rhs = 0.0;
			shape.SetSize( ndof );

			for ( int q = 0; q < ir.GetNPoints(); ++q )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( q );
				tr->SetIntPoint( &ip );
				fe->CalcShape( ip, shape );

				double const weight = ip.weight*tr->Weight();
				mfem::AddMult_a_VVt( weight, shape, mass );

				int const point = firstPoint[ e ] + q;
				double value = 0.0;
				if ( code[ point ] == 2 )
				{
					// Nothing there to read. The fallback is the Dirichlet datum
					// in every caller MEQ has, which is what a cold start would
					// have had at that point anyway.
					++misses;
					value = fallback.Eval( *tr, ip );
				}
				else
				{
					value = values( point );
					worstDist = std::max( worstDist, double( distance( point ) ) );
				}

				rhs.Add( weight*value, shape );
			}

			mfem::DenseMatrixInverse inverse( mass );
			coefficients.SetSize( ndof );
			inverse.Mult( rhs, coefficients );

			targetSpace->GetElementDofs( e, elementDofs );
			target.SetSubVector( elementDofs, coefficients );
		}

		return misses;

	}

	int FieldTransfer::missed() const
	{
		return misses;
	}

	int FieldTransfer::queried() const
	{
		return points;
	}

	double FieldTransfer::worstDistance() const
	{
		return worstDist;
	}

}
