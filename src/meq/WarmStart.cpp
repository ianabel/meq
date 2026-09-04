#include "WarmStart.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace meq
{

	FieldTransfer::FieldTransfer( mfem::Mesh &sourceMeshIn )
		: sourceMesh( sourceMeshIn ), misses( 0 ), points( 0 ), worstDist( 0.0 )
	{
		if ( sourceMesh.Dimension() != 2 )
			throw std::invalid_argument( "meq::FieldTransfer: the source mesh must be two dimensional ( r, z )" );

		// FindPointsGSLIB::Setup() requires the mesh to carry an explicit nodal
		// grid function -- "Mesh nodes are required" -- and a mesh from
		// Mesh::MakeCartesian2D or read from a .mesh file has none until it is
		// asked for one. EnsureNodes() builds it at the mesh's own order, which
		// for MEQ's straight-sided triangles is degree 1 and exact. This is why
		// the constructor takes a non-const Mesh &: setting up the search
		// modifies the mesh, and pretending otherwise would need a copy of it.
		sourceMesh.EnsureNodes();
		finder.Setup( sourceMesh );

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
		 */
		finder.SetL2AvgType( mfem::FindPointsGSLIB::NONE );
	}

	FieldTransfer::~FieldTransfer()
	{
		finder.FreeData();
	}

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

		std::vector<double> r, z;
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

			firstPoint[ e ] = static_cast<int>( r.size() );
			for ( int q = 0; q < ir.GetNPoints(); ++q )
			{
				tr->Transform( ir.IntPoint( q ), node );
				r.push_back( node( 0 ) );
				z.push_back( node( 1 ) );
			}
		}
		firstPoint[ elements ] = static_cast<int>( r.size() );

		points = static_cast<int>( r.size() );
		misses = 0;
		worstDist = 0.0;
		if ( points == 0 )
			return 0;

		// byNODES: every r, then every z. gslib's own default, passed explicitly
		// because the two orderings differ silently and a transposed point cloud
		// would give a plausible wrong answer rather than an error.
		mfem::Vector positions( 2*points );
		for ( int i = 0; i < points; ++i )
		{
			positions( i ) = r[ i ];
			positions( points + i ) = z[ i ];
		}

		finder.FindPoints( positions, mfem::Ordering::byNODES );

		mfem::Vector values;
		finder.Interpolate( source, values );

		mfem::Array<unsigned int> const &code = finder.GetCode();
		mfem::Vector const &distance = finder.GetDist();

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
					worstDist = std::max( worstDist, distance( point ) );
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
