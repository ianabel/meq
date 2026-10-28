#include "DeviceBasis.hpp"

#include "Crossings.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace meq
{

	namespace
	{

		/// `p!` for `p <= triangleBasisMaxDegree`, exactly.
		double factorialOf( int n )
		{
			double total = 1.0;
			for ( int i = 2; i <= n; ++i )
				total *= static_cast<double>( i );
			return total;
		}

		/// `p!/( i! j! k! )` -- the Bernstein normalisation of one term.
		double multinomialOf( int i, int j, int k )
		{
			return factorialOf( i + j + k )
			       /( factorialOf( i )*factorialOf( j )*factorialOf( k ) );
		}

		/**
		 * The Bernstein polynomials of degree @a degree at one point.
		 *
		 * The host's expression of the basis `evaluateTriangleBasis()`
		 * evaluates, WITH the multinomial factor that build() then folds into
		 * the change of basis. Kept here rather than shared with the device
		 * form because the two are different functions -- one is the basis and
		 * the other is a contraction against it -- and writing the shared one
		 * would mean the builder and the kernel disagreeing about whether the
		 * factor is in.
		 */
		void bernsteinAt( int degree, double xi, double eta, mfem::Vector &out )
		{
			out.SetSize( triangleBasisDofs( degree ) );
			double const lambda[ 3 ] = { 1.0 - xi - eta, xi, eta };

			int at = 0;
			for ( int j = 0; j <= degree; ++j )
				for ( int k = 0; k <= degree - j; ++k, ++at )
				{
					int const i = degree - j - k;
					double term = multinomialOf( i, j, k );
					for ( int e = 0; e < i; ++e ) term *= lambda[ 0 ];
					for ( int e = 0; e < j; ++e ) term *= lambda[ 1 ];
					for ( int e = 0; e < k; ++e ) term *= lambda[ 2 ];
					out( at ) = term;
				}
		}

		/**
		 * A lattice for the conditioning check, out to @a reach.
		 *
		 * Deliberately NOT the element's own nodes: the change of basis is
		 * solved to be exact there, so measuring it there measures the linear
		 * solve's residual and not whether the two polynomials agree. The
		 * offsets keep every point off every equispaced and Gauss-Lobatto
		 * node.
		 *
		 * **AND IT GOES OUTSIDE THE TRIANGLE, WHICH IS THE WHOLE CORRECTION.**
		 * The first version sampled the interior alone -- and the interior is
		 * the one region this file's caller does not need it to be good in.
		 * `meq::CriticalPointFinder::rootInElement()` accepts a root up to
		 * `containment` outside its own element and lets its ITERATE wander
		 * two reference units out, and the barycentric monomials are
		 * conditioned by `| lambda0 |^p` out there. A budget taken where the
		 * answer cannot come from is a budget on nothing; see
		 * `extrapolationError()`.
		 */
		void checkLattice( double reach, std::vector<double> &xi,
		                   std::vector<double> &eta )
		{
			xi.clear();
			eta.clear();
			int constexpr steps = 13;
			double const span = 1.0 + 2.0*reach;
			for ( int a = 0; a <= steps; ++a )
				for ( int b = 0; a + b <= steps + 1; ++b )
				{
					double const x =
						-reach + span*( a + 0.37 )/( steps + 1.0 );
					double const y =
						-reach + span*( b + 0.19 )/( steps + 1.0 );
					// mfem::Geometry::CheckPoint( TRIANGLE, ip, reach ): the
					// same three inequalities the consumer accepts on.
					if ( x < -reach || y < -reach || x + y > 1.0 + reach )
						continue;
					xi.push_back( x );
					eta.push_back( y );
				}
		}

	}

	bool DeviceTriangleBasis::build( mfem::FiniteElementSpace const &space )
	{
		built = false;
		dofs = 0;
		order = -1;
		worstError = 0.0;
		worstOutside = 0.0;

		mfem::Mesh const *const mesh = space.GetMesh();
		if ( mesh == nullptr || mesh->GetNE() == 0 )
			return false;

		/*
		 * ONE GEOMETRY, ASKED RATHER THAN ASSUMED. An L2_FECollection hands
		 * every element of a geometry the same element object, so element 0's
		 * shape functions are every triangle's -- and on a MIXED mesh they are
		 * every triangle's and no quadrilateral's, which is a change of basis
		 * read off the first element and applied to elements it does not
		 * describe. `GetNumGeometries( dim )` is the one question that
		 * separates them.
		 */
		if ( mesh->GetNumGeometries( mesh->Dimension() ) != 1
		     || mesh->GetElementBaseGeometry( 0 ) != mfem::Geometry::TRIANGLE )
			return false;

		mfem::FiniteElement const *const fe = space.GetFE( 0 );
		if ( fe == nullptr )
			return false;

		/*
		 * A SCALAR ELEMENT, ASKED BEFORE CalcShape IS CALLED -- BECAUSE ON A
		 * VECTOR ONE IT DOES NOT RETURN.
		 *
		 * `mfem::VectorFiniteElement::CalcShape` is a PRIVATE override whose
		 * body is `mfem_error` (fe/fe_base.hpp), and this build sets
		 * `MFEM_USE_EXCEPTIONS`, so reaching it throws `ErrorException` out of
		 * a function whose whole contract is to return false. An
		 * `RT_FECollection( 0, 2 )` on triangles reaches it: `GetDof()` is 3
		 * and `GetOrder()` is 1, so it passes the dof-count test two lines
		 * down and arrives at the one call that cannot refuse.
		 *
		 * **A vdim-2 SCALAR space is fine and is what the flux is.**
		 * `FiniteElementSpace::GetFE()` hands back the scalar element whatever
		 * the vdim, so `meq::AxisSearch` builds one basis and indexes the two
		 * components itself. The two cases look alike and are not: this tests
		 * the ELEMENT's range type, not the space's vdim.
		 */
		if ( fe->GetRangeType() != mfem::FiniteElement::SCALAR )
			return false;

		int const degreeIn = fe->GetOrder();
		int const dofsIn = fe->GetDof();
		if ( degreeIn < 0 || degreeIn > triangleBasisMaxDegree
		     || dofsIn != triangleBasisDofs( degreeIn ) )
			return false;

		/*
		 * `M( a, b ) = B_b( node_a )` and `S( a, i ) = phi_i( node_a )`, then
		 * `T = M^{-1} S`. For a nodal basis `S` is the identity and this is a
		 * plain inverse -- but it is formed from `CalcShape` rather than
		 * assumed, because "the basis is nodal" is a property of a
		 * `BasisType` argument somebody may change, and getting it wrong
		 * silently would give a change of basis that is wrong by a permutation
		 * of the dofs.
		 */
		mfem::IntegrationRule const &elementNodes = fe->GetNodes();
		if ( elementNodes.GetNPoints() != dofsIn )
			return false;

		mfem::DenseMatrix bernstein( dofsIn, dofsIn );
		mfem::DenseMatrix shape( dofsIn, dofsIn );
		mfem::Vector row( dofsIn );

		for ( int a = 0; a < dofsIn; ++a )
		{
			mfem::IntegrationPoint const &node = elementNodes.IntPoint( a );

			bernsteinAt( degreeIn, node.x, node.y, row );
			for ( int b = 0; b < dofsIn; ++b )
				bernstein( a, b ) = row( b );

			fe->CalcShape( node, row );
			for ( int i = 0; i < dofsIn; ++i )
				shape( a, i ) = row( i );
		}

		/*
		 * THE CONDITIONING OF THE INVERSE, AS A REFUSAL AND NOT AS A COMMENT.
		 *
		 * `M` is the Bernstein Vandermonde at the element's own nodes and it
		 * is what gets inverted, so `|| M ||_inf || M^-1 ||_inf` is the factor
		 * by which a round-off in `CalcShape` reaches the change of basis.
		 * Measured on this machine it is 1.0, 3.0, 11.0 and 22.6 at degrees 1
		 * to 4 -- the Bernstein basis is a partition of unity, so every row of
		 * `M` sums to 1 and `|| M ||_inf` is exactly 1 at every degree -- and
		 * the threshold is three orders above the worst of those.
		 *
		 * An earlier draft computed these two numbers and then tested only
		 * that their product was finite, which `Invert()` has already
		 * guaranteed by aborting on a singular matrix. *A quantity computed
		 * and not conditioned on is a comment with arithmetic in it*, and the
		 * header promised a refusal on ill-conditioning that could not happen.
		 */
		double normBefore = 0.0;
		for ( int a = 0; a < dofsIn; ++a )
		{
			double sum = 0.0;
			for ( int b = 0; b < dofsIn; ++b )
				sum += std::fabs( bernstein( a, b ) );
			normBefore = std::max( normBefore, sum );
		}

		bernstein.Invert();

		double normAfter = 0.0;
		for ( int a = 0; a < dofsIn; ++a )
		{
			double sum = 0.0;
			for ( int b = 0; b < dofsIn; ++b )
				sum += std::fabs( bernstein( a, b ) );
			normAfter = std::max( normAfter, sum );
		}

		double constexpr conditionCeiling = 1.0e4;
		if ( !( normBefore*normAfter <= conditionCeiling ) )
			return false;

		mfem::DenseMatrix transform( dofsIn, dofsIn );
		mfem::Mult( bernstein, shape, transform );

		/*
		 * THE MULTINOMIAL FOLDED IN. `M` was built against the Bernstein
		 * polynomials because they are what is worth inverting against; the
		 * kernel evaluates plain powers, so the factor moves here, once, and
		 * never appears in a launch.
		 */
		{
			int at = 0;
			for ( int j = 0; j <= degreeIn; ++j )
				for ( int k = 0; k <= degreeIn - j; ++k, ++at )
				{
					double const scale =
						multinomialOf( degreeIn - j - k, j, k );
					for ( int i = 0; i < dofsIn; ++i )
						transform( at, i ) *= scale;
				}
		}

		/*
		 * AND THE CHECK IS THAT THE TWO POLYNOMIALS AGREE, NOT THAT THE SOLVE
		 * CONVERGED. Each shape function is evaluated both ways over a lattice
		 * that is off every node, and the worst RELATIVE disagreement is kept.
		 *
		 * **TWO REACHES, AND ONLY ONE OF THEM IS A BUDGET.** `containment` is
		 * how far outside its own element a root may be ACCEPTED, so the error
		 * there bounds the answer and is refused on. Two reference units is
		 * where `rootInElement()` ABANDONS, so the error there bounds the
		 * iterate's path and is reported instead -- a worse evaluation out
		 * there costs a few steps taking a different route to the same root,
		 * and only becomes a correctness question when it approaches the
		 * stopping target, which is what the ceiling on `triangleBasisMaxDegree`
		 * is set by.
		 */
		double constexpr containmentReach = 0.10;
		double constexpr abandonReach = 2.0;

		auto worstAgainstCalcShape = [ & ]( double reach,
		                                    bool relative ) -> double
		{
			std::vector<double> latticeXi, latticeEta;
			checkLattice( reach, latticeXi, latticeEta );

			std::vector<mfem::real_t> coefficients(
				static_cast<std::size_t>( dofsIn ) );
			mfem::Vector exact( dofsIn );
			double worst = 0.0;

			for ( std::size_t p = 0; p < latticeXi.size(); ++p )
			{
				mfem::IntegrationPoint probe;
				probe.Set2( latticeXi[ p ], latticeEta[ p ] );
				fe->CalcShape( probe, exact );

				double scale = 0.0;
				for ( int i = 0; i < dofsIn; ++i )
					// static_cast because mfem::Vector::operator() returns
					// real_t: under MFEM_PRECISION = single that is a float
					// and std::max refuses to deduce a common type with the
					// double accumulator. That refusal is the instrument --
					// it is how the sites needing conversion get found at
					// all, rather than by an audit.
					scale = std::max(
						scale, static_cast<double>( std::fabs( exact( i ) ) ) );
				if ( !( scale > 0.0 ) )
					scale = 1.0;

				for ( int i = 0; i < dofsIn; ++i )
				{
					// The i-th shape function is the change of basis applied
					// to the i-th unit dof vector, which is column i.
					for ( int b = 0; b < dofsIn; ++b )
						coefficients[ static_cast<std::size_t>( b ) ] =
							transform( b, i );

					mfem::real_t value = 0.0;
					evaluateTriangleBasis( degreeIn, coefficients.data(),
					                       latticeXi[ p ], latticeEta[ p ],
					                       value, nullptr, nullptr );
					double const gap =
						std::fabs( static_cast<double>( value ) - exact( i ) );
					worst = std::max( worst, relative ? gap/scale : gap );
				}
			}
			return worst;
		};

		/*
		 * RELATIVE IN THE BAND AND ABSOLUTE AT THE BOUND, WHICH IS NOT A
		 * SLIP -- THE TWO NUMBERS ANSWER DIFFERENT QUESTIONS.
		 *
		 * In the containment band the question is whether the two polynomials
		 * are the same function where a root may be ACCEPTED, and that is a
		 * relative statement: the shape functions have an O( 1 ) scale there
		 * and the answer carries it. At the abandon bound the question is
		 * whether the Newton's own stopping test can still fire, and that test
		 * is `residual <= tolerance * max | q |` -- an ABSOLUTE threshold
		 * against a global scale, so a relative figure normalised by shape
		 * functions that are themselves growing would hide exactly what it is
		 * being asked about. Measured that way the abandon-bound error reads
		 * 5.5e-11 at degree 4 against 1.6e-09 at 5 and 4.8e-07 at 6, which is
		 * what put the ceiling at 4.
		 */
		worstError = worstAgainstCalcShape( containmentReach, true );
		worstOutside = worstAgainstCalcShape( abandonReach, false );

		/*
		 * THE REFUSAL, AND ITS NUMBER IS A BUDGET RATHER THAN A HOPE.
		 *
		 * The caller is a search for a STATIONARY point, where M-220 measured
		 * that a field perturbation of `d` moves the located position by about
		 * `sqrt( d )` and the value by `d`. The driver prints six decimals of
		 * a position, so a reconstruction error of 1e-12 buys a position
		 * uncertainty of 1e-06 and is the loosest thing worth taking; 1e-10
		 * would put it at 1e-05, which is a printed digit. Anything above this
		 * falls back to the host rather than being evaluated more carefully,
		 * because there is no more careful mode -- the conditioning is the
		 * basis's.
		 */
		/*
		 * **AND THE BUDGET FOLLOWS real_t, WHICH THE FIRST DRAFT'S DID NOT.**
		 *
		 * A bare 1.0e-12 is a build-configuration bug wearing a numerical
		 * costume -- M-183's phrase, and the fifth constant of that kind MEQ
		 * has found. Under `MFEM_PRECISION = single` it sits **five orders
		 * below `eps`**, so `build()` would refuse every space and every
		 * caller would silently keep the host route: a correct answer by a
		 * path nobody chose, which is the worst way for a refusal to be wrong.
		 *
		 * `std::max` against the historical value is the construction M-184
		 * settled on and it does two things at once -- the double build is
		 * **bit-unchanged**, and `std::max`'s refusal to deduce a common type
		 * makes any site that mixes a literal with a `real_t` a COMPILE error
		 * rather than a silent one.
		 *
		 * 100 eps rather than a fixed figure, because what this measures is
		 * CONDITIONING -- how much round-off the change of basis amplifies --
		 * and that is a multiple of eps by construction. The measured worst is
		 * about 24 eps at degree 4, so the budget is four times the evidence
		 * in both formats.
		 */
		double const errorBudget = std::max(
			1.0e-12,
			100.0*static_cast<double>(
				std::numeric_limits<mfem::real_t>::epsilon() ) );
		if ( !( worstError <= errorBudget ) )
			return false;

		/*
		 * UPLOADED ONCE, AND `HostWrite()` BEFORE THE FILL BECAUSE M-211 AND
		 * M-226 BOTH SAY SO. This object is built rather than refreshed, so
		 * the fill is a constructor's and the defect those measurements
		 * describe cannot reach it -- the call is here anyway, because the
		 * cost is nothing and the next person to add a refresh should find the
		 * pattern already right.
		 */
		change.SetSize( dofsIn*dofsIn );
		change.UseDevice( true );
		mfem::real_t *const host = change.HostWrite();
		for ( int b = 0; b < dofsIn; ++b )
			for ( int i = 0; i < dofsIn; ++i )
				host[ b*dofsIn + i ] =
					static_cast<mfem::real_t>( transform( b, i ) );
		countedDeviceRead( change, __LINE__ );

		nodes.SetSize( 2*dofsIn );
		nodes.UseDevice( true );
		mfem::real_t *const nodeHost = nodes.HostWrite();
		for ( int a = 0; a < dofsIn; ++a )
		{
			mfem::IntegrationPoint const &node = elementNodes.IntPoint( a );
			nodeHost[ 2*a ]     = static_cast<mfem::real_t>( node.x );
			nodeHost[ 2*a + 1 ] = static_cast<mfem::real_t>( node.y );
		}
		countedDeviceRead( nodes, __LINE__ );

		dofs = dofsIn;
		order = degreeIn;
		built = true;
		return true;
	}

	bool DeviceTriangleBasis::valid() const
	{
		return built;
	}

	TriangleBasisView DeviceTriangleBasis::hostView() const
	{
		TriangleBasisView out{};
		out.degree = -1;
		if ( !built )
			return out;
		// Counted: a host read of a phase B upload is a DOWNLOAD, and leaving
		// it out of the instrument would make the one direction of the bus
		// this class contributes to invisible.
		countedHostRead( change, __LINE__ );
		countedHostRead( nodes, __LINE__ );
		out.degree = order;
		out.dofs = dofs;
		out.change = change.HostRead();
		out.nodes = nodes.HostRead();
		return out;
	}

	TriangleBasisView DeviceTriangleBasis::view() const
	{
		TriangleBasisView out{};
		// -1 rather than 0, because 0 is a legal degree and a caller testing
		// this must be able to tell a refusal from a constant field.
		out.degree = -1;
		if ( !built )
			return out;
		out.degree = order;
		out.dofs = dofs;
		out.change = change.Read();
		out.nodes = nodes.Read();
		return out;
	}

	int DeviceTriangleBasis::dofsPerElement() const
	{
		return dofs;
	}

	int DeviceTriangleBasis::degree() const
	{
		return order;
	}

	double DeviceTriangleBasis::reconstructionError() const
	{
		return worstError;
	}

	double DeviceTriangleBasis::extrapolationError() const
	{
		return worstOutside;
	}

	CoordinateField::CoordinateField( mfem::Mesh &meshIn )
		: mesh( &meshIn )
	{
		if ( mesh->GetNodes() != nullptr )
			return;

		/*
		 * A VERTEX FIELD BY THE SPACE'S OWN VERTEX DOFS, NOT A PROJECTION.
		 * Mesh::GetNodes( GridFunction & ) would project the identity
		 * element by element; on a degree-1 H1 space that reads each vertex
		 * once per element it touches and lands the same numbers. Asking
		 * the space which vdofs a vertex owns is O( vertices ) and assumes
		 * nothing about how it numbers them.
		 */
		int const dim = mesh->Dimension();
		int const spaceDimension = mesh->SpaceDimension();
		collection = std::make_unique<mfem::H1_FECollection>( 1, dim );
		space = std::make_unique<mfem::FiniteElementSpace>(
			mesh, collection.get(), spaceDimension, mfem::Ordering::byNODES );
		vertices = std::make_unique<mfem::GridFunction>( space.get() );
		vertices->UseDevice( true );
		mfem::real_t *const host = vertices->HostWrite();
		mfem::Array<int> vdofs;
		for ( int v = 0; v < mesh->GetNV(); ++v )
		{
			space->GetVertexVDofs( v, vdofs );
			mfem::real_t const *point = mesh->GetVertex( v );
			for ( int c = 0; c < spaceDimension; ++c )
				host[ vdofs[ c ] ] = point[ c ];
		}
	}

	mfem::GridFunction &CoordinateField::field()
	{
		if ( mfem::GridFunction *nodes = mesh->GetNodes() )
			return *nodes;
		return *vertices;
	}

	mfem::GeometricFactors const &CoordinateField::factors( mfem::IntegrationRule const &rule,
	                                                         int flags )
	{
		mfem::GridFunction const &coordinates = field();
		for ( Cached const &entry : cache )
			if ( entry.field == &coordinates && entry.factors->IntRule == &rule
			     && ( entry.factors->computed_factors & flags ) == flags )
				return *entry.factors;

		cache.push_back( Cached{ &coordinates,
			std::make_unique<mfem::GeometricFactors>( coordinates, rule, flags ) } );
		// Back on the host for the host readers after this; see the class.
		coordinates.HostRead();
		return *cache.back().factors;
	}

	bool affineTrianglesOf( mfem::Mesh &mesh, mfem::Array<AffineTriangle> &out,
	                        CoordinateField *coordinates )
	{
		int const elements = mesh.GetNE();
		out.SetSize( elements );
		out.UseDevice( true );

		/*
		 * UNDER A DEVICE, ON A MESH OF STRAIGHT TRIANGLES, THE MAPS ARE BUILT
		 * THERE. The host loop below asks MFEM's transformation of every
		 * element twice, on one core under a Device: 90% of the driver's
		 * post-solve axis check on fixed-f-diiid at 102,840 triangles. MFEM's
		 * geometric factors at the reference corner are the same map -- the
		 * position there is the origin, and the Jacobian of a degree-1 map is
		 * the constant one the host checks for at the centre -- computed by
		 * MFEM's own kernel. A degenerate element still comes back with
		 * `affine` false, as the host loop leaves it.
		 *
		 * The factors come from a CoordinateField, never from the mesh, so a
		 * mesh with no nodes keeps none (see the class); and they are cached by
		 * the RULE'S ADDRESS, which is why the rule is a static.
		 */
		mfem::FiniteElementSpace const *nodal = mesh.GetNodalFESpace();
		if ( mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) && elements > 0
		     && mesh.Dimension() == 2 && mesh.SpaceDimension() == 2
		     && mesh.GetNumGeometries( 2 ) == 1
		     && mesh.GetElementBaseGeometry( 0 ) == mfem::Geometry::TRIANGLE
		     && ( nodal == nullptr || nodal->GetMaxElementOrder() == 1 ) )
		{
			static mfem::IntegrationRule const corner = []()
			{
				mfem::IntegrationRule rule( 1 );
				rule.IntPoint( 0 ).Set2w( 0.0, 0.0, 0.0 );
				return rule;
			}();
			std::unique_ptr<CoordinateField> local;
			if ( coordinates == nullptr )
			{
				local = std::make_unique<CoordinateField>( mesh );
				coordinates = local.get();
			}
			mfem::GeometricFactors const &factors = coordinates->factors(
				corner, mfem::GeometricFactors::COORDINATES
				        | mfem::GeometricFactors::JACOBIANS );
			if ( affineTrianglesKernel( factors.X, factors.J, elements, out ) )
				return true;
		}
		AffineTriangle *const host = out.HostWrite();

		bool everyOne = true;

		/*
		 * MFEM'S OWN MAP RATHER THAN THE VERTICES. Reading `v1 - v0` and
		 * `v2 - v0` off the vertex array reproduces the affine map only if the
		 * element's vertex order is the reference triangle's, which is true
		 * and is an assumption this file would then be making on MFEM's
		 * behalf. Asking the transformation costs one call per element, once
		 * per mesh, and is the same answer `rootInElement()` gets.
		 *
		 * A function-local `thread_local` to avoid constructing a
		 * `DenseMatrix` per element, NOT because of the shared-scratch trap:
		 * the two-argument overload used here is the REENTRANT one and writes
		 * only into what it is given. The hazard CLAUDE.md records is the
		 * one-argument form, which this file never calls -- said because the
		 * first draft of this comment cited it and was wrong about which
		 * overload the rule applies to.
		 */
		thread_local mfem::IsoparametricTransformation scratch;
		mfem::IntegrationPoint corner;
		corner.Set2( 0.0, 0.0 );
		mfem::IntegrationPoint centre;
		centre.Set2( 1.0/3.0, 1.0/3.0 );
		mfem::Vector position( 2 );

		for ( int e = 0; e < elements; ++e )
		{
			AffineTriangle map{};
			map.affine = false;

			if ( mesh.GetElementBaseGeometry( e ) != mfem::Geometry::TRIANGLE )
			{
				host[ e ] = map;
				everyOne = false;
				continue;
			}

			mesh.GetElementTransformation( e, &scratch );

			scratch.SetIntPoint( &corner );
			scratch.Transform( corner, position );
			map.origin[ 0 ] = static_cast<mfem::real_t>( position( 0 ) );
			map.origin[ 1 ] = static_cast<mfem::real_t>( position( 1 ) );

			mfem::DenseMatrix const &jacobian = scratch.Jacobian();
			if ( jacobian.Height() != 2 || jacobian.Width() != 2 )
			{
				host[ e ] = map;
				everyOne = false;
				continue;
			}
			map.jacobian[ 0 ] = static_cast<mfem::real_t>( jacobian( 0, 0 ) );
			map.jacobian[ 1 ] = static_cast<mfem::real_t>( jacobian( 0, 1 ) );
			map.jacobian[ 2 ] = static_cast<mfem::real_t>( jacobian( 1, 0 ) );
			map.jacobian[ 3 ] = static_cast<mfem::real_t>( jacobian( 1, 1 ) );

			/*
			 * AND WHETHER IT IS AFFINE AT ALL IS MEASURED, NOT INFERRED FROM
			 * `Mesh::GetNodes() == nullptr`. CLAUDE.md records a mesh
			 * ACQUIRING nodes without becoming curved -- `EnsureNodes()` on
			 * behalf of a point locator -- and `meq::GridSampler` paying two
			 * orders of magnitude for reading that pointer as the question.
			 * The Jacobian of an affine map is constant, so comparing it at
			 * two points is the question itself.
			 */
			scratch.SetIntPoint( &centre );
			mfem::DenseMatrix const &elsewhere = scratch.Jacobian();
			double drift = 0.0;
			double size = 0.0;
			for ( int a = 0; a < 2; ++a )
				for ( int b = 0; b < 2; ++b )
				{
					drift = std::max( drift,
						std::fabs( elsewhere( a, b )
						           - static_cast<double>( map.jacobian[ 2*a + b ] ) ) );
					size = std::max(
						size,
						static_cast<double>( std::fabs( elsewhere( a, b ) ) ) );
				}
			/*
			 * AND THE AFFINENESS TEST FOLLOWS real_t FOR THE SAME REASON. In
			 * double the drift of an affine map's Jacobian is EXACTLY zero --
			 * measured, bitwise, at three well-separated reference points --
			 * so any threshold works there. In single it is a few eps of the
			 * map's own scale, and a bare 1.0e-12 would mark every straight
			 * triangle in the mesh curved and refuse the route on a mesh that
			 * has nothing wrong with it.
			 */
			double const affineDrift = std::max(
				1.0e-12,
				10.0*static_cast<double>(
					std::numeric_limits<mfem::real_t>::epsilon() ) );
			if ( !( size > 0.0 ) || drift > affineDrift*size )
			{
				host[ e ] = map;
				everyOne = false;
				continue;
			}

			double const determinant =
				static_cast<double>( map.jacobian[ 0 ] )*map.jacobian[ 3 ]
				- static_cast<double>( map.jacobian[ 1 ] )*map.jacobian[ 2 ];
			if ( !( std::fabs( determinant ) > 0.0 ) )
			{
				host[ e ] = map;
				everyOne = false;
				continue;
			}

			map.determinant = static_cast<mfem::real_t>( determinant );
			mfem::real_t const inverse = mfem::real_t( 1.0 )
				/static_cast<mfem::real_t>( determinant );
			map.inverseJacobian[ 0 ] =  map.jacobian[ 3 ]*inverse;
			map.inverseJacobian[ 1 ] = -map.jacobian[ 1 ]*inverse;
			map.inverseJacobian[ 2 ] = -map.jacobian[ 2 ]*inverse;
			map.inverseJacobian[ 3 ] =  map.jacobian[ 0 ]*inverse;
			map.affine = true;

			host[ e ] = map;
		}

		countedDeviceRead( out, __LINE__ );
		return everyOne;
	}

}
