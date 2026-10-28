// Unit tests for meq::DeviceBasis -- evaluating a grid function at an
// arbitrary reference point without MFEM's host virtuals.
//
// THIS FILE ASSERTS A TOLERANCE, AND THAT IS A DELIBERATE DIFFERENCE FROM
// tests/unit/DeviceRowsTests.cpp, WHICH ASSERTS EQUALITY TO THE BIT.
// meq::DeviceRows reproduces a host loop that exists, so there is an ordering
// to match and matching it is the contract. A polynomial re-expressed in
// another basis has no such host loop: the change of basis is a dense matrix
// multiply, so it agrees with mfem::FiniteElement::CalcShape to round-off and
// cannot agree to the bit. Writing `==` here would be asserting something that
// is not true rather than something that is demanding.
//
// SO THE TOLERANCE IS MEASURED AND NOT CHOSEN, which is the only way a number
// like this means anything. The figures these cases print are the evidence,
// the assertion sits a decade above the worst of them, and
// meq::DeviceTriangleBasis::build() REFUSES rather than returning a basis that
// does not reach its own budget -- so a degree this machine cannot condition
// falls back to the host instead of being evaluated badly.
//
// AND IT CONFIGURES AN mfem::Device, WHICH THE FIRST VERSION DID NOT AND WHICH
// IS WHY THE FIRST VERSION WAS GREEN ON TWO REAL DEFECTS. Both siblings --
// DeviceRowsTests.cpp and DeviceQuadratureTests.cpp -- open with one, and
// without it meq::DeviceTriangleBasis::view() hands back host pointers and
// meq::affineTrianglesOf()'s output is never device-valid, so a test that
// dereferences either on the host passes for a reason that evaporates the
// moment anybody runs the code it is testing. Measured before the repair: a
// host read of the affine maps FAULTS under --device debug, and
// `view().change` SEGFAULTS under --device cuda while surviving debug -- which
// is M-206's recorded limit, the debug backend catching a raw access to the
// wrong LOCATION and being blind to a pointer into the wrong address space.
//
// WHAT THIS FILE DOES NOT COVER, SAID PLAINLY. It exercises the evaluator the
// way the HOST compiler builds it. The arithmetic is MFEM_HOST_DEVICE and
// lives in the header precisely so that one expression is compiled twice, but
// only tests/convergence/DeviceParity.cpp's --device cuda arm puts it through
// nvcc. That is the same arrangement meq::PointLocator's locateAffine() is in
// and it is recorded here so nobody reads a green run of this binary as
// evidence about a kernel.

#define BOOST_TEST_MODULE DeviceBasisTests
#include <boost/test/included/unit_test.hpp>

#include "meq/Crossings.hpp"
#include "meq/DeviceBasis.hpp"

#include "mfem.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace
{
	/**
	 * `real_t`'s machine epsilon, and the unit every tolerance here is in.
	 *
	 * **THE ASSERTIONS ARE IN EPSILONS BECAUSE THE CLAIM IS ABOUT
	 * CONDITIONING, AND CONDITIONING IS A MULTIPLE OF EPS BY CONSTRUCTION.**
	 * The first version of this file asserted fixed figures -- 1e-11, 1e-12,
	 * 1e-14 -- which are true in double and are ARTEFACTS of double there, so
	 * every case went red against `install-single` while measuring the same
	 * thing. M-219 records eleven tests repaired this way and the lesson each
	 * time: *asserting the claim directly is often a STRONGER test.* It is
	 * here. "The change of basis costs about twenty epsilons" is a statement
	 * about the algorithm; "it is below 1e-12" is a statement about `double`,
	 * and a basis three decades worse would have passed it.
	 *
	 * Measured, and the two builds agree to within a factor of two on every
	 * line: reconstruction **24 eps** in double against **23** in single, the
	 * value **4 and 7**, the gradient **19 and 15**, the affine map **2 and
	 * 2**. That agreement is the evidence the budget is a conditioning bound
	 * and not a format one.
	 */
	double epsilon()
	{
		return static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() );
	}

	/*
	 * ONE Device FOR THE BINARY, AND IT MUST BE A FUNCTION-LOCAL STATIC.
	 *
	 * mfem::Device::Get() returns a static member defined in
	 * general/device.cpp and Setup() opens MFEM_VERIFY( ngpu == -1 ); a
	 * namespace-scope Device here may be constructed before that singleton's
	 * dynamic initialisation, leaving `ngpu` at its zero-initialised 0 and
	 * producing "the mfem::Device is already configured!" on the first and
	 * only Device in the process. tests/unit/DeviceRowsTests.cpp carries the
	 * longer form of this; it is the static-initialisation-order fiasco in a
	 * device costume.
	 *
	 * "debug" by default and cuda on request:
	 *
	 *     MEQ_TEST_DEVICE=cuda build/tests/DeviceBasisTests
	 *
	 * The two answer different questions. `debug` has device memory semantics
	 * with HOST arithmetic, so it proves the bookkeeping and needs no card.
	 * **Only `cuda` can catch a pointer into the wrong address space**, which
	 * is what view() hands out -- so the cuda run is part of this file's
	 * acceptance even though CI cannot take it.
	 */
	mfem::Device &deviceOnce()
	{
		char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
		static mfem::Device device( requested ? requested : "debug" );
		return device;
	}

	/// Triangles that are genuinely differently shaped, which the first
	/// version of this file claimed and did not have: MakeCartesian2D gives
	/// twelve congruent triangles in two mirrored classes, every one of them
	/// with `jacobian[ 1 ] == 0` -- so an error confined to that one slot of
	/// the affine map would have been invisible in every case here. Shearing
	/// the vertices closes both gaps.
	mfem::Mesh triangleMesh()
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			3, 2, mfem::Element::TRIANGLE, false, 1.7, 0.9 );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			// real_t, not double: mfem::Mesh::GetVertex() returns real_t*,
			// which is float under MFEM_PRECISION = single.
			mfem::real_t *const point = mesh.GetVertex( v );
			mfem::real_t const x = point[ 0 ];
			mfem::real_t const y = point[ 1 ];
			point[ 0 ] = x + 0.21*y;
			point[ 1 ] = y + 0.13*x*x;
		}
		return mesh;
	}

	/// Reference points, including several OUTSIDE the reference triangle --
	/// which is not padding: meq::CriticalPointFinder::rootInElement() lets
	/// its Newton iterate leave by up to two units in REFERENCE coordinates,
	/// where the reference element's diameter is one, and accepts a root up to
	/// `containment` outside. An evaluator that is only right inside would be
	/// silently wrong exactly where that search spends its hardest steps.
	std::vector<std::pair<double, double>> probePoints()
	{
		return {
			{ 0.10, 0.10 }, { 0.70, 0.20 }, { 0.20, 0.70 },
			{ 0.33, 0.33 }, { 0.01, 0.98 }, { 0.49, 0.49 },
			{ -0.08, 0.20 }, { 0.20, -0.08 }, { 1.05, 0.02 },
			{ -0.40, 0.20 }, { 0.60, 0.90 }, { -0.25, -0.25 },
			{ 1.80, -0.70 }
		};
	}

	/// A dof vector with structure and several orders of magnitude in it.
	void fillDofs( mfem::GridFunction &gf, int seed )
	{
		for ( int i = 0; i < gf.Size(); ++i )
		{
			double const t = ( i*37 + seed*11 )%101/101.0;
			gf( i ) = std::sin( 6.0*t + seed )*std::pow( 10.0, 2.0*t - 1.0 );
		}
	}
}

// THE CHANGE OF BASIS REPRODUCES CalcShape, AT EVERY DEGREE IT ACCEPTS.
//
// The sweep goes to meq::triangleBasisMaxDegree rather than to MEQ's own
// k = 4, because the boundary of what build() accepts is the one place a
// conditioning budget can be wrong -- an accepted degree nothing exercises is
// a refusal nobody has checked the other side of.
//
// The claim is about a FUNCTION and not about a set of sample points, so the
// comparison is made against MFEM's own evaluation of the same grid function
// rather than against a stored table -- and at points off every node, since
// the construction is solved to be exact at the nodes and measuring it there
// would measure the linear solve's residual instead.
BOOST_AUTO_TEST_CASE( theChangeOfBasisReproducesMfemsOwnEvaluation )
{
	deviceOnce();
	mfem::Mesh mesh = triangleMesh();

	for ( int degree = 1; degree <= meq::triangleBasisMaxDegree; ++degree )
	{
		mfem::L2_FECollection collection( degree, 2,
		                                  mfem::BasisType::GaussLobatto );
		mfem::FiniteElementSpace space( &mesh, &collection );

		meq::DeviceTriangleBasis basis;
		BOOST_TEST_REQUIRE( basis.build( space ) );
		BOOST_TEST_REQUIRE( basis.degree() == degree );
		BOOST_TEST_REQUIRE( basis.dofsPerElement()
		                    == meq::triangleBasisDofs( degree ) );

		mfem::GridFunction field( &space );
		fillDofs( field, degree );

		int const dofs = basis.dofsPerElement();
		std::vector<mfem::real_t> local(
			static_cast<std::size_t>( dofs ) );
		std::vector<mfem::real_t> coefficients(
			static_cast<std::size_t>( dofs ) );

		// hostView(), not view(): this loop dereferences the pointers, and
		// view()'s are Read()'s -- device addresses under any Device.
		meq::TriangleBasisView const host = basis.hostView();
		BOOST_TEST_REQUIRE( host.change != nullptr );
		BOOST_TEST_REQUIRE( host.nodes != nullptr );

		double worstValue = 0.0;
		double scaleValue = 0.0;

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			// Blocked, because an L2 space lays its element dofs out that
			// way -- asserted rather than assumed, since the kernel indexes
			// this way and nothing else would catch a change.
			mfem::Array<int> elementDofs;
			space.GetElementDofs( e, elementDofs );
			BOOST_TEST_REQUIRE( elementDofs.Size() == dofs );
			for ( int i = 0; i < dofs; ++i )
				BOOST_TEST_REQUIRE( elementDofs[ i ] == e*dofs + i );

			for ( int i = 0; i < dofs; ++i )
				local[ static_cast<std::size_t>( i ) ] = field( e*dofs + i );
			meq::triangleBasisCoefficients( host, local.data(),
			                                coefficients.data() );

			for ( auto const &probe : probePoints() )
			{
				mfem::IntegrationPoint point;
				point.Set2( probe.first, probe.second );

				double const expected = field.GetValue( e, point );

				mfem::real_t got = 0.0;
				meq::evaluateTriangleBasis( degree, coefficients.data(),
				                            probe.first, probe.second,
				                            got, nullptr, nullptr );

				worstValue = std::max( worstValue,
				                       std::fabs( got - expected ) );
				scaleValue = std::max( scaleValue, std::fabs( expected ) );
			}
		}

		double const relative = worstValue/( scaleValue > 0.0 ? scaleValue : 1.0 );
		std::printf( "  degree %d: reconstruction %.3e (%.1f eps) in the "
		             "containment band, %.3e at the abandon bound; worst "
		             "value %.3e absolute, %.3e relative (%.1f eps)\n",
		             degree, basis.reconstructionError(),
		             basis.reconstructionError()/epsilon(),
		             basis.extrapolationError(), worstValue, relative,
		             relative/epsilon() );

		BOOST_TEST( relative < 100.0*epsilon() );
		BOOST_TEST( basis.reconstructionError() < 100.0*epsilon() );
	}
}

// AND THE GRADIENT, WHICH IS A SEPARATE CLAIM.
//
// evaluateTriangleBasis() differentiates the monomials analytically rather
// than differencing them, so a sign error in the lambda0 chain rule -- the one
// term whose derivative is -1 in BOTH directions -- would leave every value
// right and every gradient wrong. MFEM's GetGradient() is the reference, and
// because it is a PHYSICAL gradient this case also checks the stored inverse
// Jacobian, which the forward map alone cannot.
BOOST_AUTO_TEST_CASE( theReferenceGradientAgreesWithMfemsPhysicalOne )
{
	deviceOnce();
	mfem::Mesh mesh = triangleMesh();

	for ( int degree = 1; degree <= meq::triangleBasisMaxDegree; ++degree )
	{
		mfem::L2_FECollection collection( degree, 2,
		                                  mfem::BasisType::GaussLobatto );
		mfem::FiniteElementSpace space( &mesh, &collection );

		meq::DeviceTriangleBasis basis;
		BOOST_TEST_REQUIRE( basis.build( space ) );

		mfem::Array<meq::AffineTriangle> maps;
		BOOST_TEST_REQUIRE( meq::affineTrianglesOf( mesh, maps ) );
		// THE FUNNEL, AND IT IS NOT DECORATION. affineTrianglesOf() leaves its
		// output DEVICE-VALID -- it is a phase B upload -- and
		// mfem::Array<T>::operator[] is as raw an accessor as
		// mfem::Vector::operator(). Without this the loop below faults under
		// --device debug, which is how it was found.
		meq::countedHostRead( maps, __LINE__ );

		mfem::GridFunction field( &space );
		fillDofs( field, degree + 7 );

		int const dofs = basis.dofsPerElement();
		std::vector<mfem::real_t> local( static_cast<std::size_t>( dofs ) );
		std::vector<mfem::real_t> coefficients(
			static_cast<std::size_t>( dofs ) );
		meq::TriangleBasisView const host = basis.hostView();

		double worst = 0.0;
		double scale = 0.0;

		mfem::IsoparametricTransformation transformation;
		mfem::Vector physical( 2 );

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			for ( int i = 0; i < dofs; ++i )
				local[ static_cast<std::size_t>( i ) ] = field( e*dofs + i );
			meq::triangleBasisCoefficients( host, local.data(),
			                                coefficients.data() );

			meq::AffineTriangle const &map = maps[ e ];
			BOOST_TEST_REQUIRE( map.affine );

			mesh.GetElementTransformation( e, &transformation );

			for ( auto const &probe : probePoints() )
			{
				mfem::IntegrationPoint point;
				point.Set2( probe.first, probe.second );
				transformation.SetIntPoint( &point );
				field.GetGradient( transformation, physical );

				mfem::real_t value = 0.0, dXi = 0.0, dEta = 0.0;
				meq::evaluateTriangleBasis( degree, coefficients.data(),
				                            probe.first, probe.second,
				                            value, &dXi, &dEta );

				// ( d/dr, d/dz ) = J^{-T} ( d/dxi, d/deta ).
				double const dR = map.inverseJacobian[ 0 ]*dXi
				                  + map.inverseJacobian[ 2 ]*dEta;
				double const dZ = map.inverseJacobian[ 1 ]*dXi
				                  + map.inverseJacobian[ 3 ]*dEta;

				// static_cast throughout: mfem::Vector::operator() is real_t,
				// and std::max will not deduce a common type with a double
				// accumulator under MFEM_PRECISION = single.
				worst = std::max( worst, std::fabs(
					dR - static_cast<double>( physical( 0 ) ) ) );
				worst = std::max( worst, std::fabs(
					dZ - static_cast<double>( physical( 1 ) ) ) );
				scale = std::max( scale, std::fabs(
					static_cast<double>( physical( 0 ) ) ) );
				scale = std::max( scale, std::fabs(
					static_cast<double>( physical( 1 ) ) ) );
			}
		}

		double const relative = worst/( scale > 0.0 ? scale : 1.0 );
		std::printf( "  degree %d: worst gradient %.3e absolute, "
		             "%.3e relative (%.1f eps)\n",
		             degree, worst, relative, relative/epsilon() );
		BOOST_TEST( relative < 100.0*epsilon() );
	}
}

// THE ELEMENT'S OWN NODES RIDE ALONG WITH THE CHANGE OF BASIS, AND THEY ARE
// WHAT THE SEARCH SEEDS FROM.
//
// meq::CriticalPointFinder::elementSeeds() starts its second seed at the flux
// node carrying the smallest | q |, so a kernel without these would have to
// seed somewhere else and would be a different search. They are uploaded
// beside the change of basis because they are a constant of the same element;
// nothing else asserts they are the right points.
BOOST_AUTO_TEST_CASE( theViewCarriesTheElementsOwnNodes )
{
	deviceOnce();
	mfem::Mesh mesh = triangleMesh();

	for ( int degree = 1; degree <= meq::triangleBasisMaxDegree; ++degree )
	{
		mfem::L2_FECollection collection( degree, 2,
		                                  mfem::BasisType::GaussLobatto );
		mfem::FiniteElementSpace space( &mesh, &collection );

		meq::DeviceTriangleBasis basis;
		BOOST_TEST_REQUIRE( basis.build( space ) );
		meq::TriangleBasisView const host = basis.hostView();

		mfem::FiniteElement const *const fe = space.GetFE( 0 );
		mfem::IntegrationRule const &nodes = fe->GetNodes();
		BOOST_TEST_REQUIRE( nodes.GetNPoints() == host.dofs );

		for ( int i = 0; i < host.dofs; ++i )
		{
			BOOST_TEST( host.nodes[ 2*i ] == nodes.IntPoint( i ).x );
			BOOST_TEST( host.nodes[ 2*i + 1 ] == nodes.IntPoint( i ).y );
		}
	}
}

// THE AFFINE MAP IS MFEM'S OWN, AND THE ROUND TRIP SAYS SO.
//
// A reference point pushed forward and the physical point pulled back must
// return to where they started. The inverse is STORED rather than solved for,
// so an error in forming it would never be caught by the forward direction
// alone -- and on this sheared mesh no entry of the Jacobian is zero, so a
// transposition or a single wrong slot has nowhere to hide.
BOOST_AUTO_TEST_CASE( theAffineMapRoundTripsAgainstMfem )
{
	deviceOnce();
	mfem::Mesh mesh = triangleMesh();

	mfem::Array<meq::AffineTriangle> maps;
	bool const hadNodes = mesh.GetNodes() != nullptr;
	BOOST_TEST_REQUIRE( meq::affineTrianglesOf( mesh, maps ) );
	meq::countedHostRead( maps, __LINE__ );
	BOOST_TEST_REQUIRE( maps.Size() == mesh.GetNE() );
	// AND THE MESH IS AS IT WAS. Mesh::GetGeometricFactors() would have
	// given a mesh without nodes a nodal field, which moves FindPoints() by
	// round-off and a free-boundary solve by far more (meq::CoordinateField).
	BOOST_TEST( ( mesh.GetNodes() != nullptr ) == hadNodes,
	            "affineTrianglesOf() gave the mesh nodes" );

	mfem::IsoparametricTransformation transformation;
	mfem::Vector physical( 2 );
	double worstForward = 0.0;
	double worstRoundTrip = 0.0;
	double smallestOffDiagonal = 1.0e300;

	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		meq::AffineTriangle const &map = maps[ e ];
		BOOST_TEST_REQUIRE( map.affine );
		mesh.GetElementTransformation( e, &transformation );

		smallestOffDiagonal = std::min(
			smallestOffDiagonal,
			static_cast<double>( std::min( std::fabs( map.jacobian[ 1 ] ),
			                               std::fabs( map.jacobian[ 2 ] ) ) ) );

		for ( auto const &probe : probePoints() )
		{
			mfem::IntegrationPoint point;
			point.Set2( probe.first, probe.second );
			transformation.Transform( point, physical );

			mfem::real_t radius = 0.0, z = 0.0;
			meq::affineTransform( map, probe.first, probe.second, radius, z );
			worstForward = std::max( worstForward,
				static_cast<double>(
					std::max( std::fabs( radius - physical( 0 ) ),
					          std::fabs( z - physical( 1 ) ) ) ) );

			double const dr = radius - map.origin[ 0 ];
			double const dz = z - map.origin[ 1 ];
			double const xi = map.inverseJacobian[ 0 ]*dr
			                  + map.inverseJacobian[ 1 ]*dz;
			double const eta = map.inverseJacobian[ 2 ]*dr
			                   + map.inverseJacobian[ 3 ]*dz;
			worstRoundTrip = std::max( worstRoundTrip,
				std::max( std::fabs( xi - probe.first ),
				          std::fabs( eta - probe.second ) ) );
		}
	}

	std::printf( "  affine map: forward %.3e (%.1f eps), round trip %.3e "
	             "(%.1f eps), smallest off-diagonal %.3e\n",
	             worstForward, worstForward/epsilon(),
	             worstRoundTrip, worstRoundTrip/epsilon(),
	             smallestOffDiagonal );
	// The mesh must actually exercise both off-diagonal slots, or the two
	// assertions below are weaker than they look.
	BOOST_TEST( smallestOffDiagonal > 1.0e-3 );
	BOOST_TEST( worstForward < 50.0*epsilon() );
	BOOST_TEST( worstRoundTrip < 50.0*epsilon() );
}

// AND IT REFUSES WHAT IT CANNOT DO, RATHER THAN APPROXIMATING IT.
//
// The refusals are the whole of what makes a caller's fallback safe, and each
// one here is a case some mesh or space in the tree actually presents: a
// quadrilateral guess mesh (machine-f-diiid ships one), a curved mesh, a
// vector element, and a degree above the conditioning ceiling.
BOOST_AUTO_TEST_CASE( itRefusesEveryGeometryAndSpaceItCannotTake )
{
	deviceOnce();

	mfem::Mesh quads = mfem::Mesh::MakeCartesian2D(
		2, 2, mfem::Element::QUADRILATERAL );
	mfem::L2_FECollection quadCollection( 2, 2, mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace quadSpace( &quads, &quadCollection );

	meq::DeviceTriangleBasis onQuads;
	BOOST_TEST( !onQuads.build( quadSpace ) );
	BOOST_TEST( !onQuads.valid() );
	BOOST_TEST( onQuads.view().change == nullptr );
	BOOST_TEST( onQuads.view().degree == -1 );

	mfem::Array<meq::AffineTriangle> quadMaps;
	BOOST_TEST( !meq::affineTrianglesOf( quads, quadMaps ) );
	meq::countedHostRead( quadMaps, __LINE__ );
	BOOST_TEST_REQUIRE( quadMaps.Size() == quads.GetNE() );
	for ( int e = 0; e < quadMaps.Size(); ++e )
		BOOST_TEST( !quadMaps[ e ].affine );

	mfem::Mesh mesh = triangleMesh();

	// A degree above the ceiling. The ceiling is a MEASURED conditioning
	// bound, not an array size -- see meq::triangleBasisMaxDegree.
	mfem::L2_FECollection tall( meq::triangleBasisMaxDegree + 1, 2,
	                            mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace tallSpace( &mesh, &tall );
	meq::DeviceTriangleBasis tooTall;
	BOOST_TEST( !tooTall.build( tallSpace ) );

	// A VECTOR element, which refuses rather than throwing.
	// mfem::VectorFiniteElement::CalcShape is a private override whose body is
	// mfem_error, and this build sets MFEM_USE_EXCEPTIONS -- so without the
	// range-type test build() would throw out of a function whose contract is
	// to return false. RT_FECollection( 0, 2 ) reaches it: GetDof() is 3 and
	// GetOrder() is 1, which passes the dof-count test.
	mfem::RT_FECollection raviart( 0, 2 );
	mfem::FiniteElementSpace raviartSpace( &mesh, &raviart );
	meq::DeviceTriangleBasis onVectors;
	BOOST_TEST( !onVectors.build( raviartSpace ) );

	// A CURVED mesh: the change of basis is still right, because it is about
	// the reference element -- but the affine map is not, and that is the half
	// the caller refuses on. Asserted separately so the two cannot be confused.
	mfem::Mesh curved = triangleMesh();
	curved.SetCurvature( 2 );
	{
		mfem::GridFunction *const nodes = curved.GetNodes();
		BOOST_TEST_REQUIRE( nodes != nullptr );
		for ( int i = 0; i < nodes->Size(); ++i )
			( *nodes )( i ) += 0.04*std::sin( 3.0*i );
		curved.NodesUpdated();
	}
	mfem::Array<meq::AffineTriangle> curvedMaps;
	BOOST_TEST( !meq::affineTrianglesOf( curved, curvedMaps ) );

	mfem::L2_FECollection curvedCollection( 2, 2,
	                                        mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace curvedSpace( &curved, &curvedCollection );
	meq::DeviceTriangleBasis onCurved;
	BOOST_TEST( onCurved.build( curvedSpace ) );
}
