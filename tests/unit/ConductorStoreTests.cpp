#define BOOST_TEST_MODULE ConductorStoreTests
#include <boost/test/included/unit_test.hpp>

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include <netcdf>

#include "meq/ConductorStore.hpp"

/*
 * THE FILE-BACKED CONDUCTOR CACHE.
 *
 * `src/meq/ConductorStore.hpp` says what this is for; what it has to be is
 * EXACT and REFUSABLE, and those are the two things asserted here.
 *
 * Exact, because the whole reason this is a cache of points rather than an
 * mgrid-style grid of samples is that `ConductorSubtraction`'s three identities
 * read `0.000e+00` and must keep doing so. A round trip that agreed to 1e-15
 * would be a different product; so the comparisons below are `==` on doubles,
 * deliberately, and that is not a tolerance anybody forgot to add.
 *
 * Refusable, because the failure this file can cause has no symptom. A `psi_c`
 * belonging to a machine whose coils moved is a plausible field: the solve
 * converges, the rates hold, and the answer is to a problem nobody posed.
 */
namespace
{
	/// A scratch path that cleans itself up. The unit tests run with the
	/// SOURCE tree as their working directory -- tests/CMakeLists.txt pins it
	/// so ConfigTests can open examples/*.toml -- so a file written relatively
	/// would land in the checkout.
	class ScratchFile
	{
		public:
			explicit ScratchFile( char const *stem )
				: pathValue( std::filesystem::temp_directory_path()
				             / ( std::string( "meq-conductor-store-" ) + stem
				                 + ".nc" ) )
			{
				std::filesystem::remove( pathValue );
			}

			~ScratchFile()
			{
				std::error_code ignored;
				std::filesystem::remove( pathValue, ignored );
			}

			ScratchFile( ScratchFile const & ) = delete;
			ScratchFile &operator=( ScratchFile const & ) = delete;

			std::string path() const { return pathValue.string(); }

		private:
			std::filesystem::path pathValue;
	};

	/// Values chosen so that a transposition, a truncation or an off-by-one
	/// would show: no repeats, and magnitudes an eye can order.
	meq::ConductorCache sampleCache()
	{
		meq::ConductorCache cache;
		cache.signature.meshDigest = 0x0123456789abcdefULL;
		cache.signature.conductorDigest = 0xfedcba9876543210ULL;
		cache.signature.polynomialDegree = 3;
		cache.signature.elements = 4;
		cache.signature.potentialDofs = 6;

		cache.nodalPsi = { -1.5, 2.25, -3.125e-7, 4.0e12, 0.0, 6.5 };
		cache.quadraturePsi = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5 };
		cache.quadratureOffset = { 0, 2, 5, 6, 9 };

		// THE q_c HALF, over the same four elements. Two potential nodes and
		// three flux nodes per element here -- deliberately different counts,
		// because MEQ does not require the two spaces to share a node set and
		// a format that assumed they did would pass a test that used one.
		cache.criticalPotentialPsi = { 1.0, 2.0, 3.0, 4.0,
		                               5.0, 6.0, 7.0, 8.0 };
		cache.criticalPotentialOffset = { 0, 2, 4, 6, 8 };
		cache.criticalFluxOffset = { 0, 3, 6, 9, 12 };
		cache.criticalFluxQ.clear();
		for ( int i = 0; i < 12; ++i )
		{
			cache.criticalFluxQ.push_back( 0.125*( i + 1 ) );
			cache.criticalFluxQ.push_back( -0.25*( i + 1 ) );
		}
		// One screened node, which is what the axis gives.
		cache.criticalFluxUsable.assign( 12, 1 );
		cache.criticalFluxUsable[ 5 ] = 0;
		cache.criticalFluxScale = 3.0;
		return cache;
	}
}

BOOST_AUTO_TEST_CASE( the_round_trip_is_bit_for_bit )
{
	ScratchFile scratch( "roundtrip" );
	meq::ConductorCache const written = sampleCache();
	BOOST_TEST_REQUIRE( written.consistent() );

	meq::writeConductorCache( scratch.path(), written );

	std::optional< meq::ConductorCache > const read
		= meq::readConductorCache( scratch.path() );
	BOOST_TEST_REQUIRE( read.has_value(),
	                    "a file this code just wrote carries no cache" );

	BOOST_CHECK( read->signature == written.signature );
	BOOST_TEST_REQUIRE( read->nodalPsi.size() == written.nodalPsi.size() );
	BOOST_TEST_REQUIRE( read->quadraturePsi.size()
	                    == written.quadraturePsi.size() );

	// == ON DOUBLES, ON PURPOSE. See the header of this file: a cache that
	// round-tripped to within a tolerance would no longer be a substitute for
	// the recompute, and the identities that motivate the whole split would
	// stop reading exactly zero.
	for ( std::size_t i = 0; i < written.nodalPsi.size(); ++i )
		BOOST_TEST( read->nodalPsi[ i ] == written.nodalPsi[ i ],
		            "nodal value " << i << " came back as "
		            << read->nodalPsi[ i ] << " rather than "
		            << written.nodalPsi[ i ]
		            << ". This must be EXACT, not close: psi_c is what makes "
		               "the conductors resolved rather than discretised." );

	for ( std::size_t i = 0; i < written.quadraturePsi.size(); ++i )
		BOOST_TEST( read->quadraturePsi[ i ] == written.quadraturePsi[ i ] );

	BOOST_TEST( read->quadratureOffset == written.quadratureOffset,
	            boost::test_tools::per_element() );

	// THE q_c HALF, BY THE SAME EXACT COMPARISON. It is the expensive one --
	// `grad psi_c` at every flux node, a cross-section quadrature of elliptic
	// integrals per node for a rectangle -- and a gradient that came back
	// approximately would move the critical points, which is where psi_ax and
	// the plasma edge come from.
	BOOST_TEST_REQUIRE( read->hasCriticalPointTable() );
	BOOST_TEST_REQUIRE( read->criticalFluxQ.size()
	                    == written.criticalFluxQ.size() );
	for ( std::size_t i = 0; i < written.criticalFluxQ.size(); ++i )
		BOOST_TEST( read->criticalFluxQ[ i ] == written.criticalFluxQ[ i ],
		            "q_c component " << i << " came back as "
		            << read->criticalFluxQ[ i ] << " rather than "
		            << written.criticalFluxQ[ i ] );

	for ( std::size_t i = 0; i < written.criticalPotentialPsi.size(); ++i )
		BOOST_TEST( read->criticalPotentialPsi[ i ]
		            == written.criticalPotentialPsi[ i ] );

	BOOST_TEST( read->criticalFluxScale == written.criticalFluxScale );
	BOOST_TEST( read->criticalFluxOffset == written.criticalFluxOffset,
	            boost::test_tools::per_element() );
	BOOST_TEST( read->criticalPotentialOffset == written.criticalPotentialOffset,
	            boost::test_tools::per_element() );

	// THE SCREEN TRAVELS WITH THE VALUES rather than being re-derived. q_c is
	// NaN on the axis, and a reader that recomputed "usable" from the values
	// would be reading NaN to decide whether to read NaN.
	BOOST_TEST_REQUIRE( read->criticalFluxUsable.size()
	                    == written.criticalFluxUsable.size() );
	BOOST_TEST( read->criticalFluxUsable[ 5 ] == 0,
	            "the screened flux node came back usable" );
	BOOST_TEST( read->criticalFluxUsable == written.criticalFluxUsable,
	            boost::test_tools::per_element() );
}

BOOST_AUTO_TEST_CASE( the_q_c_half_is_present_as_a_whole_or_not_at_all )
{
	// Four arrays that index each other, so any one missing leaves the others
	// indexing something that is not there. The quadrature half follows the
	// same rule for the same reason; this one has more ways to be half built.
	meq::ConductorCache cache = sampleCache();
	BOOST_TEST_REQUIRE( cache.consistent() );

	cache = sampleCache();
	cache.criticalFluxUsable.pop_back();          // one screen short
	BOOST_TEST( !cache.consistent() );

	cache = sampleCache();
	cache.criticalFluxQ.pop_back();               // not two per node
	BOOST_TEST( !cache.consistent() );

	cache = sampleCache();
	cache.criticalFluxOffset = { 0, 3, 2, 9, 12 };  // not monotone
	BOOST_TEST( !cache.consistent() );

	cache = sampleCache();
	cache.criticalPotentialOffset.pop_back();     // not elements + 1
	BOOST_TEST( !cache.consistent() );

	// AND ABSENT IS FINE: a fixed-boundary run locates no critical point and
	// builds no table, which is not a truncated file.
	cache = sampleCache();
	cache.criticalPotentialPsi.clear();
	cache.criticalPotentialOffset.clear();
	cache.criticalFluxQ.clear();
	cache.criticalFluxUsable.clear();
	cache.criticalFluxOffset.clear();
	BOOST_TEST( cache.consistent() );
	BOOST_TEST( !cache.hasCriticalPointTable() );
}

BOOST_AUTO_TEST_CASE( it_embeds_in_an_existing_file_as_a_group )
{
	/*
	 * THE TWO DEPLOYMENTS ARE ONE FORMAT, which is the point of the group: a
	 * cache of its own for something that outlives a run, and a group inside
	 * the equilibrium `.nc` so one file carries the answer and the means to
	 * warm start from it. Here the "equilibrium" is a file with an unrelated
	 * variable in it, which is the part that has to survive.
	 */
	ScratchFile scratch( "embedded" );
	{
		netCDF::NcFile file( scratch.path(), netCDF::NcFile::replace );
		netCDF::NcDim dim = file.addDim( "R", 3 );
		netCDF::NcVar var = file.addVar( "psi", netCDF::ncDouble, dim );
		std::vector< double > const values{ 1.0, 2.0, 3.0 };
		var.putVar( values.data() );
	}

	meq::ConductorCache const written = sampleCache();
	meq::writeConductorCache( scratch.path(), written,
	                          meq::conductorCacheGroupName() );

	// Found without being told the group name, which is what lets one reader
	// take either deployment.
	std::optional< meq::ConductorCache > const read
		= meq::readConductorCache( scratch.path() );
	BOOST_TEST_REQUIRE( read.has_value(),
	                    "the embedded group was not found by a bare read" );
	BOOST_CHECK( read->signature == written.signature );
	BOOST_TEST( read->nodalPsi.size() == written.nodalPsi.size() );

	// And the equilibrium it was added to is untouched.
	netCDF::NcFile file( scratch.path(), netCDF::NcFile::read );
	netCDF::NcVar psi = file.getVar( "psi" );
	BOOST_TEST_REQUIRE( !psi.isNull(),
	                    "embedding the cache destroyed the file it was added to" );
	std::vector< double > back( 3, 0.0 );
	psi.getVar( back.data() );
	BOOST_TEST( back[ 2 ] == 3.0 );
}

BOOST_AUTO_TEST_CASE( a_file_with_no_cache_is_an_absence_and_not_an_error )
{
	// The ordinary case for an equilibrium written before this existed, or by
	// a run with no conductors. The caller's next move is to compute one, so
	// this must not throw.
	ScratchFile scratch( "nocache" );
	{
		netCDF::NcFile file( scratch.path(), netCDF::NcFile::replace );
		file.addDim( "R", 2 );
	}

	std::optional< meq::ConductorCache > const read
		= meq::readConductorCache( scratch.path() );
	BOOST_TEST( !read.has_value(),
	            "a file with no cache in it reported one" );
}

BOOST_AUTO_TEST_CASE( a_malformed_cache_throws_rather_than_reading_as_absent )
{
	/*
	 * THE DISTINCTION THIS ASSERTS IS THE WHOLE VALUE OF THE CHECK. A file
	 * that says it is one of these and is not is a different thing from a file
	 * that never claimed to be: reading the first as an absence sends the
	 * caller down the "compute one" path and hides a truncated or
	 * foreign-written file behind a merely slow run.
	 */
	ScratchFile scratch( "malformed" );
	{
		netCDF::NcFile file( scratch.path(), netCDF::NcFile::replace );
		file.putAtt( "format", "meq-conductor-cache" );
		file.putAtt( "format_version", netCDF::ncInt, 1 );
		file.putAtt( "mesh_digest", std::string( "0123456789abcdef" ) );
		file.putAtt( "conductor_digest", std::string( "fedcba9876543210" ) );
		file.putAtt( "polynomial_degree", netCDF::ncInt, 3 );
		file.putAtt( "elements", netCDF::ncInt, 4 );
		file.putAtt( "potential_dofs", netCDF::ncInt, 6 );

		// Six dofs claimed, three values written.
		netCDF::NcDim dim = file.addDim( "potential_dof", 3 );
		netCDF::NcVar var = file.addVar( "psi_c_dof", netCDF::ncDouble, dim );
		std::vector< double > const values{ 1.0, 2.0, 3.0 };
		var.putVar( values.data() );
	}

	BOOST_CHECK_THROW( meq::readConductorCache( scratch.path() ),
	                   std::runtime_error );
}

BOOST_AUTO_TEST_CASE( a_future_format_version_is_refused_by_name )
{
	ScratchFile scratch( "version" );
	{
		netCDF::NcFile file( scratch.path(), netCDF::NcFile::replace );
		file.putAtt( "format", "meq-conductor-cache" );
		file.putAtt( "format_version", netCDF::ncInt, 99 );
	}

	BOOST_CHECK_THROW( meq::readConductorCache( scratch.path() ),
	                   std::runtime_error );
}

BOOST_AUTO_TEST_CASE( the_signature_names_which_field_moved )
{
	/*
	 * A refusal that says only "the cache does not match" costs the reader the
	 * five-way guess this exists to answer, and the five are not equally
	 * likely: a coupled run changes the profiles and nothing else, so a
	 * conductor digest that moved means somebody edited the machine while a
	 * mesh digest that moved means the run re-meshed, which is ordinary.
	 */
	meq::ConductorCacheSignature const base = sampleCache().signature;

	BOOST_TEST( base.difference( base ).empty(),
	            "a signature differs from itself" );

	meq::ConductorCacheSignature moved = base;
	moved.conductorDigest ^= 1ULL;
	BOOST_TEST( base.difference( moved ).find( "conductors" )
	            != std::string::npos,
	            "a changed conductor digest is reported as: "
	            << base.difference( moved ) );

	moved = base;
	moved.meshDigest ^= 1ULL;
	BOOST_TEST( base.difference( moved ).find( "mesh" ) != std::string::npos );

	moved = base;
	moved.polynomialDegree += 1;
	BOOST_TEST( base.difference( moved ).find( "degree" )
	            != std::string::npos );

	moved = base;
	moved.potentialDofs += 1;
	BOOST_TEST( base.difference( moved ).find( "dof" ) != std::string::npos );

	// AND THE CONDUCTOR DIGEST IS REPORTED FIRST when several moved at once,
	// because it is the one that means a mistake rather than a re-mesh.
	moved = base;
	moved.meshDigest ^= 1ULL;
	moved.conductorDigest ^= 1ULL;
	BOOST_TEST( base.difference( moved ).find( "conductors" )
	            != std::string::npos );
}

BOOST_AUTO_TEST_CASE( an_inconsistent_cache_is_refused_at_the_writer )
{
	/*
	 * Refusing on the way OUT rather than only on the way in: a cache written
	 * in a state the reader will reject is a silent waste -- every later run
	 * recomputes and nothing says why -- where a throw names the bug at the
	 * one place that can still fix it.
	 */
	ScratchFile scratch( "inconsistent" );

	meq::ConductorCache cache = sampleCache();
	cache.quadratureOffset.back() = 8;      // the flattened size is 9
	BOOST_TEST_REQUIRE( !cache.consistent() );
	BOOST_CHECK_THROW( meq::writeConductorCache( scratch.path(), cache ),
	                   std::invalid_argument );

	cache = sampleCache();
	cache.quadratureOffset = { 0, 5, 2, 6, 9 };   // not monotone
	BOOST_TEST( !cache.consistent(),
	            "offsets that go backwards passed the consistency check" );

	cache = sampleCache();
	cache.quadratureOffset.clear();               // values without offsets
	BOOST_TEST( !cache.consistent() );

	cache = sampleCache();
	cache.nodalPsi.pop_back();                    // fewer values than dofs
	BOOST_TEST( !cache.consistent() );
}

BOOST_AUTO_TEST_CASE( the_nodal_half_may_travel_alone )
{
	// A warm start's stored guess needs psi_c at the dofs and nothing else --
	// apps/meq.cpp subtracts it node by node when [initialguess] Content is
	// "total" -- so a cache with no quadrature half is legitimate rather than
	// truncated, and must not be confused with the malformed case above.
	ScratchFile scratch( "nodalonly" );

	meq::ConductorCache cache = sampleCache();
	cache.quadraturePsi.clear();
	cache.quadratureOffset.clear();
	BOOST_TEST_REQUIRE( cache.consistent() );

	meq::writeConductorCache( scratch.path(), cache );
	std::optional< meq::ConductorCache > const read
		= meq::readConductorCache( scratch.path() );

	BOOST_TEST_REQUIRE( read.has_value() );
	BOOST_TEST( read->quadraturePsi.empty() );
	BOOST_TEST( read->nodalPsi.size() == cache.nodalPsi.size() );
}

BOOST_AUTO_TEST_CASE( the_digest_separates_what_psi_c_depends_on )
{
	/*
	 * The digest helpers are exposed so that two owners -- the solver for the
	 * mesh, meq::ConductorField for the conductors -- can mix in the same
	 * order. What they have to give is different numbers for different inputs;
	 * the only interesting case is the one a weaker mixer gets wrong.
	 */
	std::uint64_t const seed = meq::digestSeed();

	BOOST_TEST( meq::digestAppend( seed, 1.0 )
	            != meq::digestAppend( seed, 2.0 ) );

	// ORDER MATTERS, which is why both sides must mix the same way: a
	// commutative combiner would give one digest for two machines whose coils
	// are at each other's radii.
	BOOST_TEST( meq::digestAppend( meq::digestAppend( seed, 1.0 ), 2.0 )
	            != meq::digestAppend( meq::digestAppend( seed, 2.0 ), 1.0 ) );

	// A double and the integer of the same value are different inputs, so a
	// current of 3 A and three conductors must not collide.
	BOOST_TEST( meq::digestAppend( seed, 3.0 )
	            != meq::digestAppend( seed, static_cast< std::int64_t >( 3 ) ) );
}
