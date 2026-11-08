#define BOOST_TEST_MODULE QDrivenConsistency
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "meq/Dataset.hpp"
#include "meq/Run.hpp"

#include "SourceTree.hpp"

/*
 * THE q-DRIVEN SOLVE IS CONSISTENT: THE q IT DELIVERS CONVERGES TO THE q IT WAS
 * ASKED FOR AS THE MESH REFINES. RED, DELIBERATELY -- TODO, "The q-driven solve
 * is inconsistent: g^2 is a fixed-degree global fit, so the delivered q misses
 * the target by the fit's floor however fine the mesh"; MEASUREMENTS.md M-361.
 *
 * The loop carries g^2 as a polynomial of fixed degree d, fitted over the
 * surfaces, so a target whose g^2 is not such a polynomial is met only to the
 * fit's floor -- and that floor is a property of d, not of h. Every other q
 * test in the tree builds its target from a g^2 the fit represents exactly,
 * where the floor is zero and the defect cannot show.
 *
 * HERE THE TARGET IS NOT REPRESENTABLE. It is the q of a forward solve with
 *
 *     g^2( Psi ) = g0^2 + B ( e^( kappa Psi ) - 1 )/( e^kappa - 1 )
 *
 * (Psi the source's, zero on the boundary, so g = g0 there; B = 1, kappa = 2),
 * taken on a mesh resolved well past the rungs: RefinementLevels = 3 at
 * k = 3, whose q differs from RefinementLevels = 2's by 1.0e-7 relative, six
 * orders under the miss read here (M-361). The q-driven run is then posed on
 * that target, pinned by its own boundary g, on two rungs in h at k = 2, and
 * the delivered q -- the surfaces file's own column -- is compared with the
 * reference's at the same labels.
 *
 * WHAT IT ASSERTS: the max-norm miss falls by at least 4 per halving of h, a
 * rate of 2, under the k + 1 = 3 the field and the extraction reach. What it
 * reads now is the fit's floor: 1.18e-1 at degree 2 on levels 0, 1 and 2, and
 * 1.2e-2 at degree 3 on the same three -- it moves with the degree and not
 * with h (M-361). It turns green when g^2 is carried in a representation that
 * refines with the discretisation, matched to q at every knot, and the pin
 * then closes the system with no q condition giving way.
 *
 * Coarse on purpose (the q-driven box, levels 0 and 1, k = 2): MEQ owes about
 * 1e-6, and the floor this case is about sits orders above that.
 */
namespace
{
	double const g0 = 2.2;
	double const amplitude = 1.0;
	double const kappa = 2.0;

	double gSquared( double psi )
	{
		return g0*g0 + amplitude*( std::exp( kappa*psi ) - 1.0 )/( std::exp( kappa ) - 1.0 );
	}

	double ggPrime( double psi )
	{
		return 0.5*amplitude*kappa*std::exp( kappa*psi )/( std::exp( kappa ) - 1.0 );
	}

	std::string configuration( std::string const &prefix, std::string const &toroidal,
	                           int levels, int order )
	{
		return
			"[mesh]\nRMin = 0.6\nRMax = 1.4\nZMin = -0.6\nZMax = 0.6\nNR = 8\nNZ = 8\n"
			"RefinementLevels = " + std::to_string( levels ) + "\n\n"
			"[discretisation]\nPolynomialDegree = " + std::to_string( order ) + "\nTau = 1.0\n\n"
			"[source]\nType = \"mhd\"\nNormalised = true\nPsiAxis = 0.3\nMu0 = 1.0\n"
			"PPrimeFile = \"" + meq::test::sourcePath( "examples/q-driven-pprime.dat" ) + "\"\n"
			+ toroidal + "\n"
			"[boundary]\nType = \"zero\"\n\n"
			"[initialguess]\nType = \"bump\"\nAmplitude = 0.3\nCentreR = 1.0\nCentreZ = 0.0\n"
			"RadiusR = 0.4\nRadiusZ = 0.6\n\n"
			"[solver]\nNewtonMaxIterations = 60\nNewtonRelativeTolerance = 1.0e-12\n"
			"NewtonAbsoluteTolerance = 1.0e-14\n\n"
			"[output]\nDirectory = \".\"\nPrefix = \"" + prefix + "\"\nGridNR = 17\nGridNZ = 17\n"
			"FluxSurfaces = true\nFluxSurfaceCount = 24\n";
	}

	std::vector<double> variable( meq::Dataset const &data, std::string const &name )
	{
		meq::DatasetVariable const *v = data.variable( name );
		BOOST_TEST_REQUIRE( v != nullptr, "no variable " << name );
		return std::get<std::vector<double>>( v->values );
	}

	struct Surfaces
	{
		std::vector<double> label, q;
	};

	/// Run @a text and return its surfaces file's labels and q -- the q the
	/// file reports, or one built from the closed-form g for a forward run.
	Surfaces run( std::string const &directory, std::string const &prefix,
	              std::string const &text, bool forward )
	{
		meq::RunOptions options;
		options.configPath = directory + "/" + prefix + ".toml";
		options.configText = text;
		options.baseDirectory = directory;
		options.version = "QDrivenConsistency";
		meq::Run run( options );
		BOOST_TEST_REQUIRE( run.prepare() == meq::Solved );
		int const code = run.solve();
		BOOST_TEST_REQUIRE( code == meq::Solved, prefix << " exits " << code );
		std::shared_ptr<meq::Dataset const> const data =
			run.fluxSurfaces( directory + "/" + prefix + "_surfaces.nc" );
		BOOST_TEST_REQUIRE( data != nullptr );
		Surfaces out;
		out.label = variable( *data, "normalised_flux" );
		if ( forward )
		{
			std::vector<double> const vPrime = variable( *data, "V_prime" );
			std::vector<double> const inverse = variable( *data, "inverse_R_squared" );
			for ( std::size_t i = 0; i < out.label.size(); ++i )
				out.q.push_back( vPrime[ i ]*std::sqrt( gSquared( 1.0 - out.label[ i ] ) )
				                 *inverse[ i ]/( 4.0*M_PI*M_PI ) );
		}
		else
			out.q = variable( *data, "safety_factor" );
		return out;
	}

	struct Fixture
	{
		std::string directory;
		Surfaces reference;

		Fixture()
			: directory( std::filesystem::absolute( "QDrivenConsistency-files" ).string() )
		{
			std::filesystem::create_directories( directory );
			{
				// gg' as a dense Hermite table, exact slopes: its own error is
				// O( 64^-4 ), far under anything read here.
				std::ofstream table( directory + "/consistency-ggprime.dat" );
				table.precision( 17 );
				int const knots = 64;
				for ( int i = 0; i <= knots; ++i )
				{
					double const psi = static_cast<double>( i )/knots;
					table << psi << " " << ggPrime( psi ) << " " << kappa*ggPrime( psi ) << "\n";
				}
			}
			reference = run( directory, "consistency-reference",
				configuration( "consistency-reference", "GGPrimeFile = \"consistency-ggprime.dat\"\n",
				               3, 3 ), true );

			// The target in the source's Psi, one on the axis, ascending, with
			// centred slopes: q-driven-q.dat's construction.
			std::ofstream target( directory + "/consistency-q.dat" );
			target.precision( 17 );
			std::size_t const n = reference.label.size();
			std::vector<double> psi, q;
			for ( std::size_t i = n; i-- > 0; )
			{
				psi.push_back( 1.0 - reference.label[ i ] );
				q.push_back( reference.q[ i ] );
			}
			for ( std::size_t i = 0; i < n; ++i )
			{
				std::size_t const a = i == 0 ? 0 : i - 1;
				std::size_t const b = i + 1 == n ? n - 1 : i + 1;
				target << psi[ i ] << " " << q[ i ] << " "
				       << ( q[ b ] - q[ a ] )/( psi[ b ] - psi[ a ] ) << "\n";
			}
		}

		/// max | q_delivered - q_target |/| q_target | over the table's labels.
		double miss( int levels, int order, int degree )
		{
			std::string const prefix = "consistency-l" + std::to_string( levels ) + "-k"
			                           + std::to_string( order ) + "-d" + std::to_string( degree );
			Surfaces const delivered = run( directory, prefix, configuration( prefix,
				"SafetyFactorFile = \"consistency-q.dat\"\nSafetyFactorDegree = "
				+ std::to_string( degree ) + "\nR0B0 = " + std::to_string( g0 ) + "\n",
				levels, order ), false );
			BOOST_TEST_REQUIRE( delivered.label.size() == reference.label.size() );
			double worst = 0.0;
			for ( std::size_t i = 0; i < delivered.q.size(); ++i )
			{
				BOOST_TEST_REQUIRE( std::abs( delivered.label[ i ] - reference.label[ i ] ) < 1.0e-12 );
				worst = std::max( worst, std::abs( delivered.q[ i ] - reference.q[ i ] )
				                         /std::abs( reference.q[ i ] ) );
			}
			std::printf( "    level %d, k = %d, degree %d: max | q - q_target |/| q_target | %.3e\n",
			             levels, order, degree, worst );
			return worst;
		}
	};
}

BOOST_AUTO_TEST_CASE( theDeliveredQConvergesToTheTarget )
{
	Fixture fixture;
	std::printf( "\n  THE DELIVERED q AGAINST A TARGET NO POLYNOMIAL g^2 REPRESENTS (M-361)\n" );
	double const coarse = fixture.miss( 0, 2, 2 );
	double const fine = fixture.miss( 1, 2, 2 );
	std::printf( "    factor per halving of h: %.2f\n", coarse/fine );
	BOOST_TEST( fine < coarse/4.0,
		"the q the q-driven run delivers misses its target by " << coarse << " at h and "
		<< fine << " at h/2: it does not fall with the mesh. g^2 is carried as a "
		"polynomial of fixed degree fitted over the surfaces, so a target whose g^2 is "
		"not one is met only to the fit's floor, whatever the mesh. TODO, \"The q-driven "
		"solve is inconsistent: g^2 is a fixed-degree global fit\"; MEASUREMENTS.md M-361" );
}

/*
 * THE STUDY BEHIND M-361, run by name and not by the suite: the miss by degree
 * and by mesh, and the reference against a coarser one.
 */
BOOST_AUTO_TEST_CASE( theMissMovesWithTheDegreeAndNotTheMesh, * boost::unit_test::disabled() )
{
	Fixture fixture;
	std::printf( "\n  M-361: the miss by degree and mesh\n" );
	for ( int degree : { 2, 3 } )
		for ( int levels : { 0, 1, 2 } )
			fixture.miss( levels, 2, degree );

	Surfaces const coarser = run( fixture.directory, "consistency-reference-l2",
		configuration( "consistency-reference-l2", "GGPrimeFile = \"consistency-ggprime.dat\"\n",
		               2, 3 ), true );
	double worst = 0.0;
	for ( std::size_t i = 0; i < coarser.q.size(); ++i )
		worst = std::max( worst, std::abs( coarser.q[ i ] - fixture.reference.q[ i ] )
		                         /std::abs( fixture.reference.q[ i ] ) );
	std::printf( "    reference, level 2 against level 3 at k = 3: %.3e\n", worst );
}
