#define BOOST_TEST_MODULE QDrivenConsistency
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "meq/Dataset.hpp"
#include "meq/Run.hpp"
#include "meq/SafetyFactor.hpp"

#include "SourceTree.hpp"
#include "convergence/QDrivenTarget.hpp"

/*
 * THE q-DRIVEN SOLVE IS CONSISTENT: WHAT IT DELIVERS CONVERGES TO WHAT IT WAS
 * ASKED FOR AS THE MESH, THE DEGREE AND THE KNOTS REFINE (MEASUREMENTS.md
 * M-361, M-QR-3).
 *
 * THE TARGET IS NOT REPRESENTABLE at any finite knot count. It is the q of a
 * forward solve with
 *
 *     g^2( Psi ) = g0^2 + B ( e^( kappa Psi ) - 1 )/( e^kappa - 1 )
 *
 * (Psi the source's, zero on the boundary, so g = g0 there; B = 1, kappa = 2)
 * on a reference resolved past the rungs, RefinementLevels = 2 at k = 3, read
 * at each rung's own knots and three labels between each pair
 * (QDrivenTarget.hpp) -- so the target the loop reads at its knots is the
 * reference's own, and no table interpolates it there. Each rung is then posed
 * on its target, pinned by the closed form's own boundary g.
 *
 * WHAT IT MEASURES, AND WHY NOT q AT THE KNOTS. q is matched at the knots in
 * least squares, so the miss there is the discrete problem's inconsistency
 * and says little about the answer. What the run delivers is g^2, which is
 * known in closed form: the measure is max | g^2/g^2_exact - 1 | on a dense
 * grid over the knots' cut -- between the knots as well as at them -- the
 * interior, where q is finite on this box (its corners make q log-singular at
 * the edge). The part of the request outside the cut is one knot spacing
 * wide and shrinks with the mesh (meq::toroidalKnotCut()).
 *
 * WHAT IT ASSERTS: along the diagonal -- h halved, the knots following it --
 * the g^2 error falls by at least 6 a step, against the 2^3 = 8 the averages'
 * k + 1 gives at k = 2 (M-QR-3 measured 9.6 and 8.3 in the study; this case
 * reads 18 at k = 2 and 11 at k = 3, M-QR-7). M-361's degree-2 and degree-3
 * fits read 1.18e-1 and 1.2e-2 in q here, flat in h.
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
	                           int levels, int order, std::string const &surfaces = "" )
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
			+ surfaces;
	}

	struct Delivered
	{
		double gSquaredError = 0.0;
		double misfit = 0.0;
		std::size_t knots = 0;
		double seconds = 0.0;
	};

	struct Fixture
	{
		std::string directory;

		Fixture()
			: directory( std::filesystem::absolute( "QDrivenConsistency-files" ).string() )
		{
			std::filesystem::create_directories( directory );
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

		Delivered deliver( int levels, int order )
		{
			std::string const tag = "l" + std::to_string( levels ) + "-k" + std::to_string( order );
			std::string const target = "consistency-q-" + tag + ".dat";
			std::string const driven = "SafetyFactorFile = \"" + target + "\"\nR0B0 = " + std::to_string( g0 ) + "\n";
			{
				std::ofstream provisional( directory + "/" + target );
				provisional << "0.0 1.0 0.0\n1.0 1.0 0.0\n";
			}
			std::vector<double> const knots = meq::test::knotsOf(
				directory, configuration( "consistency-" + tag, driven, levels, order ) );
			std::shared_ptr<meq::Dataset const> const surfaces = meq::test::forward(
				directory, "consistency-reference-" + tag,
				configuration( "consistency-reference-" + tag,
				               "GGPrimeFile = \"consistency-ggprime.dat\"\n", 2, 3,
				               meq::test::surfacesAt( knots, 4 ) ) );
			meq::test::writeTarget( *surfaces, []( double psi ) { return std::sqrt( gSquared( psi ) ); },
			                        directory + "/" + target );

			meq::RunOptions options;
			options.configPath = directory + "/consistency-" + tag + ".toml";
			options.configText = configuration( "consistency-" + tag, driven, levels, order );
			options.baseDirectory = directory;
			options.version = "QDrivenConsistency";
			meq::Run run( options );
			BOOST_TEST_REQUIRE( run.prepare() == meq::Solved );
			auto const started = std::chrono::steady_clock::now();
			int const code = run.solve();
			Delivered out;
			out.seconds = std::chrono::duration<double>( std::chrono::steady_clock::now() - started ).count();
			BOOST_TEST_REQUIRE( code == meq::Solved, tag << " exits " << code );
			BOOST_TEST_REQUIRE( run.writeRestart( directory + "/consistency-" + tag ) == meq::Solved );
			meq::Dataset const &data = *run.dataset();
			auto numbers = [ & ]( char const *name )
			{
				meq::DatasetAttribute const *a = data.attribute( name );
				BOOST_TEST_REQUIRE( a != nullptr, "no attribute " << name );
				std::istringstream in( std::get<std::string>( a->value ) );
				std::vector<double> values;
				for ( double value = 0.0; in >> value; )
					values.push_back( value );
				return values;
			};
			meq::ToroidalFieldSpline spline( numbers( "g_squared_knots" ) );
			spline.setValues( numbers( "g_squared" ) );
			BOOST_TEST_REQUIRE( spline.size() == knots.size() );
			for ( int k = 0; k <= 2000; ++k )
			{
				double const psiN = knots.front() + ( knots.back() - knots.front() )*k/2000.0;
				out.gSquaredError = std::max( out.gSquaredError,
					std::abs( spline( 1.0 - psiN )/gSquared( 1.0 - psiN ) - 1.0 ) );
			}
			meq::DatasetAttribute const *misfit = data.attribute( "toroidal_loop_misfit" );
			BOOST_TEST_REQUIRE( misfit != nullptr );
			out.misfit = std::get<std::vector<double>>( misfit->value ).front();
			out.knots = knots.size();
			std::printf( "    level %d, k = %d, %2zu knots: max | g^2/g^2_exact - 1 | %.3e, misfit at the "
			             "knots %.3e, %.1f s\n", levels, order, out.knots, out.gSquaredError, out.misfit,
			             out.seconds );
			return out;
		}
	};
}

BOOST_AUTO_TEST_CASE( theDeliveredFieldConvergesToTheTarget )
{
	Fixture fixture;
	std::printf( "\n  THE DELIVERED g^2 AGAINST A TARGET NO FINITE KNOT SET REPRESENTS (M-QR-3)\n" );
	Delivered const coarse = fixture.deliver( 0, 2 );
	Delivered const fine = fixture.deliver( 1, 2 );
	double const factor = coarse.gSquaredError/fine.gSquaredError;
	std::printf( "    factor per halving of h, the knots following: %.2f\n", factor );
	BOOST_TEST( factor > 6.0,
		"the g^2 the q-driven run delivers misses the closed form its target was measured "
		"from by " << coarse.gSquaredError << " at h and " << fine.gSquaredError << " at h/2, "
		"a factor of " << factor << " against the 8 of k + 1 = 3: the solve does not converge "
		"to the problem posed. If it is flat, the representation does not refine with the mesh "
		"(meq::toroidalKnotCount())" );
	// THE MISFIT AT THE KNOTS IS PRINTED AND NOT ASSERTED: it is the least
	// squares' residual, and on these two rungs it reads 5.4e-6 then 1.4e-4 --
	// the innermost knot sits a knot spacing from the axis on a surface about
	// one element across at every level (TODO, *The q-driven misfit at the
	// innermost knot*). What the run delivers, g^2, is what converges.
}

BOOST_AUTO_TEST_CASE( theDeliveredFieldConvergesAtCubicDegree )
{
	Fixture fixture;
	std::printf( "\n  THE SAME AT k = 3 (M-QR-3)\n" );
	Delivered const coarse = fixture.deliver( 0, 3 );
	Delivered const fine = fixture.deliver( 1, 3 );
	double const factor = coarse.gSquaredError/fine.gSquaredError;
	std::printf( "    factor per halving of h, the knots following: %.2f\n", factor );
	BOOST_TEST( factor > 6.0,
		"at k = 3 the delivered g^2 misses by " << coarse.gSquaredError << " then "
		<< fine.gSquaredError << ", a factor of " << factor );
}
