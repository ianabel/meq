#ifndef MEQ_TESTS_QDRIVEN_TARGET_HPP
#define MEQ_TESTS_QDRIVEN_TARGET_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "meq/Dataset.hpp"
#include "meq/Run.hpp"

/*
 * A q-DRIVEN TARGET MEASURED WHERE THE LOOP READS IT.
 *
 * The loop matches q at its knots, which follow the mesh
 * (meq::Run::safetyFactorKnots()), and reads a SafetyFactorFile there through
 * meq::SplineProfile's Hermite cubic. A table measured at other labels is met
 * only to its own interpolation error, which is then what a round trip reads.
 * So a round trip measures its reference at the knots themselves and at
 * `between - 1` labels between each pair -- a family uniform in rho over the
 * knots' own cut, which contains them -- with three-point slopes, and the
 * table is exact at every knot and dense between them.
 */
namespace meq::test
{
	/// The knots of the q-driven run @a text would make, after prepare();
	/// its relative paths read against @a configPath, or @a directory.
	inline std::vector<double> knotsOf( std::string const &directory, std::string const &text,
	                                    std::string const &configPath = "" )
	{
		meq::RunOptions options;
		options.configPath = configPath.empty() ? directory + "/knots.toml" : configPath;
		options.configText = text;
		options.baseDirectory = directory;
		options.version = "QDrivenTarget";
		meq::Run run( options );
		BOOST_TEST_REQUIRE( run.prepare() == meq::Solved );
		return run.safetyFactorKnots();
	}

	/// [output] keys for a family over @a knots' cut with @a between labels
	/// to each knot interval: it contains every knot.
	inline std::string surfacesAt( std::vector<double> const &knots, int between )
	{
		std::ostringstream out;
		out.precision( 17 );
		out << "FluxSurfaces = true\nFluxSurfaceCount = "
		    << between*( static_cast<int>( knots.size() ) - 1 ) + 1
		    << "\nFluxInnerCut = " << knots.front() << "\nFluxOuterCut = " << knots.back() << "\n";
		return out.str();
	}

	inline std::vector<double> datasetVariable( meq::Dataset const &data, std::string const &name )
	{
		meq::DatasetVariable const *v = data.variable( name );
		BOOST_TEST_REQUIRE( v != nullptr, "no variable " << name );
		return std::get<std::vector<double>>( v->values );
	}

	/**
	 * q = V' g < R^-2 >/4 pi^2 of forward runs' surfaces files, their labels
	 * merged, g from @a gOfPsi in the SOURCE's Psi, written to @a path as a
	 * SafetyFactorFile: Psi ascending (one on the axis), with the slope of the
	 * parabola through each point and its neighbours.
	 */
	template<typename G>
	void writeTarget( std::vector<meq::Dataset const *> const &families, G const &gOfPsi,
	                  std::string const &path, std::string const &header = "" )
	{
		std::vector<std::pair<double, double>> rows;   // ( Psi, q )
		for ( meq::Dataset const *surfaces : families )
		{
			std::vector<double> const label = datasetVariable( *surfaces, "normalised_flux" );
			std::vector<double> const vPrime = datasetVariable( *surfaces, "V_prime" );
			std::vector<double> const inverse = datasetVariable( *surfaces, "inverse_R_squared" );
			for ( std::size_t i = 0; i < label.size(); ++i )
			{
				double const sourcePsi = 1.0 - label[ i ];
				rows.emplace_back( sourcePsi, vPrime[ i ]*gOfPsi( sourcePsi )*inverse[ i ]/( 4.0*M_PI*M_PI ) );
			}
		}
		std::sort( rows.begin(), rows.end() );
		std::vector<double> psi, q;
		for ( auto const &[ x, y ] : rows )
			if ( psi.empty() || x - psi.back() > 1e-12 )
			{
				psi.push_back( x );
				q.push_back( y );
			}
		std::size_t const n = psi.size();
		std::ofstream out( path );
		out << header;
		out.precision( 17 );
		for ( std::size_t i = 0; i < n; ++i )
		{
			std::size_t const a = i == 0 ? 0 : ( i + 1 == n ? n - 3 : i - 1 );
			double const x0 = psi[ a ], x1 = psi[ a + 1 ], x2 = psi[ a + 2 ];
			double const x = psi[ i ];
			double const slope = q[ a ]*( 2.0*x - x1 - x2 )/( ( x0 - x1 )*( x0 - x2 ) )
			                     + q[ a + 1 ]*( 2.0*x - x0 - x2 )/( ( x1 - x0 )*( x1 - x2 ) )
			                     + q[ a + 2 ]*( 2.0*x - x0 - x1 )/( ( x2 - x0 )*( x2 - x1 ) );
			out << psi[ i ] << " " << q[ i ] << " " << slope << "\n";
		}
	}

	template<typename G>
	void writeTarget( meq::Dataset const &surfaces, G const &gOfPsi, std::string const &path,
	                  std::string const &header = "" )
	{
		writeTarget( std::vector<meq::Dataset const *>{ &surfaces }, gOfPsi, path, header );
	}

	/// A forward run of @a text and its surfaces file.
	inline std::shared_ptr<meq::Dataset const> forward( std::string const &directory,
	                                                     std::string const &prefix,
	                                                     std::string const &text )
	{
		meq::RunOptions options;
		options.configPath = directory + "/" + prefix + ".toml";
		options.configText = text;
		options.baseDirectory = directory;
		options.version = "QDrivenTarget";
		meq::Run run( options );
		BOOST_TEST_REQUIRE( run.prepare() == meq::Solved );
		BOOST_TEST_REQUIRE( run.solve() == meq::Solved, prefix << " did not solve" );
		std::shared_ptr<meq::Dataset const> data = run.fluxSurfaces( directory + "/" + prefix + "_surfaces.nc" );
		BOOST_TEST_REQUIRE( data != nullptr );
		return data;
	}
}

#endif
