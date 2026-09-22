/*
 * MEQ'S ROTATING SOURCE AGAINST AN INDEPENDENT IMPLEMENTATION OF THE SAME
 * EQUATION.
 *
 * `../geq` drives `../freegs4e`'s `ProfilesCentrifugalMirror` for rotating
 * magnetic mirrors, and it implements Abel et al. 2013, Rep. Prog. Phys. 76
 * 116201, eq (136), closed by its (96) and (97) -- the same paper from the same
 * equations that `meq::RotatingSource` is built on. Two independent
 * implementations of one equation is a rarer thing than it sounds, and it is
 * what this file spends.
 *
 * NO SOLVER RUNS HERE AND NO MESH EXISTS. `tools/geq-benchmark/export_mirror.py`
 * writes a PRESCRIBED state -- a vacuum-field mirror with species, rotation and
 * the quasineutrality potential solved -- and this evaluates MEQ's own F at
 * geq's own converged state, against geq's own current density. That separates
 * the SOURCE from the discretisation and from the iteration exactly, the same
 * way `tools/freegs4e-benchmark/source_check.py` does for the MHD source. A red
 * here is a fault in one of the two sources; a green moves the search elsewhere.
 *
 * THE TWO FORMS AGREE ANALYTICALLY, WHICH IS WHY A DISAGREEMENT WOULD MEAN
 * SOMETHING. geq writes
 *
 *     mu0 R Jtor = mu0 R^2 Sum_s n_s [ T_s dlnN_s/dpsi
 *                                      + ( chi_s + 1 ) T_s dlnT_s/dpsi ]
 *                + mu0 R^2 Sum_s m_s n_s R^2 omega domega/dpsi
 *
 * and MEQ writes `F = mu0 R^2 dp/dpsi|_r + g g'` with `p = Sum_s n_s T_s`.
 * Expanding MEQ's derivative reproduces geq's bracket term for term once the
 * `Sum_s Z_s n_s d( e phi_0 )/dpsi` piece is dropped, and it drops by
 * quasineutrality -- the cancellation `docs/rotation.rst` records and the reason
 * the residual needs `phi_0` but never its derivative. So the two are the same
 * expression written by two people, and this file is the check that they were
 * both written correctly.
 *
 * WHAT THIS REACHES THAT NOTHING ELSE IN THE TREE DOES.
 * `tests/analytic/VaryingCentrifugal.hpp` closed the `C'( psi )` gap by a
 * MANUFACTURED fixture, because `C` constant is exactly what makes the equation
 * solvable in closed form and no published benchmark can vary it. That fixture
 * and MEQ share an author and a reading of the paper. This one does not: the
 * exported state drifts `C` by about 7x with both `T( psi )` and `omega( psi )`
 * varying, and nothing in it came from MEQ.
 *
 * THE GAUGE IS THE CONVERSION AND THE TEST SAYS SO OUT LOUD. MEQ pins
 * `phi_0( R_ref, psi ) = 0`, one condition per flux surface; geq pins one global
 * point at `psi_n = 0.5` on the midplane. `N_s` absorbs the difference through
 * `exp( Z_s e delta/T_s )`, a DIFFERENT factor per species, so no single
 * rescaling relates the two codes' tabulated densities. The export solves geq's
 * own quasineutrality at `R = R_ref` to do the transfer, and it also writes
 * geq's RAW `N_s` beside it -- so `theGaugeTransferIsLoadBearing` can assert
 * that handing MEQ the raw tables is WRONG. Without that assertion a transfer
 * that happened to be near-identity would look like a passing comparison.
 *
 * THE REFERENCE IS COMMITTED, so this runs without geq, without freegs4e and
 * without a Python environment. Regenerating it needs all three; the command is
 * in `tools/geq-benchmark/README.md`.
 */

#define BOOST_TEST_MODULE GeqSourceComparison
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Profiles.hpp"
#include "RotatingSource.hpp"

namespace
{
	/// The committed reference. MEQ_GEQ_REFERENCE overrides it, which is how the
	/// knot study under tools/geq-benchmark/README.md sweeps the tabulation
	/// without touching what is committed.
	std::string referenceRootPath()
	{
		char const * const override = std::getenv( "MEQ_GEQ_REFERENCE" );
		return override ? std::string( override ) + "/" : std::string( "tools/geq-benchmark/reference/" );
	}

	/// One exported case: what meta.txt says, and the tables it names.
	struct Reference
	{
		double referenceRadius = 0.0;
		double mu0 = 0.0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;

		struct SpeciesEntry
		{
			std::string name;
			double mass = 0.0;
			double charge = 0.0;
			std::string temperatureFile;
			std::string densityFile;
		};

		std::vector<SpeciesEntry> species;

		struct Point
		{
			double radius = 0.0, z = 0.0, psi = 0.0, f = 0.0;
		};

		std::vector<Point> points;       ///< geq's own phi_0, as it ships
		std::vector<Point> pointsTight;  ///< phi_0 solved pointwise; see below
		std::string directory;
	};

	Reference readReference( std::string const & caseName )
	{
		Reference reference;
		reference.directory = referenceRootPath() + caseName + "/";

		std::string const metaPath = reference.directory + "meta.txt";
		std::ifstream meta( metaPath );
		if ( !meta )
			throw std::runtime_error(
				"cannot open " + metaPath + " -- the exported geq reference is missing. "
				"It is committed, so a build tree that lacks it is running from the "
				"wrong working directory; ctest pins CMAKE_SOURCE_DIR. To regenerate "
				"it see tools/geq-benchmark/README.md." );

		std::string line;
		while ( std::getline( meta, line ) )
		{
			std::istringstream in( line );
			std::string key;
			in >> key;

			if ( key == "R_ref" )
				in >> reference.referenceRadius;
			else if ( key == "mu0" )
				in >> reference.mu0;
			else if ( key == "psi_axis" )
				in >> reference.psiAxis;
			else if ( key == "psi_boundary" )
				in >> reference.psiBoundary;
			else if ( key == "name" )
			{
				Reference::SpeciesEntry entry;
				std::string label;
				in >> entry.name
				   >> label >> entry.mass
				   >> label >> entry.charge
				   >> label >> entry.temperatureFile
				   >> label >> entry.densityFile;
				reference.species.push_back( entry );
			}
		}

		auto const readPoints = [ &reference ]( std::string const & name,
			std::vector<Reference::Point> & into )
		{
			std::string const path = reference.directory + name;
			std::ifstream points( path );
			if ( !points )
				throw std::runtime_error( "cannot open " + path );

			std::string line;
			while ( std::getline( points, line ) )
			{
				std::size_t const first = line.find_first_not_of( " \t" );
				if ( first == std::string::npos || line[ first ] == '#' )
					continue;

				Reference::Point point;
				std::istringstream in( line );
				if ( in >> point.radius >> point.z >> point.psi >> point.f )
					into.push_back( point );
			}
		};

		readPoints( "points.dat", reference.points );
		readPoints( "points_tight.dat", reference.pointsTight );

		return reference;
	}

	/// The species set in MEQ's gauge, or -- with `raw` -- in geq's, which is the
	/// mistake the export writes a second set of tables to make checkable.
	std::vector<meq::Species> speciesFrom( Reference const & reference, bool raw )
	{
		std::vector<meq::Species> species;

		for ( auto const & entry : reference.species )
		{
			std::string const density = raw
				? "Nraw_" + entry.name + ".dat"
				: entry.densityFile;

			meq::Species one;
			one.mass = entry.mass;
			one.charge = entry.charge;
			one.temperature = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( reference.directory + entry.temperatureFile ) );
			one.density = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( reference.directory + density ) );
			species.push_back( one );
		}

		return species;
	}

	std::shared_ptr<meq::Profile const> omegaFrom( Reference const & reference )
	{
		return std::make_shared<meq::SplineProfile const>(
			meq::SplineProfile::fromFile( reference.directory + "omega.dat" ) );
	}

	/// Relative L2 of MEQ's F against geq's, over the exported points, and the
	/// worst single point with it.
	struct Agreement
	{
		double relativeL2 = 0.0;
		double worstRelative = 0.0;
		double referenceNorm = 0.0;
		double referenceMax = 0.0;
	};

	Agreement compare( meq::RotatingSource const & source,
		std::vector<Reference::Point> const & points )
	{
		Agreement agreement;
		double sumDifference = 0.0;
		double sumReference = 0.0;

		for ( auto const & point : points )
		{
			double const mine = source.f( point.radius, point.z, point.psi );
			double const difference = mine - point.f;

			sumDifference += difference*difference;
			sumReference += point.f*point.f;

			agreement.referenceMax = std::max( agreement.referenceMax, std::fabs( point.f ) );

			// Against the case's own scale rather than the point's, so that a node
			// where geq's Jtor is near a zero crossing cannot dominate.
			if ( agreement.referenceMax > 0.0 )
				agreement.worstRelative = std::max(
					agreement.worstRelative, std::fabs( difference )/agreement.referenceMax );
		}

		agreement.referenceNorm = std::sqrt( sumReference );
		agreement.relativeL2 = agreement.referenceNorm > 0.0
			? std::sqrt( sumDifference )/agreement.referenceNorm
			: std::sqrt( sumDifference );

		return agreement;
	}

	void report( std::string const & caseName, Reference const & reference,
		Agreement const & own, Agreement const & tight )
	{
		std::cout << "\n  " << caseName << ": " << reference.species.size()
		          << " species, " << reference.points.size() << " points, R_ref = "
		          << reference.referenceRadius << "\n"
		          << std::scientific << std::setprecision( 6 )
		          << "    max |F_geq|                    " << own.referenceMax << "\n"
		          << "    relative L2, geq's own phi_0   " << own.relativeL2
		          << "   worst point " << own.worstRelative << "\n"
		          << "    relative L2, phi_0 tightened   " << tight.relativeL2
		          << "   worst point " << tight.worstRelative << "\n"
		          << "    sharpened by                   "
		          << own.relativeL2/tight.relativeL2 << "x\n";
	}
}

/*
 * THE COMPARISON ITSELF, ON BOTH OF MEQ'S CLOSURES AND AGAINST BOTH REFERENCES.
 *
 * `two` takes Closure::ClosedForm -- (97) is linear in phi_0 after logs at two
 * species -- and `three` takes Closure::RootFind, a safeguarded scalar Newton
 * with phi_0's psi-derivatives by implicit differentiation. Different code paths
 * through the same physics, checked by one external reference. The three-species
 * case also carries Z = 6, which is the charge weighting
 * `theChargeWeightedCombinationsAreNotPlainSums` exists for: every other
 * two-species configuration in the tree runs at Z = +-1, where the closed form's
 * Z_1 T_2 - Z_2 T_1 collapses to T_1 + T_2.
 *
 * THE TWO CASES SIT AT DIFFERENT MACH NUMBERS AND THAT IS FORCED, NOT CHOSEN.
 * With M^2 = m_i omega^2 R( z=0 )^2/T_e on the psi_n = 0.5 surface -- the number
 * a mirror is specified by -- `two` is at M = 6.0, which is the operating point
 * MIRROR-PLAN.md's exterior assumption needs, and `three` is at M = 0.76.
 * `meta.txt` records it per case.
 *
 * CARBON IS WHY, AND IT IS A PROPERTY OF MEQ'S GAUGE RATHER THAN OF THE
 * BENCHMARK. meq::Species::density is the PHYSICAL density on ONE curve
 * R = R_ref for every flux surface, so transferring a surface whose own midplane
 * radius is R_mid carries exp( m_s omega^2 ( R_ref^2 - R_mid^2 )/2T_s ) -- an
 * exponent linear in the species MASS. At M = 6 that is +31.2 to -51.8 for
 * deuterium and +187.2 to -310.5 for CARBON, a dynamic range of 1e216 in a
 * quantity stored as a double. The transfer's root find loses its bracket and
 * the three-species case cannot be posed at M = 6 at all. MIRROR-PLAN.md item 6
 * carries it, because a mirror is a high-Mach multi-species device and this is
 * MEQ's own interface rather than an artefact of the comparison.
 *
 * WHY THERE ARE TWO REFERENCES PER CASE, WHICH IS THE MEASUREMENT THIS FILE WAS
 * WRITTEN AROUND. Against geq's state exactly as it ships, MEQ agrees to
 * 1.324229e-03 at two species and 6.338398e-06 at three -- and a table-resolution
 * sweep says that number is NOT tabulation: at the low-M state it converges as
 * the knots are refined, 2.50e-06 at 129 knots to 2.223e-06 at 2049, and stops.
 * What it is instead was settled by asking geq about itself.
 * `solve_quasineutrality_global` converges the phi_0 FIELD to about 5.5e-05 rms
 * in Sum_s Z_s n_s; recomputing geq's OWN Jtor with phi_0 solved pointwise
 * instead moves it by 1.324204e-03 and 6.338347e-06, against MEQ's
 * disagreements of 1.324229e-03 and 6.338398e-06. FOUR AND SIX FIGURES, with
 * nothing of MEQ's in the first pair. So the whole of the disagreement is the
 * reference's own phi_0 convergence -- and it is 600x LARGER at M = 6, because
 * a stiffer quasineutrality problem is one geq's field solve converges less far
 * on. Against the tightened reference MEQ is three orders sharper at M = 6 and
 * five at M = 0.76.
 *
 * AND AGAINST THE TIGHTENED REFERENCE THERE IS NO FLOOR AT ALL -- what is left
 * is the interchange tables, and it converges away. Sweeping the exported knot
 * count with everything else fixed:
 *
 *     knots        129       257       513       1025      2049      4097
 *     two, M=0.76  1.11e-06  1.92e-07  2.84e-08  5.98e-09  1.45e-09  3.46e-10
 *     three,M=0.76 1.14e-06  1.97e-07  2.92e-08  6.17e-09  1.51e-09  3.66e-10
 *     two, M=6.0             --        2.84e-05  3.39e-06  4.45e-07  5.59e-08
 *
 * about 4x a doubling at M = 0.76 and a clean 8x -- third order -- at M = 6,
 * where the density table spans 22 orders and the interpolant has to work for
 * its accuracy. HIGH MACH DOES NOT BREAK THE AGREEMENT, IT MAKES THE TABLES
 * WORK HARDER: at 4097 knots M = 6 reaches 5.59e-08, which is what M = 0.76
 * reaches at 513. The committed references are 2049 knots for `two` and 1025 for
 * `three`, chosen for file size; the assertions below sit just above what they
 * measure and the sweep is what says they are a tabulation bound and not a
 * disagreement. SO THE TWO IMPLEMENTATIONS ARE THE SAME FUNCTION to whatever
 * accuracy the tables are asked to carry, and neither side has an error floor
 * of its own.
 *
 * THE TRANSFERABLE PART IS THAT A FLOOR IS NOT A TOLERANCE UNTIL SOMETHING
 * EXPLAINS IT. 2.2e-06 was stable, reproducible and resolution-independent, and
 * every one of those properties is equally consistent with a small fault in
 * MEQ. What separated them was a measurement of the REFERENCE against itself,
 * which costs one script and is the only experiment that can attribute a
 * residual to one side of a comparison.
 */
BOOST_AUTO_TEST_CASE( meqsRotatingSourceReproducesGeqsCurrentDensity )
{
	/*
	 * PER-CASE BOUNDS, BECAUSE THE TWO REFERENCES SIT AT DIFFERENT MACH
	 * NUMBERS AND THE INTERCHANGE'S ACCURACY IS A STEEP FUNCTION OF IT. See
	 * the block above: `two` is at M = 6.0 and `three` at M = 0.76, and the
	 * reason they differ is carbon, not a preference.
	 */
	struct Expectation
	{
		char const * name;
		double own;    ///< against geq as it ships -- the REFERENCE's floor
		double tight;  ///< against geq with phi_0 solved pointwise
		double worst;
	};

	for ( Expectation const & expected : {
		Expectation{ "two",   3.0e-3, 2.0e-6, 5.0e-6 },
		Expectation{ "three", 1.0e-5, 2.0e-8, 2.0e-7 } } )
	{
		std::string const caseName = expected.name;
		Reference const reference = readReference( caseName );

		BOOST_TEST_REQUIRE( reference.points.size() > 500u );
		BOOST_TEST_REQUIRE( reference.points.size() == reference.pointsTight.size() );
		BOOST_TEST_REQUIRE( reference.referenceRadius > 0.0 );

		std::vector<meq::Species> const species = speciesFrom( reference, false );

		// The transfer's own check, before any comparison. meq::RotatingSource
		// refuses a set that is not neutral on R = R_ref, and geq's
		// quasineutrality is what produced these densities -- so this failing
		// would indict the export rather than the source.
		double const midPsi = 0.5*( reference.psiAxis + reference.psiBoundary );
		double const neutrality = meq::chargeNeutralityResidual( species, midPsi );
		double scale = 0.0;
		for ( auto const & one : species )
			scale += std::fabs( one.charge )*( *one.density )( midPsi );
		BOOST_TEST( std::fabs( neutrality ) < 1.0e-12*scale );

		meq::RotatingSource const source(
			species, omegaFrom( reference ),
			std::make_shared<meq::ConstantProfile const>( 0.0 ),
			reference.referenceRadius, reference.mu0 );

		Agreement const own = compare( source, reference.points );
		Agreement const tight = compare( source, reference.pointsTight );
		report( caseName, reference, own, tight );

		// A comparison against a reference that is nearly zero says nothing.
		BOOST_TEST( own.referenceMax > 1.0e-3 );

		// Against geq as it ships. Loose, and it is the REFERENCE's floor: the
		// assertion below is the one with teeth.
		BOOST_TEST( own.relativeL2 < expected.own );

		// Against geq with the condition both codes share solved to the same
		// precision both codes can reach. This is the statement that two
		// independent implementations of RoPP (136) agree.
		BOOST_TEST( tight.relativeL2 < expected.tight );
		BOOST_TEST( tight.worstRelative < expected.worst );

		// And that the first number is the reference's rather than MEQ's. If a
		// fault ever appears in MEQ it moves BOTH, and this ratio collapses --
		// which is what stops the loose bound above quietly absorbing it.
		BOOST_TEST( own.relativeL2 > 1.0e2*tight.relativeL2 );
	}
}

/*
 * AND THE GAUGE TRANSFER IS LOAD BEARING, WHICH A GREEN COMPARISON CANNOT
 * ESTABLISH ON ITS OWN.
 *
 * geq's N_s and MEQ's n_s0 are one plasma in two gauges, differing by
 * exp( Z_s e delta( psi )/T_s ) -- per species, so no single rescaling relates
 * them. Were that factor near one, everything above would pass whether or not
 * anybody had thought about the gauge, and the next person to change the
 * exported state would find out the hard way.
 *
 * IT DOES NOT FAIL AS A WRONG NUMBER. It fails in the constructor:
 * Sum_s Z_s N_s is not zero on R = R_ref, because geq's gauge does not make it
 * so, and meq::RotatingSource checks exactly that and throws. Asserting the
 * THROW rather than a large error is the stronger statement -- MEQ cannot be
 * handed geq's tables by accident at all -- and it is also why
 * meq::chargeNeutralityResidual is public.
 */
BOOST_AUTO_TEST_CASE( theGaugeTransferIsLoadBearing )
{
	Reference const reference = readReference( "two" );

	double const midPsi = 0.5*( reference.psiAxis + reference.psiBoundary );
	std::vector<meq::Species> const raw = speciesFrom( reference, true );

	double rawScale = 0.0;
	for ( auto const & one : raw )
		rawScale += std::fabs( one.charge )*( *one.density )( midPsi );
	double const rawNeutrality = meq::chargeNeutralityResidual( raw, midPsi );

	std::cout << "\n  gauge: Sum_s Z_s N_s at R_ref in geq's own gauge is "
	          << std::scientific << std::setprecision( 6 )
	          << rawNeutrality << ", i.e. " << std::fabs( rawNeutrality )/rawScale
	          << " of the density scale -- not a rounding difference\n";

	BOOST_TEST( std::fabs( rawNeutrality ) > 1.0e-3*rawScale );

	BOOST_CHECK_THROW(
		meq::RotatingSource( raw, omegaFrom( reference ),
			std::make_shared<meq::ConstantProfile const>( 0.0 ),
			reference.referenceRadius, reference.mu0 ),
		std::invalid_argument );

	// And the transferred set, on the same tables, constructs and agrees. Both
	// halves are needed: a throw on its own would also be produced by a broken
	// reader.
	BOOST_CHECK_NO_THROW(
		meq::RotatingSource( speciesFrom( reference, false ), omegaFrom( reference ),
			std::make_shared<meq::ConstantProfile const>( 0.0 ),
			reference.referenceRadius, reference.mu0 ) );
}
