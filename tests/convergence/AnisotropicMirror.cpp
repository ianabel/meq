#define BOOST_TEST_MODULE AnisotropicMirror
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/Anisotropy.hpp"
#include "meq/Coils.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/RotatingSource.hpp"

#include "convergence/GeqReference.hpp"
#include "convergence/MirrorMachine.hpp"

/*
 * ANISOTROPY-PLAN.md AN-3a: THE ANISOTROPIC OPERATOR ON THE ROTATING MIRROR.
 *
 * MR-4's machine and plasma -- geq's rotating mirror, its two labels at their
 * points, the box support -- with the packs MESHED rather than subtracted, so
 * the anisotropic operator carries no coil load (ANISOTROPY-PLAN.md, What
 * changes section 2), and a POSITION Delta,
 *
 *     Delta( R, z ) = delta0 b( R ) h( z ),
 *     b = ( 1 - ( ( R - 0.225 )/0.075 )^2 )^3,  h = ( 1 - ( z/1.2 )^2 )^3,
 *
 * zero outside |R - 0.225| < 0.075, |z| < 1.2: inside the machine region, clear
 * of the box end at 1.5 and far from Gamma at rho = 3. Smooth to second order
 * at its edges, so it costs no order where it stops.
 *
 * WHY A POSITION Delta AND NOT ONE IN Psi. A Delta in Psi is a non-linear flux
 * law with no linear flux mass, and the mirror's curved Gamma needs the
 * extension, a boundary-face integrator ON that flux mass. Beside it
 * DarcyForm::Reconstruct() lifts the empty flux mass and drops the law, so
 * psi* would be wrong. Filed upstream as FLUX-LAW-BESIDE-AN-EXTENSION-FROM-MEQ.md,
 * and refused by the solver until it lands. A position Delta stays on the
 * linear flux mass and the extension sits beside it, which is what this
 * measures: the extension, the DtN and the meshed packs beside an anisotropy.
 *
 * Each case says what it asserts.
 */

namespace
{
	using namespace meqtest::mirror;

	meqtest::geq::GeqCase const &problem()
	{
		static meqtest::geq::GeqCase const p
			= meqtest::geq::readCase( meqtest::geq::inputDirectory( "rotating" ) );
		return p;
	}

	/// The bump above. Reentrant: it reads only its own constants.
	class MirrorBump : public meq::PositionAnisotropy
	{
		public:
			explicit MirrorBump( double amplitudeIn )
				: amplitude( amplitudeIn )
			{
			}

			double delta( double radius, double z ) const override
			{
				double const u = ( radius - 0.225 )/0.075;
				double const v = z/1.2;
				if ( std::abs( u ) >= 1.0 || std::abs( v ) >= 1.0 )
					return 0.0;
				double const b = 1.0 - u*u;
				double const h = 1.0 - v*v;
				return amplitude*b*b*b*h*h*h;
			}

		private:
			double amplitude;
	};

	/// Nonzero on Gamma_h, which the solver must refuse beside the extension.
	class EverywhereAnisotropy : public meq::PositionAnisotropy
	{
		public:
			double delta( double, double ) const override
			{
				return -0.1;
			}
	};

	std::vector<double> valuesAt( mfem::GridFunction const &field,
	                              std::vector<double> const &radius,
	                              std::vector<double> const &z )
	{
		int const count = static_cast<int>( radius.size() );
		mfem::Vector positions( 2*count );
		for ( int i = 0; i < count; ++i )
		{
			positions( i ) = radius[ static_cast<std::size_t>( i ) ];
			positions( count + i ) = z[ static_cast<std::size_t>( i ) ];
		}
		meq::PointLocator locator( *field.FESpace()->GetMesh() );
		locator.locate( positions );
		mfem::Vector values;
		locator.interpolate( field, values );
		for ( int i = 0; i < count; ++i )
			BOOST_TEST_REQUIRE( locator.codes()[ i ] != 2, "point " << i << " is off the mesh" );
		std::vector<double> out( static_cast<std::size_t>( count ) );
		for ( int i = 0; i < count; ++i )
			out[ static_cast<std::size_t>( i ) ] = values( i );
		return out;
	}

	/// Everything one solve aliases, kept alive with it.
	struct Run
	{
		std::unique_ptr<Domain> domain;
		std::unique_ptr<VacuumSolve> vacuum;
		std::shared_ptr<meq::NormalisedRotatingSource> plasma;
		std::unique_ptr<meq::CoilAugmentedNormalisedSource> source;
		std::unique_ptr<meq::ExteriorDtN> dtn;
		std::unique_ptr<meq::PositionAnisotropy> anisotropy;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<meq::GradShafranovSolver> solver;
		bool converged = false;
		std::string failure;
	};

	enum class Delta
	{
		None,
		Bump,
		Everywhere
	};

	std::unique_ptr<Run> run( int order, int level, Delta kind, double amplitude,
	                          bool tauFollows = false )
	{
		meqtest::geq::GeqCase const &p = problem();
		auto r = std::make_unique<Run>();
		PlasmaMesh mesh;
		mesh.axialEdge = p.axialEdge;
		mesh.aligned = true;
		mesh.alignConductors = true;
		r->domain = buildDomain( level, "anisotropic-mirror", mesh );

		r->vacuum = solveMeshed( *r->domain, order, defaultModes, machineCoils() );
		std::vector<double> const labels = valuesAt(
			r->vacuum->solver->potential(), { p.axisPointR, p.boundaryPointR },
			{ p.axisPointZ, p.boundaryPointZ } );

		r->plasma = std::make_shared<meq::NormalisedRotatingSource>(
			meqtest::geq::speciesFrom( p, false ), meqtest::geq::omegaFrom( p ),
			std::make_shared<meq::ConstantProfile const>( 0.0 ), p.referenceRadius,
			labels[ 0 ], machine().mu0 );
		r->plasma->setPlasmaSupportShape( meq::PlasmaSupportShape::Box, p.axialEdge );
		r->source = std::make_unique<meq::CoilAugmentedNormalisedSource>( r->plasma,
		                                                                  machineCoils() );

		r->dtn = std::make_unique<meq::ExteriorDtN>( 0.0, rhoGamma, defaultModes );
		r->solver = mirrorSolver( *r->domain->sub, order );
		meq::GradShafranovSolver &solver = *r->solver;
		if ( kind == Delta::Bump )
			r->anisotropy = std::make_unique<MirrorBump>( amplitude );
		else if ( kind == Delta::Everywhere )
			r->anisotropy = std::make_unique<EverywhereAnisotropy>();
		if ( r->anisotropy )
			solver.setAnisotropy( *r->anisotropy );
		if ( tauFollows )
			solver.setStabilisationFollowsAnisotropy( true );
		solver.setAxisFluxPoint( p.axisPointR, p.axisPointZ );
		solver.setBoundaryFluxPoint( p.boundaryPointR, p.boundaryPointZ );
		solver.setBoundaryFluxInitial( labels[ 1 ] );
		solver.setSource( *r->source, labels[ 0 ] );
		solver.setBoundaryData( r->zero );
		solver.setExtension( *r->domain->path, r->domain->gammaHMarker );
		solver.setExteriorCoupling( *r->dtn );
		solver.setInitialGuess( r->vacuum->solver->potential() );
		try
		{
			solver.solve();
			r->converged = true;
		}
		catch ( std::exception const &e )
		{
			r->failure = e.what();
		}
		return r;
	}

	/// geq's grid nodes in the machine region, coarsened: the points every
	/// level is sampled at.
	struct Cloud
	{
		std::vector<double> radius;
		std::vector<double> z;
	};

	Cloud machineCloud()
	{
		Cloud c;
		for ( double radius = 0.02; radius <= 0.62; radius += 0.04 )
			for ( double z = -1.3; z <= 1.3; z += 0.1 )
				if ( inMachineRegion( radius, z ) )
				{
					c.radius.push_back( radius );
					c.z.push_back( z );
				}
		return c;
	}

	double rms( std::vector<double> const &a, std::vector<double> const &b )
	{
		double sum = 0.0;
		for ( std::size_t i = 0; i < a.size(); ++i )
			sum += ( a[ i ] - b[ i ] )*( a[ i ] - b[ i ] );
		return std::sqrt( sum/static_cast<double>( a.size() ) );
	}
}

/// Delta = 0 through the interface is the meshed mirror to the bit: the
/// position path's weight is R times exactly one, the extension's coefficient
/// is R either way, and nothing else moves.
BOOST_AUTO_TEST_CASE( aZeroBumpIsTheMeshedMirrorBitForBit )
{
	auto const plain = run( 2, 0, Delta::None, 0.0 );
	auto const zero = run( 2, 0, Delta::Bump, 0.0 );
	BOOST_TEST_REQUIRE( plain->converged, plain->failure );
	BOOST_TEST_REQUIRE( zero->converged, zero->failure );

	double worst = 0.0;
	mfem::GridFunction const &a = plain->solver->potential();
	mfem::GridFunction const &b = zero->solver->potential();
	for ( int i = 0; i < a.Size(); ++i )
		worst = std::max( worst, std::abs( a( i ) - b( i ) ) );
	BOOST_TEST( worst == 0.0, "psi differs by " << worst );
	BOOST_TEST( plain->solver->newtonIterations() == zero->solver->newtonIterations() );
}

/// From the meshed vacuum, at two amplitudes a pickup plasma might have, the
/// bordered Newton converges in as few steps as the isotropic mirror's 3 or
/// 4. Three residuals are too few to read an order from, so none is asserted.
BOOST_AUTO_TEST_CASE( theAnisotropicMirrorConvergesFromTheVacuum )
{
	for ( double amplitude : { -0.1, -0.3 } )
	{
		auto const r = run( 2, 0, Delta::Bump, amplitude );
		BOOST_TEST_REQUIRE( r->converged, "delta0 = " << amplitude << ": " << r->failure );
		std::vector<double> const &history = r->solver->newtonResiduals();
		std::printf( "  delta0 = %+.1f: %d Newton steps, psi_ax %.6e, psi_bnd %.6e; residuals",
		             amplitude, r->solver->newtonIterations(), r->solver->psiAxis(),
		             r->solver->psiBoundary() );
		for ( double v : history )
			std::printf( " %.2e", v );
		std::printf( "\n" );
		std::fflush( stdout );
		BOOST_TEST( r->solver->newtonIterations() <= 8,
		            "delta0 = " << amplitude << ": " << r->solver->newtonIterations()
		            << " Newton steps against the isotropic mirror's 3 or 4" );
	}
}

/// MEQ AGAINST ITSELF, AND AGAINST THE ISOTROPIC MESHED MIRROR ON THE SAME LADDER:
/// k = 3, levels 0 to 3, psi* sampled over the machine region.
///
/// THE LADDER ITSELF IS NOT AT THE DESIGN RATE, AND THE ISOTROPIC CONTROL SAYS
/// SO. Its Cauchy rates read 7.28 then 2.35 against psi*'s design 4: level 0
/// pre-asymptotic, and then limited -- most likely by the meshed packs'
/// corners, r^2 log r at each rectangle's corner, which RotatingMirror's
/// subtracted filaments do not have. Not measured. So what is asserted is
/// what the anisotropy adds: on every rung its difference is within 25% of the
/// isotropic one, and falls. Measured, the two agree to 6%.
///
/// And tau( 1 - Delta ), setStabilisationFollowsAnisotropy(), beside it: the
/// measurement AN-3a owes on tau. At delta0 = -0.3 it moves the differences
/// by under 1e-3 of themselves, printed and not asserted.
BOOST_AUTO_TEST_CASE( theAnisotropicMirrorConvergesAsTheIsotropicOneDoes )
{
	Cloud const cloud = machineCloud();
	struct Variant
	{
		char const *label;
		Delta kind;
		bool tauFollows;
	};
	std::vector<std::vector<double>> differences;
	for ( Variant const &variant : { Variant{ "isotropic", Delta::None, false },
	                                 Variant{ "delta0 = -0.3", Delta::Bump, false },
	                                 Variant{ "delta0 = -0.3, tau( 1 - Delta )", Delta::Bump, true } } )
	{
		std::vector<std::vector<double>> sampled;
		for ( int level : { 0, 1, 2, 3 } )
		{
			auto r = run( 3, level, variant.kind, -0.3, variant.tauFollows );
			BOOST_TEST_REQUIRE( r->converged, "level " << level << ": " << r->failure );
			r->solver->postProcess();
			sampled.push_back( valuesAt( r->solver->postProcessedPotential(), cloud.radius,
			                             cloud.z ) );
		}
		std::vector<double> d;
		for ( std::size_t i = 1; i < sampled.size(); ++i )
			d.push_back( rms( sampled[ i - 1 ], sampled[ i ] ) );
		std::printf( "  %s, k = 3: psi* level 0-1 %.3e, 1-2 %.3e, 2-3 %.3e; rates %.3f, %.3f;"
		             " over %zu points\n", variant.label, d[ 0 ], d[ 1 ], d[ 2 ],
		             std::log2( d[ 0 ]/d[ 1 ] ), std::log2( d[ 1 ]/d[ 2 ] ),
		             cloud.radius.size() );
		std::fflush( stdout );
		differences.push_back( d );
	}

	std::vector<double> const &isotropic = differences[ 0 ];
	for ( std::size_t v = 1; v < differences.size(); ++v )
		for ( std::size_t i = 0; i < isotropic.size(); ++i )
		{
			double const ratio = differences[ v ][ i ]/isotropic[ i ];
			BOOST_TEST( ( ratio > 0.8 && ratio < 1.25 ), "variant " << v << ", rung " << i
			            << ": the anisotropic difference is " << ratio
			            << " of the isotropic one on the same ladder" );
			if ( i > 0 )
				BOOST_TEST( differences[ v ][ i ] < differences[ v ][ i - 1 ] );
		}
}

/// A Delta that does not vanish on Gamma_h is refused beside the extension:
/// the extension reads q there as grad psi/R, which it is only where Delta = 0.
BOOST_AUTO_TEST_CASE( aDeltaOnGammaIsRefusedBesideTheExtension )
{
	auto const r = run( 1, 0, Delta::Everywhere, 0.0 );
	BOOST_TEST( !r->converged );
	BOOST_TEST( r->failure.find( "Gamma_h" ) != std::string::npos, r->failure );
}
