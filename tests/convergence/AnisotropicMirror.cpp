#define BOOST_TEST_MODULE AnisotropicMirror
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "mfem.hpp"

#include "meq/Anisotropy.hpp"
#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/RotatingSource.hpp"

#include "convergence/GeqReference.hpp"
#include "convergence/MirrorMachine.hpp"
#include "convergence/RoundOffFloor.hpp"

/*
 * ANISOTROPY-PLAN.md AN-3a: THE ANISOTROPIC OPERATOR ON THE ROTATING MIRROR.
 *
 * MR-4's machine and plasma -- geq's rotating mirror, its two labels at their
 * points, the box support -- with the packs MESHED rather than subtracted, so
 * the anisotropic operator carries no coil load (ANISOTROPY-PLAN.md, What
 * changes section 2), and two Deltas, each zero well inside Gamma at rho = 3
 * and smooth to second order where it stops, so it costs no order there.
 *
 * A POSITION Delta, on the linear flux mass beside the extension's face
 * integrator:
 *
 *     Delta( R, z ) = delta0 b( R ) h( z ),
 *     b = ( 1 - ( ( R - 0.225 )/0.075 )^2 )^3,  h = ( 1 - ( z/1.2 )^2 )^3,
 *
 * zero outside |R - 0.225| < 0.075, |z| < 1.2: inside the machine region and
 * clear of the box end at 1.5.
 *
 * A Delta IN Psi, the shape AN-4 hands geq, a flux law with the extension's
 * face integrator ALONE on the flux mass form:
 *
 *     Delta( Psi, z ) = delta0 g( Psi ) h( z ),  g = ( 4 Psi ( 1 - Psi ) )^3,
 *
 * zero outside 0 < Psi < 1, between the two labels. Whether it vanishes on
 * Gamma_h depends on the solved psi, so the solver checks it there.
 *
 * A Delta THAT READS B, ANISOTROPY-PLAN.md AN-6a: the same shape as a
 * numerator, mu0 Delta_p = B_0^2 delta0 g( Psi ) h( z ), and Delta formed at
 * the iterate's B. Its delta0 is a pressure and not a Delta: the solved Delta
 * is deeper where B falls, and at fixed flux mu0 Delta_p < -|q|^2/4 has no
 * root at all -- the fold where Delta = -1 and the operator, whose
 * coefficient along grad psi is 1 + Delta, stops being elliptic.
 *
 * AND THE SAME PACKS SUBTRACTED, ANISOTROPY-PLAN.md AN-3b: psi_c of the two
 * rectangles handed to setConductorField(), the plasma source alone on the
 * mesh, and the anisotropic flux row carrying the conductors' load. Both
 * routes solve one continuous problem on one mesh, so their difference falls
 * with refinement, which is the meshed route's cross-check of the split. A
 * load missing or of the wrong sign leaves a difference of order delta0 that
 * does not.
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

	/// The Delta in Psi above: AN-4's, GeqReference.hpp's one definition of
	/// it, at the export's cut.
	class MirrorProfile : public meqtest::geq::ReferenceAnisotropy
	{
		public:
			explicit MirrorProfile( double amplitudeIn )
				: meqtest::geq::ReferenceAnisotropy( amplitudeIn, 1.2 )
			{
			}
	};

	/// The same shape as a numerator, ANISOTROPY-PLAN.md AN-6a's Delta( B ):
	/// mu0 Delta_p = B_0^2 delta0 g( Psi ) h( z ), with B_0^2 = 1.12 T^2 the
	/// isotropic meshed mirror's |B|^2 at the midplane, where B = B_z is
	/// 1.05 to 1.06 T at R = 0.22 to 0.26 (measured, level 0, k = 2). So Delta
	/// is about delta0 there and falls as 1/B^2 towards the throats, where B_z
	/// alone reaches 3.6 T at z = 1.2 -- a pressure held on a flux surface, as
	/// a kinetic closure hands it, rather than a Delta held.
	class MirrorFieldProfile : public meq::FieldAnisotropy
	{
		public:
			explicit MirrorFieldProfile( double amplitudeIn )
				: shape( amplitudeIn, 1.2 )
			{
			}

			double pressure( double radius, double z, double psi ) const override
			{
				return midplaneFieldSquared*shape.delta( radius, z, psi );
			}

			double dPressureDPsi( double radius, double z, double psi ) const override
			{
				return midplaneFieldSquared*shape.dDeltaDPsi( radius, z, psi );
			}

			Argument argument() const override
			{
				return shape.argument();
			}

			static constexpr double midplaneFieldSquared = 1.12;

		private:
			meqtest::geq::ReferenceAnisotropy shape;
	};

	/// A Delta in Psi that is nonzero everywhere, on Gamma_h included: only
	/// the solved psi can show it, and the solver must refuse that answer.
	class EverywhereInPsi : public meq::Anisotropy
	{
		public:
			double delta( double, double, double ) const override
			{
				return -0.1;
			}

			double dDeltaDPsi( double, double, double ) const override
			{
				return 0.0;
			}

			Argument argument() const override
			{
				return Argument::NormalisedFlux;
			}
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

	/// The two packs as psi_c: the subtracted route's conductors, the meshed
	/// route's machineCoils() at the same quadrature order.
	meq::ConductorField packsField()
	{
		meq::ConductorField field( machine().mu0 );
		field.add( pack( -1.0 ) );
		field.add( pack( +1.0 ) );
		field.setQuadratureOrder( referenceOrder );
		return field;
	}

	/// Everything one solve aliases, kept alive with it.
	struct Run
	{
		std::unique_ptr<meq::ConductorField> conductors;
		std::unique_ptr<Domain> domain;
		std::unique_ptr<VacuumSolve> vacuum;
		std::shared_ptr<meq::NormalisedRotatingSource> plasma;
		std::unique_ptr<meq::CoilAugmentedNormalisedSource> source;
		std::unique_ptr<meq::ExteriorDtN> dtn;
		std::unique_ptr<meq::Anisotropy> anisotropy;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<meq::GradShafranovSolver> solver;
		bool converged = false;
		std::string failure;
	};

	enum class Delta
	{
		None,
		Bump,
		Everywhere,
		Profile,
		EverywhereInPsi,
		Field
	};

	char const *kindName( Delta kind )
	{
		switch ( kind )
		{
			case Delta::None: return "isotropic";
			case Delta::Bump: return "bump";
			case Delta::Everywhere: return "everywhere";
			case Delta::Profile: return "Delta( Psi )";
			case Delta::EverywhereInPsi: return "everywhere in Psi";
			case Delta::Field: return "Delta( B )";
		}
		return "";
	}

	std::unique_ptr<Run> run( int order, int level, Delta kind, double amplitude,
	                          bool tauFollows = false, bool subtracted = false )
	{
		meqtest::geq::GeqCase const &p = problem();
		auto r = std::make_unique<Run>();
		PlasmaMesh mesh;
		mesh.axialEdge = p.axialEdge;
		mesh.aligned = true;
		mesh.alignConductors = true;
		r->domain = buildDomain( level, "anisotropic-mirror", mesh );

		std::vector<double> labels;
		if ( subtracted )
		{
			// The labels are the TOTAL flux at their points: the vacuum's
			// remainder there plus psi_c, exactly.
			r->conductors = std::make_unique<meq::ConductorField>( packsField() );
			r->vacuum = solveSubtracted( *r->domain, order, defaultModes, *r->conductors );
			labels = valuesAt( r->vacuum->solver->potential(), { p.axisPointR, p.boundaryPointR },
			                   { p.axisPointZ, p.boundaryPointZ } );
			labels[ 0 ] += r->conductors->psi( p.axisPointR, p.axisPointZ );
			labels[ 1 ] += r->conductors->psi( p.boundaryPointR, p.boundaryPointZ );
		}
		else
		{
			r->vacuum = solveMeshed( *r->domain, order, defaultModes, machineCoils() );
			labels = valuesAt( r->vacuum->solver->potential(), { p.axisPointR, p.boundaryPointR },
			                   { p.axisPointZ, p.boundaryPointZ } );
		}

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
		else if ( kind == Delta::Profile )
			r->anisotropy = std::make_unique<MirrorProfile>( amplitude );
		else if ( kind == Delta::EverywhereInPsi )
			r->anisotropy = std::make_unique<EverywhereInPsi>();
		else if ( kind == Delta::Field )
			r->anisotropy = std::make_unique<MirrorFieldProfile>( amplitude );
		if ( r->anisotropy )
			solver.setAnisotropy( *r->anisotropy );
		if ( tauFollows )
			solver.setStabilisationFollowsAnisotropy( true );
		solver.setAxisFluxPoint( p.axisPointR, p.axisPointZ );
		solver.setBoundaryFluxPoint( p.boundaryPointR, p.boundaryPointZ );
		solver.setBoundaryFluxInitial( labels[ 1 ] );
		if ( subtracted )
		{
			solver.setConductorField( *r->conductors );
			solver.setSource( *r->plasma, labels[ 0 ] );
		}
		else
			solver.setSource( *r->source, labels[ 0 ] );
		solver.setBoundaryData( r->zero );
		solver.setExtension( *r->domain->path, r->domain->gammaHMarker );
		solver.setExteriorCoupling( *r->dtn );
		solver.setInitialGuess( r->vacuum->solver->potential() );
#ifdef MFEM_USE_SINGLE
		// THE MESHED PACKS' SINGLE FLOOR IS THE PLASMA SOLVE'S TOO, and it
		// sits above the ten-eps gate the round-off stop works under, as
		// solveMeshed()'s vacuum does (M-294): measured with no target to
		// meet, every meshed plasma solve here stagnates between 1.35e-06 and
		// 2.20e-06 -- k = 2 at levels 0 and 1, k = 3 at levels 0 to 2, every
		// Delta -- drifting down by a few per cent a step, while the
		// subtracted route reaches 1.6e-08 to 6.9e-08. Whether a solve then
		// dips under the gate and stops is a coin toss: k = 3's isotropic
		// level 0 and k = 2's level 1 do not. solveMeshed()'s 64 eps, 7.6e-06, clears the worst stagnant
		// residual by 3.5 and is four orders under the step before it.
		// Double asks the solver's own default. MEASUREMENTS.md M-350.
		if ( !subtracted )
			solver.setNewtonControl( 1.0e-12,
			                         64.0*static_cast<double>(
				                         std::numeric_limits<mfem::real_t>::epsilon() ),
			                         30 );
#endif
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

	/// psi* at the cloud, the TOTAL on the subtracted route: the remainder's
	/// reconstruction plus psi_c, exactly.
	std::vector<double> totalPsiStarAt( Run &r, Cloud const &cloud )
	{
		r.solver->postProcess();
		std::vector<double> values = valuesAt( r.solver->postProcessedPotential(), cloud.radius,
		                                       cloud.z );
		if ( r.conductors )
			for ( std::size_t i = 0; i < values.size(); ++i )
				values[ i ] += r.conductors->psi( cloud.radius[ i ], cloud.z[ i ] );
		return values;
	}

	/// The solved Delta( B )'s range over the machine cloud, meshed route: q
	/// is the total flux there, and Delta is meq::fieldDelta() of the
	/// numerator at the solved Psi and |q|^2.
	std::pair<double, double> solvedFieldDelta( Run &r )
	{
		Cloud const cloud = machineCloud();
		mfem::GridFunction const &q = r.solver->flux();
		int const count = static_cast<int>( cloud.radius.size() );
		mfem::Vector positions( 2*count );
		for ( int i = 0; i < count; ++i )
		{
			positions( i ) = cloud.radius[ static_cast<std::size_t>( i ) ];
			positions( count + i ) = cloud.z[ static_cast<std::size_t>( i ) ];
		}
		// meq::PointLocator interpolates a scalar field, so each component
		// of the vector L2 flux is read through a scalar view of its block.
		mfem::FiniteElementSpace const &space = *q.FESpace();
		BOOST_TEST_REQUIRE( space.GetOrdering() == mfem::Ordering::byNODES );
		mfem::FiniteElementSpace scalar( space.GetMesh(), space.FEColl(), 1 );
		meq::PointLocator locator( *space.GetMesh() );
		locator.locate( positions );
		mfem::Vector component[ 2 ];
		for ( int c = 0; c < 2; ++c )
		{
			mfem::GridFunction view;
			view.MakeRef( &scalar, const_cast<mfem::real_t *>( q.GetData() ) + c*scalar.GetVSize() );
			locator.interpolate( view, component[ c ] );
		}
		for ( int i = 0; i < count; ++i )
			BOOST_TEST_REQUIRE( locator.codes()[ i ] != 2, "point " << i << " is off the mesh" );
		std::vector<double> const psi = valuesAt( r.solver->potential(), cloud.radius, cloud.z );
		auto const &field = dynamic_cast<meq::FieldAnisotropy const &>( *r.anisotropy );
		double const axis = r.solver->psiAxis();
		double const boundary = r.solver->psiBoundary();
		std::pair<double, double> range{ 0.0, 0.0 };
		for ( int i = 0; i < count; ++i )
		{
			std::size_t const k = static_cast<std::size_t>( i );
			double const pressure = field.pressure( cloud.radius[ k ], cloud.z[ k ],
			                                        ( psi[ k ] - boundary )/( axis - boundary ) );
			double const fluxSquared = component[ 0 ]( i )*component[ 0 ]( i )
			                           + component[ 1 ]( i )*component[ 1 ]( i );
			double const delta = meq::fieldDelta( pressure, fluxSquared ).delta;
			range.first = std::min( range.first, delta );
			range.second = std::max( range.second, delta );
		}
		return range;
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
		worst = std::max( worst, static_cast<double>( std::abs( a( i ) - b( i ) ) ) );
	BOOST_TEST( worst == 0.0, "psi differs by " << worst );
	BOOST_TEST( plain->solver->newtonIterations() == zero->solver->newtonIterations() );
}

/// From the meshed vacuum, at two amplitudes a pickup plasma might have, the
/// bordered Newton converges in as few steps as the isotropic mirror's 3 or
/// 4, for either Delta. Three residuals are too few to read an order from, so
/// none is asserted.
///
/// And Delta( B ) at delta0 = -0.1 and -0.2: a pressure, whose Delta is
/// printed over the machine region. It reaches -0.124 and -0.338, the second
/// near the -0.37 mirrorfp's GF1 table has at its midplane.
BOOST_AUTO_TEST_CASE( theAnisotropicMirrorConvergesFromTheVacuum )
{
	struct Case
	{
		Delta kind;
		double amplitude;
	};
	for ( Case const c : { Case{ Delta::Bump, -0.1 }, Case{ Delta::Bump, -0.3 },
	                       Case{ Delta::Profile, -0.1 }, Case{ Delta::Profile, -0.3 },
	                       Case{ Delta::Field, -0.1 }, Case{ Delta::Field, -0.2 } } )
	{
		Delta const kind = c.kind;
		double const amplitude = c.amplitude;
		char const *const label = kindName( kind );
		auto const r = run( 2, 0, kind, amplitude );
		BOOST_TEST_REQUIRE( r->converged, label << ", delta0 = " << amplitude << ": " << r->failure );
		std::vector<double> const &history = r->solver->newtonResiduals();
		std::printf( "  %s, delta0 = %+.1f: %d Newton steps, psi_ax %.6e, psi_bnd %.6e; residuals",
		             label, amplitude, r->solver->newtonIterations(), r->solver->psiAxis(),
		             r->solver->psiBoundary() );
		for ( double v : history )
			std::printf( " %.2e", v );
		std::printf( "\n" );
		if ( kind == Delta::Field )
		{
			std::pair<double, double> const range = solvedFieldDelta( *r );
			std::printf( "    solved Delta in [ %.3f, %.3f ] over the machine region\n",
			             range.first, range.second );
		}
		std::fflush( stdout );
		BOOST_TEST( r->solver->newtonIterations() <= 8,
		            label << ", delta0 = " << amplitude << ": " << r->solver->newtonIterations()
		            << " Newton steps against the isotropic mirror's 3 or 4" );
	}
}

/// A pressure past the fold has no equilibrium at this flux, and the solver
/// says so on the solved answer rather than returning one: at delta0 = -0.3,
/// mu0 Delta_p reaches -0.336 T^2 where |q|^2 is near 1.29 T^2, under the
/// 4|mu0 Delta_p| the root needs. -0.24 converges in 3 Newton steps; the line
/// is near -0.25, where |q|^2/4 meets the pressure.
BOOST_AUTO_TEST_CASE( aPressurePastTheFoldIsRefused )
{
	auto const r = run( 2, 0, Delta::Field, -0.3 );
	BOOST_TEST( !r->converged );
	BOOST_TEST( r->failure.find( "fold" ) != std::string::npos, r->failure );
}

/// MEQ AGAINST ITSELF, AND AGAINST THE ISOTROPIC MESHED MIRROR ON THE SAME LADDER:
/// k = 3, levels 0 to 2, psi* sampled over the machine region.
///
/// THE LADDER ITSELF IS NOT AT THE DESIGN RATE, AND THE ISOTROPIC CONTROL SAYS
/// SO. Its Cauchy rate reads 7.28 against psi*'s design 4, level 0
/// pre-asymptotic; a level 3 reads 2.35 on the next rung, limited -- most
/// likely by the meshed packs' corners, r^2 log r at each rectangle's corner,
/// which RotatingMirror's subtracted filaments do not have. Not measured.
/// Level 3 is not run: it is six k = 3 solves at 64 times level 0's elements,
/// most of the suite's wall, for a rung that adds no claim. So what is
/// asserted is what the anisotropy adds: on every rung its difference is
/// within 25% of the isotropic one, and falls. Measured, the two agree to 6%.
///
/// And tau( 1 - Delta ), setStabilisationFollowsAnisotropy(), beside it: the
/// measurement AN-3a owes on tau. At delta0 = -0.3 it moves the differences
/// by under 1e-3 of themselves, printed and not asserted.
///
/// AND THE Delta IN Psi, on the same ladder, which is the case with teeth for
/// the reconstruction: psi* is read, and beside the extension's face
/// integrator a dropped flux law moves it by order one.
///
/// AND Delta( B ), AN-6a, at a pressure either side: -0.2, whose solved Delta
/// is mirrorfp's depth, and +0.9. Its psi* is the reconstruction with Delta
/// frozen at the solved flux.
BOOST_AUTO_TEST_CASE( theAnisotropicMirrorConvergesAsTheIsotropicOneDoes )
{
	Cloud const cloud = machineCloud();
	struct Variant
	{
		char const *label;
		Delta kind;
		bool tauFollows;
		double amplitude;
	};
	std::vector<std::vector<double>> differences;
	double scale = 0.0;
	for ( Variant const &variant : { Variant{ "isotropic", Delta::None, false, 0.0 },
	                                 Variant{ "delta0 = -0.3", Delta::Bump, false, -0.3 },
	                                 Variant{ "delta0 = -0.3, tau( 1 - Delta )", Delta::Bump, true, -0.3 },
	                                 Variant{ "Delta( Psi ), delta0 = -0.3", Delta::Profile, false, -0.3 },
	                                 Variant{ "Delta( B ), delta0 = -0.2", Delta::Field, false, -0.2 },
	                                 Variant{ "Delta( B ), delta0 = +0.9", Delta::Field, false, 0.9 } } )
	{
		std::vector<std::vector<double>> sampled;
		for ( int level : { 0, 1, 2 } )
		{
			auto r = run( 3, level, variant.kind, variant.amplitude, variant.tauFollows );
			BOOST_TEST_REQUIRE( r->converged, "level " << level << ": " << r->failure );
			r->solver->postProcess();
			sampled.push_back( valuesAt( r->solver->postProcessedPotential(), cloud.radius,
			                             cloud.z ) );
		}
		std::vector<double> d;
		for ( std::size_t i = 1; i < sampled.size(); ++i )
			d.push_back( rms( sampled[ i - 1 ], sampled[ i ] ) );
		if ( differences.empty() )
			scale = rms( sampled.back(), std::vector<double>( sampled.back().size(), 0.0 ) );
		std::printf( "  %s, k = 3: psi* level 0-1 %.3e, 1-2 %.3e; rate %.3f; over %zu points\n",
		             variant.label, d[ 0 ], d[ 1 ], std::log2( d[ 0 ]/d[ 1 ] ),
		             cloud.radius.size() );
		std::fflush( stdout );
		differences.push_back( d );
	}

	/*
	 * IN SINGLE THE LADDER REACHES THE FORMAT'S FLOOR AFTER ITS FIRST RUNG.
	 * The isotropic differences read 5.22e-05, 1.42e-06 against double's
	 * 5.10e-05, 3.28e-07: rung 0 is the same number to 2.4%, rung 1 the
	 * meshed plasma solve's own round-off -- its residual stagnates at
	 * 1.4e-06 to 2.2e-06. In units of psi*'s rms over the cloud, 0.22, the floored
	 * rungs of every variant read 36 to 63 eps and rung 0 1970. A rung whose
	 * isotropic difference is under 256 eps of that rms compares two floors,
	 * so it is skipped loudly; rung 0 is asserted as in double.
	 * roundOffFloor()'s 2048 eps would skip rung 0 as well and leave the case
	 * asserting nothing. MEASUREMENTS.md M-350.
	 */
	std::vector<double> const &isotropic = differences[ 0 ];
	std::vector<bool> readable( isotropic.size(), true );
	double const ladderFloor = 256.0;
	for ( std::size_t i = 0; i < isotropic.size(); ++i )
		readable[ i ] = meq::tests::aboveTheFormatsFloor( isotropic[ i ], scale, ladderFloor );
#ifdef MFEM_USE_SINGLE
	std::printf( "  psi* rms over the cloud %.3e; this build's floor on a rung %.3e\n", scale,
	             ladderFloor*scale*static_cast<double>(
		             std::numeric_limits<mfem::real_t>::epsilon() ) );
	for ( std::size_t i = 0; i < isotropic.size(); ++i )
	{
		if ( !readable[ i ] )
			std::printf( "    *** SKIPPED rung %zu: the isotropic difference %.3e is under this"
			             " build's floor, so the rung compares two floors. M-350\n", i,
			             isotropic[ i ] );
	}
#endif
	for ( std::size_t v = 1; v < differences.size(); ++v )
		for ( std::size_t i = 0; i < isotropic.size(); ++i )
		{
			if ( !readable[ i ] )
				continue;
			double const ratio = differences[ v ][ i ]/isotropic[ i ];
			BOOST_TEST( ( ratio > 0.8 && ratio < 1.25 ), "variant " << v << ", rung " << i
			            << ": the anisotropic difference is " << ratio
			            << " of the isotropic one on the same ladder" );
			if ( i > 0 )
				BOOST_TEST( differences[ v ][ i ] < differences[ v ][ i - 1 ] );
		}
}

/// ANISOTROPY-PLAN.md AN-6b, Delta > 0 TOWARDS THE FIREHOSE LIMIT, from the
/// meshed vacuum at level 0, k = 2. MEQ against itself only: geq has no
/// converged reference there.
///
///  - the bump, a Delta of position, to delta0 = 0.99, where 1 - Delta is a
///    hundredth at its peak;
///  - Delta( B ) to a pressure of 3.0, three times B_0^2: its root never
///    reaches 1 at a positive flux, so the limit is not reachable at all;
///  - Delta( Psi ) to 0.7. Past it the Newton crawls, 12 steps at 0.9 and 23
///    at 0.95, and at 0.99 it does not converge. That is a defect, asserted
///    -- and red -- by theDeltaInPsiConvergesTowardsTheFirehoseLimitAsTheBumpDoes
///    below; MEASUREMENTS.md M-351.
///
/// Each within the isotropic mirror's 3 or 4 steps and some slack, as the
/// negative side is.
BOOST_AUTO_TEST_CASE( theAnisotropicMirrorConvergesTowardsTheFirehoseLimit )
{
	struct Case
	{
		Delta kind;
		double amplitude;
	};
	for ( Case const c : { Case{ Delta::Bump, 0.5 }, Case{ Delta::Bump, 0.9 },
	                       Case{ Delta::Bump, 0.99 }, Case{ Delta::Field, 0.9 },
	                       Case{ Delta::Field, 3.0 }, Case{ Delta::Profile, 0.5 },
	                       Case{ Delta::Profile, 0.7 } } )
	{
		auto const r = run( 2, 0, c.kind, c.amplitude );
		BOOST_TEST_REQUIRE( r->converged, kindName( c.kind ) << ", delta0 = " << c.amplitude
		                    << ": " << r->failure );
		std::printf( "  %s, delta0 = %+.2f: %d Newton steps, psi_ax %.6e, psi_bnd %.6e;"
		             " residuals", kindName( c.kind ), c.amplitude, r->solver->newtonIterations(),
		             r->solver->psiAxis(), r->solver->psiBoundary() );
		for ( double v : r->solver->newtonResiduals() )
			std::printf( " %.2e", v );
		std::printf( "\n" );
		if ( c.kind == Delta::Field )
		{
			std::pair<double, double> const range = solvedFieldDelta( *r );
			std::printf( "    solved Delta in [ %.3f, %.3f ] over the machine region\n",
			             range.first, range.second );
		}
		std::fflush( stdout );
		BOOST_TEST( r->solver->newtonIterations() <= 8, kindName( c.kind ) << ", delta0 = "
		            << c.amplitude << ": " << r->solver->newtonIterations() << " Newton steps" );
	}
}

/// DELTA( Psi ) TOWARDS THE FIREHOSE LIMIT, AS THE BUMP GOES: delta0 = 0.9
/// and 0.95 from the same meshed vacuum, level 0, k = 2, within the same
/// budget the neighbouring cases allow, as 0.7 and below converge.
///
/// RED, DELIBERATELY, WHILE THE CRAWL STANDS (MEASUREMENTS.md M-351; TODO,
/// "Anisotropic pressure"). Measured in double 12 and 23 steps, in single 10
/// and no convergence. The bordered Newton direction is exact along every
/// step and continuation in delta0 reaches 0.99 undamped, so neither the
/// Jacobian nor a fold is at fault: the Armijo search walks exact steps whose
/// curvature, through the boundary label's column, reaches 3e4 times the
/// bump's. Turns green with a fix of that, not with a looser budget.
BOOST_AUTO_TEST_CASE( theDeltaInPsiConvergesTowardsTheFirehoseLimitAsTheBumpDoes )
{
	for ( double amplitude : { 0.9, 0.95 } )
	{
		auto const r = run( 2, 0, Delta::Profile, amplitude );
		std::printf( "  Delta( Psi ), delta0 = %+.2f: %s, %d Newton steps; dampings",
		             amplitude, r->converged ? "converged" : "FAILED",
		             r->solver->newtonIterations() );
		for ( auto const &step : r->solver->borderSteps() )
			std::printf( " %.3g", step.damping );
		std::printf( "\n" );
		std::fflush( stdout );
		bool const met = r->converged && r->solver->newtonIterations() <= 8;
		BOOST_TEST( met, "Delta( Psi ), delta0 = " << amplitude << ": "
		            << ( r->converged ? "converged in " : "did not converge in " )
		            << r->solver->newtonIterations() << " Newton steps against the bump's 3"
		            " or 4 and this suite's budget of 8 -- the crawl towards the firehose"
		            " limit, MEASUREMENTS.md M-351, TODO \"Anisotropic pressure\". The"
		            " Jacobian is exact; the line search walks steps whose curvature runs"
		            " through the psi_bnd column. Candidate fixes, unmeasured: carry 1 - Delta"
		            " on the conservation row rather than as R/( 1 - Delta ) on the flux"
		            " mass, or bound the psi_bnd step relative to the label span"
		            << ( r->converged ? "" : std::string( ". " ) + r->failure ) );
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

/// A Delta in Psi that does not vanish on Gamma_h is refused too, on the
/// solved psi -- it cannot be checked before the solve -- and by name.
BOOST_AUTO_TEST_CASE( aDeltaInPsiOnGammaIsRefusedOnTheSolvedPsi )
{
	auto const r = run( 1, 0, Delta::EverywhereInPsi, 0.0 );
	BOOST_TEST( !r->converged );
	BOOST_TEST( r->failure.find( "Gamma_h" ) != std::string::npos, r->failure );
	BOOST_TEST( r->failure.find( "solved psi" ) != std::string::npos, r->failure );
}

/*
 * ---------------------------------------------------------------------------
 * AN-3b: the subtracted mirror
 * ---------------------------------------------------------------------------
 */

/// A zero bump on the subtracted route goes through the flux law -- every Delta
/// beside setConductorField() does -- and must be the isotropic subtracted
/// mirror, to the solver's tolerance rather than to the bit: the two assemble
/// the flux mass through different integrators.
BOOST_AUTO_TEST_CASE( aZeroBumpThroughTheLawIsTheSubtractedMirror )
{
	auto const plain = run( 2, 0, Delta::None, 0.0, false, true );
	auto const zero = run( 2, 0, Delta::Bump, 0.0, false, true );
	BOOST_TEST_REQUIRE( plain->converged, plain->failure );
	BOOST_TEST_REQUIRE( zero->converged, zero->failure );

	mfem::GridFunction const &a = plain->solver->potential();
	mfem::GridFunction const &b = zero->solver->potential();
	double worst = 0.0;
	for ( int i = 0; i < a.Size(); ++i )
		worst = std::max( worst, static_cast<double>( std::abs( a( i ) - b( i ) ) ) );
	double const scale = a.Normlinf();
	std::printf( "  subtracted, zero bump through the law: psi_p differs by %.3e on %.3e;"
	             " Newton %d and %d\n", worst, scale, plain->solver->newtonIterations(),
	             zero->solver->newtonIterations() );
	std::fflush( stdout );
	double bound = 1.0e-10*scale;
#ifdef MFEM_USE_SINGLE
	// psi_p IS A REMAINDER, SO ITS ROUND-OFF IS THE TOTAL FLUX'S, as
	// RotatingMirror's card-against-host pair finds: psi_p is 1e-03 at its
	// largest and the total it is taken from is 1.57. Measured, 5.4e-08, which
	// is 0.29 eps of the total and 430 of psi_p's own size; a conductors' load
	// dropped by the law is O( psi_c ). MEASUREMENTS.md M-350.
	mfem::GridFunction total;
	plain->solver->totalPotential( total );
	bound = std::max( bound, meq::tests::scalarAgreementFloor( 0.0, total.Normlinf() ) );
	std::printf( "    the total flux reaches %.3e; the bound in single is %.3e\n",
	             static_cast<double>( total.Normlinf() ), bound );
	std::fflush( stdout );
#endif
	BOOST_TEST( worst < bound, "psi_p differs by " << worst << " on " << scale );
}

/// From the subtracted vacuum, either Delta converges as the meshed route does.
BOOST_AUTO_TEST_CASE( theSubtractedAnisotropicMirrorConvergesFromTheVacuum )
{
	for ( Delta const kind : { Delta::Bump, Delta::Profile, Delta::Field } )
	{
		char const *const label = kindName( kind );
		// Delta( B )'s delta0 is a pressure, and -0.3 is past its fold.
		double const amplitude = kind == Delta::Field ? -0.2 : -0.3;
		auto const r = run( 2, 0, kind, amplitude, false, true );
		BOOST_TEST_REQUIRE( r->converged, label << ": " << r->failure );
		std::printf( "  subtracted, %s, delta0 = %+.1f: %d Newton steps, psi_ax %.6e,"
		             " psi_bnd %.6e; residuals", label, amplitude, r->solver->newtonIterations(),
		             r->solver->psiAxis(), r->solver->psiBoundary() );
		for ( double v : r->solver->newtonResiduals() )
			std::printf( " %.2e", v );
		std::printf( "\n" );
		std::fflush( stdout );
		BOOST_TEST( r->solver->newtonIterations() <= 8, label << ": "
		            << r->solver->newtonIterations() << " Newton steps" );
	}
}

/// THE SUBTRACTED MIRROR IS THE MESHED ONE: on each level the RMS of the two
/// routes' total psi* over the machine region, for the isotropic mirror and
/// for both Deltas at delta0 = -0.3. The isotropic difference is what the two
/// discretisations leave; an anisotropic one must fall as it does and stay
/// within a factor of it. A conductors' load dropped, or of the wrong sign,
/// leaves an O( delta0 ) difference that does not fall.
///
/// And Delta( B ) at delta0 = -0.2, which reads B from the TOTAL flux, the
/// remainder's u_p less q_c: one read from u_p alone sees the remainder's
/// field, not the machine's, and the two routes then solve different
/// equations.
BOOST_AUTO_TEST_CASE( theSubtractedAnisotropicMirrorIsTheMeshedOne )
{
	Cloud const cloud = machineCloud();
	struct Variant
	{
		char const *label;
		Delta kind;
		double amplitude;
	};
	std::vector<std::vector<double>> differences;
	for ( Variant const &variant : { Variant{ "isotropic", Delta::None, 0.0 },
	                                 Variant{ "bump, delta0 = -0.3", Delta::Bump, -0.3 },
	                                 Variant{ "Delta( Psi ), delta0 = -0.3", Delta::Profile, -0.3 },
	                                 Variant{ "Delta( B ), delta0 = -0.2", Delta::Field, -0.2 } } )
	{
		std::vector<double> d;
		for ( int level : { 0, 1, 2 } )
		{
			auto meshed = run( 2, level, variant.kind, variant.amplitude, false, false );
			auto subtracted = run( 2, level, variant.kind, variant.amplitude, false, true );
			BOOST_TEST_REQUIRE( meshed->converged, "meshed, level " << level << ": "
			                    << meshed->failure );
			BOOST_TEST_REQUIRE( subtracted->converged, "subtracted, level " << level << ": "
			                    << subtracted->failure );
			d.push_back( rms( totalPsiStarAt( *meshed, cloud ),
			                  totalPsiStarAt( *subtracted, cloud ) ) );
		}
		std::printf( "  %s, k = 2: subtracted - meshed psi*, level 0 %.3e, 1 %.3e, 2 %.3e;"
		             " rates %.3f, %.3f\n", variant.label, d[ 0 ], d[ 1 ], d[ 2 ],
		             std::log2( d[ 0 ]/d[ 1 ] ), std::log2( d[ 1 ]/d[ 2 ] ) );
		std::fflush( stdout );
		differences.push_back( d );
	}

	std::vector<double> const &isotropic = differences[ 0 ];
	for ( std::size_t v = 1; v < differences.size(); ++v )
		for ( std::size_t i = 0; i < isotropic.size(); ++i )
		{
			double const ratio = differences[ v ][ i ]/isotropic[ i ];
			BOOST_TEST( ( ratio > 0.5 && ratio < 2.0 ), "variant " << v << ", level " << i
			            << ": subtracted - meshed is " << ratio << " of the isotropic one" );
			if ( i > 0 )
				BOOST_TEST( differences[ v ][ i ] < differences[ v ][ i - 1 ] );
		}
}


