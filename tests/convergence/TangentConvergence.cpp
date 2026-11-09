/*
 * THE DISCRETE TANGENT CONVERGES TO THE EXACT ONE: delta psi_h = d psi_h/d lambda,
 * the derivative of the DISCRETE solution in a parameter its source reads,
 * against the closed-form d psi/d lambda, as h -> 0.
 *
 * Discretise, then differentiate. delta x = -J^{-1} ( dR/d lambda ) at the
 * converged state (GradShafranovSolver::solveTangent(), dR/d lambda from
 * assembleSourceDerivative()), J the Newton matrix. For a semilinear source
 * that is an HDG discretisation of the linearised equation
 *
 *     ( -Delta* - dF/dpsi ) delta psi = ( dF/d lambda ) delta lambda,
 *
 * with its coefficients read at psi_h rather than psi, and
 * Brezzi-Rappaz-Raviart would have it converge at psi's own rates -- if the
 * linearised operator is invertible, dF/dpsi is Lipschitz and Newton is
 * solved to round-off. This file MEASURES it, MEASUREMENTS.md M-373 and
 * M-374:
 *
 *   * the tangent against a Richardson difference of two converged solves;
 *   * the closed form: Kaltsas & Throumoulopoulos's exponential similarity
 *     solution, F = f0 e^{n psi}( 1 + eps R^2 ), whose psi( R, z; n ) is exact
 *     in the exponent n -- which F reads inside its nonlinearity, so
 *     d psi/d n is neither linear in the solution nor a constant -- over
 *     k = 1...4 and n = 4...32 cells a side, in L2 AND L-infinity, in psi, in
 *     q and in the post-processed psi*, with psi's own errors beside them;
 *     and d psi/d f0 = -1/( n f0 ), a constant the discrete tangent must
 *     reproduce to round-off, as the sanity rung;
 *   * where the conditions fail: a linearised operator near-singular, the
 *     error's constant against 1/sigma. (Where dF/dpsi jumps the tangent
 *     loses its order; a table's end is the case, M-360.)
 *
 * L-infinity here is the maximum over the points of a 2k + 6 rule on every
 * element, component-wise for q -- GridFunction::ComputeMaxError().
 */
#define BOOST_TEST_MODULE TangentConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "analytic/SimilarityExponential.hpp"
#include "convergence/ConvergenceHarness.hpp"

namespace
{
	using meq::tests::Rectangle;
	using meq::tests::standardBox;

	/*
	 * THE EXPONENTIAL FAMILY IN ITS EXPONENT n, written out because the
	 * fixture keeps its constants private. psi = ( log K + log sech^2 s )/n,
	 * K = b^2 c n/( 2 f0 ), s = m zeta, m = n sqrt( c )/4,
	 * zeta = a R^2 + 2 b z + c~, so with ds/dn = s/n
	 *
	 *     d psi/d n      = ( -psi + 1/n - 2 s tanh( s )/n )/n
	 *     grad d psi/d n = -( sqrt( c )/2 ) sech^2( s )( s/n ) grad zeta
	 *
	 * since grad psi = -( sqrt( c )/2 ) tanh( s ) grad zeta. Checked against
	 * differences of the fixture in n below.
	 */
	struct ExponentialFamily
	{
		double n = 3.0, f0 = 0.5, a = 1.0, b = 1.0, c = 1.0, cTilde = 0.0;

		meq::analytic::SimilarityExponential at( double exponent, double amplitude ) const
		{
			return meq::analytic::SimilarityExponential( exponent, amplitude, a, b, c, cTilde );
		}

		meq::analytic::SimilarityExponential fixture() const
		{
			return at( n, f0 );
		}

		double s( double radius, double z ) const
		{
			return n*std::sqrt( c )/4.0*( a*radius*radius + 2.0*b*z + cTilde );
		}

		double dPsiDn( double radius, double z ) const
		{
			double const arg = s( radius, z );
			return ( -fixture().psi( radius, z ) + 1.0/n - 2.0*arg*std::tanh( arg )/n )/n;
		}

		void dFluxDn( double radius, double z, double &qR, double &qZ ) const
		{
			double const arg = s( radius, z );
			double const sech = 1.0/std::cosh( arg );
			double const scale = -( std::sqrt( c )/2.0 )*sech*sech*( arg/n );
			qR = scale*2.0*a*radius/radius;
			qZ = scale*2.0*b/radius;
		}

		double dPsiDf0() const
		{
			return -1.0/( n*f0 );
		}
	};

	/// F( R, z, psi ) of the family, as the solver's Source; dF/dn = psi F
	/// and dF/d f0 = F/f0 at the total flux.
	using Source = meq::tests::EquilibriumSource<meq::analytic::SimilarityExponential>;

	std::vector<mfem::IntegrationRule const *> rules( int order )
	{
		std::vector<mfem::IntegrationRule const *> out( mfem::Geometry::NumGeom );
		for ( int i = 0; i < mfem::Geometry::NumGeom; ++i )
			out[ i ] = &mfem::IntRules.Get( i, 2*order + 6 );
		return out;
	}

	/// L2 and L-infinity of a scalar and a vector field's error.
	struct Errors
	{
		double l2 = 0.0, max = 0.0;
	};

	Errors scalarErrors( mfem::GridFunction const &field, mfem::Coefficient &exact, int order )
	{
		std::vector<mfem::IntegrationRule const *> irs = rules( order );
		return { field.ComputeL2Error( exact, irs.data() ), field.ComputeMaxError( exact, irs.data() ) };
	}

	Errors vectorErrors( mfem::GridFunction const &field, mfem::VectorCoefficient &exact, int order )
	{
		std::vector<mfem::IntegrationRule const *> irs = rules( order );
		return { field.ComputeL2Error( exact, irs.data() ), field.ComputeMaxError( exact, irs.data() ) };
	}

	/// One converged solve of the family and everything read off it.
	struct Solved
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::analytic::SimilarityExponential> eq;
		std::unique_ptr<Source> source;
		std::unique_ptr<mfem::FunctionCoefficient> datum;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		double h = 0.0;
	};

	std::unique_ptr<Solved> solveFamily( ExponentialFamily const &family, int order, int n )
	{
		auto out = std::make_unique<Solved>();
		Rectangle const box = standardBox();
		out->mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( box, n ) );
		out->h = box.width()/n;
		out->eq = std::make_unique<meq::analytic::SimilarityExponential>( family.fixture() );
		out->source = std::make_unique<Source>( *out->eq );
		meq::analytic::SimilarityExponential const *eq = out->eq.get();
		out->datum = std::make_unique<mfem::FunctionCoefficient>(
			[ eq ]( mfem::Vector const &x ) { return eq->psi( x( 0 ), x( 1 ) ); } );
		out->solver = std::make_unique<meq::GradShafranovSolver>( *out->mesh, order );
		out->solver->setSource( *out->source );
		out->solver->setBoundaryData( *out->datum );
		out->solver->setNewtonControl( 1.0e-14, 1.0e-15, 40 );
		out->solver->solve();
		out->solver->postProcess();
		return out;
	}

	/// The tangent of a converged solve in a parameter whose dF/dp is
	/// @a dFdp and whose datum's derivative is @a datum.
	mfem::Vector tangent( Solved &s, std::function<double( double, double, double )> const &dFdp,
	                      mfem::Coefficient *datum,
	                      int passes = meq::GradShafranovSolver::linearisedRefinementPasses )
	{
		mfem::Vector dRdp, dx;
		s.solver->assembleSourceDerivative( s.solver->state(), dFdp, dRdp );
		s.solver->solveTangent( s.solver->state(), dRdp, datum, dx, passes );
		return dx;
	}

	/// The potential and flux blocks of a state as fields -- the flux with
	/// the solver's sign, q = -( the DarcyForm block ) -- and its psi*.
	struct Fields
	{
		mfem::GridFunction psi, q, psiStar;
	};

	Fields fieldsOf( Solved &s, mfem::Vector const &state )
	{
		mfem::FiniteElementSpace &fluxSpace = s.solver->fluxSpace();
		mfem::FiniteElementSpace &potentialSpace = s.solver->potentialSpace();
		Fields out;
		out.q.SetSpace( &fluxSpace );
		out.psi.SetSpace( &potentialSpace );
		int const fluxSize = fluxSpace.GetVSize();
		for ( int i = 0; i < fluxSize; ++i )
			out.q( i ) = -state( i );
		for ( int i = 0; i < potentialSpace.GetVSize(); ++i )
			out.psi( i ) = state( fluxSize + i );
		s.solver->reconstructPotential( state, out.psiStar );
		return out;
	}

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}
}

/*
 * THE FIXTURE'S TANGENT IS ITS OWN DERIVATIVE: d psi/d n and its gradient
 * against central differences of the closed form in n, and the linearised
 * equation -Delta* dpsi - n F dpsi = psi F by differences of Delta*.
 */
BOOST_AUTO_TEST_CASE( theClosedFormTangentIsTheDerivativeOfTheClosedForm )
{
	ExponentialFamily const family;
	double const h = 1.0e-5;
	double worst = 0.0, worstOperator = 0.0;
	for ( double radius = 0.65; radius < 1.4; radius += 0.15 )
		for ( double z = -0.55; z < 0.6; z += 0.2 )
		{
			auto shifted = [ & ]( double by ) { return family.at( family.n + by, family.f0 ); };
			double const difference = ( shifted( h ).psi( radius, z ) - shifted( -h ).psi( radius, z ) )/( 2.0*h );
			double gR1, gZ1, gR2, gZ2, qR, qZ;
			shifted( h ).flux( radius, z, gR1, gZ1 );
			shifted( -h ).flux( radius, z, gR2, gZ2 );
			family.dFluxDn( radius, z, qR, qZ );
			worst = std::max( { worst, std::abs( difference - family.dPsiDn( radius, z ) ),
			                    std::abs( ( gR1 - gR2 )/( 2.0*h ) - qR ),
			                    std::abs( ( gZ1 - gZ2 )/( 2.0*h ) - qZ ) } );

			// -Delta* dpsi = n F dpsi + psi F, Delta* by differences of dPsiDn.
			double const d = 1.0e-4;
			auto u = [ & ]( double rr, double zz ) { return family.dPsiDn( rr, zz ); };
			double const inner = ( ( u( radius + d, z ) - u( radius, z ) )/d/( radius + d/2.0 )
			                       - ( u( radius, z ) - u( radius - d, z ) )/d/( radius - d/2.0 ) )/d;
			double const deltaStar = radius*inner
			                         + ( u( radius, z + d ) - 2.0*u( radius, z ) + u( radius, z - d ) )/( d*d );
			meq::analytic::SimilarityExponential const eq = family.fixture();
			double const psi = eq.psi( radius, z );
			double const f = eq.f( radius, z, psi );
			worstOperator = std::max( worstOperator,
				std::abs( -deltaStar - family.n*f*u( radius, z ) - psi*f )/( std::abs( psi*f ) + 1.0 ) );
		}
	std::printf( "  d psi/d n against differences of the closed form: %.2e; the linearised "
	             "equation's residual %.2e\n", worst, worstOperator );
	BOOST_TEST( worst < 1.0e-8 );
	BOOST_TEST( worstOperator < 1.0e-5 );
}

/*
 * THE TANGENT IS THE DERIVATIVE OF THE DISCRETE SOLUTION: against
 * Richardson's ( 4 D( h/2 ) - D( h ) )/3 of two converged solves at n +- h on
 * the same mesh, the whole state compared. Newton is closed to round-off, so
 * the difference's noise is 1e-15/h.
 */
BOOST_AUTO_TEST_CASE( theTangentIsADifferenceOfTwoConvergedSolves )
{
	int const order = 2, cells = 8;
	ExponentialFamily const family;
	auto base = solveFamily( family, order, cells );
	mfem::FunctionCoefficient datum( [ &family ]( mfem::Vector const &x )
	                                 { return family.dPsiDn( x( 0 ), x( 1 ) ); } );
	mfem::Vector const dx = tangent( *base,
		[]( double, double, double ) { return 0.0; }, &datum );

	// The source's own dF/dn, at the total flux, through the solver's.
	meq::analytic::SimilarityExponential const eq = family.fixture();
	mfem::Vector const exact = tangent( *base,
		[ &eq ]( double radius, double z, double psi ) { return psi*eq.f( radius, z, psi ); }, &datum );

	auto solvedAt = [ & ]( double by )
	{
		ExponentialFamily moved = family;
		moved.n += by;
		auto s = solveFamily( moved, order, cells );
		return mfem::Vector( s->solver->state() );
	};
	double const h = 1.0e-3;
	auto central = [ & ]( double step )
	{
		mfem::Vector d = solvedAt( step );
		d -= solvedAt( -step );
		d /= 2.0*step;
		return d;
	};
	mfem::Vector richardson = central( h/2.0 );
	richardson *= 4.0;
	richardson -= central( h );
	richardson /= 3.0;

	mfem::Vector error( richardson );
	error -= exact;
	double const relative = error.Norml2()/exact.Norml2();
	mfem::Vector datumOnly( dx );
	datumOnly -= exact;
	std::printf( "  the tangent against a Richardson difference of two converged solves: %.3e "
	             "relative (the datum's share alone misses by %.3e)\n", relative,
	             datumOnly.Norml2()/exact.Norml2() );
	// A precondition of the format, not a tolerance: in single the
	// reference, a difference of two single solves, is itself good to only
	// about 2e-2 over this step (M-378), so it cannot judge a tangent to
	// 1e-8. Single checks the tangent against the closed form instead.
#ifdef MFEM_USE_SINGLE
	BOOST_TEST_MESSAGE( "  single: a difference of two single solves cannot resolve the tangent to "
	                    "1e-8; not asserted (measured " << relative << ")" );
#else
	BOOST_TEST( relative < 1.0e-8, "the tangent misses the difference of two converged solves by "
	            << relative << " relative" );
#endif
}

/*
 * THE SANITY RUNG: d psi/d f0 = -1/( n f0 ), a constant, and the discrete
 * tangent equation reproduces a constant exactly -- the diffusion of a
 * constant is zero in HDG, and -( n F_h ) c = F_h/f0 pointwise at any psi_h.
 * So the error is round-off at every mesh and degree, whatever psi_h's is.
 */
BOOST_AUTO_TEST_CASE( theTangentInTheAmplitudeIsTheConstantItShouldBe )
{
	ExponentialFamily const family;
	mfem::ConstantCoefficient constant( family.dPsiDf0() );
	meq::analytic::SimilarityExponential const eq = family.fixture();
	for ( int order : { 1, 3 } )
		for ( int cells : { 4, 16 } )
		{
			auto s = solveFamily( family, order, cells );
			mfem::Vector const dx = tangent( *s,
				[ &eq, &family ]( double radius, double z, double psi )
				{ return eq.f( radius, z, psi )/family.f0; }, &constant );
			Fields const f = fieldsOf( *s, dx );
			Errors const e = scalarErrors( f.psi, constant, order );
			std::printf( "  d psi/d f0, k = %d, n = %2d: L2 %.2e, max %.2e against |c| = %.3f\n",
			             order, cells, e.l2, e.max, std::abs( family.dPsiDf0() ) );
			// Round-off at every mesh: 1e-11 in double, and max( 1e-11,
			// 20 eps |c| ) in single, the approved single floor (M-378).
			double const bound = std::max( 1.0e-11, 20.0*static_cast<double>(
				std::numeric_limits<mfem::real_t>::epsilon() )*std::abs( family.dPsiDf0() ) );
			BOOST_TEST( e.max < bound, "the constant tangent misses by " << e.max );
		}
}

/*
 * THE RATES. psi, q and psi* of the solution and of its tangent in n, in L2
 * and L-infinity, side by side, k = 1...4 on 4, 8, 16, 32 cells a side.
 * Asserted as the solution's are (ConvergenceHarness.hpp's checkOrder()):
 * every pair's rate at least k + 1 - 0.2 in psi and q, and k + 2 - 0.2 in
 * psi*, for the tangent in L2 -- and in L-infinity, which is what the
 * question was, from the second pair on. M-373.
 */
namespace
{
	struct Row
	{
		double h;
		Errors psi, q, star, dpsi, dq, dstar;
		/// max | d psi/d n | on the same points: what 1e-6 is of.
		double scale = 0.0;
	};

	Row measureRow( ExponentialFamily const &family, int order, int cells )
	{
		auto s = solveFamily( family, order, cells );
		meq::analytic::SimilarityExponential const eq = family.fixture();
		mfem::FunctionCoefficient dDatum( [ &family ]( mfem::Vector const &x )
		                                  { return family.dPsiDn( x( 0 ), x( 1 ) ); } );
		mfem::Vector const dx = tangent( *s,
			[ &eq ]( double radius, double z, double psi ) { return psi*eq.f( radius, z, psi ); },
			&dDatum );

		mfem::FunctionCoefficient psiExact( [ &eq ]( mfem::Vector const &x ) { return eq.psi( x( 0 ), x( 1 ) ); } );
		mfem::VectorFunctionCoefficient qExact( 2, [ &eq ]( mfem::Vector const &x, mfem::Vector &v )
		{
			double r = 0.0, z = 0.0;
			eq.flux( x( 0 ), x( 1 ), r, z );
			v( 0 ) = r;
			v( 1 ) = z;
		} );
		mfem::VectorFunctionCoefficient dqExact( 2, [ &family ]( mfem::Vector const &x, mfem::Vector &v )
		{
			double r = 0.0, z = 0.0;
			family.dFluxDn( x( 0 ), x( 1 ), r, z );
			v( 0 ) = r;
			v( 1 ) = z;
		} );

		Row row;
		row.h = s->h;
		Fields const solution = fieldsOf( *s, s->solver->state() );
		Fields const tangentFields = fieldsOf( *s, dx );
		row.psi = scalarErrors( solution.psi, psiExact, order );
		row.q = vectorErrors( solution.q, qExact, order );
		row.star = scalarErrors( solution.psiStar, psiExact, order + 1 );
		row.dpsi = scalarErrors( tangentFields.psi, dDatum, order );
		row.dq = vectorErrors( tangentFields.q, dqExact, order );
		row.dstar = scalarErrors( tangentFields.psiStar, dDatum, order + 1 );
		mfem::GridFunction zero( tangentFields.psi );
		zero = 0.0;
		row.scale = scalarErrors( zero, dDatum, order ).max;
		return row;
	}

	/*
	 * SIZED TO 1e-6, NOT TO ROUND-OFF. MEQ is accurate to about 1e-6
	 * globally and a suite that hyper-resolves slows everyone, so each
	 * degree's ladder stops near 1e-6 of max | d psi/d n | in L-infinity, or
	 * at 32 cells where the degree is too low to get there: k = 4 stops at 16
	 * cells (8.7e-08 absolute, 1.1e-07 relative), the rest at 32.
	 */
	std::vector<int> ladder( int order )
	{
		return order == 4 ? std::vector<int>{ 4, 8, 16 } : std::vector<int>{ 4, 8, 16, 32 };
	}

	std::vector<Row> rateTable( int order, std::vector<int> const &meshes, bool assert )
	{
		ExponentialFamily const family;
		std::vector<Row> rows;
		for ( int cells : meshes )
			rows.push_back( measureRow( family, order, cells ) );

		std::printf( "\n  k = %d: the solution and its tangent in n, errors (rates)\n", order );
		std::printf( "    %6s | %-21s %-21s | %-21s %-21s | %-21s %-21s\n", "h",
		             "psi L2", "psi max", "q L2", "q max", "psi* L2", "psi* max" );
		auto cell = [ & ]( std::size_t i, double Errors::*field, Errors Row::*which )
		{
			double const e = rows[ i ].*which.*field;
			if ( i == 0 )
				std::printf( " %9.3e (    ) ", e );
			else
				std::printf( " %9.3e (%4.2f) ", e,
				             rate( rows[ i - 1 ].*which.*field, e, rows[ i - 1 ].h/rows[ i ].h ) );
		};
		for ( char const *label : { "psi_h", "dpsi_h" } )
		{
			bool const tangentRow = label[ 0 ] == 'd';
			std::printf( "    %s\n", label );
			for ( std::size_t i = 0; i < rows.size(); ++i )
			{
				std::printf( "    %6.4f |", rows[ i ].h );
				cell( i, &Errors::l2, tangentRow ? &Row::dpsi : &Row::psi );
				cell( i, &Errors::max, tangentRow ? &Row::dpsi : &Row::psi );
				std::printf( "|" );
				cell( i, &Errors::l2, tangentRow ? &Row::dq : &Row::q );
				cell( i, &Errors::max, tangentRow ? &Row::dq : &Row::q );
				std::printf( "|" );
				cell( i, &Errors::l2, tangentRow ? &Row::dstar : &Row::star );
				cell( i, &Errors::max, tangentRow ? &Row::dstar : &Row::star );
				std::printf( "\n" );
			}
		}
		std::printf( "    the tangent's max error relative to max | d psi/d n | = %.4f:", rows.front().scale );
		for ( Row const &row : rows )
			std::printf( " %.2e", row.dpsi.max/row.scale );
		std::printf( "\n" );
		std::fflush( stdout );
		if ( !assert )
			return rows;

		double const wanted = order + 1.0 - meq::tests::rateSlack;
		for ( std::size_t i = 1; i < rows.size(); ++i )
		{
			double const ratio = rows[ i - 1 ].h/rows[ i ].h;
			struct Check
			{
				char const *name;
				Errors Row::*which;
				double Errors::*norm;
				double wanted;
				/// The solution's own errors in the same field.
				Errors Row::*solution;
			};
			/*
			 * L-INFINITY FROM THE SECOND PAIR ON. The coarsest pair, 4 to 8
			 * cells, is pre-asymptotic in the maximum for the SOLUTION as much
			 * as for its tangent -- psi's own max reads 3.81 at k = 3 and 4.66
			 * at k = 4 there (M-373) -- so it is skipped for both, by
			 * name, and every finer pair carries the full bound.
			 */
			bool const maxAsserted = i >= 2;
			for ( Check const &check : {
				      Check{ "d psi L2", &Row::dpsi, &Errors::l2, wanted, &Row::psi },
				      Check{ "d psi max", &Row::dpsi, &Errors::max, wanted, &Row::psi },
				      Check{ "d q L2", &Row::dq, &Errors::l2, wanted, &Row::q },
				      Check{ "d q max", &Row::dq, &Errors::max, wanted, &Row::q },
				      Check{ "d psi* L2", &Row::dstar, &Errors::l2, wanted + 1.0, &Row::star },
				      Check{ "d psi* max", &Row::dstar, &Errors::max, wanted + 1.0, &Row::star } } )
			{
				if ( check.norm == &Errors::max && !maxAsserted )
					continue;
				/*
				 * THE TANGENT IS HELD TO ITS SOLUTION'S RATE WHERE THE SOLUTION
				 * HAS ONE. A pair on which psi_h's own error in the same field
				 * and norm misses k + 1 - 0.2 is past the format's floor for
				 * the SOLUTION -- in single, q and psi* stop falling from 16
				 * cells -- and the tangent, refined against the same residual,
				 * can be no better. The precondition is the solution's, and it
				 * holds on every pair in double (M-378).
				 */
				double const solutionRate = rate( rows[ i - 1 ].*check.solution.*check.norm,
				                                  rows[ i ].*check.solution.*check.norm, ratio );
				if ( solutionRate < check.wanted )
				{
					std::printf( "    k = %d, h = %.4f: %s not asserted -- the solution's own rate there is %.2f\n",
					             order, rows[ i ].h, check.name, solutionRate );
					continue;
				}
				double const r = rate( rows[ i - 1 ].*check.which.*check.norm, rows[ i ].*check.which.*check.norm, ratio );
				BOOST_TEST( r >= check.wanted, "k = " << order << ", h = " << rows[ i ].h << ": the tangent's "
				            << check.name << " converged at " << r << ", wanted " << check.wanted );
			}
		}
		return rows;
	}
}

BOOST_AUTO_TEST_CASE( theTangentConvergesAtTheSolutionsRatesAtDegreeOne ) { rateTable( 1, ladder( 1 ), true ); }
BOOST_AUTO_TEST_CASE( theTangentConvergesAtTheSolutionsRatesAtDegreeTwo ) { rateTable( 2, ladder( 2 ), true ); }
BOOST_AUTO_TEST_CASE( theTangentConvergesAtTheSolutionsRatesAtDegreeThree ) { rateTable( 3, ladder( 3 ), true ); }
BOOST_AUTO_TEST_CASE( theTangentConvergesAtTheSolutionsRatesAtDegreeFour ) { rateTable( 4, ladder( 4 ), true ); }

/*
 * MEASUREMENT ONLY, NOT IN THE SUITE (disabled; run it by name): where the
 * low degrees reach 1e-6 in L-infinity, k = 1 and 2 to 128 cells.
 * M-373.
 */
BOOST_AUTO_TEST_CASE( whereTheTangentReachesOneInAMillion, * boost::unit_test::disabled() )
{
	rateTable( 1, { 32, 64, 128 }, false );
	rateTable( 2, { 32, 64, 128 }, false );
}

/*
 * ITERATIVE REFINEMENT, BY PASS: the tangent in n's error against the closed
 * form after 0...4 passes, beside the solution's own, k = 2...4. MEASUREMENT
 * ONLY (disabled; run by name, in each precision). M-378.
 */
BOOST_AUTO_TEST_CASE( theTangentErrorByRefinementPass, * boost::unit_test::disabled() )
{
	ExponentialFamily const family;
	meq::analytic::SimilarityExponential const eq = family.fixture();
	mfem::FunctionCoefficient psiExact( [ &eq ]( mfem::Vector const &x ) { return eq.psi( x( 0 ), x( 1 ) ); } );
	mfem::FunctionCoefficient dDatum( [ &family ]( mfem::Vector const &x )
	                                  { return family.dPsiDn( x( 0 ), x( 1 ) ); } );
	std::printf( "\n  THE TANGENT IN n BY REFINEMENT PASS: d psi max against the closed form\n" );
	for ( int order : { 2, 3, 4 } )
		for ( int cells : { 8, 16, 32 } )
		{
			if ( order == 4 && cells == 32 )
				continue;
			auto s = solveFamily( family, order, cells );
			Fields const solution = fieldsOf( *s, s->solver->state() );
			Errors const e = scalarErrors( solution.psi, psiExact, order );
			std::printf( "    k = %d, %2d cells: psi_h max %.3e | d psi max:", order, cells, e.max );
			for ( int passes = 0; passes <= 4; ++passes )
			{
				mfem::Vector const dx = tangent( *s,
					[ &eq ]( double radius, double z, double psi ) { return psi*eq.f( radius, z, psi ); },
					&dDatum, passes );
				Fields const f = fieldsOf( *s, dx );
				std::printf( " %d: %.3e", passes, scalarErrors( f.psi, dDatum, order ).max );
			}
			std::printf( "\n" );
			std::fflush( stdout );
		}
}

/*
 * NEAR-SINGULARITY: AN AFFINE SOURCE F = mu psi + lambda g, WITH mu WALKING UP
 * TO THE FIRST EIGENVALUE OF -Delta*. M-374.
 *
 * g is manufactured from phi = sin( a ( R - R_0 ) ) sin( b ( z - z_0 ) ), zero
 * on Gamma, as g = -Delta* phi - mu phi, so the exact psi is lambda phi and the
 * exact tangent in lambda is phi -- for EVERY mu. The discrete tangent solves
 * ( A_h - mu ) delta psi_h = g_h, so its error is the consistency error
 * amplified by || ( A_h - mu )^{-1} ||, which grows as mu reaches A_h's
 * least eigenvalue mu_1,h, the pole of the error's projection on phi --
 * located, not assumed, and sigma_min taken as proportional to the relative
 * distance delta from it. (An inverse-power estimate of the norm through the
 * tangent solve was tried and did not converge in twenty steps far from the
 * pole, where the gap is small; it put mu_1,h anywhere from 22.5 to 23.6.)
 */
namespace
{
	struct Affine
	{
		double mu = 0.0;
		double radiusLow = 0.6, zLow = -0.6, width = 0.8, height = 1.2;

		double a() const { return M_PI/width; }
		double b() const { return M_PI/height; }

		double phi( double radius, double z ) const
		{
			return std::sin( a()*( radius - radiusLow ) )*std::sin( b()*( z - zLow ) );
		}

		/// -Delta* phi - mu phi, Delta* = d_RR - ( 1/R ) d_R + d_zz.
		double g( double radius, double z ) const
		{
			double const sr = std::sin( a()*( radius - radiusLow ) ), cr = std::cos( a()*( radius - radiusLow ) );
			double const sz = std::sin( b()*( z - zLow ) );
			return ( a()*a() + b()*b() - mu )*sr*sz + ( a()/radius )*cr*sz;
		}
	};

	class AffineSource : public meq::Source
	{
		public:
			explicit AffineSource( Affine const &problemIn ) : problem( problemIn ) {}
			double f( double radius, double z, double psi ) const override
			{
				return problem.mu*psi + problem.g( radius, z );
			}
			double dFdPsi( double, double, double ) const override
			{
				return problem.mu;
			}
		private:
			Affine problem;
	};
}

namespace
{
	struct AffinePoint
	{
		bool solved = false;
		double l2 = 0.0, max = 0.0;
		/// < delta psi_h - phi, phi >: the error's share along the mode that
		/// goes singular, which changes sign across its eigenvalue.
		double projection = 0.0;
	};

	AffinePoint affineTangent( double mu, int order, int cells )
	{
		Affine problem;
		problem.mu = mu;
		AffineSource source( problem );
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), cells );
		mfem::ConstantCoefficient zero( 0.0 );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( zero );
		solver.setNewtonControl( 1.0e-14, 1.0e-15, 40 );
		AffinePoint out;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			return out;
		}
		mfem::Vector dRdp, dx;
		solver.assembleSourceDerivative( solver.state(),
			[ &problem ]( double radius, double z, double ) { return problem.g( radius, z ); }, dRdp );
		solver.solveTangent( solver.state(), dRdp, nullptr, dx );
		int const offset = solver.fluxSpace().GetVSize();
		mfem::GridFunction dpsi( &solver.potentialSpace() );
		for ( int i = 0; i < dpsi.Size(); ++i )
			dpsi( i ) = dx( offset + i );
		mfem::FunctionCoefficient phi( [ &problem ]( mfem::Vector const &x ) { return problem.phi( x( 0 ), x( 1 ) ); } );
		Errors const e = scalarErrors( dpsi, phi, order );
		mfem::GridFunction error( dpsi );
		mfem::GridFunction exact( &solver.potentialSpace() );
		exact.ProjectCoefficient( phi );
		error -= exact;
		mfem::LinearForm weight( &solver.potentialSpace() );
		weight.AddDomainIntegrator( new mfem::DomainLFIntegrator( phi ) );
		weight.Assemble();
		out.solved = true;
		out.l2 = e.l2;
		out.max = e.max;
		out.projection = weight( error );
		return out;
	}
}

namespace
{
	void nearSingularity( int order, int cells )
	{
		/*
		 * mu_1,h, A_h's least eigenvalue on this mesh, is where the projection's
		 * sign flips with its size growing on both sides: bracketed by a scan
		 * and bisected. Then the tangent at mu_1,h ( 1 - delta ), delta the
		 * relative distance -- sigma_min of the linearised operator is
		 * proportional to it.
		 */
		double lower = 0.0, upper = 0.0;
		double previous = 0.0, previousMu = 0.0;
		for ( double mu = 20.0; mu <= 30.0; mu += 0.05 )
		{
			AffinePoint const p = affineTangent( mu, order, cells );
			if ( !p.solved )
				continue;
			if ( previousMu > 0.0 && previous*p.projection < 0.0
			     && std::abs( p.projection ) > 1.0e-6 && std::abs( previous ) > 1.0e-6 )
			{
				lower = previousMu;
				upper = mu;
				break;
			}
			previous = p.projection;
			previousMu = mu;
		}
		BOOST_REQUIRE_MESSAGE( upper > lower, "no sign change of the tangent's error along the least mode "
		                       "between mu = 20 and 30" );
		double const lowSign = affineTangent( lower, order, cells ).projection;
		for ( int i = 0; i < 50; ++i )
		{
			double const mid = 0.5*( lower + upper );
			AffinePoint const p = affineTangent( mid, order, cells );
			if ( p.solved && p.projection*lowSign > 0.0 )
				lower = mid;
			else
				upper = mid;
		}
		/*
		 * THE SIGN FLIP IS NOT THE POLE. The projection is A/( mu_1,h - mu ) + B,
		 * which changes sign at the pole AND at the zero A/B away from it; the
		 * bisection finds one of the two. So the pole is fitted: three points
		 * 1e-2 to 1e-4 below the flip give P_i mu_1 - C + B mu_i = P_i mu_i,
		 * linear in ( mu_1, C, B ).
		 */
		double const flip = 0.5*( lower + upper );
		double m[ 3 ][ 4 ];
		int row = 0;
		for ( double delta : { 1.0e-2, 1.0e-3, 1.0e-4 } )
		{
			double const mu = flip*( 1.0 - delta );
			double const p = affineTangent( mu, order, cells ).projection;
			m[ row ][ 0 ] = p;
			m[ row ][ 1 ] = -1.0;
			m[ row ][ 2 ] = mu;
			m[ row ][ 3 ] = p*mu;
			++row;
		}
		for ( int col = 0; col < 3; ++col )
			for ( int r = col + 1; r < 3; ++r )
			{
				double const factor = m[ r ][ col ]/m[ col ][ col ];
				for ( int c = col; c < 4; ++c )
					m[ r ][ c ] -= factor*m[ col ][ c ];
			}
		double x[ 3 ];
		for ( int r = 2; r >= 0; --r )
		{
			double sum = m[ r ][ 3 ];
			for ( int c = r + 1; c < 3; ++c )
				sum -= m[ r ][ c ]*x[ c ];
			x[ r ] = sum/m[ r ][ r ];
		}
		double const eigen = x[ 0 ];
		std::printf( "\n  the sign flip at mu = %.12f; the fitted pole at %.12f\n", flip, eigen );

		std::printf( "\n  NEAR-SINGULARITY: k = %d, n = %d, F = mu psi + g, exact tangent phi for every mu\n",
		             order, cells );
		std::printf( "    the linearised operator's least eigenvalue mu_1,h = %.12f\n", eigen );
		std::printf( "    %10s %14s %14s %14s %14s %14s\n", "delta", "mu", "dpsi L2 error", "dpsi max error",
		             "L2 x delta", "amplification" );
		double far = 0.0;
		std::vector<double> errors, deltas;
		for ( double delta : { 1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 } )
		{
			double const mu = eigen*( 1.0 - delta );
			AffinePoint const p = affineTangent( mu, order, cells );
			if ( !p.solved )
			{
				std::printf( "    %10.1e %14.8f did not solve\n", delta, mu );
				continue;
			}
			errors.push_back( p.l2 );
			deltas.push_back( delta );
			if ( far == 0.0 )
				far = p.max;
			std::printf( "    %10.1e %14.8f %14.4e %14.4e %14.4e %14.2f\n", delta, mu, p.l2, p.max, p.l2*delta,
			             p.max/far );
		}
		std::fflush( stdout );
		// What follows from 1/sigma, and is asserted: close in, the error grows as
		// 1/delta -- the last three decades within a factor of three of it.
		BOOST_REQUIRE( errors.size() >= 6u );
		std::size_t const n = errors.size();
		double const growth = std::log10( errors[ n - 1 ]/errors[ n - 4 ] )/std::log10( deltas[ n - 4 ]/deltas[ n - 1 ] );
		std::printf( "    the error's growth over the last three decades of delta: delta^-%.3f\n", growth );
		BOOST_TEST( growth > 0.8 );
	}
}

/*
 * MEASUREMENT ONLY, NOT IN THE SUITE (disabled; run it by name): a scan of
 * two hundred solves and a bisection per mesh. Two meshes, so the residue's
 * fall with h is measured as well as the amplification.
 */
BOOST_AUTO_TEST_CASE( theTangentsErrorConstantGrowsAsTheOperatorNearsSingularity,
                      * boost::unit_test::disabled() )
{
	nearSingularity( 2, 8 );
	nearSingularity( 2, 16 );
}
