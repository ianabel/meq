#define BOOST_TEST_MODULE QRefiningStudy
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/ContinuousField.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/FluxExtraction.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "SourceTree.hpp"

/*
 * PHASE-1 STUDY FOR THE q-DRIVEN SOLVE'S REFINING REPRESENTATION (branch
 * q-refining-fit, MEASUREMENTS.md M-QR-1 and up). NOT A TEST: every case is
 * disabled and run by name. It prototypes, on M-361's case, the design
 * TODO "The q-driven solve is inconsistent" asks for:
 *
 *   g^2( Psi ) carried as its VALUES at the family's surfaces plus one at the
 *   plasma edge, interpolated by a not-a-knot cubic spline in the source's Psi;
 *   q matched at every surface; the edge value held by R0B0; and the square
 *   system solved by Newton on the exact tangent -- one bordered back-solve
 *   per unknown at a retained linearisation, and the extraction's response
 *   to each tangent field differenced along it (no equilibrium re-solved).
 *
 * The target is the q of a forward solve with the closed form
 * g^2 = 2.2^2 + ( e^( 2 Psi ) - 1 )/( e^2 - 1 ) on a reference mesh, extracted
 * at whatever labels a case asks for, so no table interpolates it.
 *
 * MODES, by environment (the study's knobs, not a design):
 *   default        the square system: an edge node held by R0B0 (M-QR-2: ill-
 *                  conditioned through the edge layer)
 *   QR_LS          no edge node; the layer is the end cubic continued, R0B0 a
 *                  row, Gauss-Newton on the m + 1 rows (M-QR-3)
 *   QR_AXISNODE    a node on the axis, collocated on q from the axis Hessian
 *   QR_INNER/OUTER the cut in Psi_N; QR_COUNTS, QR_DIAGONAL, QR_ORDER,
 *                  QR_MINLEVEL/MAXLEVEL the sweep; QR_START_EXACT the closed
 *                  form as the start; QR_NZ, QR_SYMMETRIC the mesh (by default
 *                  every vertex is moved off z = 0, M-QR-2)
 * The probes warmProbe, fieldProbe and innerProbe are the measurements M-QR-1
 * and M-QR-2 cite.
 */
namespace
{
	double const g0 = 2.2;
	double const amplitude = 1.0;
	double const kappa = 2.0;
	double const fourPiSquared = 4.0*M_PI*M_PI;

	/// QR_AXISNODE: a node at the axis, collocated on q there from the
	/// located axis's Hessian, q_ax = g/( R^2 sqrt( det dq/dx ) ).
	bool axisNode() { return std::getenv( "QR_AXISNODE" ) != nullptr; }
	/// QR_LS: no edge node. The layer outside the last surface is the end
	/// cubic continued, R0B0 is a ROW ( g^2( Psi = 0 ) - g0^2 )/g0^2 beside the
	/// q rows, and the overdetermined system is solved by Gauss-Newton.
	bool leastSquares() { return std::getenv( "QR_LS" ) != nullptr; }
	/// The R0B0 row under QR_LS: ( g^2 continued to the edge - g0^2 )/g0^2.
	/// Set by the caller to the representation in use.
	std::function<double( std::vector<double> const & )> pinRow;
	std::function<std::vector<double>()> pinGradient;

	double innerCut() { return std::getenv( "QR_INNER" ) ? std::atof( std::getenv( "QR_INNER" ) ) : 0.05; }
	double outerCut() { return std::getenv( "QR_OUTER" ) ? std::atof( std::getenv( "QR_OUTER" ) ) : 0.95; }

	double closedGSquared( double psi )
	{
		return g0*g0 + amplitude*( std::exp( kappa*psi ) - 1.0 )/( std::exp( kappa ) - 1.0 );
	}

	class ClosedGGPrime : public meq::Profile
	{
		public:
			double operator()( double psi ) const override
			{
				return 0.5*amplitude*kappa*std::exp( kappa*psi )/( std::exp( kappa ) - 1.0 );
			}
			double prime( double psi ) const override { return kappa*( *this )( psi ); }
			double doublePrime( double psi ) const override { return kappa*kappa*( *this )( psi ); }
	};

	/// A not-a-knot cubic spline through ( x_i, y_i ), linear in y, the end
	/// pieces continued as polynomials.
	class CubicSpline
	{
		public:
			CubicSpline() = default;

			explicit CubicSpline( std::vector<double> nodes ) : x( std::move( nodes ) )
			{
				std::size_t const n = x.size();
				if ( n < 4 )
					throw std::invalid_argument( "CubicSpline: four nodes at least" );
				// A M = D y; store A^{-1} D.
				std::vector<double> a( n*n, 0.0 ), d( n*n, 0.0 );
				auto h = [ this ]( std::size_t i ) { return x[ i + 1 ] - x[ i ]; };
				// not-a-knot at x_1: ( M1 - M0 )/h0 = ( M2 - M1 )/h1
				a[ 0*n + 0 ] = -1.0/h( 0 );
				a[ 0*n + 1 ] = 1.0/h( 0 ) + 1.0/h( 1 );
				a[ 0*n + 2 ] = -1.0/h( 1 );
				for ( std::size_t i = 1; i + 1 < n; ++i )
				{
					a[ i*n + i - 1 ] = h( i - 1 );
					a[ i*n + i ] = 2.0*( h( i - 1 ) + h( i ) );
					a[ i*n + i + 1 ] = h( i );
					d[ i*n + i + 1 ] += 6.0/h( i );
					d[ i*n + i ] += -6.0/h( i ) - 6.0/h( i - 1 );
					d[ i*n + i - 1 ] += 6.0/h( i - 1 );
				}
				std::size_t const l = n - 1;
				a[ l*n + l - 2 ] = -1.0/h( l - 2 );
				a[ l*n + l - 1 ] = 1.0/h( l - 2 ) + 1.0/h( l - 1 );
				a[ l*n + l ] = -1.0/h( l - 1 );
				// Gauss-Jordan with partial pivoting on [ A | D ].
				for ( std::size_t c = 0; c < n; ++c )
				{
					std::size_t p = c;
					for ( std::size_t r = c + 1; r < n; ++r )
						if ( std::abs( a[ r*n + c ] ) > std::abs( a[ p*n + c ] ) )
							p = r;
					for ( std::size_t k = 0; k < n; ++k )
					{
						std::swap( a[ c*n + k ], a[ p*n + k ] );
						std::swap( d[ c*n + k ], d[ p*n + k ] );
					}
					double const pivot = a[ c*n + c ];
					for ( std::size_t k = 0; k < n; ++k )
					{
						a[ c*n + k ] /= pivot;
						d[ c*n + k ] /= pivot;
					}
					for ( std::size_t r = 0; r < n; ++r )
					{
						if ( r == c )
							continue;
						double const f = a[ r*n + c ];
						if ( f == 0.0 )
							continue;
						for ( std::size_t k = 0; k < n; ++k )
						{
							a[ r*n + k ] -= f*a[ c*n + k ];
							d[ r*n + k ] -= f*d[ c*n + k ];
						}
					}
				}
				secondFromValues = d;
			}

			void setValues( std::vector<double> const &values )
			{
				std::size_t const n = x.size();
				y = values;
				m.assign( n, 0.0 );
				for ( std::size_t i = 0; i < n; ++i )
					for ( std::size_t k = 0; k < n; ++k )
						m[ i ] += secondFromValues[ i*n + k ]*y[ k ];
			}

			/// Derivative @a order ( 0 to 3 ) at @a t.
			double operator()( double t, int order = 0 ) const
			{
				std::size_t const n = x.size();
				std::size_t i = 0;
				if ( t >= x[ n - 2 ] )
					i = n - 2;
				else if ( t > x[ 0 ] )
					i = static_cast<std::size_t>( std::upper_bound( x.begin(), x.end(), t ) - x.begin() ) - 1;
				double const h = x[ i + 1 ] - x[ i ];
				double const a = x[ i + 1 ] - t, b = t - x[ i ];
				double const ca = y[ i ]/h - m[ i ]*h/6.0, cb = y[ i + 1 ]/h - m[ i + 1 ]*h/6.0;
				switch ( order )
				{
					case 0: return m[ i ]*a*a*a/( 6.0*h ) + m[ i + 1 ]*b*b*b/( 6.0*h ) + ca*a + cb*b;
					case 1: return -m[ i ]*a*a/( 2.0*h ) + m[ i + 1 ]*b*b/( 2.0*h ) - ca + cb;
					case 2: return m[ i ]*a/h + m[ i + 1 ]*b/h;
					default: return ( m[ i + 1 ] - m[ i ] )/h;
				}
			}

			std::vector<double> const &nodes() const { return x; }

		private:
			std::vector<double> x, y, m, secondFromValues;
	};

	/// gg' = ( g^2 )'/2 from the spline of g^2 in the source's Psi.
	class SplineGGPrime : public meq::Profile
	{
		public:
			explicit SplineGGPrime( CubicSpline s ) : spline( std::move( s ) ) {}
			double operator()( double psi ) const override { return 0.5*spline( psi, 1 ); }
			double prime( double psi ) const override { return 0.5*spline( psi, 2 ); }
			double doublePrime( double psi ) const override { return 0.5*spline( psi, 3 ); }
		private:
			CubicSpline spline;
	};

	std::unique_ptr<mfem::Mesh> boxMesh( int levels )
	{
		auto mesh = std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian2D(
			8, std::getenv( "QR_NZ" ) ? std::atoi( std::getenv( "QR_NZ" ) ) : 8, mfem::Element::TRIANGLE, false, 0.8, 1.2 ) );
		mesh->Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out = in;
			out( 0 ) += 0.6;
			out( 1 ) -= 0.6;
		} );
		for ( int i = 0; i < levels; ++i )
			mesh->UniformRefinement();
		// THE AXIS OFF THE MESH LINE z = 0. The box is up-down symmetric, so
		// its axis sits on z = 0, which is a row of faces at every level; the
		// bordered inner Newton then stalls as the located axis changes
		// element (observed at NZ = 8). Every vertex is moved by
		// 0.03 sin( pi ( z + 0.6 )/1.2 ), which keeps the box and takes every
		// interior row of faces off z = 0. QR_SYMMETRIC keeps the plain mesh.
		if ( std::getenv( "QR_SYMMETRIC" ) == nullptr )
			mesh->Transform( []( mfem::Vector const &in, mfem::Vector &out )
			{
				out = in;
				out( 1 ) += 0.03*std::sin( M_PI*( in( 1 ) + 0.6 )/1.2 );
			} );
		return mesh;
	}

	double seconds( std::chrono::steady_clock::time_point since )
	{
		return std::chrono::duration<double>( std::chrono::steady_clock::now() - since ).count();
	}

	/// One equilibrium machine: a mesh, a solver and its normalised source, one
	/// solver for every solve (warm-started), as the driver's loop has.
	struct Machine
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		std::unique_ptr<mfem::ConstantCoefficient> zero;
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		int order = 2;
		int levels = 0;
		bool solvedOnce = false;
		mfem::Vector lastState;

		Machine( int levelsIn, int orderIn, std::shared_ptr<meq::Profile const> ggPrime )
			: order( orderIn ), levels( levelsIn )
		{
			mesh = boxMesh( levels );
			auto pPrime = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( meq::test::sourcePath( "examples/q-driven-pprime.dat" ) ) );
			source = std::make_unique<meq::NormalisedMHDSource>( pPrime, ggPrime, 0.3, 1.0 );
			solver = std::make_unique<meq::GradShafranovSolver>( *mesh, order );
			zero = std::make_unique<mfem::ConstantCoefficient>( 0.0 );
			guess = std::make_unique<mfem::FunctionCoefficient>( []( mfem::Vector const &x )
			{
				return 0.3*std::sin( M_PI*( x( 0 ) - 0.6 )/0.8 )*std::sin( M_PI*( x( 1 ) + 0.6 )/1.2 );
			} );
			solver->setSource( *source, 0.3 );
			solver->setBoundaryData( *zero );
			solver->setInitialGuess( *guess );
			solver->setNewtonControl( 1.0e-12, 1.0e-14, 60 );
		}

		/// Solve and post-process; false if the solve failed.
		bool solve( std::shared_ptr<meq::Profile const> ggPrime )
		{
			source->setGGPrime( std::move( ggPrime ) );
			if ( solvedOnce && std::getenv( "QR_COLD" ) == nullptr )
				solver->setInitialState( lastState, {} );
			try
			{
				solver->solve();
				solver->postProcess();
				solvedOnce = true;
				lastState = solver->state();
			}
			catch ( std::exception const &error )
			{
				std::printf( "      solve failed: %s\n", error.what() );
				std::printf( "        residuals:" );
				for ( double r : solver->newtonResiduals() )
					std::printf( " %.2e", r );
				std::printf( "  psi_ax %.6e\n", solver->psiAxis() );
				return false;
			}
			return true;
		}

		int fluxSize() { return solver->fluxSpace().GetVSize(); }
		int potentialSize() { return solver->potentialSpace().GetVSize(); }
		int traceSize() { return solver->traceSpace().GetVSize(); }
	};

	/// psi*, q* and psihat of some state, ready for the continuous fields.
	struct Fields
	{
		std::unique_ptr<mfem::GridFunction> psiStar, qStar, trace;
	};

	Fields solvedFields( Machine &machine )
	{
		Fields f;
		f.psiStar = std::make_unique<mfem::GridFunction>( machine.solver->postProcessedPotential() );
		f.qStar = std::make_unique<mfem::GridFunction>( machine.solver->postProcessedFlux() );
		f.trace = std::make_unique<mfem::GridFunction>( machine.solver->trace() );
		return f;
	}

	Fields directionFields( Machine &machine, mfem::Vector const &dx )
	{
		Fields f;
		f.psiStar = std::make_unique<mfem::GridFunction>( machine.solver->postProcessedPotential() );
		f.qStar = std::make_unique<mfem::GridFunction>( machine.solver->postProcessedFlux() );
		f.trace = std::make_unique<mfem::GridFunction>( machine.solver->trace() );
		machine.solver->reconstructFields( dx, *f.psiStar, *f.qStar );
		int const offset = machine.fluxSize() + machine.potentialSize();
		for ( int i = 0; i < f.trace->Size(); ++i )
			( *f.trace )( i ) = dx( offset + i );
		return f;
	}

	/// base + s direction, field by field.
	Fields combine( Fields const &base, Fields const &direction, double s )
	{
		Fields f;
		f.psiStar = std::make_unique<mfem::GridFunction>( *base.psiStar );
		f.qStar = std::make_unique<mfem::GridFunction>( *base.qStar );
		f.trace = std::make_unique<mfem::GridFunction>( *base.trace );
		f.psiStar->Add( s, *direction.psiStar );
		f.qStar->Add( s, *direction.qStar );
		f.trace->Add( s, *direction.trace );
		return f;
	}

	meq::FluxFamilyOptions familyOptions( std::size_t surfaces, double inner = -1.0, double outer = -1.0 )
	{
		if ( inner < 0.0 )
			inner = innerCut();
		if ( outer < 0.0 )
			outer = outerCut();
		meq::FluxFamilyOptions options;
		options.surfaces = surfaces;
		options.innerCut = inner;
		options.outerCut = outer;
		return options;
	}

	meq::FluxSurfaceFamily extract( Fields const &fields, meq::CriticalPoint const &axis,
	                                meq::FluxFamilyOptions const &options )
	{
		meq::ContinuousFields const continuous =
			meq::buildContinuousFields( *fields.psiStar, *fields.qStar, *fields.trace );
		meq::ContourTracer const tracer( *continuous.potential, *continuous.flux );
		return meq::extractFluxSurfaces( tracer, axis, 0.0, options );
	}

	/// The axis the family of the machine's solved state is labelled by.
	meq::CriticalPoint solvedAxis( Machine &machine )
	{
		meq::CriticalPointFinder const finder( *machine.solver );
		meq::ContourTracer const tracer( *machine.solver );
		return meq::labelledAxis( *machine.solver, tracer, finder.findAxis() );
	}

	/// V' < R^-2 >/4 pi^2 per surface: q over g.
	std::vector<double> qOverG( meq::FluxSurfaceFamily const &family )
	{
		std::vector<double> out;
		for ( meq::FluxSurface const &s : family.surfaces )
			out.push_back( s.vPrime*s.inverseRSquared/fourPiSquared );
		return out;
	}

	/// q/g per condition: the surfaces' V' < R^-2 >/4 pi^2 and, with an axis
	/// node, 1/( R_ax^2 sqrt( det dq~/dx ) ) at the axis of q~.
	std::vector<double> geometry( Fields const &fields, meq::CriticalPoint const &axis,
	                              meq::FluxFamilyOptions const &options )
	{
		meq::ContinuousFields const continuous =
			meq::buildContinuousFields( *fields.psiStar, *fields.qStar, *fields.trace );
		meq::ContourTracer const tracer( *continuous.potential, *continuous.flux );
		std::vector<double> out = qOverG( meq::extractFluxSurfaces( tracer, axis, 0.0, options ) );
		if ( axisNode() )
		{
			meq::CriticalPointFinder const finder( *continuous.flux, *continuous.potential );
			meq::CriticalPoint const a = finder.findAxis();
			out.push_back( 1.0/( a.radius*a.radius*std::sqrt( a.determinant ) ) );
		}
		return out;
	}

	/// The reference: closed-form g^2 on a fine mesh. q at any labels.
	struct Reference
	{
		std::unique_ptr<Machine> machine;
		Fields fields;
		meq::CriticalPoint axis;

		Reference( int levels, int order )
		{
			auto closed = std::make_shared<ClosedGGPrime const>();
			machine = std::make_unique<Machine>( levels, order, closed );
			auto const start = std::chrono::steady_clock::now();
			BOOST_TEST_REQUIRE( machine->solve( closed ) );
			fields = solvedFields( *machine );
			axis = solvedAxis( *machine );
			std::printf( "  reference: level %d, k = %d, psi_ax %.10f, %.1f s\n", levels, order,
			             axis.psi, seconds( start ) );
		}

		std::vector<double> q( meq::FluxFamilyOptions const &options, std::vector<double> *labels = nullptr,
		                       bool withAxis = true )
		{
			meq::FluxSurfaceFamily const family = extract( fields, axis, options );
			std::vector<double> out = qOverG( family );
			for ( std::size_t i = 0; i < out.size(); ++i )
			{
				double const psiN = family.surfaces[ i ].normalisedFlux;
				out[ i ] *= std::sqrt( closedGSquared( 1.0 - psiN ) );
				if ( labels )
					labels->push_back( psiN );
			}
			if ( axisNode() && withAxis )
			{
				std::vector<double> const g = geometry( fields, axis, options );
				out.push_back( g.back()*std::sqrt( closedGSquared( 1.0 ) ) );
			}
			return out;
		}
	};

	/*
	 * THE REPRESENTATION. Nodes in the source's Psi, ascending: the edge Psi = 0
	 * (index 0, held by R0B0), then 1 - Psi_N of the family's surfaces,
	 * outermost first. Unknowns are the values at the surfaces, u_1 .. u_m.
	 */
	struct Representation
	{
		std::vector<double> psiNodes;   // source Psi, ascending, edge first
		std::vector<double> labels;     // family Psi_N, ascending (surface i)
		CubicSpline spline;

		explicit Representation( std::vector<double> const &familyLabels )
			: labels( familyLabels )
		{
			if ( !leastSquares() )
				psiNodes.push_back( 0.0 );
			for ( std::size_t i = labels.size(); i-- > 0; )
				psiNodes.push_back( 1.0 - labels[ i ] );
			if ( axisNode() )
				psiNodes.push_back( 1.0 );
			spline = CubicSpline( psiNodes );
		}

		std::size_t surfaces() const { return labels.size(); }
		std::size_t unknowns() const { return labels.size() + ( axisNode() ? 1 : 0 ); }

		/// Node index of surface i.
		std::size_t node( std::size_t i ) const { return labels.size() - i - ( leastSquares() ? 1 : 0 ); }

		/// All node values from the unknowns (surface order) and the edge value.
		std::vector<double> values( std::vector<double> const &u, double edge ) const
		{
			std::vector<double> v( psiNodes.size() );
			if ( !leastSquares() )
				v[ 0 ] = edge;
			for ( std::size_t i = 0; i < surfaces(); ++i )
				v[ node( i ) ] = u[ i ];
			if ( u.size() > surfaces() )
				v.back() = u[ surfaces() ];
			return v;
		}

		std::shared_ptr<meq::Profile const> ggPrime( std::vector<double> const &u, double edge ) const
		{
			CubicSpline s( spline );
			s.setValues( values( u, edge ) );
			return std::make_shared<SplineGGPrime const>( s );
		}

		/// The cardinal function of surface j's node: g^2 with u = e_j, edge 0.
		CubicSpline cardinal( std::size_t j ) const
		{
			std::vector<double> e( unknowns(), 0.0 );
			e[ j ] = 1.0;
			CubicSpline s( spline );
			s.setValues( values( e, 0.0 ) );
			return s;
		}

		CubicSpline of( std::vector<double> const &u, double edge ) const
		{
			CubicSpline s( spline );
			s.setValues( values( u, edge ) );
			return s;
		}
	};

	struct Timing
	{
		double solve = 0.0, linearise = 0.0, tangents = 0.0, reconstruct = 0.0, extract = 0.0;
		int solves = 0, extractions = 0;
	};

	/*
	 * THE EXACT-TANGENT JACOBIAN of r_i = q_i/q*_i - 1 in u at the machine's
	 * solved state: q_i = ( V' < R^-2 > )_i sqrt( u_i )/4 pi^2, so
	 *   dr_i/du_j = [ q_i/( 2 u_i ) delta_ij + sqrt( u_i )/4pi^2 D_j( V'<R^-2> )_i ]/q*_i
	 * with D_j the derivative along the bordered tangent [ dx_j ; dpsi_ax,j ]
	 * = -[ J c ; r d ]^{-1} [ dR/du_j ; 0 ], differenced in the EXTRACTION only,
	 * centrally at a relative step @a step of the field.
	 */
	std::vector<double> tangentJacobian( Machine &machine, Representation const &rep,
	                                     std::vector<double> const &u, std::vector<double> const &target,
	                                     meq::FluxFamilyOptions const &options, double step, Timing &timing,
	                                     std::vector<double> *qAtState = nullptr )
	{
		std::size_t const m = rep.unknowns();
		meq::GradShafranovSolver &solver = *machine.solver;
		mfem::Vector const state( solver.state() );
		double const span = solver.psiAxis() - solver.psiBoundary();

		Fields const base = solvedFields( machine );
		meq::CriticalPoint const axis = solvedAxis( machine );
		auto t0 = std::chrono::steady_clock::now();
		std::vector<double> const baseQG = geometry( base, axis, options );
		timing.extract += seconds( t0 );
		++timing.extractions;
		if ( qAtState )
		{
			qAtState->clear();
			for ( std::size_t i = 0; i < m; ++i )
				qAtState->push_back( baseQG[ i ]*std::sqrt( u[ i ] ) );
		}

		t0 = std::chrono::steady_clock::now();
		solver.linearise( state );
		timing.linearise += seconds( t0 );

		t0 = std::chrono::steady_clock::now();
		std::vector<mfem::Vector> dx( m );
		std::vector<double> dAxis( m );
		for ( std::size_t j = 0; j < m; ++j )
		{
			CubicSpline const phi = rep.cardinal( j );
			mfem::Vector dRdu;
			solver.assembleSourceDerivative( state,
				[ &phi, span, &solver ]( double, double, double psi )
				{
					double const x = ( psi - solver.psiBoundary() )/span;
					return 0.5*phi( x, 1 )/span;
				}, dRdu );
			std::vector<double> ds;
			solver.linearisedTangent( dRdu, {}, nullptr, dx[ j ], ds );
			dAxis[ j ] = ds.front();
		}
		timing.tangents += seconds( t0 );

		double const fieldScale = base.psiStar->Normlinf();
		std::vector<double> jac( m*m, 0.0 );
		for ( std::size_t j = 0; j < m; ++j )
		{
			t0 = std::chrono::steady_clock::now();
			Fields const direction = directionFields( machine, dx[ j ] );
			timing.reconstruct += seconds( t0 );
			double const s = step*fieldScale/std::max( direction.psiStar->Normlinf(), 1e-300 );
			t0 = std::chrono::steady_clock::now();
			std::vector<double> plus, minus;
			for ( double sign : { 1.0, -1.0 } )
			{
				meq::CriticalPoint moved = axis;
				moved.psi += sign*s*dAxis[ j ];
				( sign > 0 ? plus : minus ) = geometry( combine( base, direction, sign*s ), moved, options );
			}
			timing.extract += seconds( t0 );
			timing.extractions += 2;
			for ( std::size_t i = 0; i < m; ++i )
			{
				double const d = ( plus[ i ] - minus[ i ] )/( 2.0*s );
				double value = std::sqrt( u[ i ] )*d;
				if ( i == j )
					value += baseQG[ i ]*std::sqrt( u[ i ] )/( 2.0*u[ i ] );
				jac[ i*m + j ] = value/target[ i ];
			}
		}
		if ( leastSquares() )
		{
			std::vector<double> const row = pinGradient();
			jac.insert( jac.end(), row.begin(), row.end() );
		}
		return jac;
	}

	/// r at the machine's solved state.
	std::vector<double> residual( Machine &machine, std::vector<double> const &u,
	                              std::vector<double> const &target, meq::FluxFamilyOptions const &options,
	                              Timing &timing )
	{
		auto const t0 = std::chrono::steady_clock::now();
		std::vector<double> qg;
		try
		{
			qg = geometry( solvedFields( machine ), solvedAxis( machine ), options );
		}
		catch ( std::exception const &error )
		{
			std::printf( "      extraction refused: %.120s\n", error.what() );
			return {};
		}
		timing.extract += seconds( t0 );
		++timing.extractions;
		std::vector<double> r( u.size() );
		for ( std::size_t i = 0; i < u.size(); ++i )
			r[ i ] = qg[ i ]*std::sqrt( u[ i ] )/target[ i ] - 1.0;
		if ( leastSquares() )
			r.push_back( pinRow( u ) );
		return r;
	}

	double maxAbs( std::vector<double> const &v )
	{
		double out = 0.0;
		for ( double x : v )
			out = std::max( out, std::abs( x ) );
		return out;
	}

	/// Dense solve by Gaussian elimination with partial pivoting.
	std::vector<double> denseSolve( std::vector<double> a, std::vector<double> b )
	{
		std::size_t const n = b.size();
		for ( std::size_t c = 0; c < n; ++c )
		{
			std::size_t p = c;
			for ( std::size_t r = c + 1; r < n; ++r )
				if ( std::abs( a[ r*n + c ] ) > std::abs( a[ p*n + c ] ) )
					p = r;
			for ( std::size_t k = 0; k < n; ++k )
				std::swap( a[ c*n + k ], a[ p*n + k ] );
			std::swap( b[ c ], b[ p ] );
			for ( std::size_t r = c + 1; r < n; ++r )
			{
				double const f = a[ r*n + c ]/a[ c*n + c ];
				for ( std::size_t k = c; k < n; ++k )
					a[ r*n + k ] -= f*a[ c*n + k ];
				b[ r ] -= f*b[ c ];
			}
		}
		std::vector<double> x( n );
		for ( std::size_t i = n; i-- > 0; )
		{
			double s = b[ i ];
			for ( std::size_t k = i + 1; k < n; ++k )
				s -= a[ i*n + k ]*x[ k ];
			x[ i ] = s/a[ i*n + i ];
		}
		return x;
	}

	struct Outcome
	{
		bool converged = false;
		int iterations = 0;
		std::vector<double> u;
		std::vector<double> history;
		std::vector<double> misfit;
		int rejected = 0;
		Timing timing;
	};

	/// Newton with an Armijo backtracking search on || r ||_2^2, the exact
	/// tangent Jacobian each step.
	Outcome newton( Machine &machine, Representation const &rep, std::vector<double> u,
	                std::vector<double> const &target, meq::FluxFamilyOptions const &options,
	                double tolerance, int maxIterations, bool verbose )
	{
		Outcome out;
		auto t0 = std::chrono::steady_clock::now();
		if ( !machine.solve( rep.ggPrime( u, g0*g0 ) ) )
			return out;
		out.timing.solve += seconds( t0 );
		++out.timing.solves;
		std::vector<double> r = residual( machine, u, target, options, out.timing );
		auto norm2 = []( std::vector<double> const &v )
		{
			double s = 0.0;
			for ( double x : v )
				s += x*x;
			return s;
		};
		bool lastRejected = false;
		for ( int it = 0; it < maxIterations; ++it )
		{
			out.history.push_back( maxAbs( r ) );
			if ( verbose )
			{
				meq::CriticalPoint const a = solvedAxis( machine );
				std::printf( "      Newton %2d: max |q/q* - 1| %.3e  axis ( %.6f, %.3e ) element %d psi_ax %.6f, inner its %d\n",
				             it, maxAbs( r ), a.radius, a.z, a.element, a.psi, machine.solver->newtonIterations() );
			}
			if ( maxAbs( r ) < tolerance )
			{
				out.converged = true;
				out.iterations = it;
				out.u = u;
				return out;
			}
			if ( it > 0 && lastRejected )
			{
				// The machine holds the last rejected trial; put the iterate back.
				t0 = std::chrono::steady_clock::now();
				machine.solve( rep.ggPrime( u, g0*g0 ) );
				out.timing.solve += seconds( t0 );
				++out.timing.solves;
			}
			std::vector<double> const jac = tangentJacobian( machine, rep, u, target, options, 1.0e-5, out.timing );
			std::vector<double> minusR( r.size() );
			for ( std::size_t i = 0; i < r.size(); ++i )
				minusR[ i ] = -r[ i ];
			std::vector<double> du;
			if ( leastSquares() )
			{
				// Gauss-Newton: J^T J du = -J^T r.
				std::size_t const n = u.size(), rows = r.size();
				std::vector<double> normal( n*n, 0.0 ), rhs( n, 0.0 );
				for ( std::size_t a = 0; a < n; ++a )
				{
					for ( std::size_t b = 0; b < n; ++b )
						for ( std::size_t k = 0; k < rows; ++k )
							normal[ a*n + b ] += jac[ k*n + a ]*jac[ k*n + b ];
					for ( std::size_t k = 0; k < rows; ++k )
						rhs[ a ] -= jac[ k*n + a ]*r[ k ];
				}
				du = denseSolve( normal, rhs );
				double step = 0.0, size = 0.0;
				for ( std::size_t a = 0; a < n; ++a )
				{
					step = std::max( step, std::abs( du[ a ] ) );
					size = std::max( size, std::abs( u[ a ] ) );
				}
				if ( step < 1.0e-11*size )
				{
					out.converged = true;
					out.iterations = it;
					out.u = u;
					out.misfit = r;
					return out;
				}
			}
			else
				du = denseSolve( jac, minusR );
			if ( verbose )
			{
				std::printf( "        r:" );
				for ( double x : r )
					std::printf( " %+.3e", x );
				std::printf( "\n        u:" );
				for ( double x : u )
					std::printf( " %.4f", x );
				std::printf( "\n        du:" );
				for ( double x : du )
					std::printf( " %+.3e", x );
				std::printf( "\n" );
			}
			double alpha = 1.0;
			double const f0 = norm2( r );
			for ( int back = 0; back < 12; ++back )
			{
				std::vector<double> trial( u );
				bool admissible = true;
				for ( std::size_t i = 0; i < u.size(); ++i )
				{
					trial[ i ] += alpha*du[ i ];
					admissible = admissible && trial[ i ] > 0.0;
				}
				if ( admissible )
				{
					t0 = std::chrono::steady_clock::now();
					admissible = machine.solve( rep.ggPrime( trial, g0*g0 ) );
					out.timing.solve += seconds( t0 );
					++out.timing.solves;
				}
				if ( admissible )
				{
					std::vector<double> const rt = residual( machine, trial, target, options, out.timing );
					if ( !rt.empty() && norm2( rt ) <= ( 1.0 - 2.0e-4*alpha )*f0 )
					{
						u = trial;
						r = rt;
						lastRejected = false;
						break;
					}
				}
				if ( verbose )
					std::printf( "        step %.4f rejected\n", alpha );
				++out.rejected;
				lastRejected = true;
				alpha *= 0.5;
				if ( back == 11 )
				{
					// Under least squares the merit's floor is the misfit, and a
					// step that cannot lower it is the minimum: converged.
					out.converged = leastSquares();
					out.misfit = r;
					out.u = u;
					out.iterations = it;
					return out;
				}
			}
		}
		out.u = u;
		out.iterations = maxIterations;
		return out;
	}

	/*
	 * The measures. ( a ) g^2 against the closed form over the cut and over
	 * [ 0, 1 ]; ( b ) q off the knots: the delivered state's own family at
	 * the midpoint labels, g from the representation, against the reference
	 * at the same labels.
	 */
	struct Miss
	{
		double gSquaredCut = 0.0, gSquaredAll = 0.0, qOffKnots = 0.0, qOnKnots = 0.0;
	};

	Miss measure( Machine &machine, Representation const &rep, std::vector<double> const &u,
	              Reference &reference, std::size_t m )
	{
		Miss miss;
		CubicSpline const s = rep.of( u, g0*g0 );
		double const lo = innerCut(), hi = outerCut();
		for ( int k = 0; k <= 2000; ++k )
		{
			double const psiN = static_cast<double>( k )/2000.0;
			double const e = std::abs( s( 1.0 - psiN )/closedGSquared( 1.0 - psiN ) - 1.0 );
			miss.gSquaredAll = std::max( miss.gSquaredAll, e );
			if ( psiN >= lo && psiN <= hi )
				miss.gSquaredCut = std::max( miss.gSquaredCut, e );
		}
		double const rhoIn = std::sqrt( innerCut() ), rhoOut = std::sqrt( outerCut() );
		double const dRho = ( rhoOut - rhoIn )/static_cast<double>( m - 1 );
		double const a = rhoIn + 0.5*dRho, b = rhoOut - 0.5*dRho;
		meq::FluxFamilyOptions const mid = familyOptions( m - 1, a*a, b*b );
		std::vector<double> labels;
		std::vector<double> const qRef = reference.q( mid, &labels );
		meq::FluxSurfaceFamily const family = extract( solvedFields( machine ), solvedAxis( machine ), mid );
		std::vector<double> const qg = qOverG( family );
		for ( std::size_t i = 0; i < qg.size(); ++i )
		{
			double const q = qg[ i ]*std::sqrt( s( 1.0 - family.surfaces[ i ].normalisedFlux ) );
			miss.qOffKnots = std::max( miss.qOffKnots, std::abs( q/qRef[ i ] - 1.0 ) );
		}
		return miss;
	}
}

/*
 * M-QR-1: THE EXACT-TANGENT JACOBIAN AGAINST RE-SOLVED DIFFERENCES. Level 0,
 * k = 2, m surfaces, at the flat start g^2 = g0^2 and at a point near the
 * answer: every column dq/du_j both ways.
 */
BOOST_AUTO_TEST_CASE( theTangentJacobianIsTheDifferencedOne, * boost::unit_test::disabled() )
{
	std::size_t const m = 8;
	meq::FluxFamilyOptions const options = familyOptions( m );
	Reference reference( 2, 3 );
	std::vector<double> labels;
	std::vector<double> const target = reference.q( options, &labels );
	Representation const rep( labels );

	for ( int level : { 0, 1 } )
	for ( int point = 0; point < 2; ++point )
	{
		std::vector<double> u( m, g0*g0 );
		if ( point == 1 )
			for ( std::size_t i = 0; i < m; ++i )
				u[ i ] = closedGSquared( 1.0 - labels[ i ] );
		if ( axisNode() )
			u.push_back( point == 1 ? closedGSquared( 1.0 ) : g0*g0 );
		std::size_t const n = u.size();
		Machine machine( level, 2, rep.ggPrime( u, g0*g0 ) );
		BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( u, g0*g0 ) ) );
		Timing timing;
		std::vector<double> q0;
		std::printf( "  probe: | state before Jacobian | %.6e, potential gf max %.6e\n",
		             machine.solver->state().Normlinf(), machine.solver->potential().Normlinf() );
		std::vector<double> const exact = tangentJacobian( machine, rep, u, target, options, std::getenv( "QR_S" ) ? std::atof( std::getenv( "QR_S" ) ) : 1.0e-5, timing, &q0 );
		if ( std::getenv( "QR_PROBE" ) )
		{
			{
				mfem::Vector const now( machine.solver->state() );
				std::printf( "  probe: | state after Jacobian | %.6e, potential gf max %.6e\n", now.Normlinf(),
				             machine.solver->potential().Normlinf() );
			}
			std::printf( "  probe: re-solve at the same u after the Jacobian: %s\n",
			             machine.solve( rep.ggPrime( u, g0*g0 ) ) ? "ok" : "FAILED" );
			Machine fresh( level, 2, rep.ggPrime( u, g0*g0 ) );
			std::vector<double> up( u );
			up[ 0 ] *= 1.001;
			std::printf( "  probe: fresh machine at u + 0.1%% e_0: %s\n",
			             fresh.solve( rep.ggPrime( up, g0*g0 ) ) ? "ok" : "FAILED" );
		}

		// Re-solved differences, Richardson of central differences at h and h/2.
		std::vector<double> diff( n*n, 0.0 );
		meq::CriticalPoint const baseAxis = solvedAxis( machine );
		auto qAt = [ & ]( std::vector<double> const &v )
		{
			BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( v, g0*g0 ) ) );
			if ( std::getenv( "QR_FIXAXIS" ) )
			{
				meq::CriticalPoint axis = baseAxis;
				axis.psi = machine.solver->psiAxis();
				std::vector<double> const qg = qOverG( extract( solvedFields( machine ), axis, options ) );
				std::vector<double> r( v.size() );
				for ( std::size_t i = 0; i < v.size(); ++i )
					r[ i ] = qg[ i ]*std::sqrt( v[ i ] )/target[ i ] - 1.0;
				return r;
			}
			std::vector<double> r = residual( machine, v, target, options, timing );
			return r;
		};
		for ( std::size_t j = 0; j < n; ++j )
		{
			std::vector<double> col( n );
			std::vector<double> d1( n ), d2( n );
			for ( int pass = 0; pass < 2; ++pass )
			{
				double const h = ( pass == 0 ? 2.0 : 1.0 )*( std::getenv( "QR_RH" ) ? std::atof( std::getenv( "QR_RH" ) ) : 1.0e-5 )*u[ j ];
				std::vector<double> up( u ), um( u );
				up[ j ] += h;
				um[ j ] -= h;
				std::vector<double> const rp = qAt( up ), rm = qAt( um );
				for ( std::size_t i = 0; i < n; ++i )
					( pass == 0 ? d1 : d2 )[ i ] = ( rp[ i ] - rm[ i ] )/( 2.0*h );
			}
			for ( std::size_t i = 0; i < n; ++i )
				diff[ i*n + j ] = ( 4.0*d2[ i ] - d1[ i ] )/3.0;
		}
		double worst = 0.0, scale = 0.0, worstDiag = 0.0, offDiagScale = 0.0;
		for ( std::size_t i = 0; i < n; ++i )
			for ( std::size_t j = 0; j < n; ++j )
			{
				scale = std::max( scale, std::abs( diff[ i*n + j ] ) );
				worst = std::max( worst, std::abs( exact[ i*n + j ] - diff[ i*n + j ] ) );
				if ( i != j )
					offDiagScale = std::max( offDiagScale, std::abs( diff[ i*n + j ] ) );
				else
					worstDiag = std::max( worstDiag, std::abs( diff[ i*n + j ] ) );
			}
		std::printf( "  level %d, %s: max | J_tangent - J_resolved | %.3e against max | J | %.3e "
		             "(diagonal %.3e, off-diagonal %.3e): relative %.3e\n",
		             level, point == 0 ? "flat start" : "near the answer", worst, scale, worstDiag,
		             offDiagScale, worst/scale );
		std::printf( "    one Jacobian: linearise %.3f s, %zu tangents %.3f s, reconstruct %.3f s, "
		             "extractions %.3f s over %d\n",
		             timing.linearise, m, timing.tangents, timing.reconstruct, timing.extract,
		             timing.extractions );
		std::printf( "    row 0 tangent:  " );
		for ( std::size_t j = 0; j < m; ++j )
			std::printf( " %+.6e", exact[ 0*m + j ] );
		std::printf( "\n    row 0 resolved: " );
		for ( std::size_t j = 0; j < m; ++j )
			std::printf( " %+.6e", diff[ 0*m + j ] );
		std::printf( "\n" );
	}
}

/*
 * M-QR-2: THE MISS UNDER h AND THE KNOT COUNT. Newton on the exact tangent,
 * R0B0 = g0, from the flat start.
 */
BOOST_AUTO_TEST_CASE( warmProbe, * boost::unit_test::disabled() )
{
	std::size_t const m = 8;
	meq::FluxFamilyOptions const options = familyOptions( m );
	std::vector<double> labels = meq::fluxSurfaceLabels( options );
	for ( double &l : labels )
		l = l*l;
	Representation const rep( labels );
	std::vector<double> u( m, g0*g0 );
	for ( double f : { 1.0, 1.001, 1.0, 1.0001, 1.0 } )
	{
		Machine machine( 0, 2, rep.ggPrime( u, g0*g0 ) );
		BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( u, g0*g0 ) ) );
		std::printf( "  base: %d its, psi_ax %.8f\n", machine.solver->newtonIterations(), machine.solver->psiAxis() );
		for ( std::size_t j = 0; j < m; ++j )
		{
			std::vector<double> up( u );
			up[ j ] *= f;
			bool const ok = machine.solve( rep.ggPrime( up, g0*g0 ) );
			std::printf( "  warm, u_%zu x %.4f: %s, %d its\n", j, f, ok ? "ok" : "FAILED", machine.solver->newtonIterations() );
		}
		break;
	}
}

BOOST_AUTO_TEST_CASE( fieldProbe, * boost::unit_test::disabled() )
{
	std::size_t const m = 8;
	meq::FluxFamilyOptions const options = familyOptions( m );
	std::vector<double> labels = meq::fluxSurfaceLabels( options );
	for ( double &l : labels )
		l = l*l;
	Representation const rep( labels );
	std::vector<double> u( m, g0*g0 );
	Machine machine( 0, 2, rep.ggPrime( u, g0*g0 ) );
	BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( u, g0*g0 ) ) );
	meq::GradShafranovSolver &solver = *machine.solver;
	for ( std::size_t j : { 0ul, 1ul, 4ul } )
	{
		BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( u, g0*g0 ) ) );
		mfem::Vector const state( solver.state() );
		double const span = solver.psiAxis();
		solver.linearise( state );
		CubicSpline const phi = rep.cardinal( j );
		mfem::Vector dRdu, dx;
		solver.assembleSourceDerivative( state, [ &phi, span ]( double, double, double psi )
			{ return 0.5*phi( psi/span, 1 )/span; }, dRdu );
		std::vector<double> ds;
		solver.linearisedTangent( dRdu, {}, nullptr, dx, ds );
		Fields const lin = directionFields( machine, dx );
		double const h = ( std::getenv( "QR_H" ) ? std::atof( std::getenv( "QR_H" ) ) : 1.0e-4 )*u[ j ];
		std::vector<mfem::Vector> states;
		std::vector<Fields> fields;
		std::vector<double> axes;
		for ( double sign : { 1.0, -1.0 } )
		{
			std::vector<double> v( u );
			v[ j ] += sign*h;
			BOOST_TEST_REQUIRE( machine.solve( rep.ggPrime( v, g0*g0 ) ) );
			states.emplace_back( solver.state() );
			fields.push_back( solvedFields( machine ) );
			axes.push_back( solver.psiAxis() );
		}
		auto rel = [ h ]( mfem::Vector const &a, mfem::Vector const &b, mfem::Vector const &d )
		{
			mfem::Vector diff( a );
			diff -= b;
			diff *= 1.0/( 2.0*h );
			diff -= d;
			return diff.Normlinf()/d.Normlinf();
		};
		std::printf( "  j = %zu: state %.3e, psi* %.3e, q* %.3e, trace %.3e, psi_ax %.6e vs %.6e\n", j,
		             rel( states[ 0 ], states[ 1 ], dx ),
		             rel( *fields[ 0 ].psiStar, *fields[ 1 ].psiStar, *lin.psiStar ),
		             rel( *fields[ 0 ].qStar, *fields[ 1 ].qStar, *lin.qStar ),
		             rel( *fields[ 0 ].trace, *fields[ 1 ].trace, *lin.trace ),
		             ( axes[ 0 ] - axes[ 1 ] )/( 2.0*h ), ds.front() );
	}
}

class ScaledGGPrime : public meq::Profile
{
	public:
		explicit ScaledGGPrime( double f ) : factor( f ) {}
		double operator()( double psi ) const override { return factor*closed( psi ); }
		double prime( double psi ) const override { return factor*closed.prime( psi ); }
		double doublePrime( double psi ) const override { return factor*closed.doublePrime( psi ); }
	private:
		double factor;
		ClosedGGPrime closed;
};

BOOST_AUTO_TEST_CASE( innerProbe, * boost::unit_test::disabled() )
{
	std::size_t const m = 8;
	std::vector<double> labels = meq::fluxSurfaceLabels( familyOptions( m ) );
	for ( double &l : labels )
		l = l*l;
	Representation const rep( labels );
	std::vector<double> exact( m );
	for ( std::size_t i = 0; i < m; ++i )
		exact[ i ] = closedGSquared( 1.0 - labels[ i ] );
	auto closed = std::make_shared<ClosedGGPrime const>();
	int const level = std::getenv( "QR_LEVEL" ) ? std::atoi( std::getenv( "QR_LEVEL" ) ) : 0;
	for ( int kind = 0; kind < 3; ++kind )
		for ( double f : { 1.0, 1.01, 1.03, 1.1, 0.97, 0.9 } )
		{
			Machine machine( level, 2, closed );
			BOOST_TEST_REQUIRE( machine.solve( closed ) );
			std::shared_ptr<meq::Profile const> profile;
			char const *what = "";
			if ( kind == 0 )
			{
				profile = std::make_shared<ScaledGGPrime const>( f );
				what = "closed form x f";
			}
			else if ( kind == 1 )
			{
				std::vector<double> u( exact );
				for ( double &x : u )
					x = g0*g0 + f*( x - g0*g0 );
				profile = rep.ggPrime( u, g0*g0 );
				what = "spline of the closed form, gg' x f";
			}
			else
			{
				std::vector<double> u( exact );
				u[ 0 ] *= f;
				profile = rep.ggPrime( u, g0*g0 );
				what = "spline, innermost node x f";
			}
			bool const ok = machine.solve( profile );
			std::printf( "  %-36s f = %.2f: %s, %d its\n", what, f, ok ? "ok" : "FAILED", machine.solver->newtonIterations() );
		}
}

BOOST_AUTO_TEST_CASE( theMissFallsWithTheKnotsAndTheMesh, * boost::unit_test::disabled() )
{
	Reference reference( 3, 3 );
	int const order = std::getenv( "QR_ORDER" ) ? std::atoi( std::getenv( "QR_ORDER" ) ) : 2;
	int const maxLevel = std::getenv( "QR_MAXLEVEL" ) ? std::atoi( std::getenv( "QR_MAXLEVEL" ) ) : 2;
	std::vector<std::size_t> counts = { 6, 12, 24, 48 };
	if ( std::getenv( "QR_COUNTS" ) )
	{
		counts.clear();
		std::string list = std::getenv( "QR_COUNTS" );
		std::size_t at = 0;
		while ( at < list.size() )
		{
			std::size_t const comma = list.find( ',', at );
			counts.push_back( std::stoul( list.substr( at, comma - at ) ) );
			if ( comma == std::string::npos )
				break;
			at = comma + 1;
		}
	}
	std::printf( "\n  k = %d: level, surfaces m, Newton its, solves, max|g^2/g^2_exact - 1| on the cut, "
	             "on [0,1], max|q/q* - 1| off the knots, time\n", order );
	int const minLevel = std::getenv( "QR_MINLEVEL" ) ? std::atoi( std::getenv( "QR_MINLEVEL" ) ) : 0;
	for ( int level = minLevel; level <= maxLevel; ++level )
		for ( std::size_t mBase : counts )
		{
			// QR_DIAGONAL: the knots refine with the mesh, m = m0 2^level.
			std::size_t const m = std::getenv( "QR_DIAGONAL" ) ? mBase << level : mBase;
			meq::FluxFamilyOptions const options = familyOptions( m );
			std::vector<double> labels;
			std::vector<double> const target = reference.q( options, &labels );
			Representation const rep( labels );
			pinRow = [ &rep ]( std::vector<double> const &v )
			{
				return ( rep.of( v, g0*g0 )( 0.0 ) - g0*g0 )/( g0*g0 );
			};
			pinGradient = [ &rep ]()
			{
				std::vector<double> row;
				for ( std::size_t j = 0; j < rep.unknowns(); ++j )
					row.push_back( rep.cardinal( j )( 0.0 )/( g0*g0 ) );
				return row;
			};
			std::vector<double> start( m, g0*g0 );
			if ( axisNode() )
				start.push_back( g0*g0 );
			if ( std::getenv( "QR_START_EXACT" ) )
			{
				// The closed form at the nodes: the study of the answer, apart
				// from the globalisation that reaches it.
				for ( std::size_t i = 0; i < m; ++i )
					start[ i ] = closedGSquared( 1.0 - labels[ i ] );
				if ( axisNode() )
					start[ m ] = closedGSquared( 1.0 );
			}
			std::size_t const n = start.size();
			Machine machine( level, order, rep.ggPrime( start, g0*g0 ) );
			auto const t0 = std::chrono::steady_clock::now();
			Outcome const outcome = newton( machine, rep, start, target, options, 1.0e-11, 12,
			                                std::getenv( "QR_VERBOSE" ) != nullptr );
			double const wall = seconds( t0 );
			if ( !outcome.converged )
			{
				std::printf( "    %d  %3zu  NOT CONVERGED after %d, history", level, m, outcome.iterations );
				for ( double h : outcome.history )
					std::printf( " %.2e", h );
				std::printf( "\n" );
				continue;
			}
			Miss const miss = measure( machine, rep, outcome.u, reference, m );
			if ( leastSquares() )
			{
				double q = 0.0;
				for ( std::size_t i = 0; i + 1 < outcome.misfit.size(); ++i )
					q = std::max( q, std::abs( outcome.misfit[ i ] ) );
				std::printf( "      least squares: max | q/q* - 1 | at the knots %.3e, pin row %.3e\n", q,
				             outcome.misfit.empty() ? 0.0 : outcome.misfit.back() );
			}
			if ( std::getenv( "QR_CONDITION" ) )
			{
				// The pinned Jacobian at the answer: || J^-1 ||_inf and the
				// direction it amplifies most, in u relative to u.
				Timing scratch;
				std::vector<double> const jac = tangentJacobian( machine, rep, outcome.u, target, options, 1.0e-5, scratch );
				double jNorm = 0.0;
				for ( std::size_t i = 0; i < n; ++i )
				{
					double row = 0.0;
					for ( std::size_t j = 0; j < n; ++j )
						row += std::abs( jac[ i*n + j ] );
					jNorm = std::max( jNorm, row );
				}
				std::vector<double> inverse( n*n );
				for ( std::size_t j = 0; j < n; ++j )
				{
					std::vector<double> e( n, 0.0 );
					e[ j ] = 1.0;
					std::vector<double> const col = denseSolve( jac, e );
					for ( std::size_t i = 0; i < n; ++i )
						inverse[ i*n + j ] = col[ i ];
				}
				double iNorm = 0.0;
				for ( std::size_t i = 0; i < n; ++i )
				{
					double row = 0.0;
					for ( std::size_t j = 0; j < n; ++j )
						row += std::abs( inverse[ i*n + j ] );
					iNorm = std::max( iNorm, row );
				}
				// The response of u/u to a uniform relative change of the target q.
				std::vector<double> ones( n, 1.0 );
				std::vector<double> const response = denseSolve( jac, ones );
				std::printf( "      || J ||_inf %.3e, || J^-1 ||_inf %.3e, cond %.3e; du/u for dq/q = 1 everywhere:", jNorm,
				             iNorm, jNorm*iNorm );
				for ( std::size_t i = 0; i < n; ++i )
					std::printf( " %+.3e", response[ i ]/outcome.u[ i ] );
				std::printf( "\n" );
			}
			if ( std::getenv( "QR_VERBOSE" ) )
			{
				std::printf( "      u/g2_exact - 1 at the nodes:" );
				for ( std::size_t i = 0; i < m; ++i )
					std::printf( " %+.3e", outcome.u[ i ]/closedGSquared( 1.0 - labels[ i ] ) - 1.0 );
				if ( axisNode() )
					std::printf( "  axis %+.3e", outcome.u[ m ]/closedGSquared( 1.0 ) - 1.0 );
				std::printf( "\n      labels:" );
				for ( std::size_t i = 0; i < m; ++i )
					std::printf( " %.4f", labels[ i ] );
				std::printf( "\n" );
			}
			std::printf( "    %d  %3zu  %2d  %3d   %.3e   %.3e   %.3e   %.1f s (solve %.1f, linearise %.1f, "
			             "tangents %.1f, reconstruct %.1f, extract %.1f over %d)  history",
			             level, m, outcome.iterations, outcome.timing.solves, miss.gSquaredCut,
			             miss.gSquaredAll, miss.qOffKnots, wall, outcome.timing.solve,
			             outcome.timing.linearise, outcome.timing.tangents, outcome.timing.reconstruct,
			             outcome.timing.extract, outcome.timing.extractions );
			for ( double h : outcome.history )
				std::printf( " %.1e", h );
			std::printf( "\n" );
			std::fflush( stdout );
		}
}
