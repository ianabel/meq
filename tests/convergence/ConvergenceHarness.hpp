#ifndef MEQ_TESTS_CONVERGENCEHARNESS_HPP
#define MEQ_TESTS_CONVERGENCEHARNESS_HPP

/*
 * The machinery every convergence test in this directory shares: the rectangle
 * mesh, the meq::Source adapter over a tests/analytic fixture, one solve and one
 * measurement, the rate arithmetic, the tables, and the two kinds of assertion
 * on top of them.
 *
 * TWO KINDS, and the difference is the whole reason this header exists.
 *
 *   checkOrder()     an EXACT solution is available. Asserts a rate of k+1
 *                    against it, in psi and in the flux, plus an absolute
 *                    error ceiling -- because a rate is blind to a solution
 *                    wrong by a constant factor or a sign, which is the
 *                    failure CLAUDE.md's "Testing stance" is organised around.
 *
 *   checkSelfOrder() NO exact solution is available. Compares successive
 *                    refinements against each other,
 *                    Delta^k_j( f ) := || f_h^k - f_h^(k-1) ||_j, which is what
 *                    refs/HDG-GradShafranov.pdf does for its own Example 6, and
 *                    asserts that THAT falls at k+1. It is a weaker statement
 *                    -- the paper calls it "a pessimistic estimate of the true
 *                    error" -- and in particular it cannot see a convergent
 *                    solver converging to the wrong function. Never use it
 *                    where an exact solution exists.
 *
 * WHY THE SELF-CONVERGENCE NORM IS A SAMPLED ONE. Successive levels live on
 * different meshes, so || f_h^k - f_h^(k-1) || is not an integral either space
 * can compute: the coarse function would have to be represented on the fine
 * mesh first. What is done here instead is what the same paper does when it
 * reports its own "off the grid" errors -- both levels are evaluated on a FIXED
 * cloud of sample points, through mfem::Mesh::FindPoints, and the norms are
 * discrete:
 *
 *     Delta_2   = sqrt( ( |Omega|/N ) sum_i ( a_i - b_i )^2 )
 *     Delta_inf = max_i | a_i - b_i |
 *
 * The cloud is the same at every level and at every polynomial degree, so the
 * quadrature error it carries is a fixed bias and not noise, and it cancels out
 * of a ratio of two successive Deltas to the extent that it is the same bias.
 * It is deliberately NOT a random draw: a test that fails only sometimes is
 * worse than no test.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

namespace meq
{
namespace tests
{

	/// The axisymmetric box a benchmark is posed on. r must stay well away from
	/// zero: the operator carries a 1/r and several of the closed forms carry a
	/// log r.
	struct Rectangle
	{
		double rMin, rMax, zMin, zMax;

		double width() const
		{
			return rMax - rMin;
		}

		double height() const
		{
			return zMax - zMin;
		}

		double area() const
		{
			return width()*height();
		}
	};

	/// The box of examples/manufactured.toml, which encloses the ITER-like
	/// double-null boundary refs/HDG-GradShafranov.pdf poses its Example 5 on,
	/// and which SolovievConvergence.cpp uses as well.
	inline Rectangle standardBox()
	{
		return Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}

	/// A triangulated rectangle with n cells a side. Triangles rather than
	/// quadrilaterals, as in both papers.
	inline mfem::Mesh makeMesh( Rectangle const &box, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE,
		                                               false, box.width(),
		                                               box.height() );
		double const rMin = box.rMin;
		double const zMin = box.zMin;
		mesh.Transform( [ rMin, zMin ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	/// meq::Source is the interface the Newton path takes: F and dF/dpsi as
	/// plain doubles. The analytic fixtures in tests/analytic deliberately depend
	/// on neither MFEM nor src/meq -- Soloviev.hpp is the pattern -- so the join
	/// between them lives here rather than in either. Every fixture in that
	/// directory already spells f() and dFdPsi() the way meq::Source does, so one
	/// template covers all of them.
	///
	/// It is a thin forward and nothing else. In particular it does not apply a
	/// sign, a 1/r or a normalisation: meq::Source::f() is documented to be F as
	/// eq (2) writes it, and each fixture is documented to return F as its paper
	/// writes it, which is the same F. Anything clever here would be a convention
	/// change hidden in a test helper.
	template<typename Equilibrium>
	class EquilibriumSource : public meq::Source
	{
		public:
			explicit EquilibriumSource( Equilibrium const &eqIn )
				: eq( eqIn )
			{
			}

			double f( double r, double z, double psi ) const override
			{
				return eq.f( r, z, psi );
			}

			double dFdPsi( double r, double z, double psi ) const override
			{
				return eq.dFdPsi( r, z, psi );
			}

		private:
			Equilibrium eq;
	};

	/// The same thin forward, for a fixture whose profiles are functions of
	/// NORMALISED flux and which therefore carries a normalisation the SOLVER
	/// owns rather than the caller. setNormalisation() forwards to the fixture's
	/// setPsiAxis() and nothing else differs -- in particular there is still no
	/// sign, no 1/r and no scaling applied here.
	///
	/// The fixture is held by value, so the solver's writes do not reach the
	/// caller's copy. That is deliberate: the value that matters afterwards is
	/// GradShafranovSolver::psiAxis(), which is the unknown's converged value,
	/// and a fixture quietly left holding an intermediate normalisation would be
	/// a second and disagreeing answer to the same question.
	template<typename Equilibrium>
	class NormalisedEquilibriumSource : public meq::NormalisedSource
	{
		public:
			explicit NormalisedEquilibriumSource( Equilibrium const &eqIn )
				: eq( eqIn )
			{
			}

			double f( double r, double z, double psi ) const override
			{
				return eq.f( r, z, psi );
			}

			double dFdPsi( double r, double z, double psi ) const override
			{
				return eq.dFdPsi( r, z, psi );
			}

			void setNormalisation( double psiAxis ) override
			{
				eq.setPsiAxis( psiAxis );
			}

			double normalisation() const override
			{
				return eq.psiAxis();
			}

			/// The fixture as the solver last left it, for a caller that wants to
			/// evaluate F or dF/dpsi at the converged normalisation.
			Equilibrium const &equilibrium() const
			{
				return eq;
			}

		private:
			Equilibrium eq;
	};

	/// One point on a convergence curve measured against an exact solution.
	struct Measurement
	{
		double h;
		int traceDofs;
		double errorPsi;
		double errorFlux;
		int newtonIterations;
	};

	/// One point on a self-convergence curve: no errors, because there is
	/// nothing to compare against, but the solution sampled on the shared cloud
	/// so that successive levels can be compared with each other.
	struct SelfMeasurement
	{
		double h;
		int traceDofs;
		int newtonIterations;
		/// False if Newton ran out of iterations. GradShafranovSolver::solve()
		/// throws in that case and does not recover a solution, so the samples
		/// are empty -- but the residual history and the iteration count survive,
		/// and they are the record worth keeping. A caller must report a false
		/// here rather than retrying with a looser tolerance.
		bool converged;
		std::vector<double> residuals;
		/// psi_h at each cloud point, and the two flux components interleaved.
		std::vector<double> psiSamples;
		std::vector<double> fluxSamples;
		/// The extreme values of psi_h over the cloud, which is how a run that
		/// landed on the trivial branch, or wandered off, is recognised.
		double psiMin, psiMax;
	};

	/// One solve of a NORMALISED equilibrium, where psi on the magnetic axis is
	/// an unknown of the non-linear system rather than an input.
	struct NormalisedMeasurement
	{
		double h;
		int traceDofs;
		int newtonIterations;
		bool converged;
		std::vector<double> residuals;
		/// The converged normalisation, and the constraint residual
		/// psi_ax - max psi_h that says whether it is self consistent.
		double psiAxis;
		double constraint;
		/// The extreme nodal values of psi_h. psiMax is what psi_ax is
		/// constrained to equal, so the pair is the whole self-consistency
		/// statement.
		double psiMin, psiMax;
	};

	inline double rate( double coarseError, double fineError, double refinementRatio )
	{
		return std::log( coarseError/fineError )/std::log( refinementRatio );
	}

	/// The observed order of a Newton sequence between three consecutive
	/// residuals, log( r2/r1 )/log( r1/r0 ). Two for a quadratically convergent
	/// iteration once it is in the asymptotic regime, one for a linearly
	/// convergent one anywhere.
	inline double newtonOrder( double r0, double r1, double r2 )
	{
		return std::log( r2/r1 )/std::log( r1/r0 );
	}

	/// Four dyadic refinements from 4 cells a side, so three measured rates per
	/// quantity per order.
	inline std::vector<int> const &dyadicMeshes()
	{
		static std::vector<int> const sizes = { 4, 8, 16, 32 };
		return sizes;
	}

	/// k+1, less the slack allowed for a two-mesh rate estimate. Wider than the
	/// linear benchmark's 0.15 because there are no theoretical estimates for
	/// the non-linear case -- refs/HDG-GradShafranov.pdf says so in as many
	/// words -- and its own Table 5 wanders between 1.86 and 2.03 at k = 1.
	double const rateSlack = 0.2;

	/*
	 * ---------------------------------------------------------------------
	 * The sample cloud, for self-convergence
	 * ---------------------------------------------------------------------
	 */

	/// A fixed set of points at which every level is evaluated. Construct one
	/// per benchmark and share it across meshes and polynomial degrees.
	///
	/// WHY A LOW-DISCREPANCY SEQUENCE AND NOT A LATTICE. The function being
	/// measured, psi_h^k - psi_h^(k-1), is dominated by the coarser level's
	/// discretisation error, which oscillates on the scale of that level's
	/// elements. A regular lattice of sample points can resonate with that
	/// oscillation -- the two spacings are commensurate on a dyadic sequence --
	/// and then the sampled norm is biased differently at each refinement, which
	/// lands directly in the measured rate. A Halton sequence has no such scale,
	/// so what is left is an unbiased sampling error of order 1/sqrt( N ), which
	/// at N = 12000 is under one per cent and moves a rate by about 0.03.
	///
	/// It is a deterministic sequence, not a random draw. A test that fails only
	/// sometimes is worse than no test.
	class SampleCloud
	{
		public:
			/// @a count Halton points over @a box.
			explicit SampleCloud( Rectangle const &box, int count = 12000 )
			{
				for ( int i = 1; i <= count; ++i )
				{
					points.push_back( box.rMin + halton( i, 2 )*box.width() );
					points.push_back( box.zMin + halton( i, 3 )*box.height() );
				}
				weight = box.area()/static_cast<double>( size() );
			}

			/// An explicit list, for a domain that is not a rectangle. @a area is
			/// the measure the discrete L2 norm is scaled by.
			SampleCloud( std::vector<double> const &pointsIn, double area )
				: points( pointsIn )
			{
				weight = area/static_cast<double>( size() );
			}

			/// The i-th term of the van der Corput sequence in @a base, which is
			/// the one-dimensional building block of a Halton sequence.
			static double halton( int index, int base )
			{
				double result = 0.0;
				double fraction = 1.0/base;
				int i = index;
				while ( i > 0 )
				{
					result += fraction*( i % base );
					i /= base;
					fraction /= base;
				}
				return result;
			}

			int size() const
			{
				return static_cast<int>( points.size()/2 );
			}

			double r( int i ) const
			{
				return points[ 2*i ];
			}

			double z( int i ) const
			{
				return points[ 2*i + 1 ];
			}

			/// |Omega|/N, the weight of one point in the discrete L2 norm.
			double pointWeight() const
			{
				return weight;
			}

			/// Evaluate a scalar and a vector GridFunction at every point of the
			/// cloud. Points that fall outside @a mesh -- which happens on a
			/// polygonal subdomain, not on a fitted rectangle -- are recorded as
			/// quiet NaN and skipped by the norms.
			///
			/// mfem::Mesh::FindPoints is O( elements x points ): it compares every
			/// point against every element centre. That is what caps the cloud
			/// size here, and it is why a benchmark should not simply ask for more
			/// points when a rate looks noisy.
			void evaluate( mfem::Mesh &mesh, mfem::GridFunction const &scalar,
			               mfem::GridFunction const &vector,
			               std::vector<double> &scalarOut,
			               std::vector<double> &vectorOut ) const
			{
				int const n = size();
				mfem::DenseMatrix matrix( 2, n );
				for ( int i = 0; i < n; ++i )
				{
					matrix( 0, i ) = r( i );
					matrix( 1, i ) = z( i );
				}

				mfem::Array<int> elements;
				mfem::Array<mfem::IntegrationPoint> ips;
				mesh.FindPoints( matrix, elements, ips, false );

				scalarOut.assign( n, std::numeric_limits<double>::quiet_NaN() );
				vectorOut.assign( 2*n, std::numeric_limits<double>::quiet_NaN() );

				mfem::Vector value( 2 );
				for ( int i = 0; i < n; ++i )
				{
					if ( elements[ i ] < 0 )
						continue;
					scalarOut[ i ] = scalar.GetValue( elements[ i ], ips[ i ] );
					vector.GetVectorValue( elements[ i ], ips[ i ], value );
					vectorOut[ 2*i ] = value( 0 );
					vectorOut[ 2*i + 1 ] = value( 1 );
				}
			}

		private:
			std::vector<double> points;
			double weight;
	};

	/// The discrete L2 norm of a - b over the cloud, skipping any point where
	/// either is NaN.
	///
	/// Works for a scalar and for a vector alike: the vector sample arrays hold
	/// both components per point, so summing over entries already sums | . |^2
	/// pointwise, and the weight is |Omega|/N with N the number of POINTS -- which
	/// is what cloud.pointWeight() is.
	inline double sampledL2( SampleCloud const &cloud, std::vector<double> const &a,
	                         std::vector<double> const &b )
	{
		// min(), not a.size(): a level whose Newton iteration failed has no
		// samples at all, and the difference against it is not a number.
		std::size_t const n = std::min( a.size(), b.size() );
		if ( n == 0 )
			return std::numeric_limits<double>::quiet_NaN();

		double sum = 0.0;
		int counted = 0;
		for ( std::size_t i = 0; i < n; ++i )
		{
			if ( std::isnan( a[ i ] ) || std::isnan( b[ i ] ) )
				continue;
			double const d = a[ i ] - b[ i ];
			sum += d*d;
			++counted;
		}
		if ( counted == 0 )
			return std::numeric_limits<double>::quiet_NaN();
		return std::sqrt( cloud.pointWeight()*sum );
	}

	/// The discrete maximum norm of a - b over the cloud.
	inline double sampledMax( std::vector<double> const &a,
	                          std::vector<double> const &b )
	{
		std::size_t const n = std::min( a.size(), b.size() );
		if ( n == 0 )
			return std::numeric_limits<double>::quiet_NaN();

		double worst = 0.0;
		for ( std::size_t i = 0; i < n; ++i )
		{
			if ( std::isnan( a[ i ] ) || std::isnan( b[ i ] ) )
				continue;
			worst = std::max( worst, std::abs( a[ i ] - b[ i ] ) );
		}
		return worst;
	}

	/*
	 * ---------------------------------------------------------------------
	 * Solving and measuring
	 * ---------------------------------------------------------------------
	 */

	/// Solve once against an exact solution and measure the errors. The
	/// Dirichlet datum is the exact psi on all four sides; the source is handed
	/// over as a meq::Source, which is what selects the Newton path.
	template<typename Equilibrium>
	Measurement measure( Equilibrium const &eq, Rectangle const &box, int order,
	                     int n, std::vector<double> *residualHistory = nullptr )
	{
		mfem::Mesh mesh = makeMesh( box, n );
		EquilibriumSource<Equilibrium> source( eq );

		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                      mfem::Vector &value )
		{
			eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		if ( residualHistory )
			*residualHistory = solver.newtonResiduals();

		Measurement point;
		point.h = box.width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.newtonIterations = solver.newtonIterations();
		return point;
	}

	/// The non-linear ordering every study in this suite runs under, in ONE
	/// place so that a test which has to branch on it -- because the two
	/// orderings differ in which KINSOL strategy converges, say -- cannot drift
	/// from what measureSelf() actually asks for.
	inline meq::GradShafranovSolver::NonlinearOrdering defaultOrdering()
	{
		return meq::GradShafranovSolver::NonlinearOrdering::NPC;
	}

	/// Solve once with no exact solution to compare against, and sample the
	/// result on @a cloud. @a boundary supplies the Dirichlet datum, which for
	/// these benchmarks is a design choice rather than a restriction of a known
	/// solution -- see the callers.
	/// @param guess  optional starting point for Newton, as psi( r, z ). Null
	///               is the default and means "start from the Dirichlet data",
	///               which for a source vanishing at psi = 0 lands on the trivial
	///               branch -- see meq::GradShafranovSolver::setInitialGuess and
	///               homogeneousDataLandsOnTheTrivialBranch in
	///               PedestalConvergence.cpp.
	template<typename Equilibrium, typename BoundaryFunction>
	SelfMeasurement measureSelf( Equilibrium const &eq, Rectangle const &box,
	                             int order, int n, SampleCloud const &cloud,
	                             BoundaryFunction boundary,
	                             int maxNewtonIterations = 30,
	                             double relativeTolerance = 1.0e-12,
	                             mfem::Coefficient *guess = nullptr,
	                             meq::GradShafranovSolver::Globalisation glob =
	                                 meq::GradShafranovSolver::Globalisation::None,
	                             meq::GradShafranovSolver::LocalSolver local =
	                                 meq::GradShafranovSolver::LocalSolver::Newton,
	                             meq::GradShafranovSolver::NonlinearOrdering ordering =
	                                 meq::GradShafranovSolver::NonlinearOrdering::NPC )
	                             // Keep in step with defaultOrdering() above.
	{
		mfem::Mesh mesh = makeMesh( box, n );
		EquilibriumSource<Equilibrium> source( eq );

		mfem::FunctionCoefficient boundaryCoeff( [ &boundary ]( mfem::Vector const &x )
		{
			return boundary( x( 0 ), x( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( boundaryCoeff );
		// The relative tolerance is a parameter because it has to be: a benchmark
		// whose first residual is small -- a weak source, or homogeneous data --
		// cannot reach 1e-12 of it, because the floor the residual can actually
		// reach is set by the conditioning of the trace solve and not by the
		// stopping rule. Newton then bounces on round-off until the iteration cap
		// and is reported as a failure that is not one. MillerConvergence.cpp
		// carries the measurement: 1.9e-14 absolute against a first residual of
		// 1.55e-3.
		solver.setNewtonControl( relativeTolerance, 1.0e-14, maxNewtonIterations );
		if ( guess )
			solver.setInitialGuess( *guess );
		solver.setGlobalisation( glob );
		solver.setLocalSolver( local );
		solver.setNonlinearOrdering( ordering );

		SelfMeasurement point;
		point.h = box.width()/static_cast<double>( n );
		point.converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const &e )
		{
			// A Newton failure is a result, not an accident: it is reported and
			// the run goes on. Nothing is retried with a looser tolerance.
			point.converged = false;
		}

		point.traceDofs = solver.numTraceDofs();
		point.newtonIterations = solver.newtonIterations();
		point.residuals = solver.newtonResiduals();
		point.psiMin = std::numeric_limits<double>::infinity();
		point.psiMax = -std::numeric_limits<double>::infinity();

		if ( point.converged )
		{
			cloud.evaluate( mesh, solver.potential(), solver.flux(),
			                point.psiSamples, point.fluxSamples );
			for ( double v : point.psiSamples )
			{
				if ( std::isnan( v ) )
					continue;
				point.psiMin = std::min( point.psiMin, v );
				point.psiMax = std::max( point.psiMax, v );
			}
		}
		return point;
	}

	/**
	 * One solve of an equilibrium whose profiles are functions of NORMALISED
	 * flux, with psi on the magnetic axis carried as an unknown.
	 *
	 * @param guess          the Newton starting point, as psi( r, z ). NOT
	 *                       optional, and not an optimisation: at a fixed
	 *                       normalisation this equation has a small solution and
	 *                       a large one, only the large one can satisfy
	 *                       max psi = psi_ax, and an unguided Newton lands on the
	 *                       small one. See
	 *                       GradShafranovSolver::setSource( NormalisedSource &, double ).
	 * @param psiAxisGuess   the starting normalisation, likewise.
	 *
	 * A failure is reported and not retried, exactly as in measureSelf(): the
	 * residual history and the iteration count survive it and are the record
	 * worth keeping.
	 */
	template<typename Equilibrium, typename BoundaryFunction>
	NormalisedMeasurement measureNormalised( Equilibrium const &eq, Rectangle const &box,
	                                         int order, int n,
	                                         BoundaryFunction boundary,
	                                         double psiAxisGuess,
	                                         mfem::Coefficient &guess,
	                                         int maxNewtonIterations = 40,
	                                         double relativeTolerance = 1.0e-10 )
	{
		mfem::Mesh mesh = makeMesh( box, n );
		NormalisedEquilibriumSource<Equilibrium> source( eq );

		mfem::FunctionCoefficient boundaryCoeff( [ &boundary ]( mfem::Vector const &x )
		{
			return boundary( x( 0 ), x( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source, psiAxisGuess );
		solver.setBoundaryData( boundaryCoeff );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( relativeTolerance, 1.0e-14, maxNewtonIterations );

		NormalisedMeasurement point;
		point.h = box.width()/static_cast<double>( n );
		point.converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			point.converged = false;
		}

		point.traceDofs = solver.numTraceDofs();
		point.newtonIterations = solver.newtonIterations();
		point.residuals = solver.newtonResiduals();
		point.psiAxis = solver.psiAxis();
		point.constraint = solver.normalisationResidual();
		point.psiMin = 0.0;
		point.psiMax = 0.0;
		if ( point.converged )
		{
			point.psiMin = solver.potential().Min();
			point.psiMax = solver.potential().Max();
		}
		return point;
	}

	/*
	 * ---------------------------------------------------------------------
	 * Tables
	 * ---------------------------------------------------------------------
	 */

	inline void printTable( char const *label, int order, Rectangle const &box,
	                        std::vector<Measurement> const &points )
	{
		std::printf( "\n  %s, k = %d, triangles on [%.1f,%.1f]x[%.1f,%.1f]\n",
		             label, order, box.rMin, box.rMax, box.zMin, box.zMax );
		std::printf( "  %8s %9s %14s %7s %14s %7s %7s\n",
		             "h", "trace", "L2(psi)", "rate", "L2(q)", "rate", "Newton" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %14.6e %7s %14.6e %7s %7d\n",
				             p.h, p.traceDofs, p.errorPsi, "-", p.errorFlux, "-",
				             p.newtonIterations );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( "  %8.5f %9d %14.6e %7.3f %14.6e %7.3f %7d\n",
				             p.h, p.traceDofs,
				             p.errorPsi, rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorFlux, rate( points[ i - 1 ].errorFlux, p.errorFlux, ratio ),
				             p.newtonIterations );
			}
		}
		std::fflush( stdout );
	}

	/// One row of the self-convergence table, the difference between two
	/// consecutive levels in both norms and both unknowns.
	struct SelfDifference
	{
		double h;
		double l2Psi, maxPsi;
		double l2Flux, maxFlux;
	};

	inline std::vector<SelfDifference>
	selfDifferences( SampleCloud const &cloud,
	                 std::vector<SelfMeasurement> const &points )
	{
		std::vector<SelfDifference> out;
		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			SelfDifference d;
			d.h = points[ i ].h;
			d.l2Psi = sampledL2( cloud, points[ i ].psiSamples,
			                     points[ i - 1 ].psiSamples );
			d.maxPsi = sampledMax( points[ i ].psiSamples, points[ i - 1 ].psiSamples );
			d.l2Flux = sampledL2( cloud, points[ i ].fluxSamples,
			                      points[ i - 1 ].fluxSamples );
			d.maxFlux = sampledMax( points[ i ].fluxSamples,
			                        points[ i - 1 ].fluxSamples );
			out.push_back( d );
		}
		return out;
	}

	inline void printSelfTable( char const *label, int order,
	                            std::vector<SelfMeasurement> const &points,
	                            std::vector<SelfDifference> const &diffs )
	{
		// Both norms, as refs/HDG-GradShafranov.pdf reports them: its
		// Delta^k_j( f ) is indexed by j in { 2, infinity }. The rate columns are
		// on the L2 differences, which are the ones asserted on; the maximum
		// differences are printed because a localised feature that the L2 norm
		// averages away -- a pedestal, an internal layer -- shows up there first.
		std::printf( "\n  %s, k = %d -- self convergence, "
		             "Delta^k( f ) = || f_h^k - f_h^(k-1) ||\n", label, order );
		std::printf( "  %8s %9s %7s %11s %11s %12s %6s %12s %6s %11s %11s\n",
		             "h", "trace", "Newton", "min psi", "max psi",
		             "D2(psi)", "rate", "D2(q)", "rate", "Dinf(psi)", "Dinf(q)" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			SelfMeasurement const &p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %7d %11.4e %11.4e %12s %6s %12s %6s "
				             "%11s %11s\n",
				             p.h, p.traceDofs, p.newtonIterations, p.psiMin, p.psiMax,
				             "-", "-", "-", "-", "-", "-" );
				continue;
			}
			SelfDifference const &d = diffs[ i - 1 ];
			if ( i == 1 )
			{
				std::printf( "  %8.5f %9d %7d %11.4e %11.4e %12.4e %6s %12.4e %6s "
				             "%11.4e %11.4e\n",
				             p.h, p.traceDofs, p.newtonIterations, p.psiMin, p.psiMax,
				             d.l2Psi, "-", d.l2Flux, "-", d.maxPsi, d.maxFlux );
			}
			else
			{
				SelfDifference const &previous = diffs[ i - 2 ];
				double const ratio = previous.h/d.h;
				std::printf( "  %8.5f %9d %7d %11.4e %11.4e %12.4e %6.3f %12.4e %6.3f "
				             "%11.4e %11.4e\n",
				             p.h, p.traceDofs, p.newtonIterations, p.psiMin, p.psiMax,
				             d.l2Psi, rate( previous.l2Psi, d.l2Psi, ratio ),
				             d.l2Flux, rate( previous.l2Flux, d.l2Flux, ratio ),
				             d.maxPsi, d.maxFlux );
			}
		}
		std::fflush( stdout );
	}

	/// The Newton residual history of one run, printed with its observed order.
	/// A history that grinds down linearly means the Jacobian disagrees with the
	/// residual, and no amount of mesh refinement fixes it.
	inline void printNewtonHistory( char const *label, int order, double h,
	                                std::vector<double> const &history )
	{
		std::printf( "\n  %s, Newton residual history at k = %d, h = %.5f\n",
		             label, order, h );
		std::printf( "  %5s %16s %16s %8s\n", "it", "||r||", "||r||/||r_0||", "order" );
		for ( std::size_t i = 0; i < history.size(); ++i )
		{
			double const relative = history.front() > 0.0
			                        ? history[ i ]/history.front() : 0.0;
			if ( i >= 2 && history[ i - 1 ] > 0.0 && history[ i - 2 ] > 0.0 )
			{
				std::printf( "  %5zu %16.6e %16.6e %8.3f\n", i, history[ i ], relative,
				             newtonOrder( history[ i - 2 ], history[ i - 1 ],
				                          history[ i ] ) );
			}
			else
			{
				std::printf( "  %5zu %16.6e %16.6e %8s\n", i, history[ i ], relative,
				             "-" );
			}
		}
		std::fflush( stdout );
	}

	/// The best observed Newton order over the triples that are above the
	/// round-off floor. Below the floor the ratios are noise.
	inline double bestNewtonOrder( std::vector<double> const &history,
	                              double relativeFloor = 1.0e-12 )
	{
		if ( history.size() < 3 )
			return 0.0;
		double const floor = relativeFloor*history.front();
		double best = 0.0;
		for ( std::size_t i = 2; i < history.size(); ++i )
		{
			if ( history[ i ] < floor || history[ i - 1 ] < floor )
				continue;
			best = std::max( best, newtonOrder( history[ i - 2 ], history[ i - 1 ],
			                                    history[ i ] ) );
		}
		return best;
	}

	/*
	 * ---------------------------------------------------------------------
	 * Assertions
	 * ---------------------------------------------------------------------
	 */

	/// Measure a whole dyadic sequence against an exact solution and assert k+1
	/// in psi and in the flux, plus an absolute ceiling on the finest mesh.
	///
	/// The ceilings are not decoration. The rate is blind to a solution wrong by
	/// a constant factor or a sign, and a wrong convention converges at the right
	/// rate to the wrong function -- so the ceiling is the only assertion here
	/// that can see that. They sit at roughly three times the measured value and
	/// each caller records the measurement it was set from.
	template<typename Equilibrium>
	void checkOrder( Equilibrium const &eq, char const *label, int order,
	                 double psiCeiling, double fluxCeiling,
	                 Rectangle const &box = standardBox(),
	                 std::vector<int> const &sizes = dyadicMeshes(),
	                 double slack = rateSlack )
	{
		std::vector<Measurement> points;
		points.reserve( sizes.size() );
		for ( int n : sizes )
			points.push_back( measure( eq, box, order, n ) );

		printTable( label, order, box, points );

		double const expected = order + 1.0 - slack;

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const ratePsi = rate( points[ i - 1 ].errorPsi,
			                             points[ i ].errorPsi, ratio );
			double const rateFlux = rate( points[ i - 1 ].errorFlux,
			                              points[ i ].errorFlux, ratio );

			BOOST_TEST( ratePsi >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi converged at " << ratePsi << ", wanted " << expected );
			BOOST_TEST( rateFlux >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": q converged at " << rateFlux << ", wanted " << expected );
		}

		BOOST_TEST( points.back().errorPsi < psiCeiling,
		            "k = " << order << ": L2 error in psi is " << points.back().errorPsi
		            << ", above the ceiling " << psiCeiling );
		BOOST_TEST( points.back().errorFlux < fluxCeiling,
		            "k = " << order << ": L2 error in q is " << points.back().errorFlux
		            << ", above the ceiling " << fluxCeiling );
	}

	/// Measure a dyadic sequence with no exact solution to compare against, and
	/// assert that the difference between consecutive levels falls at k+1.
	///
	/// WHAT THIS CAN AND CANNOT SEE. It is a Cauchy statement about the sequence,
	/// not a statement about its limit: a solver converging steadily to the wrong
	/// function passes it. That is why it is used only where no closed form
	/// exists, and why every benchmark that has one uses checkOrder() instead.
	/// refs/HDG-GradShafranov.pdf makes the same point about its own Example 6 --
	/// "even if the relative difference between two successive approximations is a
	/// pessimistic estimate of the true error, the method still performs
	/// satisfactorily".
	///
	/// It does see, and this is what it is for here: a Newton iteration that fails
	/// or stalls, a run that has fallen onto a different solution branch between
	/// one refinement and the next -- which shows up as a difference that does not
	/// shrink at all -- and a spatial rate that is short of design order.
	///
	/// Returns the measurements so that a caller can assert more.
	template<typename Equilibrium, typename BoundaryFunction>
	std::vector<SelfMeasurement>
	checkSelfOrderAgainst( Equilibrium const &eq, char const *label, int order,
	                       SampleCloud const &cloud, BoundaryFunction boundary,
	                       Rectangle const &box, std::vector<int> const &sizes,
	                       double psiFloor, double fluxFloor,
	                       int maxNewtonIterations = 30,
	                       double relativeTolerance = 1.0e-12 )
	{
		std::vector<SelfMeasurement> points;
		points.reserve( sizes.size() );
		for ( int n : sizes )
			points.push_back( measureSelf( eq, box, order, n, cloud, boundary,
			                               maxNewtonIterations,
			                               relativeTolerance ) );

		std::vector<SelfDifference> diffs = selfDifferences( cloud, points );
		printSelfTable( label, order, points, diffs );

		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			BOOST_TEST( points[ i ].converged,
			            label << ", k = " << order << ", h = " << points[ i ].h
			            << ": Newton did NOT converge in " << points[ i ].newtonIterations
			            << " iterations. That is a finding about this benchmark, not a "
			            "tolerance to be relaxed -- but SINCE meq USES "
			            "NonlinearOrdering::NPC, MEASURE THE OTHER ORDERING BEFORE "
			            "BLAMING THE BENCHMARK. setNonlinearOrdering( "
			            "CondenseThenLinearise ) is one line and is kept as the "
			            "backup for exactly this. THE PARITY GAP THIS MESSAGE USED "
			            "TO NAME IS GONE WITH ITS CAUSE: it was a property of "
			            "MFEM's LineariseThenCondense, a trace-only operator that "
			            "kept the linearisation as hidden state, and upstream "
			            "deleted that mode. NPC holds no state between calls at "
			            "all, so a residual is a function of its argument and "
			            "there is no frozen-Jacobian local correction to truncate. "
			            "A failure here is therefore a NEW finding rather than a "
			            "known one -- record which orderings reach it. Do not "
			            "clear it by lowering the cap or switching the ordering "
			            "back" );
		}

		// A RATE TAKEN ACROSS A SOLVE THAT DID NOT CONVERGE IS NOT A RATE. Such a
		// pair yields NaN, NaN fails every floor, and the three assertions that
		// follow then bury the one message above that says what actually went
		// wrong -- measured, a single non-converged point produced three failures
		// where one was informative. So the pair is skipped and said to be
		// skipped. This is not a relaxation: the convergence assertion above has
		// already fired for that point, and the rate is unmeasurable rather than
		// merely disappointing.
		for ( std::size_t i = 1; i < diffs.size(); ++i )
		{
			bool measurable = true;
			for ( SelfMeasurement const &point : points )
				if ( !point.converged )
					measurable = false;

			if ( !measurable )
			{
				std::printf( "    rates not asserted: a solve in this sequence did "
				             "not converge, so the differences are not errors\n" );
				break;
			}

			double const ratio = diffs[ i - 1 ].h/diffs[ i ].h;
			double const ratePsi = rate( diffs[ i - 1 ].l2Psi, diffs[ i ].l2Psi, ratio );
			double const rateFlux = rate( diffs[ i - 1 ].l2Flux, diffs[ i ].l2Flux, ratio );

			BOOST_TEST( ratePsi >= psiFloor,
			            label << ", k = " << order << ", h = " << diffs[ i ].h
			            << ": Delta( psi ) fell at " << ratePsi << ", wanted "
			            << psiFloor );
			BOOST_TEST( rateFlux >= fluxFloor,
			            label << ", k = " << order << ", h = " << diffs[ i ].h
			            << ": Delta( q ) fell at " << rateFlux << ", wanted "
			            << fluxFloor );
		}

		return points;
	}

	/// The common case: both floors at k+1 less @a slack.
	template<typename Equilibrium, typename BoundaryFunction>
	std::vector<SelfMeasurement>
	checkSelfOrder( Equilibrium const &eq, char const *label, int order,
	                SampleCloud const &cloud, BoundaryFunction boundary,
	                Rectangle const &box = standardBox(),
	                std::vector<int> const &sizes = dyadicMeshes(),
	                double slack = rateSlack, int maxNewtonIterations = 30 )
	{
		return checkSelfOrderAgainst( eq, label, order, cloud, boundary, box, sizes,
		                              order + 1.0 - slack, order + 1.0 - slack,
		                              maxNewtonIterations );
	}

}
}

#endif // MEQ_TESTS_CONVERGENCEHARNESS_HPP
