/*
 * InversionScaling -- the cost of the solution-inversion chain, which had never
 * been measured. INVERSION-PLAN.md stage IN-P.
 *
 * This is a BENCHMARK, not a test: it is built but NOT registered with ctest,
 * because everything it reports is a timing and CLAUDE.md's standing rule is
 * that a threaded timing on this machine is a measurement about the machine.
 * What it DOES assert -- and the only thing it exits non-zero for -- are the
 * correctness properties that make the timings mean anything:
 *
 *   * a threaded extraction reproduces the serial one BIT FOR BIT, at
 *     0.000e+00 and not at a tolerance. Independent surfaces and independent
 *     rays reassociate nothing, so exactness is available, and a tolerance
 *     would be an admission that something is shared. The precedent is
 *     SolverContract::threadedAssemblyReproducesSerialAssemblyExactly;
 *   * every candidate KERNEL below reproduces the answer of the code it would
 *     replace. A faster wrong answer is not a result, and this file prices four
 *     replacements without adopting any of them.
 *
 * Everything in INVERSION-PLAN.md section 7 is an ACCURACY measurement. The one
 * cost datum that existed before this file is that FluxSurfaceConvergence runs
 * in 16.6 s, which is a test and not a workload. So the questions here are the
 * ones section 11 poses and could not answer:
 *
 *   1. where the time actually goes, per stage and as a share;
 *   2. which hand-written loops a library kernel beats, AND BY HOW MUCH --
 *      "this is three per cent of the run" is a perfectly good finding and
 *      stops somebody spending a week on it;
 *   3. what threading over surfaces and over rays is worth, against the
 *      bit-exactness assertion;
 *   4. whether continuation in the flux label -- which is in direct tension
 *      with parallelism over surfaces, since it makes them sequential -- pays
 *      for the independence it gives up.
 *
 * MKL_NUM_THREADS=1 IS NOT NEGOTIABLE HERE AND THIS FILE DOES NOT GET TO RELAX
 * IT. CLAUDE.md records ComputeH()'s element-local dense LU degrading by a
 * factor of forty at k = 3 under threaded MKL, and the variable is
 * process-wide, so the parallelism measured here is OpenMP over independent
 * work and never threaded BLAS. A run that appears to gain from raising MKL
 * threads is measuring the solve getting slower somewhere else.
 *
 * Run it through tests/performance/inversion-scan.sh, which drives one process
 * per thread count for the reason TraceSolverScaling's scan.sh does: MKL and
 * the OpenMP runtime both fix state at first use, so an in-process sweep
 * measures the first setting repeatedly.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/SurfaceAverage.hpp"
#include "meq/SurfaceFit.hpp"
#include "meq/Zernike.hpp"

#include "analytic/Soloviev.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MEQ_HAVE_EIGEN
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#endif

namespace
{

	using Clock = std::chrono::steady_clock;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	double const twoPi = 6.283185307179586476925286766559;

	double seconds( Clock::time_point start )
	{
		return std::chrono::duration<double>( Clock::now() - start ).count();
	}

	/// Best of @a repeats, which is the right statistic for a timing whose noise
	/// is one-sided: a run can be slowed by the machine and cannot be sped up by
	/// it, so the minimum is the closest thing to the cost of the work. Lifted
	/// from TraceSolverScaling.cpp, and deliberately the same statistic so the
	/// two harnesses' numbers can be read side by side.
	template <typename F>
	double bestOf( int repeats, F &&f )
	{
		double best = 1.0e300;
		for ( int i = 0; i < repeats; ++i )
		{
			auto start = Clock::now();
			f();
			best = std::min( best, seconds( start ) );
		}
		return best;
	}

	struct Box
	{
		double rMin, rMax, zMin, zMax;
	};

	/// NOT the standard benchmark rectangle, for the reason
	/// SurfaceAverageConvergence.cpp records: nstx()'s axis sits at r = 1.318
	/// and its surfaces are elongated, so Psi_N = 0.50 reaches r in
	/// [ 0.81, 1.66 ] and does not fit in [0.6,1.4]x[-0.6,0.6] at all.
	Box nstxBox()
	{
		return Box{ 0.60, 1.90, -1.10, 1.10 };
	}

	/// Duplicated from tests/convergence/ConvergenceHarness.hpp rather than
	/// included, because that header pulls in Boost.Test and this binary is not
	/// a test. TraceSolverScaling.cpp does the same and for the same reason.
	mfem::Mesh makeMesh( Box const &b, int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE,
		                                               false, b.rMax - b.rMin,
		                                               b.zMax - b.zMin );
		double const rMin = b.rMin;
		double const zMin = b.zMin;
		mesh.Transform( [ rMin, zMin ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	struct ExactAxis
	{
		double r;
		double z;
		double psi;
	};

	/// The magnetic axis of the CLOSED FORM, to round-off. Newton on
	/// grad( psi ) = 0 with the Hessian by central differences; the Hessian's
	/// accuracy does not reach the answer, the fixed point being where the
	/// analytic gradient vanishes whatever steered it there.
	ExactAxis exactAxis( Equilibrium const &eq, double rGuess, double zGuess )
	{
		double r = rGuess;
		double z = zGuess;
		double const step = 1.0e-5;

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );

			double a0 = 0.0, a1 = 0.0, b0 = 0.0, b1 = 0.0;
			double hessian[ 2 ][ 2 ];
			eq.gradPsi( r + step, z, a0, b0 );
			eq.gradPsi( r - step, z, a1, b1 );
			hessian[ 0 ][ 0 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 0 ] = ( b0 - b1 )/( 2.0*step );
			eq.gradPsi( r, z + step, a0, b0 );
			eq.gradPsi( r, z - step, a1, b1 );
			hessian[ 0 ][ 1 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 1 ] = ( b0 - b1 )/( 2.0*step );

			double const det = hessian[ 0 ][ 0 ]*hessian[ 1 ][ 1 ]
			                   - hessian[ 0 ][ 1 ]*hessian[ 1 ][ 0 ];
			r += -(  hessian[ 1 ][ 1 ]*gr - hessian[ 0 ][ 1 ]*gz )/det;
			z += -( -hessian[ 1 ][ 0 ]*gr + hessian[ 0 ][ 0 ]*gz )/det;
		}

		return ExactAxis{ r, z, eq.psi( r, z ) };
	}

	/// psi at normalised flux @a fraction, from the CLOSED FORM so that it is
	/// the same level at every mesh and every degree -- which is what lets a
	/// cost sweep in h compare like with like.
	double levelAt( ExactAxis const &axis, double fraction )
	{
		return axis.psi*( 1.0 - fraction );
	}

	/// One solve, kept alive, post-processed so that both candidate potentials
	/// exist. Member order is load bearing and the class is non-copyable because
	/// the coefficients capture `this`.
	class Solved
	{
		public:
			Solved( Equilibrium const &eqIn, Box const &boxIn, int orderIn, int n )
				: eq( eqIn ),
				  mesh( makeMesh( boxIn, n ) ),
				  sourceCoeff( [ this ]( mfem::Vector const &x )
				  {
					  return eq.f( x( 0 ), x( 1 ), 0.0 );
				  } ),
				  psiCoeff( [ this ]( mfem::Vector const &x )
				  {
					  return eq.psi( x( 0 ), x( 1 ) );
				  } ),
				  solver( mesh, orderIn )
			{
				solveSeconds = bestOf( 1, [ this ]()
				{
					solver.setSource( sourceCoeff );
					solver.setBoundaryData( psiCoeff );
					solver.solve();
				} );
				postProcessSeconds = bestOf( 1, [ this ]()
				{
					solver.postProcess();
				} );
			}

			Solved( Solved const & ) = delete;
			Solved &operator=( Solved const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			mfem::Mesh &theMesh()
			{
				return mesh;
			}

			double solveCost() const
			{
				return solveSeconds;
			}

			double postProcessCost() const
			{
				return postProcessSeconds;
			}

		private:
			Equilibrium eq;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			meq::GradShafranovSolver solver;
			double solveSeconds = 0.0;
			double postProcessSeconds = 0.0;
	};

	/// The workload this file is about: a family of nested surfaces, each traced,
	/// each populated by rays, each turned into a quadrature. It is what
	/// MANTA-COUPLING.md's consumer asks for and what IN-6 will cache.
	struct Family
	{
		std::vector<double> normalisedFlux;
		std::vector<meq::Contour> contours;
		std::vector<meq::AngleParametrisation> fits;
		std::vector<meq::SurfaceAverages> averages;
	};

	struct FamilySpec
	{
		int surfaces = 12;
		int angles = 48;
		double smallest = 0.15;
		double largest = 0.50;
	};

	std::vector<double> fluxLadder( FamilySpec const &spec )
	{
		std::vector<double> out;
		for ( int i = 0; i < spec.surfaces; ++i )
			out.push_back( spec.smallest + ( spec.largest - spec.smallest )
			                               *i/( spec.surfaces - 1.0 ) );
		return out;
	}

	/// Trace one surface, populate it by rays, and build its quadrature. The
	/// three stages the breakdown separates, in the order a consumer runs them.
	void extractOne( meq::ContourTracer const &tracer,
	                 meq::CriticalPoint const &axis, ExactAxis const &exact,
	                 double normalisedFlux, std::size_t angles,
	                 meq::Contour &contour, meq::AngleParametrisation &fit,
	                 meq::SurfaceAverages &average )
	{
		contour = tracer.traceFromAxis( levelAt( exact, normalisedFlux ), axis );
		fit = tracer.fitByAngle( contour, axis, angles );
		average = meq::surfaceAverages( tracer, fit );
	}

	Family extractFamily( meq::ContourTracer const &tracer,
	                      meq::CriticalPoint const &axis, ExactAxis const &exact,
	                      FamilySpec const &spec )
	{
		Family family;
		family.normalisedFlux = fluxLadder( spec );
		std::size_t const count = family.normalisedFlux.size();
		family.contours.resize( count );
		family.fits.resize( count );
		family.averages.resize( count );

		for ( std::size_t i = 0; i < count; ++i )
			extractOne( tracer, axis, exact, family.normalisedFlux[ i ],
			            static_cast<std::size_t>( spec.angles ),
			            family.contours[ i ], family.fits[ i ],
			            family.averages[ i ] );

		return family;
	}

	/// The sample cloud a SurfaceFit is built from: every ray node of every
	/// surface, relabelled by the axis shape -- which IN-3 measured to be worth
	/// 45x to 660x and is not optional.
	std::vector<meq::SurfaceSample> samplesOf( Family const &family,
	                                           meq::CriticalPoint const &axis )
	{
		std::vector<meq::SurfaceSample> geometric;
		for ( std::size_t i = 0; i < family.fits.size(); ++i )
		{
			meq::AngleParametrisation const &rays = family.fits[ i ];
			for ( std::size_t j = 0; j < rays.count(); ++j )
			{
				meq::SurfaceSample sample;
				sample.normalisedFlux = family.normalisedFlux[ i ];
				sample.theta = twoPi*static_cast<double>( j )/rays.count();
				sample.r = rays.pointR[ j ];
				sample.z = rays.pointZ[ j ];
				geometric.push_back( sample );
			}
		}

		meq::AxisShape const shape = meq::axisShapeFromSamples( geometric, axis.r,
		                                                        axis.z );
		return meq::relabelByAxisShape( geometric, shape );
	}

	/// The seam gaugeFreeFit() reads the solved field through: one lambda over
	/// ContourTracer::sampleAt(), with an element hint threaded through it
	/// because CLAUDE.md records Mesh::FindPoints as O( elements x points ) and
	/// this is invoked once per node per Gauss-Newton residual evaluation.
	meq::NormalisedFluxField solvedField( meq::ContourTracer const &tracer,
	                                      ExactAxis const &exact, int &hint )
	{
		meq::NormalisedFluxField field;
		field.sample = [ &tracer, &hint, &exact ]( double r, double z,
		                                           double &normalisedFlux,
		                                           double &gradientR,
		                                           double &gradientZ )
		{
			double psi = 0.0, fluxR = 0.0, fluxZ = 0.0;
			if ( !tracer.sampleAt( r, z, psi, fluxR, fluxZ, hint ) )
				return false;

			normalisedFlux = 1.0 - psi/exact.psi;
			gradientR = -r*fluxR/exact.psi;
			gradientZ = -r*fluxZ/exact.psi;
			return true;
		};
		return field;
	}

	// -----------------------------------------------------------------------
	// THE CANDIDATE KERNEL: every Zernike mode at one point, by recurrence.
	//
	// meq::zernike() is one Boost jacobi() AND one jacobi_prime() per mode --
	// radialValue() computes both and zernikeRadial() throws the second away --
	// plus a checkMode(), an integerPower() and a trigonometric call. A design
	// matrix row is that once per mode; SurfaceFit::position() is THREE of them
	// per mode, because evaluateMode() asks for the value and both derivatives
	// separately and each re-enters radialValue(). So a fit evaluated at a point
	// costs about six special-function evaluations per mode.
	//
	// The radial polynomials at one rho satisfy Kintner's three-term recurrence
	// in l at fixed m, so ALL of them cost O( modes ) arithmetic operations
	// TOGETHER. This is the standard construction and is not novel; what is
	// worth measuring is the size of the difference on meq's own workload and
	// whether the recurrence agrees with the Boost route to something better
	// than the fit's own residual.
	//
	// IT IS PRICED HERE AND NOT ADOPTED. Zernike.cpp's route is exact, is what
	// every accuracy number in INVERSION-PLAN.md section 7 was measured with,
	// and its header records that the naive factorial sum -- a DIFFERENT
	// alternative -- loses eight digits at l = 30. A recurrence is not that
	// sum, but the burden is on the replacement.
	// -----------------------------------------------------------------------
	class ZernikeTable
	{
		public:
			explicit ZernikeTable( int maxDegreeIn )
				: maxDegree( maxDegreeIn ),
				  modes( meq::zernikeModes( maxDegreeIn ) ),
				  radial( static_cast<std::size_t>( maxDegreeIn + 1 )
				          *static_cast<std::size_t>( maxDegreeIn + 1 ), 0.0 ),
				  radialPrime( radial.size(), 0.0 ),
				  cosine( static_cast<std::size_t>( maxDegreeIn + 1 ), 0.0 ),
				  sine( static_cast<std::size_t>( maxDegreeIn + 1 ), 0.0 ),
				  value( modes.size(), 0.0 ),
				  dRadius( modes.size(), 0.0 ),
				  dTheta( modes.size(), 0.0 )
			{
			}

			std::vector<meq::ZernikeMode> const &modeList() const
			{
				return modes;
			}

			/// Every mode, its radial derivative and its angular derivative, at
			/// one ( rho, theta ). O( modes ) arithmetic and no special function.
			void evaluate( double rho, double theta )
			{
				int const top = maxDegree;
				std::size_t const stride = static_cast<std::size_t>( top + 1 );

				// R_m^m = rho^m and its derivative, then Kintner in l.
				double power = 1.0;
				for ( int m = 0; m <= top; ++m )
				{
					radial[ at( m, m ) ] = power;
					radialPrime[ at( m, m ) ] =
						( m == 0 ) ? 0.0 : m*power/safeRho( rho );
					power *= rho;
				}

				for ( int m = 0; m + 2 <= top; ++m )
				{
					double const rm = radial[ at( m, m ) ];
					double const rm2 = rm*rho*rho;
					radial[ at( m + 2, m ) ] = ( m + 2 )*rm2 - ( m + 1 )*rm;
					radialPrime[ at( m + 2, m ) ] =
						( m + 2 )*( m + 2 )*rm2/safeRho( rho )
						- ( m + 1 )*m*rm/safeRho( rho );
				}

				for ( int m = 0; m <= top; ++m )
				{
					for ( int l = m + 4; l <= top; l += 2 )
					{
						double const k1 = ( l + m )*( l - m )*( l - 2 )/2.0;
						double const k2 = 2.0*l*( l - 1 )*( l - 2 );
						double const k3 = -static_cast<double>( m )*m*( l - 1 )
						                  - static_cast<double>( l )*( l - 1 )*( l - 2 );
						double const k4 = -static_cast<double>( l )*( l + m - 2 )
						                  *( l - m - 2 )/2.0;

						double const a = k2*rho*rho + k3;
						radial[ at( l, m ) ] = ( a*radial[ at( l - 2, m ) ]
						                         + k4*radial[ at( l - 4, m ) ] )/k1;
						radialPrime[ at( l, m ) ] =
							( 2.0*k2*rho*radial[ at( l - 2, m ) ]
							  + a*radialPrime[ at( l - 2, m ) ]
							  + k4*radialPrime[ at( l - 4, m ) ] )/k1;
					}
				}

				// cos( m theta ) and sin( m theta ) by the Chebyshev recurrence,
				// which is two multiplications a mode against a library call.
				double const c1 = std::cos( theta );
				double const s1 = std::sin( theta );
				cosine[ 0 ] = 1.0;
				sine[ 0 ] = 0.0;
				if ( top >= 1 )
				{
					cosine[ 1 ] = c1;
					sine[ 1 ] = s1;
				}
				for ( int m = 2; m <= top; ++m )
				{
					cosine[ m ] = 2.0*c1*cosine[ m - 1 ] - cosine[ m - 2 ];
					sine[ m ] = 2.0*c1*sine[ m - 1 ] - sine[ m - 2 ];
				}

				for ( std::size_t j = 0; j < modes.size(); ++j )
				{
					int const l = modes[ j ].l;
					int const m = modes[ j ].m;
					int const absM = ( m < 0 ) ? -m : m;
					double const rad = radial[ at( l, absM ) ];
					double const radPrime = radialPrime[ at( l, absM ) ];

					if ( m >= 0 )
					{
						value[ j ] = rad*cosine[ absM ];
						dRadius[ j ] = radPrime*cosine[ absM ];
						dTheta[ j ] = -absM*rad*sine[ absM ];
					}
					else
					{
						value[ j ] = rad*sine[ absM ];
						dRadius[ j ] = radPrime*sine[ absM ];
						dTheta[ j ] = absM*rad*cosine[ absM ];
					}
				}

				(void)stride;
			}

			std::vector<double> const &values() const
			{
				return value;
			}

			std::vector<double> const &radialDerivatives() const
			{
				return dRadius;
			}

			std::vector<double> const &angularDerivatives() const
			{
				return dTheta;
			}

		private:
			std::size_t at( int l, int m ) const
			{
				return static_cast<std::size_t>( l )
				       *static_cast<std::size_t>( maxDegree + 1 )
				       + static_cast<std::size_t>( m );
			}

			/// rho = 0 is an ordinary point of every admissible mode -- that is
			/// the whole content of the parity constraint -- but the CLOSED FORM
			/// for the base cases divides by rho. The numerator vanishes there
			/// for every m >= 1, so this only has to keep 0/0 from becoming NaN.
			/// Zernike.cpp meets the same thing and writes the m = 0 branch out
			/// explicitly for the same reason.
			static double safeRho( double rho )
			{
				return ( rho == 0.0 ) ? 1.0 : rho;
			}

			int maxDegree;
			std::vector<meq::ZernikeMode> modes;
			std::vector<double> radial;
			std::vector<double> radialPrime;
			std::vector<double> cosine;
			std::vector<double> sine;
			std::vector<double> value;
			std::vector<double> dRadius;
			std::vector<double> dTheta;
	};

	/// What the harness returns to main(): a count of correctness failures, so
	/// that a section can fail loudly without ending the run and losing the
	/// tables after it.
	struct Verdict
	{
		int failures = 0;

		void check( bool condition, char const *what )
		{
			if ( !condition )
			{
				std::printf( "    *** %s\n", what );
				++failures;
			}
		}
	};

}

// ===========================================================================
// 1. WHERE THE TIME GOES
// ===========================================================================

namespace
{

	void printCorrectorHistogram( Family const &family )
	{
		int worst = 0;
		std::vector<long> bucket( 40, 0 );
		long points = 0;
		long stalled = 0;
		long fallbacks = 0;
		double worstResidual = 0.0;

		for ( meq::Contour const &contour : family.contours )
		{
			stalled += contour.stalledCorrections;
			fallbacks += contour.fallbackLocations;
			worstResidual = std::max( worstResidual, contour.worstResidual );
			for ( meq::ContourPoint const &point : contour.points )
			{
				int const its = std::min( point.correctorIterations, 39 );
				++bucket[ its ];
				++points;
				worst = std::max( worst, point.correctorIterations );
			}
		}

		long cumulative = 0;
		double total = 0.0;
		for ( std::size_t i = 0; i < bucket.size(); ++i )
			total += static_cast<double>( i )*bucket[ i ];

		std::printf( "\n  corrector iterations per accepted point -- the"
		             " DISTRIBUTION, because a mean hides the tail\n" );
		std::printf( "    %4s  %9s  %8s  %8s\n", "its", "points", "share", "cum" );
		for ( std::size_t i = 0; i < bucket.size(); ++i )
		{
			if ( bucket[ i ] == 0 )
				continue;
			cumulative += bucket[ i ];
			std::printf( "    %4zu  %9ld  %7.2f%%  %7.2f%%\n", i, bucket[ i ],
			             100.0*bucket[ i ]/points, 100.0*cumulative/points );
		}
		std::printf( "    %ld points, mean %.3f, worst %d\n", points,
		             total/points, worst );
		std::printf( "    stalledCorrections %ld, fallbackLocations %ld,"
		             " worst residual %.3e\n", stalled, fallbacks, worstResidual );

		long rayIterations = 0;
		long rayNodes = 0;
		int worstRay = 0;
		long bisections = 0;
		long stalledRays = 0;
		long rayFallbacks = 0;
		for ( meq::AngleParametrisation const &fit : family.fits )
		{
			rayIterations += fit.totalIterations;
			rayNodes += static_cast<long>( fit.count() );
			worstRay = std::max( worstRay, fit.worstIterations );
			bisections += fit.bisections;
			stalledRays += fit.stalledRays;
			rayFallbacks += fit.fallbackLocations;
		}
		std::printf( "    ray Newton: %ld nodes, %.3f iterations each, worst %d,"
		             " %ld bisections, %ld stalled\n", rayNodes,
		             static_cast<double>( rayIterations )/rayNodes, worstRay,
		             bisections, stalledRays );
		std::printf( "    ray fallbackLocations %ld -- SEPARATELY from the"
		             " trace's, because they are different loops and only one of\n"
		             "    them is documented to be at zero\n", rayFallbacks );
	}

	void breakdown( Solved &solved, ExactAxis const &exact,
	                FamilySpec const &spec, int degree, int repeats,
	                Verdict &verdict )
	{
		meq::CriticalPointFinder finder( solved.theSolver() );

		double const axisCost = bestOf( repeats, [ & ]()
		{
			meq::CriticalPoint const found = finder.findAxis();
			(void)found;
		} );
		meq::CriticalPoint const axis = finder.findAxis();

		meq::ContourTracer tracer( solved.theSolver() );
		std::vector<double> const ladder = fluxLadder( spec );
		std::size_t const angles = static_cast<std::size_t>( spec.angles );

		// The three stages, timed one at a time on the same surfaces rather than
		// by subtracting a total: each is re-run from the state the one before it
		// left, so nothing is charged to the wrong column.
		std::vector<meq::Contour> contours( ladder.size() );
		double const traceCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < ladder.size(); ++i )
				contours[ i ] = tracer.traceFromAxis( levelAt( exact, ladder[ i ] ),
				                                      axis );
		} );

		std::vector<meq::AngleParametrisation> fits( ladder.size() );
		double const fitCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < ladder.size(); ++i )
				fits[ i ] = tracer.fitByAngle( contours[ i ], axis, angles );
		} );

		std::vector<meq::SurfaceAverages> averages( ladder.size() );
		double const averageCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < ladder.size(); ++i )
				averages[ i ] = meq::surfaceAverages( tracer, fits[ i ] );
		} );

		Family family;
		family.normalisedFlux = ladder;
		family.contours = contours;
		family.fits = fits;
		family.averages = averages;

		std::vector<meq::SurfaceSample> const samples = samplesOf( family, axis );

		meq::SurfaceFitOptions options;
		options.discEdge = spec.largest;

		double const linearFitCost = bestOf( repeats, [ & ]()
		{
			meq::SurfaceFit const fit( degree, samples, options );
			(void)fit;
		} );
		meq::SurfaceFit const linear( degree, samples, options );

		// The residual pass of the constructor, timed separately through the
		// public API, because it is the same arithmetic: one position() per
		// sample. That is what lets the QR and the Jacobi sweep be read off by
		// difference, which is stated rather than hidden -- there is no seam in
		// SurfaceFit that would let them be timed directly.
		double const positionPassCost = bestOf( repeats, [ & ]()
		{
			double sink = 0.0;
			for ( meq::SurfaceSample const &sample : samples )
			{
				double r = 0.0, z = 0.0;
				linear.position( sample.normalisedFlux, sample.theta, r, z );
				sink += r + z;
			}
			if ( sink == 1.0e300 )
				std::printf( " " );
		} );

		// And the design-matrix assembly, replicated exactly: sqrt( weight )
		// times meq::zernike() over zernikeModes( degree ), which is what
		// SurfaceFit's constructor does with modeValue() on the Zernike branch.
		std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( degree );
		std::vector<double> design( samples.size()*modes.size(), 0.0 );
		double const assemblyCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < samples.size(); ++i )
			{
				double const argument =
					std::sqrt( samples[ i ].normalisedFlux/options.discEdge );
				double const scale = std::sqrt( samples[ i ].weight );
				for ( std::size_t j = 0; j < modes.size(); ++j )
					design[ i*modes.size() + j ] =
						scale*meq::zernike( modes[ j ].l, modes[ j ].m, argument,
						                    samples[ i ].theta );
			}
		} );

		int hint = -1;
		meq::NormalisedFluxField const field = solvedField( tracer, exact, hint );
		std::vector<meq::DiscNode> const nodes = meq::discNodesFrom( samples );

		meq::GaugeFreeFitOptions gauge;
		meq::GaugeFreeFitReport report;
		double const gaugeCost = bestOf( repeats, [ & ]()
		{
			meq::GaugeFreeFitReport local;
			meq::SurfaceFit const freed = meq::gaugeFreeFit( linear, field, nodes,
			                                                 gauge, local );
			(void)freed;
		} );
		meq::SurfaceFit const freed = meq::gaugeFreeFit( linear, field, nodes,
		                                                 gauge, report );

		double const chain = traceCost + fitCost + averageCost + linearFitCost
		                     + gaugeCost + axisCost;

		std::printf( "\n=== 1. WHERE THE TIME GOES ===\n" );
		std::printf( "    %d surfaces, %d angles, Psi_N in [ %.2f, %.2f ],"
		             " fit degree L = %d\n", spec.surfaces, spec.angles,
		             spec.smallest, spec.largest, degree );
		std::printf( "    %zu samples, %zu modes, %d ray nodes in the family\n\n",
		             samples.size(), modes.size(),
		             spec.surfaces*spec.angles );

		std::printf( "  %-34s %11s %8s %14s\n", "stage", "seconds", "share",
		             "per surface" );
		auto row = []( char const *name, double cost, double whole, int surfaces )
		{
			std::printf( "  %-34s %11.5f %7.2f%% %14.6f\n", name, cost,
			             100.0*cost/whole, cost/surfaces );
		};
		row( "findAxis (once per family)", axisCost, chain, spec.surfaces );
		row( "trace  (predictor-corrector)", traceCost, chain, spec.surfaces );
		row( "fitByAngle  (the ray solves)", fitCost, chain, spec.surfaces );
		row( "surfaceAverages  (quadrature)", averageCost, chain, spec.surfaces );
		row( "SurfaceFit  (linear, L2)", linearFitCost, chain, spec.surfaces );
		row( "gaugeFreeFit  (Gauss-Newton)", gaugeCost, chain, spec.surfaces );
		std::printf( "  %-34s %11.5f %7.2f%%\n", "TOTAL extraction", chain, 100.0 );
		std::printf( "\n  for scale, and NOT part of the chain:"
		             " solve %.4f s, postProcess %.4f s\n",
		             solved.solveCost(), solved.postProcessCost() );

		std::printf( "\n  inside SurfaceFit's constructor, the two pieces that"
		             " can be timed directly and the one that cannot\n" );
		std::printf( "    design matrix, meq::zernike() per mode   %10.5f s"
		             "  %6.2f%% of the fit\n", assemblyCost,
		             100.0*assemblyCost/linearFitCost );
		std::printf( "    residual pass, position() per sample     %10.5f s"
		             "  %6.2f%%\n", positionPassCost,
		             100.0*positionPassCost/linearFitCost );
		std::printf( "    QR + one-sided Jacobi, BY DIFFERENCE     %10.5f s"
		             "  %6.2f%%\n",
		             linearFitCost - assemblyCost - positionPassCost,
		             100.0*( linearFitCost - assemblyCost - positionPassCost )
		                 /linearFitCost );

		// WHERE gaugeFreeFit's TIME GOES, and it is not the linear algebra.
		// The field arrives as one callable, so wrapping it counts and times
		// every evaluation exactly -- no sampling, no profiler, and the number
		// is attributable to a named function rather than to a stack.
		long samples2 = 0;
		double fieldSeconds = 0.0;
		int wrappedHint = -1;
		meq::NormalisedFluxField instrumented;
		instrumented.sample = [ &tracer, &wrappedHint, &exact, &samples2,
		                        &fieldSeconds ]( double r, double z,
		                                         double &normalisedFlux,
		                                         double &gradientR,
		                                         double &gradientZ )
		{
			auto const start = Clock::now();
			double psi = 0.0, fluxR = 0.0, fluxZ = 0.0;
			bool const ok = tracer.sampleAt( r, z, psi, fluxR, fluxZ,
			                                 wrappedHint );
			fieldSeconds += seconds( start );
			++samples2;
			if ( !ok )
				return false;
			normalisedFlux = 1.0 - psi/exact.psi;
			gradientR = -r*fluxR/exact.psi;
			gradientZ = -r*fluxZ/exact.psi;
			return true;
		};

		double instrumentedCost = 0.0;
		{
			meq::GaugeFreeFitReport local;
			auto const start = Clock::now();
			meq::SurfaceFit const wrapped = meq::gaugeFreeFit( linear,
			                                                   instrumented,
			                                                   nodes, gauge,
			                                                   local );
			instrumentedCost = seconds( start );
			(void)wrapped;
		}

		std::printf( "\n  gaugeFreeFit: %d Gauss-Newton iterations [%s],"
		             " %.5f s each\n", report.iterations, report.stop,
		             gaugeCost/std::max( 1, report.iterations ) );
		std::printf( "    surface error %.3e -> %.3e, %zu columns,"
		             " %zu soft directions, %d rejected steps\n",
		             report.initialSurfaceError, report.surfaceError,
		             report.columns, report.softDirections,
		             report.rejectedSteps );
		std::printf( "    %ld field evaluations, %.5f s in them -- %.1f%% of the"
		             " %.5f s that run took.\n"
		             "    THE LINEAR ALGEBRA IS NOT WHERE THE TIME IS: the QR and"
		             " the Jacobi sweep are the other %.1f%%.\n",
		             samples2, fieldSeconds, 100.0*fieldSeconds/instrumentedCost,
		             instrumentedCost,
		             100.0*( 1.0 - fieldSeconds/instrumentedCost ) );

		// AND WHAT ONE LINE BUYS, end to end. setWalkDepth() is a public knob
		// with a default of four, and section 3's K5 shows the walk missing on a
		// third of the ray nodes at that depth -- every miss an O( elements )
		// Mesh::FindPoints. This is the same chain with the depth raised, and
		// the answer has to be identical because FindPoints returns the element
		// the walk would have found.
		tracer.setWalkDepth( 12 );
		std::vector<meq::AngleParametrisation> deepFits( ladder.size() );
		double const deepFitCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < ladder.size(); ++i )
				deepFits[ i ] = tracer.fitByAngle( contours[ i ], axis, angles );
		} );
		double const deepAverageCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < ladder.size(); ++i )
			{
				meq::SurfaceAverages const a =
					meq::surfaceAverages( tracer, deepFits[ i ] );
				(void)a;
			}
		} );
		meq::GaugeFreeFitReport deepReport;
		double const deepGaugeCost = bestOf( repeats, [ & ]()
		{
			meq::GaugeFreeFitReport local;
			meq::SurfaceFit const f = meq::gaugeFreeFit( linear, field, nodes,
			                                             gauge, local );
			(void)f;
		} );
		meq::SurfaceFit const deepFreed = meq::gaugeFreeFit( linear, field, nodes,
		                                                      gauge, deepReport );
		(void)deepFreed;

		double worstRadius = 0.0;
		for ( std::size_t i = 0; i < ladder.size(); ++i )
			for ( std::size_t j = 0; j < deepFits[ i ].count(); ++j )
				worstRadius = std::max( worstRadius,
					std::fabs( deepFits[ i ].radius[ j ]
					           - fits[ i ].radius[ j ] ) );
		tracer.setWalkDepth( 4 );

		std::printf( "\n  setWalkDepth( 4 ) against setWalkDepth( 12 ), same"
		             " chain, same answer\n" );
		std::printf( "    %-22s %11s %11s %9s\n", "", "depth 4", "depth 12",
		             "speedup" );
		std::printf( "    %-22s %11.5f %11.5f %8.2fx\n", "fitByAngle", fitCost,
		             deepFitCost, fitCost/deepFitCost );
		std::printf( "    %-22s %11.5f %11.5f %8.2fx\n", "surfaceAverages",
		             averageCost, deepAverageCost, averageCost/deepAverageCost );
		std::printf( "    %-22s %11.5f %11.5f %8.2fx\n", "gaugeFreeFit",
		             gaugeCost, deepGaugeCost, gaugeCost/deepGaugeCost );
		std::printf( "    the fitted radii differ by %.3e and the surface error"
		             " by %.3e%s\n", worstRadius,
		             std::fabs( deepReport.surfaceError - report.surfaceError ),
		             ( worstRadius == 0.0 ) ? "   (the same answer)" : "" );
		verdict.check( worstRadius == 0.0,
		               "raising setWalkDepth() changed the answer. The walk and"
		               " Mesh::FindPoints must locate the same element" );

		printCorrectorHistogram( family );

		// A LOUD PRINT AND NOT A VERDICT FAILURE, deliberately.
		// FluxSurfaces.hpp is explicit that fallbackLocations "is a performance
		// statement, not a correctness one, since FindPoints returns the same
		// element", so failing on it would be asserting on a timing -- which is
		// the one thing this file is not allowed to do. It is worth shouting
		// about because Mesh::FindPoints is O( elements x points ): it is a
		// cliff rather than a slope, and no timing table would say why it moved.
		long traceFallbacks = 0;
		long fitFallbacks = 0;
		for ( meq::Contour const &contour : family.contours )
			traceFallbacks += contour.fallbackLocations;
		for ( meq::AngleParametrisation const &fit : family.fits )
			fitFallbacks += fit.fallbackLocations;
		if ( traceFallbacks + fitFallbacks > 0 )
			std::printf( "\n  !!! Mesh::FindPoints was reached %ld times from the"
			             " trace and %ld from the ray fit. It is\n"
			             "      O( elements x points ), so a count that scales"
			             " with the mesh turns a linear loop quadratic.\n",
			             traceFallbacks, fitFallbacks );
		(void)verdict;

		std::printf( "\n  the fit that came out of it: worst residual"
		             " %.3e, condition %.3e, min Jacobian %+.3e\n",
		             std::max( linear.diagnostics().worstR,
		                       linear.diagnostics().worstZ ),
		             linear.diagnostics().conditionNumber,
		             report.minimumJacobian );
		(void)freed;
	}

}

// ===========================================================================
// 2. SCALING
// ===========================================================================

namespace
{

	/// One number for the whole chain, so that a sweep has something to plot.
	double chainCost( Solved &solved, ExactAxis const &exact,
	                  FamilySpec const &spec, int degree, int repeats )
	{
		meq::CriticalPointFinder finder( solved.theSolver() );
		meq::CriticalPoint const axis = finder.findAxis();
		meq::ContourTracer tracer( solved.theSolver() );

		meq::SurfaceFitOptions options;
		options.discEdge = spec.largest;

		return bestOf( repeats, [ & ]()
		{
			Family const family = extractFamily( tracer, axis, exact, spec );
			std::vector<meq::SurfaceSample> const samples = samplesOf( family,
			                                                            axis );
			meq::SurfaceFit const fit( degree, samples, options );
			(void)fit;
		} );
	}

	void scaling( Equilibrium const &eq, ExactAxis const &exact,
	              std::vector<int> const &orders, std::vector<int> const &sizes,
	              FamilySpec const &base, int degree, int repeats )
	{
		std::printf( "\n=== 2. SCALING ===\n" );

		std::printf( "\n  (a) polynomial degree k and mesh size n, family and fit"
		             " held fixed\n" );
		std::printf( "    %2s %5s %9s %11s %11s %11s %11s %11s\n", "k", "n",
		             "elements", "solve", "trace", "rays", "averages", "chain" );

		for ( int order : orders )
		{
			for ( int n : sizes )
			{
				Solved solved( eq, nstxBox(), order, n );
				meq::CriticalPointFinder finder( solved.theSolver() );
				meq::CriticalPoint const axis = finder.findAxis();
				meq::ContourTracer tracer( solved.theSolver() );

				std::vector<double> const ladder = fluxLadder( base );
				std::size_t const angles = static_cast<std::size_t>( base.angles );

				std::vector<meq::Contour> contours( ladder.size() );
				double const traceCost = bestOf( repeats, [ & ]()
				{
					for ( std::size_t i = 0; i < ladder.size(); ++i )
						contours[ i ] = tracer.traceFromAxis(
							levelAt( exact, ladder[ i ] ), axis );
				} );

				std::vector<meq::AngleParametrisation> fits( ladder.size() );
				double const fitCost = bestOf( repeats, [ & ]()
				{
					for ( std::size_t i = 0; i < ladder.size(); ++i )
						fits[ i ] = tracer.fitByAngle( contours[ i ], axis, angles );
				} );

				double const averageCost = bestOf( repeats, [ & ]()
				{
					for ( std::size_t i = 0; i < ladder.size(); ++i )
					{
						meq::SurfaceAverages const a =
							meq::surfaceAverages( tracer, fits[ i ] );
						(void)a;
					}
				} );

				double const whole = chainCost( solved, exact, base, degree,
				                                repeats );

				std::printf( "    %2d %5d %9d %11.5f %11.5f %11.5f %11.5f %11.5f\n",
				             order, n, solved.theMesh().GetNE(),
				             solved.solveCost(), traceCost, fitCost, averageCost,
				             whole );
				std::fflush( stdout );
			}
		}

		// The three parameters a consumer chooses, against the k and n it does
		// not: one solve, reused, so what moves is the extraction alone.
		Solved solved( eq, nstxBox(), 2, 48 );
		meq::CriticalPointFinder finder( solved.theSolver() );
		meq::CriticalPoint const axis = finder.findAxis();
		meq::ContourTracer tracer( solved.theSolver() );

		std::printf( "\n  (b) points per surface N, at k = 2, n = 48, %d surfaces\n",
		             base.surfaces );
		std::printf( "    %5s %11s %11s %11s\n", "N", "trace", "rays", "averages" );
		for ( int angles : { 24, 48, 96, 192, 384 } )
		{
			FamilySpec spec = base;
			spec.angles = angles;
			std::vector<double> const ladder = fluxLadder( spec );

			std::vector<meq::Contour> contours( ladder.size() );
			double const traceCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
					contours[ i ] = tracer.traceFromAxis(
						levelAt( exact, ladder[ i ] ), axis );
			} );

			std::vector<meq::AngleParametrisation> fits( ladder.size() );
			double const fitCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
					fits[ i ] = tracer.fitByAngle(
						contours[ i ], axis,
						static_cast<std::size_t>( angles ) );
			} );

			double const averageCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
				{
					meq::SurfaceAverages const a =
						meq::surfaceAverages( tracer, fits[ i ] );
					(void)a;
				}
			} );

			std::printf( "    %5d %11.5f %11.5f %11.5f\n", angles, traceCost,
			             fitCost, averageCost );
			std::fflush( stdout );
		}

		std::printf( "\n  (c) number of surfaces, at k = 2, n = 48, N = %d\n",
		             base.angles );
		std::printf( "    %5s %11s %11s %11s\n", "n_psi", "trace", "rays",
		             "averages" );
		for ( int surfaces : { 4, 8, 16, 32 } )
		{
			FamilySpec spec = base;
			spec.surfaces = surfaces;
			std::vector<double> const ladder = fluxLadder( spec );
			std::size_t const angles = static_cast<std::size_t>( spec.angles );

			std::vector<meq::Contour> contours( ladder.size() );
			double const traceCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
					contours[ i ] = tracer.traceFromAxis(
						levelAt( exact, ladder[ i ] ), axis );
			} );

			std::vector<meq::AngleParametrisation> fits( ladder.size() );
			double const fitCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
					fits[ i ] = tracer.fitByAngle( contours[ i ], axis, angles );
			} );

			double const averageCost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < ladder.size(); ++i )
				{
					meq::SurfaceAverages const a =
						meq::surfaceAverages( tracer, fits[ i ] );
					(void)a;
				}
			} );

			std::printf( "    %5d %11.5f %11.5f %11.5f\n", surfaces, traceCost,
			             fitCost, averageCost );
			std::fflush( stdout );
		}

		std::printf( "\n  (d) fit degree L, on ONE fixed sample cloud"
		             " (%d surfaces x %d angles)\n", base.surfaces, base.angles );
		std::printf( "    %4s %8s %11s %11s %11s %13s\n", "L", "modes",
		             "SurfaceFit", "gaugeFree", "per GN it", "condition" );

		Family const family = extractFamily( tracer, axis, exact, base );
		std::vector<meq::SurfaceSample> const samples = samplesOf( family, axis );
		int hint = -1;
		meq::NormalisedFluxField const field = solvedField( tracer, exact, hint );
		std::vector<meq::DiscNode> const nodes = meq::discNodesFrom( samples );

		meq::SurfaceFitOptions options;
		options.discEdge = base.largest;

		for ( int L : { 4, 8, 12, 16, 20 } )
		{
			if ( meq::zernikeModeCount( L ) > samples.size() )
			{
				std::printf( "    %4d %8zu   (%zu samples cannot determine"
				             " them)\n", L, meq::zernikeModeCount( L ),
				             samples.size() );
				continue;
			}

			double const linearCost = bestOf( repeats, [ & ]()
			{
				meq::SurfaceFit const fit( L, samples, options );
				(void)fit;
			} );
			meq::SurfaceFit const linear( L, samples, options );

			meq::GaugeFreeFitOptions gauge;
			meq::GaugeFreeFitReport report;
			double const gaugeCost = bestOf( repeats, [ & ]()
			{
				meq::GaugeFreeFitReport local;
				meq::SurfaceFit const freed = meq::gaugeFreeFit( linear, field,
				                                                 nodes, gauge,
				                                                 local );
				(void)freed;
			} );
			meq::SurfaceFit const freed = meq::gaugeFreeFit( linear, field, nodes,
			                                                 gauge, report );
			(void)freed;

			std::printf( "    %4d %8zu %11.5f %11.5f %11.5f %13.3e\n", L,
			             meq::zernikeModeCount( L ), linearCost, gaugeCost,
			             gaugeCost/std::max( 1, report.iterations ),
			             linear.diagnostics().conditionNumber );
			std::fflush( stdout );
		}
	}

}

// ===========================================================================
// 3. LIBRARY KERNELS
// ===========================================================================

namespace
{

	void kernels( Solved &solved, ExactAxis const &exact, FamilySpec const &spec,
	              int degree, int repeats, Verdict &verdict )
	{
		meq::CriticalPointFinder finder( solved.theSolver() );
		meq::CriticalPoint const axis = finder.findAxis();
		meq::ContourTracer tracer( solved.theSolver() );

		Family const family = extractFamily( tracer, axis, exact, spec );
		std::vector<meq::SurfaceSample> const samples = samplesOf( family, axis );

		meq::SurfaceFitOptions options;
		options.discEdge = spec.largest;
		meq::SurfaceFit const linear( degree, samples, options );

		std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( degree );
		std::size_t const rows = samples.size();
		std::size_t const cols = modes.size();

		std::printf( "\n=== 3. LIBRARY KERNELS: what a hand-written loop is"
		             " costing, PRICED rather than asserted ===\n" );
		std::printf( "    L = %d, %zu modes, %zu samples\n", degree, cols, rows );

		// -------------------------------------------------------------------
		// K1. Zernike: Boost per mode, against a recurrence per point.
		// -------------------------------------------------------------------
		std::vector<double> byBoost( rows*cols, 0.0 );
		double const boostAll = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < rows; ++i )
			{
				double const rho =
					std::sqrt( samples[ i ].normalisedFlux/options.discEdge );
				for ( std::size_t j = 0; j < cols; ++j )
					byBoost[ i*cols + j ] = meq::zernike( modes[ j ].l,
					                                      modes[ j ].m, rho,
					                                      samples[ i ].theta );
			}
		} );

		ZernikeTable table( degree );
		std::vector<double> byRecurrence( rows*cols, 0.0 );
		double const recurrenceAll = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < rows; ++i )
			{
				double const rho =
					std::sqrt( samples[ i ].normalisedFlux/options.discEdge );
				table.evaluate( rho, samples[ i ].theta );
				std::vector<double> const &v = table.values();
				for ( std::size_t j = 0; j < cols; ++j )
					byRecurrence[ i*cols + j ] = v[ j ];
			}
		} );

		double worstMode = 0.0;
		for ( std::size_t i = 0; i < rows*cols; ++i )
			worstMode = std::max( worstMode,
			                      std::fabs( byBoost[ i ] - byRecurrence[ i ] ) );

		std::printf( "\n  K1  every mode at every sample, VALUES ONLY\n" );
		std::printf( "      meq::zernike() per mode (Boost jacobi)   %10.5f s\n",
		             boostAll );
		std::printf( "      Kintner recurrence, all modes per point  %10.5f s"
		             "   %6.1fx\n", recurrenceAll, boostAll/recurrenceAll );
		std::printf( "      worst absolute disagreement              %10.3e"
		             "   (modes are O(1))\n", worstMode );

		// The recurrence is a CANDIDATE and the burden is on it. 1e-11 is far
		// above the round-off of an O(1) quantity and far below anything a fit
		// residual could notice -- the discrete fits in section 7 sit at 1e-7
		// and better -- so this separates "a stable recurrence" from "a
		// recurrence that has started to cancel", which is the failure
		// Zernike.hpp records the explicit factorial sum suffering.
		verdict.check( worstMode < 1.0e-11,
		               "the Kintner recurrence disagrees with the Boost jacobi"
		               " route by more than 1e-11: it is not a usable"
		               " replacement at this degree" );

		// -------------------------------------------------------------------
		// K2. Evaluating a fit at many points: per-point loop against one GEMM.
		// -------------------------------------------------------------------
		std::size_t const evalPoints = 20000;
		std::vector<double> evalFlux( evalPoints, 0.0 );
		std::vector<double> evalTheta( evalPoints, 0.0 );
		for ( std::size_t i = 0; i < evalPoints; ++i )
		{
			evalFlux[ i ] = spec.smallest
				+ ( spec.largest - spec.smallest )*( i % 97 )/96.0;
			evalTheta[ i ] = twoPi*( i % 211 )/211.0;
		}

		std::vector<double> loopR( evalPoints, 0.0 );
		std::vector<double> loopZ( evalPoints, 0.0 );
		double const loopEval = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < evalPoints; ++i )
				linear.position( evalFlux[ i ], evalTheta[ i ], loopR[ i ],
				                 loopZ[ i ] );
		} );

		std::vector<double> const &cR = linear.majorRadiusCoefficients();
		std::vector<double> const &cZ = linear.heightCoefficients();

		std::vector<double> vandermonde( evalPoints*cols, 0.0 );
		std::vector<double> gemmR( evalPoints, 0.0 );
		std::vector<double> gemmZ( evalPoints, 0.0 );
		double const gemmEval = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < evalPoints; ++i )
			{
				double const rho = std::sqrt( evalFlux[ i ]/options.discEdge );
				table.evaluate( rho, evalTheta[ i ] );
				std::vector<double> const &v = table.values();
				for ( std::size_t j = 0; j < cols; ++j )
					vandermonde[ i*cols + j ] = v[ j ];
			}

			for ( std::size_t i = 0; i < evalPoints; ++i )
			{
				double r = 0.0;
				double z = 0.0;
				for ( std::size_t j = 0; j < cols; ++j )
				{
					r += vandermonde[ i*cols + j ]*cR[ j ];
					z += vandermonde[ i*cols + j ]*cZ[ j ];
				}
				gemmR[ i ] = r;
				gemmZ[ i ] = z;
			}
		} );

		double worstEval = 0.0;
		double scaleEval = 0.0;
		for ( std::size_t i = 0; i < evalPoints; ++i )
		{
			worstEval = std::max( worstEval,
			                      std::max( std::fabs( loopR[ i ] - gemmR[ i ] ),
			                                std::fabs( loopZ[ i ] - gemmZ[ i ] ) ) );
			scaleEval = std::max( scaleEval, std::fabs( loopR[ i ] ) );
		}

		std::printf( "\n  K2  the fit evaluated at %zu points\n", evalPoints );
		std::printf( "      SurfaceFit::position() per point         %10.5f s\n",
		             loopEval );
		std::printf( "      Vandermonde by recurrence, then a GEMM   %10.5f s"
		             "   %6.1fx\n", gemmEval, loopEval/gemmEval );
		std::printf( "      worst positional disagreement            %10.3e m"
		             "  (positions are O(1) m)\n", worstEval );
		verdict.check( worstEval < 1.0e-11*std::max( 1.0, scaleEval ),
		               "the batched evaluation does not reproduce"
		               " SurfaceFit::position()" );

		// -------------------------------------------------------------------
		// K3. The least-squares solve.
		// -------------------------------------------------------------------
		std::printf( "\n  K3  the least-squares solve on the %zu x %zu design"
		             " matrix\n", rows, cols );

		std::vector<double> rhsR( rows, 0.0 );
		std::vector<double> rhsZ( rows, 0.0 );
		for ( std::size_t i = 0; i < rows; ++i )
		{
			rhsR[ i ] = samples[ i ].r;
			rhsZ[ i ] = samples[ i ].z;
		}

#ifdef MEQ_HAVE_EIGEN
		// EIGEN IS SINGLE-THREADED HERE BY CONSTRUCTION AND THAT IS CHECKED
		// RATHER THAN ASSUMED. Eigen only ever parallelises general
		// matrix-matrix products, and Parallelizer.h bails to the sequential
		// path when omp_get_num_threads() > 1 -- so it does not nest inside an
		// OpenMP-over-surfaces region. Outside one it would use
		// omp_get_max_threads(), which on this machine is not what
		// MKL_NUM_THREADS=1 asks for, so it is pinned: a threaded GEMM here
		// would be a second thread team competing with the one over surfaces.
		Eigen::setNbThreads( 1 );

		Eigen::MatrixXd A( static_cast<Eigen::Index>( rows ),
		                   static_cast<Eigen::Index>( cols ) );
		for ( std::size_t i = 0; i < rows; ++i )
			for ( std::size_t j = 0; j < cols; ++j )
				A( static_cast<Eigen::Index>( i ),
				   static_cast<Eigen::Index>( j ) ) = byBoost[ i*cols + j ];

		Eigen::MatrixXd B( static_cast<Eigen::Index>( rows ), 2 );
		for ( std::size_t i = 0; i < rows; ++i )
		{
			B( static_cast<Eigen::Index>( i ), 0 ) = rhsR[ i ];
			B( static_cast<Eigen::Index>( i ), 1 ) = rhsZ[ i ];
		}

		Eigen::MatrixXd householder;
		double const eigenQr = bestOf( repeats, [ & ]()
		{
			Eigen::HouseholderQR<Eigen::MatrixXd> qr( A );
			householder = qr.solve( B );
		} );

		Eigen::MatrixXd complete;
		double eigenCod = 0.0;
		double codThreshold = 0.0;
		{
			Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
			eigenCod = bestOf( repeats, [ & ]()
			{
				cod.compute( A );
				cod.setThreshold( 1.0e-10 );
				complete = cod.solve( B );
			} );
			codThreshold = cod.threshold();
		}

		Eigen::MatrixXd bdc;
		double const eigenSvd = bestOf( repeats, [ & ]()
		{
			Eigen::BDCSVD<Eigen::MatrixXd> svd(
				A, Eigen::ComputeThinU | Eigen::ComputeThinV );
			svd.setThreshold( 1.0e-10 );
			bdc = svd.solve( B );
		} );

		auto worstAgainstMeq = [ & ]( Eigen::MatrixXd const &solution )
		{
			double worst = 0.0;
			double scale = 0.0;
			for ( std::size_t j = 0; j < cols; ++j )
			{
				worst = std::max(
					worst,
					std::max( std::fabs( solution( static_cast<Eigen::Index>( j ),
					                               0 ) - cR[ j ] ),
					          std::fabs( solution( static_cast<Eigen::Index>( j ),
					                               1 ) - cZ[ j ] ) ) );
				scale = std::max( scale, std::max( std::fabs( cR[ j ] ),
				                                   std::fabs( cZ[ j ] ) ) );
			}
			return worst/std::max( 1.0e-300, scale );
		};

		std::printf( "      Eigen HouseholderQR                     %10.5f s"
		             "   coefficients differ by %.3e (relative)\n",
		             eigenQr, worstAgainstMeq( householder ) );
		std::printf( "      Eigen CompleteOrthogonalDecomposition   %10.5f s"
		             "   %.3e, threshold %.1e, rank deficiency %zu\n",
		             eigenCod, worstAgainstMeq( complete ), codThreshold,
		             cols - static_cast<std::size_t>(
		                 Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>(
		                     A ).rank() ) );
		std::printf( "      Eigen BDCSVD                            %10.5f s"
		             "   %.3e\n", eigenSvd, worstAgainstMeq( bdc ) );
		std::printf( "      meq's own QR + one-sided Jacobi is not separable"
		             " from its constructor; see section 1 for the difference\n" );

		// The coefficients are a means, not the deliverable. What has to survive
		// a change of decomposition is the POSITION, so that is what is checked:
		// a rank-revealing solve may put a different vector in the soft tail and
		// still describe the same surface, which is exactly the shape IN-4
		// measured -- "a soft tail with no gap".
		double worstPosition = 0.0;
		for ( std::size_t i = 0; i < rows; ++i )
		{
			double const rho =
				std::sqrt( samples[ i ].normalisedFlux/options.discEdge );
			table.evaluate( rho, samples[ i ].theta );
			std::vector<double> const &v = table.values();
			double r = 0.0, z = 0.0, rq = 0.0, zq = 0.0;
			for ( std::size_t j = 0; j < cols; ++j )
			{
				r += v[ j ]*cR[ j ];
				z += v[ j ]*cZ[ j ];
				rq += v[ j ]*householder( static_cast<Eigen::Index>( j ), 0 );
				zq += v[ j ]*householder( static_cast<Eigen::Index>( j ), 1 );
			}
			worstPosition = std::max( worstPosition,
			                          std::hypot( r - rq, z - zq ) );
		}
		std::printf( "      the two fits as CURVES differ by %.3e m at the"
		             " sample points\n", worstPosition );
		verdict.check( worstPosition < 1.0e-9,
		               "Eigen's least-squares solution is a different curve from"
		               " meq's, not merely a different coefficient vector: the"
		               " swap would move IN-3's and IN-4's numbers" );
#else
		std::printf( "      Eigen is not available in this build, so the"
		             " comparison is absent rather than assumed\n" );
		(void)rhsR;
		(void)rhsZ;
#endif

		// -------------------------------------------------------------------
		// K4. Field evaluation in the tracer: locate against evaluate, and what
		//     batching by element is worth.
		// -------------------------------------------------------------------
		std::vector<double> pointR;
		std::vector<double> pointZ;
		for ( meq::AngleParametrisation const &fit : family.fits )
			for ( std::size_t j = 0; j < fit.count(); ++j )
			{
				pointR.push_back( fit.pointR[ j ] );
				pointZ.push_back( fit.pointZ[ j ] );
			}

		std::size_t const cloud = pointR.size();

		double const seamCost = bestOf( repeats, [ & ]()
		{
			int hint = -1;
			for ( std::size_t i = 0; i < cloud; ++i )
			{
				double psi = 0.0, qR = 0.0, qZ = 0.0;
				tracer.sampleAt( pointR[ i ], pointZ[ i ], psi, qR, qZ, hint );
			}
		} );

		// The same points, located once, so that the location and the evaluation
		// can be charged separately.
		//
		// TWO PASSES, AND THE FIRST IS NOT TIMED. Pass one resolves every point
		// exhaustively so that pass two can be handed the hint the tracer would
		// have -- the element the PREVIOUS point was located in, however it was
		// located. Without that a single miss poisons every hint after it and
		// the walk reports a hundred per cent miss rate, which is a measurement
		// of the seed and not of the walk. Measured: it does.
		//
		// AND THE WALK IS REPLICATED FAITHFULLY RATHER THAN APPROXIMATED.
		// ContourTracer::locate() is private, so this is a copy of it:
		// TransformBack into the hint, then a breadth-first widening over face
		// neighbours to four rings, then give up. A LINEAR SCAN OVER THE
		// ELEMENTS WAS WRITTEN FIRST AND IS A DIFFERENT MEASUREMENT -- it read
		// 1.20 s against sampleAt()'s own 0.054 s on the same points, twenty-two
		// times the cost of the thing it was supposed to be a PART of. That is
		// the cost of the FALLBACK and not of the walk, and the two look alike
		// in a profile while differing by three orders.
		mfem::Mesh &mesh = solved.theMesh();
		mfem::Table const &neighbours = mesh.ElementToElementTable();
		std::vector<int> element( cloud, -1 );
		std::vector<mfem::IntegrationPoint> reference( cloud );
		long misses = 0;

		{
			mfem::IsoparametricTransformation scratch;
			mfem::Vector point( 2 );
			mfem::IntegrationPoint ip;
			for ( std::size_t i = 0; i < cloud; ++i )
			{
				point( 0 ) = pointR[ i ];
				point( 1 ) = pointZ[ i ];
				for ( int e = 0; e < mesh.GetNE(); ++e )
				{
					mesh.GetElementTransformation( e, &scratch );
					if ( scratch.TransformBack( point, ip )
					     == mfem::InverseElementTransformation::Inside )
					{
						element[ i ] = e;
						reference[ i ] = ip;
						break;
					}
				}
			}
		}

		double const locateCost = bestOf( repeats, [ & ]()
		{
			mfem::IsoparametricTransformation scratch;
			mfem::Vector point( 2 );
			mfem::IntegrationPoint ip;
			std::vector<int> visited;
			misses = 0;

			auto tryElement = [ & ]( int candidate ) -> bool
			{
				if ( candidate < 0 || candidate >= mesh.GetNE() )
					return false;
				mesh.GetElementTransformation( candidate, &scratch );
				return scratch.TransformBack( point, ip )
				       == mfem::InverseElementTransformation::Inside;
			};

			for ( std::size_t i = 0; i < cloud; ++i )
			{
				point( 0 ) = pointR[ i ];
				point( 1 ) = pointZ[ i ];
				int const hint = ( i == 0 ) ? -1 : element[ i - 1 ];

				if ( hint >= 0 && tryElement( hint ) )
					continue;

				bool found = false;
				if ( hint >= 0 )
				{
					visited.clear();
					visited.push_back( hint );
					std::size_t frontier = 0;
					for ( int depth = 0; depth < 4 && !found; ++depth )
					{
						std::size_t const end = visited.size();
						if ( frontier >= end )
							break;
						for ( std::size_t k = frontier; k < end && !found; ++k )
						{
							int const *row = neighbours.GetRow( visited[ k ] );
							int const size = neighbours.RowSize( visited[ k ] );
							for ( int j = 0; j < size && !found; ++j )
							{
								int const candidate = row[ j ];
								if ( candidate < 0
								     || std::find( visited.begin(), visited.end(),
								                   candidate ) != visited.end() )
									continue;
								visited.push_back( candidate );
								found = tryElement( candidate );
							}
						}
						frontier = end;
					}
				}

				// A MISS IS COUNTED AND NOT PAID FOR HERE, so that this column
				// stays a measurement of the WALK. The tracer hands a miss to
				// Mesh::FindPoints, which is O( elements ) and would swamp
				// everything else: charging it here would report the fallback
				// and call it the walk.
				if ( !found )
					++misses;
			}
		} );

		mfem::GridFunction const &potential =
			solved.theSolver().postProcessedPotential();
		mfem::GridFunction const &flux = solved.theSolver().postProcessedFlux();

		std::vector<double> perPointPsi( cloud, 0.0 );
		double const evaluateCost = bestOf( repeats, [ & ]()
		{
			mfem::Vector q( 2 );
			for ( std::size_t i = 0; i < cloud; ++i )
			{
				perPointPsi[ i ] = potential.GetValue( element[ i ],
				                                       reference[ i ] );
				flux.GetVectorValue( element[ i ], reference[ i ], q );
			}
		} );

		// Batched: the points grouped by element, one GetValues() per element.
		// MFEM'S BATCHED CALL IS A LOOP AND NOT A GEMM, which is worth saying
		// because it is the opposite of what one would guess from the name:
		// GridFunction::GetValues( i, ir, vals ) calls CalcShape() once per
		// point and dots it with the local data. What it amortises is the dof
		// gather and the temporaries, not the arithmetic.
		std::vector<std::vector<std::size_t>> byElement( mesh.GetNE() );
		for ( std::size_t i = 0; i < cloud; ++i )
			if ( element[ i ] >= 0 )
				byElement[ element[ i ] ].push_back( i );

		std::size_t occupied = 0;
		std::size_t deepest = 0;
		for ( std::vector<std::size_t> const &list : byElement )
		{
			if ( !list.empty() )
				++occupied;
			deepest = std::max( deepest, list.size() );
		}

		std::vector<double> batchedPsi( cloud, 0.0 );
		double const batchedCost = bestOf( repeats, [ & ]()
		{
			mfem::Vector values;
			mfem::DenseMatrix vectors;
			mfem::IsoparametricTransformation scratch;
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				std::vector<std::size_t> const &list = byElement[ e ];
				if ( list.empty() )
					continue;

				mfem::IntegrationRule rule( static_cast<int>( list.size() ) );
				for ( std::size_t k = 0; k < list.size(); ++k )
					rule.IntPoint( static_cast<int>( k ) ) = reference[ list[ k ] ];

				// The ElementTransformation overload for the flux, and OUR OWN
				// transformation: the ( int, ir, vals, tr ) one reaches for the
				// mesh's shared scratch and computes the physical coordinates of
				// every point, which nothing here wants.
				mesh.GetElementTransformation( e, &scratch );
				potential.GetValues( e, rule, values );
				flux.GetVectorValues( scratch, rule, vectors );
				for ( std::size_t k = 0; k < list.size(); ++k )
					batchedPsi[ list[ k ] ] = values( static_cast<int>( k ) );
			}
		} );

		double worstPsi = 0.0;
		for ( std::size_t i = 0; i < cloud; ++i )
			worstPsi = std::max( worstPsi,
			                     std::fabs( perPointPsi[ i ] - batchedPsi[ i ] ) );

		std::printf( "\n  K4  the field at %zu points -- the tracer's own"
		             " primitive, taken apart\n", cloud );
		std::printf( "      ContourTracer::sampleAt(), hint threaded  %10.5f s\n",
		             seamCost );
		std::printf( "      the WALK alone, misses not paid for       %10.5f s"
		             "   %5.1f%% of it\n", locateCost,
		             100.0*locateCost/seamCost );
		std::printf( "      evaluating alone, per point               %10.5f s"
		             "   %5.1f%%\n", evaluateCost,
		             100.0*evaluateCost/seamCost );
		std::printf( "      the REST is Mesh::FindPoints             %10.5f s"
		             "   %5.1f%%, on %ld walk misses of %zu points\n",
		             seamCost - locateCost - evaluateCost,
		             100.0*( seamCost - locateCost - evaluateCost )/seamCost,
		             misses, cloud );
		std::printf( "      evaluating batched by element             %10.5f s"
		             "   %6.2fx over per point\n", batchedCost,
		             evaluateCost/batchedCost );
		std::printf( "      %zu points over %zu elements, %.2f per element,"
		             " deepest %zu\n", cloud, occupied,
		             static_cast<double>( cloud )/occupied, deepest );
		std::printf( "      batched vs per point, worst psi difference %.3e"
		             "\n", worstPsi );
		verdict.check( worstPsi == 0.0,
		               "GetValues() over an IntegrationRule does NOT reproduce"
		               " GetValue() point by point, which it must -- it is the"
		               " same arithmetic in the same order" );

		// -------------------------------------------------------------------
		// K5. setWalkDepth(), which is not a library kernel at all but is where
		//     K4 says the time actually is. Four rings of FACE neighbours reach
		//     about four triangles in a straight line and fewer than that
		//     diagonally, and consecutive ray nodes are one to two cells apart
		//     -- so the walk misses, and every miss is an O( elements ) scan.
		//     The knob already exists and nothing had measured it.
		// -------------------------------------------------------------------
		std::printf( "\n  K5  ContourTracer::setWalkDepth() -- the same rays,"
		             " re-fitted at each depth\n" );
		std::printf( "      %6s %11s %14s %12s\n", "depth", "seconds",
		             "fallbacks", "worst |dx|" );

		std::vector<meq::AngleParametrisation> reference4;
		for ( int depth : { 2, 4, 6, 8, 12, 16 } )
		{
			tracer.setWalkDepth( depth );
			std::vector<meq::AngleParametrisation> fits( family.contours.size() );
			double const cost = bestOf( repeats, [ & ]()
			{
				for ( std::size_t i = 0; i < family.contours.size(); ++i )
					fits[ i ] = tracer.fitByAngle(
						family.contours[ i ], axis,
						static_cast<std::size_t>( spec.angles ) );
			} );

			long fallbacks = 0;
			double worst = 0.0;
			for ( std::size_t i = 0; i < fits.size(); ++i )
			{
				fallbacks += fits[ i ].fallbackLocations;
				for ( std::size_t j = 0; j < fits[ i ].count(); ++j )
					worst = std::max( worst,
						std::fabs( fits[ i ].radius[ j ]
						           - family.fits[ i ].radius[ j ] ) );
			}

			std::printf( "      %6d %11.5f %14ld %12.3e%s\n", depth, cost,
			             fallbacks, worst,
			             worst == 0.0 ? "   (same answer)" : "" );

			// THE DEPTH MUST NOT CHANGE THE ANSWER, and that is the property
			// that makes raising it free rather than a trade. FluxSurfaces.hpp
			// says the fallback "returns the same element", so a deeper walk
			// finds by neighbour what a shallower one found by scan -- and the
			// located element, the reference point and therefore every field
			// value are identical. A difference here would mean the walk and
			// FindPoints disagree about which element a point is in, which on a
			// discontinuous field is two different answers.
			verdict.check( worst == 0.0,
			               "changing setWalkDepth() changed the fitted radii."
			               " The walk and Mesh::FindPoints must locate the same"
			               " element, so this is a real disagreement and not a"
			               " tuning question" );
		}
		tracer.setWalkDepth( 4 );
	}

}

// ===========================================================================
// 4. THREADS
// ===========================================================================

namespace
{

	/// A ray solve on the tracer's own seam: safeguarded Newton on rho along the
	/// ray theta from the axis, with the derivative from q.
	///
	/// A HARNESS REIMPLEMENTATION AND LABELLED AS ONE. ContourTracer::
	/// fitByAngle() runs the whole ray loop inside one call, so its rays cannot
	/// be threaded from outside without a library change. This is the same
	/// iteration on the same primitive, and it is here to answer the question
	/// INVERSION-PLAN.md section 11.2 item 3 asks -- what the rays are worth if
	/// they were threaded -- rather than to replace anything.
	double solveRay( meq::ContourTracer const &tracer, double axisR, double axisZ,
	                 double theta, double level, double seed, int &iterations )
	{
		double const uR = std::cos( theta );
		double const uZ = std::sin( theta );
		double rho = seed;
		iterations = 0;

		for ( int i = 0; i < 40; ++i )
		{
			double psi = 0.0, qR = 0.0, qZ = 0.0;
			int hint = -1;
			if ( !tracer.sampleAt( axisR + rho*uR, axisZ + rho*uZ, psi, qR, qZ,
			                       hint ) )
				return rho;

			++iterations;
			double const r = axisR + rho*uR;

			// d psi / d rho = grad psi . u = r q . u.
			double const slope = r*( qR*uR + qZ*uZ );
			if ( !( std::fabs( slope ) > 0.0 ) )
				return rho;

			double const stepSize = ( level - psi )/slope;
			rho += stepSize;
			if ( std::fabs( stepSize ) < 1.0e-13 )
				break;
		}

		return rho;
	}

	/*
	 * THREADING, AND THE FIRST THING TO SAY IS THAT THE CHAIN CANNOT BE THREADED
	 * AS IT STANDS. That is a measurement and not a caution: the first version
	 * of this section put one shared const ContourTracer in an
	 * omp parallel for over surfaces, which is exactly what INVERSION-PLAN.md
	 * section 11.2 item 1 proposes, and it ABORTED --
	 *
	 *     terminate called after throwing an instance of 'std::runtime_error'
	 *       what():  ContourTracer::traceFromAxis: the axis is not in the mesh
	 *
	 * -- because every entry point reaches mfem::Mesh::FindPoints, and
	 * FindPoints is not reentrant. It loops over every element calling the
	 * SHARED GetElementTransformation( int ) overload, and it builds a
	 * vertex-to-element table on the way. Section 11.3 anticipated the second
	 * half of this and got the premise wrong: it says the fallback "must stay
	 * outside, or be serialised. It is already O( elements x points ) and the
	 * tracer reports zero fallbacks, so this costs nothing to honour." The
	 * tracer reports zero fallbacks ON A TRACE. traceFromAxis() begins by
	 * sampling the axis with no hint at all, so it takes the fallback ONCE PER
	 * SURFACE unconditionally, and section 3's K5 measures fitByAngle() taking
	 * it 183 times in 576 rays at the default walk depth.
	 *
	 * So what is measured here is the parallelism that is AVAILABLE, by giving
	 * each thread its own mesh, its own solve and its own tracer -- which shares
	 * nothing at all, at the cost of one solve per thread outside the timed
	 * region. It answers section 11.2's question ("what are surfaces and rays
	 * worth threaded") without pretending the sharing is not there, and the
	 * bit-exactness assertion still does its job: with genuinely independent
	 * state the answers must be identical to the last bit, and anything left
	 * over -- a thread_local that is not, a lazily built global table -- shows up
	 * as a difference.
	 *
	 * THE SERIAL PASS RUNS FIRST ON PURPOSE. mfem::IntRules builds its rules
	 * lazily on first Get(), and Mesh::ElementToElementTable() caches on first
	 * call; both would be a data race if a parallel region met them cold. Doing
	 * the reference extraction first warms every one of them.
	 */
	void threads( Equilibrium const &eq, ExactAxis const &exact,
	              FamilySpec const &spec, int order, int n, Verdict &verdict )
	{
		int maximum = 1;
#ifdef _OPENMP
		maximum = omp_get_max_threads();
#endif

		std::printf( "\n=== 4. THREADS ===\n" );
		std::printf( "    omp_get_max_threads() = %d. OpenMP over INDEPENDENT"
		             " WORK and never threaded BLAS: MKL_NUM_THREADS stays at 1.\n",
		             maximum );
		std::printf( "    ONE MESH, ONE SOLVE AND ONE TRACER PER THREAD, because"
		             " mfem::Mesh::FindPoints is not reentrant and every\n"
		             "    entry point reaches it. See the comment on this"
		             " function; the shared-tracer version aborts.\n" );

		std::vector<std::unique_ptr<Solved>> perThread;
		std::vector<std::unique_ptr<meq::ContourTracer>> tracers;
		std::vector<meq::CriticalPoint> axes;

		double const buildCost = bestOf( 1, [ & ]()
		{
			for ( int t = 0; t < maximum; ++t )
			{
				perThread.push_back(
					std::unique_ptr<Solved>( new Solved( eq, nstxBox(), order, n ) ) );
				meq::CriticalPointFinder finder( perThread.back()->theSolver() );
				axes.push_back( finder.findAxis() );
				tracers.push_back( std::unique_ptr<meq::ContourTracer>(
					new meq::ContourTracer( perThread.back()->theSolver() ) ) );
			}
		} );

		// The per-thread solves must agree BIT FOR BIT, or the extraction below
		// cannot and the failure would be blamed on the threading. Checked
		// rather than assumed: a direct solve is deterministic at
		// MKL_NUM_THREADS=1 and this says so.
		double worstField = 0.0;
		for ( int t = 1; t < maximum; ++t )
		{
			mfem::GridFunction const &a =
				perThread[ 0 ]->theSolver().postProcessedPotential();
			mfem::GridFunction const &b =
				perThread[ t ]->theSolver().postProcessedPotential();
			for ( int i = 0; i < a.Size(); ++i )
				worstField = std::max( worstField, std::fabs( a( i ) - b( i ) ) );
		}
		std::printf( "    %d independent solves in %.2f s; they agree to %.3e"
		             "%s\n", maximum, buildCost, worstField,
		             worstField == 0.0 ? " (bit for bit)" : "" );
		verdict.check( worstField == 0.0,
		               "two independent solves of the same problem in one process"
		               " differ. Everything below compares extractions of two"
		               " different fields and means nothing" );

		std::vector<double> const ladder = fluxLadder( spec );
		std::size_t const angles = static_cast<std::size_t>( spec.angles );
		std::size_t const count = ladder.size();

		std::vector<meq::Contour> serialContours( count );
		std::vector<meq::AngleParametrisation> serialFits( count );
		std::vector<meq::SurfaceAverages> serialAverages( count );
		double const serialCost = bestOf( 2, [ & ]()
		{
			for ( std::size_t i = 0; i < count; ++i )
				extractOne( *tracers[ 0 ], axes[ 0 ], exact, ladder[ i ], angles,
				            serialContours[ i ], serialFits[ i ],
				            serialAverages[ i ] );
		} );

		std::printf( "\n  (a) over SURFACES -- trace, rays and quadrature, %zu"
		             " independent surfaces\n", count );
		std::printf( "    %8s %11s %9s %16s\n", "threads", "seconds", "speedup",
		             "worst difference" );
		std::printf( "    %8d %11.5f %9s %16s\n", 1, serialCost, "1.00x",
		             "(the reference)" );

		std::vector<int> counts;
		for ( int t = 2; t <= maximum; t *= 2 )
			counts.push_back( t );
		if ( maximum > 1 && ( counts.empty() || counts.back() != maximum ) )
			counts.push_back( maximum );

		for ( int t : counts )
		{
			std::vector<meq::Contour> contours( count );
			std::vector<meq::AngleParametrisation> fits( count );
			std::vector<meq::SurfaceAverages> averages( count );

			double const cost = bestOf( 2, [ & ]()
			{
#ifdef _OPENMP
				#pragma omp parallel for num_threads( t ) schedule( dynamic )
#endif
				for ( int i = 0; i < static_cast<int>( count ); ++i )
				{
					int which = 0;
#ifdef _OPENMP
					which = omp_get_thread_num();
#endif
					extractOne( *tracers[ which ], axes[ which ], exact,
					            ladder[ i ], angles, contours[ i ], fits[ i ],
					            averages[ i ] );
				}
			} );

			// EXACT, at 0.000e+00 and not at a tolerance. Independent surfaces
			// over independent state reassociate nothing, so exactness is
			// available -- and a tolerance here would be an admission that
			// something IS shared, which is precisely the defect being guarded
			// against.
			double worst = 0.0;
			for ( std::size_t i = 0; i < count; ++i )
			{
				if ( contours[ i ].points.size()
				     != serialContours[ i ].points.size()
				     || fits[ i ].count() != serialFits[ i ].count()
				     || averages[ i ].nodes.size()
				        != serialAverages[ i ].nodes.size() )
				{
					worst = 1.0e300;
					continue;
				}

				for ( std::size_t j = 0; j < contours[ i ].points.size(); ++j )
				{
					worst = std::max( worst,
						std::fabs( contours[ i ].points[ j ].r
						           - serialContours[ i ].points[ j ].r ) );
					worst = std::max( worst,
						std::fabs( contours[ i ].points[ j ].z
						           - serialContours[ i ].points[ j ].z ) );
				}
				for ( std::size_t j = 0; j < fits[ i ].count(); ++j )
					worst = std::max( worst,
						std::fabs( fits[ i ].radius[ j ]
						           - serialFits[ i ].radius[ j ] ) );
				for ( std::size_t j = 0; j < averages[ i ].nodes.size(); ++j )
					worst = std::max( worst,
						std::fabs( averages[ i ].nodes[ j ].weight
						           - serialAverages[ i ].nodes[ j ].weight ) );
				worst = std::max( worst, std::fabs( averages[ i ].vPrime
				                                    - serialAverages[ i ].vPrime ) );
			}

			std::printf( "    %8d %11.5f %8.2fx %16.3e%s\n", t, cost,
			             serialCost/cost, worst,
			             worst == 0.0 ? "   (bit for bit)" : "" );
			std::fflush( stdout );

			verdict.check( worst == 0.0,
			               "a threaded extraction over surfaces did NOT reproduce"
			               " the serial one bit for bit, on state that is"
			               " per-thread. Something else is shared" );
		}

		// -------------------------------------------------------------------
		// (b) over RAYS.
		// -------------------------------------------------------------------
		double const level = levelAt( exact, 0.35 );
		std::size_t const rayCount = 4096;
		std::vector<double> thetas( rayCount, 0.0 );
		for ( std::size_t i = 0; i < rayCount; ++i )
			thetas[ i ] = twoPi*static_cast<double>( i )/rayCount;

		std::vector<double> serialRadius( rayCount, 0.0 );
		double const raySerial = bestOf( 2, [ & ]()
		{
			for ( std::size_t i = 0; i < rayCount; ++i )
			{
				int its = 0;
				serialRadius[ i ] = solveRay( *tracers[ 0 ], axes[ 0 ].r,
				                              axes[ 0 ].z, thetas[ i ], level,
				                              0.25, its );
			}
		} );

		std::printf( "\n  (b) over RAYS -- %zu independent 1-D solves on"
		             " ContourTracer::sampleAt()\n", rayCount );
		std::printf( "    %8s %11s %9s %16s\n", "threads", "seconds", "speedup",
		             "worst difference" );
		std::printf( "    %8d %11.5f %9s %16s\n", 1, raySerial, "1.00x",
		             "(the reference)" );

		for ( int t : counts )
		{
			std::vector<double> radius( rayCount, 0.0 );
			double const cost = bestOf( 2, [ & ]()
			{
#ifdef _OPENMP
				#pragma omp parallel for num_threads( t ) schedule( static )
#endif
				for ( int i = 0; i < static_cast<int>( rayCount ); ++i )
				{
					int which = 0;
#ifdef _OPENMP
					which = omp_get_thread_num();
#endif
					int its = 0;
					radius[ i ] = solveRay( *tracers[ which ], axes[ which ].r,
					                        axes[ which ].z, thetas[ i ], level,
					                        0.25, its );
				}
			} );

			double worst = 0.0;
			for ( std::size_t i = 0; i < rayCount; ++i )
				worst = std::max( worst, std::fabs( radius[ i ]
				                                    - serialRadius[ i ] ) );

			std::printf( "    %8d %11.5f %8.2fx %16.3e%s\n", t, cost,
			             raySerial/cost, worst,
			             worst == 0.0 ? "   (bit for bit)" : "" );
			std::fflush( stdout );

			verdict.check( worst == 0.0,
			               "a threaded ray population did NOT reproduce the"
			               " serial one bit for bit" );
		}
	}

}

// ===========================================================================
// 5. THE TWO ALGORITHMIC LEVERS, AND THE CONSUMER'S CALL PATTERN
// ===========================================================================

namespace
{

	void levers( Solved &solved, ExactAxis const &exact, FamilySpec const &spec,
	             int repeats )
	{
		meq::CriticalPointFinder finder( solved.theSolver() );
		meq::CriticalPoint const axis = finder.findAxis();
		meq::ContourTracer tracer( solved.theSolver() );

		std::vector<double> const ladder = fluxLadder( spec );
		std::size_t const count = ladder.size();
		std::size_t const angles = static_cast<std::size_t>( spec.angles );

		std::printf( "\n=== 5. THE ALGORITHMIC LEVERS ===\n" );

		// -------------------------------------------------------------------
		// (a) Continuation in the flux label.
		// -------------------------------------------------------------------
		std::vector<meq::Contour> fresh( count );
		double const freshCost = bestOf( repeats, [ & ]()
		{
			for ( std::size_t i = 0; i < count; ++i )
				fresh[ i ] = tracer.traceFromAxis( levelAt( exact, ladder[ i ] ),
				                                   axis );
		} );

		std::vector<meq::Contour> continued( count );
		double const continuedCost = bestOf( repeats, [ & ]()
		{
			double startR = axis.r;
			double startZ = axis.z;
			for ( std::size_t i = 0; i < count; ++i )
			{
				if ( i == 0 )
					continued[ i ] = tracer.traceFromAxis(
						levelAt( exact, ladder[ i ] ), axis );
				else
					continued[ i ] = tracer.trace( levelAt( exact, ladder[ i ] ),
					                               startR, startZ );
				startR = continued[ i ].points.front().r;
				startZ = continued[ i ].points.front().z;
			}
		} );

		long freshIterations = 0;
		long continuedIterations = 0;
		double worstAgreement = 0.0;
		for ( std::size_t i = 0; i < count; ++i )
		{
			freshIterations += fresh[ i ].correctorIterationsTotal;
			continuedIterations += continued[ i ].correctorIterationsTotal;
			worstAgreement = std::max(
				worstAgreement,
				std::fabs( fresh[ i ].hermiteLength()
				           - continued[ i ].hermiteLength() ) );
		}

		std::printf( "\n  (a) continuation in the flux label, %zu surfaces\n",
		             count );
		std::printf( "    %-40s %11s %14s\n", "", "seconds",
		             "corrector its" );
		std::printf( "    %-40s %11.5f %14ld\n", "traceFromAxis, each from"
		             " scratch", freshCost, freshIterations );
		std::printf( "    %-40s %11.5f %14ld\n", "trace() from the previous"
		             " surface", continuedCost, continuedIterations );
		std::printf( "    ratio %.3fx, and the two families agree in arc length"
		             " to %.3e\n", freshCost/continuedCost, worstAgreement );
		std::printf( "    READ IT AGAINST THE STRUCTURE: the predictor for every"
		             " point after the first ALREADY comes from the previous\n"
		             "    point of the SAME surface, so continuation can only"
		             " save the SEEDING -- traceFromAxis()'s ray bracket\n"
		             "    and the first corrector. It cannot make the trace"
		             " itself cheaper, and it makes the surfaces sequential.\n" );

		// -------------------------------------------------------------------
		// (b) The MaNTA call pattern, and what a per-psi cache is worth.
		// -------------------------------------------------------------------
		// Sixty nodes against twelve surfaces, which is the shape
		// MANTA-COUPLING.md section 5 describes: the physics grid is finer
		// than any surface family a geometry cache would hold, so the ratio
		// the cache buys is nodes/surfaces and the measurement is whether
		// anything eats it.
		int const nodes = 60;
		int const residuals = 3;

		std::printf( "\n  (b) MANTA-COUPLING.md section 5's pointwise pattern:"
		             " %d physics nodes, %d residual evaluations\n",
		             nodes, residuals );

		double const naive = bestOf( 1, [ & ]()
		{
			for ( int residual = 0; residual < residuals; ++residual )
				for ( int node = 0; node < nodes; ++node )
				{
					double const psiN = spec.smallest
						+ ( spec.largest - spec.smallest )*node/( nodes - 1.0 );
					meq::Contour contour;
					meq::AngleParametrisation fit;
					meq::SurfaceAverages average;
					extractOne( tracer, axis, exact, psiN, angles, contour, fit,
					            average );
				}
		} );

		double const cached = bestOf( 1, [ & ]()
		{
			for ( int residual = 0; residual < residuals; ++residual )
			{
				Family const family = extractFamily( tracer, axis, exact, spec );
				for ( int node = 0; node < nodes; ++node )
				{
					// The cache hit: the family is already there and the node
					// reads it. Interpolation between surfaces is IN-6's work
					// and is a handful of flops whatever it turns out to be, so
					// what is timed here is the MISS rate and nothing else.
					double const psiN = spec.smallest
						+ ( spec.largest - spec.smallest )*node/( nodes - 1.0 );
					double closest = 1.0e300;
					std::size_t best = 0;
					for ( std::size_t i = 0; i < family.normalisedFlux.size(); ++i )
						if ( std::fabs( family.normalisedFlux[ i ] - psiN )
						     < closest )
						{
							closest = std::fabs( family.normalisedFlux[ i ] - psiN );
							best = i;
						}
					double const sink = family.averages[ best ].vPrime;
					if ( sink == 1.0e300 )
						std::printf( " " );
				}
			}
		} );

		std::printf( "    %-46s %11.5f s\n", "one extraction per node per"
		             " residual", naive );
		std::printf( "    %-46s %11.5f s   %6.1fx\n", "one FAMILY per residual,"
		             " nodes read it", cached, naive/cached );
		std::printf( "    misses: %d against %d, i.e. a hit rate of %.1f%%\n",
		             residuals*spec.surfaces, residuals*nodes,
		             100.0*( 1.0 - static_cast<double>( spec.surfaces )/nodes ) );

		// -------------------------------------------------------------------
		// (c) dGeometry_dpsi by differencing, which is the number to design
		//     against.
		// -------------------------------------------------------------------
		int const fieldDofs =
			solved.theSolver().postProcessedPotential().Size();
		double const oneFamily = bestOf( repeats, [ & ]()
		{
			Family const family = extractFamily( tracer, axis, exact, spec );
			(void)family;
		} );

		std::printf( "\n  (c) dGeometry_dpsi by differencing, shaped"
		             " ( nGeometry, nFieldDOF )\n" );
		std::printf( "    one family extraction   %11.5f s\n", oneFamily );
		std::printf( "    nFieldDOF               %11d   (psi* on this mesh)\n",
		             fieldDofs );
		std::printf( "    differenced Jacobian    %11.1f s   = %.2f hours,"
		             " ONE evaluation, serial\n", oneFamily*fieldDofs,
		             oneFamily*fieldDofs/3600.0 );
		std::printf( "    it is embarrassingly parallel over psi DOFs, so a"
		             " core count divides it -- and nothing else here does.\n"
		             "    INVERSION-PLAN.md section 11.4's shape derivative is"
		             " the lever that changes the exponent rather than the\n"
		             "    constant, and this is the number it has to beat.\n" );
	}

}

int main( int argc, char **argv )
{
	int repeats = 2;
	int order = 2;
	int n = 48;
	int degree = 10;
	FamilySpec spec;
	std::string sections = "1,2,3,4,5";

	auto parseList = []( char const *text )
	{
		std::vector<int> out;
		std::string item;
		for ( char const *c = text; ; ++c )
		{
			if ( *c == ',' || *c == '\0' )
			{
				if ( !item.empty() )
					out.push_back( std::atoi( item.c_str() ) );
				item.clear();
				if ( *c == '\0' )
					break;
			}
			else
				item.push_back( *c );
		}
		return out;
	};

	std::vector<int> orders = { 1, 2, 3 };
	std::vector<int> sizes = { 24, 48, 96 };

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--repeats" && i + 1 < argc )
			repeats = std::atoi( argv[ ++i ] );
		else if ( arg == "--order" && i + 1 < argc )
			order = std::atoi( argv[ ++i ] );
		else if ( arg == "--size" && i + 1 < argc )
			n = std::atoi( argv[ ++i ] );
		else if ( arg == "--degree" && i + 1 < argc )
			degree = std::atoi( argv[ ++i ] );
		else if ( arg == "--surfaces" && i + 1 < argc )
			spec.surfaces = std::atoi( argv[ ++i ] );
		else if ( arg == "--angles" && i + 1 < argc )
			spec.angles = std::atoi( argv[ ++i ] );
		else if ( arg == "--orders" && i + 1 < argc )
			orders = parseList( argv[ ++i ] );
		else if ( arg == "--sizes" && i + 1 < argc )
			sizes = parseList( argv[ ++i ] );
		else if ( arg == "--sections" && i + 1 < argc )
			sections = argv[ ++i ];
	}

	auto wanted = [ &sections ]( char which )
	{
		return sections.find( which ) != std::string::npos;
	};

	char const *ompEnv = std::getenv( "OMP_NUM_THREADS" );
	char const *mklEnv = std::getenv( "MKL_NUM_THREADS" );

	std::printf( "\n=== meq solution-inversion scaling (IN-P) ===\n" );
	std::printf( "  OMP_NUM_THREADS=%-6s MKL_NUM_THREADS=%-6s best of %d\n",
	             ompEnv ? ompEnv : "(unset)", mklEnv ? mklEnv : "(unset)",
	             repeats );
	if ( mklEnv == nullptr || std::string( mklEnv ) != "1" )
		std::printf( "  NOTE: MKL_NUM_THREADS is not 1. CLAUDE.md records"
		             " ComputeH()'s element-local dense LU degrading by a factor"
		             " of\n        forty at k = 3 under threaded MKL, so the"
		             " solve below is not comparable with any other run.\n" );
	std::printf( "  build:" );
#ifdef _OPENMP
	std::printf( " OPENMP" );
#endif
#ifdef MFEM_THREAD_SAFE
	std::printf( " MFEM_THREAD_SAFE" );
#endif
#ifdef MEQ_HAVE_EIGEN
	std::printf( " EIGEN" );
#endif
	std::printf( "\n" );

#if defined( _OPENMP ) && !defined( MFEM_THREAD_SAFE )
	std::printf( "  *** MFEM_THREAD_SAFE is OFF. FiniteElement evaluation is"
	             " then not reentrant and section 4 is measuring a race.\n" );
#endif

	Equilibrium const eq = Equilibrium::nstx();
	ExactAxis const exact = exactAxis( eq, 1.3, 0.0 );

	Verdict verdict;

	if ( wanted( '1' ) || wanted( '3' ) || wanted( '4' ) || wanted( '5' ) )
	{
		Solved solved( eq, nstxBox(), order, n );
		std::printf( "\n  the equilibrium: nstx() Solov'ev on"
		             " [%.2f,%.2f]x[%.2f,%.2f], k = %d, n = %d, %d elements,"
		             " %d psi* dofs\n",
		             nstxBox().rMin, nstxBox().rMax, nstxBox().zMin,
		             nstxBox().zMax, order, n, solved.theMesh().GetNE(),
		             solved.theSolver().postProcessedPotential().Size() );

		if ( wanted( '1' ) )
			breakdown( solved, exact, spec, degree, repeats, verdict );
		if ( wanted( '3' ) )
			kernels( solved, exact, spec, degree, repeats, verdict );
		if ( wanted( '4' ) )
			threads( eq, exact, spec, order, n, verdict );
		if ( wanted( '5' ) )
			levers( solved, exact, spec, repeats );
	}

	if ( wanted( '2' ) )
		scaling( eq, exact, orders, sizes, spec, degree, repeats );

	std::printf( "\n  correctness, which is what makes the timings mean"
	             " anything\n" );
	if ( verdict.failures == 0 )
		std::printf( "    every property held: threaded == serial bit for bit,"
		             " and every candidate kernel reproduces what it would\n"
		             "    replace. Mesh::FindPoints fallbacks are printed above"
		             " and are NOT one of these, being a cost and not an\n"
		             "    answer -- read the count there before reading the"
		             " tables.\n\n" );
	else
		std::printf( "\n    *** %d correctness propert%s failed above. The"
		             " timings are not a result until that is fixed.\n\n",
		             verdict.failures,
		             verdict.failures == 1 ? "y" : "ies" );

	return verdict.failures == 0 ? 0 : 1;
}
