#include "meq/SafetyFactorSolve.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

// EIGEN FOR THE NEWTON SYSTEM, AND IT IS HEADER ONLY, so this file's consumers
// take on nothing -- the same arrangement SurfaceFit.cpp records. It is here
// because the Jacobian is a dense n x n that can be RANK DEFICIENT, and a
// rank-revealing decomposition is the one thing that distinguishes "the step is
// zero" from "the step is infinite". See DenseJacobianSolver.
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "mfem.hpp"

namespace meq
{

	namespace
	{

		/**
		 * `sigma_min/sigma_max` of a dense matrix -- 1 when its rows are
		 * orthogonal and 0 when it is singular.
		 *
		 * THE OUTER JACOBIAN CAN BE RANK DEFICIENT AND THAT IS A REAL STATE:
		 * it means the degree `g^2` is fitted in is higher than the surface
		 * family determines, so some combination of coefficients does not
		 * change the equilibrium the loop comes back with. Nothing else here
		 * can tell that apart from a hard problem, so it is measured once and
		 * used in both places that care.
		 */
		double singularValueRatio( mfem::DenseMatrix const &matrix )
		{
			Eigen::Index const n = matrix.Height();
			Eigen::MatrixXd a( n, matrix.Width() );
			for ( Eigen::Index i = 0; i < n; ++i )
				for ( Eigen::Index j = 0; j < matrix.Width(); ++j )
					a( i, j ) = matrix( static_cast<int>( i ),
					                    static_cast<int>( j ) );

			Eigen::JacobiSVD<Eigen::MatrixXd> const svd( a );
			double const largest = svd.singularValues()( 0 );
			if ( !( largest > 0.0 ) )
				return 0.0;
			return svd.singularValues()( svd.singularValues().size() - 1 )/largest;
		}

		/// The rank test, LAPACK's own: only a matrix singular to working
		/// precision fails it, so a merely ill-conditioned Newton step is left
		/// alone rather than silently damped into a different algorithm.
		double rankFloor( int n )
		{
			return static_cast<double>( n )*std::numeric_limits<double>::epsilon();
		}

		/**
		 * `R( c ) = G( c ) - c` as an mfem::Operator, with a dense differenced
		 * gradient.
		 *
		 * KINSOL solves `oper( x ) = 0`, which is this residual's own form --
		 * so unlike GradShafranov.cpp's ShiftedResidual there is nothing to
		 * shift here. That asymmetry is worth naming: NewtonSolver::Mult forms
		 * `oper( x ) - b` and KINSolver::Mult does not, and a residual written
		 * for one and handed to the other converges to a different problem.
		 */
		class FixedPointResidual : public mfem::Operator
		{
			public:
				FixedPointResidual(
					int n,
					std::function<std::vector<double>(
						std::vector<double> const & )> const &map,
					double difference, int &evaluations )
					: mfem::Operator( n ), map( map ), difference( difference ),
					  evaluations( evaluations ), gradient( n, n )
				{
				}

				void Mult( mfem::Vector const &x, mfem::Vector &y ) const override
				{
					std::vector<double> const image = evaluate( x );
					for ( int i = 0; i < x.Size(); ++i )
						y( i ) = image[ static_cast<std::size_t>( i ) ] - x( i );
				}

				/// The conditioning of the last Jacobian formed, as
				/// `sigma_min/sigma_max`. Meaningless before one has been.
				double conditioning() const { return lastConditioning; }

				mfem::Operator &GetGradient( mfem::Vector const &x ) const override
				{
					// CACHED ON THE POINT, because solveForToroidalField forms the
					// Jacobian at the start itself -- see the rank check there --
					// and KINSOL's first act is to form it at the same point. Each
					// column is two map evaluations and each map evaluation is a
					// whole equilibrium solve, so re-deriving it would make that
					// check cost 2n solves rather than none.
					if ( haveGradient && sameAs( x, gradientPoint ) )
						return gradient;

					int const n = x.Size();
					mfem::Vector forward( n ), backward( n ), plus( n ), minus( n );
					for ( int j = 0; j < n; ++j )
					{
						// Relative to the coefficient with the option as a
						// floor: c_0 is O( g^2 ) while the top coefficient can
						// legitimately be near zero, and one absolute step
						// would be far too large for the first and far too
						// small for the last.
						double const step =
							difference*std::max( std::abs( x( j ) ), 1.0 );

						plus = x;
						plus( j ) += step;
						minus = x;
						minus( j ) -= step;
						Mult( plus, forward );
						Mult( minus, backward );

						for ( int i = 0; i < n; ++i )
							gradient( i, j ) =
								( forward( i ) - backward( i ) )/( 2.0*step );
					}

					gradientPoint = x;
					haveGradient = true;
					lastConditioning = singularValueRatio( gradient );
					return gradient;
				}

			private:
				static bool sameAs( mfem::Vector const &a, mfem::Vector const &b )
				{
					if ( a.Size() != b.Size() )
						return false;
					for ( int i = 0; i < a.Size(); ++i )
						if ( a( i ) != b( i ) )
							return false;
					return true;
				}

				std::vector<double> evaluate( mfem::Vector const &x ) const
				{
					std::vector<double> in( static_cast<std::size_t>( x.Size() ) );
					for ( int i = 0; i < x.Size(); ++i )
						in[ static_cast<std::size_t>( i ) ] = x( i );

					++evaluations;
					std::vector<double> out = map( in );
					if ( out.size() != in.size() )
						throw std::invalid_argument(
							"meq::solveForToroidalField: the map returned a "
							"different number of coefficients than it was "
							"given, so the fixed point is not on the space "
							"being iterated in" );
					return out;
				}

				std::function<std::vector<double>(
					std::vector<double> const & )> const &map;
				double difference;
				int &evaluations;
				mutable mfem::DenseMatrix gradient;
				mutable mfem::Vector gradientPoint;
				mutable bool haveGradient = false;
				mutable double lastConditioning = 1.0;
		};

		/**
		 * The Newton system, solved densely and RANK-REVEALINGLY.
		 *
		 * **THIS EXISTS BECAUSE A SINGULAR OUTER JACOBIAN IS A REAL STATE AND
		 * A KRYLOV METHOD ANSWERS IT WITH AN INFINITY.** The first version was
		 * `mfem::GMRESSolver`, on the argument that an unpreconditioned Krylov
		 * method terminates in at most `n` steps on a dense `n x n` -- which is
		 * true and is not the difficulty. `Update()`'s back-substitution is
		 * `y( i ) /= h( i, i )`, unguarded, and on a singular system
		 * `h( 0, 0 )` is zero: the correction comes back infinite, the iterate
		 * goes to NaN, and **KINSOL then never stops, because every one of its
		 * convergence and failure tests is a comparison against NaN and every
		 * comparison against NaN is false.** Measured on a two-variable map
		 * with no root at all: the map was still being evaluated, at NaN, after
		 * two million calls. Bounding GMRES's iteration count does not touch
		 * it -- the divide happens on the first one.
		 *
		 * A truncated SVD returns the MINIMUM-NORM least-squares step instead,
		 * which on a zero Jacobian is exactly zero. KINSOL then takes a step it
		 * cannot improve on and terminates by its own rules, in a state this
		 * file can report. Eigen rather than a hand-rolled decomposition is
		 * this tree's standing preference, and `SurfaceFit.cpp` already depends
		 * on it for the same reason -- a soft singular tail with no gap in it.
		 *
		 * It also RECORDS what it saw. `sigma_min/sigma_max` at every Jacobian
		 * formed is the only evidence available afterwards for whether the
		 * coefficients determined their own fixed point, and KINSOL's own flag
		 * cannot carry it: a zero step is reported as `KIN_STEP_LT_STPTOL`,
		 * which MFEM maps to **converged**.
		 */
		class DenseJacobianSolver : public mfem::Solver
		{
			public:
				explicit DenseJacobianSolver( double &conditioning )
					: mfem::Solver( 0 ), conditioning( conditioning )
				{
				}

				void SetOperator( mfem::Operator const &op ) override
				{
					// KINSOL hands this whatever GetGradient() returned, which
					// here is always the dense difference. A dynamic_cast
					// rather than a static one because a wrong answer would be
					// silent, and MFEM_VERIFY is not available on a path
					// reached from C frames -- it throws.
					jacobian = dynamic_cast<mfem::DenseMatrix const *>( &op );
					height = width = op.Height();
				}

				void Mult( mfem::Vector const &b, mfem::Vector &x ) const override
				{
					x = 0.0;
					if ( jacobian == nullptr )
						return;

					Eigen::Index const n = jacobian->Height();
					Eigen::MatrixXd a( n, n );
					Eigen::VectorXd rhs( n );
					for ( Eigen::Index i = 0; i < n; ++i )
					{
						rhs( i ) = b( static_cast<int>( i ) );
						for ( Eigen::Index j = 0; j < n; ++j )
							a( i, j ) = ( *jacobian )( static_cast<int>( i ),
							                           static_cast<int>( j ) );
					}

					Eigen::JacobiSVD<Eigen::MatrixXd> const svd(
						a, Eigen::ComputeThinU | Eigen::ComputeThinV );
					Eigen::VectorXd const &sigma = svd.singularValues();
					double const largest = sigma( 0 );

					double const floor = rankFloor( static_cast<int>( n ) )*largest;

					Eigen::VectorXd step = Eigen::VectorXd::Zero( n );
					for ( Eigen::Index k = 0; k < n; ++k )
						if ( sigma( k ) > floor )
							step += ( svd.matrixU().col( k ).dot( rhs )
							          /sigma( k ) )*svd.matrixV().col( k );

					for ( Eigen::Index i = 0; i < n; ++i )
						x( static_cast<int>( i ) ) = step( i );

					double const ratio =
						largest > 0.0 ? sigma( n - 1 )/largest : 0.0;
					conditioning = std::min( conditioning, ratio );
				}

			private:
				mfem::DenseMatrix const *jacobian = nullptr;
				double &conditioning;
		};

	}

	OuterNewtonResult solveForToroidalField(
		std::vector<double> start,
		std::function<std::vector<double>( std::vector<double> const & )> const
			&map,
		OuterNewtonOptions const &options )
	{
		if ( start.empty() )
			throw std::invalid_argument(
				"meq::solveForToroidalField: no starting coefficients" );
		if ( !map )
			throw std::invalid_argument(
				"meq::solveForToroidalField: no map was given" );
		if ( options.maxIterations < 1 || !( options.difference > 0.0 )
		     || !( options.functionTolerance > 0.0 ) )
			throw std::invalid_argument(
				"meq::solveForToroidalField: the options do not describe an "
				"iteration" );

		int const n = static_cast<int>( start.size() );
		OuterNewtonResult out;
		out.coefficients = start;

		FixedPointResidual residual( n, map, options.difference,
		                             out.mapEvaluations );

		mfem::Vector solution( n );
		for ( int i = 0; i < n; ++i )
			solution( i ) = start[ static_cast<std::size_t>( i ) ];

		/*
		 * THE RANK CHECK HAPPENS HERE AND NOT INSIDE KINSOL, AND IT IS NOT AN
		 * OPTIMISATION -- A SINGULAR JACOBIAN CANNOT BE HANDED TO
		 * KIN_LINESEARCH AT ALL.
		 *
		 * With no slope under the residual there is no Newton direction and the
		 * honest step is zero. KINSOL's line search then interpolates on
		 *
		 *     rl <- -slpi rl^2 / ( 2( fnorm( u + rl p ) - fnorm( u ) - slpi rl ) )
		 *
		 * where `slpi` is the directional derivative `< F, J p >`. A zero
		 * Jacobian makes `slpi` zero whatever `p` is, and a step that does not
		 * move makes the bracket zero as well, so that is **0/0**. The iterate
		 * goes to NaN, and from there KINSOL never stops: every one of its
		 * convergence and failure tests is a comparison against NaN and every
		 * comparison against NaN is false. Measured -- a two-variable map with
		 * no root at all was still calling the map, at NaN, after **two million
		 * evaluations**, which on the real loop would be two million
		 * equilibrium solves.
		 *
		 * **BOUNDING THE LINEAR SOLVER DOES NOT REACH IT, WHICH IS WORTH
		 * RECORDING BECAUSE IT WAS THE FIRST REPAIR TRIED.** GMRES on a
		 * singular system divides by a zero pivot in its own back-substitution
		 * (`Update()`'s `y( i ) /= h( i, i )`, unguarded) on the FIRST
		 * iteration, so an iteration cap changes nothing. Replacing it with the
		 * rank-revealing solve used below does not either -- that correctly
		 * returns a step of exactly ZERO, and it is `slpi` rather than the step
		 * that has gone to zero. The line search cannot be rescued from
		 * outside it.
		 *
		 * **AND AN ALREADY-SATISFIED FIXED POINT LOOKS IDENTICAL AND IS THE
		 * OPPOSITE ANSWER.** `G = identity` has a zero Jacobian too, and every
		 * point is a fixed point of it, so returning immediately is right there
		 * where refusing would be wrong. The RESIDUAL is what separates the two
		 * cases, so it is what is asked.
		 *
		 * The Jacobian is cached on its point and KINSOL's first act is to form
		 * it at this same point, so on a healthy problem this check costs no
		 * map evaluations at all.
		 */
		residual.GetGradient( solution );
		if ( residual.conditioning() <= rankFloor( n ) )
		{
			mfem::Vector atStart( n );
			residual.Mult( solution, atStart );
			out.jacobianConditioning = residual.conditioning();

			std::ostringstream flat;
			if ( atStart.Normlinf() <= options.functionTolerance )
			{
				out.converged = true;
				flat << "converged in 0 KINSOL iterations and "
				     << out.mapEvaluations
				     << " map evaluations: the starting coefficients are "
				        "already a fixed point, to " << atStart.Normlinf();
			}
			else
			{
				out.converged = false;
				flat << "did not converge in 0 KINSOL iterations and "
				     << out.mapEvaluations
				     << " map evaluations: the outer Jacobian at the starting "
				        "coefficients is RANK DEFICIENT, its smallest singular "
				        "value being " << out.jacobianConditioning
				     << " of its largest while the residual there is "
				     << atStart.Normlinf()
				     << ". The coefficients do not determine their own fixed "
				        "point, so no step in the undetermined directions means "
				        "anything -- fit g^2 in a lower degree, or give the "
				        "family more surfaces";
			}
			out.status = flat.str();
			return out;
		}

		mfem::KINSolver newton( options.lineSearch ? KIN_LINESEARCH : KIN_NONE,
		                        true );

		double conditioning = residual.conditioning();
		DenseJacobianSolver linear( conditioning );

		// SetOperator BEFORE SetSolver: KINSolver's own SetSolver documents
		// that it must be called after, and MEQ's other KINSOL sites keep the
		// same order for the same reason.
		newton.SetOperator( residual );
		newton.SetSolver( linear );
		newton.SetPrintLevel( 0 );
		newton.SetMaxIter( options.maxIterations );
		newton.SetRelTol( options.functionTolerance );
		if ( options.stepTolerance > 0.0 )
			newton.SetAbsTol( options.stepTolerance );

		// KINSOL solves oper( x ) = 0, so the right-hand side is not used and
		// iterative_mode carries the starting point in `solution`.
		newton.iterative_mode = true;
		mfem::Vector const unused;
		newton.Mult( unused, solution );

		for ( int i = 0; i < n; ++i )
			out.coefficients[ static_cast<std::size_t>( i ) ] = solution( i );
		out.iterations = newton.GetNumIterations();
		out.jacobianConditioning = conditioning;

		/*
		 * A JACOBIAN THAT WENT SINGULAR MID-ITERATION IS NOT A CONVERGED SOLVE,
		 * AND KINSOL'S OWN FLAG WOULD SAY IT WAS. The truncated SVD keeps the
		 * step finite there rather than infinite, so the run ends instead of
		 * hanging -- but it ends on a step it could not improve, which KINSOL
		 * reports as `KIN_STEP_LT_STPTOL` and `mfem::KINSolver` turns into
		 * `converged = ( flag >= 0 )`. The recorded conditioning is the only
		 * thing that can tell that apart from a solve.
		 */
		bool const deficient = out.jacobianConditioning <= rankFloor( n );
		out.converged = newton.GetConverged() && !deficient;

		std::ostringstream status;
		status << ( out.converged ? "converged" : "did not converge" ) << " in "
		       << out.iterations << " KINSOL iterations and "
		       << out.mapEvaluations << " map evaluations, final norm "
		       << newton.GetFinalNorm();
		if ( deficient )
			status << "; the outer Jacobian went RANK DEFICIENT during the "
			          "iteration, its smallest singular value reaching "
			       << out.jacobianConditioning << " of its largest, so the "
			          "coefficients stopped determining their own fixed point";
		out.status = status.str();
		return out;
	}

}
