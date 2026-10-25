#include "meq/SafetyFactorSolve.hpp"

#include "meq/Crossings.hpp"
#include "meq/StageChecksum.hpp"

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
		 *
		 * **IT REMEMBERS EVERY POINT IT EVALUATED, AND WHAT THE MAP SAID ITS
		 * FLOOR WAS THERE**, for two reasons that are both about the answer
		 * rather than about cost. The convergence test is on the max norm of
		 * the residual AT THE POINT KINSOL RETURNS, against the floor AT THAT
		 * POINT -- and KINSOL reports neither: its final norm is a scaled 2-norm
		 * and the floor is not its notion at all. And a point asked for twice
		 * is answered once, so KINSOL's own opening evaluation of the starting
		 * coefficients -- which the rank check below has already made -- is not
		 * a second equilibrium solve with a second draw of the map's noise.
		 */
		class FixedPointResidual : public mfem::Operator
		{
			public:
				FixedPointResidual(
					int n,
					std::function<std::vector<double>(
						std::vector<double> const & )> const &map,
					std::function<double()> const &floorOfMap,
					double difference, int &evaluations )
					: mfem::Operator( n ), map( map ), floorOfMap( floorOfMap ),
					  difference( difference ), evaluations( evaluations ),
					  gradient( n, n )
				{
				}

				/// Called with the floor after every evaluation of the map.
				/// solveForToroidalField() uses it to keep KINSOL's own
				/// function-norm test at the floor of where it now is.
				void onFloor( std::function<void( double )> hook )
				{
					floorChanged = std::move( hook );
				}

				/**
				 * The residual's floor at @a x: the map's own, as it reported
				 * it there -- or at the last point evaluated, if @a x was not
				 * -- plus the floor of holding @a x at all.
				 *
				 * **THE ITERATE IS AN mfem::Vector, SO IT IS IN real_t, AND NO
				 * MAP IS EXACTER THAN THE POINT IT IS GIVEN.** Rounding `c` to
				 * the format moves it by up to half an ulp of each coefficient,
				 * which moves the residual by the Jacobian times that:
				 * `|| J ||_inf eps max( || c ||_inf, 1 )` bounds it, with a
				 * factor of two in hand. In single it is what stops an EXACT
				 * map short of a `1e-9` tolerance -- the unstable affine map
				 * below lands at a residual of order `1e-7` and can get no
				 * closer, since `4.84` itself is `4.84000015` there. In double
				 * it is `1e-15` and nothing reaches it.
				 */
				double floorAt( mfem::Vector const &x ) const
				{
					Evaluation const *known = find( x );
					double const ofMap =
						known != nullptr ? known->floor
						: history.empty() ? 0.0 : history.back().floor;
					countedHostRead( x, __LINE__ );
					double const ofIterate =
						gradientNorm
						*static_cast<double>(
							std::numeric_limits<mfem::real_t>::epsilon() )
						*std::max( static_cast<double>( x.Normlinf() ), 1.0 );
					return ofMap + ofIterate;
				}

				/// `|| G( x ) - x ||_inf` if @a x was evaluated, and NaN if it
				/// never was: what a failed loop reports, without asking the
				/// map that just failed again.
				double rememberedNormAt( mfem::Vector const &x ) const
				{
					Evaluation const *known = find( x );
					if ( known == nullptr )
						return std::numeric_limits<double>::quiet_NaN();
					double norm = 0.0;
					for ( double value : known->residual )
						norm = std::max( norm, std::abs( value ) );
					return norm;
				}

				/// `|| G( x ) - x ||_inf`, evaluating only if @a x never was.
				double normAt( mfem::Vector const &x ) const
				{
					mfem::Vector y( x.Size() );
					Mult( x, y );
					countedHostRead( y, __LINE__ );
					return y.Normlinf();
				}

				void Mult( mfem::Vector const &x, mfem::Vector &y ) const override
				{
					if ( Evaluation const *known = find( x ) )
					{
						y.HostWrite();
						for ( int i = 0; i < x.Size(); ++i )
							y( i ) = known->residual[ static_cast<std::size_t>( i ) ];
						y.Read();
						return;
					}

					/*
					 * `y` IS KINSOL'S OWN WORK VECTOR AND UNDER A DEVICE IT IS
					 * DEVICE-VALID WHEN IT ARRIVES.
					 *
					 * Filling it through `operator()` writes the HOST copy and
					 * updates no flags, so KINSOL reads back the device copy --
					 * a fresh N_Vector, i.e. zeros -- decides `|| F || = 0` and
					 * reports convergence at iteration 0 with the starting
					 * coefficients. The whole q-driven outer loop then returns
					 * its initial guess, and every printed number about the
					 * INNER solve is correct, which is what made it survive.
					 * MEASUREMENTS.md M-206.
					 *
					 * `HostWrite()` and not `HostReadWrite()`: every entry is
					 * overwritten below, so there is nothing to bring down and
					 * a read-write would pay a d2h copy for data about to be
					 * discarded.
					 *
					 * **AND THE HOST WRITE ALONE IS NOT ENOUGH, WHICH IS THE
					 * HALF ONLY THIS CALLER CAN DO.** `SundialsNVector` binds
					 * SUNDIALS' HOST array and its DEVICE array into one
					 * `Memory` -- `data.Wrap( h_ptr, d_ptr, ... )`,
					 * `linalg/sundials.cpp` -- and KINSOL's `N_Vector_Cuda`
					 * kernels read the device one. Correct flags do not move
					 * data by themselves, so the trailing `Read()` is what
					 * copies h2d into the array KINSOL will actually look at.
					 * It leaves BOTH valid, which is right here because MFEM
					 * cannot know which of the two pointers SUNDIALS will use
					 * next. With no Device configured it returns the host
					 * pointer and does nothing.
					 */
					stageChecksum( "outerResidual.x", x );
					std::vector<double> const image = evaluate( x );
					y.HostWrite();
					for ( int i = 0; i < x.Size(); ++i )
						y( i ) = image[ static_cast<std::size_t>( i ) ] - x( i );
					y.Read();
					stageChecksum( "outerResidual.y", y );

					Evaluation record;
					record.point = x;
					record.point.HostReadWrite();
					record.residual.resize( static_cast<std::size_t>( x.Size() ) );
					for ( int i = 0; i < x.Size(); ++i )
						record.residual[ static_cast<std::size_t>( i ) ] =
							image[ static_cast<std::size_t>( i ) ] - x( i );
					record.floor = floorOfMap ? std::abs( floorOfMap() ) : 0.0;
					history.push_back( std::move( record ) );
					if ( floorChanged )
						floorChanged( floorAt( x ) );
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
					double const relative = differenceAt( x );
					mfem::Vector forward( n ), backward( n ), plus( n ), minus( n );
					for ( int j = 0; j < n; ++j )
					{
						// Relative to the coefficient with the option as a
						// floor: c_0 is O( g^2 ) while the top coefficient can
						// legitimately be near zero, and one absolute step
						// would be far too large for the first and far too
						// small for the last.
						double const step =
							relative*std::max( std::abs( x( j ) ), mfem::real_t( 1.0 ) );

						/*
						 * THE PERTURBATION HAS TO SURVIVE, AND UNDER A DEVICE
						 * IT DOES NOT BY DEFAULT.
						 *
						 * `plus = x` is device-aware, so if `x`'s live copy is
						 * the device's then `plus` becomes device-valid --
						 * and `plus( j ) += step` writes the HOST copy and
						 * updates nothing. `Mult()` then funnels `plus`,
						 * which pulls the UNPERTURBED device copy back over
						 * the step just written: `forward` and `backward`
						 * come out equal, every column of the Jacobian is
						 * zero, and the step it implies is infinite. Measured:
						 * the q-driven outer loop reaches `-nan` in three map
						 * evaluations and then spins.
						 *
						 * `HostReadWrite()` and not `HostWrite()`: the copy of
						 * `x` is what is being perturbed, so the existing
						 * contents are wanted. MEASUREMENTS.md M-206.
						 */
						plus = x;
						plus.HostReadWrite();
						plus( j ) += step;
						minus = x;
						minus.HostReadWrite();
						minus( j ) -= step;
						Mult( plus, forward );
						Mult( minus, backward );

						countedHostRead( forward, __LINE__ );
						countedHostRead( backward, __LINE__ );
						for ( int i = 0; i < n; ++i )
							gradient( i, j ) =
								( forward( i ) - backward( i ) )/( 2.0*step );
					}

					gradientPoint = x;
					haveGradient = true;
					lastConditioning = singularValueRatio( gradient );
					gradientNorm = 0.0;
					for ( int i = 0; i < n; ++i )
					{
						double row = 0.0;
						for ( int j = 0; j < n; ++j )
							row += std::abs( gradient( i, j ) );
						gradientNorm = std::max( gradientNorm, row );
					}
					return gradient;
				}

				/**
				 * The column step, relative to each coefficient, for a
				 * Jacobian formed at @a x: the option, or the least step at
				 * which the map's noise costs a column at most 1%, whichever
				 * is the larger.
				 *
				 * **A COLUMN'S NOISE IS AT MOST `floor/s`**, the two ends each
				 * carrying up to the floor and the difference dividing by `2s`.
				 * `R = G - I`, so a Jacobian is of order one, and the column
				 * is good to 1% from `s = 100 floor` up; the step below is
				 * relative and multiplied by `max( |c_j|, 1 ) >= 1`, so
				 * `100 floor` as the relative step reaches that on every
				 * column, the smallest included. Above it the truncation,
				 * `s^2 F'''/6`, is the error, and it is smallest at the
				 * smallest step -- which is why this is the LEAST step that
				 * meets the noise and not the textbook balance of the two.
				 *
				 * **THE BALANCE WAS TRIED AND IS WORSE, BECAUSE THE MAP IS NOT
				 * NOISE PLUS A SMOOTH REMAINDER.** `( floor/|G| )^( 1/3 )`,
				 * Dennis and Schnabel section 5.4, assumes both. The q-driven
				 * map is PIECEWISE smooth: a surface node whose ray crosses an
				 * element face moves to the other element's arc of the DG
				 * level set, and its weight takes the other element's flux, so
				 * the map steps by the face jump -- `2.6e-05` in the
				 * coefficients on SafetyFactorSolver's fixture, from one node
				 * of the innermost surface, with the solve smooth through it.
				 * No floor sees that: the residual at each point is right to
				 * the floor, and it is the difference across the step that is
				 * wrong. A column that straddles one is off by the step's
				 * height over `2s`; at the balance's `4.4e-4` one sits
				 * `5.3e-4` below the root's `c_1`, every Jacobian near the root
				 * straddles it, and the loop contracts linearly at 0.19: 8
				 * iterations and 57 equilibria to `7e-09`, against 4 and 29 at
				 * the option's `1e-4`, which straddles none there and finishes
				 * quadratically. Every step from `1e-3` to `3e-2` is linear
				 * too, at 0.07 down to 0.01 a step. The option is not proof
				 * against a step either -- a column at `1e-4` that straddles
				 * one is off by four times as much -- and
				 * theLoopRecoversTheFieldItsTargetWasBuiltFrom bounds the
				 * iterations so that one is seen. MEASUREMENTS.md M-308b.
				 *
				 * **IN SINGLE THE OPTION'S STEP WAS THE STALL.** `1e-4` of `c`
				 * is `5e-4` on `c_0` and `1e-4` on `c_2`, against a map whose
				 * noise is `1.3e-04`, so a column carried an error of order
				 * one: measured on `examples/q-driven.toml`, the `( 2, 2 )`
				 * entry differenced at `1e-4` read `-1.06` where steps of
				 * `1e-3` and `1e-2` agree on `-0.21` to `-0.27`. Here the floor
				 * is `3.8e-04` and the step `0.038`. In double the floor is
				 * `3e-10` and the option is the larger by three orders, so the
				 * double path differences exactly as it did.
				 */
				double differenceAt( mfem::Vector const &x ) const
				{
					return std::max( difference, 100.0*floorAt( x ) );
				}

			private:
				struct Evaluation
				{
					mfem::Vector point;
					std::vector<double> residual;
					double floor = 0.0;
				};

				Evaluation const *find( mfem::Vector const &x ) const
				{
					for ( Evaluation const &known : history )
						if ( sameAs( x, known.point ) )
							return &known;
					return nullptr;
				}

				static bool sameAs( mfem::Vector const &a, mfem::Vector const &b )
				{
					/*
					 * THIS DECIDES WHETHER A JACOBIAN IS REUSED, SO A STALE
					 * READ HERE IS AN INEXACT NEWTON RATHER THAN A WRONG
					 * NUMBER. `a` is KINSOL's iterate; comparing its stale
					 * host copy against the cached point can answer "same" at
					 * a point the gradient does not belong to, and the
					 * iteration then contracts linearly with every norm
					 * looking healthy. MEASUREMENTS.md M-206.
					 */
					countedHostRead( a, __LINE__ );
					countedHostRead( b, __LINE__ );
					if ( a.Size() != b.Size() )
						return false;
					for ( int i = 0; i < a.Size(); ++i )
						if ( a( i ) != b( i ) )
							return false;
					return true;
				}

				std::vector<double> evaluate( mfem::Vector const &x ) const
				{
					// THE READ SIDE OF THE SAME HAZARD. `x` is KINSOL's
					// iterate, so under a Device the live copy can be the
					// device's and `operator()` is as raw an accessor as the
					// write above. It was benign on the case that found this
					// only because the stale host copy happened to hold the
					// same numbers -- which is not a property anything
					// guarantees.
					countedHostRead( x, __LINE__ );
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
				std::function<double()> const &floorOfMap;
				std::function<void( double )> floorChanged;
				double difference;
				int &evaluations;
				mutable std::vector<Evaluation> history;
				mutable mfem::DenseMatrix gradient;
				mutable mfem::Vector gradientPoint;
				mutable bool haveGradient = false;
				mutable double lastConditioning = 1.0;
				mutable double gradientNorm = 0.0;
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
					// silent. It is reached from KINSOL's C frames, so the
					// throw is caught by OuterKINSolver's guard there and the
					// loop is reported as failed with this message.
					jacobian = dynamic_cast<mfem::DenseMatrix const *>( &op );
					if ( jacobian == nullptr )
						throw std::logic_error(
							"meq::solveForToroidalField: the outer Jacobian "
							"is not the dense difference the Newton system's "
							"solve is written for" );
					height = width = op.Height();
				}

				void Mult( mfem::Vector const &b, mfem::Vector &x ) const override
				{
					/*
					 * THE SAME SEAM AS FixedPointResidual::Mult, ONE LEVEL
					 * DOWN, AND KINSOL OWNS BOTH VECTORS HERE TOO.
					 *
					 * `b` is KINSOL's residual and `x` the step it will apply,
					 * so under a Device both are `SundialsNVector`s carrying a
					 * host AND a device array. Eigen is host arithmetic, so
					 * `b` is funnelled down and `x` is filled on the host and
					 * pushed back up -- `HostWrite()` for the flags,
					 * `Read()` for the copy, exactly as the residual does.
					 *
					 * `x = 0.0` is NOT the way to clear it: that operator is
					 * device-aware and would leave the zeros device-valid,
					 * after which the host fill below reaches nothing KINSOL
					 * reads. MEASUREMENTS.md M-206.
					 */
					countedHostRead( b, __LINE__ );
					x.HostWrite();
					for ( int i = 0; i < x.Size(); ++i )
						x( i ) = 0.0;
					if ( jacobian == nullptr )
					{
						x.Read();
						return;
					}

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
					// Into the array KINSOL will read; see the note above.
					x.Read();

					double const ratio =
						largest > 0.0 ? sigma( n - 1 )/largest : 0.0;
					conditioning = std::min( conditioning, ratio );
				}

			private:
				mfem::DenseMatrix const *jacobian = nullptr;
				double &conditioning;
		};

#ifdef MFEM_USE_SUNDIALS
		/**
		 * mfem::KINSolver with the outer loop's two settings that MFEM does not
		 * expose: a Jacobian-vector product taken from the Jacobian, and a
		 * function-norm tolerance that can move while KINSOL runs.
		 *
		 * **THE PRODUCT IS THE ONE THAT MATTERS, AND WITHOUT IT EVERY OUTER
		 * ITERATION SPENT AN EQUILIBRIUM ON A NUMBER IT COULD NOT TRUST.**
		 * `mfem::KINSolver` presents a custom linear solver to KINSOL as
		 * `SUNLINEARSOLVER_MATRIX_ITERATIVE`, so KINSOL treats the solve as
		 * inexact and, under its default forcing term, forms `J p` after every
		 * one -- for the line search's slope `< F, J p >` and for the forcing
		 * term's `|| J p ||`. MFEM installs a product only on the path WITHOUT
		 * a custom solver, so here KINSOL fell back to its own difference
		 * quotient: one map evaluation at `u + sigma p`, `sigma` of order
		 * `sqrt( eps )`. Measured in the trace of `examples/q-driven.toml`:
		 * one equilibrium per iteration at a point `3e-4` of a step from the
		 * iterate, its residual differing from the iterate's by less than the
		 * map's own noise. In single that is a slope made of noise, handed to
		 * an Armijo test.
		 *
		 * The dense Jacobian already formed at that iterate IS `J`, so its
		 * product is exact to the difference and costs nothing.
		 */
		class OuterKINSolver : public mfem::KINSolver
		{
			public:
				OuterKINSolver( int strategy, FixedPointResidual const &residual )
					: mfem::KINSolver( strategy, true ), residual( residual )
				{
				}

				/// After SetSolver(), which is what attaches the linear solver
				/// a product can be set on.
				void takeProductsFromTheJacobian()
				{
					// KINSOL REFUSES A PRODUCT ON A LINEAR SOLVER THAT CANNOT
					// TAKE ONE, and MFEM's wrapper declares no `setatimes`:
					// its solve is the dense one below, which needs no
					// products of its own. So the wrapper is given one that
					// accepts and ignores them -- the product is for KINSOL's
					// line search, not for the solve.
					LSA->ops->setatimes = []( SUNLinearSolver, void *,
					                          SUNATimesFn ) -> int
					{
						return 0;
					};
					flag = KINSetJacTimesVecFn( sundials_mem,
					                            OuterKINSolver::jacobianTimes );
					MFEM_VERIFY( flag == KIN_SUCCESS,
					             "KINSetJacTimesVecFn() failed" );
				}

				/**
				 * EVERY CALLBACK KINSOL MAKES INTO THIS FILE, BEHIND A CATCH.
				 *
				 * KINSOL is C. The residual, the Jacobian's setup, the Newton
				 * system's solve and the Jacobian-vector product all run C++
				 * that can throw -- the map itself, a map that changes the
				 * length of its argument, an MFEM_VERIFY -- and an exception
				 * unwinding through KINSOL's frames is undefined behaviour,
				 * which works only where the C was compiled with unwind
				 * tables. So each callback catches, records the first
				 * message, and hands KINSOL `-1`, its unrecoverable failure:
				 * KINSOL returns at once with `KIN_SYSFUNC_FAIL`,
				 * `KIN_LSETUP_FAIL` or `KIN_LSOLVE_FAIL`, and
				 * solveForToroidalField() reports the loop as not converged
				 * with @ref failure() in its status.
				 *
				 * MFEM's own callbacks are the ones replaced -- the residual
				 * through `KINSetSysFunc`, the setup through `KINSetJacFn`,
				 * the solve in the linear solver's table -- and each guard
				 * calls the one it replaces, so nothing else changes.
				 * After takeProductsFromTheJacobian().
				 */
				void guardTheCallbacks()
				{
					flag = KINSetSysFunc( sundials_mem, OuterKINSolver::residualOf );
					MFEM_VERIFY( flag == KIN_SUCCESS, "KINSetSysFunc() failed" );
					flag = KINSetJacFn( sundials_mem, OuterKINSolver::setUp );
					MFEM_VERIFY( flag == KIN_SUCCESS, "KINSetJacFn() failed" );
					LSA->ops->solve = OuterKINSolver::solveTheSystem;
				}

				/// Whether any callback threw -- true even where the message
				/// could not be recorded.
				bool callbackFailed() const { return failed; }

				/// What the first callback that threw said, and where.
				std::string const &failure() const { return failureMessage; }

				/// KINSOL's own function-norm tolerance, changed mid-solve
				/// as the map's floor is learned. It is read at each
				/// iteration's stopping test, so a change takes effect there.
				void setFunctionTolerance( double tolerance ) const
				{
					if ( sundials_mem != nullptr )
						KINSetFuncNormTol( sundials_mem,
						                   static_cast<sunrealtype>( tolerance ) );
				}

			private:
				/// The object behind a KINSOL `user_data` or a SUNDIALS
				/// `content` pointer: both are the mfem::KINSolver that
				/// KINSolver::SetOperator() and SetSolver() registered, which
				/// is this object.
				static OuterKINSolver *self( void *registered )
				{
					return static_cast<OuterKINSolver *>(
						static_cast<mfem::KINSolver *>( registered ) );
				}

				/// One callback, behind the catch guardTheCallbacks()
				/// describes. Never throws.
				template <typename Body>
				int guarded( char const *where, Body const &body ) noexcept
				{
					try
					{
						return body();
					}
					catch ( std::exception const &error )
					{
						record( where, error.what() );
					}
					catch ( ... )
					{
						record( where, "an exception that is not a "
						               "std::exception" );
					}
					return -1;
				}

				void record( char const *where, char const *what ) noexcept
				{
					if ( !failureMessage.empty() )
						return;
					try
					{
						failureMessage = std::string( "evaluating " ) + where
						                 + " threw: " + what;
					}
					catch ( ... )
					{
						// An allocation failure while describing a failure;
						// KINSOL still gets its flag, and the status says
						// that much from the fallback in solveForToroidalField.
					}
					failed = true;
				}

				static int residualOf( N_Vector u, N_Vector fu,
				                       void *userData ) noexcept
				{
					OuterKINSolver *solver = self( userData );
					return solver->guarded( "the residual G( c ) - c", [ & ]
					{
						return mfem::KINSolver::Mult( u, fu, userData );
					} );
				}

				static int setUp( N_Vector u, N_Vector fu, SUNMatrix matrix,
				                  void *userData, N_Vector work1,
				                  N_Vector work2 ) noexcept
				{
					OuterKINSolver *solver = self( matrix->content );
					return solver->guarded( "the differenced outer Jacobian", [ & ]
					{
						return mfem::KINSolver::LinSysSetup( u, fu, matrix, userData,
						                                     work1, work2 );
					} );
				}

				static int solveTheSystem( SUNLinearSolver linear, SUNMatrix matrix,
				                           N_Vector step, N_Vector rhs,
				                           sunrealtype tolerance ) noexcept
				{
					OuterKINSolver *solver = self( linear->content );
					return solver->guarded( "the Newton system's solve", [ & ]
					{
						return mfem::KINSolver::LinSysSolve( linear, matrix, step, rhs,
						                                     tolerance );
					} );
				}

				static int jacobianTimes( N_Vector v, N_Vector jv, N_Vector u,
				                          sunbooleantype *newU,
				                          void *userData ) noexcept
				{
					OuterKINSolver *solver = self( userData );
					return solver->guarded( "the Jacobian-vector product", [ & ]
					{
						return solver->multiplyByTheJacobian( v, jv, u, newU );
					} );
				}

				int multiplyByTheJacobian( N_Vector v, N_Vector jv, N_Vector u,
				                           sunbooleantype *newU )
				{
					// THE SAME DEVICE SEAMS AS FixedPointResidual::Mult: all
					// three are KINSOL's vectors, read and written on the host
					// and pushed back up. MEASUREMENTS.md M-206.
					mfem::SundialsNVector const point( u );
					mfem::SundialsNVector const direction( v );
					mfem::SundialsNVector product( jv );
					if ( *newU )
					{
						// Cached on the point, and LinSysSetup has just formed
						// it at this same iterate: no map evaluation.
						jacobian = &residual.GetGradient( point );
						*newU = SUNFALSE;
					}
					auto const *dense =
						dynamic_cast<mfem::DenseMatrix const *>( jacobian );
					if ( dense == nullptr )
						throw std::logic_error(
							"meq::solveForToroidalField: the outer Jacobian is "
							"not the dense difference its product is taken from" );

					countedHostRead( direction, __LINE__ );
					product.HostWrite();
					for ( int i = 0; i < dense->Height(); ++i )
					{
						double sum = 0.0;
						for ( int j = 0; j < dense->Width(); ++j )
							sum += ( *dense )( i, j )*direction( j );
						product( i ) = static_cast<mfem::real_t>( sum );
					}
					product.Read();
					return 0;
				}

				FixedPointResidual const &residual;
				std::string failureMessage;
				bool failed = false;
		};
#endif

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
		     || !( options.functionTolerance > 0.0 )
		     || options.stepTolerance < 0.0 )
			throw std::invalid_argument(
				"meq::solveForToroidalField: the options do not describe an "
				"iteration" );

		int const n = static_cast<int>( start.size() );
		OuterNewtonResult out;
		out.coefficients = start;

		FixedPointResidual residual( n, map, options.mapFloor, options.difference,
		                             out.mapEvaluations );

		// THE CONVERGENCE THRESHOLD, at a point the map was evaluated at. The
		// tolerance is the caller's; the floor is the map's own, and no
		// residual below it is evidence of anything.
		auto threshold = [ &options ]( double floor )
		{
			return std::max( options.functionTolerance, floor );
		};

		mfem::Vector solution( n );
		for ( int i = 0; i < n; ++i )
			solution( i ) = start[ static_cast<std::size_t>( i ) ];
		/*
		 * AND THE STARTING POINT AND THE ANSWER CROSS THE SAME SEAM.
		 * `solution` becomes a KINSOL N_Vector at `newton.Mult()` below --
		 * `iterative_mode` is true, so it carries the guess IN and the answer
		 * OUT -- and both ends are filled and read here through `operator()`.
		 * Same obligation as FixedPointResidual::Mult; MEASUREMENTS.md M-206.
		 */
		solution.Read();

		/*
		 * THE RESIDUAL AT THE START FIRST, BECAUSE THE JACOBIAN'S STEP NEEDS
		 * THE FLOOR AND THE FLOOR IS THE MAP'S ANSWER TO HAVING BEEN EVALUATED.
		 * KINSOL's own first act is the same evaluation, which the residual
		 * answers from what it remembers, so this costs nothing.
		 */
		/*
		 * A MAP THAT CANNOT BE EVALUATED IS A LOOP THAT DID NOT CONVERGE, AND
		 * IT IS REPORTED, NOT THROWN -- here, ahead of KINSOL, for the same
		 * reason as inside it, where it must be (OuterKINSolver's guard). One
		 * contract wherever the map fails: `converged` false, the coefficients
		 * the last point the loop had accepted, and the map's own message in
		 * `status`.
		 */
		auto mapFailed = [ & ]( int iterations, std::string const &why )
		{
			out.converged = false;
			out.iterations = iterations;
			out.jacobianConditioning = residual.conditioning();
			out.floor = residual.floorAt( solution );
			out.finalResidual = residual.rememberedNormAt( solution );
			for ( int i = 0; i < n; ++i )
				out.coefficients[ static_cast<std::size_t>( i ) ] = solution( i );
			std::ostringstream failed;
			failed << "did not converge in " << iterations
			       << " KINSOL iterations and " << out.mapEvaluations
			       << " map evaluations: " << why
			       << "; the coefficients returned are the last the loop "
			          "accepted, at a residual of " << out.finalResidual;
			out.status = failed.str();
			return out;
		};

		double normAtStart = 0.0;
		try
		{
			normAtStart = residual.normAt( solution );
			residual.GetGradient( solution );
		}
		catch ( std::exception const &error )
		{
			return mapFailed( 0, std::string( "evaluating the map at the "
			                                  "starting coefficients threw: " )
			                     + error.what() );
		}

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
		 * map evaluations at all. It was formed above, beside the starting
		 * residual, so that a map failing at the start is reported there.
		 */
		if ( residual.conditioning() <= rankFloor( n ) )
		{
			out.jacobianConditioning = residual.conditioning();
			out.floor = residual.floorAt( solution );
			out.finalResidual = normAtStart;

			std::ostringstream flat;
			if ( normAtStart <= threshold( out.floor ) )
			{
				out.converged = true;
				flat << "converged in 0 KINSOL iterations and "
				     << out.mapEvaluations
				     << " map evaluations: the starting coefficients are "
				        "already a fixed point, to " << normAtStart
				     << " against a tolerance of " << options.functionTolerance
				     << " and a floor of " << out.floor;
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
				     << normAtStart
				     << ". The coefficients do not determine their own fixed "
				        "point, so no step in the undetermined directions means "
				        "anything -- fit g^2 in a lower degree, or give the "
				        "family more surfaces";
			}
			out.status = flat.str();
			return out;
		}

#ifdef MFEM_USE_SUNDIALS
		OuterKINSolver newton( options.lineSearch ? KIN_LINESEARCH : KIN_NONE,
		                       residual );

		double conditioning = residual.conditioning();
		DenseJacobianSolver linear( conditioning );

		// SetOperator BEFORE SetSolver: KINSolver's own SetSolver documents
		// that it must be called after, and MEQ's other KINSOL sites keep the
		// same order for the same reason.
		newton.SetOperator( residual );
		newton.SetSolver( linear );
		newton.takeProductsFromTheJacobian();
		newton.guardTheCallbacks();
		newton.SetPrintLevel( 0 );
		newton.SetMaxIter( options.maxIterations );

		/*
		 * **THE TOLERANCE IS ABSOLUTE, AND IT IS SET, BECAUSE THE ONE KINSOL WAS
		 * GIVEN BEFORE WAS NEITHER ASKED FOR NOR SEEN.** `mfem::KINSolver::Mult`
		 * turns a relative tolerance into an absolute one by multiplying it by
		 * `|| F( c0 ) ||_inf` -- and then uses its OWN absolute tolerance
		 * instead whenever that is the larger, which it was here every time:
		 * `1e-9` times a starting residual of 3.9 is `3.9e-9`, and MFEM's
		 * default absolute tolerance is `eps^( 1/3 )`, **4.9e-03 in single and
		 * 6.1e-06 in double**. So the loop stopped on KINSOL's `KIN_SUCCESS` at
		 * a norm nobody had asked for: measured on `examples/q-driven.toml` in
		 * single, a max norm of 4.20e-03, its g^2 1-2% from the closed form,
		 * reported as "converged ... final norm 0.00502" against a tolerance
		 * that read 1e-9. And it cost an extra equilibrium, the residual MFEM
		 * evaluates to scale by.
		 *
		 * So the relative path is switched off and the absolute one carries
		 * max( tolerance, floor ) -- the threshold the answer is judged by --
		 * kept at the floor of where the iterate now is by the hook below.
		 */
		newton.SetRelTol( 0.0 );
		newton.SetAbsTol( threshold( residual.floorAt( solution ) ) );
		residual.onFloor( [ &newton, &threshold ]( double floor )
		{
			newton.setFunctionTolerance( threshold( floor ) );
		} );
		if ( options.stepTolerance > 0.0 )
			newton.SetScaledStepTol( options.stepTolerance );

		/*
		 * A FRESH JACOBIAN AT EVERY ITERATE -- A NEWTON, AND NOT THE MODIFIED
		 * NEWTON KINSOL DEFAULTS TO.
		 *
		 * KINSOL re-forms the Jacobian every ten iterations unless told
		 * otherwise, which is a saving where a Jacobian is expensive beside a
		 * residual. Here it is `2n` residuals against one, so the saving is
		 * small and the price is the convergence rate: the Jacobian KINSOL kept
		 * was the one at the START, the vacuum field `g^2 = const`, a different
		 * equilibrium from the answer. Measured on `examples/q-driven.toml` in
		 * single, the loop then contracted LINEARLY -- the max norm going
		 * 0.244, 0.191, 0.149, ..., 0.017, 0.012, 0.010, 0.008, 0.0042 over
		 * eight iterations on that one Jacobian -- which is the crawl this tree
		 * treats as a defect and not as convergence. With a fresh Jacobian at
		 * every iterate the same run reaches its floor in four.
		 */
		newton.SetMaxSetupCalls( 1 );

		// KINSOL solves oper( x ) = 0, so the right-hand side is not used and
		// iterative_mode carries the starting point in `solution`.
		newton.iterative_mode = true;
		mfem::Vector const unused;
		newton.Mult( unused, solution );
		// The answer comes back wherever KINSOL last wrote it.
		countedHostRead( solution, __LINE__ );

		if ( newton.callbackFailed() )
		{
			OuterNewtonResult failed = mapFailed(
				newton.GetNumIterations(),
				newton.failure().empty()
					? std::string( "a KINSOL callback threw" )
					: newton.failure() );
			failed.jacobianConditioning = conditioning;
			return failed;
		}

		for ( int i = 0; i < n; ++i )
			out.coefficients[ static_cast<std::size_t>( i ) ] = solution( i );
		out.iterations = newton.GetNumIterations();
		out.jacobianConditioning = conditioning;

		/*
		 * **CONVERGED MEANS THE RESIDUAL AT THE ANSWER IS AT OR BELOW
		 * max( TOLERANCE, FLOOR ), AND NOTHING ELSE DOES.**
		 *
		 * `mfem::KINSolver::GetConverged()` is `flag >= 0`, which counts
		 * `KIN_STEP_LT_STPTOL` -- the step became too small to move the
		 * iterate, which is what KINSOL reports when it can make no more
		 * progress -- as converged wherever the residual stalled. And it
		 * counts as failure a line search that could not decrease a residual
		 * already inside the map's noise, which is what a noisy map does once
		 * it has been solved. Both are the same mistake: KINSOL's flag says
		 * how it STOPPED, and the question is where it got to. So how it
		 * stopped is reported and the norm is the test.
		 *
		 * The norm is the max norm, KINSOL's own stopping norm, at the point
		 * KINSOL returned -- which it evaluated, so this is read from what the
		 * residual remembers rather than solved for again. `GetFinalNorm()`
		 * is not it: that is a scaled 2-norm, `sqrt( n )` times larger at
		 * worst, and was what the status reported.
		 *
		 * A JACOBIAN THAT WENT SINGULAR MID-ITERATION IS STILL NOT A
		 * CONVERGED SOLVE. The truncated SVD keeps the step finite there
		 * rather than infinite, so the run ends instead of hanging -- but at
		 * a point the coefficients do not determine, which no residual can
		 * vouch for.
		 */
		// KINSOL returns a point it evaluated, so this is remembered and not
		// solved again; the catch is for a point it did not.
		try
		{
			out.finalResidual = residual.normAt( solution );
		}
		catch ( std::exception const &error )
		{
			return mapFailed( out.iterations,
			                  std::string( "evaluating the map at the point "
			                               "KINSOL returned threw: " )
			                  + error.what() );
		}
		out.floor = residual.floorAt( solution );
		double const judged = threshold( out.floor );
		bool const deficient = out.jacobianConditioning <= rankFloor( n );
		bool const reached = out.finalResidual <= judged;
		out.converged = reached && !deficient;

		// HOW KINSOL STOPPED, said in words: the flag itself is overwritten
		// by MFEM before it returns, but the three cases it can be are told
		// apart by what is left.
		char const *stopped =
			reached ? "with the residual within the threshold"
			: out.iterations >= options.maxIterations ? "at the iteration cap"
			: newton.GetConverged()
				? "on a step below KINSOL's step tolerance, the residual having "
				  "stalled above the threshold"
				: "on a line search that could not decrease the residual";

		std::ostringstream status;
		status << ( out.converged ? "converged" : "did not converge" ) << " in "
		       << out.iterations << " KINSOL iterations and "
		       << out.mapEvaluations << " map evaluations, final norm "
		       << out.finalResidual << " against a threshold of " << judged
		       << " -- the larger of the tolerance " << options.functionTolerance
		       << " and the round-off floor " << out.floor
		       << "; KINSOL stopped " << stopped;
		if ( deficient )
			status << "; the outer Jacobian went RANK DEFICIENT during the "
			          "iteration, its smallest singular value reaching "
			       << out.jacobianConditioning << " of its largest, so the "
			          "coefficients stopped determining their own fixed point";
		out.status = status.str();
		return out;
#else
		out.status = "meq::safetyFactorSolve: the q-driven outer Newton needs KINSOL, and MFEM was built without MFEM_USE_SUNDIALS";
		out.converged = false;
		return out;
#endif
	}

}
