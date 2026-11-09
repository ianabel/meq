#define BOOST_TEST_MODULE RecoveryPrecision
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "mfem.hpp"

#include "ExtractionCrossing.hpp"

/*
 * WHY MFEM'S ELEMENT-LOCAL HDG RECOVERY LOSES PRECISION IN SINGLE -- the
 * instrument behind MEASUREMENTS.md's M-RD-1 onwards. BUILT, NOT REGISTERED:
 * it asserts nothing.
 *
 * M-371 measured that in single the trace psihat is solved to ~0.4 eps while
 * RecoverFEMSolution()'s psi_h and q_h carry ~200x and ~70x of that, and that
 * Reconstruct()'s psi* has a floor of its own near 5.5e-7. Three readings of
 * the recovery are separated here:
 *
 *   H1  the arithmetic of the local solves -- the stored local matrices are
 *       fine, and a double recomputation from the SAME float data is quiet,
 *       as is a refinement pass with a double residual;
 *   H2  cancellation in B A^-1 b_u - b_p, which a real_t refinement cannot
 *       recover and a double residual or a rescaling can;
 *   H3  the stored float data itself, so that even an exact solve of the
 *       stored float system is noisy.
 *
 * HOW THE LOCAL DATA IS REACHED. DarcyHybridization keeps A and the Schur
 * complement FACTORED in place, B raw, and the face blocks C^T and E raw, all
 * protected or private, and MEQ keeps its DarcyForm private. They are read
 * here through the explicit-instantiation access idiom (an explicit
 * instantiation is exempt from access checking), read only. The UNFACTORED A
 * and D are re-assembled with the very integrators the hybridization holds,
 * in the order DarcyForm::Assemble() adds them, and the study checks that
 * factoring them as ComputeElementH() does reproduces the stored factors BIT
 * FOR BIT, and that MultInv() on the local right-hand side formed here
 * reproduces RecoverFEMSolution()'s fields bit for bit, before any number is
 * read. The same is done for Reconstruct()'s local problem.
 *
 * THE NOISE MEASURE is M-371's: 21 states 1e-5 apart in the source mix A of
 * ExtractionCrossing.hpp's linear box, each its own solve; sigma is the root
 * mean square of the second differences over sqrt( 6 ), of the field's
 * largest value at the middle state, printed as root mean square over dofs /
 * worst dof. THE LOCAL MATRICES DO NOT DEPEND ON A -- only the source load
 * does -- so they are the same bits in every state, and a rounding of their
 * entries is a FIXED perturbation that a second difference cancels. H3 can
 * therefore only show in sigma through the per-state inputs (the load, the
 * trace), which (b) splits.
 *
 * Every method is applied to the same per-state inputs:
 *
 *   (a)   RecoverFEMSolution(), as MEQ ran it;
 *   (b)   the local system [ A sB^T ; B D ] solved by LU in DOUBLE from the
 *         float-stored A, B, D, C^T, E, loads and trace, the local right-hand
 *         side formed in double;
 *   (b')  the same, from the local right-hand side as MFEM formed it in real_t;
 *   (c)   (a) plus one refinement pass, residual in real_t, correction by
 *         MultInv() on the stored factors;
 *   (d)   (a) plus one pass with the residual in double (right-hand side
 *         formed in double), correction on the stored factors, sum in double;
 *   (d')  (d) against MFEM's real_t right-hand side;
 *   (d2)  two passes of (d);
 *   (e1)  MultInv()'s Schur recovery in real_t on blocks scaled per block by
 *         powers of two (the h-power rescaling);
 *   (e2)  the same, scaled by a symmetric Ruiz equilibration of the full
 *         local matrix, rounded to powers of two (exact in any precision);
 *   (f)   one LU of the full local matrix in real_t, no Schur complement;
 *   (f2)  (f) on the Ruiz-equilibrated matrix.
 *
 * (b) is also split into the noise the trace carries and the noise the load
 * carries, by freezing the other at the middle state's value.
 *
 * RECONSTRUCT() is two steps: ReconstructTotalFlux() makes u_t in RT_k from
 * ( u, p, psihat ), and ReconstructFluxAndPot() solves one local problem per
 * element -- u*, p* and an enriched trace, closed by p*'s element mean -- by
 * DenseMatrixInverse in real_t. The study re-assembles that local problem
 * (bit for bit), then solves it in double, refines it, equilibrates it,
 * feeds it fields recovered by (d), and forms its load in double: the load
 * is linear in u_t and p_c, so it is G u_t + P p_c with G and P fixed per
 * element, whose columns come from the same real_t code on unit inputs.
 */

// ------------------------------------------------------------------------
// The access idiom. Read only; nothing below writes through these.
// ------------------------------------------------------------------------
#if defined( __GNUC__ ) && !defined( __clang__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif
namespace recoveryaccess
{
	template <typename Tag, typename Tag::type Member>
	struct Rob
	{
		friend typename Tag::type get( Tag )
		{
			return Member;
		}
	};

#define MEQ_RECOVERY_ROB( NAME, TYPE, MEMBER ) \
	struct NAME \
	{ \
		using type = TYPE; \
		friend type get( NAME ); \
	}; \
	template struct Rob<NAME, MEMBER>;

	using Hyb = mfem::DarcyHybridization;
	using Base = mfem::Hybridization;
	MEQ_RECOVERY_ROB( AfData, mfem::Array<mfem::real_t> Base::*, &Base::Af_data )
	MEQ_RECOVERY_ROB( AfIpiv, mfem::Array<int> Base::*, &Base::Af_ipiv )
	MEQ_RECOVERY_ROB( AfOffsets, mfem::Array<int> Base::*, &Base::Af_offsets )
	MEQ_RECOVERY_ROB( AfFOffsets, mfem::Array<int> Base::*, &Base::Af_f_offsets )
	MEQ_RECOVERY_ROB( BfData, mfem::Vector Hyb::*, &Hyb::Bf_data )
	MEQ_RECOVERY_ROB( BfOffsets, mfem::Array<int> Hyb::*, &Hyb::Bf_offsets )
	MEQ_RECOVERY_ROB( DfData, mfem::Vector Hyb::*, &Hyb::Df_data )
	MEQ_RECOVERY_ROB( DfIpiv, mfem::Array<int> Hyb::*, &Hyb::Df_ipiv )
	MEQ_RECOVERY_ROB( DfOffsets, mfem::Array<int> Hyb::*, &Hyb::Df_offsets )
	MEQ_RECOVERY_ROB( DfFOffsets, mfem::Array<int> Hyb::*, &Hyb::Df_f_offsets )
	MEQ_RECOVERY_ROB( Bsym, bool Hyb::*, &Hyb::bsym )
	MEQ_RECOVERY_ROB( GetFDofs, void ( Hyb::* )( int, mfem::Array<int> & ) const, &Hyb::GetFDofs )
	MEQ_RECOVERY_ROB( GetElementFaces, void ( Hyb::* )( int, mfem::Array<int> & ) const,
	                  &Hyb::GetElementFaces )
	MEQ_RECOVERY_ROB( GetCtFaceMatrix, void ( Hyb::* )( int, int, mfem::DenseMatrix & ) const,
	                  &Hyb::GetCtFaceMatrix )
	MEQ_RECOVERY_ROB( GetEFaceMatrix, void ( Hyb::* )( int, int, mfem::DenseMatrix & ) const,
	                  &Hyb::GetEFaceMatrix )
	MEQ_RECOVERY_ROB( MultInv,
	                  void ( Hyb::* )( int, mfem::Vector const &, mfem::Vector const &, mfem::Vector &,
	                                   mfem::Vector &, bool, mfem::Vector * ) const,
	                  &Hyb::MultInv )

	using Form = mfem::DarcyForm;
	MEQ_RECOVERY_ROB( Reconstruction, std::unique_ptr<Form> Form::*, &Form::reconstruction )
	MEQ_RECOVERY_ROB( FormMu, std::unique_ptr<mfem::BilinearForm> Form::*, &Form::M_u )
	MEQ_RECOVERY_ROB( FormMp, std::unique_ptr<mfem::BilinearForm> Form::*, &Form::M_p )
	MEQ_RECOVERY_ROB( FormB, std::unique_ptr<mfem::MixedBilinearForm> Form::*, &Form::B )

	using Solver = meq::GradShafranovSolver;
	MEQ_RECOVERY_ROB( SolverDarcy, std::unique_ptr<mfem::DarcyForm> Solver::*, &Solver::darcy )
	MEQ_RECOVERY_ROB( SolverDarcyRhs, mfem::BlockVector Solver::*, &Solver::darcyRhs )
	MEQ_RECOVERY_ROB( SolverDarcySolution, mfem::BlockVector Solver::*, &Solver::darcySolution )
	MEQ_RECOVERY_ROB( SolverTraceX, mfem::Vector Solver::*, &Solver::traceX )
	MEQ_RECOVERY_ROB( SolverTotalFlux, mfem::GridFunction Solver::*, &Solver::totalFluxGf )
	MEQ_RECOVERY_ROB( SolverEnrichedFlux, mfem::GridFunction Solver::*, &Solver::enrichedFluxGf )
	MEQ_RECOVERY_ROB( SolverPostProcessed, mfem::GridFunction Solver::*, &Solver::postProcessedGf )
	MEQ_RECOVERY_ROB( SolverEnrichedTrace, mfem::GridFunction Solver::*, &Solver::enrichedTraceGf )
#undef MEQ_RECOVERY_ROB
}
#if defined( __GNUC__ ) && !defined( __clang__ )
#pragma GCC diagnostic pop
#endif

// LAPACK in double, whatever real_t is. MKL is on the link line through MFEM.
extern "C"
{
	void dgetrf_( int const *, int const *, double *, int const *, int *, int * );
	void dgetrs_( char const *, int const *, int const *, double const *, int const *, int const *, double *,
	              int const *, int * );
	void dgesvd_( char const *, char const *, int const *, int const *, double *, int const *, double *,
	              double *, int const *, double *, int const *, double *, int const *, int * );
}

namespace
{
	// get() is found by argument-dependent lookup on the tags.

	double const centre = 0.30;
	double const spacing = 1.0e-5;
	int const states = 21;
	double const eps = std::numeric_limits<mfem::real_t>::epsilon();

	// ---------------------------------------------------------------- dense
	/// Column-major, double.
	struct Dense
	{
		int rows = 0, cols = 0;
		std::vector<double> a;
		Dense() = default;
		Dense( int m, int n ) : rows( m ), cols( n ), a( static_cast<std::size_t>( m*n ), 0.0 ) {}
		double &operator()( int i, int j ) { return a[ static_cast<std::size_t>( i + j*rows ) ]; }
		double operator()( int i, int j ) const { return a[ static_cast<std::size_t>( i + j*rows ) ]; }
	};

	Dense promote( mfem::DenseMatrix const &m )
	{
		Dense out( m.Height(), m.Width() );
		for ( int j = 0; j < m.Width(); ++j )
			for ( int i = 0; i < m.Height(); ++i )
				out( i, j ) = static_cast<double>( m( i, j ) );
		return out;
	}

	std::vector<double> promote( mfem::Vector const &v )
	{
		std::vector<double> out( static_cast<std::size_t>( v.Size() ) );
		for ( int i = 0; i < v.Size(); ++i )
			out[ static_cast<std::size_t>( i ) ] = static_cast<double>( v( i ) );
		return out;
	}

	/// x = M^-1 b by LAPACK LU in double.
	std::vector<double> luSolve( Dense m, std::vector<double> b )
	{
		int const n = m.rows;
		std::vector<int> ipiv( static_cast<std::size_t>( n ) );
		int info = 0, one = 1;
		dgetrf_( &n, &n, m.a.data(), &n, ipiv.data(), &info );
		if ( info != 0 )
			throw std::runtime_error( "dgetrf failed" );
		char const trans = 'N';
		dgetrs_( &trans, &n, &one, m.a.data(), &n, ipiv.data(), b.data(), &n, &info );
		return b;
	}

	/// The 2-norm condition number, by SVD in double.
	double cond2( Dense m )
	{
		int const r = m.rows, c = m.cols, k = std::min( r, c );
		std::vector<double> s( static_cast<std::size_t>( k ) );
		double dummy = 0.0;
		int one = 1, info = 0, lwork = -1;
		char const no = 'N';
		double query = 0.0;
		dgesvd_( &no, &no, &r, &c, m.a.data(), &r, s.data(), &dummy, &one, &dummy, &one, &query, &lwork, &info );
		lwork = static_cast<int>( query );
		std::vector<double> work( static_cast<std::size_t>( lwork ) );
		dgesvd_( &no, &no, &r, &c, m.a.data(), &r, s.data(), &dummy, &one, &dummy, &one, work.data(), &lwork,
		         &info );
		return s.back() > 0.0 ? s.front()/s.back() : std::numeric_limits<double>::infinity();
	}

	double normInf( std::vector<double> const &v )
	{
		double m = 0.0;
		for ( double x : v )
			m = std::max( m, std::abs( x ) );
		return m;
	}

	double normInf( Dense const &m )
	{
		double best = 0.0;
		for ( int i = 0; i < m.rows; ++i )
		{
			double row = 0.0;
			for ( int j = 0; j < m.cols; ++j )
				row += std::abs( m( i, j ) );
			best = std::max( best, row );
		}
		return best;
	}

	std::vector<double> multiply( Dense const &m, std::vector<double> const &x )
	{
		std::vector<double> y( static_cast<std::size_t>( m.rows ), 0.0 );
		for ( int j = 0; j < m.cols; ++j )
			for ( int i = 0; i < m.rows; ++i )
				y[ static_cast<std::size_t>( i ) ] += m( i, j )*x[ static_cast<std::size_t>( j ) ];
		return y;
	}

	double powerOfTwo( double x )
	{
		return std::exp2( std::round( std::log2( x ) ) );
	}

	/// Symmetric Ruiz equilibration, rounded to powers of two so that the
	/// scaling itself is exact in any precision.
	std::vector<double> ruiz( Dense const &m, int iterations = 12 )
	{
		int const n = m.rows;
		std::vector<double> d( static_cast<std::size_t>( n ), 1.0 );
		for ( int it = 0; it < iterations; ++it )
		{
			std::vector<double> r( static_cast<std::size_t>( n ), 0.0 );
			for ( int j = 0; j < n; ++j )
				for ( int i = 0; i < n; ++i )
					r[ static_cast<std::size_t>( i ) ] = std::max( r[ static_cast<std::size_t>( i ) ],
						std::abs( d[ static_cast<std::size_t>( i ) ]*m( i, j )*d[ static_cast<std::size_t>( j ) ] ) );
			for ( int i = 0; i < n; ++i )
				if ( r[ static_cast<std::size_t>( i ) ] > 0.0 )
					d[ static_cast<std::size_t>( i ) ] /= std::sqrt( r[ static_cast<std::size_t>( i ) ] );
		}
		for ( double &x : d )
			x = powerOfTwo( x );
		return d;
	}

	Dense scaled( Dense const &m, std::vector<double> const &d )
	{
		Dense out = m;
		for ( int j = 0; j < m.cols; ++j )
			for ( int i = 0; i < m.rows; ++i )
				out( i, j ) = d[ static_cast<std::size_t>( i ) ]*m( i, j )*d[ static_cast<std::size_t>( j ) ];
		return out;
	}

	mfem::DenseMatrix demote( Dense const &m )
	{
		mfem::DenseMatrix out( m.rows, m.cols );
		for ( int j = 0; j < m.cols; ++j )
			for ( int i = 0; i < m.rows; ++i )
				out( i, j ) = static_cast<mfem::real_t>( m( i, j ) );
		return out;
	}

	mfem::Vector demote( std::vector<double> const &v )
	{
		mfem::Vector out( static_cast<int>( v.size() ) );
		for ( std::size_t i = 0; i < v.size(); ++i )
			out( static_cast<int>( i ) ) = static_cast<mfem::real_t>( v[ i ] );
		return out;
	}

	double median( std::vector<double> v )
	{
		if ( v.empty() )
			return 0.0;
		std::sort( v.begin(), v.end() );
		return v[ v.size()/2 ];
	}

	double maxOf( std::vector<double> const &v )
	{
		return v.empty() ? 0.0 : *std::max_element( v.begin(), v.end() );
	}

	// --------------------------------------------------------------- sigma
	struct Sigma
	{
		double rms = 0.0, worst = 0.0;
	};

	Sigma sigma( std::vector<std::vector<double>> const &series )
	{
		std::size_t const n = series.front().size();
		double const scale = std::max( normInf( series[ series.size()/2 ] ), 1.0e-300 );
		Sigma out;
		double sumSquares = 0.0;
		for ( std::size_t i = 0; i < n; ++i )
		{
			double sum = 0.0;
			for ( std::size_t j = 1; j + 1 < series.size(); ++j )
			{
				double const d2 = series[ j + 1 ][ i ] - 2.0*series[ j ][ i ] + series[ j - 1 ][ i ];
				sum += d2*d2;
			}
			double const s = std::sqrt( sum/static_cast<double>( series.size() - 2 ) )/std::sqrt( 6.0 )/scale;
			out.worst = std::max( out.worst, s );
			sumSquares += s*s;
		}
		out.rms = std::sqrt( sumSquares/static_cast<double>( n ) );
		return out;
	}

	double relativeDistance( std::vector<double> const &x, std::vector<double> const &reference )
	{
		double d = 0.0;
		for ( std::size_t i = 0; i < x.size(); ++i )
			d = std::max( d, std::abs( x[ i ] - reference[ i ] ) );
		return d/std::max( normInf( reference ), 1.0e-300 );
	}

	std::vector<double> roundedToReal( std::vector<double> v )
	{
		for ( double &x : v )
			x = static_cast<double>( static_cast<mfem::real_t>( x ) );
		return v;
	}

	// ------------------------------------------------- one element's data
	struct Element
	{
		int na = 0, nd = 0;
		mfem::Array<int> uDofs, pDofs;
		mfem::DenseMatrix A, B, D;          ///< raw, real_t, re-assembled and checked
		mfem::Vector buReal, bpReal;        ///< the local right-hand side as MFEM forms it
		std::vector<double> buDouble, bpDouble; ///< the same, formed in double from real_t pieces
		/// The pieces, for the diagnostics: the loads alone and the trace terms alone.
		std::vector<double> buLoad, bpLoad, buTrace, bpTrace;
	};

	/// The full local matrix [ A sB^T ; B D ] in double, s = -1 under bsym.
	Dense fullMatrix( Element const &e, double s )
	{
		int const n = e.na + e.nd;
		Dense m( n, n );
		for ( int j = 0; j < e.na; ++j )
			for ( int i = 0; i < e.na; ++i )
				m( i, j ) = e.A( i, j );
		for ( int j = 0; j < e.na; ++j )
			for ( int i = 0; i < e.nd; ++i )
			{
				m( e.na + i, j ) = e.B( i, j );
				m( j, e.na + i ) = s*static_cast<double>( e.B( i, j ) );
			}
		for ( int j = 0; j < e.nd; ++j )
			for ( int i = 0; i < e.nd; ++i )
				m( e.na + i, e.na + j ) = e.D( i, j );
		return m;
	}

	std::vector<double> concat( std::vector<double> a, std::vector<double> const &b )
	{
		a.insert( a.end(), b.begin(), b.end() );
		return a;
	}

	/// MultInv()'s Schur recovery transcribed, factoring its own copies of the
	/// blocks in real_t: with the raw blocks it is ComputeElementH() followed
	/// by MultInv(), so it is the vehicle for the scaled variants.
	void schurRecoveryReal( mfem::DenseMatrix const &aIn, mfem::DenseMatrix const &b, mfem::DenseMatrix const &dIn,
	                        bool bsym, mfem::Vector const &bu, mfem::Vector const &bp, mfem::Vector &u,
	                        mfem::Vector &p )
	{
		int const na = aIn.Height(), nd = dIn.Height();
		mfem::DenseMatrix a( aIn ), s( dIn );
		mfem::Array<int> ipivA( na ), ipivS( nd );
		mfem::LUFactors luA( a.Data(), ipivA.GetData() );
		luA.Factor( na );
		mfem::DenseMatrix aiBt( na, nd );
		aiBt.Transpose( b );
		if ( !bsym )
			aiBt.Neg();
		luA.Solve( na, nd, aiBt.Data() );
		mfem::AddMult( b, aiBt, s );
		mfem::LUFactors luS( s.Data(), ipivS.GetData() );
		luS.Factor( nd );

		u = bu;
		luA.Solve( na, 1, u.GetData() );
		p.SetSize( nd );
		b.Mult( u, p );
		p -= bp;
		luS.Solve( nd, 1, p.GetData() );
		p.Neg();
		mfem::Vector w( na );
		b.MultTranspose( p, w );
		luA.Solve( na, 1, w.GetData() );
		if ( bsym )
			u += w;
		else
			u -= w;
	}

	/// The same blocks scaled by diag( du, dp ) on both sides, solved by
	/// schurRecoveryReal(), and scaled back.
	void scaledSchurReal( Element const &e, bool bsym, std::vector<double> const &d, mfem::Vector &u,
	                      mfem::Vector &p )
	{
		auto du = [ & ]( int i ) { return static_cast<mfem::real_t>( d[ static_cast<std::size_t>( i ) ] ); };
		auto dp = [ & ]( int i ) { return static_cast<mfem::real_t>( d[ static_cast<std::size_t>( e.na + i ) ] ); };
		mfem::DenseMatrix a( e.A ), b( e.B ), dd( e.D );
		for ( int j = 0; j < e.na; ++j )
			for ( int i = 0; i < e.na; ++i )
				a( i, j ) *= du( i )*du( j );
		for ( int j = 0; j < e.na; ++j )
			for ( int i = 0; i < e.nd; ++i )
				b( i, j ) *= dp( i )*du( j );
		for ( int j = 0; j < e.nd; ++j )
			for ( int i = 0; i < e.nd; ++i )
				dd( i, j ) *= dp( i )*dp( j );
		mfem::Vector bu( e.buReal ), bp( e.bpReal );
		for ( int i = 0; i < e.na; ++i )
			bu( i ) *= du( i );
		for ( int i = 0; i < e.nd; ++i )
			bp( i ) *= dp( i );
		schurRecoveryReal( a, b, dd, bsym, bu, bp, u, p );
		for ( int i = 0; i < e.na; ++i )
			u( i ) *= du( i );
		for ( int i = 0; i < e.nd; ++i )
			p( i ) *= dp( i );
	}

	/// One LU of the full local matrix in real_t, optionally equilibrated.
	void fullLuReal( Element const &e, double s, std::vector<double> const *d, mfem::Vector &u, mfem::Vector &p )
	{
		int const n = e.na + e.nd;
		Dense m = fullMatrix( e, s );
		std::vector<double> scale( static_cast<std::size_t>( n ), 1.0 );
		if ( d )
			scale = *d;
		mfem::DenseMatrix mr = demote( scaled( m, scale ) );
		mfem::Vector x( n );
		for ( int i = 0; i < e.na; ++i )
			x( i ) = e.buReal( i )*static_cast<mfem::real_t>( scale[ static_cast<std::size_t>( i ) ] );
		for ( int i = 0; i < e.nd; ++i )
			x( e.na + i ) = e.bpReal( i )*static_cast<mfem::real_t>( scale[ static_cast<std::size_t>( e.na + i ) ] );
		mfem::Array<int> ipiv( n );
		mfem::LUFactors lu( mr.Data(), ipiv.GetData() );
		lu.Factor( n );
		lu.Solve( n, 1, x.GetData() );
		u.SetSize( e.na );
		p.SetSize( e.nd );
		for ( int i = 0; i < e.na; ++i )
			u( i ) = x( i )*static_cast<mfem::real_t>( scale[ static_cast<std::size_t>( i ) ] );
		for ( int i = 0; i < e.nd; ++i )
			p( i ) = x( e.na + i )*static_cast<mfem::real_t>( scale[ static_cast<std::size_t>( e.na + i ) ] );
	}

	// --------------------------------------------------- one state, read
	/// Everything this study reads from one solved state.
	struct State
	{
		std::map<std::string, std::vector<double>> u, p;   ///< global, per method
		std::vector<double> trace, buGlobal, bpGlobal;
		std::vector<double> psiStar, totalFlux, enrichedFlux;
		std::map<std::string, std::vector<double>> psiStarMethod;
		std::map<std::string, std::vector<double>> totalFluxMethod;
		std::map<std::string, std::vector<double>> reconstructionLoad;
		double totalFluxDifference = 0.0, linearisedLoadDifference = 0.0;
		// Checks, worst over elements.
		double factorADifference = 0.0, factorSDifference = 0.0;
		double recoveryUDifference = 0.0, recoveryPDifference = 0.0;
		double reconstructionDifference = 0.0;
		// Per-element diagnostics, middle state only.
		std::map<std::string, std::vector<double>> diagnostics;
		std::map<std::string, std::vector<double>> reconstructionDiagnostics;
	};

	void scatter( std::vector<double> &global, mfem::Array<int> const &dofs, mfem::Vector const &local )
	{
		for ( int i = 0; i < dofs.Size(); ++i )
			global[ static_cast<std::size_t>( dofs[ i ] ) ] = static_cast<double>( local( i ) );
	}

	void scatter( std::vector<double> &global, mfem::Array<int> const &dofs, std::vector<double> const &local,
	              std::size_t offset = 0 )
	{
		for ( int i = 0; i < dofs.Size(); ++i )
			global[ static_cast<std::size_t>( dofs[ i ] ) ] = local[ offset + static_cast<std::size_t>( i ) ];
	}

	/// Re-assembles the unfactored D of every element as
	/// DarcyForm::AssemblePotHDGFaces() adds it: interior faces in face order,
	/// then the boundary faces in boundary-element order.
	std::vector<mfem::DenseMatrix> assemblePotentialBlocks( mfem::DarcyHybridization &hyb,
	                                                       mfem::FiniteElementSpace &fesP )
	{
		mfem::Mesh &mesh = *fesP.GetMesh();
		std::vector<mfem::DenseMatrix> d( static_cast<std::size_t>( mesh.GetNE() ) );
		for ( int el = 0; el < mesh.GetNE(); ++el )
		{
			int const s = fesP.GetFE( el )->GetDof()*fesP.GetVDim();
			d[ static_cast<std::size_t>( el ) ].SetSize( s );
			d[ static_cast<std::size_t>( el ) ] = 0.0;
		}
		mfem::BilinearFormIntegrator *interior = hyb.GetPotConstraintIntegrator();
		mfem::DenseMatrix elmat, aux;
		if ( interior )
			for ( int f = 0; f < mesh.GetNumFaces(); ++f )
			{
				if ( !mesh.FaceIsInterior( f ) )
					continue;
				mfem::FaceElementTransformations *ftr = mesh.GetFaceElementTransformations( f );
				int const el1 = ftr->Elem1No, el2 = ftr->Elem2No;
				mfem::FiniteElement const *fe1 = fesP.GetFE( el1 );
				mfem::FiniteElement const *fe2 = fesP.GetFE( el2 );
				interior->AssembleHDGFaceMatrix( *hyb.TraceFE( f ), *fe1, *fe2, *ftr, elmat );
				int const n1 = fe1->GetDof()*fesP.GetVDim(), n2 = fe2->GetDof()*fesP.GetVDim();
				mfem::DenseMatrix &d1 = d[ static_cast<std::size_t>( el1 ) ];
				mfem::DenseMatrix &d2 = d[ static_cast<std::size_t>( el2 ) ];
				for ( int j = 0; j < n1; ++j )
					for ( int i = 0; i < n1; ++i )
						d1( i, j ) += elmat( i, j );
				for ( int j = 0; j < n2; ++j )
					for ( int i = 0; i < n2; ++i )
						d2( i, j ) += elmat( n1 + i, n1 + j );
			}
		int const integrators = hyb.NumBdrPotConstraintIntegrators();
		if ( integrators > 0 )
		{
			mfem::Array<int> marker( mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0 );
			marker = 0;
			for ( int k = 0; k < integrators; ++k )
			{
				mfem::Array<int> *m = hyb.GetBdrPotConstraintIntegratorMarker( k );
				if ( m == nullptr )
				{
					marker = 1;
					break;
				}
				for ( int i = 0; i < marker.Size(); ++i )
					marker[ i ] |= ( *m )[ i ];
			}
			for ( int be = 0; be < fesP.GetNBE(); ++be )
			{
				int const attribute = mesh.GetBdrAttribute( be );
				if ( marker[ attribute - 1 ] == 0 )
					continue;
				if ( !mesh.GetBdrFaceTransformations( be ) )
					continue;
				int const face = mesh.GetBdrElementFaceIndex( be );
				mfem::FaceElementTransformations *ftr = mesh.GetFaceElementTransformations( face );
				mfem::FiniteElement const *fe = fesP.GetFE( ftr->Elem1No );
				elmat.SetSize( 0 );
				for ( int k = 0; k < integrators; ++k )
				{
					mfem::Array<int> *m = hyb.GetBdrPotConstraintIntegratorMarker( k );
					if ( m && ( *m )[ attribute - 1 ] == 0 )
						continue;
					hyb.GetBdrPotConstraintIntegrator( k ).AssembleHDGFaceMatrix( *hyb.TraceFE( face ), *fe, *fe,
					                                                              *ftr, aux );
					if ( elmat.Height() > 0 )
						elmat += aux;
					else
						elmat = aux;
				}
				if ( elmat.Height() == 0 )
					continue;
				int const n = fe->GetDof()*fesP.GetVDim();
				mfem::DenseMatrix &d1 = d[ static_cast<std::size_t>( ftr->Elem1No ) ];
				for ( int j = 0; j < n; ++j )
					for ( int i = 0; i < n; ++i )
						d1( i, j ) += elmat( i, j );
			}
		}
		return d;
	}

	// ------------------------------------------- Reconstruct()'s local problem
	/// DarcyForm::ReconstructFluxAndPot()'s element problem, assembled as that
	/// routine's per-element lambda assembles it on MEQ's linear path (no
	/// lifted potential, no frozen flux law, a linear face constraint).
	struct ReconstructionProblem
	{
		mfem::DenseMatrix elmat;
		mfem::Vector rhs;
		int ndofU = 0, ndofP = 0;
		mfem::Array<int> dofsP, vdofsU;
	};

	class Reconstructor
	{
		public:
			Reconstructor( mfem::DarcyForm &darcy, mfem::DarcyHybridization &hyb, mfem::GridFunction const &totalFlux,
			               mfem::GridFunction const &enrichedTrace, mfem::BlockVector const &solution )
				: rec( *( darcy.*get( recoveryaccess::Reconstruction{} ) ) ),
				  hyb( hyb ),
				  ut( totalFlux ),
				  fesTr( *const_cast<mfem::FiniteElementSpace *>( enrichedTrace.FESpace() ) ),
				  fesPc( *const_cast<mfem::FiniteElementSpace *>( darcy.PotentialFESpace() ) ),
				  pc( const_cast<mfem::FiniteElementSpace *>( darcy.PotentialFESpace() ),
				      const_cast<mfem::Vector &>( solution.GetBlock( 1 ) ), 0 ),
				  divergence( &ut, rec.PotentialFESpace()->GetVDim() ),
				  load( divergence )
			{
			}

			void assemble( int z, ReconstructionProblem &out )
			{
				mfem::FiniteElementSpace &fesU = *rec.FluxFESpace();
				mfem::FiniteElementSpace &fesP = *rec.PotentialFESpace();
				mfem::FiniteElementSpace const &fesUt = *ut.FESpace();
				mfem::Mesh &mesh = *fesU.GetMesh();
				int const neq = fesP.GetVDim();
				mfem::BilinearForm *mU = ( rec.*get( recoveryaccess::FormMu{} ) ).get();
				mfem::BilinearForm *mP = ( rec.*get( recoveryaccess::FormMp{} ) ).get();
				mfem::MixedBilinearForm *b = ( rec.*get( recoveryaccess::FormB{} ) ).get();
				mfem::BilinearFormIntegrator *cBfi = hyb.GetFluxConstraintIntegrator();
				mfem::BilinearFormIntegrator *cBfiP = hyb.GetPotConstraintIntegrator();

				mesh.GetElementEdges( z, faces, oris );
				fesU.GetElementVDofs( z, out.vdofsU, doftrans );
				fesP.GetElementVDofs( z, out.dofsP, doftrans );
				int const ndofU = out.vdofsU.Size();
				int const ndofP = out.dofsP.Size();
				int const ndofPs = ndofP/neq;
				out.ndofU = ndofU;
				out.ndofP = ndofP;
				int ndofTr = 0;
				for ( int f : faces )
					ndofTr += neq*fesTr.GetFaceElement( f )->GetDof();
				int const n = ndofU + ndofP + ndofTr;
				out.elmat.SetSize( n, n );
				out.elmat = 0.0;
				out.rhs.SetSize( n );
				out.rhs = 0.0;

				mfem::FiniteElement const *feP = fesP.GetFE( z );
				mesh.GetElementTransformation( z, &eltrans );
				mfem::ElementTransformation *tr = &eltrans;
				fesPc.GetElementVDofs( z, dofsPc, doftrans );
				pc.GetSubVector( dofsPc, solPc );
				int const ndofPcS = dofsPc.Size()/neq;

				mU->ComputeElementMatrix( z, muZ, eltrans, ceWork );
				out.elmat.CopyMN( muZ, 0, 0 );
				b->ComputeElementMatrix( z, bZ, eltrans, ceWork );
				out.elmat.CopyMN( bZ, ndofU, 0 );
				bZ.Neg();
				out.elmat.CopyMNt( bZ, 0, ndofU );
				if ( mP )
				{
					mP->ComputeElementMatrix( z, mpZ, eltrans, ceWork );
					out.elmat.CopyMN( mpZ, ndofU, ndofU );
				}
				rhsP.MakeRef( out.rhs, ndofU, ndofP );
				load.AssembleRHSElementVect( *feP, *tr, rhsP );

				int offTr = ndofU + ndofP;
				for ( int f : faces )
				{
					mfem::FiniteElement const *feTr = fesTr.GetFaceElement( f );
					int const ndofTrFs = feTr->GetDof();
					int const ndofTrF = neq*ndofTrFs;
					mesh.GetFaceElementTransformations( f, facetrans, f1, f2 );
					mfem::FaceElementTransformations *ftr = &facetrans;
					mfem::FiniteElement const *feU1 = fesU.GetFE( ftr->Elem1No );
					mfem::FiniteElement const *feU2 = ( ftr->Elem2No >= 0 ) ? fesU.GetFE( ftr->Elem2No ) : feU1;
					cBfi->AssembleFaceMatrix( *feTr, *feU1, *feU2, *ftr, ctF );
					int const offU = ( ftr->Elem1No == z ) ? 0 : feU1->GetDof()*fesU.GetVDim();
					ctFz.CopyMN( ctF, ndofU, ndofTrF, offU, 0 );
					out.elmat.CopyMN( ctFz, 0, offTr );
					out.elmat.CopyMNt( ctFz, offTr, 0 );
					int const side = ( ftr->Elem1No == z ) ? 0 : 1;
					cBfiP->AssembleHDGFaceMatrix( side, *feTr, *feP, *ftr, deghF );
					dFz.CopyMN( deghF, ndofP, ndofP, 0, 0 );
					out.elmat.AddMatrix( dFz, ndofU, ndofU );
					out.elmat.CopyMN( deghF, ndofP, ndofTrF, 0, ndofP, ndofU, offTr );
					out.elmat.CopyMN( deghF, ndofTrF, ndofP, ndofP, 0, offTr, ndofU );
					out.elmat.CopyMN( deghF, ndofTrF, ndofTrF, ndofP, ndofP, offTr, offTr );

					mfem::FiniteElement const *feUt = fesUt.GetFaceElement( f );
					int const ndofUtF = feUt->GetDof();
					int order = feUt->GetOrder() + feTr->GetOrder();
					if ( feTr->GetMapType() != mfem::FiniteElement::VALUE )
						order += ftr->OrderW();
					mfem::IntegrationRule const &ir = mfem::IntRules.Get( feUt->GetGeomType(), order );
					shapeUt.SetSize( ndofUtF );
					shapeTr.SetSize( ndofTrFs );
					rhsF.MakeRef( out.rhs, offTr, ndofTrF );
					rhsF = 0.0;
					fesUt.GetFaceVDofs( f, vdofsUt );
					ut.GetSubVector( vdofsUt, utF );
					for ( int q = 0; q < ir.GetNPoints(); ++q )
					{
						mfem::IntegrationPoint const &ip = ir.IntPoint( q );
						feUt->CalcShape( ip, shapeUt );
						feTr->CalcShape( ip, shapeTr );
						for ( int e = 0; e < neq; ++e )
						{
							mfem::Vector const utFe( utF.GetData() + e*ndofUtF, ndofUtF );
							mfem::real_t const utQ = shapeUt*utFe;
							mfem::real_t wQ = ip.weight*utQ;
							if ( feTr->GetMapType() == mfem::FiniteElement::INTEGRAL )
							{
								ftr->SetIntPoint( &ip );
								wQ /= ftr->Weight();
							}
							mfem::Vector rhsFe( rhsF.GetData() + e*ndofTrFs, ndofTrFs );
							rhsFe.Add( wQ, shapeTr );
						}
					}
					if ( ftr->Elem1No != z )
						rhsF.Neg();
					offTr += ndofTrF;
				}

				{
					mfem::FiniteElement const *fePc = fesPc.GetFE( z );
					int const order = feP->GetOrder() + tr->OrderW();
					mfem::IntegrationRule const &ir = mfem::IntRules.Get( feP->GetGeomType(), order );
					mesh.GetElementTransformation( z, &eltrans );
					shapeP.SetSize( ndofPs );
					shapePc.SetSize( ndofPcS );
					sumPc.SetSize( neq );
					sumPc = 0.0;
					massP.SetSize( ndofPs );
					massP = 0.0;
					for ( int q = 0; q < ir.GetNPoints(); ++q )
					{
						mfem::IntegrationPoint const &ip = ir.IntPoint( q );
						tr->SetIntPoint( &ip );
						fePc->CalcShape( ip, shapePc );
						mfem::real_t const wQ = ip.weight*tr->Weight();
						for ( int e = 0; e < neq; ++e )
						{
							mfem::Vector const pcE( solPc.GetData() + e*ndofPcS, ndofPcS );
							mfem::real_t const value = shapePc*pcE;
							sumPc( e ) += value*wQ;
						}
						feP->CalcShape( ip, shapeP );
						massP.Add( wQ, shapeP );
					}
					for ( int e = 0; e < neq; ++e )
					{
						int const iP = e*ndofPs;
						out.elmat.SetRow( ndofU + iP, 0.0 );
						for ( int i = 0; i < ndofPs; ++i )
							out.elmat( ndofU + iP, ndofU + iP + i ) = massP( i );
						out.rhs( ndofU + iP ) = sumPc( e );
					}
				}
			}

			int elements() const
			{
				return rec.FluxFESpace()->GetNE();
			}

			int potentialSize() const
			{
				return rec.PotentialFESpace()->GetVSize();
			}

			int fluxSize() const
			{
				return rec.FluxFESpace()->GetVSize();
			}

		private:
			mfem::DarcyForm &rec;
			mfem::DarcyHybridization &hyb;
			mfem::GridFunction const &ut;
			mfem::FiniteElementSpace &fesTr;
			mfem::FiniteElementSpace &fesPc;
			mfem::GridFunction pc;
			mfem::VectorDivergenceGridFunctionCoefficient divergence;
			mfem::VectorDomainLFIntegrator load;

			mfem::Array<int> faces, oris, dofsPc, vdofsUt;
			mfem::DofTransformation doftrans;
			mfem::IsoparametricTransformation eltrans, f1, f2;
			mfem::FaceElementTransformations facetrans;
			mfem::DenseMatrix muZ, bZ, mpZ, ceWork, ctF, ctFz, deghF, dFz;
			mfem::Vector solPc, rhsP, rhsF, shapeUt, shapeTr, utF, shapeP, shapePc, sumPc, massP;
	};

	// ------------------------------------------------------- one state
	/// The middle state's global trace and loads, for the split of (b)'s
	/// noise into the part the trace carries and the part the load carries.
	struct Frozen
	{
		std::vector<double> trace, bu, bp;
	};

	State readState( meq::tests::crossing::Solved &solved, bool diagnostics, Frozen const *frozen )
	{
		auto &solver = const_cast<meq::GradShafranovSolver &>( solved.theSolver() );
		mfem::DarcyForm &darcy = *( solver.*get( recoveryaccess::SolverDarcy{} ) );
		mfem::BlockVector &rhs = solver.*get( recoveryaccess::SolverDarcyRhs{} );
		mfem::BlockVector &solution = solver.*get( recoveryaccess::SolverDarcySolution{} );
		mfem::Vector &lambda = solver.*get( recoveryaccess::SolverTraceX{} );
		mfem::DarcyHybridization &hyb = *darcy.GetHybridization();

		auto const &afData = hyb.*get( recoveryaccess::AfData{} );
		auto const &afIpiv = hyb.*get( recoveryaccess::AfIpiv{} );
		auto const &afOffsets = hyb.*get( recoveryaccess::AfOffsets{} );
		auto const &afFOffsets = hyb.*get( recoveryaccess::AfFOffsets{} );
		auto const &bfData = hyb.*get( recoveryaccess::BfData{} );
		auto const &bfOffsets = hyb.*get( recoveryaccess::BfOffsets{} );
		auto const &dfData = hyb.*get( recoveryaccess::DfData{} );
		auto const &dfIpiv = hyb.*get( recoveryaccess::DfIpiv{} );
		auto const &dfOffsets = hyb.*get( recoveryaccess::DfOffsets{} );
		auto const &dfFOffsets = hyb.*get( recoveryaccess::DfFOffsets{} );
		bool const bsym = hyb.*get( recoveryaccess::Bsym{} );
		double const s = bsym ? -1.0 : 1.0;
		auto getFDofs = get( recoveryaccess::GetFDofs{} );
		auto getElementFaces = get( recoveryaccess::GetElementFaces{} );
		auto getCt = get( recoveryaccess::GetCtFaceMatrix{} );
		auto getE = get( recoveryaccess::GetEFaceMatrix{} );
		auto multInv = get( recoveryaccess::MultInv{} );

		mfem::FiniteElementSpace &fesU = *darcy.FluxFESpace();
		mfem::FiniteElementSpace &fesP = *darcy.PotentialFESpace();
		mfem::Mesh &mesh = *fesU.GetMesh();
		int const elements = mesh.GetNE();
		if ( lambda.Size() != hyb.ConstraintFESpace()->GetVSize() )
			throw std::logic_error( "the trace is not the constraint space's size" );

		mfem::Vector const &bu = rhs.GetBlock( 0 );
		mfem::Vector const &bp = rhs.GetBlock( 1 );
		bu.HostRead();
		bp.HostRead();
		lambda.HostRead();
		solution.HostRead();

		State out;
		out.trace = promote( lambda );
		out.buGlobal = promote( bu );
		out.bpGlobal = promote( bp );
		std::size_t const uSize = static_cast<std::size_t>( fesU.GetVSize() );
		std::size_t const pSize = static_cast<std::size_t>( fesP.GetVSize() );
		char const *const methods[] = { "a", "a-check", "b", "b'", "c", "d", "d'", "d2", "d rounded", "e1", "e2",
		                                "f", "f2", "b, trace noise only", "b, load noise only" };
		for ( char const *m : methods )
		{
			out.u[ m ].assign( uSize, 0.0 );
			out.p[ m ].assign( pSize, 0.0 );
		}
		out.u[ "a" ] = promote( solution.GetBlock( 0 ) );
		out.p[ "a" ] = promote( solution.GetBlock( 1 ) );

		std::vector<mfem::DenseMatrix> const dRaw = assemblePotentialBlocks( hyb, fesP );
		mfem::BilinearForm &fluxMass = *darcy.GetFluxMassForm();

		mfem::Array<int> faces, cDofs;
		mfem::Vector lamF, wk, u, p, du, dp;
		mfem::DenseMatrix ct, e;
		for ( int el = 0; el < elements; ++el )
		{
			Element x;
			x.na = afFOffsets[ el + 1 ] - afFOffsets[ el ];
			x.nd = dfFOffsets[ el + 1 ] - dfFOffsets[ el ];
			( hyb.*getFDofs )( el, x.uDofs );
			fesP.GetElementVDofs( el, x.pDofs );

			// ---- the raw blocks, checked against the stored factors
			fluxMass.ComputeElementMatrix( el, x.A );
			x.B.SetSize( x.nd, x.na );
			for ( int j = 0; j < x.na; ++j )
				for ( int i = 0; i < x.nd; ++i )
					x.B( i, j ) = bfData[ bfOffsets[ el ] + i + j*x.nd ];
			x.D = dRaw[ static_cast<std::size_t>( el ) ];
			{
				mfem::DenseMatrix a( x.A ), sc( x.D );
				mfem::Array<int> ipivA( x.na ), ipivS( x.nd );
				mfem::LUFactors luA( a.Data(), ipivA.GetData() );
				luA.Factor( x.na );
				mfem::DenseMatrix aiBt( x.na, x.nd );
				aiBt.Transpose( x.B );
				if ( !bsym )
					aiBt.Neg();
				luA.Solve( x.na, x.nd, aiBt.Data() );
				mfem::AddMult( x.B, aiBt, sc );
				mfem::LUFactors luS( sc.Data(), ipivS.GetData() );
				luS.Factor( x.nd );
				for ( int i = 0; i < x.na*x.na; ++i )
					out.factorADifference = std::max( out.factorADifference,
						std::abs( static_cast<double>( a.Data()[ i ] ) - afData[ afOffsets[ el ] + i ] ) );
				for ( int i = 0; i < x.na; ++i )
					if ( ipivA[ i ] != afIpiv[ afFOffsets[ el ] + i ] )
						out.factorADifference = std::max( out.factorADifference, 1.0 );
				for ( int i = 0; i < x.nd*x.nd; ++i )
					out.factorSDifference = std::max( out.factorSDifference,
						std::abs( static_cast<double>( sc.Data()[ i ] ) - dfData[ dfOffsets[ el ] + i ] ) );
				for ( int i = 0; i < x.nd; ++i )
					if ( ipivS[ i ] != dfIpiv[ dfFOffsets[ el ] + i ] )
						out.factorSDifference = std::max( out.factorSDifference, 1.0 );
			}

			// ---- the local right-hand side, as ComputeSolution() forms it
			bu.GetSubVector( x.uDofs, x.buReal );
			bp.GetSubVector( x.pDofs, x.bpReal );
			if ( bsym )
				x.bpReal.Neg();
			x.buLoad = promote( x.buReal );
			x.bpLoad = promote( x.bpReal );
			x.buDouble = x.buLoad;
			x.bpDouble = x.bpLoad;
			x.buTrace.assign( static_cast<std::size_t>( x.na ), 0.0 );
			x.bpTrace.assign( static_cast<std::size_t>( x.nd ), 0.0 );
			// (b) split: the load frozen at the middle state (trace noise only),
			// and the trace frozen there (load noise only).
			std::vector<double> buTraceOnly, bpTraceOnly, buLoadOnly, bpLoadOnly;
			if ( frozen )
			{
				for ( int i = 0; i < x.na; ++i )
				{
					buTraceOnly.push_back( frozen->bu[ static_cast<std::size_t>( x.uDofs[ i ] ) ] );
					buLoadOnly.push_back( x.buLoad[ static_cast<std::size_t>( i ) ] );
				}
				for ( int i = 0; i < x.nd; ++i )
				{
					bpTraceOnly.push_back( ( bsym ? -1.0 : 1.0 )*frozen->bp[ static_cast<std::size_t>( x.pDofs[ i ] ) ] );
					bpLoadOnly.push_back( x.bpLoad[ static_cast<std::size_t>( i ) ] );
				}
			}
			( hyb.*getElementFaces )( el, faces );
			for ( int f = 0; f < faces.Size(); ++f )
			{
				int el1 = -1, el2 = -1;
				mesh.GetFaceElements( faces[ f ], &el1, &el2 );
				( hyb.*getCt )( faces[ f ], el1 != el, ct );
				hyb.TraceVDofs( faces[ f ], cDofs );
				lambda.GetSubVector( cDofs, lamF );
				ct.AddMult_a( -1.0, lamF, x.buReal );
				( hyb.*getE )( faces[ f ], el1 != el, e );
				e.AddMult_a( -1.0, lamF, x.bpReal );
				for ( int c = 0; c < cDofs.Size(); ++c )
				{
					double const l = static_cast<double>( lamF( c ) );
					double const lFrozen = frozen ? frozen->trace[ static_cast<std::size_t>( cDofs[ c ] ) ] : 0.0;
					if ( frozen )
					{
						for ( int i = 0; i < x.na; ++i )
						{
							buTraceOnly[ static_cast<std::size_t>( i ) ] -= ct( i, c )*l;
							buLoadOnly[ static_cast<std::size_t>( i ) ] -= ct( i, c )*lFrozen;
						}
						for ( int i = 0; i < x.nd; ++i )
						{
							bpTraceOnly[ static_cast<std::size_t>( i ) ] -= e( i, c )*l;
							bpLoadOnly[ static_cast<std::size_t>( i ) ] -= e( i, c )*lFrozen;
						}
					}
					for ( int i = 0; i < x.na; ++i )
						x.buTrace[ static_cast<std::size_t>( i ) ] -= ct( i, c )*l;
					for ( int i = 0; i < x.nd; ++i )
						x.bpTrace[ static_cast<std::size_t>( i ) ] -= e( i, c )*l;
				}
			}
			for ( int i = 0; i < x.na; ++i )
				x.buDouble[ static_cast<std::size_t>( i ) ] += x.buTrace[ static_cast<std::size_t>( i ) ];
			for ( int i = 0; i < x.nd; ++i )
				x.bpDouble[ static_cast<std::size_t>( i ) ] += x.bpTrace[ static_cast<std::size_t>( i ) ];

			// ---- (a) re-run through MultInv(): must be RecoverFEMSolution()'s bits
			( hyb.*multInv )( el, x.buReal, x.bpReal, u, p, false, &wk );
			scatter( out.u[ "a-check" ], x.uDofs, u );
			scatter( out.p[ "a-check" ], x.pDofs, p );
			for ( int i = 0; i < x.na; ++i )
				out.recoveryUDifference = std::max( out.recoveryUDifference,
					std::abs( static_cast<double>( u( i ) ) - out.u[ "a" ][ static_cast<std::size_t>( x.uDofs[ i ] ) ] ) );
			for ( int i = 0; i < x.nd; ++i )
				out.recoveryPDifference = std::max( out.recoveryPDifference,
					std::abs( static_cast<double>( p( i ) ) - out.p[ "a" ][ static_cast<std::size_t>( x.pDofs[ i ] ) ] ) );
			mfem::Vector const uA( u ), pA( p );

			// ---- (b), (b'): LU in double on the stored data
			Dense const m = fullMatrix( x, s );
			std::vector<double> const bD = concat( x.buDouble, x.bpDouble );
			std::vector<double> const bR = concat( promote( x.buReal ), promote( x.bpReal ) );
			std::vector<double> const xb = luSolve( m, bD );
			std::vector<double> const xb1 = luSolve( m, bR );
			scatter( out.u[ "b" ], x.uDofs, xb );
			scatter( out.p[ "b" ], x.pDofs, xb, static_cast<std::size_t>( x.na ) );
			scatter( out.u[ "b'" ], x.uDofs, xb1 );
			scatter( out.p[ "b'" ], x.pDofs, xb1, static_cast<std::size_t>( x.na ) );
			if ( frozen )
			{
				std::vector<double> const xt = luSolve( m, concat( buTraceOnly, bpTraceOnly ) );
				scatter( out.u[ "b, trace noise only" ], x.uDofs, xt );
				scatter( out.p[ "b, trace noise only" ], x.pDofs, xt, static_cast<std::size_t>( x.na ) );
				std::vector<double> const xl = luSolve( m, concat( buLoadOnly, bpLoadOnly ) );
				scatter( out.u[ "b, load noise only" ], x.uDofs, xl );
				scatter( out.p[ "b, load noise only" ], x.pDofs, xl, static_cast<std::size_t>( x.na ) );
			}
			// ---- (c): one pass, residual in real_t
			{
				mfem::Vector ru( x.buReal ), rp( x.bpReal );
				x.A.AddMult_a( -1.0, uA, ru );
				x.B.AddMultTranspose_a( static_cast<mfem::real_t>( -s ), pA, ru );
				x.B.AddMult_a( -1.0, uA, rp );
				x.D.AddMult_a( -1.0, pA, rp );
				( hyb.*multInv )( el, ru, rp, du, dp, false, &wk );
				mfem::Vector uc( uA ), pc( pA );
				uc += du;
				pc += dp;
				scatter( out.u[ "c" ], x.uDofs, uc );
				scatter( out.p[ "c" ], x.pDofs, pc );
			}

			// ---- (d), (d'), (d2): residual in double, correction on the stored factors
			auto refine = [ & ]( std::vector<double> const &xIn, std::vector<double> const &b )
			{
				std::vector<double> r = b;
				std::vector<double> const mx = multiply( m, xIn );
				for ( std::size_t i = 0; i < r.size(); ++i )
					r[ i ] -= mx[ i ];
				mfem::Vector ru( x.na ), rp( x.nd );
				for ( int i = 0; i < x.na; ++i )
					ru( i ) = static_cast<mfem::real_t>( r[ static_cast<std::size_t>( i ) ] );
				for ( int i = 0; i < x.nd; ++i )
					rp( i ) = static_cast<mfem::real_t>( r[ static_cast<std::size_t>( x.na + i ) ] );
				( hyb.*multInv )( el, ru, rp, du, dp, false, &wk );
				std::vector<double> xOut = xIn;
				for ( int i = 0; i < x.na; ++i )
					xOut[ static_cast<std::size_t>( i ) ] += static_cast<double>( du( i ) );
				for ( int i = 0; i < x.nd; ++i )
					xOut[ static_cast<std::size_t>( x.na + i ) ] += static_cast<double>( dp( i ) );
				return xOut;
			};
			std::vector<double> const x0 = concat( promote( uA ), promote( pA ) );
			std::vector<double> const xd = refine( x0, bD );
			std::vector<double> const xd1 = refine( x0, bR );
			std::vector<double> const xd2 = refine( xd, bD );
			std::vector<double> const xdr = roundedToReal( xd );
			scatter( out.u[ "d" ], x.uDofs, xd );
			scatter( out.p[ "d" ], x.pDofs, xd, static_cast<std::size_t>( x.na ) );
			scatter( out.u[ "d'" ], x.uDofs, xd1 );
			scatter( out.p[ "d'" ], x.pDofs, xd1, static_cast<std::size_t>( x.na ) );
			scatter( out.u[ "d2" ], x.uDofs, xd2 );
			scatter( out.p[ "d2" ], x.pDofs, xd2, static_cast<std::size_t>( x.na ) );
			scatter( out.u[ "d rounded" ], x.uDofs, xdr );
			scatter( out.p[ "d rounded" ], x.pDofs, xdr, static_cast<std::size_t>( x.na ) );

			// ---- (e1), (e2): the Schur recovery in real_t on scaled blocks
			Dense const aD = promote( x.A );
			Dense sD = promote( x.D );       // S = D - s B A^-1 B^T, in double
			{
				Dense bt( x.na, x.nd );
				for ( int j = 0; j < x.nd; ++j )
					for ( int i = 0; i < x.na; ++i )
						bt( i, j ) = x.B( j, i );
				for ( int j = 0; j < x.nd; ++j )
				{
					std::vector<double> col( bt.a.begin() + j*x.na, bt.a.begin() + ( j + 1 )*x.na );
					std::vector<double> const aiCol = luSolve( aD, col );
					for ( int i = 0; i < x.nd; ++i )
					{
						double sum = 0.0;
						for ( int k = 0; k < x.na; ++k )
							sum += x.B( i, k )*aiCol[ static_cast<std::size_t>( k ) ];
						sD( i, j ) -= s*sum;
					}
				}
			}
			std::vector<double> blockScale( static_cast<std::size_t>( x.na + x.nd ) );
			{
				double aMax = 0.0, sMax = 0.0;
				for ( int i = 0; i < x.na; ++i )
					aMax = std::max( aMax, std::abs( aD( i, i ) ) );
				for ( int i = 0; i < x.nd; ++i )
					sMax = std::max( sMax, std::abs( sD( i, i ) ) );
				for ( int i = 0; i < x.na; ++i )
					blockScale[ static_cast<std::size_t>( i ) ] = powerOfTwo( 1.0/std::sqrt( aMax ) );
				for ( int i = 0; i < x.nd; ++i )
					blockScale[ static_cast<std::size_t>( x.na + i ) ] = powerOfTwo( 1.0/std::sqrt( sMax ) );
			}
			std::vector<double> const ruizScale = ruiz( m );
			scaledSchurReal( x, bsym, blockScale, u, p );
			scatter( out.u[ "e1" ], x.uDofs, u );
			scatter( out.p[ "e1" ], x.pDofs, p );
			scaledSchurReal( x, bsym, ruizScale, u, p );
			scatter( out.u[ "e2" ], x.uDofs, u );
			scatter( out.p[ "e2" ], x.pDofs, p );

			// ---- (f), (f2): one LU of the full local matrix in real_t
			fullLuReal( x, s, nullptr, u, p );
			scatter( out.u[ "f" ], x.uDofs, u );
			scatter( out.p[ "f" ], x.pDofs, p );
			fullLuReal( x, s, &ruizScale, u, p );
			scatter( out.u[ "f2" ], x.uDofs, u );
			scatter( out.p[ "f2" ], x.pDofs, p );

			if ( !diagnostics )
				continue;

			// ---- per-element diagnostics, in double from the stored data
			auto &dg = out.diagnostics;
			dg[ "cond(A)" ].push_back( cond2( aD ) );
			dg[ "cond(S)" ].push_back( cond2( sD ) );
			dg[ "cond(M)" ].push_back( cond2( m ) );
			dg[ "cond(M), Ruiz" ].push_back( cond2( scaled( m, ruizScale ) ) );
			// cancellation inside the solve: B A^-1 b_u against b_p
			std::vector<double> const aiBu = luSolve( aD, x.buDouble );
			std::vector<double> baiBu( static_cast<std::size_t>( x.nd ), 0.0 );
			for ( int i = 0; i < x.nd; ++i )
				for ( int k = 0; k < x.na; ++k )
					baiBu[ static_cast<std::size_t>( i ) ] += x.B( i, k )*aiBu[ static_cast<std::size_t>( k ) ];
			std::vector<double> diff( static_cast<std::size_t>( x.nd ) );
			for ( int i = 0; i < x.nd; ++i )
				diff[ static_cast<std::size_t>( i ) ] = baiBu[ static_cast<std::size_t>( i ) ]
				                                        - x.bpDouble[ static_cast<std::size_t>( i ) ];
			double const nBaiBu = normInf( baiBu ), nBp = normInf( x.bpDouble ), nDiff = normInf( diff );
			{
				// the u update u = A^-1 b_u - s A^-1 B^T p, at the double answer
				std::vector<double> btp( static_cast<std::size_t>( x.na ), 0.0 );
				for ( int k = 0; k < x.na; ++k )
					for ( int i = 0; i < x.nd; ++i )
						btp[ static_cast<std::size_t>( k ) ] += s*x.B( i, k )*xb[ static_cast<std::size_t>( x.na + i ) ];
				std::vector<double> const aiBtp = luSolve( aD, btp );
				std::vector<double> ub( xb.begin(), xb.begin() + x.na );
				dg[ "|A^-1 b_u|" ].push_back( normInf( aiBu ) );
				dg[ "|A^-1 s B^T p|" ].push_back( normInf( aiBtp ) );
				dg[ "|u|" ].push_back( normInf( ub ) );
				dg[ "|p|" ].push_back( normInf( std::vector<double>( xb.begin() + x.na, xb.end() ) ) );
				dg[ "cancellation in the u update" ].push_back( normInf( aiBu )/std::max( normInf( ub ), 1.0e-300 ) );
				// ||B A^-1|| ||s B^T|| / ||S||: the growth of the Schur route
				Dense ai( x.na, x.na );
				for ( int j = 0; j < x.na; ++j )
				{
					std::vector<double> ej( static_cast<std::size_t>( x.na ), 0.0 );
					ej[ static_cast<std::size_t>( j ) ] = 1.0;
					std::vector<double> const col = luSolve( aD, ej );
					for ( int i = 0; i < x.na; ++i )
						ai( i, j ) = col[ static_cast<std::size_t>( i ) ];
				}
				Dense bai( x.nd, x.na );
				for ( int j = 0; j < x.na; ++j )
					for ( int i = 0; i < x.nd; ++i )
						for ( int k = 0; k < x.na; ++k )
							bai( i, j ) += x.B( i, k )*ai( k, j );
				Dense bD2 = promote( x.B );
				dg[ "|B A^-1| |B| / |S|" ].push_back( normInf( bai )*normInf( bD2 )/normInf( sD ) );
				dg[ "|B A^-1| |B| / |D|" ].push_back( normInf( bai )*normInf( bD2 )/normInf( promote( x.D ) ) );
				// backward error of the full-LU route (f)
				mfem::Vector uf, pf;
				fullLuReal( x, s, nullptr, uf, pf );
				std::vector<double> const xf = concat( promote( uf ), promote( pf ) );
				std::vector<double> rf = bR;
				std::vector<double> const mxf = multiply( m, xf );
				for ( std::size_t i = 0; i < rf.size(); ++i )
					rf[ i ] -= mxf[ i ];
				dg[ "backward error of (f), eps" ].push_back(
					normInf( rf )/( normInf( m )*normInf( xf ) + normInf( bR ) )/eps );
			}
			dg[ "|B A^-1 b_u|" ].push_back( nBaiBu );
			dg[ "|b_p|" ].push_back( nBp );
			dg[ "|B A^-1 b_u - b_p|" ].push_back( nDiff );
			dg[ "cancellation in the solve" ].push_back( std::max( nBaiBu, nBp )/std::max( nDiff, 1.0e-300 ) );
			// cancellation forming the right-hand side: load against C^T lambda, E lambda
			dg[ "|b_u load|" ].push_back( normInf( x.buLoad ) );
			dg[ "|C^T lambda|" ].push_back( normInf( x.buTrace ) );
			dg[ "|b_u|" ].push_back( normInf( x.buDouble ) );
			dg[ "|b_p load|" ].push_back( normInf( x.bpLoad ) );
			dg[ "|E lambda|" ].push_back( normInf( x.bpTrace ) );
			dg[ "cancellation forming b_p" ].push_back(
				std::max( normInf( x.bpLoad ), normInf( x.bpTrace ) )/std::max( normInf( x.bpDouble ), 1.0e-300 ) );
			dg[ "cancellation forming b_u" ].push_back(
				std::max( normInf( x.buLoad ), normInf( x.buTrace ) )/std::max( normInf( x.buDouble ), 1.0e-300 ) );
			// the residual of MFEM's recovered (u, p), in double, normwise backward error in eps
			{
				std::vector<double> r = bD;
				std::vector<double> const mx = multiply( m, x0 );
				for ( std::size_t i = 0; i < r.size(); ++i )
					r[ i ] -= mx[ i ];
				dg[ "backward error of (a), eps" ].push_back(
					normInf( r )/( normInf( m )*normInf( x0 ) + normInf( bD ) )/eps );
				std::vector<double> rp( r.begin() + x.na, r.end() );
				dg[ "residual of (a), potential rows / |b_p|, eps" ].push_back( normInf( rp )/nBp/eps );
				std::vector<double> ua( x0.begin(), x0.begin() + x.na ), ub( xb.begin(), xb.begin() + x.na );
				std::vector<double> pa( x0.begin() + x.na, x0.end() ), pb( xb.begin() + x.na, xb.end() );
				dg[ "forward error of (a) in p, of |p|_el, eps" ].push_back( relativeDistance( pa, pb )/eps );
				dg[ "forward error of (a) in u, of |u|_el, eps" ].push_back( relativeDistance( ua, ub )/eps );
				std::vector<double> rb = bD;
				std::vector<double> const mxb = multiply( m, xb );
				for ( std::size_t i = 0; i < rb.size(); ++i )
					rb[ i ] -= mxb[ i ];
				dg[ "backward error of (b), eps" ].push_back(
					normInf( rb )/( normInf( m )*normInf( xb ) + normInf( bD ) )/eps );
			}
		}

		// ---------------------------------------- Reconstruct()
		mfem::GridFunction &totalFlux = solver.*get( recoveryaccess::SolverTotalFlux{} );
		mfem::GridFunction &enrichedTrace = solver.*get( recoveryaccess::SolverEnrichedTrace{} );
		mfem::GridFunction &postProcessed = solver.*get( recoveryaccess::SolverPostProcessed{} );
		mfem::GridFunction &enrichedFlux = solver.*get( recoveryaccess::SolverEnrichedFlux{} );
		out.psiStar = promote( postProcessed );
		out.totalFlux = promote( totalFlux );
		out.enrichedFlux = promote( enrichedFlux );
		Reconstructor reconstructor( darcy, hyb, totalFlux, enrichedTrace, solution );
		std::size_t const psSize = static_cast<std::size_t>( reconstructor.potentialSize() );
		char const *const rMethods[] = { "a-check", "b", "c", "d", "e2" };
		for ( char const *mth : rMethods )
			out.psiStarMethod[ mth ].assign( psSize, 0.0 );
		ReconstructionProblem problem;
		mfem::DenseMatrixInverse inverse;
		mfem::Vector sol, residual, correction;
		for ( int z = 0; z < reconstructor.elements(); ++z )
		{
			reconstructor.assemble( z, problem );
			int const n = problem.elmat.Height();
			int const pOff = problem.ndofU;
			inverse.Factor( problem.elmat );
			inverse.Mult( problem.rhs, sol );
			for ( int i = 0; i < problem.ndofP; ++i )
			{
				double const mine = static_cast<double>( sol( pOff + i ) );
				out.psiStarMethod[ "a-check" ][ static_cast<std::size_t>( problem.dofsP[ i ] ) ] = mine;
				out.reconstructionDifference = std::max( out.reconstructionDifference,
					std::abs( mine - out.psiStar[ static_cast<std::size_t>( problem.dofsP[ i ] ) ] ) );
			}
			Dense const m = promote( problem.elmat );
			std::vector<double> const b = promote( problem.rhs );
			std::vector<double> const xb = luSolve( m, b );
			scatter( out.psiStarMethod[ "b" ], problem.dofsP, xb, static_cast<std::size_t>( pOff ) );
			// (c): residual in real_t, correction on the same factors
			residual = problem.rhs;
			problem.elmat.AddMult_a( -1.0, sol, residual );
			inverse.Mult( residual, correction );
			mfem::Vector solC( sol );
			solC += correction;
			for ( int i = 0; i < problem.ndofP; ++i )
				out.psiStarMethod[ "c" ][ static_cast<std::size_t>( problem.dofsP[ i ] ) ]
					= static_cast<double>( solC( pOff + i ) );
			// (d): residual in double
			std::vector<double> const x0 = promote( sol );
			std::vector<double> r = b;
			std::vector<double> const mx = multiply( m, x0 );
			for ( std::size_t i = 0; i < r.size(); ++i )
				r[ i ] -= mx[ i ];
			inverse.Mult( demote( r ), correction );
			for ( int i = 0; i < problem.ndofP; ++i )
				out.psiStarMethod[ "d" ][ static_cast<std::size_t>( problem.dofsP[ i ] ) ]
					= x0[ static_cast<std::size_t>( pOff + i ) ] + static_cast<double>( correction( pOff + i ) );
			// (e2): Ruiz-equilibrated, in real_t
			std::vector<double> const d = ruiz( m );
			mfem::DenseMatrix ms = demote( scaled( m, d ) );
			mfem::Vector bs( n );
			for ( int i = 0; i < n; ++i )
				bs( i ) = problem.rhs( i )*static_cast<mfem::real_t>( d[ static_cast<std::size_t>( i ) ] );
			mfem::DenseMatrixInverse inverseScaled( ms );
			mfem::Vector xs;
			inverseScaled.Mult( bs, xs );
			for ( int i = 0; i < problem.ndofP; ++i )
				out.psiStarMethod[ "e2" ][ static_cast<std::size_t>( problem.dofsP[ i ] ) ]
					= static_cast<double>( xs( pOff + i )*static_cast<mfem::real_t>(
						d[ static_cast<std::size_t>( pOff + i ) ] ) );
			if ( diagnostics )
			{
				auto &dg = out.reconstructionDiagnostics;
				dg[ "cond(K)" ].push_back( cond2( m ) );
				dg[ "cond(K), Ruiz" ].push_back( cond2( scaled( m, d ) ) );
				dg[ "backward error of MFEM's solve, eps" ].push_back(
					normInf( r )/( normInf( m )*normInf( x0 ) + normInf( b ) )/eps );
				std::vector<double> pa( x0.begin() + pOff, x0.begin() + pOff + problem.ndofP );
				std::vector<double> pb( xb.begin() + pOff, xb.begin() + pOff + problem.ndofP );
				dg[ "forward error in p*, of |p*|_el, eps" ].push_back( relativeDistance( pa, pb )/eps );
				double rowMin = std::numeric_limits<double>::infinity(), rowMax = 0.0;
				for ( int i = 0; i < n; ++i )
				{
					double row = 0.0;
					for ( int j = 0; j < n; ++j )
						row = std::max( row, std::abs( m( i, j ) ) );
					rowMin = std::min( rowMin, row );
					rowMax = std::max( rowMax, row );
				}
				dg[ "row scale spread" ].push_back( rowMax/rowMin );
			}
			// The local problem's load by part: the divergence rows (less the
			// row the mean constraint replaces), the face rows, the mean.
			for ( int i = 1; i < problem.ndofP; ++i )
				out.reconstructionLoad[ "divergence rows" ].push_back( b[ static_cast<std::size_t>( pOff + i ) ] );
			for ( int i = pOff + problem.ndofP; i < n; ++i )
				out.reconstructionLoad[ "face rows" ].push_back( b[ static_cast<std::size_t>( i ) ] );
			out.reconstructionLoad[ "mean row" ].push_back( b[ static_cast<std::size_t>( pOff ) ] );
		}

		/*
		 * WHICH INPUT CARRIES psi*'s FLOOR. The total flux is re-made by
		 * MFEM's own ReconstructTotalFlux() -- checked against the one
		 * Reconstruct() made -- from fields recovered another way, and the
		 * local problem is then solved on it in double.
		 */
		auto pipeline = [ & ]( std::vector<double> const &uIn, std::vector<double> const &pIn, char const *name,
		                       bool keepFlux )
		{
			mfem::BlockVector sol( darcy.GetOffsets() );
			for ( int i = 0; i < sol.GetBlock( 0 ).Size(); ++i )
				sol.GetBlock( 0 )( i ) = static_cast<mfem::real_t>( uIn[ static_cast<std::size_t>( i ) ] );
			for ( int i = 0; i < sol.GetBlock( 1 ).Size(); ++i )
				sol.GetBlock( 1 )( i ) = static_cast<mfem::real_t>( pIn[ static_cast<std::size_t>( i ) ] );
			mfem::GridFunction ut( totalFlux.FESpace() );
			darcy.ReconstructTotalFlux( sol, lambda, ut );
			if ( keepFlux )
				out.totalFluxMethod[ name ] = promote( ut );
			Reconstructor r( darcy, hyb, ut, enrichedTrace, sol );
			std::vector<double> psReal( psSize, 0.0 ), psDouble( psSize, 0.0 );
			ReconstructionProblem pr;
			mfem::DenseMatrixInverse inv;
			mfem::Vector x;
			for ( int z = 0; z < r.elements(); ++z )
			{
				r.assemble( z, pr );
				inv.Factor( pr.elmat );
				inv.Mult( pr.rhs, x );
				for ( int i = 0; i < pr.ndofP; ++i )
					psReal[ static_cast<std::size_t>( pr.dofsP[ i ] ) ] = static_cast<double>( x( pr.ndofU + i ) );
				std::vector<double> const xb = luSolve( promote( pr.elmat ), promote( pr.rhs ) );
				scatter( psDouble, pr.dofsP, xb, static_cast<std::size_t>( pr.ndofU ) );
			}
			out.psiStarMethod[ std::string( name ) + ", real_t solve" ] = psReal;
			out.psiStarMethod[ std::string( name ) + ", double solve" ] = psDouble;
		};
		{
			mfem::GridFunction utCheck( totalFlux.FESpace() );
			darcy.ReconstructTotalFlux( solution, lambda, utCheck );
			for ( int i = 0; i < utCheck.Size(); ++i )
				out.totalFluxDifference = std::max( out.totalFluxDifference,
					std::abs( static_cast<double>( utCheck( i ) - totalFlux( i ) ) ) );
		}
		pipeline( out.u[ "a" ], out.p[ "a" ], "u_t and p_c from (a)", false );
		pipeline( out.u[ "d rounded" ], out.p[ "d rounded" ], "u_t and p_c from (d) rounded", true );
		pipeline( out.u[ "a" ], out.p[ "d rounded" ], "u_t from (a), p_c from (d) rounded", false );
		pipeline( out.u[ "d rounded" ], out.p[ "a" ], "u from (d) rounded, p from (a)", true );

		/*
		 * THE LOCAL PROBLEM'S LOAD, FORMED IN DOUBLE. It is linear in the total
		 * flux's coefficients and in the potential's, with no constant part, so
		 * it is G u_t + P p_c with G and P fixed per element: their columns are
		 * assembled here by the same real_t code, from unit coefficients, and
		 * being fixed they cancel from a second difference as the local matrix
		 * does. The products are then taken in double, which removes the
		 * per-state rounding of the load's assembly and nothing else.
		 */
		{
			mfem::GridFunction utScratch( totalFlux.FESpace() );
			utScratch = 0.0;
			mfem::BlockVector solScratch( darcy.GetOffsets() );
			solScratch = 0.0;
			Reconstructor rs( darcy, hyb, utScratch, enrichedTrace, solScratch );
			mfem::FiniteElementSpace const &fesUt = *totalFlux.FESpace();
			std::vector<double> const &utD = out.totalFluxMethod.at( "u_t and p_c from (d) rounded" );
			std::vector<double> const &pcD = out.p.at( "d rounded" );
			std::vector<double> const utA = promote( totalFlux );
			std::vector<double> const &pcA = out.p.at( "a" );
			std::vector<double> psA( psSize, 0.0 ), psD( psSize, 0.0 ), psDReal( psSize, 0.0 );
			std::vector<double> psDivOnly( psSize, 0.0 ), psFaceOnly( psSize, 0.0 );
			ReconstructionProblem pr, column;
			mfem::Array<int> utDofs, pcDofs;
			mfem::DofTransformation doftrans;
			for ( int z = 0; z < reconstructor.elements(); ++z )
			{
				reconstructor.assemble( z, pr );
				int const n = pr.elmat.Height();
				const_cast<mfem::FiniteElementSpace &>( fesUt ).GetElementVDofs( z, utDofs, doftrans );
				darcy.PotentialFESpace()->GetElementVDofs( z, pcDofs );
				std::vector<double> loadA( static_cast<std::size_t>( n ), 0.0 ), loadD = loadA;
				for ( int idx : utDofs )
				{
					int const i = idx >= 0 ? idx : -1 - idx;
					utScratch( i ) = 1.0;
					rs.assemble( z, column );
					utScratch( i ) = 0.0;
					for ( int r = 0; r < n; ++r )
					{
						loadA[ static_cast<std::size_t>( r ) ] += column.rhs( r )*utA[ static_cast<std::size_t>( i ) ];
						loadD[ static_cast<std::size_t>( r ) ] += column.rhs( r )*utD[ static_cast<std::size_t>( i ) ];
					}
				}
				for ( int idx : pcDofs )
				{
					solScratch.GetBlock( 1 )( idx ) = 1.0;
					rs.assemble( z, column );
					solScratch.GetBlock( 1 )( idx ) = 0.0;
					for ( int r = 0; r < n; ++r )
					{
						loadA[ static_cast<std::size_t>( r ) ] += column.rhs( r )*pcA[ static_cast<std::size_t>( idx ) ];
						loadD[ static_cast<std::size_t>( r ) ] += column.rhs( r )*pcD[ static_cast<std::size_t>( idx ) ];
					}
				}
				double scale = 0.0, gap = 0.0;
				for ( int r = 0; r < n; ++r )
				{
					scale = std::max( scale, std::abs( static_cast<double>( pr.rhs( r ) ) ) );
					gap = std::max( gap, std::abs( static_cast<double>( pr.rhs( r ) ) - loadA[ static_cast<std::size_t>( r ) ] ) );
				}
				out.linearisedLoadDifference = std::max( out.linearisedLoadDifference, gap/scale );
				Dense const m = promote( pr.elmat );
				scatter( psA, pr.dofsP, luSolve( m, loadA ), static_cast<std::size_t>( pr.ndofU ) );
				// which rows: only the divergence rows from the double load, or only the face rows
				{
					std::vector<double> divOnly = promote( pr.rhs ), faceOnly = promote( pr.rhs );
					for ( int r = pr.ndofU; r < pr.ndofU + pr.ndofP; ++r )
						divOnly[ static_cast<std::size_t>( r ) ] = loadA[ static_cast<std::size_t>( r ) ];
					for ( int r = pr.ndofU + pr.ndofP; r < n; ++r )
						faceOnly[ static_cast<std::size_t>( r ) ] = loadA[ static_cast<std::size_t>( r ) ];
					scatter( psDivOnly, pr.dofsP, luSolve( m, divOnly ), static_cast<std::size_t>( pr.ndofU ) );
					scatter( psFaceOnly, pr.dofsP, luSolve( m, faceOnly ), static_cast<std::size_t>( pr.ndofU ) );
				}
				scatter( psD, pr.dofsP, luSolve( m, loadD ), static_cast<std::size_t>( pr.ndofU ) );
				// and with the load rounded back to real_t, solved in double
				std::vector<double> const loadDReal = roundedToReal( loadD );
				scatter( psDReal, pr.dofsP, luSolve( m, loadDReal ), static_cast<std::size_t>( pr.ndofU ) );
			}
			out.psiStarMethod[ "u_t and p_c from (a), load in double, double solve" ] = psA;
			out.psiStarMethod[ "u_t and p_c from (a), divergence rows in double, double solve" ] = psDivOnly;
			out.psiStarMethod[ "u_t and p_c from (a), face rows in double, double solve" ] = psFaceOnly;
			out.psiStarMethod[ "u_t and p_c from (d) rounded, load in double, double solve" ] = psD;
			out.psiStarMethod[ "u_t and p_c from (d) rounded, load in double rounded, double solve" ] = psDReal;
		}
		return out;
	}

	void printDiagnostics( char const *title, std::map<std::string, std::vector<double>> const &dg )
	{
		std::printf( "    %s, per element (worst / median):\n", title );
		for ( auto const &[ name, values ] : dg )
			std::printf( "      %-48s %10.3e / %10.3e\n", name.c_str(), maxOf( values ), median( values ) );
	}

	void measure( int n, int order )
	{
		Frozen frozen;
		{
			meq::tests::crossing::Solved solved( centre + ( states/2 )*spacing, n, order );
			State const s = readState( solved, false, nullptr );
			frozen.trace = s.trace;
			frozen.bu = s.buGlobal;
			frozen.bp = s.bpGlobal;
		}
		std::vector<State> series;
		for ( int j = 0; j < states; ++j )
		{
			meq::tests::crossing::Solved solved( centre + j*spacing, n, order );
			bool const middle = j == states/2;
			series.push_back( readState( solved, middle, &frozen ) );
		}

		std::printf( "\n  box, n = %d, k = %d, %s (eps = %.3e), %d states %.0e apart in A\n", n, order,
		             sizeof( mfem::real_t ) == 4 ? "single" : "double", eps, states, spacing );
		double fa = 0.0, fs = 0.0, ru = 0.0, rp = 0.0, rr = 0.0, rt = 0.0, rl = 0.0;
		for ( State const &s : series )
		{
			fa = std::max( fa, s.factorADifference );
			fs = std::max( fs, s.factorSDifference );
			ru = std::max( ru, s.recoveryUDifference );
			rp = std::max( rp, s.recoveryPDifference );
			rr = std::max( rr, s.reconstructionDifference );
			rt = std::max( rt, s.totalFluxDifference );
			rl = std::max( rl, s.linearisedLoadDifference );
		}
		std::printf( "    checks, worst over states and elements (0 is bit for bit): re-assembled A factored "
		             "%.1e, S factored %.1e; MultInv() on the re-formed load u %.1e, p %.1e; Reconstruct()'s "
		             "local problem psi* %.1e; ReconstructTotalFlux() re-run %.1e; the load as G u_t + P p_c, "
		             "of its largest entry, %.1e\n", fa, fs, ru, rp, rr, rt, rl );

		auto column = [ & ]( auto pick )
		{
			std::vector<std::vector<double>> v;
			for ( State const &s : series )
				v.push_back( pick( s ) );
			return sigma( v );
		};
		Sigma const tr = column( []( State const &s ) { return s.trace; } );
		Sigma const bpG = column( []( State const &s ) { return s.bpGlobal; } );
		std::printf( "    sigma, of the largest value, rms / worst: trace %.2e / %.2e (%.2f / %.2f eps); "
		             "load b_p %.2e / %.2e\n", tr.rms, tr.worst, tr.rms/eps, tr.worst/eps, bpG.rms, bpG.worst );
		std::printf( "    %-34s %-23s %-23s %s\n", "recovery", "p = psi_h", "u = -q_h", "of the trace's rms" );
		for ( auto const &[ name, unused ] : series.front().u )
		{
			std::string const key = name;
			Sigma const sp = column( [ & ]( State const &s ) { return s.p.at( key ); } );
			Sigma const su = column( [ & ]( State const &s ) { return s.u.at( key ); } );
			std::printf( "    %-34s %.2e / %.2e     %.2e / %.2e     p %.0fx, u %.0fx\n", name.c_str(), sp.rms,
			             sp.worst, su.rms, su.worst, sp.rms/tr.rms, su.rms/tr.rms );
		}
		Sigma const ut = column( []( State const &s ) { return s.totalFlux; } );
		Sigma const uStar = column( []( State const &s ) { return s.enrichedFlux; } );
		Sigma const ps = column( []( State const &s ) { return s.psiStar; } );
		std::printf( "    Reconstruct(): total flux u_t %.2e / %.2e, enriched flux %.2e / %.2e, psi* %.2e / %.2e\n",
		             ut.rms, ut.worst, uStar.rms, uStar.worst, ps.rms, ps.worst );
		for ( auto const &[ name, unused ] : series.front().psiStarMethod )
		{
			std::string const key = name;
			Sigma const sp = column( [ & ]( State const &s ) { return s.psiStarMethod.at( key ); } );
			std::printf( "      psi*, %-66s %.2e / %.2e\n", name.c_str(), sp.rms, sp.worst );
		}
		for ( auto const &[ name, unused ] : series.front().totalFluxMethod )
		{
			std::string const key = name;
			Sigma const sp = column( [ & ]( State const &s ) { return s.totalFluxMethod.at( key ); } );
			std::printf( "      total flux u_t, %-42s %.2e / %.2e\n", name.c_str(), sp.rms, sp.worst );
		}
		for ( auto const &[ name, unused ] : series.front().reconstructionLoad )
		{
			std::string const key = name;
			Sigma const sp = column( [ & ]( State const &s ) { return s.reconstructionLoad.at( key ); } );
			std::printf( "      the local problem's load, %-32s %.2e / %.2e\n", name.c_str(), sp.rms, sp.worst );
		}

		State const &mid = series[ states/2 ];
		printDiagnostics( "the recovery's local system, middle state", mid.diagnostics );
		printDiagnostics( "Reconstruct()'s local system, middle state", mid.reconstructionDiagnostics );

	}
}

BOOST_AUTO_TEST_CASE( n12k2 ) { measure( 12, 2 ); }
BOOST_AUTO_TEST_CASE( n16k2 ) { measure( 16, 2 ); }
BOOST_AUTO_TEST_CASE( n32k2 ) { measure( 32, 2 ); }
BOOST_AUTO_TEST_CASE( n16k3 ) { measure( 16, 3 ); }
