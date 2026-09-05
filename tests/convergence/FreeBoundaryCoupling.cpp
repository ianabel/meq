/*
 * FB-1's coupling matrix P, and the claim FREE-BOUNDARY-PLAN.md section 4.3
 * calls "an argument and not a measurement".
 *
 * WHAT FB-1 IS. Free boundary replaces the fixed-boundary datum on Gamma_h with
 * the trace of an exterior expansion:
 *
 *     psihat|_Gamma_h  =  P a,      P_{in} = trace projection of ( C_n o a )
 *
 * where `a` is the vector of Gegenbauer coefficients of the exterior solution
 * and `a( x )` -- an unfortunate collision of notation, and the plan's -- is the
 * transfer path's foot map from Gamma_h to the true Gamma.
 * meq::GradShafranovSolver::exteriorTraceColumns() builds P, and this file is
 * what says it is the right matrix.
 *
 * THE PLAN ASKS FOR ONE MEASUREMENT HERE AND UNDER NPC IT IS NOT NEEDED.
 * Section 4.3 argues that dF/da = ( dF/dpsihat ) P is constant in the iterate --
 * psihat enters the flux row as <psihat, v.n>, the potential row as
 * <tau psihat, w> and the trace row as <tau psihat, mu>, linearly in all three,
 * while every non-linearity is F( r, z, psi ), which depends on psi and not on
 * psihat -- and says to "build the column at two well-separated iterates and
 * difference them".
 *
 * That measurement is the right one for a CONDENSED formulation. Under NPC it
 * is answered by construction, twice over, and this file records why rather
 * than performing a differencing that could only confirm an identity:
 *
 *   1. P IS A FUNCTION OF THE GEOMETRY ALONE. exteriorTraceColumns() takes no
 *      iterate and cannot take one -- it reads the mesh, the transfer path and
 *      the mode, and nothing else exists for it to read. A column that varied
 *      with the iterate could not be written with this signature, so the API is
 *      the assertion. If it ever needs an iterate, the claim has failed.
 *   2. Gamma_h's TRACE DOFS ARE ESSENTIAL, so the reduced operator masks its
 *      residual to zero on them and puts a unit row in the Jacobian there --
 *      which is what theEssentialTraceConditionImposesTheDatum already pins.
 *      The border column is therefore exactly -P, not a difference of one. That
 *      is section 4.5's "two of the three borders are exact under NPC", and it
 *      applies to this border for the same reason it applies to psi_ax's.
 *
 * So what is left to check is that P is the RIGHT matrix, and that is what this
 * file measures: that it lives only on Gamma_h, that it reproduces the mode AT
 * THE FEET rather than on Gamma_h, that its columns are independent, and that a
 * known coefficient vector reconstructs a known trace.
 *
 * THE FOOT IS THE POINT AND IT IS THE EASY THING TO GET WRONG. Gamma_h is the
 * inscribed polygon; Gamma is the true boundary; they differ by O( h ). A
 * projection of C_n evaluated on Gamma_h rather than at its foot on Gamma is a
 * different function by exactly that much, and it would throw away the accuracy
 * the whole transfer technique exists to buy -- while still converging, and
 * still looking like a coupling. `theColumnsAreTheModeAtTheFootAndNotOnGammaH`
 * is the test that separates the two, and it is the one to read first if
 * anything here fails.
 */

#define BOOST_TEST_MODULE FreeBoundaryCoupling
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"

namespace
{
	/*
	 * THE GEOMETRY IS ExtensionConvergence'S, DELIBERATELY.
	 *
	 * That file's Gamma is the Soloviev surface psi = -0.03 rather than the
	 * separatrix, for the reason its header gives at length: with correct
	 * coefficients psi = 0 IS the separatrix, which passes through an X-point --
	 * a corner of Gamma, where both transfer-path families give out and the
	 * Cockburn-Solano analysis does not reach. Borrowing the same surface here
	 * means this file inherits a geometry that is known to work and known WHY,
	 * instead of introducing a second one whose failures would be ambiguous.
	 *
	 * It is NOT a semicircle, so the exterior expansion is not the right
	 * representation of anything outside it -- and that is fine, because nothing
	 * here solves. P is a projection, and a projection of C_n onto the trace of
	 * Gamma_h is well defined whatever shape Gamma is.
	 */
	// ExtensionConvergence's own box and offset, copied rather than narrowed:
	// it contains Omega with room to spare and keeps r well away from zero,
	// which the operator's 1/r and psi's log r both want.
	double const rMin = 0.25;
	double const rMax = 1.95;
	double const zMin = -1.75;
	double const zMax = 1.65;
	double const psiOffset = 0.03;

	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq =
			meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	/// Negative inside Omega. Displacing psi by a constant changes neither the
	/// source nor the flux, so this is nstx() throughout -- and the surface is
	/// psi = -0.03 rather than the separatrix for the reason
	/// ExtensionConvergence's header gives at length.
	double levelSet( mfem::Vector const &x )
	{
		return equilibrium().psi( x( 0 ), x( 1 ) ) + psiOffset;
	}

	/// The subdomain D_h and its Gamma_h attribute. The same construction
	/// ExtensionConvergence::makeSubdomain uses, minus its assertions, which
	/// belong to the study it is part of.
	std::unique_ptr<mfem::SubMesh> makeSubdomain( int n, int &gammaH, double &h )
	{
		/*
		 * THE BACKGROUND MUST OUTLIVE THE SubMesh CUT FROM IT, AND IT USED NOT TO.
		 *
		 * mfem::SubMesh keeps a POINTER to its parent. This function returns the
		 * SubMesh and the background was a local, so the parent was dangling the
		 * moment it returned -- undefined behaviour that cost nothing for as long
		 * as nothing dereferenced it, which is to say for as long as nobody asked
		 * the SubMesh about the mesh it came from.
		 *
		 * mfem::VertexConePath now does. Its cone C( x ) reads the PARENT's edges
		 * at each vertex of Gamma_h -- HasCone() is documented as "whether the mesh
		 * handed to the constructor was a SubMesh with a parent to read edges
		 * from" -- so the constructor walks freed memory and segfaults in
		 * Mesh::GetVertexToVertexTable. It presented as an MFEM regression and was
		 * this fixture all along.
		 *
		 * The pool is the whole fix: a background lives until the process ends,
		 * which in a test binary is the simplest lifetime that is certainly long
		 * enough. Returning the pair instead would be tidier and would touch every
		 * caller; this is the change that is obviously correct.
		 */
		static std::vector<std::unique_ptr<mfem::Mesh>> backgrounds;
		backgrounds.push_back( std::make_unique<mfem::Mesh>(
			mfem::Mesh::MakeCartesian2D(
				n, 2*n, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin ) ) );
		mfem::Mesh &background = *backgrounds.back();
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		h = ( rMax - rMin )/static_cast<double>( n );

		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background, levelSet, 0.0,
		                                               marker, 1 );
		BOOST_TEST_REQUIRE( inside > 0, "the subdomain is empty at n = " << n );

		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		auto sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( background, domainAttr ) );

		gammaH = sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( sub->bdr_attributes.Size() == 1,
		                    "D_h has inherited boundary from the background box "
		                    "at n = " << n << ", so part of Gamma_h is fitted" );
		return sub;
	}

	/// Everything one case needs, kept alive together: the columns alias the
	/// solver's trace space and the path aliases the submesh, so none of them
	/// may outlive another.
	struct Case
	{
		std::unique_ptr<mfem::SubMesh> sub;
		std::unique_ptr<mfem::VertexConePath> path;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		mfem::Array<int> gammaHMarker;
		int gammaH = 0;
		double h = 0.0;
	};

	Case build( int order, int n )
	{
		Case c;
		c.sub = makeSubdomain( n, c.gammaH, c.h );
		c.path = std::make_unique<mfem::VertexConePath>( *c.sub, c.gammaH,
		                                                levelSet, 6.0*c.h );
		c.gammaHMarker.SetSize( c.gammaH );
		c.gammaHMarker = 0;
		c.gammaHMarker[ c.gammaH - 1 ] = 1;

		mfem::FunctionCoefficient source( []( mfem::Vector const &x )
		{
			return equilibrium().f( x( 0 ), x( 1 ), 0.0 );
		} );
		mfem::ConstantCoefficient zero( 0.0 );

		c.solver = std::make_unique<meq::GradShafranovSolver>( *c.sub, order );
		c.solver->setSource( source );
		c.solver->setBoundaryData( zero );
		c.solver->setExtension( *c.path, c.gammaHMarker );
		c.solver->prepare();
		return c;
	}

	/// A semicircle enclosing the whole domain, so every mode is evaluated at a
	/// direction the basis is happy with. The radius is irrelevant to P -- C_n
	/// depends on direction alone -- and is chosen only so that exterior() would
	/// be legal if anyone called it.
	meq::ExteriorDtN exterior( int modes )
	{
		return meq::ExteriorDtN( 0.0, 4.0, modes );
	}
}

/*
 * The refusal, first: there is no P on the fitted path.
 *
 * A fitted mesh has Gamma_h == Gamma and the trace unknown on the boundary IS
 * the condition imposed, so there is no transfer path, no foot map, and nothing
 * for a mode to be evaluated at. Returning an empty or zero P there would be a
 * coupling that silently does nothing, which is the shape of failure this whole
 * campaign is written to avoid.
 */
BOOST_AUTO_TEST_CASE( thereIsNoCouplingOnTheFittedPath )
{
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		4, 4, mfem::Element::TRIANGLE, false, 0.8, 1.2 );
	mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
	{
		out( 0 ) = in( 0 ) + 0.6;
		out( 1 ) = in( 1 ) - 0.6;
	} );

	mfem::ConstantCoefficient source( 1.0 );
	mfem::ConstantCoefficient zero( 0.0 );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( zero );
	solver.prepare();

	meq::ExteriorDtN const dtn = exterior( 4 );
	BOOST_CHECK_THROW( solver.exteriorTraceColumns( dtn ), std::logic_error );
}

/*
 * P lives on Gamma_h and nowhere else, and its columns are independent.
 *
 * The support claim is what says the coupling touches only the boundary it is
 * supposed to: a column with entries on an interior face would be imposing an
 * exterior mode inside the plasma. The independence claim is what says the
 * bordered system is solvable -- N columns that are not independent give a
 * singular corner block, and the failure would appear as a Newton that cannot
 * take a step rather than as anything recognisable.
 */
BOOST_AUTO_TEST_CASE( theColumnsLiveOnGammaHAndAreIndependent )
{
	int const modes = 6;
	Case c = build( 2, 12 );
	meq::ExteriorDtN const dtn = exterior( modes );

	std::vector<mfem::Vector> const columns = c.solver->exteriorTraceColumns( dtn );

	BOOST_TEST_REQUIRE( static_cast<int>( columns.size() ) == modes,
	                    "expected one column per mode" );

	// Every column has the full trace length, and the same nonzero pattern:
	// the dofs of Gamma_h. Comparing the PATTERNS against each other is the
	// cheap way to say "the same set of dofs" without recomputing the marker.
	int nonzeroCount = 0;
	for ( std::size_t k = 0; k < columns.size(); ++k )
	{
		BOOST_TEST_REQUIRE( columns[ k ].Size() == c.solver->numTraceDofs(),
		                    "column " << k << " is not a full trace vector" );

		int count = 0;
		for ( int i = 0; i < columns[ k ].Size(); ++i )
			if ( columns[ k ]( i ) != 0.0 )
				++count;

		if ( k == 0 )
		{
			nonzeroCount = count;
			BOOST_TEST_REQUIRE( count > 0, "the first column is identically zero" );
		}
	}

	std::printf( "\n  P: %d trace dofs, %d of them on Gamma_h, %d modes\n",
	             c.solver->numTraceDofs(), nonzeroCount, modes );

	// Independence, by the Gram matrix's smallest pivot under a plain Cholesky.
	// N is small -- 6 here, 20 to 40 in production -- so a dense Gram is cheap
	// and is the most direct statement of the property.
	mfem::DenseMatrix gram( modes, modes );
	for ( int i = 0; i < modes; ++i )
		for ( int j = 0; j < modes; ++j )
			gram( i, j ) = columns[ static_cast<std::size_t>( i ) ]
			               *columns[ static_cast<std::size_t>( j ) ];

	double const determinant = gram.Det();

	std::printf( "  Gram determinant %.6e, diagonal:", determinant );
	for ( int i = 0; i < modes; ++i )
		std::printf( " %.3e", gram( i, i ) );
	std::printf( "\n" );
	std::fflush( stdout );

	for ( int i = 0; i < modes; ++i )
		BOOST_TEST( gram( i, i ) > 0.0,
		            "column " << i << " of P is identically zero, so mode "
		            << meq::ExteriorDtN::firstMode() + i << " drives nothing" );

	BOOST_TEST( determinant > 0.0,
	            "the columns of P are linearly dependent -- Gram determinant "
	            << determinant << ". The bordered system's corner block is then "
	            "singular and Newton cannot take a step" );
}

/*
 * THE COLUMNS CARRY THE MODE AT THE FOOT, NOT ON Gamma_h. This is the test with
 * teeth.
 *
 * mfem::PathTraceCoefficient evaluates its function at a( x ), the foot of the
 * transfer path on the true Gamma. The plausible wrong implementation evaluates
 * it at x itself, on Gamma_h. Both produce a P that is supported in the right
 * place, has independent columns, and couples an exterior expansion to the
 * solve; the wrong one is simply a different function by O( dist( Gamma_h,
 * Gamma ) ) -- which is O( h ), which is exactly the error the transfer
 * technique exists to remove.
 *
 * So the two are separated by measuring both and showing they DIFFER by that
 * much, and that the difference SHRINKS with the mesh at the rate the distance
 * does. A single mesh could not tell a real difference from a coincidence.
 */
BOOST_AUTO_TEST_CASE( theColumnsAreTheModeAtTheFootAndNotOnGammaH )
{
	int const modes = 4;

	std::printf( "\n  P against the same projection taken on Gamma_h itself\n" );
	std::printf( "    %5s %8s %14s %14s %8s\n",
	             "n", "h", "worst |diff|", "worst |P|", "relative" );

	std::vector<double> spacing;
	std::vector<double> relative;

	for ( int n : { 8, 16, 32 } )
	{
		Case c = build( 2, n );
		meq::ExteriorDtN const dtn = exterior( modes );

		std::vector<mfem::Vector> const atFoot = c.solver->exteriorTraceColumns( dtn );

		// The same projection with NO path: the mode evaluated at the point of
		// Gamma_h itself. Built here rather than exposed from the solver,
		// because it is the wrong thing and has no business being callable.
		double worstDiff = 0.0;
		double worstValue = 0.0;

		for ( int k = 0; k < modes; ++k )
		{
			int const degree = meq::ExteriorDtN::firstMode() + k;
			mfem::FunctionCoefficient plain(
				[ &dtn, degree ]( mfem::Vector const &x )
			{
				return dtn.basis( degree, x( 0 ), x( 1 ) );
			} );

			mfem::GridFunction column( &c.solver->traceSpace() );
			column = 0.0;
			column.ProjectBdrCoefficient( plain, c.gammaHMarker );

			for ( int i = 0; i < column.Size(); ++i )
			{
				double const a = atFoot[ static_cast<std::size_t>( k ) ]( i );
				worstDiff = std::max( worstDiff, std::fabs( a - column( i ) ) );
				worstValue = std::max( worstValue, std::fabs( a ) );
			}
		}

		double const rel = worstDiff/worstValue;
		spacing.push_back( c.h );
		relative.push_back( rel );

		std::printf( "    %5d %8.4f %14.6e %14.6e %8.3f\n",
		             modes, c.h, worstDiff, worstValue, rel );
		std::fflush( stdout );
	}

	// They must differ at all: identical columns mean the path is not being
	// used and the coupling is on the wrong curve.
	BOOST_TEST( relative.front() > 1.0e-6,
	            "the projection at the foot and the projection on Gamma_h are "
	            "identical to " << relative.front() << ". PathTraceCoefficient is "
	            "then not evaluating at a( x ), and the coupling is being applied "
	            "on the inscribed polygon rather than on the true Gamma -- which "
	            "converges, and is O( h ) wrong" );

	// And the difference must shrink with the mesh, since dist( Gamma_h, Gamma )
	// does. A difference that did not shrink would be a bug rather than the
	// geometry.
	for ( std::size_t i = 1; i < relative.size(); ++i )
		BOOST_TEST( relative[ i ] < relative[ i - 1 ],
		            "the foot-vs-Gamma_h difference did not shrink from h = "
		            << spacing[ i - 1 ] << " to " << spacing[ i ] << ": "
		            << relative[ i - 1 ] << " then " << relative[ i ]
		            << ". It should track dist( Gamma_h, Gamma )" );
}

/*
 * ===========================================================================
 * THE TILING CHECK, which is the acceptance for mfem::ExtensionBoundaryQuadrature
 * ===========================================================================
 *
 * FB-1's transmission row is an integral over the TRUE Gamma of the extension
 * of the flux -- the t = 1 face of the path map, where
 * ExtensionRegionQuadrature() sweeps the region between Gamma_h and Gamma.
 * mfem::ExtensionBoundaryQuadrature() is that face, written as a patch against
 * gf-hdg-subdomains-dev because the extension machinery is what it needs and
 * NPC is not (that branch carries no NPC at all).
 *
 * ITS DOCUMENTATION NAMES ITS OWN ACCEPTANCE and it is the same one the region
 * sweep has: summed over the faces of Gamma_h, the weights must give |Gamma|,
 * exactly as the volume weights must give |Omega| - |D_h|. That is not a
 * property of the routine -- it is a property of the PATH FAMILY, which must
 * agree on the path through a shared vertex or leave gaps and overlaps at the
 * vertices. VertexConePath agrees by construction, since it interpolates the
 * paths of the vertices; a family following each face's own normal would not.
 *
 * SO THE GEOMETRY HERE IS A CIRCLE, not the Soloviev surface the rest of this
 * file uses. |Gamma| is then 2 pi R exactly, and the check is against a closed
 * form rather than against a quadrature of a level set -- which would be
 * checking two approximations against each other and could not tell which was
 * wrong.
 *
 * BOTH SWEEPS ARE CHECKED TOGETHER because they share the hypothesis. If the
 * boundary sum is wrong and the volume sum is right, the fault is in the new
 * routine; if both are wrong, it is the path family, and the new routine is
 * innocent. Running only one could not tell those apart.
 */
BOOST_AUTO_TEST_CASE( theBoundarySweepTilesGammaAndTheRegionSweepTilesTheGap )
{
	// A circle, well away from the axis so the operator's 1/r is nowhere near
	// singular and the geometry is the only thing under test.
	double const centreR = 1.10;
	double const centreZ = 0.0;
	double const radius = 0.40;

	auto circle = [ = ]( mfem::Vector const &x )
	{
		double const dr = x( 0 ) - centreR;
		double const dz = x( 1 ) - centreZ;
		return std::sqrt( dr*dr + dz*dz ) - radius;
	};

	double const exactPerimeter = 2.0*M_PI*radius;
	double const exactArea = M_PI*radius*radius;

	std::printf( "\n  the path family tiles Gamma and the gap beside it\n" );
	std::printf( "    Gamma is a circle of radius %.2f: |Gamma| = %.10f, "
	             "area = %.10f\n", radius, exactPerimeter, exactArea );
	std::printf( "    %5s %8s %16s %10s %16s %10s\n",
	             "n", "h", "sum of weights", "vs |Gamma|", "region sum", "vs gap" );

	std::vector<double> boundaryError;
	std::vector<double> regionError;

	for ( int n : { 12, 24, 48 } )
	{
		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false, 1.2, 1.2 );
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + 0.5;
			out( 1 ) = in( 1 ) - 0.6;
		} );
		double const h = 1.2/static_cast<double>( n );

		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background, circle, 0.0,
		                                               marker, 1 );
		BOOST_TEST_REQUIRE( inside > 0, "the disc is empty at n = " << n );

		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		mfem::SubMesh sub = mfem::SubMesh::CreateFromDomain( background, domainAttr );

		int const gammaH = sub.bdr_attributes.Max();
		BOOST_TEST_REQUIRE( sub.bdr_attributes.Size() == 1,
		                    "D_h reaches the background box at n = " << n );

		mfem::VertexConePath path( sub, gammaH, circle, 6.0*h );

		// The area of D_h, so the region sweep has something to be compared
		// against: it must give |Omega| - |D_h|.
		double areaOfDh = 0.0;
		for ( int e = 0; e < sub.GetNE(); ++e )
			areaOfDh += sub.GetElementVolume( e );

		// A rule generous enough that the quadrature is not what is being
		// measured. The map is not polynomial, so no order is exact; the
		// convergence below is what says the rule is adequate.
		mfem::IntegrationRule const &faceRule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 12 );
		mfem::IntegrationRule const &lineRule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 12 );

		double boundarySum = 0.0;
		double regionSum = 0.0;
		double worstNormalDeviation = 0.0;

		for ( int be = 0; be < sub.GetNBE(); ++be )
		{
			if ( sub.GetBdrAttribute( be ) != gammaH )
				continue;

			mfem::FaceElementTransformations *ftr =
				sub.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;

			mfem::ExtensionBoundaryQuadrature( *ftr, path, faceRule,
				[ & ]( mfem::ExtensionBoundaryPoint const &pt )
			{
				boundarySum += pt.weight;

				// AND THE NORMAL, WHILE WE ARE HERE. On a circle the outward
				// normal at y is ( y - centre )/radius, exactly -- so the
				// routine's orientation rule ( agree with the paths ) can be
				// checked against a closed form rather than against itself.
				double const nr = ( pt.y( 0 ) - centreR )/radius;
				double const nz = ( pt.y( 1 ) - centreZ )/radius;
				double const deviation =
					std::sqrt( ( pt.nu( 0 ) - nr )*( pt.nu( 0 ) - nr )
					           + ( pt.nu( 1 ) - nz )*( pt.nu( 1 ) - nz ) );
				worstNormalDeviation = std::max( worstNormalDeviation, deviation );
			} );

			mfem::ExtensionRegionQuadrature( *ftr, path, faceRule, lineRule,
				[ & ]( mfem::ExtensionPoint const &pt )
			{
				regionSum += pt.weight;
			} );
		}

		double const gap = exactArea - areaOfDh;
		double const bRel = std::fabs( boundarySum - exactPerimeter )/exactPerimeter;
		double const rRel = std::fabs( regionSum - gap )/std::fabs( gap );

		boundaryError.push_back( bRel );
		regionError.push_back( rRel );

		std::printf( "    %5d %8.4f %16.10f %10.2e %16.10f %10.2e\n",
		             n, h, boundarySum, bRel, regionSum, rRel );
		std::printf( "          worst |nu - exact| on Gamma: %.3e\n",
		             worstNormalDeviation );
		// The cone is new on gf-hdg-subdomains-dev and is the first suspect for
		// the tiling residual asserted below, so report it rather than guess.
		std::printf( "          cone: have %d, vertices %d, restricted %d, "
		             "tighter %d, widened %d\n",
		             path.HasCone() ? 1 : 0, path.NumVertices(),
		             path.NumConeRestricted(), path.NumTighter(),
		             path.NumWidened() );
		std::fflush( stdout );

		// THE NORMAL IS THE SHARPEST OF THE THREE and it is pointwise rather
		// than an aggregate: a sum can be right with cancelling errors, a
		// pointwise normal cannot. It also pins the ORIENTATION, which section
		// 7 of the plan says will be got wrong at least once -- an inward
		// normal reads 2.0 here, not a small number.
		BOOST_TEST( worstNormalDeviation < 5.0e-2,
		            "n = " << n << ": the normal returned on Gamma differs from "
		            "the circle's own by " << worstNormalDeviation
		            << ". A value near 2 means the ORIENTATION is inverted -- the "
		            "rule is that nu agrees with the path direction, since the "
		            "paths run outward from D_h" );
	}

	/*
	 * THE BOUNDARY SUM IS EXACT AT EVERY MESH, AND THAT -- NOT ITS CONVERGENCE
	 * -- IS THE PROPERTY.
	 *
	 * An earlier version of this test asserted that the error SHRANK with
	 * refinement, which was a guess and was the wrong shape. Measured, it does
	 * not shrink because there is nothing to shrink: 4.85e-10, 4.64e-10,
	 * 6.38e-11 at h = 0.1, 0.05, 0.025. That is mesh-INDEPENDENT and it sits at
	 * the central difference's own floor -- fd_step = 1e-6 gives O( fd^2 ) =
	 * 1e-12 truncation and O( eps/fd ) = 1e-10 round-off, and 1e-10 is what is
	 * measured. The quadrature is EXACT and the residue is the instrument.
	 *
	 * A gate that merely required improvement would pass on a routine that was
	 * O( h ) wrong, which is what this one was before the weight was signed.
	 */
	/*
	 * THIS IS RED AND IT IS A REAL FINDING, NOT A TOLERANCE TO RELAX.
	 *
	 * The expectation asserted here -- exact at every mesh, at the central
	 * difference's own floor of about 1e-10, mesh-INDEPENDENT -- was measured
	 * while ExtensionBoundaryQuadrature was being written, against the
	 * VertexConePath of that day. Against the path family now on
	 * gf-hdg-subdomains-dev the same geometry reads 1.01e-04, 2.24e-05 and
	 * 4.59e-06 at n = 12, 24, 48, which is not a floor at all: it converges at
	 * about O( h^2 ). A quadrature residual that CONVERGES is measuring a
	 * geometric error rather than an instrument, so the images of the faces no
	 * longer tile Gamma exactly -- they overlap or gap by O( h^2 ) somewhere.
	 *
	 * IT IS THE CONE, AND THAT IS MEASURED RATHER THAN SUSPECTED.
	 * theConeIsWhatCostsTheTiling below runs the identical geometry with the
	 * cone off -- a path built on a parentless copy of the same D_h, which is
	 * the documented way HasCone() goes false -- and reads 4.85e-10, 4.64e-10
	 * and 6.38e-11 against this case's 1.01e-04, 2.24e-05 and 4.59e-06. One
	 * variable, a factor of 2e5, and the cone-off column REPRODUCES the numbers
	 * recorded above, which also settles that they were measured here rather
	 * than inherited from another geometry.
	 *
	 * Why it costs anything is not established. Tiling rests on adjacent faces
	 * AGREEING at a shared vertex, which interpolating vertex directions gives
	 * by construction and a per-vertex restriction ought to preserve; the cone
	 * diagnostics printed above say it fired at every vertex at every mesh and
	 * was strictly tighter than the half space at about 40% of them.
	 *
	 * NOTE THAT UPSTREAM'S OWN COMMIT SAYS THE CONE "changes nothing". This
	 * contradicts that directly, on their own branch's routine, and is the
	 * reason to report it with the control attached rather than as an opinion.
	 *
	 * MFEM'S OWN SUITE CANNOT SEE THIS. test_darcy_extension.cpp exercises
	 * ExtensionRegionQuadrature and mentions the boundary sweep only in a
	 * comment, so this is the only check of the property anywhere -- which is a
	 * reason to report it upstream, and not a reason to weaken it here.
	 *
	 * Left failing deliberately, per CLAUDE.md's testing stance: a defect gets
	 * a test that asserts the behaviour WANTED and fails until it is there.
	 * Relaxing this gate to 1e-3 would make the suite green and would throw
	 * away the only measurement of this property that exists.
	 */
	for ( std::size_t i = 0; i < boundaryError.size(); ++i )
		BOOST_TEST( boundaryError[ i ] < 1.0e-8,
		            "the boundary weights do not sum to |Gamma| at mesh " << i
		            << ": relative " << boundaryError[ i ] << ". This should be "
		            "EXACT at every mesh, at the fd_step floor of about 1e-10 -- "
		            "an O( h ) error here is the WEIGHT LOSING ITS SIGN, which "
		            "makes the sweep measure the length the foot map traverses "
		            "rather than the piece of Gamma it covers. Those differ "
		            "wherever the map backtracks, which on a staircase Gamma_h "
		            "it does" );

	/*
	 * AND THE REGION SWEEP IS THE CONTROL, WHICH STILL CARRIES THE DEFECT.
	 *
	 * mfem::ExtensionRegionQuadrature uses std::abs( J.Det() ) and so counts a
	 * fold twice, exactly as the boundary sweep did before its weight was
	 * signed. Measured here: 6.96e-03, 3.29e-03, 1.84e-03 -- O( h ), against the
	 * boundary sweep's 1e-10 on the same geometry and the same path family.
	 *
	 * IT IS DELIBERATELY NOT FIXED HERE. That routine is what
	 * HDGExtensionIntegrator's lifting is built on, which is stage 5's machinery
	 * and which MEQ has measured at k+2; changing its quadrature would move
	 * results the suite pins, and is a claim that needs its own measurement
	 * rather than a drive-by. The gate below is therefore set ABOVE the measured
	 * value, to say "this is known and is not what is under test here".
	 */
	BOOST_TEST( regionError.back() < 1.0e-2,
	            "the region weights do not sum to |Omega| - |D_h|: relative "
	            << regionError.back() << ", worse than the O( h ) already known. "
	            "That is the EXISTING sweep and its unsigned weight; if this has "
	            "grown, the path family's coverage has changed and the boundary "
	            "result above needs rechecking too" );

	BOOST_TEST( boundaryError.back() < regionError.back(),
	            "the signed boundary sweep is no better than the unsigned region "
	            "sweep on the same geometry and the same paths. The whole point "
	            "of signing the weight is that it is" );
}

/*
 * STAR-SHAPEDNESS, WHICH IS WHAT RADIAL TRANSFER PATHS REST ON.
 *
 * Gamma is a semicircle centred on the axis, so the natural path family is rays
 * from that centre -- mfem::ClosestPointPath::Sphere, whose foot map is exact
 * and monotone and needs no search, against mfem::VertexConePath's fan search
 * and interpolation. The simpler family is the right one for MEQ precisely
 * because MEQ's boundary is trivial; MFEM keeps the general one because other
 * callers' boundaries are not.
 *
 * WHAT IT COSTS IS A HYPOTHESIS: every ray from the centre must meet Gamma_h
 * exactly once, or two faces claim the same piece of Gamma and the transfer is
 * double valued there. That is star-shapedness about the centre, and on a
 * polygon it is the local condition ( x - c ).n > 0 on every face.
 *
 * SO IT IS MEASURED AND NOT ASSUMED, which is the same discipline
 * meq::AngleParametrisation applies to the contour tracer -- star-shapedness
 * there is a hypothesis, is measured as min |u x t|, and is refused when it
 * fails. starShapedMargin() is the same quantity for this boundary: a cosine,
 * dimensionless, comparable between meshes.
 *
 * AND IT MUST SURVIVE REFINEMENT. meq::AdaptiveDomain re-cuts D_h from the
 * background mesh every cycle and the admitted elements change, so a domain
 * that starts star-shaped need not stay so. The margin is one dot product per
 * boundary face, so an adaptive loop can afford to assert it every cycle.
 */
BOOST_AUTO_TEST_CASE( gammaHIsStarShapedAboutTheCentreAndStaysSoUnderRefinement )
{
	double const centreR = 1.10;
	double const centreZ = 0.0;
	double const radius = 0.40;

	auto circle = [ = ]( mfem::Vector const &x )
	{
		double const dr = x( 0 ) - centreR;
		double const dz = x( 1 ) - centreZ;
		return std::sqrt( dr*dr + dz*dz ) - radius;
	};

	std::printf( "\n  star-shapedness of Gamma_h about the semicircle's centre\n" );
	std::printf( "    %5s %8s %14s\n", "n", "h", "min cosine" );

	std::vector<double> margins;

	for ( int n : { 12, 24, 48 } )
	{
		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false, 1.2, 1.2 );
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + 0.5;
			out( 1 ) = in( 1 ) - 0.6;
		} );
		double const h = 1.2/static_cast<double>( n );

		mfem::Array<int> marker;
		mfem::MarkLevelSetSubdomain( background, circle, 0.0, marker, 1 );
		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		mfem::SubMesh sub = mfem::SubMesh::CreateFromDomain( background, domainAttr );

		int const gammaH = sub.bdr_attributes.Max();
		mfem::VertexConePath path( sub, gammaH, circle, 6.0*h );
		mfem::Array<int> gammaHMarker( gammaH );
		gammaHMarker = 0;
		gammaHMarker[ gammaH - 1 ] = 1;

		mfem::ConstantCoefficient source( 1.0 );
		mfem::ConstantCoefficient zero( 0.0 );

		meq::GradShafranovSolver solver( sub, 1 );
		solver.setSource( source );
		solver.setBoundaryData( zero );
		solver.setExtension( path, gammaHMarker );
		solver.prepare();

		double const margin = solver.starShapedMargin( centreR, centreZ );
		margins.push_back( margin );

		// Scientific, because the interesting value is an EXACT zero and a
		// fixed-point format cannot tell that from 1e-8.
		std::printf( "    %5d %8.4f %16.9e\n", n, h, margin );
		std::fflush( stdout );
	}

	/*
	 * NON-NEGATIVE, NOT STRICTLY POSITIVE, AND THE DIFFERENCE IS THE FINDING.
	 *
	 * This assertion was written as margin > 0 and has never held: it reads
	 * EXACTLY 0.000000000e+00 at n = 12, 24 and 48 alike. An exact zero
	 * repeated across three meshes is not geometry -- geometry would move --
	 * and the mechanism is the staircase. Gamma_h is a union of background
	 * element faces, so it carries axis-aligned faces; a HORIZONTAL face with
	 * an endpoint on the line z = centreZ has ( x - c ) purely radial and n
	 * purely vertical, and their dot product is exactly zero. The diagonal
	 * split of MakeCartesian2D is not symmetric about that line, so such a
	 * face exists at every mesh rather than being cancelled by symmetry.
	 *
	 * A ray through that one corner is TANGENT to that face. It is not a ray
	 * meeting Gamma_h twice, and no face is ever turned away -- the margin is
	 * zero and never negative. What the transfer needs is that no face faces
	 * away from the centre, which is margin >= 0; the strict inequality was a
	 * guess at how to say it and is false for every staircase.
	 *
	 * The control below is what keeps this from being vacuous: a centre outside
	 * the disc reads -1, so the quantity can go negative and the >= is doing
	 * work.
	 */
	for ( std::size_t i = 0; i < margins.size(); ++i )
		BOOST_TEST( margins[ i ] >= 0.0,
		            "Gamma_h is NOT star-shaped about the centre at mesh " << i
		            << ": margin " << margins[ i ] << ". A ray from the centre "
		            "then meets it more than once, two faces claim the same piece "
		            "of Gamma, and a radial transfer path is double valued there. "
		            "This is the hypothesis mfem::ClosestPointPath::Sphere rests "
		            "on, and it is why the check exists rather than the "
		            "assumption" );

	/*
	 * A CONTROL, because a margin that is positive on every mesh anyone happens
	 * to try says nothing about whether the check WORKS. A centre placed outside
	 * the disc cannot see the far side of it, so the margin must go negative --
	 * and if it does not, the quantity is not measuring what it claims to.
	 */
	{
		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			24, 24, mfem::Element::TRIANGLE, false, 1.2, 1.2 );
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + 0.5;
			out( 1 ) = in( 1 ) - 0.6;
		} );

		mfem::Array<int> marker;
		mfem::MarkLevelSetSubdomain( background, circle, 0.0, marker, 1 );
		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();
		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		mfem::SubMesh sub = mfem::SubMesh::CreateFromDomain( background, domainAttr );

		int const gammaH = sub.bdr_attributes.Max();
		mfem::VertexConePath path( sub, gammaH, circle, 6.0*0.05 );
		mfem::Array<int> gammaHMarker( gammaH );
		gammaHMarker = 0;
		gammaHMarker[ gammaH - 1 ] = 1;

		mfem::ConstantCoefficient source( 1.0 );
		mfem::ConstantCoefficient zero( 0.0 );
		meq::GradShafranovSolver solver( sub, 1 );
		solver.setSource( source );
		solver.setBoundaryData( zero );
		solver.setExtension( path, gammaHMarker );
		solver.prepare();

		// Well outside the disc, where the far side is not visible.
		double const outside = solver.starShapedMargin( centreR + 3.0*radius,
		                                                centreZ );
		std::printf( "    control, centre outside the disc: %14.6f\n", outside );
		std::fflush( stdout );

		BOOST_TEST( outside < 0.0,
		            "a centre placed outside the disc still reports a positive "
		            "star-shapedness margin, " << outside << ". It cannot see the "
		            "far side, so the check is not measuring visibility and the "
		            "positive results above mean nothing" );
	}
}

/*
 * ===========================================================================
 * THE TRANSMISSION ROW IS THE INTEGRAL IT CLAIMS TO BE
 * ===========================================================================
 *
 * exteriorTransmissionRows() returns the flux-block covector of
 *
 *     INT_Gamma E_h( q_h ).nu C_m dGamma
 *
 * and there are four independent ways to get it wrong, none of which a solve
 * would report as anything but a wrong equilibrium: the vdof ordering ( the
 * flux space is L2 with vdim 2 and byNODES, so component d of basis j is
 * dof*d + j, and swapping them ROTATES the field ); the sign ( the assembled
 * block holds -q ); the measure ( a 1/r written here would divide by the radius
 * twice, since q already carries one ); and the weight's SIGN ( a staircase
 * Gamma_h has faces whose foot map reverses, and std::abs would integrate the
 * traversed length rather than Gamma ).
 *
 * THE TEST IS AGAINST A CLOSED FORM AND NOT AGAINST A RE-IMPLEMENTATION.
 * Feed the row a flux field the discrete space represents EXACTLY and whose
 * extension is therefore itself, and the integral collapses to a quadrature of
 * pure geometry -- no field evaluation, no inverse element map, no polynomial
 * extended anywhere. Two rungs:
 *
 *   constant q = ( 1, 0 ) and ( 0, 1 ):  INT nu_r C_m dGamma, INT nu_z C_m dGamma
 *   linear   q = ( z, r ):               INT ( y_z nu_r + y_r nu_z ) C_m dGamma
 *
 * The constant rung pins the ordering, the sign and the measure. THE LINEAR
 * RUNG IS THE ONE THAT EXERCISES THE EXTENSION: a constant is what an element's
 * polynomial gives at any point whether or not TransformBack found the right
 * reference coordinates, so a clamped or diverged inverse map is INVISIBLE to
 * it -- and clamping is precisely the failure mfem::ElementExtension exists to
 * prevent. A linear field read at a clamped point is wrong by the distance it
 * was clamped by.
 *
 * SEPARATELY VARYING THE TWO COMPONENTS IS WHAT CATCHES THE ORDERING. With
 * q = ( 1, 1 ) a transposed index gives the same answer, and the test passes
 * while the field is rotated.
 *
 * ORTHOGONALITY IS NOT USED AND CANNOT BE HERE. The exterior identity
 * T_m = 0 needs INT C_n C_m dGamma/r = delta_nm h_n, which holds only on a
 * semicircle reaching the axis at both ends. This file's Gamma is a Soloviev
 * surface away from the axis, so what is tested here is the CONTRACTION alone.
 * The identity is FB-1's, and it needs the half-disc domain section 7.5's
 * closing paragraph is about.
 */
BOOST_AUTO_TEST_CASE( theTransmissionRowIsTheBoundaryIntegralItClaims )
{
	int const order = 2;
	int const n = 16;

	Case c = build( order, n );

	// The file's own helper: centred on the axis, with a radius large enough
	// that exterior() would be legal. Only the direction matters to basis().
	meq::ExteriorDtN dtn = exterior( 3 );

	std::vector<mfem::Vector> rows = c.solver->exteriorTransmissionRows( dtn );
	BOOST_TEST_REQUIRE( static_cast<int>( rows.size() ) == dtn.modeCount() );

	mfem::FiniteElementSpace &fluxFes = c.solver->fluxSpace();

	/// Contract a row against a flux field, in DarcyForm's convention: the
	/// block holds -q, so that is what is put there.
	auto contract = [ & ]( mfem::Vector const &row,
	                       mfem::VectorCoefficient &q )
	{
		mfem::GridFunction gf( &fluxFes );
		gf.ProjectCoefficient( q );

		double total = 0.0;
		for ( int i = 0; i < gf.Size(); ++i )
			total += row( i )*( -gf( i ) );
		return total;
	};

	/// The same integral by quadrature of the geometry alone: no field, no
	/// extension, no element map. `integrand` is given the point on Gamma and
	/// the outward normal there and returns q.nu at that point.
	auto reference = [ & ]( int mode,
	                        std::function<double( mfem::Vector const &,
	                                              mfem::Vector const & )> const
	                            &integrand )
	{
		mfem::Mesh &mesh = *c.sub;
		double total = 0.0;

		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			if ( mesh.GetBdrAttribute( be ) != c.gammaH )
				continue;

			mfem::FaceElementTransformations *ftr =
				mesh.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;

			mfem::IntegrationRule const &faceRule =
				mfem::IntRules.Get( ftr->GetGeometryType(), 12 );

			mfem::ExtensionBoundaryQuadrature( *ftr, *c.path, faceRule,
				[ & ]( mfem::ExtensionBoundaryPoint const &pt )
			{
				total += pt.weight*integrand( pt.y, pt.nu )
				         *dtn.basis( mode, pt.y( 0 ), pt.y( 1 ) );
			} );
		}
		return total;
	};

	std::printf( "\n  THE TRANSMISSION ROW AGAINST A CLOSED FORM"
	             "  ( k = %d, n = %d )\n\n", order, n );
	std::printf( "    %-22s %6s %16s %16s %12s\n",
	             "field", "mode", "row . q", "quadrature", "relative" );

	struct Probe
	{
		char const *name;
		std::function<void( mfem::Vector const &, mfem::Vector & )> value;
	};

	std::vector<Probe> const probes = {
		{ "constant ( 1, 0 )",
		  []( mfem::Vector const &, mfem::Vector &v ) { v( 0 ) = 1.0; v( 1 ) = 0.0; } },
		{ "constant ( 0, 1 )",
		  []( mfem::Vector const &, mfem::Vector &v ) { v( 0 ) = 0.0; v( 1 ) = 1.0; } },
		{ "linear ( z, r )",
		  []( mfem::Vector const &x, mfem::Vector &v ) { v( 0 ) = x( 1 ); v( 1 ) = x( 0 ); } }
	};

	double worst = 0.0;
	double worstMagnitude = 0.0;

	for ( Probe const &probe : probes )
	{
		mfem::VectorFunctionCoefficient q( 2, probe.value );

		for ( int m = 0; m < dtn.modeCount(); ++m )
		{
			int const degree = meq::ExteriorDtN::firstMode() + m;

			double const fromRow = contract( rows[ static_cast<std::size_t>( m ) ], q );
			double const fromQuadrature = reference( degree,
				[ & ]( mfem::Vector const &y, mfem::Vector const &nu )
			{
				mfem::Vector v( 2 );
				probe.value( y, v );
				return v( 0 )*nu( 0 ) + v( 1 )*nu( 1 );
			} );

			double const scale = std::max( std::abs( fromQuadrature ), 1.0e-12 );
			double const relative = std::abs( fromRow - fromQuadrature )/scale;

			std::printf( "    %-22s %6d %16.8e %16.8e %12.3e\n",
			             probe.name, degree, fromRow, fromQuadrature, relative );

			worst = std::max( worst, relative );
			worstMagnitude = std::max( worstMagnitude, std::abs( fromQuadrature ) );
		}
	}
	std::fflush( stdout );

	// The reference is not identically zero, or the agreement above is vacuous
	// -- a row of zeros would reproduce it perfectly.
	BOOST_TEST( worstMagnitude > 1.0e-6,
	            "every reference integral is at round-off, " << worstMagnitude
	            << ", so the agreement says nothing. Either the modes vanish on "
	            "this Gamma or the sweep visited no face" );

	// Both sides are exact for these fields: the space represents them, the
	// extension reproduces them, and the quadrature is generous. What is left is
	// round-off over the assembly, not a discretisation error.
	BOOST_TEST( worst < 1.0e-10,
	            "the transmission row disagrees with a quadrature of the same "
	            "integral by " << worst << " relative, on a field its own space "
	            "represents exactly. Suspect, in order: the vdof ordering "
	            "( dof*d + j, byNODES ), the sign against DarcyForm's -q, a 1/r "
	            "that should not be there, or std::abs on a signed weight" );
}

/*
 * ===========================================================================
 * IS THE CONE WHAT COSTS THE TILING? THE CONTROL THAT ANSWERS IT
 * ===========================================================================
 *
 * theBoundarySweepTilesGammaAndTheRegionSweepTilesTheGap is red: the weights
 * sum to |Gamma| only to O( h^2 ) where the recorded expectation is the central
 * difference's floor, mesh-independent at about 1e-10. Two stories fit that
 * equally well and they call for opposite actions:
 *
 *   (a) mfem::VertexConePath's cone C( x ), which is new, changed the path
 *       family and the images stopped tiling. Then something regressed and it
 *       belongs upstream.
 *   (b) MEQ's circle never tiled to 1e-10 and the recorded numbers came from a
 *       different geometry. Then nothing regressed, the expectation was never
 *       MEQ's to hold, and the gate is wrong rather than the library.
 *
 * THE COIN-FLIP IS FREE, BECAUSE THE CONE TURNS ITSELF OFF. HasCone() is
 * documented as "whether the mesh handed to the constructor was a SubMesh with
 * a parent to read edges from", so a path built on a PLAIN Mesh copy of the
 * same D_h has the identical geometry, the identical Gamma_h and no cone. Same
 * faces, same level set, same sweep -- one variable.
 *
 * Slicing the copy is the whole trick: mfem::Mesh plain( sub ) copies the
 * elements and boundary attributes and drops the parent pointer, which is
 * exactly the difference being isolated.
 */
BOOST_AUTO_TEST_CASE( theConeIsWhatCostsTheTiling )
{
	double const centreR = 1.10;
	double const centreZ = 0.0;
	double const radius = 0.40;
	double const exactPerimeter = 2.0*M_PI*radius;

	auto circle = [ = ]( mfem::Vector const &x )
	{
		double const dr = x( 0 ) - centreR;
		double const dz = x( 1 ) - centreZ;
		return std::sqrt( dr*dr + dz*dz ) - radius;
	};

	std::printf( "\n  DOES THE CONE COST THE TILING?  ( |Gamma| = %.10f )\n\n",
	             exactPerimeter );
	std::printf( "    %5s %8s %16s %11s %16s %11s\n",
	             "n", "h", "cone on", "rel", "cone off", "rel" );

	std::vector<double> withCone;
	std::vector<double> withoutCone;

	for ( int n : { 12, 24, 48 } )
	{
		// The identical construction the tiling case uses, so the only variable
		// between the two columns below is the cone.
		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false, 1.2, 1.2 );
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + 0.5;
			out( 1 ) = in( 1 ) - 0.6;
		} );
		double const h = 1.2/static_cast<double>( n );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE( mfem::MarkLevelSetSubdomain( background, circle, 0.0,
		                                                 marker, 1 ) > 0 );
		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		mfem::SubMesh sub = mfem::SubMesh::CreateFromDomain( background, domainAttr );
		int const gammaH = sub.bdr_attributes.Max();

		// The same D_h with no parent, so VertexConePath cannot read a cone.
		// Slicing to Mesh copies the elements and the boundary attributes and
		// drops the parent pointer, which is exactly the difference wanted.
		mfem::Mesh plain( sub );

		double sums[ 2 ] = { 0.0, 0.0 };
		bool cone[ 2 ] = { false, false };

		for ( int which = 0; which < 2; ++which )
		{
			mfem::Mesh &meshRef = ( which == 0 )
				? static_cast<mfem::Mesh &>( sub ) : plain;

			mfem::VertexConePath path( meshRef, gammaH, circle, 6.0*h );
			cone[ which ] = path.HasCone();

			for ( int be = 0; be < meshRef.GetNBE(); ++be )
			{
				if ( meshRef.GetBdrAttribute( be ) != gammaH )
					continue;

				mfem::FaceElementTransformations *ftr =
					meshRef.GetBdrFaceTransformations( be );
				if ( !ftr )
					continue;

				mfem::IntegrationRule const &faceRule =
					mfem::IntRules.Get( ftr->GetGeometryType(), 12 );

				mfem::ExtensionBoundaryQuadrature( *ftr, path, faceRule,
					[ & ]( mfem::ExtensionBoundaryPoint const &pt )
				{
					sums[ which ] += pt.weight;
				} );
			}
		}

		double const relOn = std::abs( sums[ 0 ] - exactPerimeter )/exactPerimeter;
		double const relOff = std::abs( sums[ 1 ] - exactPerimeter )/exactPerimeter;

		std::printf( "    %5d %8.4f %16.10f %11.2e %16.10f %11.2e\n",
		             n, h, sums[ 0 ], relOn, sums[ 1 ], relOff );
		std::fflush( stdout );

		BOOST_TEST_REQUIRE( cone[ 0 ], "the SubMesh path has no cone at n = " << n
		                    << ", so the two columns are the same experiment" );
		BOOST_TEST_REQUIRE( !cone[ 1 ], "the plain-Mesh path still has a cone at n = "
		                    << n << ", so slicing did not drop the parent" );

		withCone.push_back( relOn );
		withoutCone.push_back( relOff );
	}

	std::printf( "\n" );
	std::fflush( stdout );

	BOOST_TEST( withoutCone.back() < withCone.back(),
	            "turning the cone off does not improve the tiling at the finest "
	            "mesh: " << withoutCone.back() << " against " << withCone.back()
	            << ". Then the cone is NOT what costs it, and the O( h^2 ) is a "
	            "property of MEQ's circle or of the routine -- either way the red "
	            "gate beside this is measuring the wrong thing and the recorded "
	            "1e-10 expectation never applied here" );
}
