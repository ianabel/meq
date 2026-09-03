#define BOOST_TEST_MODULE FluxSurfaceConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"
#include "ConvergenceHarness.hpp"

/*
 * The acceptance test for INVERSION-PLAN.md stages IN-0 and IN-1: the
 * predictor-corrector contour tracer, and the poloidal-angle parametrisation
 * and metric built on it. ON THE FITTED PATH, where there is no band between
 * Gamma_h and Gamma and every rate is therefore attributable to the tracer
 * alone.
 *
 * SECTION 2's THREE-WAY ERROR SPLIT, MADE INTO SEPARATE NUMBERS. That split is
 * the point of the stage and conflating its three parts is how this problem
 * gets done badly, so each is measured on its own and none of them is inferred
 * from another:
 *
 *   ( a ) FIELD, { psi_h = c } against { psi = c }. Measured as
 *         | psi( x_i ) - c | / | grad psi( x_i ) | at the traced points against
 *         the CLOSED FORM, which is the first-order distance to the exact level
 *         set. Converges at k+1 in h.
 *   ( b ) POINT LOCATION, how far an accepted point is from { psi_h = c }.
 *         Measured as | psi_h( x_i ) - c |, which sits at the corrector
 *         tolerance and -- the property the corrector exists for -- DOES NOT
 *         GROW WITH PATH LENGTH. One, five and ten circuits.
 *   ( c ) REPRESENTATION, how far the interpolant BETWEEN accepted points is
 *         from the curve. Measured at segment midpoints against the DISCRETE
 *         level set, cubic Hermite against straight chords, with the chords
 *         kept as the control in the way ExtensionConvergence.cpp keeps its
 *         pinned-zero column. Fourth order against second, and at fixed
 *         Delta_s independent of h -- which is what says ( b ) and ( c ) have
 *         been separated.
 *
 * IN-1 IS THE TRAP OF SECTION 3.2 MADE INTO A MEASUREMENT, and its acceptance
 * is three columns of one table: the periodic trapezoid rule with rho' taken
 * pointwise from q, the same rule with rho' obtained by differencing
 * neighbouring node radii, and a chord sum. The middle column is the
 * deliverable -- a spectrally accurate rule reduced to second order by its
 * Jacobian, with nothing in its own output saying so. If either control ever
 * converged as fast as the first the comparison would be empty and this test
 * would be worthless, and its failure message says exactly that.
 *
 * WHERE THE FIRST COLUMN FLOORS, AND WHY THAT IS NOT A DEFECT. psi_h is
 * discontinuous across faces, so rho( theta ) for the DISCRETE contour is
 * piecewise analytic with O( h^(k+1) ) jumps at every face crossing, and a
 * trapezoid rule cannot be spectral on a function with jumps. The first column
 * therefore converges very fast and then stops at the jump level rather than at
 * round-off. That is measured, and the control that says it is the FIELD and
 * not the RULE is theMetricIdentityIsSpectralOnTheClosedForm below, which runs
 * the identical rule on the analytic equilibrium -- where rho is analytic --
 * and reaches round-off.
 */

namespace
{

	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	double const twoPi = 6.283185307179586476925286766559;

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}

	/// The magnetic axis of the CLOSED FORM, to round-off. Newton on
	/// grad( psi ) = 0 with the Hessian by central differences; the Hessian's
	/// accuracy does not reach the answer, the fixed point being where the
	/// analytic gradient vanishes whatever steered it there. Lifted from
	/// CriticalPointConvergence.cpp, which is where the argument for it is.
	struct ExactAxis
	{
		double r;
		double z;
		double psi;
	};

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

			double a0 = 0.0;
			double a1 = 0.0;
			double b0 = 0.0;
			double b1 = 0.0;
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

		ExactAxis axis;
		axis.r = r;
		axis.z = z;
		axis.psi = eq.psi( r, z );
		return axis;
	}

	/// The levels to trace, chosen from the CLOSED FORM and therefore the same
	/// at every mesh and every degree.
	///
	/// That matters: a level read off the discrete solution would move with h,
	/// and the exact contour it is compared against would move with it, so the
	/// field error ( a ) would be measured against a moving target and its rate
	/// would mean nothing. The axis is an interior MINIMUM for every Solov'ev
	/// fixture -- F is single-signed negative, so psi is a subsolution and its
	/// maximum is on the boundary; CriticalPointConvergence.cpp measures that --
	/// so psi increases outward and the largest closed contour inside the box is
	/// bounded by the smallest boundary value.
	struct LevelSet
	{
		ExactAxis axis;
		double edge;
		std::vector<double> levels;
		std::vector<double> fractions;
	};

	LevelSet levelsFor( Equilibrium const &eq, Rectangle const &box,
	                    std::vector<double> const &fractions )
	{
		LevelSet set;
		set.axis = exactAxis( eq, 0.5*( box.rMin + box.rMax ),
		                      0.5*( box.zMin + box.zMax ) );
		set.fractions = fractions;

		set.edge = std::numeric_limits<double>::infinity();
		int const samples = 400;
		for ( int i = 0; i <= samples; ++i )
		{
			double const s = static_cast<double>( i )/samples;
			double const r = box.rMin + s*box.width();
			double const z = box.zMin + s*box.height();
			set.edge = std::min( set.edge, eq.psi( r, box.zMin ) );
			set.edge = std::min( set.edge, eq.psi( r, box.zMax ) );
			set.edge = std::min( set.edge, eq.psi( box.rMin, z ) );
			set.edge = std::min( set.edge, eq.psi( box.rMax, z ) );
		}

		for ( std::size_t i = 0; i < fractions.size(); ++i )
			set.levels.push_back( set.axis.psi
			                      + fractions[ i ]*( set.edge - set.axis.psi ) );
		return set;
	}

	/// One solve, kept alive, with the post-processing done so that both
	/// candidate potentials are available. The harness's measure() destroys its
	/// solver on the way out, which is right for an error norm and no use here:
	/// the tracer borrows the potential and the flux and needs them to outlive
	/// the measurement. Member order is load bearing and the class is
	/// non-copyable because the coefficients capture `this`. Lifted from
	/// CriticalPointConvergence.cpp.
	class SolvedEquilibrium
	{
		public:
			SolvedEquilibrium( Equilibrium const &eqIn, Rectangle const &boxIn,
			                   int orderIn, int n )
				: eq( eqIn ),
				  mesh( meq::tests::makeMesh( boxIn, n ) ),
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
				solver.setSource( sourceCoeff );
				solver.setBoundaryData( psiCoeff );
				solver.solve();
				solver.postProcess();
				h = boxIn.width()/static_cast<double>( n );
			}

			SolvedEquilibrium( SolvedEquilibrium const & ) = delete;
			SolvedEquilibrium &operator=( SolvedEquilibrium const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			double meshSize() const
			{
				return h;
			}

			meq::CriticalPoint axis() const
			{
				meq::CriticalPointFinder finder( solver );
				return finder.findAxis();
			}

		private:
			Equilibrium eq;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			meq::GradShafranovSolver solver;
			double h = 0.0;
	};

	/// The first-order distance from @a point to the EXACT level set { psi = c }.
	/// This is error ( a ), and note its 1/|grad psi| weighting: the geometric
	/// error of a contour is worst exactly where the gradient is small, which is
	/// section 2's point about the axis.
	double distanceToExactContour( Equilibrium const &eq, double level, double r,
	                               double z )
	{
		double gr = 0.0;
		double gz = 0.0;
		eq.gradPsi( r, z, gr, gz );
		double const magnitude = std::sqrt( gr*gr + gz*gz );
		return magnitude > 0.0 ? std::abs( eq.psi( r, z ) - level )/magnitude
		                       : std::numeric_limits<double>::infinity();
	}

	/// The worst of that over a traced contour.
	double fieldError( Equilibrium const &eq, meq::Contour const &contour )
	{
		double worst = 0.0;
		for ( std::size_t i = 0; i < contour.points.size(); ++i )
			worst = std::max( worst, distanceToExactContour(
				eq, contour.level, contour.points[ i ].r, contour.points[ i ].z ) );
		return worst;
	}

	/// Measurement ( c ): the distance from the interpolant, sampled between the
	/// accepted points, to the DISCRETE level set. | psi_h( x ) - c | / | grad
	/// psi_h( x ) | with grad_bar( psi_h ) = r q, which is the first-order
	/// normal distance and is evaluated INSIDE whichever element holds the
	/// sample -- nothing is ever read outside an element.
	///
	/// TWO COLUMNS, AND THE SECOND IS WHY. A segment whose two endpoints were
	/// corrected in DIFFERENT elements has its endpoints on two different arcs
	/// of { psi_h = c }, and no interpolant through them can be closer to either
	/// arc than the arcs are to each other. That distance is the DG jump over
	/// | grad psi | and it is error ( a ). So the error is reported over every
	/// segment AND over the CLEAN ones -- both endpoints and the sample in one
	/// element -- where the interpolation is measured against the single
	/// polynomial arc it is actually interpolating and there is no floor at all.
	struct RepresentationError
	{
		double hermite = 0.0;
		double chord = 0.0;
		double cleanHermite = 0.0;
		double cleanChord = 0.0;
		int clean = 0;
		int total = 0;
	};

	RepresentationError representationError( meq::ContourTracer const &tracer,
	                                         meq::Contour const &contour,
	                                         int samplesPerSegment = 1 )
	{
		RepresentationError worst;
		int hint = contour.points.empty() ? -1 : contour.points.front().element;

		for ( std::size_t i = 0; i + 1 < contour.points.size(); ++i )
		{
			// A CLEAN SEGMENT IS ONE WHOSE INTERPOLATION IS MEASURABLE, and
			// there are two ways for it not to be. Both endpoints and the sample
			// must be in one element, or the arc has a jump in the middle of it;
			// and both endpoints must have REACHED the corrector tolerance, or
			// the interpolant is being asked to pass through a point that is
			// not on the arc. The second is the rarer and was the one that had
			// to be measured to be believed: a single endpoint accepted at the
			// DG jump rather than at 1e-12 puts the whole segment's error at the
			// jump, and on a coarse mesh that is a hundred times the
			// interpolation error being looked for.
			bool const sameElement =
				contour.points[ i ].element == contour.points[ i + 1 ].element
				&& contour.points[ i ].residual <= contour.correctorTarget
				&& contour.points[ i + 1 ].residual <= contour.correctorTarget;

			for ( int j = 0; j < samplesPerSegment; ++j )
			{
				double const t = ( j + 1.0 )/( samplesPerSegment + 1.0 );

				double r = 0.0;
				double z = 0.0;
				double psi = 0.0;
				double qR = 0.0;
				double qZ = 0.0;

				++worst.total;
				bool clean = sameElement;

				contour.pointOnSegment( i, t, r, z );
				if ( tracer.sampleAt( r, z, psi, qR, qZ, hint ) )
				{
					clean = clean && hint == contour.points[ i ].element;
					double const magnitude = r*std::sqrt( qR*qR + qZ*qZ );
					if ( magnitude > 0.0 )
					{
						double const error = std::abs( psi - contour.level )/magnitude;
						worst.hermite = std::max( worst.hermite, error );
						if ( clean )
							worst.cleanHermite = std::max( worst.cleanHermite, error );
					}
				}

				contour.chordOnSegment( i, t, r, z );
				if ( tracer.sampleAt( r, z, psi, qR, qZ, hint ) )
				{
					clean = clean && hint == contour.points[ i ].element;
					double const magnitude = r*std::sqrt( qR*qR + qZ*qZ );
					if ( magnitude > 0.0 )
					{
						double const error = std::abs( psi - contour.level )/magnitude;
						worst.chord = std::max( worst.chord, error );
						if ( clean )
							worst.cleanChord = std::max( worst.cleanChord, error );
					}
				}

				if ( clean )
					++worst.clean;
			}
		}
		return worst;
	}

	/// The SAME traced points, with the unit tangents replaced by those of
	/// grad( psi_h ) rather than those of q.
	///
	/// IT EXISTS BECAUSE THE OBVIOUS MEASUREMENT IS SUBTLY THE WRONG ONE, AND
	/// THIS IS THE FINDING OF THE STAGE.
	///
	/// The tangent from q is ( -q_z, +q_r )/|q|. The tangent of the curve
	/// { psi_h = c } is ( -d_z psi_h, +d_r psi_h )/| grad psi_h |. THESE ARE NOT
	/// THE SAME VECTOR: q_h and grad( psi_h )/r are separate objects in a mixed
	/// method and differ by O( h^k ) -- q converges at k+1 and a differentiated
	/// L2 potential at k, which is the entire reason meq solves for q at all.
	///
	/// So an interpolant built on q's tangents is NOT interpolating
	/// { psi_h = c }: it passes through points on that curve with tangents
	/// belonging to a slightly different one. Measured against { psi_h = c } it
	/// is fourth order in Delta_s until the mismatch takes over and then second,
	/// with a constant proportional to the discrepancy. Against the TRUE contour
	/// -- which is what anybody actually wants -- q's tangent is the better of
	/// the two, being an approximation of one order higher.
	///
	/// The representation error of INVERSION-PLAN.md section 2 ( c ) is a
	/// property of the interpolation and is asked of the curve the interpolant
	/// is measured against, so it is this consistent pairing that carries the
	/// Delta_s^4 rate. The q pairing is reported beside it, along with the
	/// measured discrepancy, so that the tail has a number attached to it rather
	/// than being a mystery.
	struct TangentComparison
	{
		meq::Contour contour;
		double worstAngle = 0.0;
	};

	TangentComparison withGradientTangents( mfem::GridFunction const &potential,
	                                        meq::Contour const &in )
	{
		TangentComparison out;
		out.contour = in;

		mfem::Mesh &mesh = *potential.FESpace()->GetMesh();

		for ( std::size_t i = 0; i < out.contour.points.size(); ++i )
		{
			meq::ContourPoint &p = out.contour.points[ i ];

			mfem::Vector point( 2 );
			point( 0 ) = p.r;
			point( 1 ) = p.z;

			mfem::ElementTransformation *trans =
				mesh.GetElementTransformation( p.element );
			mfem::IntegrationPoint ip;
			trans->TransformBack( point, ip );
			trans->SetIntPoint( &ip );

			mfem::Vector gradient( 2 );
			potential.GetGradient( *trans, gradient );

			double const magnitude = std::sqrt( gradient( 0 )*gradient( 0 )
			                                    + gradient( 1 )*gradient( 1 ) );
			if ( !( magnitude > 0.0 ) )
				continue;

			double tangentR = -gradient( 1 )/magnitude;
			double tangentZ = gradient( 0 )/magnitude;

			// grad_bar( psi ) = r q with r > 0, so the two agree in orientation
			// as well as nearly in direction. The guard is here so that a sign
			// convention changing under this file fails loudly rather than
			// silently reversing half the tangents.
			double const alignment = tangentR*p.tangentR + tangentZ*p.tangentZ;
			BOOST_TEST( alignment > 0.0,
			            "the tangent from grad( psi_h ) points the opposite way "
			            << "to the tangent from q at ( " << p.r << ", " << p.z
			            << " ), which means one of the two sign conventions has "
			            << "moved" );

			double const cross = p.tangentR*tangentZ - p.tangentZ*tangentR;
			out.worstAngle = std::max( out.worstAngle, std::abs( cross ) );

			p.tangentR = tangentR;
			p.tangentZ = tangentZ;
		}
		return out;
	}

	/// A tracer set up for a fixed step: curvature control off and the local
	/// ceiling off, so that Delta_s is exactly what was asked for and the mesh
	/// cannot silently override the quantity a study is sweeping.
	void useFixedStep( meq::ContourTracer &tracer, double deltaS )
	{
		tracer.setTargetTurn( 0.0 );
		tracer.setLocalStepCeiling( 0.0 );
		tracer.setStep( deltaS );
	}

	/// Slack on a rate asserted across a whole sequence. Pointwise quantities
	/// here are noisier than the integrated norms SolovievConvergence.cpp
	/// asserts on -- a midpoint or a traced point sits somewhere different
	/// inside its element at every refinement, so the constant in front of the
	/// power is not the same constant at every level. The two-tier pattern of
	/// ExtensionConvergence.cpp: loose per pair, tight across the sequence.
	double const sequenceSlack = 0.30;
	double const pairSlack = 0.70;

	/// The standard benchmark box, and the fixture every rate below is measured
	/// on. nstx() is the corrected coefficient set;
	/// SolovievGeometryConvergence.cpp is what checks it.
	Rectangle benchmarkBox()
	{
		return meq::tests::standardBox();
	}

	/*
	 * ---------------------------------------------------------------------
	 * THE CURVED PATH AND THE BAND: IN-0's SECOND ACCEPTANCE.
	 * ---------------------------------------------------------------------
	 *
	 * Everything above is on the FITTED path, where Gamma is the mesh boundary,
	 * every rate is attributable to the tracer alone and there is no band. Below
	 * is the other half: Omega_h is the union of background elements lying
	 * INSIDE Gamma, so Gamma_h is inscribed and there is a band O( h ) wide that
	 * is inside the plasma and outside the mesh. The outermost closed surface is
	 * psi = 0, which IS Gamma, so the surfaces the band affects are exactly
	 * those a q( psi ) profile most wants.
	 *
	 * THE FIXTURE IS ExtensionConvergence.cpp's, DELIBERATELY, so that the two
	 * files measure different things about the same geometry and the numbers
	 * here can be read beside the rates there. Gamma is the closed flux surface
	 * psi_nstx = -0.03 written as the zero set of psi := psi_nstx + 0.03: the
	 * published separatrix is not closed with these coefficients -- there is a
	 * saddle at psi = -8.7e-3 -- and it would pass through an X-point, which is
	 * a CORNER of Gamma where both path families give out. That whole argument
	 * is in ExtensionConvergence.cpp's header and is not repeated.
	 *
	 * WHAT IS MEASURED, AND WHY IT IS TWO POPULATIONS AND NOT ONE. A contour
	 * near Gamma weaves in and out of the staircase Gamma_h, so its points split
	 * into those inside Omega_h -- limited by the discretisation, at psi_h's own
	 * k+1 -- and those in the band, limited by the extension. Quoting one rate
	 * over both would hide exactly the thing this stage exists to measure, which
	 * is the .nc's own lesson: a band node holds real data and is
	 * indistinguishable from a solved one without being told. So the band error
	 * and the interior error are separate columns at every point of every sweep.
	 */

	/// The additive constant that makes Gamma a closed flux surface. See above.
	double const curvedOffset = 0.03;

	Equilibrium const &curvedEquilibrium()
	{
		static Equilibrium const eq = Equilibrium::nstx();
		return eq;
	}

	double curvedPsi( double r, double z )
	{
		return curvedEquilibrium().psi( r, z ) + curvedOffset;
	}

	/// The level set the subdomain and the paths are built from: negative inside
	/// Omega. ExtensionConvergence.cpp measures that sign rather than assuming
	/// it, and this file inherits the answer.
	double curvedLevelSet( mfem::Vector const &x )
	{
		return curvedPsi( x( 0 ), x( 1 ) );
	}

	/// The background box: it contains Omega with room to spare, keeps r well
	/// away from zero and has sides in the ratio 1:2, so n by 2n cells are
	/// square.
	Rectangle curvedBox()
	{
		return Rectangle{ 0.25, 1.95, -1.75, 1.65 };
	}

	/// D_h: the elements of a uniform background triangulation that lie entirely
	/// inside Omega, as a SubMesh whose generated boundary is Gamma_h.
	std::unique_ptr<mfem::SubMesh> makeCurvedSubdomain( int n, int &gammaH,
	                                                    double &h )
	{
		Rectangle const box = curvedBox();

		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			n, 2*n, mfem::Element::TRIANGLE, false, box.width(), box.height() );
		background.Transform( [ box ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + box.rMin;
			out( 1 ) = in( 1 ) + box.zMin;
		} );
		h = box.width()/static_cast<double>( n );

		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background, curvedLevelSet,
		                                               0.0, marker, 1 );
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
		                    "D_h has boundary inherited from the background box "
		                    "at n = " << n << ", so part of Gamma_h is fitted" );
		return sub;
	}

	/// One curved solve, kept alive with its subdomain and its path family, so
	/// that a tracer can borrow all three. Member order is load bearing: the
	/// solver borrows the path, the path borrows the mesh, and destruction runs
	/// the other way.
	class CurvedSolve
	{
		public:
			CurvedSolve( int orderIn, int n )
				: sub( makeCurvedSubdomain( n, gammaH, meshSize ) ),
				  sourceCoeff( []( mfem::Vector const &x )
				  {
					  return curvedEquilibrium().f( x( 0 ), x( 1 ), 0.0 );
				  } ),
				  zero( 0.0 ),
				  // Six h of search length: the paths are about 1.3 h long, so
				  // this is a factor of four of slack, and it is what
				  // ExtensionConvergence.cpp uses.
				  path( *sub, gammaH, curvedLevelSet, 6.0*meshSize ),
				  solver( *sub, orderIn )
			{
				marker.SetSize( gammaH );
				marker = 0;
				marker[ gammaH - 1 ] = 1;

				solver.setSource( sourceCoeff );
				solver.setBoundaryData( zero );
				solver.setExtension( path, marker );
				solver.solve();
			}

			CurvedSolve( CurvedSolve const & ) = delete;
			CurvedSolve &operator=( CurvedSolve const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			mfem::SubMesh &mesh()
			{
				return *sub;
			}

			mfem::VertexConePath const &paths() const
			{
				return path;
			}

			mfem::Array<int> const &gammaHMarker() const
			{
				return marker;
			}

			double h() const
			{
				return meshSize;
			}

			meq::CriticalPoint axis()
			{
				meq::CriticalPointFinder finder( solver );
				return finder.findAxis();
			}

		private:
			int gammaH = 0;
			double meshSize = 0.0;
			std::unique_ptr<mfem::SubMesh> sub;
			mfem::Array<int> marker;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::ConstantCoefficient zero;
			mfem::VertexConePath path;
			meq::GradShafranovSolver solver;
	};

	/// Configure @a tracer for one of the two band extensions. One place, so
	/// that the side-by-side comparison differs in the extension and in nothing
	/// else -- which is what makes it a comparison.
	void useBandExtension( meq::ContourTracer &tracer, CurvedSolve &solved,
	                       meq::BandExtension which )
	{
		tracer.setBandExtension( which, solved.gammaHMarker(), &solved.paths() );
	}

	/// The first-order distance from a point to the exact contour of the CURVED
	/// fixture. The offset moves psi and not its gradient, so gradPsi is the
	/// equilibrium's own.
	double distanceToCurvedContour( double level, double r, double z )
	{
		double gr = 0.0;
		double gz = 0.0;
		curvedEquilibrium().gradPsi( r, z, gr, gz );
		double const magnitude = std::sqrt( gr*gr + gz*gz );
		return magnitude > 0.0 ? std::abs( curvedPsi( r, z ) - level )/magnitude
		                       : std::numeric_limits<double>::infinity();
	}

	/// The two populations of INVERSION-PLAN.md section 4.3, kept apart.
	struct BandMeasurement
	{
		double bandWorst = 0.0;
		double bandRms = 0.0;
		double interiorWorst = 0.0;
		int bandPoints = 0;
		int interiorPoints = 0;
		double deepest = 0.0;
		int stalled = 0;
		meq::ContourStatus status = meq::ContourStatus::Stalled;
	};

	BandMeasurement splitByBand( meq::Contour const &contour )
	{
		BandMeasurement out;
		out.status = contour.status;
		out.deepest = contour.deepestBandPoint;
		out.stalled = contour.stalledCorrections;

		double squares = 0.0;

		// The last point of a closed contour IS the first, exactly, so counting
		// it would weight one point twice. It changes no rate and it would make
		// the counts disagree with Contour::extendedPoints, which is the number
		// a reader checks these against.
		std::size_t const n = contour.closed() && contour.points.size() > 1
			? contour.points.size() - 1 : contour.points.size();

		for ( std::size_t i = 0; i < n; ++i )
		{
			meq::ContourPoint const &p = contour.points[ i ];
			double const error = distanceToCurvedContour( contour.level, p.r, p.z );
			if ( p.extended )
			{
				++out.bandPoints;
				out.bandWorst = std::max( out.bandWorst, error );
				squares += error*error;
			}
			else
			{
				++out.interiorPoints;
				out.interiorWorst = std::max( out.interiorWorst, error );
			}
		}

		out.bandRms = out.bandPoints > 0
			? std::sqrt( squares/out.bandPoints ) : 0.0;
		return out;
	}

	/// The flux label of a contour, as a fraction of the way from the exact axis
	/// to Gamma. Gamma is psi = 0 by construction, so 1 IS the boundary and the
	/// interesting levels are the ones just short of it.
	double curvedLevel( double fraction )
	{
		static ExactAxis const axis = exactAxis( curvedEquilibrium(), 1.322, 0.008 );
		double const axisPsi = axis.psi + curvedOffset;
		return axisPsi + fraction*( 0.0 - axisPsi );
	}

}

/*
 * THE INTERPOLANT ON ITS OWN, BEFORE ANY SOLVE.
 *
 * Contour's Hermite evaluators are on the struct rather than on the tracer
 * precisely so that this is possible: a Contour built by hand from an EXACT
 * circle, with exact unit tangents, isolates the representation from the
 * discretisation, from the corrector and from the step controller. If this
 * fails, nothing below means anything; if it passes and a later table does not,
 * the interpolant is not the suspect.
 *
 * The expected constant is not free either. The header derives the midpoint
 * deviation of the chord-scaled cubic on a unit circle as
 * 1 - cos( theta/2 ) - ( 1/2 ) sin^2( theta/2 ) = theta^4/128 + O( theta^6 ),
 * and that closed form is asserted against here as well as the rate, because a
 * rate alone cannot see a scheme that is fourth order about the wrong curve.
 */
BOOST_AUTO_TEST_CASE( theCubicHermiteIsFourthOrderAndTheChordsAreSecond )
{
	double const radius = 0.3;

	std::printf( "\n  Cubic Hermite from exact points and exact unit tangents, "
	             "circle of radius %.2f\n", radius );
	std::printf( "  %6s %10s %14s %7s %14s %7s %14s\n",
	             "n", "Delta_s", "Hermite", "rate", "chord", "rate", "theta^4/128" );

	std::vector<int> const counts = { 8, 16, 32, 64, 128 };
	std::vector<double> hermite;
	std::vector<double> chord;

	for ( std::size_t c = 0; c < counts.size(); ++c )
	{
		int const n = counts[ c ];
		double const dTheta = twoPi/n;

		meq::Contour circle;
		circle.level = 0.0;
		circle.status = meq::ContourStatus::Closed;
		for ( int i = 0; i <= n; ++i )
		{
			double const theta = i*dTheta;
			meq::ContourPoint p;
			p.r = 1.0 + radius*std::cos( theta );
			p.z = radius*std::sin( theta );
			p.tangentR = -std::sin( theta );
			p.tangentZ = std::cos( theta );
			p.fluxMagnitude = 1.0;
			p.arcLength = radius*theta;
			p.element = 0;
			circle.points.push_back( p );
		}

		double worstHermite = 0.0;
		double worstChord = 0.0;
		for ( std::size_t i = 0; i + 1 < circle.points.size(); ++i )
		{
			double r = 0.0;
			double z = 0.0;

			circle.pointOnSegment( i, 0.5, r, z );
			worstHermite = std::max( worstHermite,
				std::abs( std::sqrt( ( r - 1.0 )*( r - 1.0 ) + z*z ) - radius ) );

			circle.chordOnSegment( i, 0.5, r, z );
			worstChord = std::max( worstChord,
				std::abs( std::sqrt( ( r - 1.0 )*( r - 1.0 ) + z*z ) - radius ) );
		}

		hermite.push_back( worstHermite );
		chord.push_back( worstChord );

		double const predicted = radius*dTheta*dTheta*dTheta*dTheta/128.0;
		if ( c == 0 )
			std::printf( "  %6d %10.5f %14.6e %7s %14.6e %7s %14.6e\n",
			             n, radius*dTheta, worstHermite, "-", worstChord, "-", predicted );
		else
			std::printf( "  %6d %10.5f %14.6e %7.3f %14.6e %7.3f %14.6e\n",
			             n, radius*dTheta, worstHermite,
			             rate( hermite[ c - 1 ], worstHermite, 2.0 ), worstChord,
			             rate( chord[ c - 1 ], worstChord, 2.0 ), predicted );

		// The closed form, not merely the order. Ten percent, because the O( theta^6 )
		// term is not negligible at n = 8 and is what the loose end of this is for.
		BOOST_TEST( std::abs( worstHermite - predicted ) < 0.10*predicted + 1.0e-15,
		            "n = " << n << ": the Hermite midpoint deviation is "
		            << worstHermite << " where the closed form of the header says "
		            << predicted );
	}
	std::fflush( stdout );

	double const hermiteRate = rate( hermite.front(), hermite.back(),
	                                 static_cast<double>( counts.back() )/counts.front() );
	double const chordRate = rate( chord.front(), chord.back(),
	                               static_cast<double>( counts.back() )/counts.front() );

	BOOST_TEST( hermiteRate > 4.0 - sequenceSlack,
	            "the cubic Hermite converges at " << hermiteRate
	            << ", not the fourth order the chord-scaled tangents give" );
	bool const chordsAreSecondOrder = chordRate > 2.0 - sequenceSlack
	                                  && chordRate < 2.0 + sequenceSlack;
	BOOST_TEST( chordsAreSecondOrder,
	            "the straight chords converge at " << chordRate
	            << " rather than at second order. THE CHORDS ARE THE CONTROL: if "
	            << "they ever converge as fast as the Hermite the comparison is "
	            << "empty and this test is worthless" );
}

/*
 * THE TRACER RUNS, CLOSES, AND DOES NOT REACH FOR Mesh::FindPoints.
 *
 * Three separate properties, and the third is the one nothing else would
 * notice. CLAUDE.md records FindPoints as O( elements x points ), so a tracer
 * that falls back on it once per corrector iteration is quadratic in the mesh
 * and still returns the right answer -- there is no wrong number to see. The
 * count is the only thing that says the element walk is working.
 */
BOOST_AUTO_TEST_CASE( theTracerClosesAndTheElementWalkDoesNotFallBack )
{
	struct Case
	{
		char const *name;
		Equilibrium eq;
	};

	// TWO FIXTURES, AND THE SECOND ONE IS NOT DECORATION. nstx() puts its axis
	// at r = 1.318 on a box that ends at 1.4, so an outer surface reaches the
	// edge of the mesh along +r while sitting comfortably inside it everywhere
	// else -- which is what made traceFromAxis() search eight rays rather than
	// one. iterExample2() is the well-centred one, at r = 1.051, and is what the
	// rate sweeps below are measured on.
	std::vector<Case> cases = {
		{ "iterExample2", Equilibrium::iterExample2() },
		{ "nstx",         Equilibrium::nstx() }
	};

	Rectangle const box = benchmarkBox();

	for ( std::size_t c = 0; c < cases.size(); ++c )
	{
		Equilibrium const &eq = cases[ c ].eq;
		LevelSet const set = levelsFor( eq, box, { 0.05, 0.4, 0.8 } );

		std::printf( "\n  Traced contours, Solov'ev %s, k = 2, n = 24\n",
		             cases[ c ].name );
		std::printf( "  exact axis ( %.10f, %.10f ) psi %+.6e, edge psi %+.6e\n",
		             set.axis.r, set.axis.z, set.axis.psi, set.edge );
		std::printf( "  %6s %12s %6s %9s %10s %9s %8s %8s %9s %9s\n",
		             "Psi_N", "level", "points", "status", "length", "turning",
		             "circuits", "faces", "worst |r|", "fallback" );

		SolvedEquilibrium solved( eq, box, 2, 24 );
		meq::CriticalPoint const axis = solved.axis();

		for ( std::size_t i = 0; i < set.levels.size(); ++i )
		{
			meq::ContourTracer tracer( solved.theSolver() );
			meq::Contour const contour = tracer.traceFromAxis( set.levels[ i ], axis );

			std::printf( "  %6.2f %12.5e %6d %9s %10.6f %9.5f %8d %8d %9.2e %9d\n",
			             set.fractions[ i ], set.levels[ i ],
			             static_cast<int>( contour.points.size() ),
			             meq::contourStatusName( contour.status ),
			             contour.hermiteLength(), contour.turning, contour.circuits,
			             contour.faceCrossings, contour.worstResidual,
			             contour.fallbackLocations );

			BOOST_TEST( contour.closed(),
			            cases[ c ].name << " at Psi_N = " << set.fractions[ i ]
			            << " ended as '" << meq::contourStatusName( contour.status )
			            << "' rather than closing. A partial curve labelled as a "
			            << "contour is exactly what v0-legacy:FluxSurfaces.cpp "
			            << "returned and what ContourStatus exists to prevent" );

			BOOST_TEST( contour.circuits == 1,
			            cases[ c ].name << " at Psi_N = " << set.fractions[ i ]
			            << " turned " << contour.turning << " rad, which is "
			            << contour.circuits << " circuits rather than one" );

			BOOST_TEST( contour.fallbackLocations == 0,
			            cases[ c ].name << ": the element walk fell back on "
			            << "Mesh::FindPoints " << contour.fallbackLocations
			            << " times at Psi_N = " << set.fractions[ i ]
			            << ". That is not a wrong answer, it is the trace going "
			            << "quadratic in the mesh, and nothing but this count "
			            << "would say so" );

			// The last point IS the first point, exactly: the final step is
			// shortened to land on it rather than leaving a stub.
			BOOST_TEST( contour.points.front().r == contour.points.back().r,
			            "the closing point is not the start point" );
			BOOST_TEST( contour.points.front().z == contour.points.back().z,
			            "the closing point is not the start point" );
		}
	}
	std::fflush( stdout );
}

/*
 * MEASUREMENT ( c ): THE REPRESENTATION ERROR, AGAINST THE DISCRETE LEVEL SET.
 *
 * Delta_s is swept on ONE fixed mesh, so h and k are held and what is left is
 * the interpolation alone. Three columns, and reading them in the right order
 * is the whole content of the test.
 *
 * THE BRIEF FOR THIS STAGE EXPECTED TWO COLUMNS AND THE MEASUREMENT NEEDED
 * THREE. The expectation was that a cubic Hermite through traced points, with
 * the tangents taken from q, would be O( Delta_s^4 ) against { psi_h = c }. It
 * is not, and the reason is the mixed method rather than the interpolation:
 * q_h and grad( psi_h )/r are DIFFERENT OBJECTS, converging at k+1 and at k
 * respectively, and the tangent of { psi_h = c } is built from the second. So
 * the q interpolant passes through points of { psi_h = c } carrying the
 * tangents of a slightly different curve, and measured against { psi_h = c } it
 * is fourth order until that O( h^k ) mismatch takes over and second order
 * afterwards. See withGradientTangents() above.
 *
 * The representation error is a property of an interpolation and is asked
 * relative to the curve it is measured against, so the CONSISTENT pairing is
 * what carries the rate: the same traced points with tangents from
 * grad( psi_h ). The q pairing is printed beside it with the measured tangent
 * discrepancy, so that its tail has a number attached to it, and the chords are
 * printed beside both as the control.
 *
 * NONE OF THIS ARGUES AGAINST TAKING THE TANGENT FROM q. Against the TRUE
 * contour -- which is what a flux-surface average is actually taken over -- q's
 * tangent is the better of the two by a full order, and the field error ( a )
 * measured two tests below is larger than either interpolant's deviation at
 * every mesh in the sweep. What has been found is that the DISCRETE level set
 * is not quite the right thing to measure a q-tangent interpolant against, and
 * that is worth a column rather than a footnote.
 *
 * IT ALSO HAS A FLOOR AND THE FLOOR IS PRINTED BESIDE IT. psi_h is
 * discontinuous across faces, so a segment whose two endpoints are corrected
 * onto the arcs of two different elements cannot be closer to either arc than
 * those arcs are to each other. That distance is the DG jump over | grad psi |,
 * it converges at k+1 in h and it is a component of error ( a ). The clean
 * columns exclude those segments, and representationError() above says what
 * else "clean" has to exclude and why.
 */
BOOST_AUTO_TEST_CASE( theRepresentationConvergesAtFourthOrderWithHermite )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	int const order = 3;
	int const cells = 32;

	SolvedEquilibrium solved( eq, box, order, cells );
	meq::CriticalPoint const axis = solved.axis();

	struct Pairing
	{
		char const *name;
		meq::Potential which;
		mfem::GridFunction const *potential;
	};

	std::vector<Pairing> pairings = {
		{ "psi_h with q_h", meq::Potential::Raw,
		  &solved.theSolver().potential() },
		{ "psi* with q*",   meq::Potential::PostProcessed,
		  &solved.theSolver().postProcessedPotential() }
	};

	// A first trace at the default step, only to learn the contour's length so
	// that the sweep can be expressed in points per circuit rather than in an
	// arbitrary number.
	double contourLength = 0.0;
	{
		meq::ContourTracer tracer( solved.theSolver() );
		contourLength = tracer.traceFromAxis( set.levels[ 0 ], axis ).hermiteLength();
	}

	std::vector<int> const perTurn = { 128, 256, 512, 1024 };
	double const ratio = static_cast<double>( perTurn.back() )/perTurn.front();

	std::vector<double> fluxRates;
	std::vector<double> worstTilt;

	for ( std::size_t p = 0; p < pairings.size(); ++p )
	{
		std::printf( "\n  Representation error ( c ) against the level set of "
		             "the field being rooted,\n  %s, k = %d, n = %d, "
		             "Psi_N = %.2f, contour length %.6f, h = %.5f\n",
		             pairings[ p ].name, order, cells, set.fractions[ 0 ],
		             contourLength, solved.meshSize() );
		std::printf( "  worst over segments that stay inside one element with "
		             "both ends at tolerance\n" );
		std::printf( "  %8s %11s %8s %14s %7s %14s %7s %14s %7s %10s %10s\n",
		             "per turn", "Delta_s", "points", "grad-tangent", "rate",
		             "flux-tangent", "rate", "chord", "rate", "|t_q x t_g|",
		             "jump" );

		std::vector<double> gradient;
		std::vector<double> flux;
		std::vector<double> chord;
		double tilt = 0.0;

		for ( std::size_t i = 0; i < perTurn.size(); ++i )
		{
			double const deltaS = contourLength/perTurn[ i ];

			meq::ContourTracer tracer( solved.theSolver(), pairings[ p ].which );
			useFixedStep( tracer, deltaS );
			meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );

			BOOST_TEST( contour.closed(),
			            pairings[ p ].name << ": the contour did not close at "
			            << "Delta_s = " << deltaS );

			TangentComparison const consistent =
				withGradientTangents( *pairings[ p ].potential, contour );

			RepresentationError const fromFlux = representationError( tracer, contour );
			RepresentationError const fromGradient =
				representationError( tracer, consistent.contour );

			gradient.push_back( fromGradient.cleanHermite );
			flux.push_back( fromFlux.cleanHermite );
			chord.push_back( fromFlux.cleanChord );
			tilt = std::max( tilt, consistent.worstAngle );

			double const scale = contour.points.front().r
			                     *contour.points.front().fluxMagnitude;
			double const floor = scale > 0.0 ? contour.worstFaceJump/scale : 0.0;

			if ( i == 0 )
				std::printf( "  %8d %11.3e %8d %14.6e %7s %14.6e %7s %14.6e %7s "
				             "%10.2e %10.2e\n",
				             perTurn[ i ], deltaS,
				             static_cast<int>( contour.points.size() ),
				             gradient[ i ], "-", flux[ i ], "-", chord[ i ], "-",
				             consistent.worstAngle, floor );
			else
				std::printf( "  %8d %11.3e %8d %14.6e %7.3f %14.6e %7.3f %14.6e "
				             "%7.3f %10.2e %10.2e\n",
				             perTurn[ i ], deltaS,
				             static_cast<int>( contour.points.size() ),
				             gradient[ i ], rate( gradient[ i - 1 ], gradient[ i ], 2.0 ),
				             flux[ i ], rate( flux[ i - 1 ], flux[ i ], 2.0 ),
				             chord[ i ], rate( chord[ i - 1 ], chord[ i ], 2.0 ),
				             consistent.worstAngle, floor );
		}
		std::fflush( stdout );

		double const gradientRate = rate( gradient.front(), gradient.back(), ratio );
		double const chordRate = rate( chord.front(), chord.back(), ratio );
		fluxRates.push_back( rate( flux.front(), flux.back(), ratio ) );
		worstTilt.push_back( tilt );

		// THE GRAD-TANGENT COLUMN IS THE REPRESENTATION ERROR PROPER, and it
		// must be fourth order for BOTH pairings: it is the interpolation of a
		// smooth arc by two points and two of its own tangents, and the field
		// the arc came from does not enter.
		BOOST_TEST( gradientRate > 4.0 - sequenceSlack,
		            pairings[ p ].name << ": the cubic Hermite converges at "
		            << gradientRate << " in Delta_s against the curve whose "
		            << "tangents it was given, not the fourth order two points "
		            << "and two tangents buy" );

		bool const chordsAreSecondOrder = chordRate > 2.0 - sequenceSlack
		                                  && chordRate < 2.0 + sequenceSlack;
		BOOST_TEST( chordsAreSecondOrder,
		            pairings[ p ].name << ": the straight chords converge at "
		            << chordRate << " rather than at second order. THE CHORDS "
		            << "ARE THE CONTROL: if they ever converge as fast as the "
		            << "Hermite, the fourth-order claim has nothing to be "
		            << "measured against and this test is worthless" );

		BOOST_TEST( gradient.back() < 0.005*chord.back(),
		            pairings[ p ].name << ": at the finest spacing the Hermite "
		            << "error is " << gradient.back() << " against the chords' "
		            << chord.back() << ", which is not the separation two orders "
		            << "should give" );
	}

	std::printf( "\n  the flux-tangent interpolant: %s reaches %.3f with a "
	             "tangent tilt of %.2e,\n  %s reaches %.3f with %.2e\n",
	             pairings[ 0 ].name, fluxRates[ 0 ], worstTilt[ 0 ],
	             pairings[ 1 ].name, fluxRates[ 1 ], worstTilt[ 1 ] );
	std::fflush( stdout );

	// AND THE PAIRING IS WHAT DECIDES WHETHER THE TRACER'S OWN INTERPOLANT --
	// the one built on q, which is what a caller receives -- KEEPS THE FOURTH
	// ORDER. With psi* it does; with psi_h it does not, and the tilt column says
	// by how much and why.
	BOOST_TEST( fluxRates[ 1 ] > 4.0 - sequenceSlack,
	            "with psi* the Hermite from q converges at " << fluxRates[ 1 ]
	            << ". grad( psi* ) and q* agree to " << worstTilt[ 1 ]
	            << ", which is a full order better than grad( psi_h ) against "
	            << "q_h, and that is why the default pairing keeps the order "
	            << "that the tangent from q was chosen for" );

	BOOST_TEST( worstTilt[ 0 ] > 10.0*worstTilt[ 1 ],
	            "grad( psi_h ) tilts from q_h by " << worstTilt[ 0 ]
	            << " and grad( psi* ) from q* by " << worstTilt[ 1 ]
	            << ". THAT SEPARATION IS THE CONTROL for the row above: if the "
	            << "two pairings were equally consistent there would be nothing "
	            << "to explain the difference in their flux-tangent rates and "
	            << "this table would be measuring something else" );

	BOOST_TEST( fluxRates[ 0 ] < fluxRates[ 1 ] - 0.5,
	            "the flux-tangent interpolant converges at " << fluxRates[ 0 ]
	            << " on psi_h and " << fluxRates[ 1 ] << " on psi*. The whole "
	            << "point of the pair of tables is that the first is worse; if "
	            << "it stops being worse the tangent-consistency account in "
	            << "FluxSurfaces.hpp is wrong" );
}

/*
 * AND AT FIXED Delta_s IT DOES NOT KNOW ABOUT h, WHICH IS WHAT SAYS ( b ) AND
 * ( c ) HAVE BEEN SEPARATED.
 *
 * The representation error is a property of the point spacing and the
 * interpolation order and of nothing else -- INVERSION-PLAN.md section 2 calls
 * that the single most important realisation in its survey, because it is what
 * defeats the obvious objection that exact points joined by chords still give a
 * second-order line integral. Hold Delta_s and refine the mesh: if the number
 * moved with h, the interpolant would be inheriting the mesh through the point
 * PLACEMENT and the two errors would not be separate at all.
 *
 * THE q COLUMN IS PRINTED AND IS NOT FLAT, AND THAT IS THE SAME FINDING AS THE
 * TEST ABOVE RATHER THAN A SECOND ONE. Its extra term is the O( h^k ) gap
 * between q_h and grad( psi_h )/r, which does move with the mesh -- it falls,
 * at about one order here -- because it is a discretisation error and not a
 * representation error. The consistent column is the representation error and
 * it is the one that must not move.
 */
BOOST_AUTO_TEST_CASE( theRepresentationErrorDoesNotKnowAboutTheMesh )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	int const order = 2;
	std::vector<int> const meshes = { 12, 24, 48 };
	double const deltaS = 0.004;

	std::printf( "\n  Representation error ( c ) at FIXED Delta_s = %.4f, k = %d\n",
	             deltaS, order );
	std::printf( "  %6s %9s %8s %14s %14s %14s %10s %10s %9s %7s\n",
	             "n", "h", "points", "grad psi*", "from q*", "chord",
	             "|t_q x t_g|", "jump", "clean", "stalled" );

	std::vector<double> gradient;
	std::vector<double> chord;

	for ( std::size_t i = 0; i < meshes.size(); ++i )
	{
		SolvedEquilibrium solved( eq, box, order, meshes[ i ] );
		meq::CriticalPoint const axis = solved.axis();

		meq::ContourTracer tracer( solved.theSolver() );
		useFixedStep( tracer, deltaS );
		meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );

		BOOST_TEST( contour.closed(), "the contour did not close at n = " << meshes[ i ] );

		TangentComparison const consistent =
			withGradientTangents( solved.theSolver().postProcessedPotential(),
			                      contour );

		RepresentationError const fromFlux = representationError( tracer, contour );
		RepresentationError const fromGradient =
			representationError( tracer, consistent.contour );

		gradient.push_back( fromGradient.cleanHermite );
		chord.push_back( fromFlux.cleanChord );

		double const scale = contour.points.front().r
		                     *contour.points.front().fluxMagnitude;
		double const floor = scale > 0.0 ? contour.worstFaceJump/scale : 0.0;

		char clean[ 32 ];
		std::snprintf( clean, sizeof( clean ), "%d/%d", fromFlux.clean,
		               fromFlux.total );

		std::printf( "  %6d %9.5f %8d %14.6e %14.6e %14.6e %10.2e %10.2e %9s %7d\n",
		             meshes[ i ], solved.meshSize(),
		             static_cast<int>( contour.points.size() ),
		             fromGradient.cleanHermite, fromFlux.cleanHermite,
		             fromFlux.cleanChord, consistent.worstAngle, floor, clean,
		             contour.stalledCorrections );
	}
	std::fflush( stdout );

	double const worstHermite = *std::max_element( gradient.begin(), gradient.end() );
	double const bestHermite = *std::min_element( gradient.begin(), gradient.end() );
	double const worstChord = *std::max_element( chord.begin(), chord.end() );
	double const bestChord = *std::min_element( chord.begin(), chord.end() );

	// A factor of two across a sixteenfold change in dof count. A quantity that
	// tracked h would move by 4^3 = 64 at k = 2 over this sweep.
	BOOST_TEST( worstHermite < 2.0*bestHermite,
	            "the Hermite representation error moved from " << bestHermite
	            << " to " << worstHermite << " over a fourfold refinement at "
	            << "fixed Delta_s. It should not know about h at all: if it "
	            << "does, ( b ) and ( c ) are not separate and the whole error "
	            << "split of INVERSION-PLAN.md section 2 is not what this code "
	            << "implements" );

	BOOST_TEST( worstChord < 2.0*bestChord,
	            "the chord representation error moved from " << bestChord
	            << " to " << worstChord << " at fixed Delta_s" );
}

/*
 * MEASUREMENT ( a ): THE FIELD ERROR, WHICH IS THE DISCRETISATION AND NOTHING
 * TO DO WITH THE TRACER.
 *
 * The traced points are on { psi_h = c } to the corrector tolerance, so their
 * distance from { psi = c } IS the pointwise error of psi_h divided by
 * | grad psi | -- the implicit function theorem, exactly as
 * CriticalPointConvergence.cpp uses it for the axis. It converges at k+1.
 *
 * Delta_s is held fixed across the sweep so that the representation error
 * cannot contaminate it, and it is held SMALL, so that the quantity being
 * measured is the position of the accepted points rather than the curve between
 * them.
 */
BOOST_AUTO_TEST_CASE( theTracedContourConvergesAtTheFieldsOwnOrder )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.05, 0.4, 0.8 } );

	std::vector<int> const orders = { 1, 2, 3 };

	// NOT { 8, 16, 32 }, AND THE COARSEST POINT IS WHY. At k = 1 on h = 0.1 the
	// innermost contour of this equilibrium is not resolved at all -- the trace
	// runs out of the mesh and the field error reads 6.3e-1 on a contour 0.15
	// across. That is the standing objection this project keeps rediscovering
	// and CLAUDE.md records against its own self-convergence studies: a sweep
	// that starts outside the asymptotic regime measures the approach to it
	// rather than the rate.
	std::vector<int> const meshes = { 16, 32, 64 };
	double const deltaS = 0.01;

	for ( std::size_t o = 0; o < orders.size(); ++o )
	{
		int const order = orders[ o ];
		std::printf( "\n  Field error ( a ): traced { psi_h = c } against "
		             "{ psi = c }, k = %d, Delta_s = %.3f\n", order, deltaS );
		std::printf( "  %6s %9s", "n", "h" );
		for ( std::size_t i = 0; i < set.levels.size(); ++i )
			std::printf( "  %13s %6s", "dist", "rate" );
		std::printf( "\n" );

		std::vector<std::vector<double>> errors( set.levels.size() );

		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			SolvedEquilibrium solved( eq, box, order, meshes[ m ] );
			meq::CriticalPoint const axis = solved.axis();

			std::printf( "  %6d %9.5f", meshes[ m ], solved.meshSize() );

			for ( std::size_t i = 0; i < set.levels.size(); ++i )
			{
				// EXPLICITLY THE RAW PAIRING. k+1 is psi_h's order, and the
				// tracer's default is psi*, which carries k+2 into the contour --
				// measured in theRawPotentialAndThePostProcessedOneAreCompared.
				// Measuring ( a ) here against psi_h's order is what makes the
				// implicit-function-theorem argument checkable at all.
				meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
				useFixedStep( tracer, deltaS );
				meq::Contour const contour = tracer.traceFromAxis( set.levels[ i ], axis );

				BOOST_TEST( contour.closed(),
				            "k = " << order << ", n = " << meshes[ m ]
				            << ", Psi_N = " << set.fractions[ i ]
				            << ": the contour ended as '"
				            << meq::contourStatusName( contour.status ) << "'" );

				double const error = fieldError( eq, contour );
				errors[ i ].push_back( error );

				if ( m == 0 )
					std::printf( "  %13.6e %6s", error, "-" );
				else
					std::printf( "  %13.6e %6.3f", error,
					             rate( errors[ i ][ m - 1 ], error, 2.0 ) );
			}
			std::printf( "\n" );
		}
		std::fflush( stdout );

		for ( std::size_t i = 0; i < set.levels.size(); ++i )
		{
			double const whole = rate( errors[ i ].front(), errors[ i ].back(), 4.0 );
			BOOST_TEST( whole > order + 1.0 - sequenceSlack,
			            "k = " << order << ", Psi_N = " << set.fractions[ i ]
			            << ": the traced contour converges on the exact one at "
			            << whole << " rather than at the k+1 = " << order + 1
			            << " that psi_h's own order gives it through the "
			            << "implicit function theorem" );

			for ( std::size_t m = 1; m < meshes.size(); ++m )
				BOOST_TEST( rate( errors[ i ][ m - 1 ], errors[ i ][ m ], 2.0 )
				            > order + 1.0 - pairSlack,
				            "k = " << order << ", Psi_N = " << set.fractions[ i ]
				            << ": the pair " << meshes[ m - 1 ] << " -> "
				            << meshes[ m ] << " converges at "
				            << rate( errors[ i ][ m - 1 ], errors[ i ][ m ], 2.0 ) );
		}
	}
}

/*
 * MEASUREMENT ( b ): THE POINT ERROR DOES NOT ACCUMULATE, AND THAT IS THE WHOLE
 * REASON THIS FAMILY OF METHODS IS THE RIGHT ONE.
 *
 * Without a corrector, error off the contour accumulates with path length: the
 * curve drifts onto a neighbouring contour and does not close. With one, each
 * point is independently on the level set and there is no accumulation at all.
 * The measurement is one, five and ten circuits of the same contour --
 * ten times the path length, sixty-three radians of turning -- and neither the
 * worst residual nor the closure error may know about it.
 *
 * AND "CLOSURE AT MACHINE PRECISION" NEEDED A REFINEMENT, WHICH IS WHY THE
 * NORMAL COMPONENT IS WHAT IS ASSERTED. psi_h is discontinuous across faces, so
 * the level set is a union of per-element arcs and the curve is closed as a set
 * without being one analytic curve. Returning to the start element lands on the
 * SAME arc, so the NORMAL closure error is at the corrector tolerance;
 * tangentially it is whatever the final step leaves, which is why the final
 * step is shortened rather than taken at full length. Both are printed and the
 * normal one is asserted.
 */
BOOST_AUTO_TEST_CASE( thePointErrorDoesNotAccumulateOverManyCircuits )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	SolvedEquilibrium solved( eq, box, 2, 24 );
	meq::CriticalPoint const axis = solved.axis();

	std::printf( "\n  Error ( b ) against path length, k = 2, n = 24, Psi_N = %.2f\n",
	             set.fractions[ 0 ] );
	std::printf( "  %9s %8s %10s %9s %14s %8s %14s %14s\n",
	             "circuits", "points", "length", "turning", "worst |psi-c|",
	             "worst it", "closure normal", "closure tang." );

	std::vector<int> const circuits = { 1, 5, 10 };
	std::vector<double> residuals;
	std::vector<double> normals;
	std::vector<double> tangentials;
	std::vector<double> curvatures;

	for ( std::size_t i = 0; i < circuits.size(); ++i )
	{
		meq::ContourTracer tracer( solved.theSolver() );
		tracer.setCircuits( circuits[ i ] );
		meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );

		BOOST_TEST( contour.closed(),
		            circuits[ i ] << " circuits ended as '"
		            << meq::contourStatusName( contour.status ) << "'" );
		BOOST_TEST( contour.circuits == circuits[ i ],
		            "asked for " << circuits[ i ] << " circuits and turned "
		            << contour.turning << " rad, which is " << contour.circuits );

		residuals.push_back( contour.worstResidual );
		normals.push_back( contour.closureNormal );
		tangentials.push_back( contour.closureTangential );
		curvatures.push_back( twoPi*contour.circuits/contour.length() );

		std::printf( "  %9d %8d %10.5f %9.4f %14.6e %8d %14.6e %14.6e\n",
		             circuits[ i ], static_cast<int>( contour.points.size() ),
		             contour.length(), contour.turning, contour.worstResidual,
		             contour.worstCorrectorIterations, contour.closureNormal,
		             contour.closureTangential );
	}
	std::fflush( stdout );

	// Not "small": FLAT. A drifting tracer's residual would grow with the path,
	// and the whole argument for a corrector is that this one cannot.
	BOOST_TEST( residuals.back() < 4.0*residuals.front() + 1.0e-30,
	            "the worst | psi_h - c | went from " << residuals.front()
	            << " over one circuit to " << residuals.back() << " over "
	            << circuits.back() << ". It must not know about path length: "
	            << "that is the property the corrector exists for, and its "
	            << "absence would mean the tracer is integrating rather than "
	            << "projecting" );

	// AND THE NORMAL CLOSURE ERROR, WHICH NEEDED ITS OWN REFINEMENT AND IS THE
	// PLACE THIS FILE'S FIRST ACCOUNT WAS WRONG.
	//
	// It is NOT at the corrector tolerance and it cannot be. The final corrected
	// point and the start point are both on { psi_h = c } to 1e-13, but they are
	// separated ALONG the curve by whatever the shortened final step leaves, and
	// an arc of curvature kappa departs from its own tangent line by
	// kappa g_t^2 / 2 over a tangential gap g_t. So the normal component of the
	// gap is bounded by the curve BENDING, not by drift -- measured, it sits at
	// about a quarter of that bound.
	//
	// That is the discriminating statement, and it is what is asserted: a
	// tracer that was integrating rather than projecting would have a normal
	// closure error growing with path length and unrelated to g_t.
	for ( std::size_t i = 0; i < circuits.size(); ++i )
		BOOST_TEST( normals[ i ] <= 2.0*curvatures[ i ]*tangentials[ i ]
		                            *tangentials[ i ] + 1.0e-14,
		            circuits[ i ] << " circuits: the normal closure error is "
		            << normals[ i ] << " against the " << curvatures[ i ]
		            << " * " << tangentials[ i ] << "^2 that the curve's own "
		            << "bending over the tangential gap accounts for. An excess "
		            << "over that is drift, which is exactly what the corrector "
		            << "is there to prevent" );

	// AND IT DOES NOT KNOW ABOUT PATH LENGTH. Five circuits against ten -- the
	// path doubled and nothing else changed. One circuit is deliberately not the
	// comparison: the step controller starts cold there, so the final step aims
	// differently and the tangential gap it leaves is four times smaller, which
	// moves the bound above with it. Comparing against it would be comparing two
	// different final steps rather than two different path lengths.
	BOOST_TEST( normals.back() < 2.0*normals[ 1 ] + 1.0e-14,
	            "the NORMAL closure error went from " << normals[ 1 ]
	            << " over " << circuits[ 1 ] << " circuits to " << normals.back()
	            << " over " << circuits.back() << ". Doubling the path length "
	            << "must not move it: that is the property the corrector exists "
	            << "for" );
}

/*
 * THE FACE JUMP IS ERROR ( a ), MEASURED, AND IT IS NOT A TRACER DEFECT.
 *
 * A traced contour crosses faces, and the potential evaluated from the two
 * sides at the crossing disagrees by the DG jump. That number is worth
 * reporting for two reasons: it is the floor under measurement ( c ) above, so
 * a reader who does not know it is there will read a flattening Hermite column
 * as a defect in the interpolant; and it converges at the ORDER OF THE FIELD
 * BEING ROOTED, which is what says it is the discretisation. That field is
 * psi* by default, so the rate here is k+2 and not k+1 -- the same fact the
 * potential comparison below establishes from the direction of the contour,
 * arriving here from the direction of the jump.
 *
 * It is measured by BISECTING for the crossing point and evaluating from each
 * side INSIDE its own element -- never by extrapolating one element's
 * polynomial across the face, which is the trap CLAUDE.md records the .nc band
 * falling into.
 */
BOOST_AUTO_TEST_CASE( theFaceJumpIsTheDiscretisationAndConvergesAtItsOrder )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	std::vector<int> const orders = { 1, 2, 3 };
	std::vector<int> const meshes = { 8, 16, 32 };

	std::printf( "\n  The DG jump of psi* across the faces a contour crosses, "
	             "Psi_N = %.2f -- psi* because\n  that is the field the default "
	             "tracer roots, and the jump is the jump of THAT field\n",
	             set.fractions[ 0 ] );
	std::printf( "  THE MEAN IS WHAT CARRIES THE RATE: the worst jump over a "
	             "contour is a maximum over\n  a handful of crossings, so it is "
	             "a different sample at every refinement\n" );
	std::printf( "  %3s %6s %9s %8s %14s %7s %14s %7s\n",
	             "k", "n", "h", "faces", "mean jump", "rate", "worst jump", "rate" );

	for ( std::size_t o = 0; o < orders.size(); ++o )
	{
		std::vector<double> mean;
		std::vector<double> worst;
		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			SolvedEquilibrium solved( eq, box, orders[ o ], meshes[ m ] );
			meq::CriticalPoint const axis = solved.axis();

			meq::ContourTracer tracer( solved.theSolver() );
			useFixedStep( tracer, 0.004 );
			meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );

			BOOST_TEST( contour.closed(),
			            "k = " << orders[ o ] << ", n = " << meshes[ m ]
			            << ": the contour ended as '"
			            << meq::contourStatusName( contour.status ) << "'" );

			mean.push_back( contour.meanFaceJump );
			worst.push_back( contour.worstFaceJump );

			if ( m == 0 )
				std::printf( "  %3d %6d %9.5f %8d %14.6e %7s %14.6e %7s\n",
				             orders[ o ], meshes[ m ], solved.meshSize(),
				             contour.faceCrossings, contour.meanFaceJump, "-",
				             contour.worstFaceJump, "-" );
			else
				std::printf( "  %3d %6d %9.5f %8d %14.6e %7.3f %14.6e %7.3f\n",
				             orders[ o ], meshes[ m ], solved.meshSize(),
				             contour.faceCrossings, contour.meanFaceJump,
				             rate( mean[ m - 1 ], contour.meanFaceJump, 2.0 ),
				             contour.worstFaceJump,
				             rate( worst[ m - 1 ], contour.worstFaceJump, 2.0 ) );

			BOOST_TEST( contour.faceCrossings > 0,
			            "k = " << orders[ o ] << ", n = " << meshes[ m ]
			            << ": a contour of length " << contour.length()
			            << " on a mesh of size " << solved.meshSize()
			            << " crossed no faces at all, so the jump measurement is "
			            << "measuring nothing" );
		}

		// k+2, NOT k+1, because the field being rooted is psi* -- which is the
		// same statement as the potential comparison two tests down, arriving
		// from the direction of the jump rather than of the contour.
		double const whole = rate( mean.front(), mean.back(), 4.0 );
		BOOST_TEST( whole > orders[ o ] + 2.0 - pairSlack,
		            "k = " << orders[ o ] << ": the DG jump of psi* converges at "
		            << whole << " rather than at the k+2 that is psi*'s own "
		            << "order. It is a component of error ( a ) and must fall "
		            << "with the discretisation" );
	}
	std::fflush( stdout );
}

/*
 * WHICH POTENTIAL TO ROOT: THE MEASUREMENT, NOT THE PREFERENCE.
 *
 * The corrector roots a scalar field and takes its direction from the flux, and
 * meq has two candidate pairs: psi_h with q_h, which converge at k+1, and
 * psi*_h with q*_h from DarcyForm::Reconstruct(), where psi* converges at k+2
 * and is what meq reports everywhere else.
 *
 * The reason to expect psi* to win is that the local post-processing is built
 * so that grad( psi* ) matches the flux, which would make the corrector's
 * direction consistent with the function being rooted. CHECK THAT BEFORE
 * WRITING IT DOWN AS A REASON -- and it does not hold as stated: Stenberg's
 * local problem is driven by the reconstructed TOTAL flux qhat_h in RT_k, not
 * by q_h and not by q*_h, so neither pairing is the consistent one. What is
 * reported here is therefore the measurement and not the argument: the field
 * error ( a ) of the contour each pair traces, and the corrector iterations
 * each costs.
 */
BOOST_AUTO_TEST_CASE( theRawPotentialAndThePostProcessedOneAreCompared )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	std::vector<int> const orders = { 1, 2, 3 };
	std::vector<int> const meshes = { 8, 16, 32 };
	double const deltaS = 0.01;

	std::printf( "\n  Which potential the corrector should root, Psi_N = %.2f, "
	             "Delta_s = %.3f\n", set.fractions[ 0 ], deltaS );
	std::printf( "  %3s %6s %14s %6s %8s %10s %14s %6s %8s %10s %9s\n",
	             "k", "n", "psi_h dist", "rate", "it/point", "|t_q x t_g|",
	             "psi* dist", "rate", "it/point", "|t_q x t_g|", "ratio" );

	for ( std::size_t o = 0; o < orders.size(); ++o )
	{
		std::vector<double> raw;
		std::vector<double> post;

		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			SolvedEquilibrium solved( eq, box, orders[ o ], meshes[ m ] );
			meq::CriticalPoint const axis = solved.axis();

			meq::ContourTracer rawTracer( solved.theSolver(), meq::Potential::Raw );
			useFixedStep( rawTracer, deltaS );
			meq::Contour const rawContour =
				rawTracer.traceFromAxis( set.levels[ 0 ], axis );

			meq::ContourTracer postTracer( solved.theSolver(),
			                               meq::Potential::PostProcessed );
			useFixedStep( postTracer, deltaS );
			meq::Contour const postContour =
				postTracer.traceFromAxis( set.levels[ 0 ], axis );

			double const rawError = fieldError( eq, rawContour );
			double const postError = fieldError( eq, postContour );
			raw.push_back( rawError );
			post.push_back( postError );

			double const rawCost = static_cast<double>( rawContour.correctorIterationsTotal )
			                       /rawContour.points.size();
			double const postCost = static_cast<double>( postContour.correctorIterationsTotal )
			                        /postContour.points.size();

			// AND THE CONSISTENCY CLAIM, CHECKED RATHER THAN REPEATED. The
			// reason to EXPECT psi* to win is that its local post-processing is
			// built so that its gradient matches the flux, which would make the
			// corrector's direction consistent with the function being rooted.
			// This is that claim as a number, for both pairings: the sine of the
			// angle between the tangent from the flux and the tangent from the
			// potential's own gradient.
			double const rawTilt =
				withGradientTangents( solved.theSolver().potential(),
				                      rawContour ).worstAngle;
			double const postTilt =
				withGradientTangents( solved.theSolver().postProcessedPotential(),
				                      postContour ).worstAngle;

			std::printf( "  %3d %6d %14.6e %6s %8.2f %10.2e %14.6e %6s %8.2f "
			             "%10.2e %9.2f\n",
			             orders[ o ], meshes[ m ], rawError,
			             m == 0 ? "-" : "", rawCost, rawTilt, postError,
			             m == 0 ? "-" : "", postCost, postTilt,
			             rawError/postError );

			if ( m > 0 )
			{
				std::printf( "  %3s %6s %14s %6.3f %8s %10s %14s %6.3f %8s %10s "
				             "%9s\n",
				             "", "", "", rate( raw[ m - 1 ], rawError, 2.0 ), "",
				             "", "", rate( post[ m - 1 ], postError, 2.0 ), "",
				             "", "" );
			}

			bool const bothClosed = rawContour.closed() && postContour.closed();
			BOOST_TEST( bothClosed,
			            "k = " << orders[ o ] << ", n = " << meshes[ m ]
			            << ": one of the two pairings did not close" );
		}

		// Both must converge: the question this test answers is WHICH IS
		// BETTER, and a pairing that did not converge at all would not be a
		// candidate. k+1 is the bar for the raw pair; the post-processed one is
		// asserted at the same bar rather than at k+2, because what is being
		// traced is the level set of psi* and the accuracy of a level set is
		// set by BOTH the function and the direction the corrector moves in.
		double const rawRate = rate( raw.front(), raw.back(), 4.0 );
		double const postRate = rate( post.front(), post.back(), 4.0 );

		BOOST_TEST( rawRate > orders[ o ] + 1.0 - sequenceSlack,
		            "k = " << orders[ o ] << ": rooting psi_h converges at "
		            << rawRate << " rather than at the k+1 that psi_h's own "
		            << "order gives it" );

		// AND psi* CARRIES ITS OWN ORDER INTO THE CONTOUR, which is the
		// measurement that settles the default. psi* converges at k+2, the
		// corrector puts every point on { psi* = c }, and the implicit function
		// theorem hands the contour the same order.
		BOOST_TEST( postRate > orders[ o ] + 2.0 - pairSlack,
		            "k = " << orders[ o ] << ": rooting psi* converges at "
		            << postRate << " rather than at the k+2 that psi*'s own "
		            << "order gives it. If it has fallen back to k+1, the "
		            << "post-processed field is not the one being rooted" );

		BOOST_TEST( post.back() < 0.2*raw.back(),
		            "k = " << orders[ o ] << ": on the finest mesh rooting psi* "
		            << "gives " << post.back() << " against psi_h's "
		            << raw.back() << ". psi* IS THE DEFAULT and this is the "
		            << "measurement it is the default on" );
	}
	std::fflush( stdout );
}

/*
 * IN-1, THE CONTROL FIRST: THE METRIC IDENTITY ON THE CLOSED FORM.
 *
 * rho' = rho ( u . t ) / ( u' . t ) is derived in FluxSurfaces.hpp from
 * dx/dtheta being parallel to the tangent, and a derivation is not a
 * measurement. Two things are checked here and both are independent of the
 * solver:
 *
 *   * rho' from the identity against a SPECTRAL derivative of the sampled
 *     rho_j, which needs only the radii and knows nothing of the tangent;
 *   * the periodic trapezoid rule on sqrt( rho'^2 + rho^2 ) against the same
 *     rule at high resolution, which on an ANALYTIC contour must converge
 *     geometrically.
 *
 * It exists as the control for the discrete table below. If the discrete
 * first column floors, this is what says the floor is the FIELD -- psi_h's
 * jumps across faces make rho piecewise analytic rather than analytic -- and
 * not the rule.
 */
BOOST_AUTO_TEST_CASE( theMetricIdentityIsSpectralOnTheClosedForm )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );
	ExactAxis const axis = set.axis;
	double const level = set.levels[ 0 ];

	/// One ray of the exact contour: rho by Newton on psi( a + rho u ) = c with
	/// the analytic gradient, then rho' by the identity.
	auto ray = [ & ]( double theta, double seed, double &rho, double &rhoPrime )
	{
		double const cosine = std::cos( theta );
		double const sine = std::sin( theta );

		rho = seed;
		for ( int i = 0; i < 100; ++i )
		{
			double const r = axis.r + rho*cosine;
			double const z = axis.z + rho*sine;
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );
			double const slope = gr*cosine + gz*sine;
			double const residual = eq.psi( r, z ) - level;
			if ( std::abs( residual ) < 1.0e-15 )
				break;
			rho -= residual/slope;
		}

		double const r = axis.r + rho*cosine;
		double const z = axis.z + rho*sine;
		double gr = 0.0;
		double gz = 0.0;
		eq.gradPsi( r, z, gr, gz );
		double const magnitude = std::sqrt( gr*gr + gz*gz );
		double const tangentR = -gz/magnitude;
		double const tangentZ = gr/magnitude;

		double const uDotT = cosine*tangentR + sine*tangentZ;
		double const cross = cosine*tangentZ - sine*tangentR;
		rhoPrime = rho*uDotT/cross;
	};

	std::printf( "\n  The metric identity on the CLOSED FORM, Psi_N = %.2f\n",
	             set.fractions[ 0 ] );
	std::printf( "  %6s %16s %14s %14s %14s\n",
	             "N", "length", "|err| vs ref", "rho' vs FFT", "chord err" );

	// The reference: the same rule at high resolution. It is legitimate BECAUSE
	// the first column is the one being shown to converge geometrically -- an
	// algebraic column could not be its own reference.
	std::size_t const reference = 2048;
	double referenceLength = 0.0;
	{
		double seed = 0.1;
		double total = 0.0;
		for ( std::size_t j = 0; j < reference; ++j )
		{
			double rho = 0.0;
			double rhoPrime = 0.0;
			ray( twoPi*j/reference, seed, rho, rhoPrime );
			seed = rho;
			total += std::sqrt( rhoPrime*rhoPrime + rho*rho );
		}
		referenceLength = ( twoPi/reference )*total;
	}

	std::vector<std::size_t> const counts = { 16, 32, 64, 128 };
	std::vector<double> errors;

	for ( std::size_t c = 0; c < counts.size(); ++c )
	{
		std::size_t const n = counts[ c ];
		std::vector<double> rho( n, 0.0 );
		std::vector<double> rhoPrime( n, 0.0 );

		double seed = 0.1;
		for ( std::size_t j = 0; j < n; ++j )
		{
			ray( twoPi*j/n, seed, rho[ j ], rhoPrime[ j ] );
			seed = rho[ j ];
		}

		double total = 0.0;
		double chords = 0.0;
		for ( std::size_t j = 0; j < n; ++j )
		{
			total += std::sqrt( rhoPrime[ j ]*rhoPrime[ j ] + rho[ j ]*rho[ j ] );
			std::size_t const next = ( j + 1 )%n;
			double const theta = twoPi*j/n;
			double const nextTheta = twoPi*next/n;
			double const dr = rho[ next ]*std::cos( nextTheta ) - rho[ j ]*std::cos( theta );
			double const dz = rho[ next ]*std::sin( nextTheta ) - rho[ j ]*std::sin( theta );
			chords += std::sqrt( dr*dr + dz*dz );
		}
		double const length = ( twoPi/n )*total;

		// The spectral derivative of the SAMPLED radii -- an independent route
		// to rho' that uses no tangent at all. Trefethen's differentiation
		// matrix, cotangent form for even N.
		double worstDerivative = 0.0;
		for ( std::size_t j = 0; j < n; ++j )
		{
			double sum = 0.0;
			for ( std::size_t k = 0; k < n; ++k )
			{
				if ( k == j )
					continue;
				double const half = M_PI*( static_cast<double>( j )
				                           - static_cast<double>( k ) )/n;
				double const sign = ( ( j + k )%2 == 0 ) ? 1.0 : -1.0;
				sum += 0.5*sign/std::tan( half )*rho[ k ];
			}
			worstDerivative = std::max( worstDerivative, std::abs( sum - rhoPrime[ j ] ) );
		}

		double const error = std::abs( length - referenceLength );
		errors.push_back( error );
		std::printf( "  %6d %16.12f %14.6e %14.6e %14.6e\n",
		             static_cast<int>( n ), length, error, worstDerivative,
		             std::abs( chords - referenceLength ) );

		// AT THE FINEST N ONLY, and the coarser rows say why: the spectral
		// derivative of the sampled radii is itself converging -- 2.3e-2,
		// 9.3e-4, 3.8e-6, 6.4e-11 -- so at N = 16 the disagreement is the
		// CHECK's error and not the identity's. Asserting it at every N would be
		// asserting that the spectral derivative is already converged, which is
		// a different claim and a false one.
		if ( c + 1 == counts.size() )
			BOOST_TEST( worstDerivative < 1.0e-9,
			            "N = " << n << ": rho' from the identity disagrees with a "
			            << "spectral derivative of the same radii by "
			            << worstDerivative << ". The identity is derived in "
			            << "FluxSurfaces.hpp and this is the only check on the "
			            << "derivation that does not use the tangent" );
	}
	std::fflush( stdout );

	BOOST_TEST( errors.back() < 1.0e-11,
	            "on an ANALYTIC contour the trapezoid rule with rho' from the "
	            << "gradient reaches only " << errors.back()
	            << " at N = " << counts.back()
	            << ". If this is not at round-off, the rule or the identity is "
	            << "wrong and the discrete table below cannot be read" );
}

/*
 * IN-1'S ACCEPTANCE: THREE COLUMNS OF ONE TABLE, AND THE MIDDLE ONE IS THE
 * DELIVERABLE.
 *
 *   the answer     periodic trapezoid in theta, rho' pointwise from q
 *   THE TRAP       periodic trapezoid in theta, rho' by central difference
 *                  of the rho_j                                    O( N^-2 )
 *   the control    chord sum                                       O( N^-2 )
 *
 * The middle column is a spectrally accurate rule reduced to second order by
 * its Jacobian, with nothing in its own output saying so. IF EITHER CONTROL
 * EVER CONVERGES AS FAST AS THE FIRST THE COMPARISON IS EMPTY AND THIS TEST IS
 * WORTHLESS, and the failure messages say so.
 *
 * The first column is NOT spectral on a discrete contour and that is measured
 * rather than hidden: psi_h jumps across faces, so rho( theta ) is piecewise
 * analytic and the trapezoid rule cannot be geometric on it. It converges very
 * fast and then floors at the jump level.
 * theMetricIdentityIsSpectralOnTheClosedForm is the control that says so.
 */
BOOST_AUTO_TEST_CASE( theSpectralRuleIsOnlyAsGoodAsItsJacobian )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	SolvedEquilibrium solved( eq, box, 3, 32 );
	meq::CriticalPoint const axis = solved.axis();

	meq::ContourTracer tracer( solved.theSolver() );
	tracer.setTargetTurn( 0.05 );
	meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );
	BOOST_TEST( contour.closed(), "the contour to be re-parametrised did not close" );

	// AND rho' CHECKED A THIRD WAY ON THE DISCRETE CONTOUR, by the spectral
	// derivative of the sampled radii -- which uses no tangent, no flux and no
	// identity, only the rho_j. It is the same check
	// theMetricIdentityIsSpectralOnTheClosedForm makes on the analytic contour,
	// made here on the field the solver actually produced, and it is what
	// exercises AngleParametrisation::spectralRadiusPrime().
	{
		meq::AngleParametrisation const fit = tracer.fitByAngle( contour, axis, 256 );
		std::vector<double> const spectral = fit.spectralRadiusPrime();

		double worst = 0.0;
		double scale = 0.0;
		for ( std::size_t j = 0; j < spectral.size(); ++j )
		{
			worst = std::max( worst, std::abs( spectral[ j ] - fit.radiusPrime[ j ] ) );
			scale = std::max( scale, std::abs( fit.radiusPrime[ j ] ) );
		}

		std::printf( "  rho' from the identity against a spectral derivative of "
		             "the same radii: %.3e\n  on a scale of %.3e, that is %.2e "
		             "relative; fit fallbacks %d\n",
		             worst, scale, worst/scale, fit.fallbackLocations );

		BOOST_TEST( worst < 1.0e-4*scale,
		            "rho' from q and a spectral derivative of the radii disagree "
		            << "by " << worst << " on a scale of " << scale
		            << ". The two share the traced points and nothing else, so a "
		            << "disagreement is the identity of FluxSurfaces.hpp being "
		            << "wrong or the level set being far rougher than the DG "
		            << "jump accounts for" );
	}

	std::size_t const reference = 4096;
	meq::AngleParametrisation const referenceFit =
		tracer.fitByAngle( contour, axis, reference );
	double const referenceLength = referenceFit.length();

	std::printf( "\n  Arc length of { psi_h = c } three ways, k = 3, n = 32, "
	             "Psi_N = %.2f\n", set.fractions[ 0 ] );
	std::printf( "  Hermite arc length of the traced curve %.12f, "
	             "reference (N = %d, from q) %.12f\n",
	             contour.hermiteLength( 12 ), static_cast<int>( reference ),
	             referenceLength );
	double const jumpScale = contour.points.front().r
	                         *contour.points.front().fluxMagnitude;
	std::printf( "  transversality min |u x t| = %.4f, bisections %d, "
	             "worst |psi - c| %.2e\n",
	             referenceFit.transversality, referenceFit.bisections,
	             referenceFit.worstResidual );
	std::printf( "  the DG jump on this contour is %.2e in psi, %.2e as a "
	             "distance -- THAT is where\n  the first column floors, and "
	             "theMetricIdentityIsSpectralOnTheClosedForm reaches 4e-15\n  "
	             "on an analytic contour with the identical rule\n",
	             contour.worstFaceJump,
	             jumpScale > 0.0 ? contour.worstFaceJump/jumpScale : 0.0 );
	std::printf( "  %6s %14s %7s %14s %7s %14s %7s\n",
	             "N", "from q", "rate", "differenced", "rate", "chords", "rate" );

	std::vector<std::size_t> const counts = { 16, 32, 64, 128 };
	std::vector<double> fromFlux;
	std::vector<double> differenced;
	std::vector<double> chords;

	for ( std::size_t c = 0; c < counts.size(); ++c )
	{
		meq::AngleParametrisation const fit =
			tracer.fitByAngle( contour, axis, counts[ c ] );

		double const a = std::abs( fit.length() - referenceLength );
		double const b = std::abs( fit.differencedLength() - referenceLength );
		double const d = std::abs( fit.chordLength() - referenceLength );

		fromFlux.push_back( a );
		differenced.push_back( b );
		chords.push_back( d );

		if ( c == 0 )
			std::printf( "  %6d %14.6e %7s %14.6e %7s %14.6e %7s\n",
			             static_cast<int>( counts[ c ] ), a, "-", b, "-", d, "-" );
		else
			std::printf( "  %6d %14.6e %7.3f %14.6e %7.3f %14.6e %7.3f\n",
			             static_cast<int>( counts[ c ] ), a,
			             rate( fromFlux[ c - 1 ], a, 2.0 ), b,
			             rate( differenced[ c - 1 ], b, 2.0 ), d,
			             rate( chords[ c - 1 ], d, 2.0 ) );
	}
	std::fflush( stdout );

	double const ratio = static_cast<double>( counts.back() )/counts.front();
	double const fluxRate = rate( fromFlux.front(), fromFlux.back(), ratio );
	double const differencedRate = rate( differenced.front(), differenced.back(), ratio );
	double const chordRate = rate( chords.front(), chords.back(), ratio );

	BOOST_TEST( differencedRate < 3.0,
	            "the differenced Jacobian converges at " << differencedRate
	            << ". THAT COLUMN IS THE CONTROL: it is the trap of "
	            << "INVERSION-PLAN.md section 3.2, a spectrally accurate rule "
	            << "fed a second-order metric, and if it ever converges as fast "
	            << "as the pointwise one the comparison is empty and this test "
	            << "is worthless" );

	BOOST_TEST( chordRate < 3.0,
	            "the chord sum converges at " << chordRate
	            << ". Same objection: it is a control and it must stay second "
	            << "order or there is nothing to compare against" );

	BOOST_TEST( differencedRate > 2.0 - pairSlack,
	            "the differenced Jacobian converges at " << differencedRate
	            << " rather than the second order a central difference gives" );

	BOOST_TEST( fluxRate > differencedRate + 1.0,
	            "the pointwise metric from q converges at " << fluxRate
	            << " against the differenced one's " << differencedRate
	            << ". The whole content of IN-1 is that taking the Jacobian from "
	            << "q rather than from a difference is worth more than a "
	            << "constant" );

	// AT THE FINEST N, NOT THE COARSEST. At N = 16 neither column is in its
	// asymptotic regime and they sit within a factor of ten of each other; the
	// separation the table exists to show is what the two rates open up over the
	// sweep, and it is five orders by the end.
	BOOST_TEST( fromFlux.back() < 1.0e-4*differenced.back(),
	            "at N = " << counts.back() << " the pointwise metric is "
	            << fromFlux.back() << " against the differenced one's "
	            << differenced.back() << ", which is not the separation the "
	            << "table is supposed to demonstrate" );
}

/*
 * THE SAME METRIC A SECOND, CHEAP WAY -- TWO INDEPENDENT ROUTES TO ONE NUMBER.
 *
 * The fit's | dx/dtheta | integrates to the arc length of the surface; so does
 * the Hermite interpolant's own speed, integrated over the traced segments.
 * They share the field and the contour and NOTHING else: one is a ray
 * parametrisation with a metric from the identity, the other is a cubic through
 * the traced points. Agreement is evidence for both; disagreement would not say
 * which.
 */
BOOST_AUTO_TEST_CASE( theMetricAndTheHermiteArcLengthAgree )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.05, 0.4, 0.8 } );

	SolvedEquilibrium solved( eq, box, 3, 32 );
	meq::CriticalPoint const axis = solved.axis();

	std::printf( "\n  Two routes to the arc length of { psi_h = c }, k = 3, n = 32\n" );
	std::printf( "  %6s %8s %16s %16s %12s %10s %9s\n",
	             "Psi_N", "points", "Hermite", "metric (N=1024)", "relative",
	             "min|u x t|", "bisect" );

	for ( std::size_t i = 0; i < set.levels.size(); ++i )
	{
		meq::ContourTracer tracer( solved.theSolver() );
		tracer.setTargetTurn( 0.05 );
		meq::Contour const contour = tracer.traceFromAxis( set.levels[ i ], axis );
		BOOST_TEST( contour.closed(),
		            "Psi_N = " << set.fractions[ i ] << ": the contour did not close" );

		meq::AngleParametrisation const fit = tracer.fitByAngle( contour, axis, 1024 );

		double const hermite = contour.hermiteLength( 12 );
		double const metric = fit.length();
		double const relative = std::abs( hermite - metric )/metric;

		std::printf( "  %6.2f %8d %16.10f %16.10f %12.3e %10.4f %9d\n",
		             set.fractions[ i ], static_cast<int>( contour.points.size() ),
		             hermite, metric, relative, fit.transversality, fit.bisections );

		BOOST_TEST( fit.transverse,
		            "Psi_N = " << set.fractions[ i ]
		            << ": the fit reports min | u x t | = " << fit.transversality );

		BOOST_TEST( relative < 1.0e-6,
		            "Psi_N = " << set.fractions[ i ] << ": the Hermite arc length "
		            << hermite << " and the metric's " << metric << " disagree by "
		            << relative << " relative. They share the field and the "
		            << "contour and nothing else, so a disagreement is a defect "
		            << "in one of the two and this test does not say which" );

		// AND THE INTERPOLANT IS A CURVE AND NOT A TABLE, which is what
		// Contour::pointAtArcLength() is for and what a consumer of this will
		// actually call. Sixty-four points at equal spacing along it, each
		// measured against the level set it is supposed to lie on -- the same
		// quantity representationError() measures at midpoints, asked of an
		// arbitrary parameter instead.
		double worstOff = 0.0;
		int hint = contour.points.front().element;
		for ( int j = 0; j < 64; ++j )
		{
			double r = 0.0;
			double z = 0.0;
			double psi = 0.0;
			double qR = 0.0;
			double qZ = 0.0;
			contour.pointAtArcLength( j*contour.length()/64.0, r, z );
			if ( !tracer.sampleAt( r, z, psi, qR, qZ, hint ) )
				continue;
			double const magnitude = r*std::sqrt( qR*qR + qZ*qZ );
			if ( magnitude > 0.0 )
				worstOff = std::max( worstOff,
				                     std::abs( psi - contour.level )/magnitude );
		}

		std::printf( "  %6s %8s %16s %16s   worst off-curve at 64 equal arc "
		             "lengths %.3e\n", "", "", "", "", worstOff );

		BOOST_TEST( worstOff < 1.0e-5,
		            "Psi_N = " << set.fractions[ i ] << ": a point taken from the "
		            << "interpolant at an arbitrary arc length is " << worstOff
		            << " off the level set, which is far worse than the segment "
		            << "midpoints are. pointAtArcLength() is what a consumer "
		            << "calls, so it has to be as good as the curve it indexes" );
	}
	std::fflush( stdout );
}

/*
 * STAR-SHAPEDNESS IS A HYPOTHESIS, AND THE FIT REFUSES RATHER THAN RETURNING A
 * NUMBER WHEN IT FAILS.
 *
 * INVERSION-PLAN.md section 3.4 keeps ray bisection as a cross-check and not as
 * the primary route precisely because it assumes star-shapedness, which fails
 * on indented cross-sections and near the axis where the bracket degenerates.
 * The denominator of the metric identity, u x t, vanishes exactly when the ray
 * is TANGENT to the curve -- so the hypothesis is not a qualitative worry, it is
 * a number, and IndexAudit::transversality is the model for reporting it.
 *
 * The demonstration is a fit about a point that is NOT the axis: far enough off
 * centre and some ray grazes the surface. It must throw, and the message must
 * say which hypothesis failed.
 */
BOOST_AUTO_TEST_CASE( theFitRefusesWhenTheRaysStopBeingTransverse )
{
	Equilibrium const eq = Equilibrium::iterExample2();
	Rectangle const box = benchmarkBox();
	LevelSet const set = levelsFor( eq, box, { 0.4 } );

	SolvedEquilibrium solved( eq, box, 2, 24 );
	meq::CriticalPoint const axis = solved.axis();

	meq::ContourTracer tracer( solved.theSolver() );
	tracer.setTargetTurn( 0.05 );
	meq::Contour const contour = tracer.traceFromAxis( set.levels[ 0 ], axis );
	BOOST_TEST( contour.closed(), "the contour did not close" );

	// How far off centre the surface is, so that the displaced origin below can
	// be reported as a fraction of it rather than as a bare length.
	double smallest = std::numeric_limits<double>::infinity();
	double largest = 0.0;
	for ( std::size_t i = 0; i < contour.points.size(); ++i )
	{
		double const dr = contour.points[ i ].r - axis.r;
		double const dz = contour.points[ i ].z - axis.z;
		double const radius = std::sqrt( dr*dr + dz*dz );
		smallest = std::min( smallest, radius );
		largest = std::max( largest, radius );
	}

	std::printf( "\n  Star-shapedness as a measurement, Psi_N = %.2f\n",
	             set.fractions[ 0 ] );
	std::printf( "  the surface spans rho in [ %.5f, %.5f ] about the axis\n",
	             smallest, largest );
	std::printf( "  %14s %14s %10s %s\n",
	             "origin offset", "min |u x t|", "bisect", "outcome" );

	meq::AngleParametrisation const good = tracer.fitByAngle( contour, axis, 128 );
	std::printf( "  %14.5f %14.4f %10d %s\n", 0.0, good.transversality,
	             good.bisections, "accepted" );

	BOOST_TEST( good.transversality > 0.5,
	            "about the axis itself the rays should be nearly normal to the "
	            << "surface and min | u x t | reads " << good.transversality );

	// Displaced well outside the surface: every ray from there either misses it
	// or grazes it, so the fit must refuse.
	meq::CriticalPoint offCentre = axis;
	offCentre.r = axis.r + 3.0*largest;

	bool refused = false;
	try
	{
		meq::AngleParametrisation const bad = tracer.fitByAngle( contour, offCentre, 128 );
		std::printf( "  %14.5f %14.4f %10d %s\n", 3.0*largest, bad.transversality,
		             bad.bisections, "ACCEPTED -- it should not have been" );
	}
	catch ( std::exception const &error )
	{
		refused = true;
		std::printf( "  %14.5f %14s %10s refused: %.90s\n", 3.0*largest, "-", "-",
		             error.what() );
	}
	std::fflush( stdout );

	BOOST_TEST( refused,
	            "a fit about a point outside the surface was accepted. Ray "
	            << "methods assume star-shapedness and this one must report the "
	            << "hypothesis rather than survive its failure -- a number "
	            << "returned from a degenerate fit is worse than no number" );
}

/*
 * IN-0's SECOND ACCEPTANCE, PART ONE: A CONTOUR THAT CROSSES THE BAND COMPLETES,
 * AND SAYS SO PER POINT.
 *
 * The prior art is what this is written against. v0-legacy:FluxSurfaces.cpp
 * printed "Terminating because curve left domain" and returned a partial curve,
 * which for an outer surface is an ARC LABELLED AS A CLOSED CONTOUR -- and a
 * flux-surface average over it is wrong by an amount nothing reports. So the
 * first thing asserted is the control: with no extension configured the trace
 * does NOT close and returns ContourStatus::LeftMesh. If that ever closes on
 * its own the rest of this file is measuring nothing.
 *
 * AND THE FLAG IS A MASK AND NOT A COUNT. CLAUDE.md records extrapolated_nodes
 * having been a count rather than a mask as the other half of a real defect in
 * the .nc: a band node holds real data, so the `inside` mask says 1 there, and
 * nothing downstream could tell WHICH nodes were continued. The same is true of
 * a contour, so ContourPoint::extended is checked point by point against an
 * INDEPENDENT statement of the same fact -- mfem::Mesh::FindPoints, which knows
 * nothing of the band machinery.
 */
BOOST_AUTO_TEST_CASE( theBandIsCrossedOnlyWhenAnExtensionSaysHow )
{
	int const order = 2;
	int const n = 24;
	double const fraction = 0.90;
	double const level = curvedLevel( fraction );

	CurvedSolve solved( order, n );
	meq::CriticalPoint const axis = solved.axis();

	std::printf( "\n  The band, k = %d, n = %d, h = %.5f, Psi_N = %.2f, "
	             "psi = %.6e\n", order, n, solved.h(), fraction, level );
	std::printf( "  %14s %9s %8s %8s %8s %9s %9s\n", "extension", "status",
	             "points", "band", "stalled", "deepest", "deep/h" );

	meq::BandExtension const methods[] = { meq::BandExtension::None,
	                                       meq::BandExtension::FluxTaylor,
	                                       meq::BandExtension::TransferLift };

	for ( meq::BandExtension method : methods )
	{
		meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
		if ( method != meq::BandExtension::None )
			useBandExtension( tracer, solved, method );

		BOOST_TEST( ( tracer.bandExtension() == method ),
		            "setBandExtension() did not take" );

		meq::Contour const contour = tracer.traceFromAxis( level, axis );

		std::printf( "  %14s %9s %8zu %8d %8d %9.3e %9.3f\n",
		             meq::bandExtensionName( method ),
		             meq::contourStatusName( contour.status ),
		             contour.points.size(), contour.extendedPoints,
		             contour.stalledCorrections, contour.deepestBandPoint,
		             contour.deepestBandPoint/solved.h() );
		std::fflush( stdout );

		if ( method == meq::BandExtension::None )
		{
			// THE CONTROL. A contour at Psi_N = 0.9 on an inscribed Gamma_h has
			// to leave the mesh, and without an extension the tracer must say so
			// rather than return the arc.
			BOOST_TEST( !contour.closed(),
			            "with no band extension the contour at Psi_N = " << fraction
			            << " closed anyway, so it never met Gamma_h and this "
			            "whole comparison is empty -- pick a level nearer Gamma" );
			BOOST_TEST( ( contour.status == meq::ContourStatus::LeftMesh ),
			            "with no band extension the trace ended as '"
			            << meq::contourStatusName( contour.status )
			            << "' rather than 'left mesh'" );
			BOOST_TEST( contour.extendedPoints == 0,
			            "a trace with no extension reported "
			            << contour.extendedPoints << " band points" );
			BOOST_TEST( !contour.crossesBand(),
			            "crossesBand() is true with no extension configured" );
			continue;
		}

		BOOST_TEST( contour.closed(),
		            meq::bandExtensionName( method ) << ": the contour ended as '"
		            << meq::contourStatusName( contour.status )
		            << "' rather than closing across the band" );
		BOOST_TEST( contour.extendedPoints > 0,
		            meq::bandExtensionName( method )
		            << ": the contour closed without a single point in the band, "
		            "so it did not cross it" );
		BOOST_TEST( contour.crossesBand(),
		            meq::bandExtensionName( method )
		            << ": extendedPoints is positive but crossesBand() is false" );
		BOOST_TEST( ( contour.bandExtension == method ),
		            "the contour does not record which extension answered for it" );

		// THE MASK, POINT BY POINT, AGAINST AN INDEPENDENT LOCATOR.
		int flaggedButInside = 0;
		int outsideButUnflagged = 0;
		int counted = 0;
		std::size_t const last = contour.points.size() - 1;

		for ( std::size_t i = 0; i < last; ++i )
		{
			meq::ContourPoint const &p = contour.points[ i ];

			mfem::DenseMatrix matrix( 2, 1 );
			matrix( 0, 0 ) = p.r;
			matrix( 1, 0 ) = p.z;
			mfem::Array<int> found;
			mfem::Array<mfem::IntegrationPoint> ips;
			solved.mesh().FindPoints( matrix, found, ips, false );

			bool const inside = found.Size() > 0 && found[ 0 ] >= 0;
			if ( p.extended && inside )
				++flaggedButInside;
			if ( !p.extended && !inside )
				++outsideButUnflagged;
			if ( p.extended )
			{
				++counted;
				BOOST_TEST( p.bandDepth > 0.0,
				            "a band point reports a depth of " << p.bandDepth );
			}
			else
			{
				BOOST_TEST( p.bandDepth == 0.0,
				            "an interior point reports a band depth of "
				            << p.bandDepth );
			}
		}

		BOOST_TEST( flaggedButInside == 0,
		            meq::bandExtensionName( method ) << ": " << flaggedButInside
		            << " points are flagged as extended while Mesh::FindPoints "
		            "puts them inside the mesh" );
		BOOST_TEST( outsideButUnflagged == 0,
		            meq::bandExtensionName( method ) << ": " << outsideButUnflagged
		            << " points are outside the mesh and NOT flagged, which is "
		            "the defect the mask exists to prevent -- a consumer cannot "
		            "drop what it is not told about" );
		BOOST_TEST( counted == contour.extendedPoints,
		            "the per-point mask counts " << counted << " band points and "
		            "Contour::extendedPoints says " << contour.extendedPoints
		            << "; a count that disagrees with its own mask is worse than "
		            "no count" );
	}

	// AND THE REACH IS A HARD LIMIT THAT SAYS SO WHEN IT BITES. An extension is
	// only analysed within O( h ) of Gamma_h, so setBandExtension() takes a
	// reach and refuses beyond it -- and the refusal has to surface as a STATUS
	// rather than as a short curve, which is the whole complaint against the
	// prior art. Swept on this benchmark: Gamma's deepest point sits at 0.92 to
	// 0.98 of a face length outside Gamma_h, so a reach of one works with two
	// per cent to spare and 0.75 truncates at every mesh. Half is well inside
	// that and is asserted here because a limit that silently did nothing would
	// pass every other test in this file.
	{
		meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
		tracer.setBandExtension( meq::BandExtension::TransferLift,
		                         solved.gammaHMarker(), &solved.paths(),
		                         mfem::PositionFunction(), 0.5 );

		meq::Contour const contour = tracer.traceFromAxis( level, axis );
		std::printf( "  %14s %9s %8zu %8d\n", "reach = 0.5",
		             meq::contourStatusName( contour.status ),
		             contour.points.size(), contour.extendedPoints );
		std::fflush( stdout );

		BOOST_TEST( !contour.closed(),
		            "a reach of half a face length let the contour cross a band "
		            "that is nearly a whole face deep, so the limit is not being "
		            "applied" );
		BOOST_TEST( ( contour.status == meq::ContourStatus::LeftMesh ),
		            "a trace stopped by the reach ended as '"
		            << meq::contourStatusName( contour.status )
		            << "' rather than 'left mesh', so a consumer is not told why "
		            "it is holding an arc" );
	}
}

/*
 * IN-0's SECOND ACCEPTANCE, PART TWO, AND THE DELIVERABLE OF THIS HALF OF THE
 * STAGE: THE TWO BAND EXTENSIONS MEASURED SIDE BY SIDE ON THE SAME CONTOURS.
 *
 * INVERSION-PLAN.md section 4.3 lists three candidates and says the choice is
 * to be MEASURED, not argued. Two of them are fields and are measured here; the
 * third -- psi = 0 exactly on Gamma, which is known analytically -- is not a
 * field at all but an ANCHOR, and it is what the lift's line integral starts
 * from. It is deliberately not used as a correction: CLAUDE.md records that
 * blending an extrapolated value towards a known zero does not rescue a bad
 * extension, because ( 1 - t ) v scales a positive value down and never changes
 * its sign. The error was in WHERE the field was evaluated.
 *
 * WHAT SEPARATES THE TWO, AND IT IS ARITHMETIC RATHER THAN TASTE.
 *
 *   FluxTaylor    psi( x0 ) + r0 q( x0 ) . ( p - x0 ). Its error is the error
 *                 of psi_h at the foot, O( h^(k+1) ), PLUS the Taylor remainder
 *                 of a second-order expansion carried over a band of width
 *                 O( h ), which is O( h^2 ) whatever k is. So it is exact in
 *                 order at k = 1 and CAPS THE BAND AT SECOND ORDER from k = 2
 *                 on. It is the control.
 *
 *   TransferLift  g( a( x0 ) ) + the line integral of -r q back from Gamma,
 *                 with q outside the mesh supplied by the method's own
 *                 extension operator. Its error is the error of q integrated
 *                 over a path of length O( h ) -- so O( h^(k+2) ), a full order
 *                 BETTER than psi_h's own, and the band should not limit the
 *                 contour at all.
 *
 * THE TWO POPULATIONS ARE KEPT APART, per section 4.3. The interior column is
 * the same contour's points inside Omega_h and converges at psi_h's k+1
 * whichever extension is configured -- it is the control on the control, since
 * an extension that damaged the interior would show up there and nowhere else.
 *
 * AND THE CONTROL MUST LOSE, OR THE COMPARISON IS EMPTY. That is the shape
 * ExtensionConvergence.cpp uses for its pinned-zero column and the reason is
 * the same: a measurement against something that also works measures nothing.
 */
BOOST_AUTO_TEST_CASE( theTwoBandExtensionsAreMeasuredSideBySide )
{
	std::vector<int> const orders = { 1, 2, 3 };

	// NOT { 12, 24, 48 }, AND THE COARSEST POINT IS WHY -- the same objection
	// theTracedContourConvergesAtTheFieldsOwnOrder records against its own
	// sweep. At k = 1 on h = 0.142 a contour at Psi_N = 0.9 is not resolved at
	// all: the interior error reads 8.1e-2 on a contour about 2.5 across, the
	// Taylor extension's own O( h^2 ) is then 2e-2 against a level of 2.4e-2,
	// and the trace WANDERS IN THE BAND -- 187,000 band points and
	// ContourStatus::TooLong. That is the extension being useless at that
	// resolution and being honest about it, not a defect, but a sweep that
	// starts there measures the approach to the asymptotic regime rather than
	// the rate.
	std::vector<int> const meshes = { 16, 32, 64 };

	double const fraction = 0.90;
	double const level = curvedLevel( fraction );
	double const deltaS = 0.02;

	// A contour of this surface is about 5 long, so 250 points at Delta_s. The
	// cap is eighty times that, and it is here so that a trace which DOES wander
	// costs a second rather than a minute -- finding out cheaply that something
	// does not work is worth having.
	int const pointCap = 20000;

	meq::BandExtension const methods[] = { meq::BandExtension::FluxTaylor,
	                                       meq::BandExtension::TransferLift };

	for ( std::size_t o = 0; o < orders.size(); ++o )
	{
		int const order = orders[ o ];

		std::vector<std::vector<BandMeasurement>> measured( 2 );
		std::vector<double> meshWidth;

		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			CurvedSolve solved( order, meshes[ m ] );
			meq::CriticalPoint const axis = solved.axis();
			meshWidth.push_back( solved.h() );

			for ( int which = 0; which < 2; ++which )
			{
				meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
				useBandExtension( tracer, solved, methods[ which ] );
				useFixedStep( tracer, deltaS );
				tracer.setMaxPoints( pointCap );

				meq::Contour const contour = tracer.traceFromAxis( level, axis );

				BOOST_TEST( contour.closed(),
				            meq::bandExtensionName( methods[ which ] ) << ", k = "
				            << order << ", n = " << meshes[ m ]
				            << ": the contour ended as '"
				            << meq::contourStatusName( contour.status ) << "'" );
				BOOST_TEST( contour.crossesBand(),
				            meq::bandExtensionName( methods[ which ] ) << ", k = "
				            << order << ", n = " << meshes[ m ]
				            << ": the contour did not reach the band, so there is "
				            "nothing here to compare" );

				measured[ which ].push_back( splitByBand( contour ) );
			}
		}

		std::printf( "\n  Band and interior error against the exact contour, "
		             "Psi_N = %.2f, psi_h ( raw ), k = %d, Delta_s = %.3f\n",
		             fraction, order, deltaS );
		std::printf( "  %14s %5s %8s %7s %6s %8s %13s %6s %13s %6s\n",
		             "extension", "n", "h", "points", "band", "deep/h",
		             "band error", "rate", "interior", "rate" );

		for ( int which = 0; which < 2; ++which )
		{
			for ( std::size_t m = 0; m < meshes.size(); ++m )
			{
				BandMeasurement const &p = measured[ which ][ m ];
				std::printf( "  %14s %5d %8.5f %7d %6d %8.3f %13.6e",
				             meq::bandExtensionName( methods[ which ] ),
				             meshes[ m ], meshWidth[ m ],
				             p.bandPoints + p.interiorPoints, p.bandPoints,
				             p.deepest/meshWidth[ m ], p.bandWorst );
				if ( m == 0 )
					std::printf( " %6s %13.6e %6s\n", "-", p.interiorWorst, "-" );
				else
					std::printf( " %6.3f %13.6e %6.3f\n",
					             rate( measured[ which ][ m - 1 ].bandWorst,
					                   p.bandWorst, 2.0 ),
					             p.interiorWorst,
					             rate( measured[ which ][ m - 1 ].interiorWorst,
					                   p.interiorWorst, 2.0 ) );
			}
		}

		double const span = meshWidth.front()/meshWidth.back();
		double bandRate[ 2 ];
		double interiorRate[ 2 ];
		for ( int which = 0; which < 2; ++which )
		{
			bandRate[ which ] = rate( measured[ which ].front().bandWorst,
			                          measured[ which ].back().bandWorst, span );
			interiorRate[ which ] = rate( measured[ which ].front().interiorWorst,
			                              measured[ which ].back().interiorWorst,
			                              span );
		}

		double const ratio = measured[ 0 ].back().bandWorst
		                     /measured[ 1 ].back().bandWorst;

		std::printf( "    across the sequence, k = %d:  band  flux Taylor %.3f, "
		             "transfer lift %.3f\n", order, bandRate[ 0 ], bandRate[ 1 ] );
		std::printf( "    %31s interior  flux Taylor %.3f, transfer lift %.3f\n",
		             "", interiorRate[ 0 ], interiorRate[ 1 ] );
		std::printf( "    on the finest mesh the lift is %.1f times closer to the "
		             "exact contour in the band\n", ratio );
		std::fflush( stdout );

		// THE ANSWER. The lift's error is the error of q integrated over a path
		// of length O( h ), so it is a full order better than psi_h's own k+1 and
		// the band does not limit the contour at all. Asserted at k+1 rather than
		// at k+2, which is what it measures: the sequence slack here has to
		// absorb the shape of D_h changing with h, which is not a smooth function
		// of h -- ExtensionConvergence.cpp's reason for its own two-tier
		// assertion.
		BOOST_TEST( bandRate[ 1 ] > order + 1.0 - sequenceSlack,
		            "k = " << order << ": the transfer lift converges in the band "
		            "at " << bandRate[ 1 ] << " rather than at the k+1 = "
		            << order + 1 << " that the flux's own order gives it" );

		// THE CONTROL MUST LOSE, or the comparison is empty and this test is
		// worthless -- the shape ExtensionConvergence.cpp uses for its
		// pinned-zero column.
		//
		// ASSERTED ON THE MARGIN AND NOT ON THE RATE, AND THE deep/h COLUMN IS
		// WHY. A contour at a fixed Psi_N sits a fixed distance inside Gamma
		// while Gamma_h climbs towards Gamma, so its excursion into the band
		// shrinks faster than h -- 136 band points become 9 across this sweep --
		// and BOTH columns then converge faster than the extension they are
		// built on. Neither rate here is the extension's own order and reading
		// them as one would be wrong.
		// theOutermostSurfaceIsTracedEntirelyInTheBand is where the orders are
		// measured, on a curve that stays in the band at every mesh.
		if ( order >= 2 )
			BOOST_TEST( ratio > 100.0,
			            "k = " << order << ": the transfer lift is only " << ratio
			            << " times closer than the flux Taylor step on the finest "
			            "mesh. The Taylor step's remainder is second order in the "
			            "band width whatever k is and the lift's is the flux's "
			            "own, so if these two are close the control has stopped "
			            "being a control and the comparison is empty" );

		// AND NEITHER EXTENSION MAY DAMAGE THE INTERIOR. The same contour's
		// points inside Omega_h are the discretisation's own and must converge at
		// k+1 whichever extension answered for the rest of the curve; an
		// extension that perturbed them would show up here and nowhere else.
		for ( int which = 0; which < 2; ++which )
			BOOST_TEST( interiorRate[ which ] > order + 1.0 - sequenceSlack,
			            meq::bandExtensionName( methods[ which ] ) << ", k = "
			            << order << ": the points INSIDE the mesh converge at "
			            << interiorRate[ which ] << " rather than at k+1 = "
			            << order + 1 << ", so the extension has reached back into "
			            "the discretisation" );
	}
}

/*
 * IN-0's SECOND ACCEPTANCE, PART THREE: THE OUTERMOST SURFACE, WHICH IS Gamma
 * ITSELF AND LIES ENTIRELY IN THE BAND.
 *
 * WHY THIS MEASUREMENT EXISTS BESIDE THE ONE ABOVE, AND IT IS NOT A REFINEMENT
 * OF IT. Read the deep/h column of that table: 0.92, 0.66, 0.26. A contour at a
 * FIXED Psi_N sits a fixed distance inside Gamma, and Gamma_h climbs towards
 * Gamma as h falls, so the excursion the curve makes into the band SHRINKS
 * FASTER THAN h -- 136 band points become 9. Both columns there therefore
 * converge faster than the extension's own order, and reading either as that
 * order would be wrong. It is a realistic measurement of what the band costs a
 * given surface and it is not a measurement of the extension.
 *
 * The level psi = 0 is. It is Gamma, it is attained nowhere inside Omega_h --
 * psi_h is strictly negative there -- so EVERY point of it is answered by the
 * extension at every mesh, and the exact answer is known in closed form because
 * Gamma is the level set the subdomain was cut from. Nothing shrinks and there
 * is nothing to disentangle.
 *
 * AND IT IS THE THIRD CANDIDATE OF INVERSION-PLAN.md SECTION 4.3 MADE INTO A
 * NUMBER. That candidate is "psi = 0 exactly on Gamma, and Gamma is known
 * analytically, so the outermost contour is known for free with no tracing at
 * all". It is an ANCHOR and not a field: TransferLift starts its line integral
 * from g( a( x0 ) ) and so carries it, FluxTaylor knows nothing about Gamma and
 * so does not. This table is what the anchor is worth, and the lift is
 * SUPPOSED to win it -- that is the design and not an artefact. What the table
 * adds is how much, and that the Taylor step converges at all.
 */
BOOST_AUTO_TEST_CASE( theOutermostSurfaceIsTracedEntirelyInTheBand )
{
	std::vector<int> const orders = { 1, 2, 3 };
	std::vector<int> const meshes = { 16, 32, 64 };
	double const deltaS = 0.02;
	int const pointCap = 20000;

	meq::BandExtension const methods[] = { meq::BandExtension::FluxTaylor,
	                                       meq::BandExtension::TransferLift };

	for ( std::size_t o = 0; o < orders.size(); ++o )
	{
		int const order = orders[ o ];

		std::vector<std::vector<BandMeasurement>> measured( 2 );
		std::vector<double> meshWidth;

		for ( std::size_t m = 0; m < meshes.size(); ++m )
		{
			CurvedSolve solved( order, meshes[ m ] );
			meq::CriticalPoint const axis = solved.axis();
			meshWidth.push_back( solved.h() );

			for ( int which = 0; which < 2; ++which )
			{
				meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
				useBandExtension( tracer, solved, methods[ which ] );
				useFixedStep( tracer, deltaS );
				tracer.setMaxPoints( pointCap );

				meq::Contour const contour = tracer.traceFromAxis( 0.0, axis );

				BOOST_TEST( contour.closed(),
				            meq::bandExtensionName( methods[ which ] ) << ", k = "
				            << order << ", n = " << meshes[ m ]
				            << ": Gamma itself ended as '"
				            << meq::contourStatusName( contour.status ) << "'" );

				BandMeasurement const split = splitByBand( contour );

				// ESSENTIALLY ALL OF IT, AND "ESSENTIALLY" IS EARNED RATHER THAN
				// HEDGED. psi_h is negative throughout Omega_h to the accuracy of
				// the solve, so { psi_h = 0 } is outside it -- but only to that
				// accuracy: near Gamma_h the imposed datum is O( h ) below zero,
				// and a solution that creeps a little above it there puts a stray
				// point of the traced curve back inside. Measured, that is ONE
				// point of 322 at k = 1, n = 64 with the Taylor step and none
				// anywhere else. It is the discretisation and not the tracer, and
				// the property this test needs is that the curve is a band
				// measurement rather than a mixture.
				int const total = split.bandPoints + split.interiorPoints;
				BOOST_TEST( split.bandPoints > 0.98*total,
				            meq::bandExtensionName( methods[ which ] ) << ", k = "
				            << order << ", n = " << meshes[ m ] << ": only "
				            << split.bandPoints << " of " << total
				            << " points of Gamma are in the band, so this is not "
				            "the pure band measurement it is written to be" );

				measured[ which ].push_back( split );
			}
		}

		std::printf( "\n  Gamma itself, traced wholly in the band: distance from "
		             "the exact Gamma, k = %d, Delta_s = %.3f\n", order, deltaS );
		std::printf( "  %14s %5s %8s %7s %8s %13s %6s %13s\n", "extension", "n",
		             "h", "points", "deep/h", "worst", "rate", "rms" );

		for ( int which = 0; which < 2; ++which )
			for ( std::size_t m = 0; m < meshes.size(); ++m )
			{
				BandMeasurement const &p = measured[ which ][ m ];
				std::printf( "  %14s %5d %8.5f %7d %8.3f %13.6e",
				             meq::bandExtensionName( methods[ which ] ),
				             meshes[ m ], meshWidth[ m ], p.bandPoints,
				             p.deepest/meshWidth[ m ], p.bandWorst );
				if ( m == 0 )
					std::printf( " %6s %13.6e\n", "-", p.bandRms );
				else
					std::printf( " %6.3f %13.6e\n",
					             rate( measured[ which ][ m - 1 ].bandWorst,
					                   p.bandWorst, 2.0 ), p.bandRms );
			}

		double const span = meshWidth.front()/meshWidth.back();
		double const taylorRate = rate( measured[ 0 ].front().bandWorst,
		                                measured[ 0 ].back().bandWorst, span );
		double const liftRate = rate( measured[ 1 ].front().bandWorst,
		                              measured[ 1 ].back().bandWorst, span );
		double const ratio = measured[ 0 ].back().bandWorst
		                     /measured[ 1 ].back().bandWorst;

		std::printf( "    across the sequence, k = %d: flux Taylor %.3f, "
		             "transfer lift %.3f, and the lift is %.1f times closer on "
		             "the finest mesh\n", order, taylorRate, liftRate, ratio );
		std::fflush( stdout );

		// THE BAND MUST NOT LIMIT THE CONTOUR. The lift's error is the error of
		// q integrated over a path of length O( h ), so it is at worst psi_h's
		// own k+1 and measures nearer k+2 from k = 2 on -- 3.995 and 4.848 here.
		// Asserted at k+1 because that is the property IN-2 needs: a surface in
		// the band converges at the rate a surface inside Omega_h does.
		BOOST_TEST( liftRate > order + 1.0 - sequenceSlack,
		            "k = " << order << ": the transfer lift converges on Gamma at "
		            << liftRate << " rather than at the k+1 = " << order + 1
		            << " that the flux's own order gives it, so the band IS "
		            "limiting the outermost surface" );

		// AND THE CONTROL MUST LOSE. Here it can be asserted on the RATE, which
		// it cannot be on a contour at fixed Psi_N: this curve is wholly in the
		// band at every mesh, so nothing shrinks and 2.138 is the Taylor step's
		// own second order rather than an artefact of a vanishing excursion. If
		// it ever reaches k+1 the two extensions are the same thing and the
		// comparison is empty.
		BOOST_TEST( taylorRate > 1.5,
		            "k = " << order << ": the flux Taylor step converges on Gamma "
		            "at " << taylorRate << ", which is not the second order its "
		            "remainder gives it -- the control is broken rather than "
		            "merely losing, and a comparison against something broken "
		            "measures nothing" );

		if ( order >= 2 )
		{
			BOOST_TEST( taylorRate < order + 1.0 - sequenceSlack,
			            "k = " << order << ": the flux Taylor step converges on "
			            "Gamma at " << taylorRate << ", which is the k+1 = "
			            << order + 1 << " it cannot have. Its remainder is second "
			            "order in the band width at every k, so if this is real "
			            "the control has stopped being a control" );
			BOOST_TEST( liftRate > taylorRate + 1.0,
			            "k = " << order << ": the transfer lift converges at "
			            << liftRate << " against the Taylor step's " << taylorRate
			            << ", less than the full order that separates an "
			            "O( h^2 ) remainder from an integrated O( h^(k+1) ) flux" );
			BOOST_TEST( ratio > 100.0,
			            "k = " << order << ": the transfer lift is only " << ratio
			            << " times closer to Gamma than the Taylor step on the "
			            "finest mesh" );
		}
	}
}

/*
 * IN-0's SECOND ACCEPTANCE, PART FOUR: THE FLAG SURVIVES THE REPARAMETRISATION,
 * BECAUSE IN-2 IS WHAT NEEDS IT.
 *
 * INVERSION-PLAN.md section 4.3 asks that "every flux-surface quantity computed
 * on a surface that crosses the band be flagged as such", so that IN-2 reports
 * the two populations separately rather than averaging them in. A flux-surface
 * average is not taken over the traced points: it is taken over a
 * REPARAMETRISATION of them, at equispaced poloidal angles found by their own
 * ray Newton. Those rays go through the same seam and land where they land, so
 * the band does not respect the angular grid and a surface can be inside
 * Omega_h at one theta and outside it at the next.
 *
 * So AngleParametrisation carries the same pair the Contour does -- a per-node
 * mask and a count -- and the mask is checked here against
 * mfem::Mesh::FindPoints, which knows nothing of the band machinery. A count
 * alone would say a quantity is affected and not WHERE, which is the .nc's own
 * defect and the reason section 4.3 says so explicitly.
 */
BOOST_AUTO_TEST_CASE( theBandFlagReachesTheAngleParametrisation )
{
	int const order = 2;
	int const n = 24;
	std::size_t const count = 128;
	double const fraction = 0.90;
	double const level = curvedLevel( fraction );

	CurvedSolve solved( order, n );
	meq::CriticalPoint const axis = solved.axis();

	meq::ContourTracer tracer( solved.theSolver(), meq::Potential::Raw );
	useBandExtension( tracer, solved, meq::BandExtension::TransferLift );

	meq::Contour const contour = tracer.traceFromAxis( level, axis );
	BOOST_TEST_REQUIRE( contour.closed(),
	                    "the contour ended as '"
	                    << meq::contourStatusName( contour.status ) << "'" );
	BOOST_TEST_REQUIRE( contour.crossesBand(),
	                    "the contour does not cross the band, so there is no flag "
	                    "to carry" );

	meq::AngleParametrisation const fit =
		tracer.fitByAngle( contour, axis, count );

	BOOST_TEST( fit.extended.size() == count,
	            "the fit carries " << fit.extended.size() << " flags for "
	            << count << " nodes" );
	BOOST_TEST( fit.extendedNodes > 0,
	            "the contour crosses the band with " << contour.extendedPoints
	            << " points and the fit reports no band node at all, so the flag "
	            "is lost exactly where IN-2 would read it" );
	BOOST_TEST( fit.crossesBand(),
	            "extendedNodes is positive but crossesBand() is false" );
	BOOST_TEST( fit.deepestBandNode > 0.0,
	            "a band node reports a depth of " << fit.deepestBandNode );

	int flaggedButInside = 0;
	int outsideButUnflagged = 0;
	for ( std::size_t j = 0; j < count; ++j )
	{
		mfem::DenseMatrix matrix( 2, 1 );
		matrix( 0, 0 ) = fit.pointR[ j ];
		matrix( 1, 0 ) = fit.pointZ[ j ];
		mfem::Array<int> found;
		mfem::Array<mfem::IntegrationPoint> ips;
		solved.mesh().FindPoints( matrix, found, ips, false );

		bool const inside = found.Size() > 0 && found[ 0 ] >= 0;
		if ( fit.extended[ j ] && inside )
			++flaggedButInside;
		if ( !fit.extended[ j ] && !inside )
			++outsideButUnflagged;
	}

	std::printf( "\n  The band through the angle fit, k = %d, n = %d, "
	             "Psi_N = %.2f: %d of %zu nodes extended, deepest %.3e "
	             "( %.3f h ), transversality %.3f\n", order, n, fraction,
	             fit.extendedNodes, count, fit.deepestBandNode,
	             fit.deepestBandNode/solved.h(), fit.transversality );
	std::fflush( stdout );

	BOOST_TEST( flaggedButInside == 0,
	            flaggedButInside << " nodes are flagged as extended while "
	            "Mesh::FindPoints puts them inside the mesh" );
	BOOST_TEST( outsideButUnflagged == 0,
	            outsideButUnflagged << " nodes are outside the mesh and NOT "
	            "flagged, which is the defect the mask exists to prevent: IN-2 "
	            "would average them in with the solved ones and report one rate "
	            "over two populations" );
}

/*
 * AND THE BAND MUST WORK WITH THE TRACER'S OWN DEFAULT PAIRING, WHICH IS NOT
 * THE ONE THE RATES ABOVE ARE MEASURED ON.
 *
 * Every band table in this file roots psi_h and steps with q_h, deliberately:
 * k+1 and k+2 are the two orders being told apart and a study of the band has
 * to say which field it is the band of. But Potential::PostProcessed is the
 * DEFAULT, and both extensions read whichever pair the tracer was constructed
 * with -- the lift integrates q* and the Taylor step steps with it. A band
 * written for one pair and quietly wrong on the other would be invisible in
 * every table above, so it is checked here: the same surface, both pairings,
 * and psi* must be at least as close to Gamma as psi_h is.
 */
BOOST_AUTO_TEST_CASE( theBandWorksOnThePostProcessedPairingToo )
{
	int const order = 2;
	int const n = 32;

	CurvedSolve solved( order, n );
	solved.theSolver().postProcess();
	meq::CriticalPoint const axis = solved.axis();

	double error[ 2 ] = { 0.0, 0.0 };
	int extended[ 2 ] = { 0, 0 };
	meq::Potential const pairing[ 2 ] = { meq::Potential::Raw,
	                                      meq::Potential::PostProcessed };

	for ( int which = 0; which < 2; ++which )
	{
		meq::ContourTracer tracer( solved.theSolver(), pairing[ which ] );
		useBandExtension( tracer, solved, meq::BandExtension::TransferLift );
		useFixedStep( tracer, 0.02 );
		tracer.setMaxPoints( 20000 );

		meq::Contour const contour = tracer.traceFromAxis( 0.0, axis );
		BOOST_TEST_REQUIRE( contour.closed(),
		                    ( which == 0 ? "psi_h" : "psi*" )
		                    << ": Gamma ended as '"
		                    << meq::contourStatusName( contour.status ) << "'" );

		BandMeasurement const split = splitByBand( contour );
		error[ which ] = split.bandWorst;
		extended[ which ] = split.bandPoints;
	}

	std::printf( "\n  Gamma in the band, both pairings, k = %d, n = %d: "
	             "psi_h %.6e over %d band points, psi* %.6e over %d\n",
	             order, n, error[ 0 ], extended[ 0 ], error[ 1 ], extended[ 1 ] );
	std::fflush( stdout );

	BOOST_TEST( extended[ 1 ] > 0,
	            "the post-processed pairing traced Gamma without entering the "
	            "band, so the extension is not being reached at all on the "
	            "tracer's own default" );
	BOOST_TEST( error[ 1 ] <= error[ 0 ],
	            "the post-processed pairing is " << error[ 1 ] << " from Gamma "
	            "against psi_h's " << error[ 0 ] << ". psi* carries a higher "
	            "order and q* agrees with its gradient an order better, so the "
	            "band cannot be costing it accuracy unless the extension is "
	            "reading the wrong field" );
}
