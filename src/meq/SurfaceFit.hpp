#ifndef MEQ_SURFACEFIT_HPP
#define MEQ_SURFACEFIT_HPP

#include <cstddef>
#include <functional>
#include <vector>

#include "Zernike.hpp"

/*
 * The representation: R and z as truncated Zernike expansions on the disc,
 * fitted to traced surface POINTS. INVERSION-PLAN.md stage IN-3.
 *
 * ---------------------------------------------------------------------------
 * A NAME COLLISION THAT IS SETTLED BEFORE ANYTHING ELSE, BECAUSE BOTH
 * QUANTITIES MEET IN THIS FILE AND ONE IS A FUNCTION OF THE OTHER.
 * ---------------------------------------------------------------------------
 *
 * In src/meq/Zernike.hpp and in INVERSION-PLAN.md section 4.1, `rho` is the
 * DISC RADIAL COORDINATE, rho = sqrt( Psi_N ): a FLUX label, dimensionless,
 * zero at the magnetic axis and one at the edge of the disc.
 *
 * In meq::AngleParametrisation -- which is one of this file's two natural
 * inputs -- `rho` is the GEOMETRIC DISTANCE in metres from the magnetic axis
 * to the surface along a ray, and `radiusPrime` is its theta-derivative.
 *
 * THEY ARE NOT THE SAME QUANTITY AND NEITHER IS A REPARAMETRISATION OF THE
 * OTHER: the first labels which surface, the second says how far away that
 * surface is in a given direction. Conflating them produces a fit that is
 * plausible-looking and wrong everywhere except on the one surface where the
 * two happen to agree numerically, and nothing in a coefficient table says so.
 *
 * So NOTHING in this file is called rho. The disc coordinate is
 * `discRadius` and the flux label it is built from is `normalisedFlux`; the
 * geometric distance does not appear here at all, because a sample is given as
 * a POSITION ( r, z ) and never as a distance. The tree has settled one such
 * collision already -- q is the flux and the safety factor is safetyFactor,
 * see src/meq/SurfaceAverage.hpp -- and this is the second.
 *
 * ---------------------------------------------------------------------------
 * IT FITS POINTS. IT IS NOT A PROJECTION.
 * ---------------------------------------------------------------------------
 *
 * Zernike.hpp deliberately carries no fitting routine and no quadrature, and
 * says why: IN-3 fits traced surface POINTS, which is a linear least-squares
 * problem against a point cloud, not an integral against a callable. The input
 * here is therefore a set of SurfaceSample -- each one a position ( r, z ) at a
 * known normalised flux and a known poloidal angle -- and that is exactly what
 * ContourTracer::fitByAngle() produces per surface, and what
 * meq::analytic::surfaceQuadrature() produces on an exact field.
 *
 * A caller with an AngleParametrisation writes the loop itself:
 *
 *     for ( std::size_t j = 0; j < fit.count(); ++j )
 *         samples.push_back( { psiN, 2*pi*j/fit.count(),
 *                              fit.pointR[ j ], fit.pointZ[ j ] } );
 *     ...
 *     samples = relabelByAxisShape(
 *         samples, axisShapeFromSamples( samples, axis.r, axis.z ) );
 *
 * rather than this file taking a dependency on FluxSurfaces.hpp, and therefore
 * on MFEM, for a loop. THE LAST LINE IS NOT OPTIONAL and the next section is
 * why: without it the fit converges algebraically. NOTHING HERE INCLUDES MFEM, for the same two reasons
 * Zernike.hpp, Profiles.hpp and Source.hpp do not: the fit is pure arithmetic
 * on plain doubles, so it is testable without the finite element library, and
 * continuous integration -- which cannot obtain the MFEM branch MEQ needs --
 * can build and test it.
 *
 * ---------------------------------------------------------------------------
 * THE ANGLE IS A LABEL AND THE CHOICE OF LABEL DECIDES WHETHER THE BASIS
 * CONVERGES AT ALL. THE OBVIOUS CHOICE IS THE WRONG ONE.
 * ---------------------------------------------------------------------------
 *
 * SurfaceSample::theta is whatever the caller says it is. It is NOT required
 * to be the geometric poloidal angle, and the headline finding of stage IN-3 is
 * that it had better not be.
 *
 * The geometric poloidal angle about the magnetic axis is the obvious choice --
 * it is what ContourTracer::fitByAngle() and meq::analytic::surfaceQuadrature()
 * both produce, it is consistent across surfaces because the axis is one point,
 * and it needs no extra information. It is also, for any surface that is not a
 * CIRCLE, a parametrisation in which the disc map is NOT SMOOTH AT THE AXIS,
 * and no basis that is smooth there can converge against it.
 *
 * The argument is one line. A function that is smooth in Cartesian ( x, y ) has
 * an expansion f( 0 ) + a x + b y + O( rho^2 ), so its coefficient of rho^1 is
 * exactly alpha cos( theta ) + beta sin( theta ) -- ONE angular harmonic. Take
 * a family of nested similar ELLIPSES with semi-axes ( A, B ) and parametrise
 * each by the geometric angle: then
 *
 *     R - R_axis = discRadius * A B cos( theta )
 *                  / sqrt( A^2 sin^2 theta + B^2 cos^2 theta )
 *
 * whose coefficient of discRadius^1 carries cos( theta ), cos( 3 theta ),
 * cos( 5 theta ) and so on for ever. Every one of those beyond the first is a
 * mode the Zernike index constraint EXCLUDES at radial degree one, precisely
 * because it is not smooth at the origin -- so what is excluded is exactly what
 * the parametrisation has put there.
 *
 * MEASURED, ON THE ELLIPSES THEMSELVES, WHERE THE ANSWER IS KNOWN EXACTLY.
 * Parametrised by the ellipse's own parameter the family is a POLYNOMIAL of
 * degree one in ( x, y ) and the fit reproduces it to 2.9e-15 at maxDegree = 2.
 * The IDENTICAL points relabelled by their geometric angle are wrong by
 * 1.5e-01 at degree 2 and still wrong by 1.1e-02 at degree 20, decaying like
 * L^-1.2. Same points, same basis, same solver; only the label differs.
 *
 * SO THIS FILE SUPPLIES THE RELABELLING, AND IT IS COMPUTED FROM THE SURFACES
 * RATHER THAN ASSUMED. Near the axis every equilibrium's surfaces are ellipses
 * -- psi has a non-degenerate critical point there, so its level sets are the
 * level sets of a quadratic form -- and that form is available either from the
 * Hessian of psi at the axis or, with no field at all, by fitting a quadratic
 * form to the innermost traced surface. AxisShape is that ellipse,
 * axisShapeFromHessian() and axisShapeFromSamples() are the two routes to it,
 * and relabelByAxisShape() turns geometric angles into the label that makes the
 * leading term smooth. On SolovievEquilibrium::nstx() it is worth between 45x
 * and 660x in the worst fit error, at no cost in conditioning.
 *
 * IT IS AN EXACT REPARAMETRISATION AND NOT AN APPROXIMATION. The map from the
 * geometric angle to the label is a fixed bijection of the circle, so a caller
 * asking "where is this surface at geometric angle theta" evaluates the fit at
 * shapedPoloidalAngle( shape, theta ) and gets the right point, at every
 * radius. What is approximate is only the CLAIM that this label makes the map
 * smooth, and that claim is exact at leading order and no better.
 *
 * AND THAT LIMIT IS STRUCTURAL, WHICH IS THIS STAGE'S SECOND FINDING. A
 * relabelling that does not depend on discRadius has exactly one function's
 * worth of freedom, and it can therefore fix exactly one order: it makes the
 * discRadius^1 coefficient a single harmonic and leaves discRadius^2,
 * discRadius^3 and the rest carrying whatever harmonics the surface shaping
 * puts there. A parametrisation that is smooth to ALL orders has to vary with
 * the surface, which is what a code that SOLVES for its coefficients gets for
 * free and a post-hoc fit does not. Measured, the consequence is exactly where
 * it should be: with the innermost surface at Psi_N = 0.10 the coefficient
 * envelope decays geometrically, and pulling it in to Psi_N = 0.02 leaves an
 * algebraic tail that no degree removes. That is INVERSION-PLAN.md section
 * 4.4's tension -- axis regularity against everything else -- as a table rather
 * than an argument.
 *
 * AND IT IS RESOLVED, LOWER DOWN THIS FILE, BY TAKING THE ANGLE OUT OF THE
 * PROBLEM STATEMENT ALTOGETHER. meq::gaugeFreeFit() asks each disc node only to
 * LAND ON ITS SURFACE and lets the truncated basis choose where along it to sit,
 * which is what a code that SOLVES for its coefficients gets for free. Measured
 * on nstx() at degree 16 with the innermost surface at Psi_N = 0.02, the worst
 * distance from a node to the surface it belongs on falls from 3.78e-04 to
 * 4.16e-09 -- ninety thousand times -- and the coefficient envelope resumes
 * geometric decay. See the gauge-free section below; everything above this line
 * is still the warm start it is built on, and is not superseded by it.
 *
 * THE OTHER COST IS NAMED TOO: a ray angle of any kind assumes the surface is
 * STAR-SHAPED about the axis, which is the hypothesis
 * AngleParametrisation::transversality measures and refuses on. This file does
 * not re-check it -- it never sees a field -- and a caller that took its
 * samples from fitByAngle() has already had it checked.
 *
 * ---------------------------------------------------------------------------
 * THE RADIAL EXTENT: WHAT ROD THE DISC IS CUT TO, AND WHY THERE IS AN OPTION.
 * ---------------------------------------------------------------------------
 *
 * The disc is discRadius in [ 0, 1 ]. Section 4.1 defines discRadius =
 * sqrt( Psi_N ), so discRadius = 1 is the SEPARATRIX -- and a fit almost never
 * has data there, because the separatrix carries an X-point where the surface
 * develops a corner and every flux-surface quantity diverges logarithmically.
 * So the sample set covers an ANNULUS and there are two ways to place it:
 *
 *   discEdge = 1        FIT ON A PARTIAL DISC. discRadius = sqrt( Psi_N )
 *                       literally, so a coefficient means what section 4.1 says
 *                       it means and two fits with different outermost surfaces
 *                       are directly comparable. THE DEFAULT, on that ground and
 *                       on no other.
 *
 *   discEdge = Psi_max  RESCALE ONTO THE UNIT DISC. discRadius =
 *                       sqrt( Psi_N / Psi_max ), so the outermost fitted surface
 *                       sits at discRadius = 1 and the samples reach the edge.
 *                       The coefficients are then an expansion in a DIFFERENT
 *                       coordinate: they are not comparable against a fit with
 *                       another Psi_max, and a consumer reading them has to be
 *                       told which -- which is why majorRadiusExpansion() hands
 *                       back an object that does not record it only for the
 *                       default.
 *
 * AND THE TWO FIT THE SAME FUNCTION, WHICH IS THE PART THAT IS EASY TO GET
 * WRONG IN EITHER DIRECTION. A Zernike expansion of degree L spans exactly the
 * polynomials of degree L in ( x, y ), and that space is closed under scaling
 * ( x, y ) -- so rescaling the disc is a change of BASIS and not a change of
 * model, and the two fits agree to every digit. Measured: identical worst
 * errors, to better than 1e-6 relative, at every inner limit tried.
 *
 * SO THE CHOICE IS PURELY CONDITIONING, AND IT IS WORTH UP TO ELEVEN THOUSAND.
 * With the innermost surface at Psi_N = 0.02 the partial disc conditions at
 * 8.9e4 and the rescaled one at 7.8. RESCALE UNLESS YOU NEED THE COEFFICIENTS
 * TO BE COMPARABLE ACROSS FITS; it costs nothing and the default is what it is
 * because a silent change of coordinate is worse than a conditioning number a
 * caller can read. SurfaceFitDiagnostics::conditionNumber is that number, and
 * it is the condition number of the DESIGN MATRIX rather than of the normal
 * matrix -- see the solve below; the two differ by a square.
 *
 * ---------------------------------------------------------------------------
 * THE SOLVE: QR AND A SINGULAR VALUE DECOMPOSITION, NEVER NORMAL EQUATIONS.
 * ---------------------------------------------------------------------------
 *
 * Zernike is orthogonal on the disc under the uniform measure, so a
 * well-distributed sample gives a well-conditioned problem -- but the sample is
 * whatever the caller traced, which is equispaced in theta and in whatever
 * levels were asked for, and that is NOT the uniform measure on the disc. The
 * conditioning is therefore a measurement and not an assumption.
 *
 * The design matrix is reduced by a column-pivoted Householder QR and its
 * triangular factor is decomposed by a one-sided Jacobi rotation sweep, which
 * gives the singular values of the design matrix itself to high relative
 * accuracy and a solve that squares nothing. Normal equations would square the
 * condition number, so a design matrix at 1e8 -- which an ill-placed sample set
 * reaches easily -- would lose every digit through a route whose only symptom is
 * a fit that is merely poor. Modes whose singular value falls below a relative
 * floor are DISCARDED rather than inverted, and the count is reported.
 *
 * THAT PAIR IS Eigen::JacobiSVD AND NOT A HAND-ROLLED ONE, since 2026-09-03.
 * SurfaceFit.cpp used to carry both by hand, about a hundred and fifty lines,
 * and the replacement was NOT a performance decision -- measured, the two are
 * within a factor of two on MEQ's own workload and the whole solve is under half
 * a per cent of the extraction chain. It is the standing preference for a
 * maintained library over a hand-rolled algorithm. The reason JacobiSVD rather
 * than the faster BDCSVD or CompleteOrthogonalDecomposition is in the comment on
 * decompose() in the .cpp, and it turns on IN-4's finding that the gauge is a
 * soft tail with no gap -- which needs the small singular values to relative
 * rather than absolute accuracy. Eigen is header only, so this file's promise of
 * plain doubles in and coefficients out is unaffected: nothing here mentions it.
 *
 * AND THE PLAUSIBLE ADVICE ABOUT WHERE TO PUT THE SURFACES IS WRONG, WHICH IS
 * WORTH RECORDING BECAUSE IT SURVIVED INTO A DRAFT OF THIS HEADER. The disc
 * measure is discRadius d( discRadius ) = d( Psi_N )/2, so GAUSS-LEGENDRE
 * LEVELS IN Psi_N with equispaced angles ought to make the discrete inner
 * product the continuous one and the fit a projection at a condition number of
 * order one. Measured, they do not: over equispaced-in-Psi_N,
 * equispaced-in-discRadius and Gauss-Legendre-in-Psi_N layouts on the same
 * surfaces the three condition numbers agree to within twenty per cent at every
 * inner limit tried, and Gauss is sometimes the WORST of the three. The
 * orthogonality argument needs the nodes to span the whole disc, and a sample
 * set has a HOLE in the middle -- nobody traces the surface at Psi_N = 0,
 * because it is a point.
 *
 * WHAT ACTUALLY DECIDES THE CONDITIONING IS THAT HOLE AND THE DISC EDGE, AND
 * BOTH BY ORDERS. On SolovievEquilibrium::nstx() at degree 16, with the edge
 * rescaled onto the outermost surface, the condition number reads 7.8, 37, 319
 * and 7.3e4 as the innermost surface is pulled out from Psi_N = 0.02 to 0.05,
 * 0.10 and 0.25 -- four orders for a hole of 0.18 to 0.65 of the radius. And
 * leaving the edge at Psi_N = 1 with the innermost surface at 0.02 costs a
 * factor of ELEVEN THOUSAND, 8.9e4 against 7.8, for a fit that is the same
 * function to every digit. tests/convergence/SurfaceFitConvergence.cpp is that
 * table.
 *
 * ---------------------------------------------------------------------------
 * THE CONTROLS ARE FIRST-CLASS MEMBERS AND NOT TEST SCAFFOLDING.
 * ---------------------------------------------------------------------------
 *
 * Two of this file's design decisions are exactly the kind CLAUDE.md records as
 * "invisible to a convergence table" -- they change the RATE and nothing in the
 * output says why -- so each has a losing alternative that is buildable from
 * here, in the way meq::SurfaceAverages keeps its differenced metric and
 * meq::BandExtension keeps its flux Taylor step.
 *
 *   FitRadialCoordinate::NormalisedFlux   feeds the basis Psi_N instead of
 *                                         sqrt( Psi_N ). Zernike.hpp argues
 *                                         that this puts a square-root branch
 *                                         point at the axis and that every
 *                                         basis then converges ALGEBRAICALLY
 *                                         against it. That argument is now a
 *                                         measurement and this enum is what
 *                                         makes it one.
 *
 *   FitBasis::TensorProduct               any power of the radial coordinate
 *                                         times any Fourier mode in theta, with
 *                                         no parity constraint -- so it admits
 *                                         discRadius^2 cos( theta ), which is
 *                                         C^1 and not C^2 at the origin, AND it
 *                                         admits cos( theta ) itself, which does
 *                                         not even have a limit there. The
 *                                         second is what makes the control
 *                                         decisive at the axis: under it the
 *                                         fitted magnetic axis is a CURVE
 *                                         rather than a point.
 *
 * NEITHER IS AN ANSWER. Both are refused by nothing and returned by nothing
 * except a caller that asked for them by name.
 *
 * ---------------------------------------------------------------------------
 * THE AXIS COMES OUT RIGHT FOR FREE, AND "FOR FREE" IS A STRUCTURAL STATEMENT.
 * ---------------------------------------------------------------------------
 *
 * R_l^m contains only the powers discRadius^|m| ... discRadius^l, so every mode
 * with m != 0 VANISHES at discRadius = 0 and what is left there is
 * sum_l c_{l,0} R_l^0( 0 ) -- one number, with no theta in it. So a Zernike fit
 * puts every surface's collapse point at the same place BY CONSTRUCTION, to
 * round-off, whatever the coefficients are and however badly the fit was
 * conditioned. There is no special case, no limit to take and no ray to choose.
 *
 * WHAT IS NOT FREE, AND IS THE MEASUREMENT WORTH TAKING, IS WHETHER THAT POINT
 * IS THE MAGNETIC AXIS. It is an extrapolation: the sample set has a hole in
 * the middle, since a traced surface at Psi_N = 0 is a point and nobody traces
 * one. axis() is that extrapolation and the suite measures how far it lands
 * from the axis the field's own critical-point finder reports.
 *
 * ---------------------------------------------------------------------------
 * WHAT THIS IS NOT.
 * ---------------------------------------------------------------------------
 *
 * It is not a fit of OPEN surfaces. Every mode here is periodic in theta and
 * regular at the disc centre; an open surface terminating on the domain
 * boundary has endpoints and wants Chebyshev, which is INVERSION-PLAN.md IN-5
 * and is deferred with free boundary.
 *
 * It is not the psi-varying element of section 4.4. This is the GLOBAL Zernike
 * candidate of that section -- the simplest of the three -- and IN-4 is where it
 * is measured against the other two. What this file supplies to that decision is
 * the global column and the machinery to build the others'.
 *
 * It does not know where its samples came from, so it cannot flag a surface as
 * band-limited. A caller mixing surfaces that crossed the band between Gamma_h
 * and Gamma with surfaces that did not is fitting two different accuracies into
 * one expansion, and INVERSION-PLAN.md section 4.3 says that is a decision to be
 * taken deliberately. Contour::crossesBand() and
 * AngleParametrisation::crossesBand() are where the caller finds out.
 */

namespace meq
{

	/// Which radial coordinate the basis is fed. See the header: this is a
	/// CONTROL as much as an option, and the losing branch is the one
	/// Zernike.hpp argues against.
	enum class FitRadialCoordinate
	{
		/// discRadius = sqrt( Psi_N / discEdge ). THE ANSWER.
		DiscRadius,

		/// Psi_N / discEdge, fed to the same basis. THE CONTROL: it puts a
		/// square-root branch point at the axis, and the decay goes algebraic.
		NormalisedFlux
	};

	/// Which set of modes. See the header.
	enum class FitBasis
	{
		/// l >= |m|, l - |m| even. THE ANSWER, and the parity constraint is the
		/// whole reason the disc centre is an ordinary point.
		Zernike,

		/// a^p cos( m theta ) and a^p sin( |m| theta ) for every p + |m| <= L
		/// with NO parity constraint and no p >= |m| constraint -- the naive
		/// tensor product of "any polynomial in the radius" with "any Fourier
		/// mode in the angle". THE CONTROL.
		TensorProduct
	};

	/// "Zernike", "tensor product", "disc radius", "normalised flux". For
	/// printing a table that says which column is which.
	char const *fitBasisName( FitBasis basis );
	char const *fitRadialCoordinateName( FitRadialCoordinate coordinate );

	/**
	 * One traced point: a position, the surface it is on, and the poloidal angle
	 * it was found at.
	 *
	 * A plain aggregate, in the style of meq::ZernikeMode and meq::SurfaceNode.
	 * The angle is the geometric poloidal angle about the magnetic axis; the
	 * flux label is Psi_N, zero on the axis and one on the separatrix, which is
	 * the normalisation tests/analytic/FluxSurfaceReference.hpp uses and the one
	 * INVERSION-PLAN.md is written in.
	 *
	 * THE ANGLE IS NOT REQUIRED TO BE EQUISPACED and the samples are not
	 * required to lie on a tensor grid: the fit is a least-squares problem
	 * against a cloud and the cloud may be anything. What the conditioning
	 * actually turns on is measured under THE SOLVE above, and it is not the
	 * layout of the levels.
	 *
	 * AND theta IS A LABEL RATHER THAN NECESSARILY THE GEOMETRIC ANGLE. See the
	 * angle section of this header: a fit whose samples carry the geometric
	 * poloidal angle converges algebraically and one whose samples have been put
	 * through relabelByAxisShape() does not.
	 */
	struct SurfaceSample
	{
		double normalisedFlux = 0.0;  ///< Psi_N of the surface this point is on
		double theta = 0.0;           ///< poloidal angle about the magnetic axis
		double r = 0.0;
		double z = 0.0;

		/// Least-squares weight. One unless a caller has a reason; a surface
		/// known only to the band extension's order is the obvious one.
		double weight = 1.0;
	};

	/**
	 * The ellipse the flux surfaces collapse onto at the magnetic axis.
	 *
	 * psi has a non-degenerate critical point at the axis, so its level sets
	 * there are the level sets of a quadratic form and every equilibrium's
	 * innermost surfaces are ellipses. Two numbers describe one up to scale --
	 * how it is turned and how flat it is -- and those two are exactly what the
	 * angular relabelling needs. See the header on why the relabelling is the
	 * headline finding of this stage.
	 */
	struct AxisShape
	{
		/// The angle from the +r direction to the ellipse's SHORT semi-axis,
		/// which is the eigenvector of the LARGER eigenvalue of the quadratic
		/// form. In radians.
		double tilt = 0.0;

		/// short semi-axis over long semi-axis, so it lies in ( 0, 1 ] and is
		/// exactly 1 for a circular axis -- where the relabelling below is the
		/// identity, which is the cheapest check that it is the right way up.
		/// It is the reciprocal of the near-axis ELONGATION.
		double semiAxisRatio = 1.0;
	};

	/// The axis ellipse from the Hessian of psi at the axis. The Hessian is
	/// symmetric, so three numbers describe it.
	/// @throws std::invalid_argument if the form is degenerate or indefinite --
	///         at a genuine magnetic axis it is neither, and a caller that has
	///         handed this a SADDLE should be told rather than given an angle.
	AxisShape axisShapeFromHessian( double dRR, double dRZ, double dZZ );

	/// The same ellipse WITHOUT any field at all: the quadratic form that best
	/// fits the innermost surface of @a samples about ( axisR, axisZ ), by
	/// linear least squares in its three coefficients.
	///
	/// THIS IS THE ROUTE TO PREFER, and not because it is more accurate -- the
	/// Hessian route is exact and this one is the innermost surface's own shape,
	/// which differs from it at O( Psi_N ). It is the route to prefer because it
	/// needs nothing but the points already in hand: no field, no derivative, no
	/// MFEM, and no second opinion about where the axis is.
	///
	/// @throws std::invalid_argument if the innermost surface carries fewer than
	///         three points, or if its best form is degenerate.
	AxisShape axisShapeFromSamples( std::vector<SurfaceSample> const &samples,
	                                double axisR, double axisZ );

	/// The label that makes the leading-order disc map smooth, from the
	/// geometric poloidal angle. The identity when the axis is circular.
	double shapedPoloidalAngle( AxisShape const &shape, double geometricAngle );

	/// Its exact inverse: the geometric poloidal angle a label corresponds to.
	double geometricPoloidalAngle( AxisShape const &shape, double label );

	/// A copy of @a samples with every theta relabelled by
	/// shapedPoloidalAngle(). The positions are untouched -- this is a change of
	/// parametrisation and not of geometry, which is why it can be applied after
	/// the tracing rather than during it.
	std::vector<SurfaceSample> relabelByAxisShape(
		std::vector<SurfaceSample> const &samples, AxisShape const &shape );

	// -----------------------------------------------------------------------
	// THE GAUGE-FREE FIT. INVERSION-PLAN.md section 4.4, stage IN-4.
	// -----------------------------------------------------------------------
	//
	// EVERYTHING ABOVE PINS THE ANGLE AS AN INPUT, AND THAT IS THE WHOLE OF
	// IN-3's SECOND FINDING. A SurfaceSample carries a theta, so the fit is
	// asked to put a particular point at a particular label; a relabelling that
	// does not depend on discRadius straightens exactly one order, and below
	// Psi_N ~ 0.1 an algebraic tail appears that no degree removes.
	//
	// DESC does not have that problem, and Panici et al. say why in one
	// sentence: "DESC, while not explicitly enforcing any poloidal angle
	// constraints, ends up finding an optimal representation through the course
	// of the optimization procedure." A SOLVE is free to slide the angle to
	// whatever its truncated basis represents best. A FIT is not -- unless it is
	// asked for less.
	//
	// So ask for less. Require each disc node only to LAND ON THE RIGHT SURFACE:
	//
	//     minimise over the coefficients of R and z
	//         sum_j [ Psi_N( x( discRadius_j, theta_j ) ) - Psi_j ]^2
	//
	// with x the expansion being solved for. Nothing says where along its
	// surface a node must sit, so the angle is free and the truncated basis may
	// choose it. It needs no force balance and no second physics solver, because
	// MEQ already has psi; and the Jacobian needs grad Psi_N, which is the
	// SOLVED FLUX -- grad_bar psi = r q, so grad Psi_N = r q / ( psi_bnd -
	// psi_ax ). That is INVERSION-PLAN.md section 3.2's "the gradient is a
	// solved unknown" paying off a third time.
	//
	// IT IS A GEOMETRIC GAUSS-NEWTON WARM-STARTED FROM THE LINEAR FIT ABOVE, and
	// the two are not alternatives: the linear fit is what puts the iterate in
	// the right basin, and the gauge-free step is what lets it leave the
	// parametrisation it arrived in.
	//
	// THE RESIDUAL IS SCALED INTO METRES AND THAT IS NOT COSMETIC. Psi_N is
	// dimensionless and its gradient varies by orders across the disc, so a
	// least-squares problem in the flux residual weights the outer surfaces
	// against the inner ones by whatever | grad Psi_N | happens to be. Dividing
	// each row by | grad Psi_N | makes the residual the NORMAL DISTANCE from the
	// node to its surface, in metres, which is the quantity a reader of this
	// class actually cares about and the one every acceptance below is stated
	// in.
	//
	// ---------------------------------------------------------------------------
	// THE GAUGE, WHICH IS THE ONE REAL DESIGN DECISION, AND THE SPECTRUM SAYS
	// THE OBVIOUS PICTURE OF IT IS WRONG.
	// ---------------------------------------------------------------------------
	//
	// The freedom is theta -> theta + Lambda( discRadius, theta ), one function's
	// worth, and the linearised system cannot see it: a tangential displacement
	// of a node leaves it on its own surface, so it changes no residual. The
	// expectation going in was therefore that the Gauss-Newton matrix is RANK
	// DEFICIENT along the gauge, by about half its columns.
	//
	// MEASURED, IT IS NOT, AND THE ANSWER DEPENDS ON THE FIELD. On a family of
	// SIMILAR nested ellipses the exactly-null directions are a handful -- three,
	// at every degree from 2 to 16 -- and on the nstx() Solov'ev equilibrium
	// there are NONE AT ALL. What both have instead is a long SOFT TAIL with no
	// gap in it: at degree 20 on the ellipse family the singular values run
	// smoothly from 44 down to 5e-05 and then fall off a cliff to 1e-14, and on
	// nstx() at degree 16 they reach 8e-08 of the largest with no cliff beneath
	// them. The reason is that a tangential slide is only exactly null when the
	// displaced map is still IN the truncated space, and Lambda times
	// d x / d theta generally is not. So the gauge is not a subspace to be
	// projected out; it is a direction in which the problem is merely very soft.
	//
	// AND THE SOFT TAIL IS ENOUGH ON ITS OWN, which is the part that matters:
	// removing the gauge on nstx(), where nothing is exactly singular, still
	// leaves the map FOLDED and the fit five orders worse. A field with no exact
	// null space is not a field that can do without a gauge.
	//
	// THAT MATTERS BECAUSE IT MOVES WHAT THE GAUGE IS. Three treatments, and the
	// first two are the ones INVERSION-PLAN.md section 4.4 names:
	//
	//   MinimumNorm     the pseudo-inverse step: directions below a relative
	//                   singular-value floor are NOT MOVED. Among the ways to
	//                   satisfy the surface constraint it prefers the one with
	//                   the smallest step, which is a spectral-condensation
	//                   preference arrived at for free and with no weights to
	//                   choose. THE DEFAULT.
	//
	//                   AND IT IS MIN-NORM RELATIVE TO THE WARM START, NOT
	//                   ABSOLUTELY. Each STEP avoids the null directions, so the
	//                   answer stays as close as it can to the coefficients the
	//                   linear fit arrived with. That is a property to like --
	//                   the linear fit's label is a sensible place to be pinned
	//                   -- but it is not "the smallest coefficients", and a
	//                   reader who assumed it was would be wrong.
	//
	//   SpectralWidth   an explicit penalty, lambda sum |m|^spectralExponent
	//                   ( cR^2 + cz^2 ), which is Hirshman & Breslau's idea and
	//                   the constraint VMEC uses where DESC lets the solve find
	//                   the angle. It fixes the gauge outright rather than by
	//                   preference -- and it BIASES the surface residual, by
	//                   O( lambda / sigma^2 ) in the directions the constraint
	//                   does determine, so lambda is a real tuning parameter
	//                   where the floor above is a threshold on a cliff.
	//
	//                   IT LOSES ON ITS OWN METRIC AS WELL AS ON ACCURACY, which
	//                   was not expected. Panici's spectral width M( p, q ) is a
	//                   RATIO of two weighted sums of the same coefficients, and
	//                   a quadratic penalty shrinks both -- so driving the
	//                   coefficients down does not drive the ratio down.
	//                   Measured on nstx() at degree 12, twelve decades of
	//                   lambda move the width by 1.6% and in the WRONG
	//                   direction, and cost 44x in surface error; the
	//                   minimum-norm step, which asks nothing about the
	//                   spectrum, already sits below every penalised value.
	//                   Hirshman & Breslau MINIMISE M itself, which is not a
	//                   quadratic problem, and a quadratic surrogate for it is
	//                   not the same thing. Kept because it is the alternative
	//                   INVERSION-PLAN.md section 4.4 names and because a
	//                   comparison with no losing column is not a comparison.
	//
	//   None            THE CONTROL. Every direction inverted however small its
	//                   singular value, and no trust region either. See below for
	//                   why the trust region has to go with it.
	//
	// A THIRD TREATMENT EXISTS AND IT IS THE ONE THAT ACTUALLY MADE THIS ROBUST:
	// A TRUST REGION. The step is damped Levenberg-Marquardt, sigma/( sigma^2 +
	// mu ), with mu adapted -- lowered on an accepted step, raised until one is.
	// As mu falls the step tends to the plain pseudo-inverse, so on an easy
	// problem the damping is inert; where it is not, it is what stops a soft
	// direction of singular value 5e-05 being inverted into a step of 1e5.
	// Measured on the ellipse family, the undamped pseudo-inverse reaches
	// round-off at degree 12 and FAILS at degree 16 and 20 (4.6e-03 and 3.8e-03,
	// against a linear fit at 1.6e-02); damped, every degree from 2 to 16 lands
	// between 8e-16 and 2e-09.
	//
	// AND THAT IS WHY SurfaceGauge::None DISABLES THE DAMPING TOO. A trust region
	// IS a Tikhonov gauge -- damping the step is choosing among the directions
	// the constraint does not determine -- so a "no gauge" control that kept it
	// would be a control of nothing, and it would quietly pass.
	//
	// ---------------------------------------------------------------------------
	// THE SAFEGUARD, WHICH IS NOT OPTIONAL.
	// ---------------------------------------------------------------------------
	//
	// Minimising a surface residual alone admits degenerate answers that a
	// residual column cannot see: nodes can bunch, and the disc map can fold. The
	// residual of a folded map is as good as the residual of a sound one, because
	// every node is still on its surface -- it is the same species of defect as
	// the tensor-product control above, which fits the sample cloud eight times
	// better than Zernike and puts the magnetic axis on a curve. So
	// SurfaceFit::mapJacobian() is reported over the whole fitted annulus and its
	// SIGN is asserted. DESC does the equivalent, correcting boundary orientation
	// to keep the Jacobian positive.
	//
	// ---------------------------------------------------------------------------
	// WHAT IT IS NOT.
	// ---------------------------------------------------------------------------
	//
	// It is not an equilibrium solve. Nothing here enforces force balance; psi is
	// an input and the only thing being solved for is where the disc coordinates
	// put their points.
	//
	// It does NOT constrain the axis. A node at discRadius = 0 would ask for
	// Psi_N( x ) = 0, which is true at the magnetic axis and would pin it for
	// free -- except that grad Psi_N VANISHES there, so the row scaling divides
	// by zero and the linearisation carries no information about where to move.
	// That is INVERSION-PLAN.md section 2's error ( a ) and its 1/| grad psi |
	// met a third time. discNodesFrom() therefore refuses a node at Psi_N = 0,
	// and the axis stays what section "THE AXIS COMES OUT RIGHT FOR FREE" above
	// says it is: an extrapolation into the hole, whose accuracy is measured
	// rather than imposed.
	//
	// And it is MFEM-free like everything else here: the field arrives as a
	// callable, so the caller supplies either a closed form or a lambda over the
	// tracer, and this file never learns which.

	/**
	 * Psi_N and its gradient at a point, as ONE callable.
	 *
	 * ONE RATHER THAN TWO, WHICH IS A DEPARTURE FROM THE OBVIOUS DESIGN AND IS
	 * DELIBERATE. A pair of std::function<double( double, double )> -- one for
	 * the value and one for each gradient component -- reads better and costs
	 * the consumer THREE POINT LOCATIONS per node per Gauss-Newton iteration.
	 * On the discrete leg a location is meq::ContourTracer::sampleAt(), and
	 * CLAUDE.md records mfem::Mesh::FindPoints as O( elements x points ); the
	 * tracer's own seam hands back psi and q together for exactly this reason.
	 *
	 * The bool is the second half of it: a discrete field cannot answer
	 * everywhere, and a node pushed outside the mesh by a trial step must be a
	 * REFUSAL rather than a value. A refusal on the starting iterate is an
	 * error and throws; a refusal on a trial step is a rejected step, which is
	 * what the trust region is for.
	 *
	 * @param normalisedFlux out: Psi_N, zero on the axis and one on the
	 *                       separatrix, in the same normalisation the nodes
	 *                       carry.
	 * @param gradientR      out: d Psi_N / d r. For a solved field this is
	 *                       r q_r / ( psi_bnd - psi_ax ), which is the flux and
	 *                       not a differentiated potential.
	 */
	struct NormalisedFluxField
	{
		std::function<bool( double r, double z, double &normalisedFlux,
		                    double &gradientR, double &gradientZ )> sample;
	};

	/**
	 * One point of the disc the gauge-free fit constrains.
	 *
	 * IT CARRIES NO POSITION, AND THAT IS THE ENTIRE DIFFERENCE FROM
	 * meq::SurfaceSample. A sample says "the surface passes through HERE at
	 * THIS angle"; a node says only "evaluate the map at these coordinates and
	 * require the answer to be on this surface". The angle is a parameter value
	 * rather than a claim, which is what gives the fit its gauge freedom.
	 */
	struct DiscNode
	{
		double normalisedFlux = 0.0;  ///< the surface this node must land on
		double theta = 0.0;           ///< where on the disc the map is evaluated

		/// Least-squares weight, in the same sense SurfaceSample::weight has --
		/// but applied to a residual already scaled into metres.
		double weight = 1.0;
	};

	/// The ( Psi_N, theta ) labels of @a samples, with their weights, and their
	/// positions dropped. The ordinary way to build a node set: the gauge-free
	/// fit is asked to hold the surfaces the linear fit was given, and nothing
	/// about where along them its points sat.
	/// @throws std::invalid_argument if a sample sits at Psi_N = 0, where
	///         grad Psi_N vanishes and the constraint carries no information.
	std::vector<DiscNode> discNodesFrom(
		std::vector<SurfaceSample> const &samples );

	/// How the tangential freedom is fixed. See the header.
	enum class SurfaceGauge
	{
		/// Pseudo-inverse steps under a relative singular-value floor, damped by
		/// an adaptive trust region. THE ANSWER.
		MinimumNorm,

		/// Hirshman & Breslau's spectral condensation as an explicit quadratic
		/// penalty on the coefficients. Fixes the gauge outright and biases the
		/// residual; the weight is a real tuning parameter.
		SpectralWidth,

		/// Every direction inverted however small its singular value, and NO
		/// trust region. THE CONTROL, and it disables the damping as well as the
		/// floor because a trust region is itself a gauge.
		None
	};

	char const *surfaceGaugeName( SurfaceGauge gauge );

	/// What the gauge-free fit was asked for. A plain aggregate; every default
	/// here is the measured one and the header says which measurement.
	struct GaugeFreeFitOptions
	{
		SurfaceGauge gauge = SurfaceGauge::MinimumNorm;

		/// Gauss-Newton steps. A cap and not a target: the stopping rules below
		/// normally fire first, at eight to fifteen steps from a linear warm
		/// start.
		int maxIterations = 20;

		/// Stop when the worst nodal distance to its surface falls below this,
		/// IN METRES.
		///
		/// 1e-10 METRES IS A TENTH OF A NANOMETRE and is orders below anything a
		/// tokamak geometry means, below the discretisation error of any solved
		/// field, and below the tracing error of any traced one. Chasing 1e-13
		/// instead is available and costs roughly a doubling of the iteration
		/// count -- measured on the ellipse family, roughly a doubling -- for a
		/// number nothing downstream can use. It is a default and not a limit;
		/// the acceptance suite lowers it where it is measuring the ITERATION
		/// rather than the geometry.
		double tolerance = 1.0e-10;

		/// A step that improves the worst nodal distance by less than this
		/// fraction counts as buying nothing. The iteration grinds for tens of
		/// steps past the point where it is buying anything, so this and the
		/// count below are what actually terminate it.
		///
		/// Ten per cent rather than one: a Gauss-Newton step on this problem
		/// either buys orders or buys nothing, so the threshold sits in a gap
		/// and not on a slope.
		double improvement = 0.10;

		/// How many CONSECUTIVE such steps before giving up. Three rather than
		/// one, and that is a measurement: the trust region raises its damping on
		/// a rejected step and lowers it on an accepted one, so a run that has
		/// just been damped hard takes a step or two to recover its stride --
		/// stopping at the first poor step left degree 6 of the ellipse family at
		/// 3.7e-04 where it reaches 1e-13 given three. It is the same shape as
		/// meq::ContourTracer's corrector, which accepts after four
		/// non-improving steps rather than after one.
		int stallLimit = 3;

		/// Singular values below this multiple of the largest are NOT INVERTED.
		/// 1e-10 is a threshold on a cliff and not a tuning parameter: the
		/// exactly-null directions sit at 1e-14 of the largest and the softest
		/// non-null one measured is 8e-08, so anything between them gives the
		/// same answer. It is the trust region below, not this floor, that
		/// decides what happens to the soft tail. Ignored by
		/// SurfaceGauge::None.
		double singularValueFloor = 1.0e-10;

		/// lambda of the spectral-width penalty, relative to the largest squared
		/// singular value so that it is scale free. Used by
		/// SurfaceGauge::SpectralWidth alone.
		double spectralWeight = 1.0e-6;

		/// p of the penalty weight ( |m| / maxDegree )^p. NORMALISED BY THE
		/// DEGREE, which Hirshman & Breslau's bare |m|^p is not: at p = 4 and
		/// L = 16 a bare power spans five orders, so a weight tuned at one
		/// exponent means something else at another and a sensitivity sweep in p
		/// would be measuring the normalisation. Four rather than two because
		/// Panici's
		/// spectral width M( p, q ) at p = q = 2 has |m|^4 in its numerator, so
		/// this is the quadratic form the diagnostic is built on rather than a
		/// second convention. m = 0 modes carry weight zero at any p, which is
		/// Hirshman & Breslau's structure and not an omission.
		double spectralExponent = 4.0;

		/// The trust region's initial damping, relative to the largest squared
		/// singular value. Small enough that the first step is the plain
		/// pseudo-inverse wherever that works.
		double damping = 1.0e-12;

		/// How many times the damping may be raised before a step is abandoned.
		int dampingRetries = 30;
	};

	/// What the gauge-free fit did. Printed, and asserted on where the header
	/// says a property rather than a number is at stake.
	struct GaugeFreeFitReport
	{
		int iterations = 0;
		bool converged = false;

		/// Why it stopped: "tolerance", "no improvement", "iteration cap" or
		/// "step rejected". A string rather than an enum because it is printed
		/// and never branched on.
		char const *stop = "";

		/// The worst distance from a node to the surface it is supposed to be
		/// on, IN METRES, before and after. This is the quantity the whole
		/// exercise is about and it is GAUGE INVARIANT -- it says nothing about
		/// where along its surface a node sits, which is exactly the freedom
		/// being granted.
		double initialSurfaceError = 0.0;
		double surfaceError = 0.0;

		/// The WORST that quantity reached at any accepted iterate, including
		/// the first. An iteration that converged reads this equal to
		/// initialSurfaceError; one that WANDERED reads it far above, which is
		/// the single number the ungauged control exists to produce. A final
		/// error alone cannot tell the two apart, because an iteration may
		/// return from an excursion.
		double worstExcursion = 0.0;

		/// Panici eq. ( 6 ) at p = q = 2, before and after. DESC's own
		/// diagnostic and the one their Figure 5 is drawn in.
		double initialSpectralWidth = 0.0;
		double spectralWidth = 0.0;

		/// Euclidean norm of the coefficient vector, both components together,
		/// before and after. The blunt instrument that catches a gauge that has
		/// let the coefficients run away.
		double initialCoefficientNorm = 0.0;
		double coefficientNorm = 0.0;

		/// det d( R, z )/d( discRadius, theta ) over the fitted annulus, at the
		/// end. THE SAFEGUARD: these must share a sign, or the map has folded
		/// and the residual will not say so.
		double minimumJacobian = 0.0;
		double maximumJacobian = 0.0;

		/// The singular spectrum of the FIRST Gauss-Newton matrix, which is what
		/// says how rank deficient the problem really is. `nullDirections`
		/// counts singular values below 1e-10 of the largest -- the exact gauge
		/// -- and `softDirections` those below 1e-04 of it, which is the tail
		/// with no gap in it.
		std::size_t columns = 0;
		std::size_t nullDirections = 0;
		std::size_t softDirections = 0;
		double largestSingularValue = 0.0;
		double smallestSingularValue = 0.0;

		/// Where the trust region ended up, relative to the largest squared
		/// singular value OF THE FIRST Gauss-Newton matrix -- the same scale
		/// GaugeFreeFitOptions::damping is given in, so the two are directly
		/// comparable. A number that stayed at its initial value says the
		/// damping was inert and the step was the plain pseudo-inverse
		/// throughout.
		double damping = 0.0;

		/// Trial steps the trust region rejected, over the whole run.
		int rejectedSteps = 0;

		/// The Euclidean norm of the FIRST trial step, undamped. Read against
		/// initialCoefficientNorm: a gauge that inverts a singular value the
		/// data does not determine produces a correction orders larger than the
		/// thing being corrected, and this is where that shows.
		double firstStepNorm = 0.0;
	};

	/// What the fit was asked for. A plain aggregate so that a caller taking the
	/// defaults writes nothing.
	struct SurfaceFitOptions
	{
		FitBasis basis = FitBasis::Zernike;
		FitRadialCoordinate coordinate = FitRadialCoordinate::DiscRadius;

		/// The normalised flux that maps to discRadius = 1. One means
		/// discRadius = sqrt( Psi_N ) literally, which is section 4.1's
		/// definition; see the header for what the alternative buys and costs.
		double discEdge = 1.0;

		/// Singular values below this multiple of the largest are treated as
		/// zero and their directions dropped from the solution rather than
		/// inverted. 1e-10 rather than machine epsilon because the fit is
		/// against traced data whose own accuracy is far above it, so inverting
		/// a direction the data does not determine buys noise.
		double singularValueFloor = 1.0e-10;
	};

	/// What the fit cost and how trustworthy it is. Printed rather than
	/// asserted on, except the condition number, which the header says is a
	/// measurement and not an assumption.
	struct SurfaceFitDiagnostics
	{
		std::size_t samples = 0;
		std::size_t modes = 0;

		/// Of the DESIGN MATRIX, not of the normal matrix -- the two differ by a
		/// square, and quoting the wrong one understates the difficulty by half
		/// its digits.
		double largestSingularValue = 0.0;
		double smallestSingularValue = 0.0;
		double conditionNumber = 0.0;

		/// Directions the singular value floor discarded. Non-zero means the
		/// sample set does not determine the expansion at this degree, which is
		/// a statement about the samples and not about the basis.
		std::size_t discardedModes = 0;

		/// Root mean square and worst residual of the fit against its own
		/// samples, per component, in the units of r and z.
		double residualR = 0.0;
		double residualZ = 0.0;
		double worstR = 0.0;
		double worstZ = 0.0;

		/// The extent of the sample set in the disc coordinate. The first is the
		/// radius of the HOLE in the middle -- the fit is an extrapolation
		/// inside it, which is what axis() is doing.
		double smallestSampledRadius = 0.0;
		double largestSampledRadius = 0.0;
	};

	/**
	 * R( discRadius, theta ) and z( discRadius, theta ), fitted.
	 *
	 * A value type: copyable, movable, holding coefficients and nothing else. It
	 * caches no scratch, so one fit may be evaluated from several threads at
	 * once -- which is what MANTA-COUPLING.md section 5's pointwise call pattern
	 * will want.
	 */
	class SurfaceFit
	{
		public:
			/**
			 * Fit @a samples at degree @a maxDegree.
			 *
			 * @throws std::invalid_argument if the degree is out of range, if
			 *         there are fewer samples than modes, if a sample carries a
			 *         negative or non-finite normalised flux, if a sample lies
			 *         outside the disc ( Psi_N > discEdge ), or if discEdge is
			 *         not positive.
			 * @throws std::runtime_error if the design matrix has no non-zero
			 *         singular value at all, which means the samples determine
			 *         nothing.
			 */
			SurfaceFit( int maxDegree, std::vector<SurfaceSample> const &samples,
			            SurfaceFitOptions const &options = SurfaceFitOptions() );

			/// The position the fit puts on surface @a normalisedFlux at angle
			/// @a theta.
			///
			/// @a theta IS THE LABEL THE SAMPLES CARRIED, which is the geometric
			/// poloidal angle only if the samples were not relabelled. A caller
			/// that fitted relabelByAxisShape()'d samples and then asks for a
			/// geometric angle here gets the wrong point on the right surface --
			/// so ask at shapedPoloidalAngle( shape, theta ), with the same shape
			/// the samples were relabelled by.
			void position( double normalisedFlux, double theta,
			               double &r, double &z ) const;

			/// d( r, z )/d( discRadius ) at fixed theta. FINITE EVERYWHERE,
			/// including on the axis, which is why INVERSION-PLAN.md IN-3 states
			/// its acceptance in this derivative and not in the next one.
			///
			/// It is d/d( discRadius ) for BOTH coordinates: under
			/// FitRadialCoordinate::NormalisedFlux the chain factor 2 discRadius
			/// is applied here, so the two settings are directly comparable and
			/// a control column is a control of the coordinate rather than of
			/// what is being differentiated.
			void radialDerivative( double normalisedFlux, double theta,
			                       double &r, double &z ) const;

			/// d( r, z )/d( Psi_N ) at fixed theta.
			///
			/// IT DIVERGES LIKE 1/( 2 sqrt( Psi_N ) ) AT THE AXIS AND THAT IS THE
			/// COORDINATE RATHER THAN A DEFECT. INVERSION-PLAN.md section 4.4:
			/// the m = 1 coefficient of R is the minor radius, which grows like
			/// sqrt( Psi ) in EVERY equilibrium, so the derivative of a surface's
			/// POSITION with respect to flux is unbounded there universally.
			/// MaNTA's own answer is that its fluxes vanish on axis and its nodes
			/// need not sit at Psi = 0; what is left is a conditioning question,
			/// and the growth of this factor as the innermost node approaches the
			/// axis is measured in tests/convergence/SurfaceFitConvergence.cpp.
			/// Prefer radialDerivative() wherever a caller can work in
			/// discRadius.
			void fluxDerivative( double normalisedFlux, double theta,
			                     double &r, double &z ) const;

			/// d( r, z )/d( theta ) at fixed surface.
			void angularDerivative( double normalisedFlux, double theta,
			                        double &r, double &z ) const;

			/// | dx/dtheta | at fixed surface: the in-surface metric of IN-1,
			/// arrived at from the REPRESENTATION rather than from the field.
			/// The field's own route is AngleParametrisation::speed, which is
			/// sqrt( radiusPrime^2 + radius^2 ) built pointwise from q, and the
			/// two agreeing is a check on both.
			double angularSpeed( double normalisedFlux, double theta ) const;

			/// det d( R, z )/d( discRadius, theta ): whether the disc map is a
			/// map at all.
			///
			/// A FIT WITH A BEAUTIFUL RESIDUAL AND A FOLDED MAP IS THE QUIET
			/// WRONG ANSWER THIS CLASS IS MOST EXPOSED TO, and it is the
			/// gauge-free fit above that exposes it: requiring only that nodes
			/// LAND on their surfaces says nothing about their order along one,
			/// so nodes may bunch, cross, and turn the map over while every
			/// residual stays perfect. Sample this over the fitted annulus and
			/// check that it keeps ONE SIGN.
			///
			/// WHICH sign is not fixed and must not be assumed: it is positive
			/// for a theta that runs anticlockwise in ( r, z ) and negative for
			/// one that runs the other way, and both are legitimate labels. It
			/// vanishes like discRadius at the centre for the ordinary reason
			/// that polar coordinates do, so a check that includes
			/// discRadius = 0 measures that and not a fold.
			double mapJacobian( double normalisedFlux, double theta ) const;

			/// The spectral width of Panici et al. eq. ( 6 ), which is DESC's
			/// own diagnostic and the one their Figure 5 is drawn in:
			///
			///     M( p, q ) = sum |m|^( p + q ) ( cR^2 + cz^2 )
			///                 / sum |m|^p ( cR^2 + cz^2 )
			///
			/// with p >= 0 and q > 0. It is a weighted mean of |m|^q, so it has
			/// the units of a mode number and a SMALLER value means a more
			/// condensed angular spectrum -- which is what Hirshman & Breslau's
			/// constraint minimises in VMEC and what DESC reaches without being
			/// asked. Zero when every coefficient sits at m = 0, since both sums
			/// are then empty; that is a degenerate spectrum rather than an
			/// infinitely good one.
			///
			/// @throws std::invalid_argument if p < 0 or q <= 0.
			double spectralWidth( double p = 2.0, double q = 2.0 ) const;

			/// The point the fit puts at discRadius = 0, which is where every
			/// surface collapses.
			///
			/// FOR THE ZERNIKE BASIS THIS IS INDEPENDENT OF theta BY
			/// CONSTRUCTION -- see the header -- so the one-argument form is the
			/// answer and the two-argument form exists to MEASURE that, and to
			/// measure what the tensor-product control does instead, where it is
			/// a curve rather than a point.
			void axis( double &r, double &z ) const;
			void axisAtAngle( double theta, double &r, double &z ) const;

			/// max | coefficient | over every mode of radial degree @a degree,
			/// taken over both components. This is the coefficient ENVELOPE, and
			/// its decay with degree is IN-3's first acceptance: geometric for a
			/// map that is analytic on the disc, algebraic against a branch
			/// point. Zero for a degree above maxDegree().
			double coefficientEnvelope( int degree ) const;

			/// The coefficients of the r and z components, in the fit's own mode
			/// order -- which is zernikeModes() order for FitBasis::Zernike.
			std::vector<double> const &majorRadiusCoefficients() const;
			std::vector<double> const &heightCoefficients() const;

			/// The r component as a meq::ZernikeExpansion, for a consumer that
			/// wants the basis object rather than this class.
			/// @throws std::runtime_error unless the basis is FitBasis::Zernike,
			///         the coordinate is FitRadialCoordinate::DiscRadius AND
			///         discEdge is one -- a ZernikeExpansion carries no record of
			///         any of the three, so handing one back from a control fit,
			///         or from a fit on a rescaled disc, would be handing back an
			///         object whose argument means something other than
			///         sqrt( Psi_N ) with nothing to say so.
			ZernikeExpansion majorRadiusExpansion() const;
			ZernikeExpansion heightExpansion() const;

			int maxDegree() const;
			SurfaceFitOptions const &options() const;
			SurfaceFitDiagnostics const &diagnostics() const;

			/// The disc coordinate a normalised flux maps to, and the inverse.
			/// Exposed because a caller printing a table wants to label its rows
			/// with the coordinate the fit is actually in, and because the
			/// mapping carries discEdge, which is easy to forget.
			double discRadiusOf( double normalisedFlux ) const;
			double normalisedFluxOf( double discRadius ) const;

		private:
			/// One basis function, indexed the same way for both bases: the
			/// radial index is l for Zernike and the power p for the tensor
			/// product, and the angular index carries the sign convention of
			/// meq::ZernikeMode -- m >= 0 is cos( m theta ), m < 0 is
			/// sin( |m| theta ).
			struct FitMode
			{
				int radial;
				int angular;
			};

			/// A fit whose coefficients came from somewhere other than a linear
			/// least-squares solve. PRIVATE, with gaugeFreeFit() as its only
			/// caller, because the diagnostics of such a fit are not the
			/// diagnostics of a fit at all -- there is no design matrix, so no
			/// condition number and no residual against samples -- and a public
			/// constructor handing back an object with half its report unset
			/// would be exactly the mislabelling checkExpansionIsPlain() exists
			/// to refuse.
			SurfaceFit( int maxDegree, SurfaceFitOptions const &options,
			            std::vector<FitMode> modes,
			            std::vector<double> majorRadius,
			            std::vector<double> height );

			friend SurfaceFit gaugeFreeFit( SurfaceFit const &start,
			                                NormalisedFluxField const &field,
			                                std::vector<DiscNode> const &nodes,
			                                GaugeFreeFitOptions const &options,
			                                GaugeFreeFitReport &report );

			/// Refuse to expose the coefficients as a bare ZernikeExpansion
			/// unless they really are one. See the accessors above.
			void checkExpansionIsPlain( char const *where ) const;

			/// The value of one mode alone, which is all the design matrix
			/// wants. Kept separate from evaluateMode() because on the Zernike
			/// branch the derivatives cost three further special function
			/// evaluations per entry and the matrix has samples x modes of them.
			double modeValue( FitMode mode, double argument, double theta ) const;

			/// value and both derivatives of one mode, in the basis's own
			/// argument. Written as one call because the caller always wants at
			/// least two of the three and the trigonometry is shared.
			void evaluateMode( FitMode mode, double argument, double theta,
			                   double &value, double &dArgument,
			                   double &dTheta ) const;

			/// The basis argument for a normalised flux, and its derivative with
			/// respect to Psi_N. The second is 1/( 2 discEdge argument ) for the
			/// disc radius and 1/discEdge for the control, and it is the factor
			/// Zernike.hpp says is the one that gets dropped.
			double argumentOf( double normalisedFlux ) const;

			/// The three sums, at a given basis argument.
			void evaluateAll( double argument, double theta,
			                  double &r, double &z,
			                  double &dArgumentR, double &dArgumentZ,
			                  double &dThetaR, double &dThetaZ ) const;

			int degree = 0;
			SurfaceFitOptions option;
			SurfaceFitDiagnostics diagnostic;
			std::vector<FitMode> modeList;
			std::vector<double> coefficientR;
			std::vector<double> coefficientZ;
	};

	/**
	 * Refit @a start so that every node LANDS ON ITS SURFACE, with the poloidal
	 * angle free. INVERSION-PLAN.md section 4.4; see the gauge-free section of
	 * this header for the whole argument.
	 *
	 * @a start supplies the degree, the options, the mode list and the initial
	 * coefficients, and is not otherwise consulted -- in particular the samples
	 * it was built from are gone by then, which is the point: this solve knows
	 * nothing about where any surface point was, only which surfaces exist and
	 * what field defines them.
	 *
	 * THE WARM START IS NOT AN OPTIMISATION. The surface residual is minimised
	 * by any map onto the right surfaces, including folded ones and ones that
	 * traverse a surface twice, so the basin matters; the linear fit is what
	 * puts the iterate in the right one. Starting from zero coefficients is not
	 * supported and would not be meaningful.
	 *
	 * @throws std::invalid_argument if the node set is empty, if it carries
	 *         fewer rows than the 2 x modes unknowns, if a node lies outside the
	 *         disc, or if an option is out of range.
	 * @throws std::runtime_error if the field refuses a node at the STARTING
	 *         iterate -- which means the warm start is not close enough for the
	 *         field to answer, and is a fault in the caller's field or start
	 *         rather than in the iteration.
	 */
	SurfaceFit gaugeFreeFit( SurfaceFit const &start,
	                         NormalisedFluxField const &field,
	                         std::vector<DiscNode> const &nodes,
	                         GaugeFreeFitOptions const &options,
	                         GaugeFreeFitReport &report );

	/// Every mode of the tensor-product control at degree at most @a maxDegree:
	/// every ( p, m ) with p >= 0 and p + |m| <= maxDegree, in the same
	/// degree-then-m order zernikeModes() uses. Exposed so that a test can count
	/// them -- the control carries ( L + 1 )^2 modes against Zernike's
	/// ( L + 1 )( L + 2 )/2, so it is a GENEROUS control and anything it does
	/// worse it does worse with more freedom, not less.
	std::size_t tensorProductModeCount( int maxDegree );

}

#endif // MEQ_SURFACEFIT_HPP
