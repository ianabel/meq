#ifndef MEQ_SURFACEFIT_HPP
#define MEQ_SURFACEFIT_HPP

#include <cstddef>
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
 * continuous integration -- which cannot obtain the MFEM branch meq needs --
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
 * than an argument, and it is IN-4's to resolve.
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
 * The design matrix is reduced by Householder QR and its triangular factor is
 * decomposed by a one-sided Jacobi rotation sweep, which gives the singular
 * values of the design matrix itself to high relative accuracy and a solve that
 * squares nothing. Normal equations would square the condition number, so a
 * design matrix at 1e8 -- which an ill-placed sample set reaches easily -- would
 * lose every digit through a route whose only symptom is a fit that is merely
 * poor. Modes whose singular value falls below a relative floor are DISCARDED
 * rather than inverted, and the count is reported.
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

	/// Every mode of the tensor-product control at degree at most @a maxDegree:
	/// every ( p, m ) with p >= 0 and p + |m| <= maxDegree, in the same
	/// degree-then-m order zernikeModes() uses. Exposed so that a test can count
	/// them -- the control carries ( L + 1 )^2 modes against Zernike's
	/// ( L + 1 )( L + 2 )/2, so it is a GENEROUS control and anything it does
	/// worse it does worse with more freedom, not less.
	std::size_t tensorProductModeCount( int maxDegree );

}

#endif // MEQ_SURFACEFIT_HPP
