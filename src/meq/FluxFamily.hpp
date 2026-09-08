#ifndef MEQ_FLUXFAMILY_HPP
#define MEQ_FLUXFAMILY_HPP

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/*
 * The flux-surface family, the flux label it is expressed against, and the
 * per-psi cache: INVERSION-PLAN.md stage IN-6.
 *
 * MFEM-FREE, DELIBERATELY, AND THE REASON IS THE CONTRACT RATHER THAN THE
 * ARITHMETIC. What is in this file is a container of doubles, an interpolation
 * in one variable, and a cache whose invalidation rule is the whole of
 * MANTA-COUPLING.md section 8. That rule is the part of IN-6 with a stated
 * obligation attached -- a served answer must equal a cold one BIT FOR BIT and
 * a stale cache must be impossible rather than unlikely -- and it is testable
 * with no finite element library at all. Continuous integration cannot obtain
 * the MFEM branch meq needs, so anything that can be gated there should be, in
 * the same way and for the same reasons as Profiles, Source, Zernike,
 * SurfaceFit, Coils, ExteriorDtN and PlasmaComponent.
 *
 * The half that DOES need MFEM -- turning a solved psi_h and q_h into the
 * surfaces below -- is meq::extractFluxSurfaces() in FluxExtraction.hpp.
 *
 *
 * 1. THE FLUX LABEL IS rho = sqrt( Psi_N ), AND THAT IS A MEASUREMENT.
 *
 * Psi_N = ( psi_ax - psi ) / ( psi_ax - psi_bnd ), zero on the magnetic axis
 * and one on the plasma boundary. The geometry is NOT smooth in it: psi has a
 * quadratic maximum at the axis, so Psi_N behaves like ( distance )^2 there and
 * a surface's minor radius grows like sqrt( Psi_N ). Parametrise by Psi_N
 * directly and every quantity acquires a square-root branch point at the axis.
 *
 * Zernike.hpp argues that from the geometry and IN-3 measured it: on nstx() at
 * L = 20 the coefficient envelope reads 4.89e-05 against 8.62e-03, a ratio of
 * 176, and the worst fit error 3.44e-05 against 7.01e-03, a factor of 204 --
 * with the conditioning UNTOUCHED by the choice, which is what says the gain is
 * the branch point and not the algebra.
 *
 * So the surfaces are laid out in rho, the interpolation below is in rho, and
 * the file written by meq::FluxGridWriter carries both columns. A consumer
 * declaring a coordinate -- MANTA-COUPLING.md section 5's FieldModelSpec::label
 * -- is being handed rho.
 *
 *
 * 2. THE SEPARATRIX CUT IS A DECISION AND IT IS RECORDED HERE.
 *
 * INVERSION-PLAN.md section 8's second risk asks for exactly this: everything
 * degrades approaching the separatrix -- 1/| grad psi | in every weight, the
 * corner in the surface, the logarithm no polynomial basis catches -- and
 * production codes simply stop. LIUQE cuts at Psi_N = 0.95; FreeGS extrapolates
 * outside [ 0.01, 0.99 ]. The instruction is to decide meq's cut deliberately
 * rather than discover it as a convergence failure.
 *
 * meq's default is Psi_N in [ 0.05, 0.95 ], and the two ends are cut for
 * DIFFERENT reasons, which is why they are two numbers and not one:
 *
 *   the outer end   the surface approaches a corner and V' diverges
 *                   logarithmically. On the curved path it also crosses the
 *                   band between Gamma_h and Gamma, where the field is
 *                   continued rather than solved.
 *
 *   the inner end   the surface shrinks to a point, V' -> 0, and
 *                   d( geometry )/d( Psi_N ) diverges like 1/( 2 sqrt( Psi_N ) )
 *                   -- which belongs to the COORDINATE and not to the fit, and
 *                   is measured as such: IN-3 reads the product with
 *                   2 sqrt( Psi_N ) settling at 1.148 down to Psi_N = 0.005.
 *                   What actually fails first is the trace: traceFromAxis()
 *                   brackets the level along a ray, and near the axis the
 *                   bracket is a fraction of an element wide.
 *
 * The numbers themselves are measured on meq's own fixture rather than
 * inherited -- tests/convergence/FluxGridConvergence.cpp sweeps both ends and
 * prints where each quantity gives out, and CLAUDE.md's IN-6 section carries
 * the table. They are configurable, because the right cut depends on the mesh
 * and on what the consumer needs, and because a cut chosen by the library and
 * not sayable by the caller is a cut nobody can measure.
 *
 * A QUERY OUTSIDE THE CUT THROWS AND IS NOT EXTRAPOLATED. That is the opposite
 * of FreeGS's choice and it is MANTA-COUPLING.md section 8's first contract:
 * a state meq cannot evaluate at must throw, so the consumer's integrator
 * treats it as recoverable and retries with a smaller step. An extrapolated
 * V' past the separatrix is exactly the "plausible-looking nonsense value" that
 * converts a recoverable step into a wrong answer.
 *
 *
 * 3. THE INTERPOLATION IS A MONOTONE CUBIC, AND MONOTONE IS THE LOAD-BEARING
 *    WORD.
 *
 * A family is extracted at a few tens of levels and a consumer asks at its own
 * nodes, so something has to interpolate. An ordinary cubic through four
 * neighbours overshoots, and V' -> 0 at the axis means an overshoot there is a
 * NEGATIVE volume derivative -- a nonsense value that no consumer would flag,
 * arriving through the one route section 8 forbids. Fritsch-Carlson limited
 * slopes cannot overshoot and reproduce a cubic where the data is smooth.
 *
 * It is built on meq::HermiteCubicSpline rather than beside it: the standing
 * preference is to take the maintained implementation, and what is added here
 * is the twenty-line slope rule and nothing else.
 *
 *
 * 4. THE CACHE IS KEYED ON THE psi VECTOR, COMPARED BITWISE, AND A HASH WAS
 *    REJECTED.
 *
 * MANTA-COUPLING.md section 5: Geometry is POINTWISE -- called once per physics
 * node, per residual evaluation, handed the whole psi vector each time. Done
 * naively that locates and integrates a surface per node per residual, which
 * INVERSION-PLAN.md section 11.1 calls a requirement on this stage rather than
 * an optimisation. IN-P measured what the cache is worth: nodes / surfaces,
 * 5.1x on a sixty-node case, 6.55 s against 1.29 s.
 *
 * The key is the psi vector itself, compared entry by entry with std::memcmp.
 * A hash would be cheaper and is WRONG HERE: section 8 asks for a stale cache
 * to be impossible, and a hash makes it unlikely. The comparison is O( nDOF )
 * against an extraction that traces, fits and integrates a whole family, so it
 * is free at every size worth caching at; tests/unit/FluxFamilyTests.cpp prints
 * the cost of a served query and CLAUDE.md's IN-6 section carries the number.
 *
 * memcmp AND NOT ==, and the difference is in the safe direction. Two vectors
 * differing only in the sign of a zero compare EQUAL under == and UNEQUAL under
 * memcmp, so memcmp recomputes where == would serve; it never claims equality
 * where the bits differ. A NaN compares unequal to itself under == -- which
 * would defeat the cache outright -- and equal to itself under memcmp, which is
 * correct, the input being the same input. The non-finite guard below is what
 * stops a NaN psi being served a geometry rather than a refusal.
 *
 *
 * 5. resetForRun() IS NOT REDUNDANT WITH THE KEY, AND THAT IS THE SUBTLE HALF.
 *
 * The bitwise key catches a changed psi. It does NOT catch a changed EXTRACTOR:
 * a second run whose psi happens to start where the first one ended, against a
 * different mesh or a different solver, matches the key and is served the first
 * run's geometry. That is MANTA-COUPLING.md section 8's "a cache keyed on the
 * object rather than the run", and the document records what it cost the last
 * time: a second run that completed, looked plausible, and was wrong in the
 * eleventh digit. MaNTA pins the property with a test asserting that a reused
 * solver matches a fresh one BIT FOR BIT, zero tolerance; the corresponding
 * assertions here are in tests/unit/FluxFamilyTests.cpp.
 *
 * So the owner calls resetForRun() whenever anything the extractor closes over
 * changes. The cache cannot detect that for itself and does not pretend to.
 *
 *
 * 6. A FAILED EXTRACTION COMMITS NOTHING.
 *
 * If the extractor throws, the cache is left holding NOTHING -- not the
 * previous family under the previous key, which would serve a different state's
 * geometry the moment the caller retried at the old psi, and not the previous
 * family under the new key, which would serve it immediately. Both are fudges.
 * The exception propagates, which is what section 8 asks for and what makes the
 * consumer's integrator retry with a smaller step.
 */

namespace meq
{

	/// One flux surface: where it is, and the flux-surface averages over it.
	///
	/// The averages are computed by meq::surfaceAverages() and copied here, so
	/// that a family is a value a consumer can hold, write out and interpolate
	/// without keeping a mesh, a solver or a tracer alive behind it. That is
	/// what makes the cache below MFEM-free and it is what a NetCDF file is
	/// anyway.
	struct FluxSurface
	{
		/// psi at the surface, in meq's own units, and the two labels derived
		/// from it. Psi_N is zero on the axis and one on the plasma boundary;
		/// radial is sqrt( Psi_N ), which is the label everything here is
		/// expressed against. See the header, section 1.
		double level = 0.0;
		double normalisedFlux = 0.0;
		double radial = 0.0;

		/// The surface, at angles equispaced in the poloidal angle about the
		/// magnetic axis: theta_j = 2 pi j / N, not stored because it is implied
		/// by the index.
		std::vector<double> r;
		std::vector<double> z;

		/// 1 where the node's field came from the band extension outside the
		/// mesh rather than from an element, 0 otherwise.
		///
		/// PER NODE AND NOT A COUNT. INVERSION-PLAN.md section 4.3 is explicit
		/// that a consumer must be able to tell WHICH data came from an
		/// extension, and CLAUDE.md records that the count-not-mask version of
		/// exactly this was half of a real defect in the ( R, Z ) grid file --
		/// nothing downstream could tell which of 1667 nodes were continued.
		/// extendedNodes below is the count and it is a summary of this, never
		/// a substitute for it.
		///
		/// char rather than bool because std::vector<bool> is a bitset whose
		/// elements have no address, and this is data that is written to a file
		/// and passed around beside the columns above.
		std::vector<char> extended;

		/// V' = dV/dpsi = closed-integral 2 pi R dl / | grad psi |, with V the
		/// volume enclosed by the surface. The convention is SurfaceAverage.hpp's
		/// and the 2 pi R is IN.
		double vPrime = 0.0;

		/// The enclosed volume and the enclosed poloidal cross-section area,
		/// by Green's theorem on the same nodes: V = closed-integral pi R^2 dz
		/// and A = closed-integral R dz.
		///
		/// AN INDEPENDENT ROUTE TO V' AND THAT IS WHY THEY ARE HERE. The coarea
		/// formula gives V'( psi ) = -dV/dpsi, so the family carries two
		/// completely different integrals -- a weighted line integral with
		/// 1/| grad psi | in it, and a Green's-theorem area with no gradient
		/// anywhere -- which must agree. They share the nodes and nothing else,
		/// and the check needs no reference value. It is asserted in
		/// tests/convergence/FluxGridConvergence.cpp, Richardson-extrapolated
		/// for the reason SurfaceAverage.hpp gives at length.
		///
		/// dz/dtheta is taken POINTWISE from rho' -- which comes from the solved
		/// flux -- and never by differencing neighbouring z_j. That is
		/// INVERSION-PLAN.md section 3.2's metric trap, and it is the same trap
		/// on a new integrand.
		double volume = 0.0;
		double crossSectionArea = 0.0;

		/// closed-integral dl, the poloidal circumference, and
		/// closed-integral 2 pi R dl, the area of the toroidal surface.
		double arcLength = 0.0;
		double surfaceArea = 0.0;

		/// < R^-2 >, < | grad psi |^2 / R^2 >, < | grad psi | > and
		/// < | grad psi |^2 >. The first two are what RoPP (142)'s safety factor
		/// and the averaged Grad-Shafranov identity are built from; the last two
		/// are what < | grad rho | > and < | grad rho |^2 > are, times the
		/// analytic drho/dpsi below.
		double inverseRSquared = 0.0;
		double gradPsiSquaredOverRSquared = 0.0;
		double absGradPsi = 0.0;
		double gradPsiSquared = 0.0;

		/// V' g < R^-2 > / 4 pi^2, RoPP (142). Zero, and NOT to be read, unless
		/// the family's safetyFactorAvailable is true: g( psi ) = R B_toroidal
		/// is the caller's to supply from its own profile, a meq::Source not
		/// carrying it. NEVER named q -- see SurfaceAverage.hpp on the collision.
		double safetyFactor = 0.0;

		/// How many nodes carried band data, whether any did, and how far
		/// outside Gamma_h the deepest sat. A surface with crossesBand true is
		/// limited by the extension rather than by the discretisation and must
		/// be reported separately.
		int extendedNodes = 0;
		bool crossesBand = false;
		double deepestBandNode = 0.0;

		/// The worst | psi_h( x_j ) - level | over the nodes: INVERSION-PLAN.md
		/// section 2's error ( b ) for this surface, and the thing to read
		/// before believing any of the rest.
		double worstResidual = 0.0;

		/// min | u x t | over the fit, and whether it cleared the floor.
		/// Star-shapedness is a hypothesis and this is it, measured.
		double transversality = 0.0;
		bool transverse = true;

		/// Rays accepted at the closest point they could reach rather than at
		/// the tolerance, and times the element walk fell back on
		/// Mesh::FindPoints. Both are carried from the producing fit.
		int stalledRays = 0;
		int fallbackLocations = 0;

		std::size_t count() const
		{
			return r.size();
		}
	};

	/// A quantity of a surface: what the interpolation below is applied to.
	///
	/// A CALLABLE AND NOT AN ENUMERATION, for the reason SurfaceAverage.hpp
	/// gives: MANTA-COUPLING.md says the geometry slot list "is negotiated with
	/// the transport physics case, not fixed by MaNTA" and that its illustrative
	/// set is being revised. A list that is still moving must not be baked in.
	/// The named wrappers below are one line each.
	using SurfaceQuantity = std::function<double( FluxSurface const & )>;

	/// A family of flux surfaces at one psi: what the cache holds and what the
	/// ( Psi, theta ) file is written from.
	struct FluxSurfaceFamily
	{
		/// The surfaces, ordered by strictly increasing radial label.
		std::vector<FluxSurface> surfaces;

		/// The magnetic axis the rays were drawn from, and the two flux values
		/// the label is normalised by.
		double axisR = 0.0;
		double axisZ = 0.0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;

		/// The cut, in Psi_N, that this family was extracted between. Stored so
		/// that a family says what it was cut at rather than leaving it to be
		/// remembered -- and so that a query outside it can refuse in the terms
		/// the caller asked in. See the header, section 2.
		double innerCut = 0.0;
		double outerCut = 0.0;

		/// Nodes per surface.
		std::size_t angles = 0;

		/// Whether safetyFactor was computed. False when the caller supplied no
		/// g( psi ), in which case the field is zero and must not be read.
		bool safetyFactorAvailable = false;

		std::size_t size() const
		{
			return surfaces.size();
		}

		bool empty() const
		{
			return surfaces.empty();
		}

		/// The label range actually covered, which is sqrt() of the cut.
		double innerLabel() const;
		double outerLabel() const;

		/// dPsi_N/dpsi = -1 / ( psi_ax - psi_bnd ), and
		/// drho/dpsi = dPsi_N/dpsi / ( 2 rho ). BOTH ANALYTIC: the label is an
		/// affine function of psi and a square root of that, so neither needs a
		/// difference, and < | grad rho | > is | drho/dpsi | < | grad psi | >
		/// exactly.
		///
		/// @throws std::runtime_error from radialDerivative() at rho = 0, where
		///         drho/dpsi is unbounded. That is the coordinate and not a
		///         defect, and it is one of the reasons for the inner cut.
		double normalisedFluxDerivative() const;
		double radialDerivative( double label ) const;

		/// The quantity @a of at the label @a label, by the monotone cubic of
		/// the header's section 3.
		///
		/// @param label rho = sqrt( Psi_N ).
		/// @throws std::domain_error if @a label is outside [ innerLabel(),
		///         outerLabel() ] -- REFUSED AND NOT EXTRAPOLATED, per the
		///         header's section 2 -- or std::runtime_error if the family has
		///         fewer than two surfaces to interpolate between.
		double at( double label, SurfaceQuantity const &of ) const;

		/// The same against normalised flux, for a caller that has Psi_N in
		/// hand. Exactly at( std::sqrt( normalisedFlux ), of ), and it exists so
		/// that the square root is written once.
		double atNormalisedFlux( double normalisedFlux,
		                         SurfaceQuantity const &of ) const;

		/// Whether @a label is inside the cut. A caller that would rather test
		/// than catch; at() refuses on exactly this predicate.
		bool covers( double label ) const;

		/// The named slots, one line each, so that a consumer does not write the
		/// lambda and a reader does not have to check that it wrote the right
		/// one. Every one of them is at() with a different callable.
		double vPrimeAt( double label ) const;
		double volumeAt( double label ) const;
		double inverseRSquaredAt( double label ) const;
		double gradPsiSquaredOverRSquaredAt( double label ) const;
		double arcLengthAt( double label ) const;
		double surfaceAreaAt( double label ) const;

		/// < | grad rho | > and < | grad rho |^2 >, which are what a 1-D
		/// transport metric on rho actually reads. The chain factor is the
		/// analytic radialDerivative() above and NOT a difference.
		double absGradLabelAt( double label ) const;
		double gradLabelSquaredAt( double label ) const;

		/// Whether any surface of the family crossed the band, and how many
		/// nodes over the whole family did. The per-node mask is on each
		/// surface; these are the summary and never a substitute for it.
		bool crossesBand() const;
		int extendedNodes() const;

		/// The worst | psi_h - level | over every node of every surface.
		double worstResidual() const;
	};

	/**
	 * The per-psi geometry cache: MANTA-COUPLING.md section 5's requirement and
	 * section 8's contract, in one object.
	 *
	 * It holds at most one family. A query whose psi vector is bitwise identical
	 * to the held one is served from it; anything else re-extracts. See the
	 * header for why the key is a bitwise comparison rather than a hash, why
	 * resetForRun() is not redundant with it, and why a failed extraction
	 * commits nothing.
	 *
	 * IT IS NOT THREAD SAFE AND IS NOT MEANT TO BE. The consumer's call pattern
	 * is a serial loop over physics nodes inside a residual evaluation; a cache
	 * shared between threads would need a lock around the extraction and would
	 * then serialise the only expensive part anyway.
	 */
	class GeometryCache
	{
		public:
			/// Extract a family from a psi vector. Throws when it cannot, which
			/// is propagated -- see the header's section 6.
			using Extractor =
				std::function<FluxSurfaceFamily( std::vector<double> const & )>;

			/// @param extractorIn what to call on a miss. Must not be empty.
			/// @throws std::invalid_argument if @a extractorIn is empty, since a
			///         cache that cannot fill itself would report every query as
			///         a miss and then fail with something unrelated.
			explicit GeometryCache( Extractor extractorIn );

			/**
			 * The family at @a psi, extracted if it is not the held one.
			 *
			 * @throws std::invalid_argument if @a psi is empty or holds a
			 *         non-finite entry. THE SECOND IS A GUARD AND NOT
			 *         FASTIDIOUSNESS: MANTA-COUPLING.md section 1 guarantees
			 *         meq is called at states far from equilibrium as a normal
			 *         part of the consumer's Newton iteration, and a NaN in psi
			 *         would otherwise reach the tracer and come back as a
			 *         geometry rather than as a refusal.
			 * @throws whatever the extractor throws, with nothing committed.
			 */
			FluxSurfaceFamily const &family( std::vector<double> const &psi );

			/// One pointwise query: the shape MANTA-COUPLING.md section 5
			/// describes. Exactly family( psi ).at( label, of ).
			double at( std::vector<double> const &psi, double label,
			           SurfaceQuantity const &of );

			/// Drop whatever is held.
			///
			/// CALL IT WHENEVER ANYTHING THE EXTRACTOR CLOSES OVER CHANGES --
			/// the mesh, the solver, the options. The key cannot see any of
			/// that, and MANTA-COUPLING.md section 8 records what a cache keyed
			/// on the object rather than the run cost the last time it was got
			/// wrong. Cheap and idempotent; call it on any doubt.
			void resetForRun();

			/// Whether a family is held. False after construction, after
			/// resetForRun(), and after an extraction that threw.
			bool held() const;

			/// Attempted extractions, queries served, and queries served from
			/// the held family. COUNTS AND NOT TIMINGS, deliberately: what the
			/// cache promises is that a family is extracted once per distinct
			/// psi, and that is a count. A timing here would be a measurement
			/// about the machine, which is what tests/performance is for.
			///
			/// extractions() counts ATTEMPTS, so an extraction that threw is in
			/// it: the work was done and the caller paid for it.
			std::size_t extractions() const;
			std::size_t queries() const;
			std::size_t hits() const;

		private:
			Extractor extractor;
			std::vector<double> key;
			FluxSurfaceFamily heldFamily;
			bool valid = false;
			std::size_t extractionCount = 0;
			std::size_t queryCount = 0;
			std::size_t hitCount = 0;
	};

	/// Psi_N from psi, and psi from Psi_N, with the convention stated once:
	/// Psi_N = ( psi_ax - psi ) / ( psi_ax - psi_bnd ), zero on the axis.
	///
	/// @throws std::invalid_argument if the span psi_ax - psi_bnd is zero, which
	///         is not a normalisation at all.
	double normalisedFlux( double psi, double psiAxis, double psiBoundary );
	double fluxAtNormalised( double normalisedFluxIn, double psiAxis,
	                         double psiBoundary );

}

#endif // MEQ_FLUXFAMILY_HPP
