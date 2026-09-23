#ifndef MEQ_CONDUCTORSTORE_HPP
#define MEQ_CONDUCTORSTORE_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

/**
 * `psi_c` ON DISK: A FILE-BACKED VERSION OF THE CONDUCTOR CACHES.
 *
 * `COIL-SUBTRACTION-PLAN.md` §1 splits the flux as `psi = psi_c + psi_p` and
 * licenses a cache of the conductors' own field, because the conductors do not
 * move during a forward solve. MEQ builds it ONCE PER MESH at three sets of
 * points -- `psi_c` at every potential dof and at every quadrature point
 * meq::SourceIntegrator visits, and **`q_c = ( 1/R ) grad_bar psi_c` at every
 * flux-space element node**, which meq::CriticalPointFinder needs because the
 * critical points are roots of the PHYSICAL flux. That precompute is the whole
 * reason the split pays for a machine of rectangles, where one field point
 * costs about 1e4 Carlson evaluations against 23 for the same machine as
 * filaments ( **[M-142](MEASUREMENTS.md#m-142)** ).
 *
 * **`q_c` IS THE EXPENSIVE ONE AND IT IS NOT `psi_c`.** A first version of this
 * file cached `psi_c` alone, measured 3.51x, and left the warm arm seventeen
 * times slower than the same machine as filaments -- both of them reading
 * `psi_c` from the same file. `perf` on the WARM run put the gradient kernels
 * at 61% of what remained. A gradient is a different quantity and no `psi_c`
 * cache covers it; **[M-170](MEASUREMENTS.md#m-170)** has the numbers and the
 * lesson.
 *
 * **THE COST IS PER MESH, AND A COUPLED RUN PAYS IT PER SOLVE.** Driven from
 * MaNTA the profiles move and the machine does not: `p'` and `g g'` change,
 * `psi_c` does not, and every one of those solves rebuilds a cache that could
 * not have changed. This is the file that lets the second solve read it.
 *
 * **IT IS A CACHE AND NOT A FIELD, WHICH DECIDES EVERY QUESTION ABOUT THE
 * FORMAT.** VMEC's `mgrid` samples the vacuum field on a uniform `( R, phi, Z )`
 * grid and interpolates it back, which makes it mesh-free and makes it
 * approximate. This is the other trade: the values are stored **at exactly the
 * points they are used at**, so a reload is bit-for-bit what the recompute
 * would have produced and `ConductorSubtraction`'s three identities keep
 * reading `0.000e+00`. Nothing here is interpolated and nothing here is
 * resampled. The price is that the file is **bound to one mesh, one degree and
 * one quadrature rule**, and is refused rather than adapted when any of them
 * moves.
 *
 * **SO THE SIGNATURE IS THE LOAD-BEARING PART.** A stale cache is not a slow
 * run, it is `psi_c` from a different machine added to a solve that converges;
 * the failure has no symptom. readConductorCache() therefore returns the
 * signature with the data and the caller must compare it, and
 * meq::GradShafranovSolver::adoptConductorCache() refuses on any mismatch and
 * says which field differed.
 *
 * **IT IS A NETCDF-4 GROUP EITHER WAY**, which is the one thing the two
 * deployments need in common: written as the root group of a file of its own
 * for a cache that outlives a run, or added as a named group to the
 * equilibrium `.nc` so that one file carries the answer and the means to warm
 * start from it. MEQ's `.nc` is already netCDF-4, so embedding costs no format
 * change; see `tools/README.md` for which reader takes which file.
 */
namespace meq
{
	/// The identity a cached `psi_c` is valid for, and nothing else. Every
	/// field is cheap to compute without building the cache, which is the
	/// point: a caller decides whether to trust a file before paying for the
	/// thing the file would save.
	struct ConductorCacheSignature
	{
		/// The mesh's vertices, connectivity and attributes. Positions decide
		/// where `psi_c` was evaluated, so a mesh that moved by one vertex is a
		/// different cache.
		std::uint64_t meshDigest = 0;

		/// Every conductor's geometry and current, the permeability, and the
		/// cross-section quadrature order -- everything `psi_c` is a function
		/// of besides the field point.
		std::uint64_t conductorDigest = 0;

		/// The potential space's degree. The dof positions are the element's
		/// `GaussLobatto` nodes, so this changes them.
		int polynomialDegree = -1;

		int elements = 0;
		int potentialDofs = 0;

		bool operator==( ConductorCacheSignature const &other ) const;
		bool operator!=( ConductorCacheSignature const &other ) const;

		/// The first field that differs, as a sentence naming both values, or
		/// an empty string when they agree. This is what a refusal prints: a
		/// bare "the cache does not match" sends the reader to guess which of
		/// five things moved.
		std::string difference( ConductorCacheSignature const &other ) const;
	};

	/// `psi_c` at the two sets of points the solver caches it at. Both may be
	/// empty -- a solve with no conductor field caches nothing -- and a file
	/// may legitimately carry only the nodal half, since that is the half a
	/// warm start's stored guess needs.
	struct ConductorCache
	{
		ConductorCacheSignature signature;

		/// `psi_c` at every potential dof, indexed by dof. The potential space
		/// is L2, so each dof belongs to exactly one element.
		std::vector< double > nodalPsi;

		/// `psi_c` at every source quadrature point, flattened, with
		/// `quadratureOffset[ e ]` the first entry of element `e`.
		std::vector< double > quadraturePsi;

		/// `elements + 1` entries, a prefix sum. Stored rather than recomputed
		/// so a reader that is not MEQ can index the flattened array, and so a
		/// reload can check the rule it would have used against the rule that
		/// was used.
		std::vector< int > quadratureOffset;

		/*
		 * AND `q_c`, WHICH IS THE EXPENSIVE HALF AND IS NOT `psi_c` AT ALL.
		 *
		 * meq::CriticalPointFinder keeps its own table over the element nodes
		 * of the potential and flux spaces, because the critical points are
		 * roots of the PHYSICAL flux and this solver holds the remainder. The
		 * flux half is `q_c = ( 1/R ) grad_bar psi_c`, and measured on a
		 * machine of rectangles it dominates a warm solve outright: caching
		 * `psi_c` alone left `filamentGradKernel` and `ellint_rd` at 61% of the
		 * profile, because a gradient is a different quantity and the solver's
		 * caches never held it.
		 *
		 * It is indexed by ( element, local node ) rather than by dof, which is
		 * why it cannot share `nodalPsi` above: the finder walks elements and
		 * MEQ does not require the two spaces to share a node set.
		 */
		std::vector< double > criticalPotentialPsi;
		std::vector< int > criticalPotentialOffset;

		/// `( qR, qZ )` interleaved, per flux-space element node.
		std::vector< double > criticalFluxQ;

		/// One byte per flux node. `q_c` is NaN on the axis -- that is
		/// meq::ConductorField::flux()'s contract -- so the screen travels with
		/// the values rather than being re-derived from them.
		std::vector< unsigned char > criticalFluxUsable;
		std::vector< int > criticalFluxOffset;

		/// `max |q_c|` over the flux nodes, reduced at build time.
		double criticalFluxScale = 0.0;

		/// Whether the critical-point half is present. It travels as a whole:
		/// the offsets index the values and the screen indexes the nodes, so a
		/// partial one indexes nothing.
		bool hasCriticalPointTable() const;

		bool empty() const;

		/// Whether the two arrays agree with the signature and with each other
		/// -- offsets monotone, the last offset the flattened size, the nodal
		/// count the dof count. A file is checked on the way in, because
		/// everything downstream indexes these without bounds checks.
		bool consistent() const;

		private:
			/// consistent()'s critical-point half, split out because it has
			/// four arrays that must agree with each other and with the
			/// element count, and because "absent is fine, partial is not" is
			/// a different rule from the one the quadrature half follows.
			bool criticalPointTableIsConsistent() const;
	};

	/// FNV-1a over the IEEE bytes. Exposed because the mesh and conductor
	/// digests are built by callers that own those objects -- the solver knows
	/// its mesh, meq::ConductorField knows its conductors -- and a digest is
	/// only comparable if both sides mix in the same order.
	std::uint64_t digestSeed();
	std::uint64_t digestAppend( std::uint64_t digest, double value );
	std::uint64_t digestAppend( std::uint64_t digest, std::int64_t value );

	/// The group name MEQ writes when embedding, and looks for when reading a
	/// file it was not told the group of.
	char const *conductorCacheGroupName();

	/// Write `cache`.
	///
	/// With an empty `group` this REPLACES `path` with a new netCDF-4 file
	/// whose root group is the cache. With a group name it OPENS `path` for
	/// writing and adds the group, which is how the equilibrium `.nc` comes to
	/// carry it -- so the equilibrium must already have been written.
	///
	/// @throws std::invalid_argument if the cache is not consistent(), because
	///         writing a cache that will be refused on the way back in is a
	///         silent waste rather than an error anyone sees.
	/// @throws std::runtime_error naming `path` if netCDF refuses.
	void writeConductorCache( std::string const &path,
	                          ConductorCache const &cache,
	                          std::string const &group = std::string() );

	/// Read a cache from `path`.
	///
	/// With an empty `group` this takes the root group if it carries the
	/// format attribute and otherwise looks for conductorCacheGroupName(), so
	/// one call reads both deployments. Returns nothing -- rather than
	/// throwing -- when the file exists but carries no cache, since "this `.nc`
	/// has no cache in it" is an ordinary answer and the caller's next move is
	/// to compute one.
	///
	/// @throws std::runtime_error naming `path` if the file cannot be opened,
	///         or if a cache is present but malformed. A malformed cache IS an
	///         error: it is a file claiming to be one of these.
	std::optional< ConductorCache >
	readConductorCache( std::string const &path,
	                    std::string const &group = std::string() );
}

#endif
