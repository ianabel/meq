#include "SourceView.hpp"

#include "Coils.hpp"
#include "Profiles.hpp"

#include <cmath>

namespace meq
{

	namespace
	{

		/// A meq::Profile as a ProfileView, or false.
		///
		/// meq::ScaledProfile is folded into `scale` rather than given a kind,
		/// and the fold is RECURSIVE because nothing stops a caller nesting
		/// them -- makeNormalisedSource wraps a profile once, the I_p border
		/// can wrap the result again, and a view built from only the outermost
		/// factor would be quietly off by the inner one.
		bool viewOfProfile( Profile const &profile, ProfileView &view,
		                    std::vector<KnotView> &knots,
		                    mfem::real_t scale = mfem::real_t( 1 ) )
		{
			if ( auto const *scaled = dynamic_cast<ScaledProfile const *>( &profile ) )
				return viewOfProfile( scaled->unscaled(), view, knots,
				                      scale*static_cast<mfem::real_t>( scaled->scale() ) );

			if ( auto const *constant = dynamic_cast<ConstantProfile const *>( &profile ) )
			{
				view.kind = static_cast<int>( ProfileKind::Constant );
				view.scale = scale;
				view.constantValue = static_cast<mfem::real_t>( constant->value() );
				view.offset = 0;
				view.count = 0;
				return true;
			}

			if ( auto const *spline = dynamic_cast<SplineProfile const *>( &profile ) )
			{
				std::vector<Knot> const &data = spline->knots();
				// A Hermite cubic needs two knots to be a cubic at all, and
				// evaluateProfile() indexes [ i + 1 ] unconditionally.
				if ( data.size() < 2 )
					return false;

				view.kind = static_cast<int>( ProfileKind::Spline );
				view.scale = scale;
				view.constantValue = mfem::real_t( 0 );
				view.offset = static_cast<int>( knots.size() );
				view.count = static_cast<int>( data.size() );
				for ( Knot const &k : data )
					knots.push_back( KnotView{
						static_cast<mfem::real_t>( k.psi ),
						static_cast<mfem::real_t>( k.value ),
						static_cast<mfem::real_t>( k.derivative ) } );

				/*
				 * UNIFORMITY IS MEASURED FROM THE STORED ABSCISSAE, NOT
				 * ASSUMED FROM THE FILE. A profile written on a uniform grid
				 * and one written on a nearly-uniform one look identical to a
				 * reader and behave differently in the kernel, so the test is
				 * on the data: every spacing within a relative 1e-9 of the
				 * first. That tolerance is far looser than the ulp-level
				 * wobble of decimal abscissae read back as binary, and far
				 * tighter than any grading a profile would deliberately carry.
				 *
				 * AND GETTING IT WRONG EITHER WAY IS SAFE, which is why a
				 * tolerance is acceptable here at all: a uniform profile
				 * misjudged non-uniform is merely slower, and a non-uniform
				 * one misjudged uniform is CORRECTED back by
				 * profileInterval()'s two comparisons. The flag is a
				 * performance hint whose failure mode is performance.
				 */
				view.uniform = 0;
				view.firstKnot = static_cast<mfem::real_t>( data.front().psi );
				view.inverseSpacing = mfem::real_t( 1 );
				if ( data.size() >= 2 )
				{
					double const h = data[ 1 ].psi - data[ 0 ].psi;
					bool uniform = h > 0.0;
					for ( std::size_t i = 1; uniform && i + 1 < data.size(); ++i )
					{
						double const step = data[ i + 1 ].psi - data[ i ].psi;
						if ( std::abs( step - h ) > 1.0e-9*std::abs( h ) )
							uniform = false;
					}
					if ( uniform )
					{
						view.uniform = 1;
						view.inverseSpacing = static_cast<mfem::real_t>( 1.0/h );
					}
				}
				return true;
			}

			// Anything else -- a profile kind added later, or one read through
			// a wrapper this function does not know -- has no POD form. Refuse.
			return false;
		}

	}

	bool makeSourceView( Source const &source, SourceView &view,
	                     std::vector<KnotView> &knots )
	{
		view = SourceView{};

		if ( auto const *soloviev = dynamic_cast<SolovievSource const *>( &source ) )
		{
			view.kind = static_cast<int>( SourceKind::Soloviev );
			view.solovievA = static_cast<mfem::real_t>( soloviev->a() );
			return true;
		}

		auto const *mhd = dynamic_cast<NormalisedMHDSource const *>( &source );
		if ( mhd == nullptr )
			return false;

		/*
		 * THE WRAPPERS ARE REFUSED BY THE dynamic_cast ABOVE AND THAT IS
		 * DELIBERATE RATHER THAN INCIDENTAL. A coil-augmented source derives
		 * from NormalisedSource and adds a term this view does not carry, and
		 * meq::RotatingSource closes its species balance by a root find -- so
		 * a cast that accepted a derived class would build a view describing
		 * the plasma term alone and silently drop the rest. The cast is to the
		 * EXACT production class for that reason; a new kind is a new
		 * SourceKind, never a widened cast.
		 */

		view.kind = static_cast<int>( SourceKind::NormalisedMHD );
		view.permeability = static_cast<mfem::real_t>( mhd->mu0() );
		view.psiAxis = static_cast<mfem::real_t>( mhd->normalisation() );
		view.psiBoundary = static_cast<mfem::real_t>( mhd->boundaryNormalisation() );
		view.currentScale = static_cast<mfem::real_t>( mhd->currentScale() );
		view.confineToPlasma = mhd->plasmaSupport() ? 1 : 0;
		view.supportAxis = static_cast<mfem::real_t>( mhd->supportAxis() );
		view.supportBoundary = static_cast<mfem::real_t>( mhd->supportBoundary() );

		return viewOfProfile( mhd->pPrime(), view.pPrime, knots )
		       && viewOfProfile( mhd->ggPrime(), view.ggPrime, knots );
	}

	bool makeCoilViews( Source const &source, std::vector<CoilView> &coils,
	                    mfem::real_t &permeability )
	{
		coils.clear();
		permeability = mfem::real_t( 0 );

		/*
		 * ASKED OF meq::Source RATHER THAN OF A WRAPPER TYPE. `conductors()`
		 * is the virtual every source answers -- null when it has none -- so
		 * this needs no cast and cannot go stale against a new wrapper, which
		 * is the failure makeSourceView()'s exact-class cast guards against
		 * from the other side.
		 */
		CoilSet const *const set = source.conductors();
		if ( set == nullptr )
			return true;   // no conductors: an empty list, not a refusal

		permeability = static_cast<mfem::real_t>( set->mu0() );
		coils.reserve( set->size() );
		for ( std::size_t i = 0; i < set->size(); ++i )
		{
			Coil const &c = set->coil( i );
			coils.push_back( CoilView{
				static_cast<mfem::real_t>( c.minRadius() ),
				static_cast<mfem::real_t>( c.maxRadius() ),
				static_cast<mfem::real_t>( c.zMin() ),
				static_cast<mfem::real_t>( c.zMax() ),
				static_cast<mfem::real_t>( c.currentDensity() ) } );
		}
		return true;
	}

	bool makePlasmaSourceView( Source const &source, SourceView &view,
	                           std::vector<KnotView> &knots )
	{
		/*
		 * UNWRAP, THEN DELEGATE. Both wrappers hold the plasma source and both
		 * expose it as plasma(); the normalised one additionally overrides
		 * scaledF() to return it, which is the function this view is a POD of.
		 * Delegating rather than duplicating is what keeps the ONE expression
		 * of the profile walk -- viewOfProfile() -- and it means a new
		 * SourceKind is added in exactly one place.
		 *
		 * The normalised wrapper is checked FIRST because it derives from
		 * NormalisedSource and the other does not, so the two casts are
		 * disjoint; the order is for a reader rather than for the compiler.
		 */
		if ( auto const *normalised =
		         dynamic_cast<CoilAugmentedNormalisedSource const *>( &source ) )
			return makeSourceView( normalised->plasma(), view, knots );

		if ( auto const *augmented =
		         dynamic_cast<CoilAugmentedSource const *>( &source ) )
			return makeSourceView( augmented->plasma(), view, knots );

		// Not wrapped: scaledF() is f() by meq::Source's default, so the
		// ordinary view already describes it and there is nothing to strip.
		return makeSourceView( source, view, knots );
	}

}
