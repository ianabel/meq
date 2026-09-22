#ifndef MEQ_FIELDVIEWS_HPP
#define MEQ_FIELDVIEWS_HPP

/*
 * EVALUATING THE PHYSICAL FIELD WHEN THE SOLVED ONE IS ONLY PART OF IT.
 *
 * Under `COIL-SUBTRACTION-PLAN.md`'s split the solver holds `psi_p` and the
 * physical flux is `psi_c + psi_p`. `COIL-SUBTRACTION-PLAN.md` §8.1 counts
 * about sixty-eight places that read the solved field, and §3 names the hazard
 * plainly: every consumer that forgets is a silent wrong answer rather than a
 * failure. This file exists so that the forgetting is hard.
 *
 * WHAT MAKES IT HARD IS THE NAMED CONSTRUCTOR AND NOT THE DOCUMENTATION.
 * There is no default constructor and no implicit conversion from a
 * `GridFunction`; a caller writes either
 *
 *     PotentialView::physical( solver )     psi_c + psi_p, the real flux
 *     PotentialView::remainder( solver )    psi_p alone, deliberately
 *
 * so the choice is made in a word at the call site rather than by knowing what
 * this class does. `remainder()` is spelled out because there ARE consumers
 * that want it -- the residual estimator measures the discretisation of the
 * problem actually solved -- and a reviewer can grep for it.
 *
 * AND IT IS FREE WHEN THE SPLIT IS NOT IN USE, WHICH IS WHY IT IS A HEADER.
 * Every method is inline and the conductor pointer is null on every existing
 * path, so what the compiler sees is the same `GetValue()` call it saw before
 * behind one predictable branch. `carriesConductors()` is false, the addition
 * never happens, and nothing that does not use the split moves by a bit --
 * which is the property the acceptance asserts rather than assumes.
 *
 * THE COST WHEN IT IS IN USE IS REAL AND IS NOT HIDDEN. `psi_c` is a Carlson
 * elliptic integral per filament and a cross-section quadrature per rectangle,
 * evaluated at every point a consumer asks about. `COIL-SUBTRACTION-PLAN.md`
 * §1 says `psi_c` should be computed once per mesh rather than once per
 * evaluation, and **this class is where that cache belongs when it is built**:
 * one type to change rather than sixty-eight call sites.
 *
 * A TRAP THIS FILE HAS TO AVOID AND EVERY CONSUMER WOULD MEET SEPARATELY.
 * Adding `psi_c` needs the point's `( R, z )`, which needs the element's
 * transformation -- and `mfem::Mesh::GetElementTransformation( int )` hands out
 * SHARED SCRATCH, which `CLAUDE.md` records as a silent wrong answer under
 * threading and which cost this project six call sites once already. The
 * overloads below take a caller-supplied transformation where one is in hand,
 * and use a function-local `thread_local` where it is not.
 */

#include "mfem.hpp"

#include "ConductorField.hpp"

namespace meq
{
	namespace detail
	{
		/// The physical point of an integration point, through the REENTRANT
		/// transformation overload. See the trap note at the top of this file.
		inline void viewPoint( mfem::Mesh &mesh, int element,
		                       mfem::IntegrationPoint const &ip,
		                       double &radius, double &z )
		{
			// Function-local rather than a member, so a transformation held
			// live across a call cannot be reset underneath it -- the same fix
			// the six sites in FluxSurfaces and CriticalPoints took.
			thread_local mfem::IsoparametricTransformation transformation;
			mesh.GetElementTransformation( element, &transformation );
			transformation.SetIntPoint( &ip );

			double coordinates[ 3 ] = { 0.0, 0.0, 0.0 };
			mfem::Vector position( coordinates, 3 );
			transformation.Transform( ip, position );
			radius = position( 0 );
			z = position( 1 );
		}

		inline void viewPoint( mfem::ElementTransformation &transformation,
		                       mfem::IntegrationPoint const &ip,
		                       double &radius, double &z )
		{
			double coordinates[ 3 ] = { 0.0, 0.0, 0.0 };
			mfem::Vector position( coordinates, 3 );
			transformation.Transform( ip, position );
			radius = position( 0 );
			z = position( 1 );
		}
	}

	/**
	 * `psi` at a point: the solved potential, plus `psi_c` when the split is in
	 * use. Built by name — `physical()` or `remainder()` — and by nothing else.
	 */
	class PotentialView
	{
		public:
			/// The PHYSICAL flux, `psi_c + psi_p`. Pass null conductors and it
			/// is the solved field, which is what every existing path is.
			static PotentialView physical( mfem::GridFunction const &solved,
			                               ConductorField const *conductors )
			{
				return PotentialView( solved, conductors );
			}

			/// The SOLVED field alone, deliberately. For consumers that mean
			/// the discretisation of the problem actually solved rather than
			/// the physical equilibrium — the residual estimator is the one
			/// that does. Greppable on purpose.
			static PotentialView remainder( mfem::GridFunction const &solved )
			{
				return PotentialView( solved, nullptr );
			}

			bool carriesConductors() const
			{
				return conductorField != nullptr;
			}

			/// @param element  the mesh element index the point lies in.
			double value( int element,
			              mfem::IntegrationPoint const &ip ) const
			{
				double const solvedValue = field->GetValue( element, ip );
				if ( !conductorField )
					return solvedValue;

				double radius = 0.0;
				double z = 0.0;
				detail::viewPoint( *field->FESpace()->GetMesh(), element, ip,
				                   radius, z );
				return solvedValue + conductorField->psi( radius, z );
			}

			/// The overload for a caller that already holds the element's
			/// transformation, which is every assembly loop. Cheaper, and it
			/// cannot meet the shared-scratch trap at all.
			double value( mfem::ElementTransformation &transformation,
			              mfem::IntegrationPoint const &ip ) const
			{
				double const solvedValue = field->GetValue( transformation,
				                                            ip );
				if ( !conductorField )
					return solvedValue;

				double radius = 0.0;
				double z = 0.0;
				detail::viewPoint( transformation, ip, radius, z );
				return solvedValue + conductorField->psi( radius, z );
			}

			/// The solved field itself, for the handful of places that must
			/// pass a GridFunction on to something else. **Reading this is
			/// reading the REMAINDER**; it is named so that a call site says so.
			mfem::GridFunction const &solvedField() const
			{
				return *field;
			}

		private:
			PotentialView( mfem::GridFunction const &solvedIn,
			               ConductorField const *conductorsIn )
				: field( &solvedIn ), conductorField( conductorsIn )
			{
			}

			mfem::GridFunction const *field;
			ConductorField const *conductorField;
	};

	/**
	 * `q = ( 1/R ) grad_bar( psi )` at a point, on the same contract.
	 *
	 * **`q_c` ADDS, BECAUSE `Δ*` IS LINEAR AND `q` IS A RELABELLING.** MEQ
	 * solves for `q` rather than differentiating `psi`, so the physical flux is
	 * `q_c + q_p` with `q_c` from meq::ConductorField::flux() — and `B` follows
	 * by `B_R = −q_z`, `B_Z = +q_r`, which is why meq::Field needs nothing of
	 * its own here.
	 */
	class FluxVectorView
	{
		public:
			static FluxVectorView physical( mfem::GridFunction const &solved,
			                                ConductorField const *conductors )
			{
				return FluxVectorView( solved, conductors );
			}

			static FluxVectorView remainder( mfem::GridFunction const &solved )
			{
				return FluxVectorView( solved, nullptr );
			}

			bool carriesConductors() const
			{
				return conductorField != nullptr;
			}

			/// @note NaN on the axis when conductors are carried, exactly as
			///       meq::filamentFlux() is and for the same reason: `q` is
			///       `psi`'s gradient over `R`.
			void value( int element, mfem::IntegrationPoint const &ip,
			            mfem::Vector &out ) const
			{
				field->GetVectorValue( element, ip, out );
				if ( !conductorField )
					return;

				double radius = 0.0;
				double z = 0.0;
				detail::viewPoint( *field->FESpace()->GetMesh(), element, ip,
				                   radius, z );

				double qR = 0.0;
				double qZ = 0.0;
				conductorField->flux( radius, z, qR, qZ );
				out( 0 ) += qR;
				out( 1 ) += qZ;
			}

			mfem::GridFunction const &solvedField() const
			{
				return *field;
			}

		private:
			FluxVectorView( mfem::GridFunction const &solvedIn,
			                ConductorField const *conductorsIn )
				: field( &solvedIn ), conductorField( conductorsIn )
			{
			}

			mfem::GridFunction const *field;
			ConductorField const *conductorField;
	};
}

#endif
