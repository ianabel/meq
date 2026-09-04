#ifndef MEQ_FIELD_HPP
#define MEQ_FIELD_HPP

/*
 * The poloidal magnetic field, from the solved flux.
 *
 * This file is four lines of arithmetic and a page of justification, which is
 * the right ratio: the arithmetic is a relabelling and the justification is the
 * whole argument for the discretisation.
 *
 * WHERE B COMES FROM. refs/Miller.pdf eq (1) writes the field as
 *
 *     B = grad(phi) x grad(psi) + f(psi) grad(phi),
 *
 * phi the toroidal angle, so grad(phi) = e_phi / r. With
 * grad(psi) = ( d_r psi ) e_r + ( d_z psi ) e_z and
 * e_phi x e_r = e_z, e_phi x e_z = -e_r,
 *
 *     B_pol = ( 1/r )[ ( d_r psi ) e_z - ( d_z psi ) e_r ],
 *
 * that is
 *
 *     B_R = -( 1/r ) d_z psi,      B_Z = +( 1/r ) d_r psi.
 *
 * AND MEQ ALREADY SOLVES FOR THAT. The HDG flux variable is
 * q = ( 1/r ) grad_bar( psi ) = ( q_r, q_z ), so
 *
 *     B_R = -q_z,      B_Z = +q_r
 *
 * pointwise, with NO DERIVATIVE TAKEN ANYWHERE. That is the whole reason the
 * mixed method is worth its extra unknowns: README.md claims the physically
 * interesting output is the field and that a mixed method resolves it at the
 * same order as psi rather than one lower, and this is where that is cashed in.
 * q converges at k+1, so B does.
 *
 * A DIFFERENTIATED psi WOULD LOSE AN ORDER, and worse than that on this problem:
 * tests/convergence/EstimatorConvergence.cpp measures the estimator's eta_2 both
 * ways and building it on a differentiated potential costs exactly one order at
 * every k. The same would apply to a field obtained that way.
 *
 * THE SIGNS ARE MEASURED, NOT ARGUED. The derivation above is a derivation, and
 * of this project's four convention questions three went against the
 * derivation -- the Solov'ev source sign, tau in eq (8e), and DarcyForm holding
 * -q. tests/convergence/FieldConvergence.cpp checks B against central
 * differences of an EXACT psi and against the fixture's analytic flux, because a
 * field pointing the wrong way is invisible to every convergence rate in the
 * suite.
 *
 * NOTE ON WHICH q. GradShafranovSolver::flux() has already undone DarcyForm's
 * sign convention -- the block vector holds -q, see the file comment there --
 * so the mapping above applies to flux() and NOT to the raw block. Getting that
 * wrong flips the field, which is exactly what the test exists to catch.
 */

#include "mfem.hpp"

namespace meq
{
	/// B_pol from the HDG flux q, as a vdim-2 GridFunction on q's own space:
	/// B_R = -q_z, B_Z = +q_r. @a q must be what
	/// GradShafranovSolver::flux() returns, in MEQ's sign convention.
	///
	/// @a field is resized and given q's space if it does not have one. The
	/// spaces must match if it does.
	void poloidalField( mfem::GridFunction const &q, mfem::GridFunction &field );

	/// The same relabelling as a Coefficient, for sampling B without building a
	/// GridFunction for it. @a component is 0 for B_R and 1 for B_Z.
	class PoloidalFieldCoefficient : public mfem::Coefficient
	{
		public:
			PoloidalFieldCoefficient( mfem::GridFunction const &qIn, int componentIn );

			/// MFEM's spelling, from mfem::Coefficient.
			double Eval( mfem::ElementTransformation &tr, // NOLINT(readability-identifier-naming)
			             mfem::IntegrationPoint const &ip ) override;

		private:
			mfem::GridFunction const &q;
			int component;
			mfem::Vector value;
	};
}

#endif // MEQ_FIELD_HPP
