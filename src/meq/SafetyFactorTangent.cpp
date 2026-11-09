#include "meq/SafetyFactorTangent.hpp"

#include <chrono>
#include <cmath>
#include <stdexcept>

#include "meq/ContinuousField.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/ResidentIterate.hpp"

namespace meq
{

	namespace
	{
		double secondsSince( std::chrono::steady_clock::time_point since )
		{
			return std::chrono::duration<double>( std::chrono::steady_clock::now() - since ).count();
		}

		/*
		 * A STRETCH OF A VECTOR, OR ITS DOT WITH ANOTHER, ACCUMULATED IN DOUBLE:
		 * on the device where one is configured and the kernels are built
		 * (meq::doubleDotKernel()), else on the host through HostRead(). What
		 * the pin rows' quadratures are contracted with, so that in single a
		 * sum over every potential dof does not round in real_t.
		 */
		double sumInDouble( mfem::Vector const &a, int aBegin, mfem::Vector const *b, int bBegin, int count )
		{
			double out = 0.0;
			if ( mfem::Device::Allows( mfem::Backend::DEVICE_MASK )
			     && doubleDotKernel( a, aBegin, b, bBegin, count, out ) )
				return out;
			mfem::real_t const *x = a.HostRead();
			mfem::real_t const *y = b != nullptr ? b->HostRead() : nullptr;
			out = 0.0;
			for ( int k = 0; k < count; ++k )
				out += static_cast<double>( x[ aBegin + k ] )
				       *( y != nullptr ? static_cast<double>( y[ bBegin + k ] ) : 1.0 );
			return out;
		}
	}

	double plasmaCells( mfem::Mesh const &mesh )
	{
		double cells = 0.0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
			cells += mesh.GetElementGeometry( e ) == mfem::Geometry::TRIANGLE ? 0.5 : 1.0;
		return cells;
	}

	SafetyFactorEquilibrium::SafetyFactorEquilibrium( GradShafranovSolver &solverIn,
	                                                  NormalisedMHDSource &sourceIn, Options optionsIn )
		: solver( solverIn ), source( sourceIn ), options( std::move( optionsIn ) )
	{
		if ( !( options.step > 0.0 ) )
			throw std::invalid_argument( "meq::SafetyFactorEquilibrium: the difference step must be positive" );
		if ( ( options.gammaHMarker == nullptr ) != ( options.path == nullptr ) )
			throw std::invalid_argument( "meq::SafetyFactorEquilibrium: the band needs both its "
			                             "Gamma_h marker and its transfer path" );
	}

	std::vector<double> SafetyFactorEquilibrium::labels() const
	{
		std::vector<double> out = fluxSurfaceLabels( options.family );
		for ( double &rho : out )
			rho *= rho;
		return out;
	}

	SafetyFactorEquilibrium::Fields SafetyFactorEquilibrium::solvedFields() const
	{
		Fields f;
		f.psiStar = std::make_unique<mfem::GridFunction>( solver.postProcessedPotential() );
		f.qStar = std::make_unique<mfem::GridFunction>( solver.postProcessedFlux() );
		f.trace = std::make_unique<mfem::GridFunction>( solver.trace() );
		f.psi = std::make_unique<mfem::GridFunction>( solver.potential() );
		f.flux = std::make_unique<mfem::GridFunction>( solver.flux() );
		return f;
	}

	/*
	 * THE TANGENT'S FIELDS. The state is ( q, psi, psihat ) in DarcyForm's
	 * sign, which holds -q: flux() is its negation, so the tangent's flux is
	 * too. psi* and q* by reconstructFields(), which is linear in the state.
	 */
	SafetyFactorEquilibrium::Fields SafetyFactorEquilibrium::along( mfem::Vector const &dx ) const
	{
		Fields f = solvedFields();
		solver.reconstructFields( dx, *f.psiStar, *f.qStar );
		int const fluxSize = f.flux->Size();
		int const potentialSize = f.psi->Size();
		// Views of dx's blocks, copied by mfem::Vector's own device-aware
		// operations: no host loop over the dofs.
		mfem::Vector &blocks = const_cast<mfem::Vector &>( dx );
		mfem::Vector const fluxBlock( blocks, 0, fluxSize );
		mfem::Vector const potentialBlock( blocks, fluxSize, potentialSize );
		mfem::Vector const traceBlock( blocks, fluxSize + potentialSize, f.trace->Size() );
		f.flux->Set( -1.0, fluxBlock );
		f.psi->Set( 1.0, potentialBlock );
		f.trace->Set( 1.0, traceBlock );
		return f;
	}

	SafetyFactorEquilibrium::Fields SafetyFactorEquilibrium::combine( Fields const &base,
	                                                                  Fields const &direction, double s )
	{
		Fields f;
		auto sum = []( mfem::GridFunction const &a, mfem::GridFunction const &b, double t )
		{
			auto out = std::make_unique<mfem::GridFunction>( a );
			out->Add( t, b );
			return out;
		};
		f.psiStar = sum( *base.psiStar, *direction.psiStar, s );
		f.qStar = sum( *base.qStar, *direction.qStar, s );
		f.trace = sum( *base.trace, *direction.trace, s );
		f.psi = sum( *base.psi, *direction.psi, s );
		f.flux = sum( *base.flux, *direction.flux, s );
		return f;
	}

	/*
	 * q/g AT THE KNOTS OF A FIELD: the family on the continuous pair, rooted
	 * from the axis the raw fields locate -- CriticalPointFinder over q_h and
	 * psi_h, as the write phase's family is -- and labelled by @a psiAxis,
	 * the solve's own psi_ax (meq::labelledAxis()).
	 */
	CriticalPoint SafetyFactorEquilibrium::locatedAxis( Fields const &fields ) const
	{
		CriticalPointFinder const finder( *fields.flux, *fields.psi );
		return finder.findAxis();
	}

	std::vector<double> SafetyFactorEquilibrium::qOverG( Fields const &fields, double psiAxis ) const
	{
		CriticalPoint axis = locatedAxis( fields );
		axis.psi = psiAxis;
		return qOverG( fields, axis );
	}

	std::vector<double> SafetyFactorEquilibrium::qOverG( Fields const &fields, CriticalPoint const &axis ) const
	{
		auto const started = std::chrono::steady_clock::now();
		ContinuousFields const continuous =
			buildContinuousFields( *fields.psiStar, *fields.qStar, *fields.trace );
		ContourTracer tracer( *continuous.potential, *continuous.flux );
		if ( options.gammaHMarker != nullptr )
			tracer.setBandExtension( BandExtension::TransferLift, *options.gammaHMarker, options.path );
		FluxSurfaceFamily const family = extractFluxSurfaces( tracer, axis, solver.psiBoundary(),
		                                                      options.family );
		std::vector<double> out;
		out.reserve( family.surfaces.size() );
		for ( FluxSurface const &surface : family.surfaces )
			out.push_back( surface.vPrime*surface.inverseRSquared/( 4.0*M_PI*M_PI ) );
		secondsExtracting += secondsSince( started );
		extractionCount += 1;
		return out;
	}

	double SafetyFactorEquilibrium::measure( ToroidalFieldSpline const &spline ) const
	{
		switch ( options.pin )
		{
			case ScalePin::EdgeField:
				return std::sqrt( std::max( spline( 0.0 ), 0.0 ) );
			case ScalePin::PlasmaCurrent:
				return solver.deliveredPlasmaCurrent();
			case ScalePin::ToroidalFlux:
				break;
		}
		double const boundary = solver.psiBoundary();
		double const span = solver.psiAxis() - boundary;
		return options.sign*solver.integrateOverDomain(
			[ & ]( double radius, double, double psi )
			{
				if ( !( radius > 0.0 ) || span == 0.0 )
					return 0.0;
				double const normalised = ( psi - boundary )/span;
				if ( !( normalised > 0.0 ) )
					return 0.0;
				double const value = spline( normalised );
				return value > 0.0 ? std::sqrt( value )/radius : 0.0;
			} );
	}

	bool SafetyFactorEquilibrium::evaluate( ToroidalFieldSpline const &spline, SafetyFactorEvaluation &out )
	{
		auto const started = std::chrono::steady_clock::now();
		try
		{
			source.setGGPrime( spline.ggPrime() );
			if ( warm.Size() > 0 )
				solver.setInitialState( warm, {} );
			if ( options.solve )
				options.solve();
			else
				solver.solve();
			newtonIterations += solver.newtonIterations();
			solver.postProcess();
			secondsSolving += secondsSince( started );
			out.qOverG = qOverG( solvedFields(), solver.psiAxis() );
			out.measured = measure( spline );
		}
		catch ( std::exception const & )
		{
			secondsSolving += secondsSince( started );
			admittedValues.clear();
			return false;
		}
		warm = solver.state();
		admitted = warm;
		admittedValues = spline.values();
		return true;
	}

	bool SafetyFactorEquilibrium::differentiate( ToroidalFieldSpline const &spline,
	                                             SafetyFactorDerivative &out )
	{
		if ( admittedValues.empty() || spline.values() != admittedValues )
			throw std::logic_error( "meq::SafetyFactorEquilibrium::differentiate: the solver does "
			                        "not hold the state of this spline" );
		auto started = std::chrono::steady_clock::now();
		std::size_t const m = spline.size();
		mfem::Vector const state( solver.state() );
		solver.linearise( state );

		std::vector<mfem::Vector> dx( m ), sources( m );
		std::vector<double> dAxis( m, 0.0 );
		std::vector<std::shared_ptr<Profile const>> directions( m );
		for ( std::size_t j = 0; j < m; ++j )
		{
			directions[ j ] = spline.cardinalGGPrime( j );
			Profile const &direction = *directions[ j ];
			mfem::Vector &dRdu = sources[ j ];
			solver.assembleSourceDerivative( state,
				[ this, &direction ]( double radius, double z, double psi )
				{
					return source.ggPrimeDerivative( radius, z, psi, direction );
				}, dRdu );
			std::vector<double> ds;
			solver.linearisedTangent( dRdu, {}, nullptr, dx[ j ], ds );
			dAxis[ j ] = ds.empty() ? 0.0 : ds.front();
		}
		secondsTangents += secondsSince( started );
		lastAxisDerivative = dAxis;

		Fields const base = solvedFields();
		double const scale = base.psiStar->Normlinf();
		double const psiAxis = solver.psiAxis();
		CriticalPoint const axis = locatedAxis( base );

		/*
		 * dq_h/dx at the axis, in its own element: the root's derivative is
		 * -( dq_h/dx )^-1 dq along a tangent. The point by an unrestricted
		 * Newton inverse -- a root a little outside its element is the
		 * ordinary case on a face (CriticalPoint::outside), and the default
		 * solver would clamp it onto the element (CLAUDE.md, MFEM's API).
		 */
		mfem::Mesh &mesh = *base.flux->FESpace()->GetMesh();
		mfem::IsoparametricTransformation transformation;
		mesh.GetElementTransformation( axis.element, &transformation );
		mfem::InverseElementTransformation inverse( &transformation );
		inverse.SetSolverType( mfem::InverseElementTransformation::Newton );
		mfem::IntegrationPoint at;
		mfem::Vector physical( 2 );
		physical( 0 ) = axis.radius;
		physical( 1 ) = axis.z;
		inverse.Transform( physical, at );
		transformation.SetIntPoint( &at );
		mfem::DenseMatrix gradient;
		base.flux->GetVectorGradient( transformation, gradient );
		double const det = gradient( 0, 0 )*gradient( 1, 1 ) - gradient( 0, 1 )*gradient( 1, 0 );
		if ( !( std::abs( det ) > 0.0 ) )
			return false;

		out.qOverG.assign( m*m, 0.0 );
		for ( std::size_t j = 0; j < m; ++j )
		{
			Fields const direction = along( dx[ j ] );
			double const size = direction.psiStar->Normlinf();
			if ( !( size > 0.0 ) )
				continue;
			mfem::Vector dq( 2 );
			direction.flux->GetVectorValue( transformation, at, dq );
			double const dR = -( gradient( 1, 1 )*dq( 0 ) - gradient( 0, 1 )*dq( 1 ) )/det;
			double const dZ = -( -gradient( 1, 0 )*dq( 0 ) + gradient( 0, 0 )*dq( 1 ) )/det;
			double const s = options.step*scale/size;
			auto moved = [ & ]( double t )
			{
				CriticalPoint point = axis;
				point.radius += t*dR;
				point.z += t*dZ;
				point.psi = psiAxis + t*dAxis[ j ];
				return point;
			};
			std::vector<double> const plus = qOverG( combine( base, direction, s ), moved( s ) );
			std::vector<double> const minus = qOverG( combine( base, direction, -s ), moved( -s ) );
			for ( std::size_t i = 0; i < m; ++i )
				out.qOverG[ i*m + j ] = ( plus[ i ] - minus[ i ] )/( 2.0*s );
		}

		/*
		 * THE PIN ROWS, LINEARISED AS QUADRATURES of the source rule through
		 * the solver's own assemblers (their device routes where the state is
		 * resident), contracted in double. assembleSourceDerivative() with
		 * `f` is `-int f phi_i/R` over the potential dofs, so its potential
		 * block summed is `-int f/R` and dotted with a tangent's potential
		 * block `-int f dpsi/R`.
		 */
		out.measured.assign( m, 0.0 );
		int const potentialStart = solver.fluxSpace().GetVSize();
		int const potentialSize = solver.potentialSpace().GetVSize();
		if ( options.pin == ScalePin::ToroidalFlux )
		{
			// Phi = sign int g( Psi )/R, Psi = ( psi - psi_b )/span:
			// dPhi = sign int [ phi_j/( 2g ) + g^2'( Psi ) ( dpsi - Psi dpsi_ax )/( 2g span ) ]/R.
			double const boundary = solver.psiBoundary();
			double const span = psiAxis - boundary;
			auto inside = [ &spline, boundary, span ]( double psi, double &x, double &g )
			{
				x = ( psi - boundary )/span;
				if ( !( x > 0.0 ) )
					return false;
				double const value = spline( x );
				if ( !( value > 0.0 ) )
					return false;
				g = std::sqrt( value );
				return true;
			};
			mfem::Vector slope, slopeTimesPsi;
			solver.assembleSourceDerivative( state,
				[ &spline, &inside, span ]( double, double, double psi )
				{
					double x = 0.0, g = 0.0;
					return inside( psi, x, g ) ? spline( x, 1 )/( 2.0*g*span ) : 0.0;
				}, slope );
			solver.assembleSourceDerivative( state,
				[ &spline, &inside, span ]( double, double, double psi )
				{
					double x = 0.0, g = 0.0;
					return inside( psi, x, g ) ? x*spline( x, 1 )/( 2.0*g*span ) : 0.0;
				}, slopeTimesPsi );
			double const againstAxis = -sumInDouble( slopeTimesPsi, potentialStart, nullptr, 0, potentialSize );
			for ( std::size_t j = 0; j < m; ++j )
			{
				mfem::Vector direct;
				solver.assembleSourceDerivative( state,
					[ &spline, &inside, j ]( double, double, double psi )
					{
						double x = 0.0, g = 0.0;
						return inside( psi, x, g ) ? spline.cardinal( j, x )/( 2.0*g ) : 0.0;
					}, direct );
				double const valueTerm = -sumInDouble( direct, potentialStart, nullptr, 0, potentialSize );
				double const fieldTerm = -sumInDouble( slope, potentialStart, &dx[ j ], potentialStart,
				                                       potentialSize );
				out.measured[ j ] = options.sign*( valueTerm + fieldTerm - dAxis[ j ]*againstAxis );
			}
		}
		else if ( options.pin == ScalePin::PlasmaCurrent )
		{
			mfem::Vector row;
			double againstAxis = 0.0;
			solver.plasmaCurrentLinearisation( state, row, againstAxis );
			for ( std::size_t j = 0; j < m; ++j )
				out.measured[ j ] = sumInDouble( row, potentialStart, &dx[ j ], potentialStart, potentialSize )
				                    + againstAxis*dAxis[ j ]
				                    - sumInDouble( sources[ j ], potentialStart, nullptr, 0, potentialSize );
		}
		return true;
	}

}
