#include "meq/SafetyFactorTangent.hpp"

#include <chrono>
#include <functional>
#include <cmath>
#include <stdexcept>

#include "meq/ConductorField.hpp"
#include "meq/ContinuousField.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/ResidentIterate.hpp"
#include "meq/StageChecksum.hpp"

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

		/// A spline function of the source's Psi as a Profile, for
		/// hermiteKnots(): @a f( psi, derivative ).
		class SplineFunction : public Profile
		{
			public:
				explicit SplineFunction( std::function<double( double, int )> fIn ) : f( std::move( fIn ) ) {}
				double operator()( double psi ) const override { return f( psi, 0 ); }
				double prime( double psi ) const override { return f( psi, 1 ); }
				double doublePrime( double psi ) const override { return f( psi, 2 ); }

			private:
				std::function<double( double, int )> f;
		};

		/*
		 * A PIECEWISE CUBIC OF THE SPLINE'S AS A SplineProfile, EXACTLY: its
		 * Hermite knots at 0, the spline's own knots in Psi and 1, between which
		 * it is a cubic (hermiteKnots()). A SplineProfile has a device form, so
		 * the integrand it enters is a PsiIntegrand the plasma loops' kernel
		 * evaluates (GradShafranovSolver::assemblePsiIntegrandColumn(), M-386).
		 * Past [ 0, 1 ] it continues as a table does; the pins read Psi > 0, and
		 * Psi > 1 only by a located axis's overshoot.
		 */
		std::shared_ptr<Profile const> hermiteOf( ToroidalFieldSpline const &spline,
		                                          std::function<double( double, int )> f )
		{
			std::vector<double> at{ 0.0 };
			std::vector<double> const &labels = spline.labels();
			for ( auto label = labels.rbegin(); label != labels.rend(); ++label )
				at.push_back( 1.0 - *label );
			at.push_back( 1.0 );
			return std::make_shared<SplineProfile const>( hermiteKnots( SplineFunction( std::move( f ) ), at ) );
		}

		std::shared_ptr<Profile const> squareOf( ToroidalFieldSpline const &spline )
		{
			return hermiteOf( spline, [ &spline ]( double psi, int d ) { return spline( psi, d ); } );
		}

		/// `-int f phi_i/R` of @a f, summed over the potential block: `-int f/R`.
		double columnSum( GradShafranovSolver &solver, mfem::Vector const &state, PsiIntegrand const &f,
		                  mfem::Vector &column )
		{
			solver.assemblePsiIntegrandColumn( state, f, column );
			return sumInDouble( column, solver.fluxSpace().GetVSize(), nullptr, 0,
			                    solver.potentialSpace().GetVSize() );
		}
	}

	double plasmaCells( mfem::Mesh const &mesh )
	{
		double cells = 0.0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
			cells += mesh.GetElementGeometry( e ) == mfem::Geometry::TRIANGLE ? 0.5 : 1.0;
		return cells;
	}

	double plasmaCells( GradShafranovSolver const &solver )
	{
		mfem::GridFunction const &psi = solver.potential();
		mfem::Mesh const &mesh = *psi.FESpace()->GetMesh();
		auto weight = [ &mesh ]( int e )
		{
			return mesh.GetElementGeometry( e ) == mfem::Geometry::TRIANGLE ? 0.5 : 1.0;
		};
		double cells = 0.0;
		if ( solver.plasmaComponentElements() > 0 )
		{
			for ( int e = 0; e < mesh.GetNE(); ++e )
				if ( solver.elementInPlasma( e ) )
					cells += weight( e );
			return cells;
		}
		/*
		 * A FREE EDGE WITH NO SUPPORT: the plasma is Psi > 0 of the solved
		 * state, the total field's, read at each element's centre -- a
		 * limited or diverted plasma inside a box it does not fill. Counted on
		 * the mesh, its knots outrun what the plasma's cells resolve, and the
		 * rows go nearly dependent (M-405).
		 */
		double const boundary = solver.psiBoundary();
		double const span = solver.psiAxis() - boundary;
		if ( !solver.boundaryFluxIsAnUnknown() || span == 0.0 )
			return plasmaCells( mesh );
		ConductorField const *conductors = solver.conductorField();
		mfem::IsoparametricTransformation tr;
		mfem::Vector centre( 2 );
		psi.HostRead();
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::IntegrationPoint const &ip = mfem::Geometries.GetCenter( mesh.GetElementBaseGeometry( e ) );
			const_cast<mfem::Mesh &>( mesh ).GetElementTransformation( e, &tr );
			tr.Transform( ip, centre );
			double value = psi.GetValue( e, ip );
			if ( conductors )
				value += conductors->psi( centre( 0 ), centre( 1 ) );
			if ( ( value - boundary )/span > 0.0 )
				cells += weight( e );
		}
		return cells > 0.0 ? cells : plasmaCells( mesh );
	}

	double deliveredToroidalFlux( GradShafranovSolver &solver, ToroidalFieldSpline const &field,
	                              double sign )
	{
		// g( Psi ) = sqrt( G ), G the spline's Hermite form: a PsiIntegrand, on
		// the device under one.
		std::shared_ptr<Profile const> const square = squareOf( field );
		PsiIntegrand f;
		f.kind = PsiIntegrandKind::Root;
		f.square = square.get();
		mfem::Vector const state( solver.state() );
		mfem::Vector assembled;
		return -sign*columnSum( solver, state, f, assembled );
	}

	SafetyFactorEquilibrium::SafetyFactorEquilibrium( GradShafranovSolver &solverIn,
	                                                  NormalisedMHDSource &sourceIn, Options optionsIn )
		: solver( solverIn ), source( sourceIn ), options( std::move( optionsIn ) )
	{
		if ( !( options.step > 0.0 ) )
			throw std::invalid_argument( "meq::SafetyFactorEquilibrium: the difference step must be positive" );
		/*
		 * THE FLUX PIN BESIDE A FREE EDGE IS REFUSED. g( 0 ) = R0B0 is not
		 * zero, so Phi = int_{Psi > 0} g/R jumps by g/R times a weight each
		 * time a quadrature point crosses the edge psi_bnd moves, and is not
		 * differentiable in the state: re-solved differences of it read the
		 * crossings (QDrivenFreeEdge). A fixed boundary's edge is Gamma_h and
		 * does not move. TODO, *The flux pin beside a free plasma edge*.
		 */
		if ( options.pin == ScalePin::ToroidalFlux && solver.boundaryFluxIsAnUnknown() )
			throw std::invalid_argument( "meq::SafetyFactorEquilibrium: the ToroidalFlux pin is refused "
			                             "beside a free plasma edge (a limiter or an X-point): the delivered "
			                             "flux integrates g( 0 ) != 0 up to an edge that moves, so it jumps "
			                             "as quadrature points cross it. Pin R0B0 or PlasmaCurrent" );
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

	CriticalPoint SafetyFactorEquilibrium::locatedAxis( Fields const & ) const
	{
		// The total field's, of psi_ax's sense, seeded at the solve's own axis.
		return solvedAxis( solver );
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
			solver.continuousInterpolant().apply( *fields.psiStar, *fields.qStar, *fields.trace );
		secondsExtracting += secondsSince( started );
		return qOverG( *continuous.potential, *continuous.flux, axis, solver.psiBoundary() );
	}

	std::vector<double> SafetyFactorEquilibrium::qOverG( mfem::GridFunction const &potential,
	                                                     mfem::GridFunction const &flux,
	                                                     CriticalPoint const &axis, double boundary ) const
	{
		auto const started = std::chrono::steady_clock::now();
		FluxSurfaceFamily const family = extractFamily( potential, flux, axis, boundary );
		std::vector<double> out;
		out.reserve( family.surfaces.size() );
		for ( FluxSurface const &surface : family.surfaces )
			out.push_back( surface.vPrime*surface.inverseRSquared/( 4.0*M_PI*M_PI ) );
		secondsExtracting += secondsSince( started );
		extractionCount += 1;
		return out;
	}

	FluxSurfaceFamily SafetyFactorEquilibrium::extractFamily( mfem::GridFunction const &potential,
	                                                          mfem::GridFunction const &flux,
	                                                          CriticalPoint const &axis, double boundary ) const
	{
		ContourTracer tracer( potential, flux );
		// The conductors' field, constant along every tangent, on top of the
		// plasma's -- the surfaces are the total field's.
		tracer.setConductorField( solver.conductorField() );
		if ( options.gammaHMarker != nullptr )
			tracer.setBandExtension( BandExtension::TransferLift, *options.gammaHMarker, options.path );
		return extractFluxSurfaces( tracer, axis, boundary, options.family );
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
		return deliveredToroidalFlux( solver, spline, options.sign );
	}

	void SafetyFactorEquilibrium::warmFromSolver()
	{
		warm = solver.state();
		warmExterior = solver.exteriorCoefficients();
		warmBoundary = solver.psiBoundary();
		warmAxis = solver.psiAxis();
	}

	bool SafetyFactorEquilibrium::evaluate( ToroidalFieldSpline const &spline, SafetyFactorEvaluation &out )
	{
		auto const started = std::chrono::steady_clock::now();
		try
		{
			source.setGGPrime( spline.ggPrime() );
			if ( warm.Size() > 0 )
			{
				// The whole iterate: the field, the exterior modes, and psi_bnd
				// where it is free -- else a refused trial's would be the start.
				solver.setInitialState( warm, warmExterior );
				solver.setBoundaryFluxInitial( warmBoundary );
				// And psi_ax's label: the source keeps the last solve's, a refused
				// trial's included, and a start whose label is that trial's is
				// not the warm state's (M-406).
				if ( warmAxis != 0.0 )
					solver.setAxisFluxInitial( warmAxis );
			}
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
		catch ( std::exception const &error )
		{
			lastRefusal = error.what();
			secondsSolving += secondsSince( started );
			admittedValues.clear();
			return false;
		}
		warm = solver.state();
		warmExterior = solver.exteriorCoefficients();
		warmBoundary = solver.psiBoundary();
		warmAxis = solver.psiAxis();
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

		/*
		 * EVERY BORDER THE SOLVE CARRIES IS IN THE TANGENT: psi_ax always,
		 * psi_bnd where the plasma edge is free (a limiter, the X-point),
		 * the X-point's position and the exterior modes. The extraction reads
		 * the first two as its labels, so both move with the tangent.
		 */
		GradShafranovSolver::LinearisedBorders const layout = solver.linearisedLayout();
		std::vector<mfem::Vector> dx( m ), sources( m );
		std::vector<double> dAxis( m, 0.0 ), dBoundary( m, 0.0 );
		std::vector<std::shared_ptr<Profile const>> directions( m );
		for ( std::size_t j = 0; j < m; ++j )
		{
			directions[ j ] = spline.cardinalGGPrime( j );
			Profile const &direction = *directions[ j ];
			mfem::Vector &dRdu = sources[ j ];
			// The source's own kernel on a substituted view under a Device;
			// ggPrimeDerivative() through the host loop otherwise (M-409).
			solver.assembleGGPrimeDerivative( state, direction, dRdu );
			std::vector<double> ds;
			auto checksum = []( char const *stage, mfem::Vector const &v )
			{
				if ( !stageChecksumWanted() )
					return;
				mfem::real_t const *host = v.HostRead();
				stageChecksum( stage, std::vector<double>( host, host + v.Size() ) );
			};
			checksum( "toroidalTangent.source", dRdu );
			solver.linearisedTangent( dRdu, {}, nullptr, dx[ j ], ds );
			checksum( "toroidalTangent.state", dx[ j ] );
			if ( layout.axis >= 0 )
				dAxis[ j ] = ds[ static_cast<std::size_t>( layout.axis ) ];
			if ( layout.boundary >= 0 )
				dBoundary[ j ] = ds[ static_cast<std::size_t>( layout.boundary ) ];
		}
		secondsTangents += secondsSince( started );
		lastAxisDerivative = dAxis;
		stageChecksum( "toroidalTangent.borders", dAxis );
		lastBoundaryDerivative = dBoundary;

		Fields const base = solvedFields();
		double const scale = base.psiStar->Normlinf();
		double const psiAxis = solver.psiAxis();
		double const psiBoundary = solver.psiBoundary();
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
		if ( ConductorField const *conductors = solver.conductorField() )
		{
			// The root is the total q's: the conductors' gradient, centrally
			// differenced as the bordered Newton's X-point corner takes it.
			double const step = 1.0e-6*std::max( 1.0, std::abs( axis.radius ) );
			double aR = 0.0, aZ = 0.0, bR = 0.0, bZ = 0.0;
			conductors->flux( axis.radius + step, axis.z, aR, aZ );
			conductors->flux( axis.radius - step, axis.z, bR, bZ );
			gradient( 0, 0 ) += ( aR - bR )/( 2.0*step );
			gradient( 1, 0 ) += ( aZ - bZ )/( 2.0*step );
			conductors->flux( axis.radius, axis.z + step, aR, aZ );
			conductors->flux( axis.radius, axis.z - step, bR, bZ );
			gradient( 0, 1 ) += ( aR - bR )/( 2.0*step );
			gradient( 1, 1 ) += ( aZ - bZ )/( 2.0*step );
		}
		double const det = gradient( 0, 0 )*gradient( 1, 1 ) - gradient( 0, 1 )*gradient( 1, 0 );
		if ( !( std::abs( det ) > 0.0 ) )
			return false;

		/*
		 * psi~ AND q~ ARE LINEAR IN THE STATE (ContinuousField.hpp): the
		 * solver's kept interpolant applied once for the base and once per
		 * tangent (meq::ContinuousInterpolant, its quadratures done once per
		 * mesh, on the device under one).
		 */
		auto const building = std::chrono::steady_clock::now();
		ContinuousFields const baseContinuous =
			solver.continuousInterpolant().apply( *base.psiStar, *base.qStar, *base.trace );
		secondsExtracting += secondsSince( building );
		CriticalPoint labelled = axis;
		labelled.psi = psiAxis;

		// The axis root's motion along tangent j, in its own element.
		auto originMotion = [ & ]( Fields const &direction, double &dR, double &dZ )
		{
			mfem::Vector dq( 2 );
			direction.flux->GetVectorValue( transformation, at, dq );
			dR = -( gradient( 1, 1 )*dq( 0 ) - gradient( 0, 1 )*dq( 1 ) )/det;
			dZ = -( -gradient( 1, 0 )*dq( 0 ) + gradient( 0, 0 )*dq( 1 ) )/det;
		};

		out.qOverG.assign( m*m, 0.0 );
		if ( !options.differenced )
		{
			std::vector<double> originR( m, 0.0 ), originZ( m, 0.0 );
			for ( std::size_t j = 0; j < m; ++j )
				originMotion( along( dx[ j ] ), originR[ j ], originZ[ j ] );
			auto const shaping = std::chrono::steady_clock::now();
			bool const shaped = shapeDerivative( baseContinuous, labelled, psiBoundary, dx, dAxis, dBoundary,
			                                     originR, originZ, out.qOverG );
			secondsExtracting += secondsSince( shaping );
			if ( !shaped )
				return false;
		}
		else
		{
			/*
			 * THE DIFFERENCED EXTRACTION (Options::differenced): the perturbed
			 * pair is the base's plus s times the tangent's, extracted centrally
			 * by default (Options::step, M-399).
			 */
			std::vector<double> const atBase = options.central ? std::vector<double>()
				: qOverG( *baseContinuous.potential, *baseContinuous.flux, labelled, psiBoundary );
			mfem::GridFunction potential( baseContinuous.potential->FESpace() );
			mfem::GridFunction flux( baseContinuous.flux->FESpace() );

			for ( std::size_t j = 0; j < m; ++j )
			{
				Fields const direction = along( dx[ j ] );
				double const size = direction.psiStar->Normlinf();
				if ( !( size > 0.0 ) )
					continue;
				double dR = 0.0, dZ = 0.0;
				originMotion( direction, dR, dZ );
				double const s = options.step*scale/size;
				auto moved = [ & ]( double t )
				{
					CriticalPoint point = axis;
					point.radius += t*dR;
					point.z += t*dZ;
					point.psi = psiAxis + t*dAxis[ j ];
					return point;
				};
				auto const interpolating = std::chrono::steady_clock::now();
				ContinuousFields const tangent =
					solver.continuousInterpolant().apply( *direction.psiStar, *direction.qStar, *direction.trace );
				secondsExtracting += secondsSince( interpolating );
				auto perturbed = [ & ]( double t )
				{
					potential = *baseContinuous.potential;
					potential.Add( t, *tangent.potential );
					flux = *baseContinuous.flux;
					flux.Add( t, *tangent.flux );
					return qOverG( potential, flux, moved( t ), psiBoundary + t*dBoundary[ j ] );
				};
				std::vector<double> const plus = perturbed( s );
				if ( options.central )
				{
					std::vector<double> const minus = perturbed( -s );
					for ( std::size_t i = 0; i < m; ++i )
						out.qOverG[ i*m + j ] = ( plus[ i ] - minus[ i ] )/( 2.0*s );
				}
				else
					for ( std::size_t i = 0; i < m; ++i )
						out.qOverG[ i*m + j ] = ( plus[ i ] - atBase[ i ] )/s;
			}
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
			// dPhi = sign int [ phi_j/( 2g ) + g^2'( Psi ) ( dpsi - dpsi_b
			//                   - Psi ( dpsi_ax - dpsi_b ) )/( 2g span ) ]/R.
			// The interior's change: where psi_b is free the edge Psi = 0
			// moves too, and g( 0 ) = R0B0 is not zero there, so the
			// quadrature's own edge is what this does not carry (TODO, *The
			// flux pin beside a free plasma edge*).
			/*
			 * THE FOUR INTEGRANDS AS DATA (PsiIntegrandKind): G'/( 2 g span ),
			 * Psi G'/( 2 g span ) and each cardinal H_j/( 2 g ), G and H_j the
			 * spline's Hermite forms, so each column is the plasma loops' kernel
			 * under a Device (assemblePsiIntegrandColumn(), M-386).
			 */
			std::shared_ptr<Profile const> const square = squareOf( spline );
			PsiIntegrand slopeIntegrand;
			slopeIntegrand.kind = PsiIntegrandKind::SlopeOverRoot;
			slopeIntegrand.square = square.get();
			PsiIntegrand psiSlopeIntegrand = slopeIntegrand;
			psiSlopeIntegrand.kind = PsiIntegrandKind::PsiSlopeOverRoot;
			mfem::Vector slope, slopeTimesPsi;
			double const slopeSum = columnSum( solver, state, slopeIntegrand, slope );
			double const againstAxis = -columnSum( solver, state, psiSlopeIntegrand, slopeTimesPsi );
			double const againstBoundary = layout.boundary >= 0 ? -slopeSum - againstAxis : 0.0;
			for ( std::size_t j = 0; j < m; ++j )
			{
				std::shared_ptr<Profile const> const cardinal = hermiteOf( spline,
					[ &spline, j ]( double psi, int d ) { return spline.cardinal( j, psi, d ); } );
				PsiIntegrand ratio;
				ratio.kind = PsiIntegrandKind::RatioOverRoot;
				ratio.square = square.get();
				ratio.numerator = cardinal.get();
				mfem::Vector direct;
				double const valueTerm = -columnSum( solver, state, ratio, direct );
				double const fieldTerm = -sumInDouble( slope, potentialStart, &dx[ j ], potentialStart,
				                                       potentialSize );
				out.measured[ j ] = options.sign*( valueTerm + fieldTerm - dAxis[ j ]*againstAxis
				                                   - dBoundary[ j ]*againstBoundary );
			}
		}
		else if ( options.pin == ScalePin::PlasmaCurrent )
		{
			mfem::Vector row;
			double againstAxis = 0.0, againstBoundary = 0.0;
			solver.plasmaCurrentLinearisation( state, row, againstAxis, &againstBoundary );
			for ( std::size_t j = 0; j < m; ++j )
				out.measured[ j ] = sumInDouble( row, potentialStart, &dx[ j ], potentialStart, potentialSize )
				                    + againstAxis*dAxis[ j ] + againstBoundary*dBoundary[ j ]
				                    - sumInDouble( sources[ j ], potentialStart, nullptr, 0, potentialSize );
		}
		return true;
	}

}
