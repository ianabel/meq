#include "Estimator.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace meq
{

	namespace
	{
		/// The diameter of a set of mesh vertices: the largest distance between
		/// any two of them. See elementDiameter() in the header for why this and
		/// not mfem::Mesh::GetElementSize().
		double diameter( mfem::Mesh &mesh, mfem::Array<int> const &vertices )
		{
			int const dim = mesh.SpaceDimension();
			double largest = 0.0;
			for ( int i = 0; i < vertices.Size(); ++i )
			{
				for ( int j = i + 1; j < vertices.Size(); ++j )
				{
					largest = std::max( largest,
					                    mfem::Distance( mesh.GetVertex( vertices[ i ] ),
					                                    mesh.GetVertex( vertices[ j ] ),
					                                    dim ) );
				}
			}
			return largest;
		}
	}

	double elementDiameter( mfem::Mesh &mesh, int element )
	{
		mfem::Array<int> vertices;
		mesh.GetElementVertices( element, vertices );
		return diameter( mesh, vertices );
	}

	double faceDiameter( mfem::Mesh &mesh, int face )
	{
		mfem::Array<int> vertices;
		mesh.GetFaceVertices( face, vertices );
		return diameter( mesh, vertices );
	}

	ResidualEstimator::ResidualEstimator( GradShafranovSolver &solverIn,
	                                      Source const &sourceIn )
		: solver( &solverIn ),
		  source( &sourceIn ),
		  sourceCoeff( nullptr ),
		  potentialChoice( Potential::PostProcessed ),
		  comparisonChoice( TraceComparison::Projected ),
		  extraQuadratureOrder( 4 ),
		  sequence( -1 )
	{
		sums.fill( 0.0 );
	}

	ResidualEstimator::ResidualEstimator( GradShafranovSolver &solverIn,
	                                      mfem::Coefficient &sourceIn )
		: solver( &solverIn ),
		  source( nullptr ),
		  sourceCoeff( &sourceIn ),
		  potentialChoice( Potential::PostProcessed ),
		  comparisonChoice( TraceComparison::Projected ),
		  extraQuadratureOrder( 4 ),
		  sequence( -1 )
	{
		sums.fill( 0.0 );
	}

	void ResidualEstimator::setPotential( Potential potentialIn )
	{
		if ( potentialIn != potentialChoice )
			Reset();
		potentialChoice = potentialIn;
	}

	ResidualEstimator::Potential ResidualEstimator::potential() const
	{
		return potentialChoice;
	}

	void ResidualEstimator::setTraceComparison( TraceComparison comparisonIn )
	{
		if ( comparisonIn != comparisonChoice )
			Reset();
		comparisonChoice = comparisonIn;
	}

	ResidualEstimator::TraceComparison ResidualEstimator::traceComparison() const
	{
		return comparisonChoice;
	}

	void ResidualEstimator::setTransferredBoundary( mfem::Array<int> const &markerIn,
	                                                mfem::Coefficient *datumIn )
	{
		markerIn.Copy( transferredBoundary );
		transferredDatumCoefficient = datumIn;
		Reset();
	}

	void ResidualEstimator::setExtraQuadratureOrder( int extraIn )
	{
		if ( extraIn != extraQuadratureOrder )
			Reset();
		extraQuadratureOrder = extraIn;
	}

	void ResidualEstimator::Reset()
	{
		sequence = -1;
	}

	char const *ResidualEstimator::name( Term term )
	{
		switch ( term )
		{
			case Term::Divergence:    return "eta_1";
			case Term::Constitutive:  return "eta_2";
			case Term::FluxJump:      return "eta_3";
			case Term::PotentialJump: return "eta_4";
			case Term::TraceMismatch: return "eta_5";
		}
		return "eta_?";
	}

	double ResidualEstimator::sourceValue( mfem::ElementTransformation &tr,
	                                       mfem::IntegrationPoint const &ip,
	                                       double r, double z, double psi ) const
	{
		if ( source )
			return source->f( r, z, psi );
		return sourceCoeff->Eval( tr, ip );
	}

	mfem::Vector const &ResidualEstimator::GetLocalErrors()
	{
		compute();
		return errors;
	}

	mfem::real_t ResidualEstimator::GetTotalError() const
	{
		compute();
		double sum = 0.0;
		for ( double value : sums )
			sum += value;
		return std::sqrt( sum );
	}

	double ResidualEstimator::component( Term term ) const
	{
		compute();
		return std::sqrt( sums[ static_cast<int>( term ) ] );
	}

	mfem::Vector const &ResidualEstimator::localSquares( Term term ) const
	{
		compute();
		return squares[ static_cast<int>( term ) ];
	}

	void ResidualEstimator::compute() const
	{
		// psi* is the whole reason this class needs the solver rather than a bag
		// of GridFunctions, and a zero psi* would make every term look like a
		// slightly wrong version of the psi_h estimator instead of failing. So it
		// is refused rather than silently accepted.
		if ( potentialChoice == Potential::PostProcessed && !solver->isPostProcessed() )
			throw std::logic_error( "meq::ResidualEstimator: the estimator uses the post-processed potential, but GradShafranovSolver::postProcess() has not been called" );

		mfem::GridFunction const &fluxGf = solver->flux();
		mfem::GridFunction const &traceGf = solver->trace();
		mfem::GridFunction const &psiGf =
			( potentialChoice == Potential::PostProcessed )
				? solver->postProcessedPotential()
				: solver->potential();

		mfem::FiniteElementSpace &traceFes = solver->traceSpace();
		mfem::Mesh &mesh = *solver->potentialSpace().GetMesh();

		if ( mesh.GetSequence() == sequence && errors.Size() == mesh.GetNE() )
			return;

		int const numElements = mesh.GetNE();
		int const dim = mesh.Dimension();

		for ( mfem::Vector &term : squares )
		{
			term.SetSize( numElements );
			term = 0.0;
		}
		errors.SetSize( numElements );
		errors = 0.0;

		// Twice the degree of the potential in use, plus room for the fact that
		// neither 1/r nor F is a polynomial. Whichever potential is selected, the
		// two branches are integrated on the SAME rule, so the psi*-against-psi_h
		// comparison in EstimatorConvergence.cpp is not a comparison of two
		// quadratures.
		int const psiOrder = psiGf.FESpace()->GetMaxElementOrder();
		int const quadratureOrder = 2*psiOrder + extraQuadratureOrder;

		mfem::Vector point;
		mfem::Vector fluxValue;
		mfem::Vector gradient;

		// ---------------------------------------------------------------------
		// The element terms, eta_1 and eta_2.
		// ---------------------------------------------------------------------
		for ( int e = 0; e < numElements; ++e )
		{
			mfem::ElementTransformation *tr = mesh.GetElementTransformation( e );
			double const hK = elementDiameter( mesh, e );
			mfem::IntegrationRule const &ir =
				mfem::IntRules.Get( mesh.GetElementGeometry( e ), quadratureOrder );

			double divergenceTerm = 0.0;
			double constitutiveTerm = 0.0;

			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				tr->SetIntPoint( &ip );
				tr->Transform( ip, point );

				double const r = point( 0 );
				double const z = point( 1 );
				double const weight = ip.weight*tr->Weight();
				double const psi = psiGf.GetValue( e, ip );

				// eta_1: the residual of -div_bar q = F/r, so F/r + div_bar q.
				// F is evaluated at the potential in use -- psi* for the published
				// estimator -- which is what makes this the residual of the
				// semi-linear equation.
				double const residual = sourceValue( *tr, ip, r, z, psi )/r
				                        + fluxGf.GetDivergence( *tr );
				divergenceTerm += weight*residual*residual;

				// eta_2: the residual of q = ( 1/r ) grad_bar psi. This is the
				// term that differentiates the potential, and the term the paper
				// says loses an order if the potential is psi_h.
				fluxGf.GetVectorValue( e, ip, fluxValue );
				psiGf.GetGradient( *tr, gradient );
				gradient /= r;
				fluxValue -= gradient;
				constitutiveTerm += weight*( fluxValue*fluxValue );
			}

			squares[ static_cast<int>( Term::Divergence ) ]( e ) = hK*hK*divergenceTerm;
			squares[ static_cast<int>( Term::Constitutive ) ]( e ) = constitutiveTerm;
		}

		// ---------------------------------------------------------------------
		// The edge terms, eta_3 to eta_5.
		//
		// The paper writes these as a sum over elements of a sum over the edges
		// of each, which counts every interior edge twice; the 1/2 in front of
		// eta_3 and eta_4 is what undoes that. Iterating over FACES instead and
		// giving each of the two neighbours half is the same arithmetic done once
		// -- and it has to be a face loop anyway, because a jump needs both sides.
		//
		// eta_5 carries no 1/2 and is genuinely two different numbers on an
		// interior edge: psi* has a different trace from each side, and it is
		// exactly that difference from psihat_h that the term measures.
		// ---------------------------------------------------------------------
		mfem::Array<int> traceVdofs;
		mfem::Vector traceLocal;
		mfem::Vector traceShape;
		mfem::Vector normal( dim );
		mfem::Vector fluxOne;
		mfem::Vector fluxTwo;
		mfem::DenseMatrix faceMass;
		mfem::Vector faceLoad;
		mfem::Vector datumLoad;
		mfem::Vector datumCoefficients;
		mfem::Array<int> const faceToBdr = mesh.GetFaceToBdrElMap();
		mfem::Vector projected;
		mfem::Vector residualCoefficients;
		mfem::Vector massTimesResidual;

		for ( int f = 0; f < mesh.GetNumFaces(); ++f )
		{
			mfem::FaceElementTransformations *ftr = mesh.GetFaceElementTransformations( f );
			if ( !ftr )
				continue;

			int const elementOne = ftr->Elem1No;
			int const elementTwo = ftr->Elem2No;
			bool const interior = ( elementTwo >= 0 );
			double const he = faceDiameter( mesh, f );

			// Whether eta_5 leaves this face out. Only a boundary face can be
			// excluded, and only when setTransferredBoundary() has named its
			// attribute -- see that method for what psihat_h is worth on such a
			// face and what leaving it out costs.
			bool transferred = false;
			if ( !interior && transferredBoundary.Size() > 0 )
			{
				int const bdrElement = faceToBdr[ f ];
				if ( bdrElement >= 0 )
				{
					int const attribute = mesh.GetBdrAttribute( bdrElement );
					transferred = attribute >= 1 && attribute <= transferredBoundary.Size()
					              && transferredBoundary[ attribute - 1 ];
				}
			}

			// psihat_h on this edge. The trace collection is
			// DG_Interface_FECollection with the default VALUE map type, so the
			// shape functions are values and shape . dofs is psihat_h -- there is
			// no quadrature weight to divide out, as there would be for an
			// INTEGRAL face element. And every face carries its own dofs, so no
			// orientation question arises: GetFaceVDofs and GetFaceElement agree
			// on the face's own reference frame, which is the frame the face
			// transformation below uses too.
			//
			// This is the route the deleted GridFunction::GetValueFacet used to
			// provide. It does not exist in MFEM 4.9.1, and the same pattern is
			// what fem/darcy/estimators_hdg.cpp and darcyform.cpp's reconstruction
			// use to reach a trace value.
			mfem::FiniteElement const *traceFe = traceFes.GetFaceElement( f );
			int const traceDof = traceFe->GetDof();
			traceFes.GetFaceVDofs( f, traceVdofs );
			traceGf.GetSubVector( traceVdofs, traceLocal );
			traceShape.SetSize( traceDof );

			mfem::IntegrationRule const &ir =
				mfem::IntRules.Get( ftr->GetGeometryType(), quadratureOrder );

			double fluxJump = 0.0;
			double potentialJump = 0.0;

			// eta_5, one number per side of the face -- see the block comment
			// above -- and each of them needed in both forms: as printed, and with
			// the difference taken inside M_h. The face mass matrix and the load
			// vector of psi* against the trace basis are what the second needs, and
			// they are accumulated on the same quadrature pass as the first.
			double mismatch[ 2 ] = { 0.0, 0.0 };
			double projectedMismatch[ 2 ] = { 0.0, 0.0 };
			int const sides = interior ? 2 : 1;
			int const element[ 2 ] = { elementOne, elementTwo };

			// On a transferred face with a datum to hand, psihat_h is NOT the
			// condition imposed -- it is a zero standing in for one -- so eta_5
			// reads phi_h instead. Without a datum the face is left out, which is
			// what this class did unconditionally before.
			bool const useDatum = transferred && transferredDatumCoefficient;

			for ( int side = 0; side < sides; ++side )
			{
				faceMass.SetSize( traceDof );
				faceMass = 0.0;
				faceLoad.SetSize( traceDof );
				faceLoad = 0.0;
				datumLoad.SetSize( traceDof );
				datumLoad = 0.0;

				for ( int i = 0; i < ir.GetNPoints(); ++i )
				{
					mfem::IntegrationPoint const &ip = ir.IntPoint( i );

					// Sets the face point and, through Loc1/Loc2, the
					// corresponding points in both neighbouring elements.
					ftr->SetAllIntPoints( &ip );

					// The face measure. FaceElementTransformations IS the face
					// transformation, so Weight() here is |dx/dxi| along the edge
					// and ip.weight sums to one over the reference segment: the
					// product is ds. Leaving it out -- which the pre-modernisation
					// estimator did -- turns every edge term into an unweighted
					// average and makes h_e the only geometry in them.
					double const weight = ip.weight*ftr->Weight();

					traceFe->CalcShape( ip, traceShape );

					// Both halves of phi_h want the FACE transformation and the
					// FACE integration point -- mfem::PathLiftCoefficient
					// MFEM_VERIFYs the dynamic_cast, since a lifting is defined
					// along a path issuing from a face and a path family may need
					// the outward normal. ftr->SetAllIntPoints( &ip ) above is
					// what makes both usable here.
					double const hatValue = useDatum
						? transferredDatumCoefficient->Eval( *ftr, ip )
						: traceShape*traceLocal;

					mfem::IntegrationPoint const &eip =
						side ? ftr->GetElement2IntPoint() : ftr->GetElement1IntPoint();
					double const psiValue = psiGf.GetValue( element[ side ], eip );

					mismatch[ side ] += weight*( hatValue - psiValue )
					                          *( hatValue - psiValue );

					mfem::AddMult_a_VVt( weight, traceShape, faceMass );
					faceLoad.Add( weight*psiValue, traceShape );
					if ( useDatum )
						datumLoad.Add( weight*hatValue, traceShape );

					if ( side != 0 )
						continue;

					// eta_3 and eta_4 need both sides at once, so they are done on
					// the first pass, where the second element's integration point
					// is available from the same SetAllIntPoints().
					if ( !interior )
						continue;

					mfem::IntegrationPoint const &ipOne = ftr->GetElement1IntPoint();
					mfem::IntegrationPoint const &ipTwo = ftr->GetElement2IntPoint();
					double const psiOne = psiGf.GetValue( elementOne, ipOne );
					double const psiTwo = psiGf.GetValue( elementTwo, ipTwo );

					// eta_4: [[ psi* ]] = psi*+ - psi*-.
					potentialJump += weight*( psiOne - psiTwo )*( psiOne - psiTwo );

					// eta_3: [[ q ]] = q+ . n+ + q- . n- = ( q+ - q- ) . n+, a
					// scalar. CalcOrtho returns the normal scaled by the face
					// Jacobian, so it must be divided by the same weight to be a
					// unit normal -- omit that and the jump is multiplied by the
					// edge length, which is a second h_e nobody asked for.
					mfem::CalcOrtho( ftr->Jacobian(), normal );
					normal /= ftr->Weight();

					fluxGf.GetVectorValue( elementOne, ipOne, fluxOne );
					fluxGf.GetVectorValue( elementTwo, ipTwo, fluxTwo );
					fluxOne -= fluxTwo;
					double const jump = fluxOne*normal;
					fluxJump += weight*jump*jump;
				}

				// || psihat_h - P_M psi* ||_e^2 = ( b - c )^T M ( b - c ), with
				// M c = load the L2( e ) projection of psi*|_K onto the trace
				// space. b - c is what is left of the difference once the part of
				// psi* that M_h cannot represent is taken out -- and that part is
				// O( h^(k+1) ) even for the EXACT solution, which is why the two
				// forms converge at different rates. See the file comment.
				mfem::DenseMatrixInverse inverse( faceMass );
				projected.SetSize( traceDof );
				inverse.Mult( faceLoad, projected );

				residualCoefficients.SetSize( traceDof );
				if ( useDatum )
				{
					// phi_h is a general function on the face, not an element of
					// M_h, so its coefficients have to be found the same way
					// psi*'s were rather than read off a GridFunction. Where the
					// datum is NOT in use this stays traceLocal, which is EXACT
					// -- psihat_h is already in M_h -- so every number this
					// estimator produced before is reproduced bit for bit.
					datumCoefficients.SetSize( traceDof );
					inverse.Mult( datumLoad, datumCoefficients );
					residualCoefficients = datumCoefficients;
				}
				else
					residualCoefficients = traceLocal;
				residualCoefficients -= projected;
				massTimesResidual.SetSize( traceDof );
				faceMass.Mult( residualCoefficients, massTimesResidual );
				projectedMismatch[ side ] = residualCoefficients*massTimesResidual;
			}

			if ( !transferred || useDatum )
			{
				for ( int side = 0; side < sides; ++side )
				{
					double const chosen = ( comparisonChoice == TraceComparison::Projected )
						? projectedMismatch[ side ] : mismatch[ side ];
					squares[ static_cast<int>( Term::TraceMismatch ) ]( element[ side ] )
						+= chosen/he;
				}
			}

			if ( interior )
			{
				squares[ static_cast<int>( Term::FluxJump ) ]( elementOne )
					+= 0.5*he*fluxJump;
				squares[ static_cast<int>( Term::FluxJump ) ]( elementTwo )
					+= 0.5*he*fluxJump;
				squares[ static_cast<int>( Term::PotentialJump ) ]( elementOne )
					+= 0.5*potentialJump/he;
				squares[ static_cast<int>( Term::PotentialJump ) ]( elementTwo )
					+= 0.5*potentialJump/he;
			}
		}

		for ( int t = 0; t < termCount; ++t )
		{
			sums[ t ] = squares[ t ].Sum();
			for ( int e = 0; e < numElements; ++e )
				errors( e ) += squares[ t ]( e );
		}
		for ( int e = 0; e < numElements; ++e )
			errors( e ) = std::sqrt( errors( e ) );

		sequence = mesh.GetSequence();
	}

	void markDoerfler( mfem::Vector const &localErrors, double gamma,
	                   mfem::Array<int> &marked )
	{
		if ( gamma <= 0.0 || gamma > 1.0 )
			throw std::invalid_argument( "meq::markDoerfler: the marking parameter must lie in ( 0, 1 ]" );

		int const count = localErrors.Size();
		std::vector<int> order( count );
		std::iota( order.begin(), order.end(), 0 );
		std::sort( order.begin(), order.end(), [ &localErrors ]( int a, int b )
		{
			return localErrors( a ) > localErrors( b );
		} );

		double total = 0.0;
		for ( int e = 0; e < count; ++e )
			total += localErrors( e )*localErrors( e );

		// A minimal set, so the accumulation stops at the first index that
		// reaches the bulk rather than after it. With a total of zero there is
		// nothing to refine and nothing is marked -- which is the honest answer,
		// and is also the condition GS-2 section 3.3 names as the one under which
		// the computational domains stop exhausting Omega.
		marked.SetSize( 0 );
		double accumulated = 0.0;
		double const wanted = gamma*total;
		for ( int i = 0; i < count && accumulated < wanted; ++i )
		{
			marked.Append( order[ i ] );
			accumulated += localErrors( order[ i ] )*localErrors( order[ i ] );
		}
	}

	void markMaximum( mfem::Vector const &localErrors, double gamma,
	                  mfem::Array<int> &marked )
	{
		if ( gamma < 0.0 || gamma > 1.0 )
			throw std::invalid_argument( "meq::markMaximum: the marking parameter must lie in [ 0, 1 ]" );

		double largest = 0.0;
		for ( int e = 0; e < localErrors.Size(); ++e )
			largest = std::max( largest, localErrors( e ) );

		marked.SetSize( 0 );
		if ( largest <= 0.0 )
			return;

		double const threshold = gamma*largest;
		for ( int e = 0; e < localErrors.Size(); ++e )
		{
			if ( localErrors( e ) >= threshold )
				marked.Append( e );
		}
	}

	namespace
	{
		/// Whether an element meets Omega: some sampled point has phi <= 0.
		///
		/// Sampled on the same lattice mfem::MarkLevelSetSubdomain uses for the
		/// inside test -- the vertices, plus the refined reference lattice when
		/// extraRefine >= 1 -- so the two markings are decided from the same
		/// evidence and cannot disagree about an element. An element the inside
		/// test rejects because one lattice point had phi > 0 is caught here by
		/// the points that had phi <= 0, which is exactly the band Gamma cuts.
		bool meetsDomain( mfem::Mesh &mesh, int element,
		                  mfem::PositionFunction const &phi, int extraRefine )
		{
			int const dim = mesh.SpaceDimension();
			mfem::Vector x( dim );
			mfem::Array<int> vertices;

			mesh.GetElementVertices( element, vertices );
			for ( int v = 0; v < vertices.Size(); ++v )
			{
				double const *c = mesh.GetVertex( vertices[ v ] );
				for ( int d = 0; d < dim; ++d )
					x( d ) = c[ d ];
				if ( phi( x ) <= 0.0 )
					return true;
			}

			if ( extraRefine < 1 )
				return false;

			mfem::RefinedGeometry *refined = mfem::GlobGeometryRefiner.Refine(
				mesh.GetElementGeometry( element ), extraRefine );
			mfem::IsoparametricTransformation tr;
			mesh.GetElementTransformation( element, &tr );
			for ( int q = 0; q < refined->RefPts.GetNPoints(); ++q )
			{
				tr.Transform( refined->RefPts.IntPoint( q ), x );
				if ( phi( x ) <= 0.0 )
					return true;
			}
			return false;
		}
	}

	AdaptiveDomain::AdaptiveDomain( mfem::Mesh const &backgroundIn,
	                                mfem::PositionFunction levelSetIn,
	                                int extraRefineIn )
		: backgroundMesh( backgroundIn, true ),
		  levelSet( std::move( levelSetIn ) ),
		  extraRefine( extraRefineIn ),
		  gammaH( 0 ),
		  refinementCount( 0 ),
		  proximityAdditions( 0 )
	{
		if ( backgroundMesh.Dimension() != 2 )
			throw std::invalid_argument( "meq::AdaptiveDomain: the background mesh must be two dimensional ( r, z )" );
		select();
	}

	void AdaptiveDomain::select()
	{
		// Step 5, and step 1 on the first call: T_h is the elements lying
		// entirely inside Omega and T_c^h is the minimal cover, the elements that
		// meet Omega at all. Re-selected from scratch after every refinement
		// rather than tracked through it, which is what makes the bookkeeping
		// trivial -- the element numbering changes under GeneralRefinement and
		// nothing here depends on it surviving.
		int const inside = mfem::MarkLevelSetSubdomain( backgroundMesh, levelSet, 0.0,
		                                                insideMarker, extraRefine );
		if ( inside <= 0 )
			throw std::runtime_error( "meq::AdaptiveDomain: no background element lies inside Omega" );

		companionMarker.SetSize( backgroundMesh.GetNE() );
		for ( int e = 0; e < backgroundMesh.GetNE(); ++e )
		{
			companionMarker[ e ] = insideMarker[ e ]
				|| meetsDomain( backgroundMesh, e, levelSet, extraRefine ) ? 1 : 0;
		}

		// Attribute 1 is T_h, 2 the rest of the companion, 3 the background that
		// takes part in neither. Only the first is used to cut the SubMesh; the
		// other two are set so that a mesh written out for a picture shows the
		// three sets apart, which is the one way to check Figure 4 by eye.
		for ( int e = 0; e < backgroundMesh.GetNE(); ++e )
		{
			int const attribute = insideMarker[ e ] ? 1 : ( companionMarker[ e ] ? 2 : 3 );
			backgroundMesh.SetAttribute( e, attribute );
		}
		backgroundMesh.SetAttributes();

		mfem::Array<int> domainAttribute( 1 );
		domainAttribute[ 0 ] = 1;
		computationalMesh = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( backgroundMesh, domainAttribute ) );

		// SubMesh gives the boundary it had to generate one new attribute, one
		// past whatever the parent already used, and leaves inherited boundary
		// with the attributes it had. Omega is strictly inside the box here, so
		// the whole of Gamma_h is generated and there is exactly one attribute.
		gammaH = computationalMesh->bdr_attributes.Max();
		if ( computationalMesh->bdr_attributes.Size() != 1 )
			throw std::runtime_error( "meq::AdaptiveDomain: the computational mesh has boundary inherited from the background box, so Omega is not strictly inside it" );

		gammaHMarkerValue.SetSize( gammaH );
		gammaHMarkerValue = 0;
		gammaHMarkerValue[ gammaH - 1 ] = 1;
	}

	mfem::SubMesh &AdaptiveDomain::computational()
	{
		return *computationalMesh;
	}

	int AdaptiveDomain::gammaHAttribute() const
	{
		return gammaH;
	}

	mfem::Array<int> const &AdaptiveDomain::gammaHMarker() const
	{
		return gammaHMarkerValue;
	}

	int AdaptiveDomain::numComputational() const
	{
		return computationalMesh->GetNE();
	}

	int AdaptiveDomain::numCompanion() const
	{
		int count = 0;
		for ( int e = 0; e < companionMarker.Size(); ++e )
			count += companionMarker[ e ];
		return count;
	}

	int AdaptiveDomain::numBackground() const
	{
		return backgroundMesh.GetNE();
	}

	int AdaptiveDomain::refinements() const
	{
		return refinementCount;
	}

	int AdaptiveDomain::lastProximityAdditions() const
	{
		return proximityAdditions;
	}

	double AdaptiveDomain::largestElement() const
	{
		double largest = 0.0;
		for ( int e = 0; e < computationalMesh->GetNE(); ++e )
			largest = std::max( largest, elementDiameter( *computationalMesh, e ) );
		return largest;
	}

	double AdaptiveDomain::smallestElement() const
	{
		double smallest = 0.0;
		for ( int e = 0; e < computationalMesh->GetNE(); ++e )
		{
			double const h = elementDiameter( *computationalMesh, e );
			smallest = ( e == 0 ) ? h : std::min( smallest, h );
		}
		return smallest;
	}

	void AdaptiveDomain::refine( mfem::Array<int> const &marked )
	{
		refineBackground( marked, true );
	}

	void AdaptiveDomain::refineWithoutCompanion( mfem::Array<int> const &marked )
	{
		refineBackground( marked, false );
	}

	void AdaptiveDomain::refineBackground( mfem::Array<int> const &marked,
	                                       bool useCompanion )
	{
		// Step 3, first half: carry M from T_h to the background, which is the
		// only translation needed -- T_c^h is a set of background elements, not a
		// mesh of its own, so "the elements in T_c^h which correspond to those
		// marked on T_h" is just the parent map.
		mfem::Array<int> const &parent = computationalMesh->GetParentElementIDMap();

		std::vector<char> refineFlag( backgroundMesh.GetNE(), 0 );
		for ( int i = 0; i < marked.Size(); ++i )
		{
			int const local = marked[ i ];
			if ( local < 0 || local >= parent.Size() )
				throw std::out_of_range( "meq::AdaptiveDomain::refine: a marked index is not an element of the computational mesh" );
			refineFlag[ parent[ local ] ] = 1;
		}

		proximityAdditions = 0;
		if ( useCompanion )
		{
			// Step 3, second half, and the whole reason the companion mesh
			// exists: every element that Gamma cuts -- in the companion but not
			// inside Omega -- and that shares an edge with something in M is
			// marked too. Its children may fall inside Omega, and that is how
			// T_h grows towards Gamma instead of standing still while h_loc
			// shrinks. Computed against the ORIGINAL flags rather than in place,
			// so that one cut element cannot recruit the next along the band and
			// walk the refinement round the whole boundary.
			std::vector<char> const seed = refineFlag;
			mfem::Table const &neighbours = backgroundMesh.ElementToElementTable();
			for ( int e = 0; e < backgroundMesh.GetNE(); ++e )
			{
				if ( insideMarker[ e ] || !companionMarker[ e ] )
					continue;
				int const *row = neighbours.GetRow( e );
				for ( int i = 0; i < neighbours.RowSize( e ); ++i )
				{
					int const other = row[ i ];
					if ( other >= 0 && seed[ other ] )
					{
						refineFlag[ e ] = 1;
						proximityAdditions++;
						break;
					}
				}
			}
		}

		mfem::Array<int> list;
		for ( int e = 0; e < backgroundMesh.GetNE(); ++e )
		{
			if ( refineFlag[ e ] )
				list.Append( e );
		}
		if ( list.Size() == 0 )
			return;

		// Step 4. The default nonconforming = -1 picks conforming refinement for
		// triangles, which is what section 3.3 asks for ("avoiding the creation of
		// hanging nodes") and what the newest-vertex-bisection it names does. It
		// propagates into neighbours that were not marked, which is why the
		// background box is kept whole rather than trimmed to the companion: a
		// propagation has to have somewhere to go.
		backgroundMesh.GeneralRefinement( list );

		// Step 5.
		select();
		refinementCount++;
	}

}
