/*
 * MR-3: THE VACUUM MIRROR, THROUGH MEQ. MIRROR-PLAN.md, stage MR-3.
 *
 * No plasma. geq's two-pack mirror machine, a half-disc made by
 * tools/mesh/halfdisc.py, and the exterior coupling live -- the rung between
 * the analytic fixtures and a full machine solve, and the one that validates
 * the geometry, the mesh, the coil model and the exterior TOGETHER and without
 * the plasma. Three acceptances, one per case:
 *
 *   (a) psi against meq::coilPsi() -- the exact vacuum field inside Gamma -- at
 *       the rate the discretisation promises, k+1, over an h-ladder at three
 *       orders;
 *   (b) the mode count N, swept on MEQ's own solve until the error stops
 *       falling, which is the number MR-0 bounded and this stage measures;
 *   (c) psi against geq's own vacuum field, with the disagreement ATTRIBUTED by
 *       measuring geq's field against meq::coilPsi() directly.
 *
 *
 * WHERE THE COILS ARE DECIDES THE ROUTE
 * -------------------------------------
 *
 * FB-7's exterior conductors cannot carry a mirror's coils. MEASURED, geq's
 * packs are R in [ 0.7855, 0.8561 ], |z| in [ 1.0, 2.0 ], so the copper spans
 * rho = 1.272 .. 2.176 -- a band, not a radius. A semicircle centred on the
 * axis can only put the packs wholly on one side of it, and:
 *
 *   * OUTSIDE Gamma ( FB-7, meq::ExteriorCoilSet ) needs rho_Gamma < 1.272. The
 *     plasma reaches the throat at |z| ~ 1.5 -- rho ~ 1.51 for the flux tube
 *     through R = 0.35 at the midplane, mirror ratio about 3.5 -- and M-169's
 *     whole argument is that Gamma must enclose it. So this route cannot carry
 *     MR-4, and a vacuum fixture built on it would validate a geometry the
 *     machine cannot use.
 *   * STRADDLING it -- M-169's rho_Gamma = 1.7 cuts both packs -- is refused by
 *     both routes, ExteriorCoilSet::clearance() and
 *     ConductorField::containment() measuring to the nearest and farthest point
 *     respectively for exactly this reason.
 *   * INSIDE Gamma is the only placement that holds the plasma and the coils on
 *     one side, so Gamma goes outside the copper: rho_Gamma = 3.0, and the
 *     background mesh's arc at 3.4 so that Gamma_h is cut from it and bridged by
 *     the transfer, as the driver does.
 *
 * Inside Gamma a rectangle can be MESHED ( a domain source through
 * meq::CoilAugmentedSource, the default and the route COIL-SUBTRACTION-PLAN.md
 * section 0a-pre says must never stop working ) or SUBTRACTED ( CS-3 ). This
 * file takes the MESHED route for (a) and (b), because a subtracted vacuum has
 * psi_p = 0 identically and nothing to converge -- section 10.1 of that plan
 * says a rate there would be wrong in both directions -- while a meshed one
 * makes the coil field the solve's own answer, which is the claim. Subtraction
 * appears in (c), where it is the instrument rather than the claim.
 *
 * AND IT MOVES THE MODE COUNT FROM THE PLASMA TO THE COILS. A source inside
 * Gamma reaching rho_max puts mode n on Gamma at about ( rho_max/rho_Gamma )^n,
 * and the copper reaches 2.176 where the plasma stops near 1.6. So N is set by
 * the machine, not by the Mach number -- see (b).
 *
 *
 * THE MACHINE IS geq'S, READ FROM THE COMMITTED REFERENCE
 * -------------------------------------------------------
 *
 * tools/geq-benchmark/export_vacuum.py writes the pack geometry and current as
 * geq computed them ( make_mirror_machine_from_field( field = 4, L = 1,
 * ri = 0.7855, dz = 3, Nr = 1, Nz = 10 ) ) into meta.txt, and this file reads
 * them rather than re-deriving geq's tape-and-turn arithmetic: one machine,
 * written down once. A pack is a meq::Coil over [ ri, ro ] x [ zc - L/2,
 * zc + L/2 ] carrying the pack's whole current at uniform density -- the
 * winding as built -- where geq carries ten point filaments at the radial
 * centre. That difference is the coil model, and (c) measures it.
 */

#define BOOST_TEST_MODULE VacuumMirror
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "RoundOffFloor.hpp"

#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/PointLocator.hpp"
#include "meq/Sampler.hpp"

#include "convergence/ConvergenceHarness.hpp"

namespace
{
	/// The committed geq reference, or MEQ_GEQ_VACUUM_REFERENCE -- the same
	/// override GeqSourceComparison takes, for a regenerated one.
	std::string referenceDirectory()
	{
		char const *override = std::getenv( "MEQ_GEQ_VACUUM_REFERENCE" );
		return override ? std::string( override )
		                : std::string( "tools/geq-benchmark/reference/vacuum" );
	}

	/// geq's machine, as geq computed it. See the file comment.
	struct MirrorMachine
	{
		double mu0 = 0.0;
		double innerRadius = 0.0;
		double outerRadius = 0.0;
		double length = 0.0;
		double centreZ = 0.0;
		double packCurrent = 0.0;
		int filamentsPerPack = 0;
	};

	MirrorMachine readMachine()
	{
		std::string const path = referenceDirectory() + "/meta.txt";
		std::ifstream in( path );
		BOOST_TEST_REQUIRE( in.good(), "cannot open " << path
			<< " -- the machine is read from the committed geq reference; "
			"tools/geq-benchmark/README.md says how to regenerate it" );
		MirrorMachine m;
		std::string key;
		while ( in >> key )
		{
			if ( key == "mu0" )
				in >> m.mu0;
			else if ( key == "coil_inner_radius" )
				in >> m.innerRadius;
			else if ( key == "coil_outer_radius" )
				in >> m.outerRadius;
			else if ( key == "coil_length" )
				in >> m.length;
			else if ( key == "coil_centre_z" )
				in >> m.centreZ;
			else if ( key == "pack_current" )
				in >> m.packCurrent;
			else if ( key == "filaments_per_pack" )
				in >> m.filamentsPerPack;
			else
			{
				std::string rest;
				std::getline( in, rest );
			}
		}
		BOOST_TEST_REQUIRE( m.mu0 > 0.0 );
		BOOST_TEST_REQUIRE( m.outerRadius > m.innerRadius );
		BOOST_TEST_REQUIRE( m.length > 0.0 );
		BOOST_TEST_REQUIRE( m.packCurrent != 0.0 );
		BOOST_TEST_REQUIRE( m.filamentsPerPack > 0 );
		return m;
	}

	MirrorMachine const &machine()
	{
		static MirrorMachine const m = readMachine();
		return m;
	}

	/// One pack, at z = sign*zc. Both carry the SAME current: a mirror, not a
	/// cusp.
	meq::Coil pack( double sign )
	{
		MirrorMachine const &m = machine();
		return meq::Coil( 0.5*( m.innerRadius + m.outerRadius ), sign*m.centreZ,
		                  0.5*( m.outerRadius - m.innerRadius ), 0.5*m.length,
		                  m.packCurrent );
	}

	/// The quadrature order of every reference evaluation of coilPsi() here.
	/// 16 rather than the default 32 because the error integrals below evaluate
	/// it at up to 700,000 points: Coils.hpp's table puts 16 at 1e-9 relative
	/// INSIDE a coil and at round-off 5 cm clear of one, against a smallest
	/// solve error of 3e-8 relative over the machine region and 5e-7 globally.
	int const referenceOrder = 16;

	std::shared_ptr<meq::CoilSet> machineCoils()
	{
		auto set = std::make_shared<meq::CoilSet>( machine().mu0 );
		set->add( pack( -1.0 ) );
		set->add( pack( +1.0 ) );
		set->setQuadratureOrder( referenceOrder );
		return set;
	}

	/// Gamma, and the background mesh's own arc outside it. See the file
	/// comment for why Gamma is outside the copper.
	double const rhoGamma = 3.0;
	double const meshRadius = 3.4;

	/// The mode count every case uses unless it is sweeping it -- (b)'s
	/// plateau with room to spare, so the truncation stays below the
	/// discretisation error at every mesh in (a).
	int const defaultModes = 32;

	/// geq's grid box, where MR-4's plasma lives. The error over it is the one
	/// asserted at k+1 for every order: it is 0.15 m clear of the copper, so the
	/// conductor corners' r^2 log r is not in it.
	bool inMachineRegion( double radius, double z )
	{
		return radius <= 0.63 && std::abs( z ) <= 1.5;
	}

	/**
	 * The half-disc, made by halfdisc.py once per process into the build tree.
	 *
	 * The packs are MESHED TO -- FB-2's aligned rates against cut ones are
	 * 1.99/2.88/3.01 against 1.33/1.27/1.09 -- and the machine region is refined
	 * as a plasma box would be, so this is the mesh MR-4 starts from with the
	 * plasma still absent. Coarse vacuum out to the arc, because Gamma is far
	 * out and the field there is smooth.
	 *
	 * The ladder below REFINES this mesh rather than re-running halfdisc.py at
	 * smaller sizes, for DivertedMachine.hpp's reason: a re-meshed geometry
	 * moves which elements the conductors occupy, a second variable in a study
	 * of the first, and uniform refinement carries the coil edges down exactly.
	 */
	std::string const &meshFile()
	{
		static std::string const file = []()
		{
			MirrorMachine const &m = machine();
			std::string const out = std::string( MEQ_TEST_OUTPUT_DIR )
			                        + "/VacuumMirror.msh";
			double const width = m.outerRadius - m.innerRadius;
			char buffer[ 2048 ];
			std::snprintf( buffer, sizeof buffer,
				"python3 %s --rho %.17g --size 0.4"
				" --coil %.17g %.17g %.17g %.17g"
				" --coil %.17g %.17g %.17g %.17g --coil-size 0.07"
				" --plasma 0.0 -1.7 0.75 3.4 --plasma-size 0.1"
				" --check -o %s > /dev/null",
				MEQ_HALFDISC_PATH, meshRadius,
				m.innerRadius, -m.centreZ - 0.5*m.length, width, m.length,
				m.innerRadius, m.centreZ - 0.5*m.length, width, m.length,
				out.c_str() );
			BOOST_TEST_REQUIRE( std::system( buffer ) == 0,
				"halfdisc.py failed, so there is no mesh. It needs python3 with "
				"gmsh; the command was: " << buffer );
			return out;
		}();
		return file;
	}

	/// Negative inside Gamma.
	double semicircle( mfem::Vector const &x )
	{
		return std::hypot( x( 0 ), x( 1 ) ) - rhoGamma;
	}

	/// D_h cut from the background at Gamma. The background owns the parent of
	/// the SubMesh and so lives as long as it does.
	struct Domain
	{
		std::unique_ptr<mfem::Mesh> background;
		std::unique_ptr<mfem::SubMesh> sub;
		std::unique_ptr<mfem::VertexConePath> path;
		mfem::Array<int> gammaHMarker;
		int gammaH = 0;
		double h = 0.0;
	};

	std::unique_ptr<Domain> buildDomain( int refinements )
	{
		auto d = std::make_unique<Domain>();
		d->background = std::make_unique<mfem::Mesh>( meshFile().c_str(), 1, 1 );
		for ( int level = 0; level < refinements; ++level )
			d->background->UniformRefinement();

		// A MESH READ FROM A FILE HAS ITS h MEASURED, and a zero search length
		// is mfem::VertexConePath aborting on the first vertex of Gamma_h.
		for ( int e = 0; e < d->background->GetNE(); ++e )
			d->h = std::max( d->h, meq::elementDiameter( *d->background, e ) );

		mfem::Array<int> marker;
		BOOST_TEST_REQUIRE( mfem::MarkLevelSetSubdomain(
			*d->background, semicircle, meq::markingOffset( 1.0, 1.5 ), marker,
			1 ) > 0, "the subdomain inside Gamma is empty" );
		int const parentBoundaryMax = d->background->bdr_attributes.Max();
		for ( int e = 0; e < d->background->GetNE(); ++e )
			d->background->SetAttribute( e, marker[ e ] ? 1 : 2 );
		d->background->SetAttributes();

		mfem::Array<int> domainAttribute( 1 );
		domainAttribute[ 0 ] = 1;
		d->sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( *d->background, domainAttribute ) );
		d->gammaH = d->sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( d->gammaH > parentBoundaryMax,
			"no boundary was generated cutting at Gamma, so Gamma_h names the "
			"background's own arc and nothing is transferred" );
		d->gammaHMarker.SetSize( d->gammaH );
		d->gammaHMarker = 0;
		d->gammaHMarker[ d->gammaH - 1 ] = 1;
		d->path = std::make_unique<mfem::VertexConePath>(
			*d->sub, d->gammaH, semicircle, 6.0*d->h );
		return d;
	}

	/// No plasma. A meq::Source and not a Coefficient because the coupled path
	/// is NPC and refuses a Coefficient.
	struct EmptyInterior : public meq::Source
	{
		double f( double, double, double ) const override { return 0.0; }
		double dFdPsi( double, double, double ) const override { return 0.0; }
	};

	/// Everything a solve aliases, kept alive with it: the solver holds the
	/// DtN and the source by reference.
	struct VacuumSolve
	{
		std::shared_ptr<meq::CoilSet> coils;
		std::unique_ptr<meq::Source> source;
		std::unique_ptr<meq::ConductorField> subtracted;
		std::unique_ptr<meq::ExteriorDtN> dtn;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<meq::GradShafranovSolver> solver;
	};

	/// The meshed route: the packs are a domain source, the [[coils]] route.
	/// @param coupled false drops the exterior coupling, leaving psi = 0 on
	///        Gamma_h -- the control in (a).
	std::unique_ptr<VacuumSolve> solveMeshed( Domain &d, int order, int modes,
	                                          std::shared_ptr<meq::CoilSet> coils,
	                                          bool coupled = true )
	{
		auto s = std::make_unique<VacuumSolve>();
		s->coils = std::move( coils );
		s->source = std::make_unique<meq::CoilAugmentedSource>(
			std::make_shared<EmptyInterior>(), s->coils );
		s->dtn = std::make_unique<meq::ExteriorDtN>( 0.0, rhoGamma, modes );
		s->solver = std::make_unique<meq::GradShafranovSolver>( *d.sub, order );
		s->solver->setSource( *s->source );
		s->solver->setBoundaryData( s->zero );
		s->solver->setExtension( *d.path, d.gammaHMarker );
		if ( coupled )
			s->solver->setExteriorCoupling( *s->dtn );
		s->solver->solve();
		return s;
	}

	/// The subtracted route, CS-3: the conductors are taken off the right-hand
	/// side and handed to the solver as psi_c, so what it solves for is the
	/// remainder psi_p and the interior equation has no source at all.
	std::unique_ptr<VacuumSolve> solveSubtracted( Domain &d, int order, int modes,
	                                              meq::ConductorField const &field )
	{
		auto s = std::make_unique<VacuumSolve>();
		s->source = std::make_unique<EmptyInterior>();
		s->subtracted = std::make_unique<meq::ConductorField>( field );
		s->dtn = std::make_unique<meq::ExteriorDtN>( 0.0, rhoGamma, modes );
		s->solver = std::make_unique<meq::GradShafranovSolver>( *d.sub, order );
		s->solver->setSource( *s->source );
		s->solver->setBoundaryData( s->zero );
		s->solver->setExtension( *d.path, d.gammaHMarker );
		s->solver->setConductorField( *s->subtracted );
		s->solver->setExteriorCoupling( *s->dtn );
		s->solver->solve();
		return s;
	}

	struct Errors
	{
		double global = 0.0;
		double region = 0.0;
		double globalNorm = 0.0;
		double regionNorm = 0.0;
	};

	/**
	 * L2( psi_h - psi_coil ) over D_h and over the machine region, from ONE
	 * evaluation of the reference.
	 *
	 * The reference is what costs -- coilPsi() at every quadrature point, a
	 * Gauss rule over each pack per point -- so it is tabulated once per mesh
	 * and order and then charged against as many solves as share them, which
	 * is what lets (b) sweep N for the price of one table. The rule is
	 * potentialError()'s, 2k + 4. The element loop stays serial, with a
	 * caller-allocated transformation rather than the mesh's shared scratch,
	 * and only the closed-form evaluation is threaded.
	 */
	class ReferenceTable
	{
		public:
			ReferenceTable( mfem::Mesh &meshIn, int orderIn,
			                meq::CoilSet const &reference )
				: mesh( meshIn ), order( orderIn )
			{
				mfem::IsoparametricTransformation transformation;
				mfem::Vector x( 2 );
				mfem::Vector centre( 2 );
				for ( int e = 0; e < mesh.GetNE(); ++e )
				{
					mesh.GetElementTransformation( e, &transformation );
					mesh.GetElementCenter( e, centre );
					bool const inside = inMachineRegion( centre( 0 ), centre( 1 ) );
					mfem::IntegrationRule const &ir = rule( e );
					for ( int q = 0; q < ir.GetNPoints(); ++q )
					{
						mfem::IntegrationPoint const &ip = ir.IntPoint( q );
						transformation.SetIntPoint( &ip );
						transformation.Transform( ip, x );
						radius.push_back( x( 0 ) );
						height.push_back( x( 1 ) );
						weight.push_back( ip.weight*transformation.Weight() );
						region.push_back( inside ? 1 : 0 );
					}
				}
				exact.resize( radius.size() );
				long const n = static_cast<long>( radius.size() );
				#pragma omp parallel for schedule( dynamic, 256 )
				for ( long i = 0; i < n; ++i )
					exact[ static_cast<std::size_t>( i ) ] = reference.psi(
						radius[ static_cast<std::size_t>( i ) ],
						height[ static_cast<std::size_t>( i ) ] );
			}

			/// @param psi a potential on the mesh this table was built on, at
			///        the order it was built for.
			Errors errors( mfem::GridFunction const &psi ) const
			{
				BOOST_TEST_REQUIRE( psi.FESpace()->GetMesh() == &mesh );
				Errors out;
				mfem::IsoparametricTransformation transformation;
				std::size_t i = 0;
				for ( int e = 0; e < mesh.GetNE(); ++e )
				{
					mesh.GetElementTransformation( e, &transformation );
					mfem::IntegrationRule const &ir = rule( e );
					for ( int q = 0; q < ir.GetNPoints(); ++q, ++i )
					{
						mfem::IntegrationPoint const &ip = ir.IntPoint( q );
						transformation.SetIntPoint( &ip );
						double const difference =
							psi.GetValue( transformation, ip ) - exact[ i ];
						double const d2 = weight[ i ]*difference*difference;
						double const n2 = weight[ i ]*exact[ i ]*exact[ i ];
						out.global += d2;
						out.globalNorm += n2;
						if ( region[ i ] )
						{
							out.region += d2;
							out.regionNorm += n2;
						}
					}
				}
				out.global = std::sqrt( out.global );
				out.region = std::sqrt( out.region );
				out.globalNorm = std::sqrt( out.globalNorm );
				out.regionNorm = std::sqrt( out.regionNorm );
				return out;
			}

		private:
			mfem::IntegrationRule const &rule( int e ) const
			{
				return mfem::IntRules.Get( mesh.GetElementBaseGeometry( e ),
				                           2*order + 4 );
			}

			mfem::Mesh &mesh;
			int order;
			std::vector<double> radius;
			std::vector<double> height;
			std::vector<double> weight;
			std::vector<double> exact;
			std::vector<char> region;
	};

	/**
	 * The relative tail of psi_coil's own trace on Gamma past N retained
	 * modes, in the operator's weight dGamma/R = dmu/( 1 - mu^2 ) -- what M-169
	 * computed for the plasma, computed for the copper. It is the best any
	 * N-mode coupling can do ON Gamma, and entry N is the tail after degrees
	 * 2 .. N + 1.
	 *
	 * **NOT ExteriorDtN::coefficients(), WHICH ALIASES PAST ABOUT THIRTY
	 * MODES.** That method spends a fixed 30-point Gauss rule, exact for mode
	 * pairs with m + n <= 61, and its own comment says nobody will use forty
	 * modes. Measured here at 96 it returns a tail of 0.82 flat out to N = 48,
	 * which is the quadrature and not the trace. Nothing in the solve calls it
	 * -- the coupling's coefficients come out of the bordered Newton -- so this
	 * is a diagnostic's limit and not a solver's; a 400-point rule in mu is
	 * spent here instead, which resolves the highest mode with room to spare.
	 * The tail is summed DOWN from the last mode rather than subtracted from
	 * the norm, which would floor it at the square root of round-off.
	 */
	std::vector<double> traceTail( meq::CoilSet const &coils, int maxModes )
	{
		meq::ExteriorDtN const wide( 0.0, rhoGamma, maxModes );
		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( mfem::Geometry::SEGMENT, 799 );

		std::vector<double> a( static_cast<std::size_t>( maxModes ), 0.0 );
		double norm2 = 0.0;
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			double const mu = 2.0*ip.x - 1.0;
			double const w = 2.0*ip.weight/( ( 1.0 - mu )*( 1.0 + mu ) );
			double const radius = rhoGamma*std::sqrt( ( 1.0 - mu )*( 1.0 + mu ) );
			double const z = rhoGamma*mu;
			double const psi = coils.psi( radius, z );
			norm2 += w*psi*psi;
			for ( int i = 0; i < maxModes; ++i )
				a[ static_cast<std::size_t>( i ) ] +=
					w*psi*wide.basis( wide.firstMode() + i, radius, z );
		}

		std::vector<double> tail( static_cast<std::size_t>( maxModes ) + 1, 0.0 );
		double rest = 0.0;
		for ( int i = maxModes; i >= 0; --i )
		{
			tail[ static_cast<std::size_t>( i ) ] = std::sqrt( rest/norm2 );
			if ( i > 0 )
			{
				int const n = wide.firstMode() + i - 1;
				double const coefficient = a[ static_cast<std::size_t>( i - 1 ) ]
				                           /wide.mass( n );
				rest += coefficient*coefficient*wide.mass( n );
			}
		}
		return tail;
	}

	/// geq's vacuum field on its own grid, R z psi per row.
	struct GeqPoint
	{
		double radius;
		double z;
		double psi;
	};

	std::vector<GeqPoint> readGeqPoints()
	{
		std::string const path = referenceDirectory() + "/points.dat";
		std::ifstream in( path );
		BOOST_TEST_REQUIRE( in.good(), "cannot open " << path );
		std::vector<GeqPoint> points;
		std::string line;
		while ( std::getline( in, line ) )
		{
			if ( line.empty() || line[ 0 ] == '#' )
				continue;
			std::istringstream row( line );
			GeqPoint p{};
			BOOST_TEST_REQUIRE( static_cast<bool>( row >> p.radius >> p.z >> p.psi ),
			                    "malformed row in " << path << ": " << line );
			points.push_back( p );
		}
		BOOST_TEST_REQUIRE( !points.empty(), path << " has no points" );
		return points;
	}

	double rms( std::vector<double> const &v )
	{
		double sum = 0.0;
		for ( double x : v )
			sum += x*x;
		return std::sqrt( sum/static_cast<double>( v.size() ) );
	}
}


/*
 * (a) MEQ'S psi AGAINST meq::coilPsi(), AT k+1.
 *
 * With no plasma and every conductor inside Gamma, the exact answer is the
 * packs' own field: Delta* psi_coil = -mu0 R J inside the copper and 0 outside,
 * psi_coil = 0 on the axis, and its exterior a decaying Gegenbauer series that
 * the coupling represents to the mode truncation. So psi_h converges on
 * coilPsi() and the order is the discretisation's.
 *
 * TWO ROUTES TO ONE FIELD, AND THAT IS WHAT MAKES IT A TEST. The solve sees the
 * packs only as a current density, CoilSet::f() = mu0 R I/area, sampled at
 * quadrature points of the elements meshed to them; the reference integrates
 * the filament Green's function over each rectangle. A factor, a sign or a
 * misplaced conductor in either is a disagreement that does not fall with h:
 * each of M-260's four mutations -- one pack moved by its own half-width, one
 * pack's current reversed, 0.1% on CoilSet::f(), 1% on the DtN's block entry --
 * turns every order red here, the smallest of them flattening the k = 3 rate
 * to 0.03.
 *
 * WHERE IT IS ASSERTED, AND WHY TWO NORMS. Over the MACHINE REGION -- geq's own
 * box, R <= 0.63 and |z| <= 1.5, where MR-4's plasma lives and 0.15 m clear of
 * the copper -- with the two-tier rule CLAUDE_HDGGS.md gives for the extension
 * path: 0.30 of slack on every pair and 0.15 over the sequence. Over the WHOLE
 * of D_h on the finest pair only: the background is 0.4 m out to the arc, so
 * level 0 resolves the vacuum near Gamma with three elements across a metre,
 * and k = 1's first pair there reads 1.27 before reaching 2.14. Both norms reach
 * k+1 at k = 3 as well, which says the r^2 log r at the packs' corners -- which
 * caps a cut conductor's rate -- is not limiting a conductor meshed to.
 */
BOOST_AUTO_TEST_CASE( theVacuumMirrorConvergesOnTheCoilFieldAtKPlusOne )
{
	std::shared_ptr<meq::CoilSet> const reference = machineCoils();

	// THE GEOMETRY DECISION, ASKED OF THE TYPES THAT OWN IT. Every pack wholly
	// inside Gamma, measured to its FARTHEST point -- the precondition the
	// subtracted route checks and the meshed route needs just as much, since
	// a source beyond Gamma is a current the exterior series cannot hold. And
	// at M-169's rho_Gamma = 1.7 the copper straddles Gamma, so NEITHER route
	// admits it: that radius was chosen for the plasma alone.
	meq::ConductorField copper( machine().mu0 );
	meq::ExteriorCoilSet outside( machine().mu0 );
	for ( double sign : { -1.0, 1.0 } )
	{
		copper.add( pack( sign ) );
		outside.add( pack( sign ) );
	}
	double const containment = copper.containment( 0.0, rhoGamma );
	BOOST_TEST_REQUIRE( containment > 0.0,
		"a pack reaches beyond Gamma = " << rhoGamma << " by " << -containment
		<< " m, so its field is not the exterior series' to represent" );
	BOOST_TEST( copper.containment( 0.0, 1.7 ) < 0.0 );
	BOOST_TEST( outside.clearance( 0.0, 1.7 ) < 0.0 );

	std::printf( "\n  MR-3 (a): geq's mirror, no plasma, packs MESHED inside "
	             "Gamma = %.2f ( %.4f m clear ), N = %d\n",
	             rhoGamma, containment, defaultModes );
	std::printf( "    pack R [%.6f, %.6f], z +-[%.3f, %.3f], I = %.6e A each\n\n",
	             machine().innerRadius, machine().outerRadius,
	             machine().centreZ - 0.5*machine().length,
	             machine().centreZ + 0.5*machine().length, machine().packCurrent );
	std::printf( "    %2s %5s %7s %8s %13s %7s %13s %7s\n", "k", "level", "NE",
	             "h", "L2 D_h", "rate", "L2 machine", "rate" );

	int const levels = 3;
	int const orders = 3;
	std::vector<std::vector<double>> global( orders );
	std::vector<std::vector<double>> region( orders );
	double normGlobal = 0.0;
	double normRegion = 0.0;
	double controlRatio = 0.0;

	for ( int level = 0; level < levels; ++level )
	{
		std::unique_ptr<Domain> const d = buildDomain( level );
		for ( int order = 1; order <= orders; ++order )
		{
			std::unique_ptr<VacuumSolve> const s =
				solveMeshed( *d, order, defaultModes, machineCoils() );
			ReferenceTable const table( *d->sub, order, *reference );
			Errors const e = table.errors( s->solver->potential() );
			std::size_t const k = static_cast<std::size_t>( order - 1 );
			global[ k ].push_back( e.global );
			region[ k ].push_back( e.region );
			normGlobal = e.globalNorm;
			normRegion = e.regionNorm;

			double gRate = std::numeric_limits<double>::quiet_NaN();
			double rRate = std::numeric_limits<double>::quiet_NaN();
			if ( level > 0 )
			{
				gRate = meq::tests::rate( global[ k ][ level - 1 ], e.global, 2.0 );
				rRate = meq::tests::rate( region[ k ][ level - 1 ], e.region, 2.0 );
			}
			std::printf( "    %2d %5d %7d %8.4f %13.4e %7.3f %13.4e %7.3f\n",
			             order, level, d->sub->GetNE(), d->h, e.global, gRate,
			             e.region, rRate );
			std::fflush( stdout );

			// AFFINE, SO ONE STEP REACHES THE FLOOR. The coils enter as data
			// and the interior source is empty, so dF/dpsi = 0 and the
			// residual is affine in ( x, a ). Asserted on the residual history
			// rather than the count, which is M-219's lesson: whether a
			// second step is taken to notice is a property of the target.
			BOOST_TEST( meq::tests::oneExactStepReachedTheFloor(
			                s->solver->newtonResiduals() ),
				"k = " << order << ", level " << level << ": "
				<< meq::tests::stepReductionNotice( s->solver->newtonResiduals() ) );

			// THE INSTRUMENT, ONCE: the table's own L2 against the solver's.
			// It integrates the same rule with the same reference, so a
			// disagreement is the table -- a transposed point, a dropped
			// weight -- and every rate below would inherit it.
			if ( level == 0 && order == 1 )
			{
				mfem::FunctionCoefficient exact( [ &reference ]( mfem::Vector const &x )
				{
					return reference->psi( x( 0 ), x( 1 ) );
				} );
				double const direct = s->solver->potentialError( exact );
				BOOST_TEST_REQUIRE( std::abs( direct - e.global ) <= 1.0e-10*direct,
					"the error table reads " << e.global << " where "
					"potentialError() reads " << direct );
			}

			// THE CONTROL, ONCE: the same solve with the coupling REMOVED, so
			// Gamma_h carries psi = 0. A coupling that did nothing would leave
			// the packs' field truncated at the boundary and every rate above
			// would still be a rate -- towards a different function.
			if ( level == 0 && order == 2 )
			{
				std::unique_ptr<VacuumSolve> const bare =
					solveMeshed( *d, order, defaultModes, machineCoils(), false );
				Errors const c = table.errors( bare->solver->potential() );
				controlRatio = c.region/e.region;
				std::printf( "       control, coupling removed: L2 D_h %.4e, "
				             "machine %.4e ( %.0fx this row )\n",
				             c.global, c.region, controlRatio );
				std::fflush( stdout );
			}
		}
	}

	std::printf( "\n    |psi_coil|: %.4e over D_h, %.4e over the machine region\n",
	             normGlobal, normRegion );

	for ( int order = 1; order <= orders; ++order )
	{
		std::size_t const k = static_cast<std::size_t>( order - 1 );
		double const wanted = order + 1.0;
		double const overall = meq::tests::rate( region[ k ].front(),
		                                         region[ k ].back(),
		                                         std::pow( 2.0, levels - 1 ) );
		double const finestGlobal = meq::tests::rate(
			global[ k ][ levels - 2 ], global[ k ][ levels - 1 ], 2.0 );
		std::printf( "    k = %d: machine region %.3f over the sequence, "
		             "D_h %.3f on the finest pair, k+1 = %.0f\n",
		             order, overall, finestGlobal, wanted );

		if ( !meq::tests::rateOfOrderIsMeasurable( wanted ) )
		{
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice(
				order, "psi against the packs' own field", wanted ) );
			continue;
		}
		for ( int level = 1; level < levels; ++level )
		{
			double const pair = meq::tests::rate( region[ k ][ level - 1 ],
			                                      region[ k ][ level ], 2.0 );
			BOOST_TEST( pair > wanted - 0.30,
				"k = " << order << ": psi converges on the coil field at " << pair
				<< " over the machine region between levels " << level - 1
				<< " and " << level << ", against k+1 = " << wanted << ". The "
				"exact answer is coilPsi() and the region is clear of the "
				"copper, so a short rate is the source term, the coupling or the "
				"transfer -- not a singularity." );
		}
		BOOST_TEST( overall > wanted - 0.15,
			"k = " << order << ": the machine-region rate over the whole ladder "
			"is " << overall << " against k+1 = " << wanted );
		BOOST_TEST( finestGlobal > wanted - 0.30,
			"k = " << order << ": over the whole of D_h the finest pair reads "
			<< finestGlobal << " against k+1 = " << wanted );
	}

	BOOST_TEST( controlRatio > 50.0,
		"removing the exterior coupling moved the machine-region error by only "
		<< controlRatio << "x, so the coupling is not what is carrying the "
		"packs' field across Gamma" );
}

/*
 * (b) THE MODE COUNT, SWEPT ON MEQ'S OWN SOLVE -- AND IT IS THE COPPER'S.
 *
 * MR-0 bounded N from the PLASMA's current: M-169 needed 32 modes at
 * rho_Gamma = 1.7 for 2.65e-05 at M = 2. With the packs inside Gamma the
 * plasma is not what sets it. The copper reaches rho = 2.176, so the packs'
 * trace on Gamma = 3.0 carries mode n at about ( 2.176/3.0 )^n = 0.725^n, and
 * that is the spectrum the coupling must hold. The trace tail column is that
 * spectrum, computed here from coilPsi() alone with no solve.
 *
 * WHAT IS ASSERTED. That the error over D_h falls by more than two decades and
 * then STOPS -- the last two rows agree to 2% -- because past that point the
 * truncation is below the discretisation error and N buys nothing; that
 * defaultModes is already on the plateau, so (a) and (c) are not measuring the
 * truncation; and that the machine region feels the truncation too, at N = 4,
 * but plateaus far earlier than D_h does, the modes past N being harmonic
 * fields that decay inward from Gamma as ( rho/rho_Gamma )^n.
 */
BOOST_AUTO_TEST_CASE( theModeCountIsSetByTheCopperAndNotByThePlasma )
{
	int const order = 3;
	int const level = 1;
	std::vector<int> const sweep = { 4, 8, 12, 16, 20, 24, 32, 48 };

	std::shared_ptr<meq::CoilSet> const reference = machineCoils();
	std::vector<double> const tail = traceTail( *reference, 96 );

	std::unique_ptr<Domain> const d = buildDomain( level );
	ReferenceTable const table( *d->sub, order, *reference );

	std::printf( "\n  MR-3 (b): the mode count, k = %d, level %d ( %d elements ), "
	             "Gamma = %.2f\n", order, level, d->sub->GetNE(), rhoGamma );
	std::printf( "    copper to rho = %.4f, so mode n falls as about %.3f^n on Gamma\n\n",
	             std::hypot( machine().outerRadius,
	                         machine().centreZ + 0.5*machine().length ),
	             std::hypot( machine().outerRadius,
	                         machine().centreZ + 0.5*machine().length )/rhoGamma );
	std::printf( "    %3s %13s %13s %13s %13s %13s\n", "N", "trace tail",
	             "L2 D_h", "relative", "L2 machine", "relative" );

	std::vector<double> global;
	std::vector<double> region;
	for ( int modes : sweep )
	{
		std::unique_ptr<VacuumSolve> const s =
			solveMeshed( *d, order, modes, machineCoils() );
		Errors const e = table.errors( s->solver->potential() );
		global.push_back( e.global );
		region.push_back( e.region );
		std::printf( "    %3d %13.4e %13.4e %13.4e %13.4e %13.4e\n", modes,
		             tail[ static_cast<std::size_t>( modes ) ], e.global,
		             e.global/e.globalNorm, e.region, e.region/e.regionNorm );
		std::fflush( stdout );
	}

	std::size_t const last = sweep.size() - 1;
	std::size_t plateau = last;
	for ( std::size_t i = 0; i <= last; ++i )
		if ( global[ i ] <= 1.05*global[ last ] )
		{
			plateau = i;
			break;
		}
	std::printf( "\n    the error over D_h is within 5%% of N = %d's from N = %d "
	             "on; the trace tail there is %.3e\n\n", sweep[ last ],
	             sweep[ plateau ], tail[ static_cast<std::size_t>( sweep[ plateau ] ) ] );

	auto at = [ &sweep ]( std::vector<double> const &v, int modes )
	{
		std::size_t const i = static_cast<std::size_t>(
			std::find( sweep.begin(), sweep.end(), modes ) - sweep.begin() );
		BOOST_TEST_REQUIRE( i < v.size() );
		return v[ i ];
	};

	// IT FALLS: the modes are what carries the field across Gamma.
	BOOST_TEST( global.front() > 100.0*at( global, defaultModes ),
		"the error over D_h falls only from " << global.front() << " at N = "
		<< sweep.front() << " to " << at( global, defaultModes ) << " at N = "
		<< defaultModes << " -- the coupling is not resolving the packs' trace" );
	BOOST_TEST( region.front() > 100.0*at( region, defaultModes ),
		"the machine region does not feel the truncation at N = "
		<< sweep.front() << ": " << region.front() << " against "
		<< at( region, defaultModes ) );

	// AND IT STOPS: past the plateau N buys nothing, so the error left is the
	// discretisation's. A rising last row would be the coupling degrading with
	// N; a still-falling one would mean (a) is reading the truncation.
	BOOST_TEST( std::abs( global[ last ]/at( global, defaultModes ) - 1.0 ) < 0.02,
		"between N = " << defaultModes << " and N = " << sweep[ last ]
		<< " the error over D_h moved from " << at( global, defaultModes )
		<< " to " << global[ last ] << ", so it has not stopped falling" );
	BOOST_TEST( sweep[ plateau ] <= defaultModes,
		"the plateau starts at N = " << sweep[ plateau ] << ", beyond the "
		<< defaultModes << " every other case here runs at" );

	// AND THE MACHINE REGION PLATEAUS FIRST, which is the statement MR-4 can
	// spend: the plasma sees the truncation through ( rho/rho_Gamma )^n.
	BOOST_TEST( region[ 1 ] < 1.2*region[ last ],
		"the machine region is still " << region[ 1 ]/region[ last ]
		<< "x its plateau at N = " << sweep[ 1 ] );
}

/*
 * (c) AGAINST geq'S OWN VACUUM FIELD, WITH THE DISAGREEMENT ATTRIBUTED.
 *
 * tools/geq-benchmark/reference/vacuum/points.dat is geq's field on every
 * other node of its 65 x 65 grid. With no plasma that field is a closed-form
 * sum over geq's twenty point filaments -- export_vacuum.py asserts so to
 * 5e-16 -- so MEQ is being compared with a coil MODEL, and four numbers say
 * how much of the difference is whose:
 *
 *   1. geq against meq::filamentPsi() over meq::filamentStack( pack, 1, 10 ):
 *      the same twenty filaments in MEQ's own kernel. Round-off, which settles
 *      every convention at once -- amperes per filament, turns, the sign of
 *      each pack, weber per radian against weber.
 *   2. geq against coilPsi(): THE MODEL FLOOR. Ten filaments at the radial
 *      centre are the midpoint rule for the rectangle, so this is geq's
 *      quadrature of its own winding, measured with no solve anywhere -- the
 *      reference measured against itself, M-167's lesson.
 *   3. MEQ, packs meshed, against geq: equal to 2, because
 *   4. MEQ, packs meshed, against coilPsi() is hundreds of times smaller.
 *
 * AND THEN THE MODEL IS REMOVED. Handing MEQ geq's own twenty filaments,
 * SUBTRACTED ( CS-3; a filament has no other route ), the floor is gone and
 * what is left is the mode truncation of a trace with a logarithm 0.88 m from
 * Gamma. So the residual is attributed in both directions: MEQ reproduces geq
 * when given geq's machine, and differs from it by exactly geq's model when
 * given the winding.
 *
 * Points outside D_h are DROPPED, and so is any inside a conductor; geq's box
 * is wholly inside both Gamma and clear of the copper, so none is, and that is
 * asserted rather than assumed.
 */
BOOST_AUTO_TEST_CASE( theVacuumMirrorAgreesWithGeqToGeqsOwnCoilModel )
{
	int const order = 3;
	int const level = 1;

	std::shared_ptr<meq::CoilSet> const reference = machineCoils();
	std::vector<GeqPoint> const points = readGeqPoints();
	int const count = static_cast<int>( points.size() );

	meq::ConductorField geqModel( machine().mu0 );
	for ( double sign : { -1.0, 1.0 } )
		for ( meq::CurrentFilament const &filament :
		      meq::filamentStack( pack( sign ), 1, machine().filamentsPerPack ) )
			geqModel.add( filament );
	BOOST_TEST_REQUIRE( geqModel.containment( 0.0, rhoGamma ) > 0.0 );

	std::unique_ptr<Domain> const d = buildDomain( level );
	std::unique_ptr<VacuumSolve> const meshed =
		solveMeshed( *d, order, defaultModes, machineCoils() );
	std::unique_ptr<VacuumSolve> const subtracted =
		solveSubtracted( *d, order, defaultModes, geqModel );

	mfem::Vector positions( 2*count );
	for ( int i = 0; i < count; ++i )
	{
		positions( i ) = points[ static_cast<std::size_t>( i ) ].radius;
		positions( count + i ) = points[ static_cast<std::size_t>( i ) ].z;
	}
	meq::PointLocator locator( *d->sub );
	locator.locate( positions );
	mfem::Vector meshedValues;
	mfem::Vector remainder;
	locator.interpolate( meshed->solver->potential(), meshedValues );
	locator.interpolate( subtracted->solver->potential(), remainder );

	std::vector<double> self;
	std::vector<double> model;
	std::vector<double> meqGeq;
	std::vector<double> meqExact;
	std::vector<double> matched;
	std::vector<double> scale;
	int dropped = 0;
	for ( int i = 0; i < count; ++i )
	{
		GeqPoint const &p = points[ static_cast<std::size_t>( i ) ];
		if ( locator.codes()[ i ] == 2 || reference->indexContaining( p.radius, p.z ) >= 0 )
		{
			++dropped;
			continue;
		}
		double const filaments = geqModel.psi( p.radius, p.z );
		double const winding = reference->psi( p.radius, p.z );
		self.push_back( filaments - p.psi );
		model.push_back( winding - p.psi );
		meqGeq.push_back( meshedValues( i ) - p.psi );
		meqExact.push_back( meshedValues( i ) - winding );
		matched.push_back( remainder( i ) + filaments - p.psi );
		scale.push_back( p.psi );
	}
	BOOST_TEST_REQUIRE( dropped == 0,
		dropped << " of geq's " << count << " points are outside D_h or inside a "
		"pack, where its box should be inside Gamma and clear of the copper" );

	double const norm = rms( scale );
	double const conventions = rms( self )/norm;
	double const floor = rms( model )/norm;
	double const meqAgainstGeq = rms( meqGeq )/norm;
	double const meqAgainstExact = rms( meqExact )/norm;
	double const geqMachine = rms( matched )/norm;

	std::printf( "\n  MR-3 (c): against geq's vacuum field, %d points, k = %d, "
	             "level %d, N = %d -- rms relative to rms psi = %.4e\n\n",
	             count - dropped, order, level, defaultModes, norm );
	std::printf( "    1. geq's filaments in MEQ's kernel  vs geq   %11.3e\n",
	             conventions );
	std::printf( "    2. the winding, coilPsi()          vs geq   %11.3e  "
	             "<- geq's coil model\n", floor );
	std::printf( "    3. MEQ, packs meshed               vs geq   %11.3e\n",
	             meqAgainstGeq );
	std::printf( "    4. MEQ, packs meshed               vs 2.    %11.3e\n",
	             meqAgainstExact );
	std::printf( "    5. MEQ, geq's filaments subtracted vs geq   %11.3e\n\n",
	             geqMachine );
	std::fflush( stdout );

	// What this build can resolve, relative: the double build never reaches it.
	double const resolvable = meq::tests::roundOffFloor( norm )/norm;

	// 1. THE CONVENTIONS. Two implementations of one Green's function, scipy's
	// elliptic integrals against Carlson's, so round-off and not a tolerance.
	BOOST_TEST( conventions < 1.0e-12,
		"geq's twenty filaments evaluated by MEQ's own kernel differ from geq "
		"by " << conventions << " relative -- a convention is wrong: current per "
		"filament, turns, a pack's sign, or Wb against Wb/rad" );

	// 2 AND 3. THE DISAGREEMENT IS geq'S COIL MODEL.
	BOOST_TEST( floor > 1.0e-5,
		"geq's filaments and the winding agree to " << floor << ", so there is "
		"no model floor to attribute and this case has stopped measuring one" );
	BOOST_TEST( std::abs( meqAgainstGeq/floor - 1.0 ) < 0.05,
		"MEQ differs from geq by " << meqAgainstGeq << " where geq's own coil "
		"model differs from the winding by " << floor << " -- the disagreement "
		"is not the model, so some of it is MEQ's" );

	// 4. AND IT IS NOT MEQ'S.
	BOOST_TEST( meqAgainstExact < std::max( 0.02*floor, resolvable ),
		"MEQ's own error at geq's points is " << meqAgainstExact << ", not small "
		"against the model floor " << floor );

	// 5. AND IT GOES AWAY WHEN MEQ IS GIVEN geq'S MACHINE.
	BOOST_TEST( geqMachine < std::max( 1.0e-3*floor, resolvable ),
		"handed geq's own filaments, MEQ still differs from geq by "
		<< geqMachine << " against a model floor of " << floor );
}
