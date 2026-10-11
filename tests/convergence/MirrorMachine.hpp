#ifndef MEQ_TESTS_MIRROR_MACHINE_HPP
#define MEQ_TESTS_MIRROR_MACHINE_HPP

/*
 * geq'S MIRROR MACHINE, ASSEMBLED ONCE FOR EVERY CASE THAT NEEDS ONE.
 * MIRROR-PLAN.md MR-3 and MR-4.
 *
 * The pack geometry and current, the half-disc mesh made by
 * tools/mesh/halfdisc.py, D_h cut from it at Gamma, the coil field as a
 * reference, geq's own twenty filaments, and the vacuum solve on the meshed
 * route. VacuumMirror is the acceptance of all of it without a plasma, and
 * RotatingMirror solves the plasma on top of it.
 *
 * IT LIVES IN A HEADER FOR DivertedMachine.hpp's REASON. RotatingMirror's
 * vacuum guess is VacuumMirror's answer and its coils are VacuumMirror's
 * machine, so a second copy would make any disagreement between the two
 * ambiguous between the plasma and the fixture.
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
 * section 0a-pre says must never stop working ) or SUBTRACTED ( CS-3 ).
 * solveMeshed() and solveSubtracted() are the two.
 *
 * AND IT MOVES THE MODE COUNT FROM THE PLASMA TO THE COILS. A source inside
 * Gamma reaching rho_max puts mode n on Gamma at about ( rho_max/rho_Gamma )^n,
 * and the copper reaches 2.176 where the plasma stops near 1.6. So N is set by
 * the machine, not by the Mach number -- VacuumMirror's (b).
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
 * centre. That difference is the coil model, and VacuumMirror's (c) measures
 * it; geqFilaments() is geq's model in MEQ's kernel, for a case that wants the
 * coil model out of a comparison.
 *
 * Include it AFTER <boost/test/unit_test.hpp>: the builders report a broken
 * fixture through BOOST_TEST_REQUIRE, so that a missing mesh or an empty
 * subdomain names itself rather than surfacing as a solver failure. A target
 * including it defines MEQ_HALFDISC_PATH and MEQ_TEST_OUTPUT_DIR, as
 * tests/CMakeLists.txt does for VacuumMirror.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Sampler.hpp"

namespace meqtest
{
namespace mirror
{
	/// The committed geq reference, or MEQ_GEQ_VACUUM_REFERENCE -- the same
	/// override GeqSourceComparison takes, for a regenerated one.
	inline std::string referenceDirectory()
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

	inline MirrorMachine readMachine()
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

	inline MirrorMachine const &machine()
	{
		static MirrorMachine const m = readMachine();
		return m;
	}

	/// One pack, at z = sign*zc. Both carry the SAME current: a mirror, not a
	/// cusp.
	inline meq::Coil pack( double sign )
	{
		MirrorMachine const &m = machine();
		return meq::Coil( 0.5*( m.innerRadius + m.outerRadius ), sign*m.centreZ,
		                  0.5*( m.outerRadius - m.innerRadius ), 0.5*m.length,
		                  m.packCurrent );
	}

	/// The quadrature order of every reference evaluation of coilPsi() here.
	/// 16 rather than the default 32 because the error integrals evaluate it at
	/// up to 700,000 points: Coils.hpp's table puts 16 at 1e-9 relative INSIDE
	/// a coil and at round-off 5 cm clear of one, against a smallest solve
	/// error of 3e-8 relative over the machine region and 5e-7 globally.
	inline int const referenceOrder = 16;

	inline std::shared_ptr<meq::CoilSet> machineCoils()
	{
		auto set = std::make_shared<meq::CoilSet>( machine().mu0 );
		set->add( pack( -1.0 ) );
		set->add( pack( +1.0 ) );
		set->setQuadratureOrder( referenceOrder );
		return set;
	}

	/// geq's own coil model in MEQ's kernel: each pack as ten filaments at its
	/// radial centre, meq::filamentStack( pack, 1, filamentsPerPack ). A
	/// filament has no route but subtraction, so this is a ConductorField.
	inline meq::ConductorField geqFilaments()
	{
		meq::ConductorField model( machine().mu0 );
		for ( double sign : { -1.0, 1.0 } )
			for ( meq::CurrentFilament const &filament :
			      meq::filamentStack( pack( sign ), 1, machine().filamentsPerPack ) )
				model.add( filament );
		return model;
	}

	/// Gamma, and the background mesh's own arc outside it. See the file
	/// comment for why Gamma is outside the copper.
	inline double const rhoGamma = 3.0;
	inline double const meshRadius = 3.4;

	/// The mode count every case uses unless it is sweeping it --
	/// VacuumMirror (b)'s plateau with room to spare, so the truncation stays
	/// below the discretisation error at every mesh in (a).
	inline int const defaultModes = 32;

	/// geq's grid box, where MR-4's plasma lives. 0.15 m clear of the copper,
	/// so the conductor corners' r^2 log r is not in it.
	inline bool inMachineRegion( double radius, double z )
	{
		return radius <= 0.63 && std::abs( z ) <= 1.5;
	}

	/**
	 * WHAT A PLASMA CASE ASKS OF THE MESH, beyond the vacuum's.
	 *
	 * THE PACKS ARE NOT MESHED TO. A plasma case subtracts geq's filaments
	 * ( CS-3 ), so the packs are no source of the solve's and their outlines
	 * buy nothing -- and they cost one: uniform refinement of the pack's
	 * elements puts a vertex EXACTLY on a filament at level 2, where psi_c is
	 * infinite and setConductorField() refuses the pair. Measured, vertex 419
	 * at ( 0.8208, -1.35 ), filament 6.
	 *
	 * THE CUTOFF MAY BE MESHED TO. A box support's source JUMPS at
	 * |z| = axialEdge, and a jump across a line that cuts elements caps the
	 * order the way a source stepping at a plasma edge does. The line is
	 * prescribed, so it can be aligned to -- the coils' argument, one line
	 * further on: halfdisc.py's --align-plasma makes z = +-axialEdge element
	 * faces, and uniform refinement keeps them so. With @a aligned false the
	 * mesh is otherwise the same, and the line cuts elements.
	 */
	struct PlasmaMesh
	{
		double axialEdge = 0.0;  ///< > 0 asks for a plasma mesh, refined to |z| <= this
		bool aligned = false;    ///< z = +-axialEdge on element faces
	};

	/**
	 * The half-disc, made by halfdisc.py once per process and @a stem into the
	 * build tree. A stem per test, because ctest runs tests in parallel and
	 * two writing one file is a race.
	 *
	 * On the vacuum's mesh the packs are MESHED TO -- FB-2's aligned rates
	 * against cut ones are 1.99/2.88/3.01 against 1.33/1.27/1.09 -- and the
	 * machine region is refined as a plasma box would be; a plasma's mesh
	 * differs as PlasmaMesh says. Coarse vacuum out to the arc on both,
	 * because Gamma is far out and the field there is smooth.
	 *
	 * A ladder REFINES this mesh rather than re-running halfdisc.py at smaller
	 * sizes, for DivertedMachine.hpp's reason: a re-meshed geometry moves which
	 * elements the conductors occupy, a second variable in a study of the
	 * first, and uniform refinement carries the coil edges down exactly.
	 *
	 * @param plasma what a plasma case asks beyond the vacuum's mesh; see
	 *        PlasmaMesh. The default is the vacuum's.
	 */
	inline std::string const &meshFile( std::string const &stem,
	                                    PlasmaMesh const &plasma = PlasmaMesh() )
	{
		static std::map<std::string, std::string> made;
		auto const found = made.find( stem );
		if ( found != made.end() )
			return found->second;

		MirrorMachine const &m = machine();
		std::string const out = std::string( MEQ_TEST_OUTPUT_DIR ) + "/" + stem
		                        + ".msh";
		double const width = m.outerRadius - m.innerRadius;
		char geometry[ 512 ];
		if ( plasma.axialEdge > 0.0 )
			std::snprintf( geometry, sizeof geometry,
				"--plasma 0.0 %.17g 0.75 %.17g%s", -plasma.axialEdge,
				2.0*plasma.axialEdge, plasma.aligned ? " --align-plasma" : "" );
		else
			std::snprintf( geometry, sizeof geometry,
				"--coil %.17g %.17g %.17g %.17g"
				" --coil %.17g %.17g %.17g %.17g --coil-size 0.07"
				" --plasma 0.0 -1.7 0.75 3.4",
				m.innerRadius, -m.centreZ - 0.5*m.length, width, m.length,
				m.innerRadius, m.centreZ - 0.5*m.length, width, m.length );
		char buffer[ 2048 ];
		std::snprintf( buffer, sizeof buffer,
			"python3 %s --rho %.17g --size 0.4 %s --plasma-size 0.1"
			" --check -o %s > /dev/null",
			MEQ_HALFDISC_PATH, meshRadius, geometry, out.c_str() );
		BOOST_TEST_REQUIRE( std::system( buffer ) == 0,
			"halfdisc.py failed, so there is no mesh. It needs python3 with "
			"gmsh; the command was: " << buffer );
		return made.emplace( stem, out ).first->second;
	}

	/// Negative inside Gamma.
	inline double semicircle( mfem::Vector const &x )
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

	/// @param refinements uniform refinements of the gmsh mesh.
	/// @param stem        meshFile()'s, one per test -- and one per mesh, so a
	///                    test making two meshes names two.
	/// @param plasma      meshFile()'s.
	inline std::unique_ptr<Domain> buildDomain( int refinements,
	                                            std::string const &stem,
	                                            PlasmaMesh const &plasma = PlasmaMesh() )
	{
		auto d = std::make_unique<Domain>();
		d->background = std::make_unique<mfem::Mesh>(
			meshFile( stem, plasma ).c_str(), 1, 1 );
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
	///        Gamma_h -- VacuumMirror (a)'s control.
	inline std::unique_ptr<VacuumSolve> solveMeshed( Domain &d, int order, int modes,
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
	inline std::unique_ptr<VacuumSolve> solveSubtracted( Domain &d, int order,
	                                                     int modes,
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
	 * is what lets VacuumMirror (b) sweep N for the price of one table. The
	 * rule is potentialError()'s, 2k + 4. The element loop stays serial, with
	 * a caller-allocated transformation rather than the mesh's shared scratch,
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
}
}

#endif
