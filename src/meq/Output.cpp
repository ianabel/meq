#include "Output.hpp"

#include "Crossings.hpp"
#include "Materials.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <set>
#include <limits>
#include <stdexcept>

#ifdef MEQ_USE_NETCDF
#include <netcdf>
#endif

namespace meq
{
	bool hasNetCDF()
	{
#ifdef MEQ_USE_NETCDF
		return true;
#else
		return false;
#endif
	}

	std::string meshText( mfem::Mesh &mesh )
	{
		// Precision 17 for writeMfem()'s reason, below: the vertices are read
		// back for an exact restart.
		std::ostringstream out;
		out.precision( 17 );
		mesh.Print( out );
		return out.str();
	}

	void writeMfem( std::string const &stem, mfem::Mesh &mesh,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux )
	{
		writeMfem( stem, meshText( mesh ), potential, flux );
	}

	void writeMfem( std::string const &stem, std::string const &text,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux )
	{
		// PRECISION 17, NOT MFEM'S DEFAULT 8 AND NOT 16. These files are read
		// back for an exact restart, and eight digits would put a 1e-8
		// perturbation into a warm start whose whole point is to be the
		// previous answer.
		//
		// AND 16 IS NOT EXACT EITHER, WHICH IS EASY TO BELIEVE IT IS.
		// std::numeric_limits<double>::max_digits10 is 17, and 16 is the number
		// of digits a double is accurate TO rather than the number needed to
		// RECOVER it -- the two differ, and the difference is the last bit.
		// Measured on 200,000 doubles in [ -1, 1 ]: at precision 16, 50,204 of
		// them -- 25.1% -- do not parse back to the value written; at 17, none.
		// So an "exact restart" at 16 was exact for three coefficients in four
		// and one bit out for the rest, silently, which is the shape of defect
		// this tree keeps meeting rather than a rounding nicety.
		auto open = [ &stem ]( char const *suffix )
		{
			std::ofstream out( stem + suffix );
			if ( !out )
				throw std::runtime_error( "meq::writeMfem: cannot write " + stem + suffix );
			out.precision( 17 );
			return out;
		};

		{ std::ofstream out = open( ".mesh" );          out << text; }
		{ std::ofstream out = open( "_psi.gf" );        potential.Save( out ); }
		{ std::ofstream out = open( "_grad_psi.gf" );   flux.Save( out ); }
	}

	void writePostProcessed( std::string const &stem,
	                         mfem::GridFunction const &postProcessed )
	{
		// Precision 17, for writeMfem()'s reason and with its measurement: these
		// are coefficients of a discrete field, eight digits is not what was
		// computed, and 16 does not recover a double -- 25.1% of them come back
		// one bit out.
		std::ofstream out( stem + "_psistar.gf" );
		if ( !out )
			throw std::runtime_error( "meq::writePostProcessed: cannot write "
			                          + stem + "_psistar.gf" );
		out.precision( 17 );
		postProcessed.Save( out );
	}

	ElementCellData::ElementCellData( mfem::Mesh &mesh, MaterialTable const &table,
	                                  ElementMaterials const &record )
		: collection( 0, mesh.Dimension() ),
		  space( &mesh, &collection ),
		  region( &space ), plasmaFraction( &space ),
		  conductorMaterial( &space ), conductorFraction( &space )
	{
		if ( record.elementCount() != mesh.GetNE() )
			throw std::invalid_argument(
				"meq::ElementCellData: the record holds "
				+ std::to_string( record.elementCount() ) + " elements and the mesh "
				+ std::to_string( mesh.GetNE() )
				+ "; it was classified on another mesh" );
		if ( record.materialArea.size() != static_cast<std::size_t>( table.size() ) )
			throw std::invalid_argument(
				"meq::ElementCellData: the record holds "
				+ std::to_string( record.materialArea.size() )
				+ " materials and the table " + std::to_string( table.size() )
				+ "; it was made with another table" );

		// One dof per element, asked of the space rather than assumed to be
		// the element's index.
		mfem::Array<int> dofs;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			space.GetElementDofs( e, dofs );
			int const at = dofs[ 0 ] >= 0 ? dofs[ 0 ] : -1 - dofs[ 0 ];
			std::size_t const element = static_cast<std::size_t>( e );
			MaterialFraction const conductor = record.dominantConductor( e, table );
			region( at ) = static_cast<double>( static_cast<int>( record.region[ element ] ) );
			plasmaFraction( at ) = record.plasmaFraction[ element ];
			conductorMaterial( at ) = static_cast<double>( conductor.material );
			conductorFraction( at ) = conductor.fraction;
		}
	}

	ElementCellData::ElementCellData( mfem::Mesh &mesh,
	                                  std::vector<double> const &regionColumn,
	                                  std::vector<double> const &plasmaFractionColumn,
	                                  std::vector<double> const &conductorMaterialColumn,
	                                  std::vector<double> const &conductorFractionColumn )
		: collection( 0, mesh.Dimension() ),
		  space( &mesh, &collection ),
		  region( &space ), plasmaFraction( &space ),
		  conductorMaterial( &space ), conductorFraction( &space )
	{
		std::size_t const elements = static_cast<std::size_t>( mesh.GetNE() );
		for ( std::vector<double> const *column : { &regionColumn, &plasmaFractionColumn,
		                                            &conductorMaterialColumn,
		                                            &conductorFractionColumn } )
			if ( column->size() != elements )
				throw std::invalid_argument(
					"meq::ElementCellData: a column holds "
					+ std::to_string( column->size() ) + " elements and the mesh "
					+ std::to_string( elements ) + "; it is of another mesh" );

		mfem::Array<int> dofs;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			space.GetElementDofs( e, dofs );
			int const at = dofs[ 0 ] >= 0 ? dofs[ 0 ] : -1 - dofs[ 0 ];
			std::size_t const element = static_cast<std::size_t>( e );
			region( at ) = regionColumn[ element ];
			plasmaFraction( at ) = plasmaFractionColumn[ element ];
			conductorMaterial( at ) = conductorMaterialColumn[ element ];
			conductorFraction( at ) = conductorFractionColumn[ element ];
		}
	}

	std::vector< std::pair< std::string, mfem::GridFunction * > >
	ElementCellData::fields() const
	{
		return { { "region", &region },
		         { "plasma_fraction", &plasmaFraction },
		         { "conductor_material", &conductorMaterial },
		         { "conductor_fraction", &conductorFraction } };
	}

	void writeVtu( std::string const &stem, mfem::Mesh &mesh,
	               mfem::GridFunction const &potential,
	               mfem::GridFunction const &field,
	               int levelsOfDetail,
	               ElementCellData const *cells )
	{
		// ParaViewDataCollection takes a name and a prefix path and joins them
		// itself, so the stem has to be split rather than handed over whole --
		// given "out/run" it would otherwise write "./out/run/..." relative to
		// the working directory under a collection literally named "out/run".
		std::string directory = ".";
		std::string name = stem;
		std::size_t const slash = stem.find_last_of( '/' );
		if ( slash != std::string::npos )
		{
			directory = stem.substr( 0, slash );
			name = stem.substr( slash + 1 );
		}

		mfem::ParaViewDataCollection collection( name, &mesh );
		collection.SetPrefixPath( directory );

		// See the header: without this a P_k field is drawn as P_1.
		collection.SetHighOrderOutput( true );
		collection.SetLevelsOfDetail( std::max( 1, levelsOfDetail ) );

		// Base64 rather than ASCII. A k = 3 solution on a refined mesh is tens
		// of megabytes written out as text, and nothing reads these by eye --
		// the .vtu header stays plain XML either way, so the file is still
		// greppable for its field names.
		collection.SetDataFormat( mfem::VTKFormat::BINARY );

		// RegisterField does not copy, and it does not take a const pointer.
		// Nothing here mutates them; the const_cast is the collection's
		// interface, not a licence.
		collection.RegisterField( "psi",
			const_cast<mfem::GridFunction *>( &potential ) );
		collection.RegisterField( "B",
			const_cast<mfem::GridFunction *>( &field ) );
		if ( cells != nullptr )
			for ( auto const &named : cells->fields() )
				collection.RegisterField( named.first, named.second );

		// A single steady state, not a time series. Cycle 0 and time 0 are what
		// ParaView shows for a file with no time axis; an adaptive run could
		// write one cycle per refinement instead, which is not done.
		collection.SetCycle( 0 );
		collection.SetTime( 0.0 );
		collection.Save();
	}

	struct VtuSeries::State
	{
		std::string name;
		std::string directory;
		int detail;
		int written;
		std::unique_ptr<mfem::ParaViewDataCollection> collection;
	};

	VtuSeries::VtuSeries( std::string const &stem, int levelsOfDetail )
		: state( new State{ stem, ".", std::max( 1, levelsOfDetail ), 0, nullptr } )
	{
		std::size_t const slash = stem.find_last_of( '/' );
		if ( slash != std::string::npos )
		{
			state->directory = stem.substr( 0, slash );
			state->name = stem.substr( slash + 1 );
		}
		state->name += "_cycles";
	}

	VtuSeries::~VtuSeries() = default;

	int VtuSeries::frames() const { return state->written; }

	void VtuSeries::append( mfem::Mesh &mesh,
	                        mfem::GridFunction const &potential,
	                        mfem::GridFunction const &field,
	                        int cycle, double time,
	                        ElementCellData const *cells )
	{
		// ONE COLLECTION, KEPT ALIVE ACROSS THE WHOLE SERIES. Rebuilding it per
		// frame was tried first and is wrong in a way that looks right: the
		// Cycle directories all appear, with the correct meshes, and the .pvd
		// index then lists only the LAST of them -- so ParaView opens the file
		// and shows one frame. ParaViewDataCollection appends to its .pvd as
		// cycles are saved and does not scan the directory, so the collection
		// has to survive between frames.
		//
		// SetMesh() is what makes that safe. The adaptive loop hands over a
		// different mesh each cycle -- refined on the fitted path, a rebuilt
		// SubMesh on the curved one -- and rebinding at the top of every append
		// means the collection never holds a mesh the loop has destroyed.
		if ( !state->collection )
		{
			state->collection = std::make_unique<mfem::ParaViewDataCollection>(
				state->name, &mesh );
			state->collection->SetPrefixPath( state->directory );
			state->collection->SetHighOrderOutput( true );
			state->collection->SetLevelsOfDetail( state->detail );
			state->collection->SetDataFormat( mfem::VTKFormat::BINARY );
		}
		else
		{
			state->collection->SetMesh( &mesh );
		}

		state->collection->RegisterField( "psi",
			const_cast<mfem::GridFunction *>( &potential ) );
		state->collection->RegisterField( "B",
			const_cast<mfem::GridFunction *>( &field ) );
		if ( cells != nullptr )
			for ( auto const &named : cells->fields() )
				state->collection->RegisterField( named.first, named.second );

		state->collection->SetCycle( cycle );
		state->collection->SetTime( time );
		state->collection->Save();

		// The FIELDS are dropped immediately, so nothing here outlives the
		// GridFunctions it was handed and the caller may destroy the solver on
		// the next line. The collection itself stays, holding only the mesh --
		// which the next append rebinds before touching.
		state->collection->DeregisterField( "psi" );
		state->collection->DeregisterField( "B" );
		if ( cells != nullptr )
			for ( auto const &named : cells->fields() )
				state->collection->DeregisterField( named.first );
		++state->written;
	}

	namespace
	{
		/// The elements that have been turned inside out -- a non-positive
		/// Jacobian determinant somewhere. Checked at the nodes of a rule two
		/// orders above the geometry, which is where a curved element folds
		/// first.
		///
		/// Returns the LIST rather than a verdict, because the displacement is
		/// backed off only where it has to be. A global verdict costs every
		/// face the worst face's limit, which on a real Miller boundary means
		/// halving all of them for the sake of one.
		std::vector<int> tangledElements( mfem::Mesh &mesh )
		{
			std::vector<int> tangled;
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				mfem::ElementTransformation *transformation =
					mesh.GetElementTransformation( e );
				mfem::IntegrationRule const &rule = mfem::IntRules.Get(
					mesh.GetElementBaseGeometry( e ),
					2*transformation->Order() + 2 );
				for ( int q = 0; q < rule.GetNPoints(); ++q )
				{
					transformation->SetIntPoint( &rule.IntPoint( q ) );
					if ( transformation->Jacobian().Det() <= 0.0 )
					{
						tangled.push_back( e );
						break;
					}
				}
			}
			return tangled;
		}
	}

	int curveBoundaryOnto( mfem::Mesh &mesh, int order,
	                       std::function<void( double, double,
	                                           double &, double & )> const &project,
	                       double &applied )
	{
		applied = 0.0;
		if ( !project )
			return 0;

		mesh.SetCurvature( std::max( 1, order ) );
		mfem::GridFunction *nodes = mesh.GetNodes();
		if ( nodes == nullptr )
			return 0;
		mfem::FiniteElementSpace const *space = nodes->FESpace();

		/*
		 * THE NODES ARE FUNNELLED BEFORE THEY ARE COPIED AND EDITED, because
		 * under an mfem::Device the solve's own kernels may have read them --
		 * and a COPY of a vector whose live copy is on the device is born
		 * device-valid (MEASUREMENTS.md M-205), so every raw read of `original`
		 * below would read a stale page, or fault under `Device( "debug" )`,
		 * which is how this was found: on `fixed-a-testtokamak` the solve
		 * finished and the output stage died here (M-238). A read-WRITE,
		 * because the loop below moves the nodes on the host.
		 */
		countedHostReadWrite( *nodes, __LINE__ );

		// The undisplaced geometry, kept so that every trial fraction starts
		// from the same place and so that failure can restore it exactly.
		mfem::GridFunction const original( *nodes );

		// One displacement per node of every boundary face. Collected rather
		// than applied, because the tangling check needs to scale them
		// together -- and because a node shared by two faces must not be
		// projected twice.
		std::map<int, double> displacement;   // vdof -> how far to move it
		int moved = 0;
		for ( int b = 0; b < mesh.GetNBE(); ++b )
		{
			mfem::Array<int> vdofs;
			space->GetBdrElementVDofs( b, vdofs );
			// vdofs are laid out component by component: the first half is R,
			// the second Z, for the byNodes ordering SetCurvature installs.
			int const perComponent = vdofs.Size()/2;
			for ( int n = 0; n < perComponent; ++n )
			{
				int const rDof = vdofs[ n ];
				int const zDof = vdofs[ perComponent + n ];
				if ( displacement.count( rDof ) )
					continue;

				double const radius = original( rDof ), z = original( zDof );
				double projectedR = radius, projectedZ = z;
				project( radius, z, projectedR, projectedZ );

				displacement[ rDof ] = projectedR - radius;
				displacement[ zDof ] = projectedZ - z;
				++moved;
			}
		}
		if ( displacement.empty() )
			return 0;

		// SPREAD THE DISPLACEMENT INTO THE INTERIOR, without which this mostly
		// does not work. Moving a boundary face by O( h ) while the vertex
		// opposite it stays put compresses that element by O( h ) -- which is
		// its whole size -- so on a real case the tangling check fires and
		// backs the displacement off to a half. The boundary is then at neither
		// Gamma_h nor Gamma, which is worse than either.
		//
		// A few Jacobi sweeps of the nodal graph turn the boundary
		// displacement into a smooth field that decays inward, so the
		// compression is shared over several elements instead of falling on
		// one. This is the cheap version of solving a Laplace problem for the
		// displacement; for a band one element deep it is enough, and the
		// tangling check below is still what decides.
		std::vector<double> spread( nodes->Size(), 0.0 );
		std::vector<char> fixed( nodes->Size(), 0 );
		for ( auto const &entry : displacement )
		{
			spread[ entry.first ] = entry.second;
			fixed[ entry.first ] = 1;
		}

		// Adjacency on the SCALAR dofs, applied to each component separately.
		// Building it on the vdofs instead couples R to Z -- they share one
		// index range -- and averages a radial displacement against a vertical
		// one. Measured, that made the tangling worse rather than better: 25%
		// of the displacement survived where moving the boundary alone had
		// managed 50%.
		std::vector<std::vector<int>> neighbourDofs( space->GetNDofs() );
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::Array<int> dofs;
			space->GetElementDofs( e, dofs );
			for ( int a = 0; a < dofs.Size(); ++a )
				for ( int b = 0; b < dofs.Size(); ++b )
					if ( a != b )
						neighbourDofs[ dofs[ a ] ].push_back( dofs[ b ] );
		}

		std::vector<double> next( spread );
		for ( int sweep = 0; sweep < 20; ++sweep )
		{
			for ( int component = 0; component < 2; ++component )
				for ( int d = 0; d < space->GetNDofs(); ++d )
				{
					int const vdof = space->DofToVDof( d, component );
					if ( fixed[ vdof ] || neighbourDofs[ d ].empty() )
						continue;
					double sum = 0.0;
					for ( int n : neighbourDofs[ d ] )
						sum += spread[ space->DofToVDof( n, component ) ];
					next[ vdof ] = sum/static_cast<double>( neighbourDofs[ d ].size() );
				}
			spread.swap( next );
		}

		// PER-NODE BACKOFF. Everything goes to Gamma first; whatever tangles
		// has the displacement of its own nodes halved, and the check repeats.
		// The gap can exceed an element's own size -- a background element cut
		// by Gamma may have its far corner a full diagonal inside -- so some
		// faces genuinely cannot reach, and a global limit would hold every
		// other face back with them.
		std::vector<double> scale( spread.size(), 1.0 );
		bool clean = false;

		/*
		 * EACH ATTEMPT STARTS FROM THE ORIGINAL NODES, COPIED ON THE HOST.
		 *
		 * `*nodes = original` is device-aware: with `UseDevice()` set it
		 * uploads `original`, copies on the device and leaves the device copy
		 * the live one -- so the raw `+=` below needed a funnel after it, and
		 * the funnel then brought the nodes straight back. Eight attempts,
		 * 25.2 MB down on `fixed-f-diiid`, for a copy whose every reader is
		 * this host loop. `original` was copied from host-valid nodes, so its
		 * host copy is current; `HostWrite()` makes the nodes' host copy the
		 * live one without moving anything.
		 */
		auto restoreOriginal = [ & ]()
		{
			mfem::real_t const *const from = original.HostRead();
			mfem::real_t *const to = nodes->HostWrite();
			std::copy( from, from + original.Size(), to );
		};

		for ( int attempt = 0; attempt < 8; ++attempt )
		{
			restoreOriginal();
			for ( std::size_t d = 0; d < spread.size(); ++d )
				( *nodes )( static_cast<int>( d ) ) += scale[ d ]*spread[ d ];

			std::vector<int> const tangled = tangledElements( mesh );
			if ( tangled.empty() )
			{
				clean = true;
				break;
			}
			for ( int e : tangled )
			{
				mfem::Array<int> vdofs;
				space->GetElementVDofs( e, vdofs );
				for ( int n = 0; n < vdofs.Size(); ++n )
					scale[ vdofs[ n ] ] *= 0.5;
			}
		}

		if ( !clean )
		{
			// Still folded. Leave the mesh exactly as it was found: a faceted
			// boundary is honest, and a folded element renders as a black spike
			// that looks like a solver failure.
			restoreOriginal();
			applied = 0.0;
			return 0;
		}

		// `applied` is the fraction of boundary nodes that reached Gamma
		// exactly, which is what a caller wants to know -- not the smallest
		// scale anywhere, which one awkward corner would dominate.
		int full = 0;
		for ( auto const &entry : displacement )
			if ( scale[ entry.first ] >= 1.0 )
				++full;
		applied = displacement.empty()
			? 0.0
			: static_cast<double>( full )/static_cast<double>( displacement.size() );
		return moved;
	}

	std::function<void( double, double, double &, double & )>
	radialProjection( BoundaryShape const &shape )
	{
		BoundaryShape const *curve = &shape;
		return [ curve ]( double radius, double z, double &outR, double &outZ )
		{
			double const centreR = curve->majorRadius();
			double const centreZ = curve->centreHeight();
			double const vR = radius - centreR, vZ = z - centreZ;
			double const rho = std::hypot( vR, vZ );
			outR = radius;
			outZ = z;
			if ( rho <= 0.0 )
				return;
			double const gap = curve->levelSet( radius, z );
			outR = radius - vR*gap/rho;
			outZ = z - vZ*gap/rho;
		};
	}

	void writeEquilibriumVtu( std::string const &stem, mfem::Mesh &mesh,
	                          BoundaryShape const *shape, int geometryOrder,
	                          mfem::GridFunction const &potential,
	                          mfem::GridFunction const &field, int degree,
	                          ElementCellData const *cells )
	{
		// On the curved path the boundary is bent out onto the true Gamma, so
		// the drawn domain is Omega rather than the inscribed polygon Omega_h.
		// writeVtu() already emits VTK Lagrange cells, so a curvilinear mesh
		// needs nothing further from the format.
		if ( shape )
		{
			double curved = 0.0;
			int const curvedNodes = curveBoundaryOnto( mesh, geometryOrder,
			                                           radialProjection( *shape ),
			                                           curved );

			// Not a warning below 100%: a face that cannot reach Gamma without
			// folding its element is a statement about the mesh being coarse
			// there, and the rest of the boundary still reaches. It is reported
			// because the VTK boundary is then not uniformly Gamma, which a
			// reader comparing it against the .nc boundary_R/Z would notice.
			if ( curvedNodes > 0 && curved < 1.0 )
				std::printf( "MEQ: %.0f%% of the boundary reached Gamma; the "
				             "rest was held back to keep its elements from "
				             "folding\n", 100.0*curved );
			else if ( curvedNodes == 0 )
				std::fprintf( stderr,
					"MEQ: warning: could not bend the boundary onto Gamma "
					"without tangling an element. The VTK boundary is the "
					"polygon Gamma_h.\n" );
		}
		writeVtu( stem, mesh, potential, field, degree, cells );
	}

	void boundaryPolyline( mfem::Mesh &mesh,
	                       std::vector<double> &radius, std::vector<double> &z,
	                       int &unreached )
	{
		radius.clear();
		z.clear();
		unreached = 0;

		// Vertex -> the boundary vertices it shares a boundary element with.
		// A closed polygonal loop gives every vertex exactly two, which is what
		// makes the walk below unambiguous.
		std::map<int, std::vector<int>> neighbours;
		for ( int b = 0; b < mesh.GetNBE(); ++b )
		{
			mfem::Array<int> vertices;
			mesh.GetBdrElementVertices( b, vertices );
			// Segments in 2D. Anything else is a mesh this function was not
			// written for, and saying so beats drawing nonsense.
			if ( vertices.Size() != 2 )
				throw std::runtime_error(
					"meq::boundaryPolyline: boundary element is not a segment; "
					"this expects a 2D mesh" );
			neighbours[ vertices[ 0 ] ].push_back( vertices[ 1 ] );
			neighbours[ vertices[ 1 ] ].push_back( vertices[ 0 ] );
		}
		if ( neighbours.empty() )
			return;

		// Walk from the first boundary vertex, always leaving by the neighbour
		// that is not where we came from.
		int const start = neighbours.begin()->first;
		int current = start;
		int previous = -1;
		std::set<int> visited;

		while ( true )
		{
			visited.insert( current );
			auto const *point = mesh.GetVertex( current );
			radius.push_back( point[ 0 ] );
			z.push_back( point[ 1 ] );

			int next = -1;
			for ( int candidate : neighbours[ current ] )
				if ( candidate != previous )
				{
					next = candidate;
					break;
				}

			// A dangling end, which a closed loop does not have. Stop rather
			// than loop forever; `unreached` below reports the shortfall.
			if ( next < 0 || next == start )
				break;
			// Already seen and not the start: the boundary is not a simple
			// loop. Same treatment.
			if ( visited.count( next ) )
				break;

			previous = current;
			current = next;
		}

		unreached = static_cast<int>( neighbours.size() - visited.size() );
	}

#ifdef MEQ_USE_NETCDF

	struct NetCDFWriter::State
	{
		netCDF::NcFile file;
		netCDF::NcDim rDim;
		netCDF::NcDim zDim;
		int nR;
		int nZ;
		bool closed;
	};

	NetCDFWriter::NetCDFWriter( std::string const &path, GridSampler const &sampler )
		: state( new State{ {}, {}, {}, sampler.nodesR(), sampler.nodesZ(), false } )
	{
		try
		{
			state->file.open( path, netCDF::NcFile::replace );
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error( "meq::NetCDFWriter: cannot create " + path
			                          + ": " + error.what() );
		}

		state->rDim = state->file.addDim( "R", static_cast<std::size_t>( state->nR ) );
		state->zDim = state->file.addDim( "Z", static_cast<std::size_t>( state->nZ ) );

		// The coordinates first, so the file is self-describing even if a later
		// write fails.
		std::vector<double> coordinate;

		coordinate.resize( static_cast<std::size_t>( state->nR ) );
		for ( int i = 0; i < state->nR; ++i ) coordinate[ i ] = sampler.rAt( i );
		netCDF::NcVar radiusVar = state->file.addVar( "R", netCDF::ncDouble, state->rDim );
		radiusVar.putAtt( "long_name", "Major radius" );
		radiusVar.putAtt( "units", "m" );
		radiusVar.putVar( coordinate.data() );

		coordinate.resize( static_cast<std::size_t>( state->nZ ) );
		for ( int j = 0; j < state->nZ; ++j ) coordinate[ j ] = sampler.zAt( j );
		netCDF::NcVar zVar = state->file.addVar( "Z", netCDF::ncDouble, state->zDim );
		zVar.putAtt( "long_name", "Height" );
		zVar.putAtt( "units", "m" );
		zVar.putVar( coordinate.data() );

		// The mask, which is the thing a reader can always trust. Written as
		// signed char because NetCDF's byte is one, and 0/1 rather than a
		// boolean because the format has no boolean.
		std::vector<signed char> mask( static_cast<std::size_t>( state->nR )*state->nZ, 0 );
		for ( int j = 0; j < state->nZ; ++j )
			for ( int i = 0; i < state->nR; ++i )
				mask[ static_cast<std::size_t>( j )*state->nR + i ] =
					sampler.located( i, j ) ? 1 : 0;

		std::vector<netCDF::NcDim> const shape{ state->zDim, state->rDim };
		netCDF::NcVar maskVar = state->file.addVar( "inside", netCDF::ncByte, shape );
		maskVar.putAtt( "long_name",
		                "1 where the node lies inside the computational domain" );
		maskVar.putVar( mask.data() );

		// AND A SECOND MASK, BECAUSE `inside` CANNOT ANSWER THE QUESTION A
		// READER ACTUALLY HAS. A node in the band between Gamma_h and Gamma is
		// located, carries real data and is genuinely inside the plasma, so
		// `inside` says 1 there and should -- but that node was CONTINUED from
		// the mesh boundary rather than solved on, and its error is an order
		// worse than its neighbours'. Until now the file said only how many such
		// nodes there were, as the `extrapolated_nodes` attribute, which is a
		// count and tells nobody WHICH. Anyone computing an error norm, fitting
		// a flux surface or differencing two runs wants to be able to drop them.
		for ( int j = 0; j < state->nZ; ++j )
			for ( int i = 0; i < state->nR; ++i )
				mask[ static_cast<std::size_t>( j )*state->nR + i ] =
					sampler.wasExtended( i, j ) ? 1 : 0;

		netCDF::NcVar bandVar =
			state->file.addVar( "extrapolated", netCDF::ncByte, shape );
		bandVar.putAtt( "long_name",
		                "1 where the node was continued across the Gamma_h-to-"
		                "Gamma band rather than solved on" );
		bandVar.putVar( mask.data() );
	}

	NetCDFWriter::~NetCDFWriter()
	{
		// Swallowed on purpose: a throw from a destructor terminates, and the
		// caller who wanted the error had close() available.
		try { close(); } catch ( ... ) { }
	}

	void NetCDFWriter::close()
	{
		if ( state && !state->closed )
		{
			state->file.close();
			state->closed = true;
		}
	}

	void NetCDFWriter::attribute( std::string const &name, std::string const &value )
	{
		state->file.putAtt( name, value );
	}

	void NetCDFWriter::attribute( std::string const &name, double value )
	{
		state->file.putAtt( name, netCDF::ncDouble, value );
	}

	void NetCDFWriter::text( std::string const &name, std::string const &value,
	                         std::string const &longName )
	{
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncString );
		var.putAtt( "long_name", longName );
		char const *pointer = value.c_str();
		var.putVar( &pointer );
	}

	void NetCDFWriter::attribute( std::string const &name, int value )
	{
		state->file.putAtt( name, netCDF::ncInt, value );
	}

	void NetCDFWriter::field( std::string const &name,
	                          std::vector<double> const &values,
	                          std::string const &longName,
	                          std::string const &units )
	{
		std::size_t const expected = static_cast<std::size_t>( state->nR )*state->nZ;
		if ( values.size() != expected )
			throw std::runtime_error( "meq::NetCDFWriter::field: " + name + " has "
			                          + std::to_string( values.size() )
			                          + " values where the grid has "
			                          + std::to_string( expected ) );

		std::vector<netCDF::NcDim> const shape{ state->zDim, state->rDim };
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncDouble, shape );
		var.putAtt( "long_name", longName );
		var.putAtt( "units", units );
		// Both the fill attribute and the `inside` mask describe the same
		// absence. See the header for why both.
		var.putAtt( "_FillValue", netCDF::ncDouble,
		            std::numeric_limits<double>::quiet_NaN() );
		var.putVar( values.data() );
	}

	void NetCDFWriter::nodeMaterial( std::vector<int> const &rows )
	{
		std::size_t const expected = static_cast<std::size_t>( state->nR )*state->nZ;
		if ( rows.size() != expected )
			throw std::runtime_error( "meq::NetCDFWriter::nodeMaterial: "
			                          + std::to_string( rows.size() )
			                          + " values where the grid has "
			                          + std::to_string( expected ) );

		std::vector<netCDF::NcDim> const shape{ state->zDim, state->rDim };
		netCDF::NcVar var = state->file.addVar( "node_material", netCDF::ncInt, shape );
		var.putAtt( "long_name", "the material table's row at the node: the "
		                         "classifier's test at a point, plasma judged on "
		                         "psi as written here" );
		var.putAtt( "comment", "an index into the material dimension "
		                       "(material_name, material_kind); -1 outside the "
		                       "domain, where inside = 0" );
		var.putAtt( "_FillValue", netCDF::ncInt, -1 );
		var.putVar( rows.data() );
	}

	void NetCDFWriter::boundary( std::vector<double> const &radius,
	                             std::vector<double> const &z )
	{
		if ( radius.size() != z.size() )
			throw std::runtime_error( "meq::NetCDFWriter::boundary: the two coordinate arrays differ in length" );
		if ( radius.empty() )
			return;

		netCDF::NcDim dim = state->file.addDim( "boundary", radius.size() );
		netCDF::NcVar radiusVar = state->file.addVar( "boundary_R", netCDF::ncDouble, dim );
		radiusVar.putAtt( "long_name", "Major radius of the prescribed boundary" );
		radiusVar.putAtt( "units", "m" );
		radiusVar.putVar( radius.data() );

		netCDF::NcVar zVar = state->file.addVar( "boundary_Z", netCDF::ncDouble, dim );
		zVar.putAtt( "long_name", "Height of the prescribed boundary" );
		zVar.putAtt( "units", "m" );
		zVar.putVar( z.data() );
	}


	namespace
	{
		/// The dimension of that name at that length, creating it if it is not
		/// there and refusing it if it is there at another length.
		netCDF::NcDim sharedDimension( netCDF::NcFile &file,
		                               std::string const &name,
		                               std::size_t length )
		{
			netCDF::NcDim existing = file.getDim( name );
			if ( existing.isNull() )
				return file.addDim( name, length );
			if ( existing.getSize() != length )
				throw std::runtime_error(
					"meq::NetCDFWriter::vector: the dimension \"" + name
					+ "\" is already " + std::to_string( existing.getSize() )
					+ " long and this array is " + std::to_string( length )
					+ ". Two arrays sharing a dimension name must share its "
					  "length" );
			return existing;
		}
	}

	void NetCDFWriter::vector( std::string const &dimension,
	                           std::string const &name,
	                           std::vector<double> const &values,
	                           std::string const &longName,
	                           std::string const &units )
	{
		if ( values.empty() )
			throw std::runtime_error( "meq::NetCDFWriter::vector: " + name
			                          + " is empty; a zero-length NetCDF "
			                            "dimension is UNLIMITED, which is not "
			                            "what an empty array means. The caller "
			                            "should not write the variable at all" );

		netCDF::NcDim dim = sharedDimension( state->file, dimension,
		                                     values.size() );
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncDouble, dim );
		var.putAtt( "long_name", longName );
		if ( !units.empty() )
			var.putAtt( "units", units );
		var.putVar( values.data() );
	}

	void NetCDFWriter::vector( std::string const &dimension,
	                           std::string const &name,
	                           std::vector<int> const &values,
	                           std::string const &longName )
	{
		if ( values.empty() )
			throw std::runtime_error( "meq::NetCDFWriter::vector: " + name
			                          + " is empty; see the double overload" );

		netCDF::NcDim dim = sharedDimension( state->file, dimension,
		                                     values.size() );
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncInt, dim );
		var.putAtt( "long_name", longName );
		var.putVar( values.data() );
	}

	void NetCDFWriter::vector( std::string const &dimension,
	                           std::string const &name,
	                           std::vector<std::string> const &values,
	                           std::string const &longName )
	{
		if ( values.empty() )
			throw std::runtime_error( "meq::NetCDFWriter::vector: " + name
			                          + " is empty; see the double overload" );

		netCDF::NcDim dim = sharedDimension( state->file, dimension,
		                                     values.size() );
		netCDF::NcVar var = state->file.addVar( name, netCDF::ncString, dim );
		var.putAtt( "long_name", longName );
		std::vector<char const *> pointers;
		pointers.reserve( values.size() );
		for ( std::string const &one : values )
			pointers.push_back( one.c_str() );
		var.putVar( pointers.data() );
	}

	namespace
	{
		/// A one-dimensional variable, written whole. The _FillValue goes on
		/// BEFORE the data, which NetCDF-4 requires.
		template <typename T>
		netCDF::NcVar column( netCDF::NcFile &file, std::string const &name,
		                      netCDF::NcType const &type, netCDF::NcDim const &dim,
		                      std::vector<T> const &values,
		                      std::string const &longName,
		                      std::string const &units = std::string(),
		                      bool fillNaN = false )
		{
			netCDF::NcVar var = file.addVar( name, type, dim );
			var.putAtt( "long_name", longName );
			if ( !units.empty() )
				var.putAtt( "units", units );
			if ( fillNaN )
				var.putAtt( "_FillValue", netCDF::ncDouble,
				            std::numeric_limits<double>::quiet_NaN() );
			var.putVar( values.data() );
			return var;
		}

		/// A name made a legal NetCDF attribute suffix: letters, digits and
		/// underscores, the driver's variableName() rule.
		std::string attributeSuffix( std::string const &name )
		{
			std::string safe;
			for ( char const c : name )
				safe += ( std::isalnum( static_cast<unsigned char>( c ) ) || c == '_' )
					? c : '_';
			return safe.empty() ? std::string( "unnamed" ) : safe;
		}

		/// CF's `flag_values` / `flag_meanings` pair for an enumeration of
		/// @a count values named by @a nameOf.
		template <typename T, typename Name>
		void flags( netCDF::NcVar &var, netCDF::NcType const &type, int count,
		            Name const &nameOf )
		{
			std::vector<T> values;
			std::string meanings;
			for ( int i = 0; i < count; ++i )
			{
				values.push_back( static_cast<T>( i ) );
				meanings += ( i == 0 ? "" : " " ) + std::string( nameOf( i ) );
			}
			var.putAtt( "flag_values", type, values.size(), values.data() );
			var.putAtt( "flag_meanings", meanings );
		}
	}

	void NetCDFWriter::materials( MaterialTable const &table,
	                              ElementMaterials const &record )
	{
		int const elements = record.elementCount();
		if ( elements <= 0 )
			throw std::invalid_argument( "meq::NetCDFWriter::materials: the record has no elements" );
		if ( record.materialArea.size() != static_cast<std::size_t>( table.size() ) )
			throw std::invalid_argument( "meq::NetCDFWriter::materials: the record holds "
			                             + std::to_string( record.materialArea.size() )
			                             + " materials and the table "
			                             + std::to_string( table.size() )
			                             + "; it was made with another table" );

		netCDF::NcFile &file = state->file;
		std::size_t const n = static_cast<std::size_t>( elements );

		// ---- the ragged record, CF's contiguous ragged array ----
		netCDF::NcDim const elementDim = sharedDimension( file, "element", n );
		netCDF::NcDim const entryDim =
			sharedDimension( file, "element_entry", record.entries.size() );

		std::vector<int> counts( n );
		for ( std::size_t e = 0; e < n; ++e )
			counts[ e ] = record.offsets[ e + 1 ] - record.offsets[ e ];
		netCDF::NcVar countVar = column( file, "element_material_count", netCDF::ncInt,
		                                 elementDim, counts,
		                                 "number of materials in the element, indexed "
		                                 "as the elements of the mesh variable" );
		countVar.putAtt( "sample_dimension", "element_entry" );

		std::vector<int> entryMaterial;
		std::vector<double> entryFraction;
		for ( MaterialFraction const &entry : record.entries )
		{
			entryMaterial.push_back( entry.material );
			entryFraction.push_back( entry.fraction );
		}
		column( file, "element_material", netCDF::ncInt, entryDim, entryMaterial,
		        "row of the material table" );
		column( file, "element_fraction", netCDF::ncDouble, entryDim, entryFraction,
		        "fraction of the element's area the material occupies, at the "
		        "solve's quadrature points; an element's fractions sum to 1", "1" );

		// ---- dense, one per element ----
		std::vector<signed char> regions( n );
		for ( std::size_t e = 0; e < n; ++e )
			regions[ e ] = static_cast<signed char>( record.region[ e ] );
		netCDF::NcVar regionVar = column( file, "element_region", netCDF::ncByte,
		                                  elementDim, regions,
		                                  "region derived from the element's fractions, "
		                                  "for plotting" );
		flags<signed char>( regionVar, netCDF::ncByte, elementRegionCount,
		                    []( int i ) { return elementRegionName( static_cast<ElementRegion>( i ) ); } );

		column( file, "element_plasma_fraction", netCDF::ncDouble, elementDim,
		        record.plasmaFraction,
		        "fraction of the element's area that carries plasma current", "1" );

		// THE DOMINANT CONDUCTOR, which the VTK's cells draw. Derivable from the
		// ragged entries and material_kind, and written anyway so that
		// meq-convert draws exactly the cells the run would have drawn
		// without a second implementation of the rule that picks it.
		std::vector<int> conductorRows( n );
		std::vector<double> conductorShares( n );
		for ( std::size_t e = 0; e < n; ++e )
		{
			MaterialFraction const conductor =
				record.dominantConductor( static_cast<int>( e ), table );
			conductorRows[ e ] = conductor.material;
			conductorShares[ e ] = conductor.fraction;
		}
		column( file, "element_conductor_material", netCDF::ncInt, elementDim,
		        conductorRows,
		        "row of the material table of the element's largest conductor, "
		        "-1 for none" );
		column( file, "element_conductor_fraction", netCDF::ncDouble, elementDim,
		        conductorShares, "that conductor's fraction of the element, 0 for none",
		        "1" );

		std::vector<signed char> excluded( record.excluded.begin(), record.excluded.end() );
		column( file, "element_excluded", netCDF::ncByte, elementDim, excluded,
		        "1 where [source] ExcludeAttributes names the element's attribute" );
		column( file, "element_attribute", netCDF::ncInt, elementDim, record.attribute,
		        "the element's attribute in the mesh variable" );

		// ---- the material table, every column of section 2.1 ----
		std::vector<Material> const &rows = table.materials();
		std::vector<std::string> names;
		std::vector<int> kinds, circuits, curveIndex;
		std::vector<double> currents, resistivities, permeabilities;
		for ( Material const &row : rows )
		{
			names.push_back( row.name );
			kinds.push_back( static_cast<int>( row.kind ) );
			circuits.push_back( row.circuit );
			curveIndex.push_back( row.permeabilityCurve );
			currents.push_back( row.current );
			resistivities.push_back( row.resistivity );
			permeabilities.push_back( row.relativePermeability );
		}
		vector( "material", "material_name", names, "material name" );
		netCDF::NcDim const materialDim = file.getDim( "material" );
		netCDF::NcVar kindVar = column( file, "material_kind", netCDF::ncInt,
		                                materialDim, kinds, "kind of material" );
		flags<int>( kindVar, netCDF::ncInt, 5,
		            []( int i ) { return materialKindName( static_cast<MaterialKind>( i ) ); } );
		column( file, "material_circuit", netCDF::ncInt, materialDim, circuits,
		        "row of the circuit table, -1 for none" );
		column( file, "material_current", netCDF::ncDouble, materialDim, currents,
		        "the coil's resolved total current", "A", true );
		column( file, "material_resistivity", netCDF::ncDouble, materialDim,
		        resistivities, "resistivity of a passive structure", "ohm m", true );
		column( file, "material_relative_permeability", netCDF::ncDouble, materialDim,
		        permeabilities,
		        "constant relative permeability; the fill where a curve holds it", "1",
		        true );
		column( file, "material_permeability_curve", netCDF::ncInt, materialDim,
		        curveIndex, "row of the permeability curve table, -1 for none" );

		// ---- the circuits, when there are any ----
		if ( !table.circuits().empty() )
		{
			std::vector<std::string> circuitNames;
			std::vector<double> circuitCurrents;
			std::vector<signed char> supplied;
			for ( Circuit const &row : table.circuits() )
			{
				circuitNames.push_back( row.name );
				circuitCurrents.push_back( row.current );
				supplied.push_back( row.supplied ? 1 : 0 );
			}
			vector( "circuit", "circuit_name", circuitNames,
			        "circuit name; \"default\" groups the coils in no [[circuits]] block" );
			netCDF::NcDim const circuitDim = file.getDim( "circuit" );
			column( file, "circuit_current", netCDF::ncDouble, circuitDim,
			        circuitCurrents,
			        "the supply's current, each member carrying its multiplier "
			        "times it; the fill for the unsupplied default group", "A", true );
			column( file, "circuit_supplied", netCDF::ncByte, circuitDim, supplied,
			        "1 for a [[circuits]] block, 0 for a grouping with no supply" );
		}

		// ---- the permeability curves, a second ragged array ----
		std::vector<PermeabilityCurve> const &curves = table.permeabilityCurves();
		if ( !curves.empty() )
		{
			std::vector<std::string> owners;
			std::vector<int> pointCounts;
			std::vector<double> h, b;
			for ( PermeabilityCurve const &curve : curves )
			{
				owners.push_back( curve.material );
				pointCounts.push_back( static_cast<int>( curve.fieldStrength.size() ) );
				h.insert( h.end(), curve.fieldStrength.begin(), curve.fieldStrength.end() );
				b.insert( b.end(), curve.fluxDensity.begin(), curve.fluxDensity.end() );
			}
			vector( "permeability_curve", "permeability_curve_material", owners,
			        "material_name of the row this curve belongs to" );
			netCDF::NcDim const curveDim = file.getDim( "permeability_curve" );
			netCDF::NcVar pointVar = column( file, "permeability_curve_point_count",
			                                 netCDF::ncInt, curveDim, pointCounts,
			                                 "number of ( H, B ) points in the curve" );
			pointVar.putAtt( "sample_dimension", "permeability_curve_point" );
			netCDF::NcDim const pointDim =
				sharedDimension( file, "permeability_curve_point", h.size() );
			column( file, "permeability_curve_h", netCDF::ncDouble, pointDim, h,
			        "magnetic field strength H, increasing from 0", "A/m" );
			netCDF::NcVar bVar = column( file, "permeability_curve_b", netCDF::ncDouble,
			                             pointDim, b, "flux density B at H", "T" );
			bVar.putAtt( "interpolation", permeabilityInterpolation );
		}

		// ---- the global summary ----
		for ( int r = 0; r < elementRegionCount; ++r )
		{
			ElementRegion const which = static_cast<ElementRegion>( r );
			file.putAtt( std::string( "region_elements_" ) + elementRegionName( which ),
			             netCDF::ncInt, record.regionCount( which ) );
		}
		std::set<std::string> used;
		for ( std::size_t m = 0; m < rows.size(); ++m )
		{
			// A COIL MAY SHARE A NAME, with another coil or with "vacuum", and
			// a second putAtt would overwrite the first without a word.
			std::string const base = "material_area_" + attributeSuffix( rows[ m ].name );
			std::string key = base;
			for ( int suffix = 0; !used.insert( key ).second; ++suffix )
				key = base + "_" + std::to_string( m ) + ( suffix == 0 ? "" : "_" + std::to_string( suffix ) );
			file.putAtt( key, netCDF::ncDouble, record.materialArea[ m ] );
		}
		file.putAtt( "region_faults", netCDF::ncInt, record.faults );
	}

	struct FluxGridWriter::State
	{
		netCDF::NcFile file;
		bool closed;
	};

	namespace
	{
		/// The surfaces' common node count, or a refusal.
		///
		/// A NetCDF variable is rectangular and a family is not obliged to be,
		/// so this is checked rather than assumed. meq::extractFluxSurfaces()
		/// fits every surface at the same angle count, so a ragged family means
		/// something assembled one by hand and got it wrong -- which is worth a
		/// message rather than a truncated file.
		std::size_t commonNodeCount( FluxSurfaceFamily const &family )
		{
			if ( family.empty() )
				throw std::invalid_argument(
					"meq::FluxGridWriter: the family is empty, so there is no "
					"( Psi, theta ) grid to write" );

			std::size_t const angles = family.surfaces.front().count();
			if ( angles < 3 )
				throw std::invalid_argument(
					"meq::FluxGridWriter: the first surface carries "
					+ std::to_string( angles )
					+ " nodes, which encloses nothing" );

			for ( std::size_t i = 0; i < family.size(); ++i )
			{
				FluxSurface const &surface = family.surfaces[ i ];

				// EVERY COLUMN, NOT JUST R. count() reports R.size(), so a
				// family whose z or whose band mask is short would pass a check
				// on count() alone and then be read past its end -- which is
				// the one failure here that would not announce itself.
				if ( surface.count() == angles
				     && surface.z.size() == angles
				     && surface.extended.size() == angles )
					continue;

				throw std::invalid_argument(
					"meq::FluxGridWriter: surface " + std::to_string( i )
					+ " carries " + std::to_string( surface.count() ) + " R, "
					+ std::to_string( surface.z.size() ) + " Z and "
					+ std::to_string( surface.extended.size() )
					+ " mask entries where the first surface carries "
					+ std::to_string( angles )
					+ " of each, so the family is ragged and cannot be written "
					"as a rectangular ( flux, theta ) array" );
			}

			return angles;
		}
	}

	FluxGridWriter::FluxGridWriter( std::string const &path,
	                                FluxSurfaceFamily const &family )
		: state( new State{ {}, false } )
	{
		std::size_t const angles = commonNodeCount( family );
		std::size_t const count = family.size();

		try
		{
			state->file.open( path, netCDF::NcFile::replace );
		}
		catch ( netCDF::exceptions::NcException const &error )
		{
			throw std::runtime_error( "meq::FluxGridWriter: cannot create " + path
			                          + ": " + error.what() );
		}

		netCDF::NcDim const fluxDim = state->file.addDim( "flux", count );
		netCDF::NcDim const thetaDim = state->file.addDim( "theta", angles );

		// The coordinates first, so the file is self-describing even if a later
		// write fails.
		auto column = [ & ]( char const *name, char const *longName,
		                     char const *units,
		                     double ( *of )( FluxSurface const & ) )
		{
			std::vector<double> values( count, 0.0 );
			for ( std::size_t i = 0; i < count; ++i )
				values[ i ] = of( family.surfaces[ i ] );

			netCDF::NcVar var = state->file.addVar( name, netCDF::ncDouble,
			                                       fluxDim );
			var.putAtt( "long_name", longName );
			var.putAtt( "units", units );
			var.putVar( values.data() );
		};

		column( "rho", "Flux label, the square root of the normalised flux", "1",
		        []( FluxSurface const &s ) { return s.radial; } );
		column( "normalised_flux",
		        "Normalised poloidal flux, 0 on the magnetic axis and 1 on the "
		        "plasma boundary", "1",
		        []( FluxSurface const &s ) { return s.normalisedFlux; } );
		column( "psi", "Poloidal flux per radian at the surface", "Wb/rad",
		        []( FluxSurface const &s ) { return s.level; } );

		std::vector<double> theta( angles, 0.0 );
		double const twoPi = 6.283185307179586476925286766559;
		for ( std::size_t j = 0; j < angles; ++j )
			theta[ j ] = twoPi*static_cast<double>( j )
			             /static_cast<double>( angles );

		netCDF::NcVar thetaVar = state->file.addVar( "theta", netCDF::ncDouble,
		                                            thetaDim );
		thetaVar.putAtt( "long_name",
		                 "Geometric poloidal angle about the magnetic axis" );
		thetaVar.putAtt( "units", "rad" );
		thetaVar.putVar( theta.data() );

		// The geometry: ( flux, theta ), theta fastest.
		std::vector<netCDF::NcDim> const shape{ fluxDim, thetaDim };
		std::vector<double> plane( count*angles, 0.0 );

		auto surfaceColumn = [ & ]( char const *name, char const *longName,
		                            char const *units,
		                            std::vector<double> const
		                                FluxSurface::*member )
		{
			for ( std::size_t i = 0; i < count; ++i )
				for ( std::size_t j = 0; j < angles; ++j )
					plane[ i*angles + j ] = ( family.surfaces[ i ].*member )[ j ];

			netCDF::NcVar var = state->file.addVar( name, netCDF::ncDouble,
			                                        shape );
			var.putAtt( "long_name", longName );
			var.putAtt( "units", units );
			var.putVar( plane.data() );
		};

		surfaceColumn( "R", "Major radius of the flux surface", "m",
		               &FluxSurface::radius );
		surfaceColumn( "Z", "Height of the flux surface", "m",
		               &FluxSurface::z );

		// THE BAND MASK, PER NODE. A count is not a mask: CLAUDE.md records that
		// carrying only `extrapolated_nodes` on the ( R, Z ) file meant nothing
		// downstream could tell WHICH nodes had been continued outward from
		// Gamma_h. The same obligation applies here, one dimension up, and it is
		// sharper: a flux surface can be inside Omega_h at one theta and outside
		// it at the next, so a per-SURFACE flag under-reports in exactly the
		// place a q( psi ) profile cares about.
		std::vector<signed char> mask( count*angles, 0 );
		for ( std::size_t i = 0; i < count; ++i )
			for ( std::size_t j = 0; j < angles; ++j )
				mask[ i*angles + j ] =
					family.surfaces[ i ].extended[ j ] != 0 ? 1 : 0;

		netCDF::NcVar bandVar = state->file.addVar( "extrapolated",
		                                            netCDF::ncByte, shape );
		bandVar.putAtt( "long_name",
		                "1 where the node's field was continued across the "
		                "Gamma_h-to-Gamma band rather than solved on" );
		bandVar.putVar( mask.data() );

		column( "V_prime", "dV/dpsi, the volume derivative", "m^3 rad/Wb",
		        []( FluxSurface const &s ) { return s.vPrime; } );
		column( "volume", "Volume enclosed by the surface", "m^3",
		        []( FluxSurface const &s ) { return s.volume; } );
		column( "cross_section_area",
		        "Poloidal cross-section area enclosed by the surface", "m^2",
		        []( FluxSurface const &s ) { return s.crossSectionArea; } );
		column( "arc_length", "Poloidal circumference of the surface", "m",
		        []( FluxSurface const &s ) { return s.arcLength; } );
		column( "surface_area", "Area of the toroidal flux surface", "m^2",
		        []( FluxSurface const &s ) { return s.surfaceArea; } );
		column( "inverse_R_squared", "Flux-surface average of R^-2", "m^-2",
		        []( FluxSurface const &s ) { return s.inverseRSquared; } );
		column( "grad_psi_squared_over_R_squared",
		        "Flux-surface average of | grad psi |^2 / R^2",
		        "Wb^2 rad^-2 m^-4",
		        []( FluxSurface const &s )
		        { return s.gradPsiSquaredOverRSquared; } );
		column( "abs_grad_psi", "Flux-surface average of | grad psi |",
		        "Wb rad^-1 m^-1",
		        []( FluxSurface const &s ) { return s.absGradPsi; } );
		column( "grad_psi_squared", "Flux-surface average of | grad psi |^2",
		        "Wb^2 rad^-2 m^-2",
		        []( FluxSurface const &s ) { return s.gradPsiSquared; } );

		// PRESENT ONLY WHEN IT MEANS SOMETHING. g( psi ) = R B_toroidal is the
		// caller's to supply -- a meq::Source carries g g' and not g -- so a
		// family without one has no safety factor, and writing a column of
		// zeroes would be indistinguishable from a machine with none.
		if ( family.safetyFactorAvailable )
			column( "safety_factor",
			        "Safety factor, V' g < R^-2 > / 4 pi^2", "1",
			        []( FluxSurface const &s ) { return s.safetyFactor; } );

		column( "worst_residual",
		        "Worst | psi_h - level | over the nodes of the surface",
		        "Wb/rad",
		        []( FluxSurface const &s ) { return s.worstResidual; } );
		column( "transversality",
		        "min | u x t | over the angle fit: how close a ray came to being "
		        "tangent to the surface", "1",
		        []( FluxSurface const &s ) { return s.transversality; } );

		std::vector<signed char> band( count, 0 );
		for ( std::size_t i = 0; i < count; ++i )
			band[ i ] = family.surfaces[ i ].crossesBand ? 1 : 0;

		netCDF::NcVar surfaceBandVar = state->file.addVar( "band",
		                                                   netCDF::ncByte,
		                                                   fluxDim );
		surfaceBandVar.putAtt( "long_name",
		                       "1 where any node of the surface is band data" );
		surfaceBandVar.putVar( band.data() );

		// The family's own provenance. It goes in unconditionally rather than
		// being left to the caller, because a file that does not say which axis
		// its label is normalised against cannot be compared with another one.
		state->file.putAtt( "flux_label", "rho = sqrt( Psi_N )" );
		state->file.putAtt( "axis_r", netCDF::ncDouble, family.axisR );
		state->file.putAtt( "axis_z", netCDF::ncDouble, family.axisZ );
		state->file.putAtt( "psi_axis", netCDF::ncDouble, family.psiAxis );
		state->file.putAtt( "psi_boundary", netCDF::ncDouble,
		                    family.psiBoundary );
		state->file.putAtt( "inner_cut", netCDF::ncDouble, family.innerCut );
		state->file.putAtt( "outer_cut", netCDF::ncDouble, family.outerCut );
		state->file.putAtt( "extrapolated_nodes", netCDF::ncInt,
		                    family.extendedNodes() );
		state->file.putAtt( "worst_residual", netCDF::ncDouble,
		                    family.worstResidual() );
	}

	FluxGridWriter::~FluxGridWriter()
	{
		// Swallowed for the reason NetCDFWriter's is: a throw from a destructor
		// terminates, and a caller who wants to see the error calls close().
		try
		{
			close();
		}
		catch ( ... )
		{
		}
	}

	void FluxGridWriter::close()
	{
		if ( state->closed )
			return;
		state->file.close();
		state->closed = true;
	}

	void FluxGridWriter::attribute( std::string const &name,
	                                std::string const &value )
	{
		state->file.putAtt( name, value );
	}

	void FluxGridWriter::attribute( std::string const &name, double value )
	{
		state->file.putAtt( name, netCDF::ncDouble, value );
	}

	void FluxGridWriter::attribute( std::string const &name, int value )
	{
		state->file.putAtt( name, netCDF::ncInt, value );
	}

#else   // MEQ_USE_NETCDF

	struct NetCDFWriter::State { };

	namespace
	{
		[[noreturn]] void unavailable()
		{
			throw std::runtime_error( "meq::NetCDFWriter: MEQ was built without netcdf-cxx4" );
		}
	}

	NetCDFWriter::NetCDFWriter( std::string const &, GridSampler const & ) { unavailable(); }
	NetCDFWriter::~NetCDFWriter() = default;
	void NetCDFWriter::close() { }
	void NetCDFWriter::vector( std::string const &, std::string const &,
	                           std::vector<double> const &,
	                           std::string const &, std::string const & )
	{ unavailable(); }
	void NetCDFWriter::vector( std::string const &, std::string const &,
	                           std::vector<int> const &, std::string const & )
	{ unavailable(); }
	void NetCDFWriter::vector( std::string const &, std::string const &,
	                           std::vector<std::string> const &, std::string const & )
	{ unavailable(); }
	void NetCDFWriter::materials( MaterialTable const &, ElementMaterials const & )
	{ unavailable(); }
	void NetCDFWriter::nodeMaterial( std::vector<int> const & ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, double ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, int ) { unavailable(); }
	void NetCDFWriter::field( std::string const &, std::vector<double> const &,
	                          std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::boundary( std::vector<double> const &,
	                             std::vector<double> const & ) { unavailable(); }
	void NetCDFWriter::text( std::string const &, std::string const &,
	                         std::string const & ) { unavailable(); }


	struct FluxGridWriter::State { };

	namespace
	{
		[[noreturn]] void fluxGridUnavailable()
		{
			throw std::runtime_error( "meq::FluxGridWriter: MEQ was built without netcdf-cxx4" );
		}
	}

	FluxGridWriter::FluxGridWriter( std::string const &,
	                                FluxSurfaceFamily const & ) { fluxGridUnavailable(); }
	FluxGridWriter::~FluxGridWriter() = default;
	void FluxGridWriter::close() { }
	void FluxGridWriter::attribute( std::string const &, std::string const & ) { fluxGridUnavailable(); }
	void FluxGridWriter::attribute( std::string const &, double ) { fluxGridUnavailable(); }
	void FluxGridWriter::attribute( std::string const &, int ) { fluxGridUnavailable(); }
#endif  // MEQ_USE_NETCDF
}
