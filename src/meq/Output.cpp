#include "Output.hpp"

#include <algorithm>
#include <fstream>
#include <map>
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

	void writeMfem( std::string const &stem, mfem::Mesh &mesh,
	                mfem::GridFunction const &potential,
	                mfem::GridFunction const &flux )
	{
		// Precision 16, not MFEM's default 8. These files are read back for an
		// exact restart, and eight digits is not exact -- it would put a 1e-8
		// perturbation into a warm start whose whole point is to be the previous
		// answer.
		auto open = [ &stem ]( char const *suffix )
		{
			std::ofstream out( stem + suffix );
			if ( !out )
				throw std::runtime_error( "meq::writeMfem: cannot write " + stem + suffix );
			out.precision( 16 );
			return out;
		};

		{ std::ofstream out = open( ".mesh" );          mesh.Print( out ); }
		{ std::ofstream out = open( "_psi.gf" );        potential.Save( out ); }
		{ std::ofstream out = open( "_grad_psi.gf" );   flux.Save( out ); }
	}

	void writeVtu( std::string const &stem, mfem::Mesh &mesh,
	               mfem::GridFunction const &potential,
	               mfem::GridFunction const &field,
	               int levelsOfDetail )
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

		// A single steady state, not a time series. Cycle 0 and time 0 are what
		// ParaView shows for a file with no time axis; an adaptive run could
		// write one cycle per refinement instead, which is not done.
		collection.SetCycle( 0 );
		collection.SetTime( 0.0 );
		collection.Save();
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

				double const r = original( rDof ), z = original( zDof );
				double projectedR = r, projectedZ = z;
				project( r, z, projectedR, projectedZ );

				displacement[ rDof ] = projectedR - r;
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
		for ( int attempt = 0; attempt < 8; ++attempt )
		{
			*nodes = original;
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
			*nodes = original;
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

	void boundaryPolyline( mfem::Mesh &mesh,
	                       std::vector<double> &r, std::vector<double> &z,
	                       int &unreached )
	{
		r.clear();
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
			double const *point = mesh.GetVertex( current );
			r.push_back( point[ 0 ] );
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
		netCDF::NcVar rVar = state->file.addVar( "R", netCDF::ncDouble, state->rDim );
		rVar.putAtt( "long_name", "Major radius" );
		rVar.putAtt( "units", "m" );
		rVar.putVar( coordinate.data() );

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

	void NetCDFWriter::boundary( std::vector<double> const &r,
	                             std::vector<double> const &z )
	{
		if ( r.size() != z.size() )
			throw std::runtime_error( "meq::NetCDFWriter::boundary: the two coordinate arrays differ in length" );
		if ( r.empty() )
			return;

		netCDF::NcDim dim = state->file.addDim( "boundary", r.size() );
		netCDF::NcVar rVar = state->file.addVar( "boundary_R", netCDF::ncDouble, dim );
		rVar.putAtt( "long_name", "Major radius of the prescribed boundary" );
		rVar.putAtt( "units", "m" );
		rVar.putVar( r.data() );

		netCDF::NcVar zVar = state->file.addVar( "boundary_Z", netCDF::ncDouble, dim );
		zVar.putAtt( "long_name", "Height of the prescribed boundary" );
		zVar.putAtt( "units", "m" );
		zVar.putVar( z.data() );
	}

#else   // MEQ_USE_NETCDF

	struct NetCDFWriter::State { };

	namespace
	{
		[[noreturn]] void unavailable()
		{
			throw std::runtime_error( "meq::NetCDFWriter: meq was built without netcdf-cxx4" );
		}
	}

	NetCDFWriter::NetCDFWriter( std::string const &, GridSampler const & ) { unavailable(); }
	NetCDFWriter::~NetCDFWriter() = default;
	void NetCDFWriter::close() { }
	void NetCDFWriter::attribute( std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, double ) { unavailable(); }
	void NetCDFWriter::attribute( std::string const &, int ) { unavailable(); }
	void NetCDFWriter::field( std::string const &, std::vector<double> const &,
	                          std::string const &, std::string const & ) { unavailable(); }
	void NetCDFWriter::boundary( std::vector<double> const &,
	                             std::vector<double> const & ) { unavailable(); }

#endif  // MEQ_USE_NETCDF
}
