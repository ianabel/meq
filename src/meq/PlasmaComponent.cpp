#include "PlasmaComponent.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace meq
{
	void PlasmaComponent::setAdjacency( std::vector< int > offsets,
	                                    std::vector< int > neighbours )
	{
		// A malformed graph is checked here rather than trusted, because the
		// symptom downstream is a MASK -- a plausible list of element numbers --
		// and nothing that consumes a mask can tell a wrong one from a right
		// one. The same argument this tree makes for every closed-form check.
		if ( offsets.empty() )
			throw std::invalid_argument( "meq::PlasmaComponent::setAdjacency: the offset array must carry at least the leading zero" );
		if ( offsets.front() != 0 )
			throw std::invalid_argument( "meq::PlasmaComponent::setAdjacency: the offsets must start at zero" );

		int const nodes = static_cast< int >( offsets.size() ) - 1;
		for ( int e = 0; e < nodes; ++e )
			if ( offsets[ e + 1 ] < offsets[ e ] )
				throw std::invalid_argument( "meq::PlasmaComponent::setAdjacency: the offsets must be non-decreasing" );

		if ( offsets.back() != static_cast< int >( neighbours.size() ) )
			throw std::invalid_argument( "meq::PlasmaComponent::setAdjacency: the last offset must be the neighbour count" );

		for ( int n : neighbours )
			if ( n < 0 || n >= nodes )
				throw std::invalid_argument( "meq::PlasmaComponent::setAdjacency: a neighbour index is outside the graph" );

		rowOffsets = std::move( offsets );
		neighbourList = std::move( neighbours );
		clear();
	}

	int PlasmaComponent::nodeCount() const
	{
		return rowOffsets.empty() ? 0 : static_cast< int >( rowOffsets.size() ) - 1;
	}

	void PlasmaComponent::fill( std::vector< char > const &carriesPlasma, int seed,
	                            int blockNode )
	{
		int const nodes = nodeCount();
		if ( static_cast< int >( carriesPlasma.size() ) != nodes )
			throw std::invalid_argument( "meq::PlasmaComponent::fill: the candidate mask must carry one entry per node, and carries "
			                             + std::to_string( carriesPlasma.size() )
			                             + " against " + std::to_string( nodes ) );
		if ( seed < 0 || seed >= nodes )
			throw std::invalid_argument( "meq::PlasmaComponent::fill: the seed is outside the graph" );

		labels.assign( static_cast< std::size_t >( nodes ), -1 );
		ring.clear();
		ringSize = 0;
		labelCount = 0;
		seedLabelValue = -1;
		componentSize = 0;
		candidateSize = 0;
		filled = true;

		auto candidate = [ & ]( int e )
		{
			return carriesPlasma[ static_cast< std::size_t >( e ) ] != 0
			       && e != blockNode;
		};

		for ( int e = 0; e < nodes; ++e )
			if ( carriesPlasma[ static_cast< std::size_t >( e ) ] != 0 )
				++candidateSize;

		// EVERY component, not only the seed's, and the whole labelling costs
		// the same single pass the one fill would. componentCount() is what says
		// whether the connectivity test found anything to do; see the header.
		std::vector< int > stack;
		stack.reserve( static_cast< std::size_t >( nodes ) );

		for ( int start = 0; start < nodes; ++start )
		{
			if ( !candidate( start ) || labels[ static_cast< std::size_t >( start ) ] >= 0 )
				continue;

			int const mark = labelCount++;
			int size = 0;

			labels[ static_cast< std::size_t >( start ) ] = mark;
			stack.push_back( start );

			// Depth first over an explicit stack rather than by recursion: a
			// component here is a mesh region and can be every element in the
			// mesh, so the recursion depth would be the element count.
			while ( !stack.empty() )
			{
				int const e = stack.back();
				stack.pop_back();
				++size;

				for ( int i = rowOffsets[ static_cast< std::size_t >( e ) ];
				      i < rowOffsets[ static_cast< std::size_t >( e ) + 1 ]; ++i )
				{
					int const n = neighbourList[ static_cast< std::size_t >( i ) ];
					if ( !candidate( n ) || labels[ static_cast< std::size_t >( n ) ] >= 0 )
						continue;
					labels[ static_cast< std::size_t >( n ) ] = mark;
					stack.push_back( n );
				}
			}

			if ( mark == labels[ static_cast< std::size_t >( seed ) ] )
			{
				seedLabelValue = mark;
				componentSize = size;
			}
		}
	}

	void PlasmaComponent::fill( std::vector< char > const &interior,
	                            std::vector< char > const &carriesPlasma,
	                            int seed, int blockNode )
	{
		int const nodes = nodeCount();
		if ( static_cast< int >( carriesPlasma.size() ) != nodes )
			throw std::invalid_argument( "meq::PlasmaComponent::fill: the candidate mask must carry one entry per node" );

		// A band taken from a set that does not contain the interior would leave
		// holes -- kept elements with dropped elements inside them -- and a mask
		// with holes in it is a plausible answer nothing downstream can check.
		for ( int e = 0; e < nodes; ++e )
			if ( interior[ static_cast< std::size_t >( e ) ] != 0
			     && carriesPlasma[ static_cast< std::size_t >( e ) ] == 0 )
				throw std::invalid_argument( "meq::PlasmaComponent::fill: the interior mask is not a subset of the candidate mask, so the two rules disagree about which is stricter" );

		// The traversal is over the STRICT set, which is what separates two
		// lobes meeting at a point: the straddling band that joins them is not
		// walked at all.
		fill( interior, seed, blockNode );

		// candidateNodes() must report the inclusive rule, since that is what
		// the pointwise test would have selected and is what the fill is being
		// compared against.
		candidateSize = 0;
		for ( int e = 0; e < nodes; ++e )
			candidateSize += carriesPlasma[ static_cast< std::size_t >( e ) ] ? 1 : 0;

		if ( seedLabelValue < 0 )
			return;

		/*
		 * THE WATERSHED. One breadth-first wave from EVERY interior component at
		 * once, through the straddling band, each element going to whichever
		 * wave reaches it first.
		 *
		 * That is what makes it parameter-free and what makes it right in both
		 * directions at once. With one interior component the whole band is
		 * assigned to it, so nothing of the plasma edge is dropped; with two,
		 * the band divides where the two waves meet, which near a saddle is the
		 * pinch -- the place the two lobes actually touch.
		 *
		 * The queue is seeded in element order and drained first-in-first-out,
		 * so an element equidistant from two components goes to the lower one
		 * and the answer does not depend on anything but the graph.
		 */
		ring.assign( static_cast< std::size_t >( nodes ), 0 );
		ringSize = 0;

		std::vector< int > assignment( static_cast< std::size_t >( nodes ), -1 );
		std::vector< int > queue;
		queue.reserve( static_cast< std::size_t >( nodes ) );

		for ( int e = 0; e < nodes; ++e )
			if ( labels[ static_cast< std::size_t >( e ) ] >= 0 )
			{
				assignment[ static_cast< std::size_t >( e ) ] =
					labels[ static_cast< std::size_t >( e ) ];
				queue.push_back( e );
			}

		for ( std::size_t head = 0; head < queue.size(); ++head )
		{
			int const e = queue[ head ];
			int const mark = assignment[ static_cast< std::size_t >( e ) ];

			for ( int i = rowOffsets[ static_cast< std::size_t >( e ) ];
			      i < rowOffsets[ static_cast< std::size_t >( e ) + 1 ]; ++i )
			{
				int const n = neighbourList[ static_cast< std::size_t >( i ) ];
				std::size_t const index = static_cast< std::size_t >( n );

				// Only the band: an interior element already carries its own
				// component's label and must never be re-assigned, which is what
				// keeps the two lobes apart however the wave travels.
				if ( carriesPlasma[ index ] == 0 || interior[ index ] != 0
				     || assignment[ index ] >= 0 || n == blockNode )
					continue;

				assignment[ index ] = mark;
				queue.push_back( n );
			}
		}

		for ( int e = 0; e < nodes; ++e )
		{
			std::size_t const index = static_cast< std::size_t >( e );
			if ( interior[ index ] == 0
			     && assignment[ index ] == seedLabelValue )
			{
				ring[ index ] = 1;
				++ringSize;
			}
		}
	}

	void PlasmaComponent::setExcluded( std::vector< char > excluded )
	{
		if ( !excluded.empty() && !labels.empty()
		     && excluded.size() != labels.size() )
			throw std::invalid_argument( "meq::PlasmaComponent::setExcluded: one entry per node, or empty" );
		if ( !excluded.empty() && !rowOffsets.empty()
		     && excluded.size() + 1 != rowOffsets.size() )
			throw std::invalid_argument( "meq::PlasmaComponent::setExcluded: one entry per node, or empty" );

		excludedNode = std::move( excluded );
		excludedSize = 0;
		for ( char const e : excludedNode )
			if ( e != 0 )
				++excludedSize;
	}

	bool PlasmaComponent::isExcluded( int element ) const
	{
		if ( excludedNode.empty() || element < 0
		     || static_cast< std::size_t >( element ) >= excludedNode.size() )
			return false;
		return excludedNode[ static_cast< std::size_t >( element ) ] != 0;
	}

	int PlasmaComponent::excludedNodes() const
	{
		return excludedSize;
	}

	bool PlasmaComponent::holds( int element ) const
	{
		// AN EXCLUSION IS A STATEMENT ABOUT THE DEVICE AND OUTRANKS EVERYTHING
		// BELOW IT, including the constant-true shortcut: a caller that never
		// asked for connectivity but did name a region outside the vessel gets
		// the exclusion and nothing else.
		if ( isExcluded( element ) )
			return false;

		// The constant true before fill(), so a caller tests unconditionally and
		// a solver that was never given a connectivity test is bit-unchanged.
		if ( !filled )
			return true;
		if ( element < 0 || element >= nodeCount() )
			return false;
		if ( !ring.empty() && ring[ static_cast< std::size_t >( element ) ] != 0 )
			return true;
		return labels[ static_cast< std::size_t >( element ) ] == seedLabelValue
		       && seedLabelValue >= 0;
	}

	bool PlasmaComponent::active() const
	{
		return filled;
	}

	int PlasmaComponent::componentNodes() const
	{
		return componentSize + ringSize;
	}

	int PlasmaComponent::ringNodes() const
	{
		return ringSize;
	}

	int PlasmaComponent::candidateNodes() const
	{
		return candidateSize;
	}

	int PlasmaComponent::componentCount() const
	{
		return labelCount;
	}

	int PlasmaComponent::label( int element ) const
	{
		if ( !filled || element < 0 || element >= nodeCount() )
			return -1;
		return labels[ static_cast< std::size_t >( element ) ];
	}

	int PlasmaComponent::seedLabel() const
	{
		return seedLabelValue;
	}

	void PlasmaComponent::clear()
	{
		// THE EXCLUSION SURVIVES A clear(), which is what makes it a property of
		// the device rather than of an iterate: clear() is called between solves
		// and whenever the adjacency is rebuilt, and a wall does not move.
		labels.clear();
		ring.clear();
		labelCount = 0;
		seedLabelValue = -1;
		componentSize = 0;
		ringSize = 0;
		candidateSize = 0;
		filled = false;
	}
}
