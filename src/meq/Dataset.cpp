#include "Dataset.hpp"

#include <algorithm>
#include <stdexcept>

namespace meq
{
	DatasetAttribute const *DatasetVariable::attribute( std::string const &wanted ) const
	{
		for ( DatasetAttribute const &one : attributes )
			if ( one.name == wanted )
				return &one;
		return nullptr;
	}

	void Dataset::addDimension( std::string const &name, std::size_t length )
	{
		if ( hasDimension( name ) )
			throw std::invalid_argument( "meq::Dataset: the dimension \"" + name
			                             + "\" is already defined" );
		dimensionList.emplace_back( name, length );
	}

	DatasetVariable &Dataset::addVariable( std::string const &name,
	                                       std::vector<std::string> const &dimensions,
	                                       DatasetValues values )
	{
		if ( variable( name ) )
			throw std::invalid_argument( "meq::Dataset: the variable \"" + name
			                             + "\" is already defined" );
		for ( std::string const &one : dimensions )
			if ( !hasDimension( one ) )
				throw std::invalid_argument( "meq::Dataset: the variable \"" + name
				                             + "\" is on \"" + one + "\", which is not a dimension" );
		variableList.push_back( DatasetVariable{ name, dimensions, std::move( values ), {} } );
		return variableList.back();
	}

	void Dataset::setAttribute( std::vector<DatasetAttribute> &on, DatasetAttribute attribute )
	{
		for ( DatasetAttribute &one : on )
			if ( one.name == attribute.name )
			{
				one = std::move( attribute );
				return;
			}
		on.push_back( std::move( attribute ) );
	}

	void Dataset::setAttribute( DatasetAttribute attribute )
	{
		setAttribute( attributeList, std::move( attribute ) );
	}

	bool Dataset::hasDimension( std::string const &name ) const
	{
		return std::any_of( dimensionList.begin(), dimensionList.end(),
		                    [ & ]( auto const &one ) { return one.first == name; } );
	}

	std::size_t Dataset::dimension( std::string const &name ) const
	{
		for ( auto const &one : dimensionList )
			if ( one.first == name )
				return one.second;
		throw std::out_of_range( "meq::Dataset: no dimension \"" + name + "\"" );
	}

	std::size_t Dataset::size( std::vector<std::string> const &dimensions ) const
	{
		std::size_t count = 1;
		for ( std::string const &one : dimensions )
			count *= dimension( one );
		return count;
	}

	DatasetVariable const *Dataset::variable( std::string const &name ) const
	{
		for ( DatasetVariable const &one : variableList )
			if ( one.name == name )
				return &one;
		return nullptr;
	}

	DatasetVariable *Dataset::variable( std::string const &name )
	{
		for ( DatasetVariable &one : variableList )
			if ( one.name == name )
				return &one;
		return nullptr;
	}

	DatasetAttribute const *Dataset::attribute( std::string const &name ) const
	{
		for ( DatasetAttribute const &one : attributeList )
			if ( one.name == name )
				return &one;
		return nullptr;
	}
}
