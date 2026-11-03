// meq::Dataset: a NetCDF file's contents held in memory.
//
// What meq::NetCDFWriter writes, recorded as it writes it, so that a front end
// holding the run -- Python, PYTHON-PLAN.md PY-4 -- reads the answer the .nc
// carries without reading the .nc: the same dimensions, variables, attributes
// and values, from the same calls. Not a second sampler and not a reader: the
// writer fills it, and nothing else does.

#ifndef MEQ_DATASET_HPP
#define MEQ_DATASET_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace meq
{
	/// The values of a variable or an attribute, in the NetCDF type they are
	/// written as: double, int, byte (signed char) or string.
	using DatasetValues = std::variant<std::vector<double>, std::vector<int>,
	                                   std::vector<signed char>, std::vector<std::string>>;

	/// An attribute: a text, or a numeric array -- a scalar attribute is an
	/// array of one, as NetCDF holds it.
	struct DatasetAttribute
	{
		std::string name;
		std::variant<std::string, std::vector<double>, std::vector<int>,
		             std::vector<signed char>> value;
	};

	/// A variable on named dimensions, row-major, the last fastest; no
	/// dimensions is a scalar, one value.
	struct DatasetVariable
	{
		std::string name;
		std::vector<std::string> dimensions;
		DatasetValues values;
		std::vector<DatasetAttribute> attributes;

		/// The attribute of that name, or null.
		DatasetAttribute const *attribute( std::string const &name ) const;
	};

	/**
	 * A FILE'S CONTENTS, IN THE ORDER WRITTEN: dimensions, variables and global
	 * attributes. Each set* replaces an existing entry of the same name, as
	 * NetCDF's putAtt does; a dimension or a variable is defined once.
	 */
	class Dataset
	{
		public:
			/// @throws std::invalid_argument for a name already defined.
			void addDimension( std::string const &name, std::size_t length );
			/// @throws std::invalid_argument for a name already defined, or a
			///         dimension not.
			DatasetVariable &addVariable( std::string const &name,
			                              std::vector<std::string> const &dimensions,
			                              DatasetValues values );
			void setAttribute( DatasetAttribute attribute );
			static void setAttribute( std::vector<DatasetAttribute> &on,
			                          DatasetAttribute attribute );

			/// The length of that dimension. @throws std::out_of_range.
			std::size_t dimension( std::string const &name ) const;
			bool hasDimension( std::string const &name ) const;
			/// The values a variable of these dimensions holds.
			std::size_t size( std::vector<std::string> const &dimensions ) const;

			/// The variable of that name, or null; and the global attribute.
			DatasetVariable const *variable( std::string const &name ) const;
			DatasetVariable *variable( std::string const &name );
			DatasetAttribute const *attribute( std::string const &name ) const;

			std::vector<std::pair<std::string, std::size_t>> const &dimensions() const
			{
				return dimensionList;
			}
			std::vector<DatasetVariable> const &variables() const { return variableList; }
			std::vector<DatasetAttribute> const &attributes() const { return attributeList; }

		private:
			std::vector<std::pair<std::string, std::size_t>> dimensionList;
			std::vector<DatasetVariable> variableList;
			std::vector<DatasetAttribute> attributeList;
	};
}

#endif // MEQ_DATASET_HPP
