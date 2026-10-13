#ifndef MEQ_TESTS_VTU_READER_HPP
#define MEQ_TESTS_VTU_READER_HPP

/*
 * READING A `.vtu` BACK, WITHOUT VTK AND WITHOUT MFEM.
 *
 * `ELEMENT-REGIONS-PLAN.md` ER-2 asserts what the VTK carries by reading the
 * file, so the reader must share no code with the writer: this decodes the
 * XML, the base64 and the zlib by hand. It reads what
 * `mfem::ParaViewDataCollection` writes -- inline `format="binary"`
 * DataArrays, each two base64 strings (a UInt32 header, then the bytes),
 * zlib-compressed when the file names `vtkZLibDataCompressor` -- and refuses
 * anything else rather than misreading it.
 *
 * Header-only and test-only: DriverAcceptance and OutputTests both include it.
 */

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#if __has_include( <zlib.h> )
#include <zlib.h>
#define MEQ_TESTS_VTU_ZLIB 1
#endif

namespace meq
{
	namespace tests
	{
		namespace vtu
		{
			inline std::vector< unsigned char > decodeBase64( std::string const &text )
			{
				auto value = []( char c ) -> int
				{
					if ( c >= 'A' && c <= 'Z' ) return c - 'A';
					if ( c >= 'a' && c <= 'z' ) return c - 'a' + 26;
					if ( c >= '0' && c <= '9' ) return c - '0' + 52;
					if ( c == '+' ) return 62;
					if ( c == '/' ) return 63;
					return -1;
				};
				std::vector< unsigned char > out;
				unsigned int buffer = 0;
				int bits = 0;
				for ( char const c : text )
				{
					if ( c == '=' )
						break;
					int const v = value( c );
					if ( v < 0 )
						throw std::runtime_error( "vtu: not base64" );
					buffer = ( buffer << 6 ) | static_cast< unsigned int >( v );
					bits += 6;
					if ( bits >= 8 )
					{
						bits -= 8;
						out.push_back( static_cast< unsigned char >( ( buffer >> bits ) & 0xFFu ) );
					}
				}
				return out;
			}

			/// The length in characters of the padded base64 of @a bytes.
			inline std::size_t encodedLength( std::size_t bytes )
			{
				return 4*( ( bytes + 2 )/3 );
			}

			inline std::string attributeOf( std::string const &tag, std::string const &name )
			{
				std::string const key = " " + name + "=\"";
				std::size_t const at = tag.find( key );
				if ( at == std::string::npos )
					return "";
				std::size_t const start = at + key.size();
				return tag.substr( start, tag.find( '"', start ) - start );
			}

			inline std::uint32_t word( std::vector< unsigned char > const &bytes, std::size_t i )
			{
				std::uint32_t value = 0;
				std::memcpy( &value, bytes.data() + 4*i, 4 );
				return value;
			}

			/// The raw bytes of one DataArray's payload.
			inline std::vector< unsigned char > payload( std::string text, bool compressed )
			{
				std::string compact;
				for ( char const c : text )
					if ( c != ' ' && c != '\n' && c != '\r' && c != '\t' )
						compact.push_back( c );

				if ( !compressed )
				{
					std::size_t const head = encodedLength( 4 );
					std::vector< unsigned char > const size =
						decodeBase64( compact.substr( 0, head ) );
					std::vector< unsigned char > data = decodeBase64( compact.substr( head ) );
					if ( size.size() < 4 || data.size() < word( size, 0 ) )
						throw std::runtime_error( "vtu: short uncompressed array" );
					data.resize( word( size, 0 ) );
					return data;
				}

#ifdef MEQ_TESTS_VTU_ZLIB
				// The block count first, from the first word, which fixes how
				// long the header string is.
				std::vector< unsigned char > const first =
					decodeBase64( compact.substr( 0, encodedLength( 4 ) ) );
				std::size_t const blocks = word( first, 0 );
				std::size_t const headerBytes = 4*( 3 + blocks );
				std::size_t const head = encodedLength( headerBytes );
				std::vector< unsigned char > const header =
					decodeBase64( compact.substr( 0, head ) );
				std::vector< unsigned char > const packed = decodeBase64( compact.substr( head ) );

				std::size_t const blockSize = word( header, 1 );
				std::size_t const lastSize = word( header, 2 );
				std::vector< unsigned char > out;
				std::size_t offset = 0;
				for ( std::size_t b = 0; b < blocks; ++b )
				{
					std::size_t const packedSize = word( header, 3 + b );
					uLongf size = static_cast< uLongf >(
						( b + 1 == blocks && lastSize != 0 ) ? lastSize : blockSize );
					std::vector< unsigned char > block( size );
					if ( offset + packedSize > packed.size()
					     || uncompress( block.data(), &size, packed.data() + offset,
					                    static_cast< uLong >( packedSize ) ) != Z_OK )
						throw std::runtime_error( "vtu: a zlib block did not inflate" );
					block.resize( size );
					out.insert( out.end(), block.begin(), block.end() );
					offset += packedSize;
				}
				return out;
#else
				throw std::runtime_error( "vtu: the file is compressed and this "
				                          "test was built without zlib.h" );
#endif
			}

			template< typename T >
			void append( std::vector< unsigned char > const &bytes, std::vector< double > &out )
			{
				std::size_t const n = bytes.size()/sizeof( T );
				for ( std::size_t i = 0; i < n; ++i )
				{
					T value;
					std::memcpy( &value, bytes.data() + i*sizeof( T ), sizeof( T ) );
					out.push_back( static_cast< double >( value ) );
				}
			}
		}

		/// One DataArray of a `.vtu`, by its `Name`, as doubles in file order.
		/// The Points array has no name and is read as "Points". Empty when
		/// the file has no such array.
		inline std::vector< double > readVtuArray( std::string const &path,
		                                           std::string const &name )
		{
			std::ifstream stream( path, std::ios::binary );
			if ( !stream )
				throw std::runtime_error( "vtu: cannot open " + path );
			std::string const text( ( std::istreambuf_iterator< char >( stream ) ),
			                        std::istreambuf_iterator< char >() );
			bool const compressed =
				text.find( "compressor=\"vtkZLibDataCompressor\"" ) != std::string::npos;
			if ( text.find( "header_type=\"UInt64\"" ) != std::string::npos )
				throw std::runtime_error( "vtu: UInt64 headers are not read here" );

			std::size_t at = 0;
			while ( ( at = text.find( "<DataArray", at ) ) != std::string::npos )
			{
				std::size_t const close = text.find( '>', at );
				std::string const tag = text.substr( at, close - at + 1 );
				std::size_t const end = text.find( "</DataArray>", close );
				std::string arrayName = vtu::attributeOf( tag, "Name" );
				std::size_t const opened = text.rfind( "<Points>", at );
				std::size_t const closed = text.rfind( "</Points>", at );
				if ( arrayName.empty() && opened != std::string::npos
				     && ( closed == std::string::npos || closed < opened ) )
					arrayName = "Points";
				if ( arrayName == name )
				{
					if ( vtu::attributeOf( tag, "format" ) != "binary" )
						throw std::runtime_error( "vtu: " + name + " is not inline binary" );
					std::vector< unsigned char > const bytes =
						vtu::payload( text.substr( close + 1, end - close - 1 ), compressed );
					std::string const type = vtu::attributeOf( tag, "type" );
					std::vector< double > out;
					if ( type == "Float64" ) vtu::append< double >( bytes, out );
					else if ( type == "Float32" ) vtu::append< float >( bytes, out );
					else if ( type == "Int32" ) vtu::append< std::int32_t >( bytes, out );
					else if ( type == "Int64" ) vtu::append< std::int64_t >( bytes, out );
					else if ( type == "UInt32" ) vtu::append< std::uint32_t >( bytes, out );
					else if ( type == "UInt8" ) vtu::append< std::uint8_t >( bytes, out );
					else throw std::runtime_error( "vtu: type " + type + " is not read here" );
					return out;
				}
				at = end;
			}
			return {};
		}

		/// Whether @a name is written as POINT data (true) or CELL data.
		inline bool vtuArrayIsPointData( std::string const &path, std::string const &name )
		{
			std::ifstream stream( path, std::ios::binary );
			std::string const text( ( std::istreambuf_iterator< char >( stream ) ),
			                        std::istreambuf_iterator< char >() );
			std::size_t const at = text.find( "Name=\"" + name + "\"" );
			if ( at == std::string::npos )
				throw std::runtime_error( "vtu: no array " + name );
			std::size_t const point = text.rfind( "<PointData", at );
			std::size_t const cell = text.rfind( "<CellData", at );
			return point != std::string::npos
			       && ( cell == std::string::npos || point > cell );
		}
	}
}

#endif // MEQ_TESTS_VTU_READER_HPP
