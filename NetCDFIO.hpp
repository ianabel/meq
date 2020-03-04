
/*
 * Class for storing psi & grad psi in NetCDF file
 */

class NetCDFConfig
{
	public:
		NetCDFConfig( const std::string &file, std::vector<double> const& R_in, std::vector<double> const& Z_in ) :
			filename( file ),
			R_dim( R_in ),
			Z_dim( Z_in )
		{
		}

		std::string filename;
		std::vector<double> R_dim,Z_dim;
}
