#ifndef MEQCONF_HPP
#define MEQCONF_HPP

/*
 * Configuration Wrapper Class for MEQ
 */

#include <vector>

#include "toml.hpp"

#include "utility.hpp"
#include "model.hpp"

namespace meq {

class Configuration {
	private:
		std::shared_ptr<Domain> ConfiguredDomain;
		std::vector<Coil> CoilSet;
		std::shared_ptr<PlasmaModel> plasma;
		std::string MeshFile,FinalMeshFile,PsiFile,GradPsiFile,NetCDFFile;

	public:
		Configuration( std::string f_name ) {
			const toml::value config = toml::parse( f_name );

			try {
				const auto options = toml::find< toml::table >( config, "options" );
				MeshFile = options.at( "MeshFile" ).as_string();
				FinalMeshFile = options.at( "FinalMeshFile" ).as_string();
				PsiFile = options.at( "PsiFile" ).as_string();
				GradPsiFile = options.at( "GradPsiFile" ).as_string();
				NetCDFFile = options.at( "NetCDFFile" ).as_string();
			} catch ( std::out_of_range &err ) {
				std::cerr << err.what() << std::endl;
				throw err;
			}
			
			const auto domainConfig = toml::find< toml::table >( config, "domain" );

			std::string BoundaryType;

			BoundaryType = domainConfig.at( "BoundaryCondition" ).as_string();
			std::map<std::string,BoundaryConditionType> BCTrans =
			{
				{"VonHagenow",VonHagenow},
				{"Prescribed",Prescribed},
				{"Zero",Zero}
			};

			BoundaryConditionType bcType;
			auto needle = BCTrans.find( BoundaryType );

			if ( needle != BCTrans.end() )
				bcType = needle->second;
			else
				bcType = Unknown;


			ConfiguredDomain = std::make_shared<Domain>( 
					domainConfig.at( "RMin" ).as_floating(), 
					domainConfig.at( "RMax" ).as_floating(), 
					domainConfig.at( "ZMin" ).as_floating(), 
					domainConfig.at( "ZMax" ).as_floating(), 
					domainConfig.at( "CellSize" ).as_floating(), 
					bcType );

			plasma = nullptr;

			const auto coils = toml::find< std::vector< toml::table > >( config, "coils" );
			for ( const auto & coil : coils )
			{
				double R = coil.at( "R" ).as_floating();
				double Z = coil.at( "Z" ).as_floating();
				double w = coil.at( "Width" ).as_floating();
				double h = coil.at( "Height" ).as_floating();
				double I = coil.at( "Current" ).as_floating();
				double coil_area = w*h;
				double j_coil = I/coil_area;
				CoilSet.emplace_back( R, Z, h, w, j_coil );
			}
		};

		~Configuration() { };

		std::shared_ptr<Domain>  GetDomain() const {return ConfiguredDomain;};
		std::vector<Coil> const& GetCoils() const {return CoilSet;};
		std::shared_ptr<PlasmaModel>  GetPlasmaModel() const {return plasma;};
		std::string const& GetMeshFile() const { return MeshFile;};
		std::string const& GetFinalMeshFile() const { return FinalMeshFile;};
		std::string const& GetPsiFile() const { return PsiFile;};
		std::string const& GetGradPsiFile() const { return GradPsiFile;};
		std::string const& GetNetCDFFile() const { return NetCDFFile;};

};

};


#endif // MEQCONF_HPP
