
/*
 * Configuration Wrapper Class for MEQ
 */

#include <vector>
#include "meq.hpp"

class Configuration {
	private:

	Domain ConfiguredDomain;
	std::vector<Coils> CoilSet;

	PlasmaModel plasma;
	
	std::string MeshFile,FinalMeshFile,PsiFile,GradPsiFile;
	
	public:
		Configuration( std::istream & is );
		Configuration( std::string const& f_name );
		Configuration( const char[] f_name );


		Domain const& GetDomain() {return Domain;};
		std::vector<Coils> const& GetCoils() {return CoilSet;};
		PlasmaModel const& GetPlasmaModel() {return plasma;};
		std::string const& GetMeshFile() { return MeshFile;};
		std::string const& GetFinalMeshFile() { return FinalMeshFile;};
		std::string const& GetPsiFile() { return PsiFile;};
		std::string const& GetGradPsiFile() { return GradPsiFile;};

}
