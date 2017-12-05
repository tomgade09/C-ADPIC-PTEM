#include "include\Simulation170925.h"

void Simulation170925::loadDataFilesIntoParticleArray()
{//not super generic, but ok - the function is virtual
	std::string importdir{ rootdir_m };
	importdir = importdir + "\\in\\data\\";
	std::vector<std::vector<std::string>> files;
	files.resize(2);
	files[0].resize(3);
	files[1].resize(3);
	files[0] = { "e_vpara.bin", "e_vperp.bin", "e_z.bin" };
	files[1] = { "i_vpara.bin", "i_vperp.bin", "i_z.bin" };

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			loadFileIntoParticleAttribute(particles_m[iii][jjj], numberOfParticlesPerType_m, importdir.c_str(), files[iii][jjj].c_str());
	}
}

void Simulation170925::prepareResults()
{
	convertMuToVPerp();
	resultsPrepared_m = true;
}