#include "include\Simulation170925.h"

void Simulation170925::loadDataFilesIntoParticleArray()//legacy - not used and may not work under new structure - namely particles_m isn't copied to the GPU first anymore
{//not super generic, but ok - the function is virtual
	std::string importdir{ rootdir_m };
	importdir = importdir + "\\in\\data\\";
	std::vector<std::vector<std::string>> files;
	files.resize(2);
	files[0].resize(3);
	files[1].resize(3);
	files[0] = { "e_vpara.bin", "e_vperp.bin", "e_z.bin" };
	files[1] = { "i_vpara.bin", "i_vperp.bin", "i_z.bin" };

	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, loadFileIntoParticleAttribute(particles_m[iii][jjj], numberOfParticlesPerType_m, importdir.c_str(), files[iii][jjj].c_str());)
}

void Simulation170925::prepareResults()
{
	convertMuToVPerp();

	//normalizes m to Re
	if (normalizedToRe_m)
	{
		LOOP_OVER_3D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, numberOfParticlesPerType_m, particles_m[iii][jjj][kk] /= RADIUS_EARTH;)
		LOOP_OVER_3D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, numberOfParticlesPerType_m, particlesorig_m[iii][jjj][kk] /= RADIUS_EARTH;)
		for (int lll = 0; lll < satelliteData_m.size(); lll++)
		{//loop over number of measurements
			LOOP_OVER_3D_ARRAY(satellites_m.size(), numberOfAttributesTracked_m, numberOfParticlesPerType_m, satelliteData_m[lll][iii][jjj][kk] /= RADIUS_EARTH;)
		}
	}
	
	resultsPrepared_m = true;
}