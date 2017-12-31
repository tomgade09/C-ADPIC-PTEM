#include "include\Simulation170925.h"

/*void Simulation170925::loadDataFilesIntoParticleArray()//legacy - not used and may not work under new structure - namely particles_m isn't copied to the GPU first anymore
{//not super generic, but ok - the function is virtual
	std::string importdir{ rootdir_m };
	importdir = importdir + "\\in\\data\\";
	std::vector<std::vector<std::string>> files;
	files.resize(2);
	files.at(0).resize(3);
	files.at(1).resize(3);
	files.at(0) = { "e_vpara.bin", "e_vperp.bin", "e_z.bin" };
	files.at(1) = { "i_vpara.bin", "i_vperp.bin", "i_z.bin" };

	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, fileIO::readDblBin(particles_m.at(iii).at(jjj), importdir + files.at(iii).at(jjj), numberOfParticlesPerType_m);)
}*/

void Simulation170925::prepareResults()
{
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(),\
		convertMuToVPerp(particleTypes_m.at(iii)->getCurrData().at(1), particleTypes_m.at(iii)->getCurrData().at(2), particleTypes_m.at(iii)->getMass());)

	//normalizes m to Re
	if (normalizedToRe_m)
	{
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->normalizeParticles(true, true);)
		
		for (int lll = 0; lll < satelliteData_m.size(); lll++)
		{//loop over number of measurements
			LOOP_OVER_3D_ARRAY(satellites_m.size(), satellites_m.at(iii)->getNumOfAttr(), satellites_m.at(iii)->getNumOfParts(), satelliteData_m.at(lll).at(iii).at(jjj).at(kk) /= RADIUS_EARTH;)
		}
	}
	
	resultsPrepared_m = true;
}