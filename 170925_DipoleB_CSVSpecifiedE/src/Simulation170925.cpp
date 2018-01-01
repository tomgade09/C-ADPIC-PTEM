#include "include\Simulation170925.h"

void Simulation170925::prepareResults()
{
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), convertMuToVPerp(particleTypes_m.at(iii));)

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_init/", true);)
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_final/", false);)

	//normalizes m to Re
	if (normalizedToRe_m)
	{
		//std::cout << "prepareResults: vpara: " << particleTypes_m.at(0)->getOrigData().at(0).at(0) << " vperp: " << particleTypes_m.at(0)->getOrigData().at(1).at(0) << " z: " << particleTypes_m.at(0)->getOrigData().at(2).at(0) << "\n";
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->normalizeParticles(true, true);)
		//std::cout << "                vpara: " << particleTypes_m.at(0)->getOrigData().at(0).at(0) << " vperp: " << particleTypes_m.at(0)->getOrigData().at(1).at(0) << " z: " << particleTypes_m.at(0)->getOrigData().at(2).at(0) << "\n";

		for (int lll = 0; lll < satelliteData_m.size(); lll++)
		{//loop over number of measurements
			LOOP_OVER_1D_ARRAY(satellites_m.size(), convertMuToVPerp(satelliteData_m.at(lll).at(iii).at(1), satelliteData_m.at(lll).at(iii).at(2), particleTypes_m.at(iii%2)->getMass());)
			LOOP_OVER_3D_ARRAY(satellites_m.size(), satellites_m.at(iii)->getNumOfAttr(), satellites_m.at(iii)->getNumOfParts(), satelliteData_m.at(lll).at(iii).at(jjj).at(kk) /= RADIUS_EARTH;)
		}
	}
	
	resultsPrepared_m = true;
}