#include "include\_simulationvariables.h"
#include "SimulationClass\SimulationAPI.h"
#include "include\AlfvenLUT.h"

DLLEXPORT double* getPointerToElectricFieldDataAPI(AlfvenLUT* simulation, int index) {
	return simulation->getPointerToElectricFieldData(index); }

//
//
//
//need to restructure to create/delete depending on specified class/subclass type
DLLEXPORT Simulation* createSimulation170925(const char* rootdir) {
	Simulation* ret = new AlfvenLUT(DT, rootdir, "ez.out");
	ret->createParticleType("elec", { "vpara", "vperp", "z" }, MASS_ELECTRON, -1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	ret->createParticleType("ions", { "vpara", "vperp", "z" }, MASS_PROTON, 1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	
	Particle* elec{ ret->getParticlePointer(0) };
	Particle* ions{ ret->getParticlePointer(1) };
	elec->loadFilesToArray("./../../in/data/");
	ions->loadFilesToArray("./../../in/data/");

	return ret; }
//
//
//

//
//
//
//need to restructure to create/delete depending on specified class/subclass type
DLLEXPORT void terminateSimulation170925(AlfvenLUT* simulation) { //change
	delete simulation; }
//
//
//

#ifndef DLLFILE
int main()//defined in fileIO.h
{
	Simulation* sim;
	sim = createSimulation170925("./../../../");

	sim->initializeSimulation();
	sim->copyDataToGPU();
	sim->iterateSimulation(5000);
	sim->copyDataToHost();
	sim->freeGPUMemory();
	sim->prepareResults();

	delete sim;

	return 0;
}
#endif