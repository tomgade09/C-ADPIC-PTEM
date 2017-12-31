#include "include\_simulationvariables.h"
#include "SimulationClass\SimulationAPI.h"
#include "include\Simulation170925.h"

DLLEXPORT double* getPointerToElectricFieldDataAPI(Simulation170925* simulation, int index) {
	return simulation->getPointerToElectricFieldData(index); }

DLLEXPORT Simulation* createSimulation170925(const char* rootdir) {
	Simulation* ret = new Simulation170925(DT, rootdir, "ez.out");
	ret->createParticleType("elec", { "vpara", "vperp", "z" }, MASS_ELECTRON, -1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	ret->createParticleType("ions", { "vpara", "vperp", "z" }, MASS_ELECTRON, -1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	ret->getParticlePointer(0)->loadFilesToArray("./../../in/data/", 100352);
	ret->getParticlePointer(1)->loadFilesToArray("./../../in/data/", 100352);
	return ret; }

DLLEXPORT void terminateSimulation170925(Simulation170925* simulation) {
	delete simulation; }

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