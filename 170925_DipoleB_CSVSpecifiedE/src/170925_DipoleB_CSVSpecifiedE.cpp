#include "include\_simulationvariables.h"
#include "SimulationClass\SimulationAPI.h"
#include "include\Simulation170925.h"

DLLEXPORT double* getPointerToElectricFieldDataAPI(Simulation170925* simulation, int index) {
	return simulation->getPointerToElectricFieldData(index); }

DLLEXPORT Simulation* createSimulation170925(const char* rootdir) {
	Simulation* ret = new Simulation170925(2, NUMPARTICLES, 3, DT, rootdir, "ez.out");
	return ret; }

DLLEXPORT void terminateSimulation170925(Simulation170925* simulation) {
	delete simulation; }

#ifndef DLLFILE
int main()//defined in SimulationAPI.h and fileIO.h
{
	Simulation* sim;
	sim = createSimulation170925("./../../../");

	sim->initializeSimulation();
	sim->copyDataToGPU();
	sim->iterateSimulation(25000);
	sim->copyDataToHost();
	sim->freeGPUMemory();
	sim->prepareResults();

	delete sim;

	return 0;
}
#endif