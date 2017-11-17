#include "include\_simulationvariables.h"
#include "include\SimulationAPI.h"
#include "include\Simulation170925.h"

DLLEXPORT void resetParticlesEscapedCountAPI(Simulation170925* simulation) {
	simulation->resetParticlesEscapedCount(); }

DLLEXPORT double* getPointerToElectricFieldDataAPI(Simulation170925* simulation, int index) {
	return simulation->getPointerToElectricFieldData(index); }

DLLEXPORT Simulation* createSimulation170925(const char* rootdir) {
	Simulation* ret = new Simulation170925(2, NUMPARTICLES, 3, DT, rootdir, "ez.out");
	return ret; }

DLLEXPORT void terminateSimulation170925(Simulation170925* simulation) {
	delete simulation; }