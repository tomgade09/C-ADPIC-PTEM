#include "include\_simulationvariables.h"
#include "SimulationClass\SimulationAPI.h"
#include "include\AlfvenLUT.h"

DLLEXPORT double* getPointerToElectricFieldDataAPI(AlfvenLUT* simulation, int index) {
	return simulation->getPointerToElectricFieldData(index); }

//
//
//
//need to restructure to create/delete depending on specified class/subclass type
DLLEXPORT Simulation* createSimulation(double dt, double simMin, double simMax, double ionT, double magT, const char* rootdir, double constEQSPS=0.0, const char* fnLUT="") {
	Simulation* ret{ nullptr };

	if (fnLUT == "")
		ret = new Simulation(dt, simMin, simMax, ionT, magT, rootdir);
	else
		ret = new AlfvenLUT(dt, simMin, simMax, ionT, magT, rootdir, fnLUT);
	
	if (constEQSPS != 0.0)
		ret->setQSPS(constEQSPS);
	
	//Simulation* ret = new AlfvenLUT(dt, rootdir, "ez.out");
	//ret->createParticleType("elec", { "vpara", "vperp", "z" }, MASS_ELECTRON, -1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	//ret->createParticleType("ions", { "vpara", "vperp", "z" }, MASS_PROTON, 1.0 * CHARGE_ELEM, 100352, 1, 2, RADIUS_EARTH);
	
	//Particle* elec{ ret->getParticlePointer(0) };
	//Particle* ions{ ret->getParticlePointer(1) };
	//elec->loadFilesToArray("./../../in/data/");
	//ions->loadFilesToArray("./../../in/data/");

	return ret; }
//
//
//

//
//
//
//need to restructure to create/delete depending on specified class/subclass type
DLLEXPORT void terminateSimulation(Simulation* simulation) { //change
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