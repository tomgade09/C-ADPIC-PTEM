#include "_simulationvariables.h"
#include "SimulationClass\SimulationAPI.h"
//#include "SimulationClass\AlfvenLUT.h"

//DLLEXPORT double* getPointerToElectricFieldDataAPI(AlfvenLUT* simulation, int index) {
	//return simulation->getPointerToElectricFieldData(index); }

DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, double ionT, double magT, const char* rootdir, double constEQSPS=0.0, const char* fnLUT="")
{
	std::string fnLUTstr{ fnLUT };

	Simulation* ret{ nullptr };

	ret = new Simulation(dt, simMin, simMax, ionT, magT, rootdir);
	
	if (constEQSPS != 0.0)
	{
		std::cout << "createSimulation: QSPS\n";
		ret->setQSPS(constEQSPS);
	}

	return ret;
}

DLLEXPORT void terminateSimulationAPI(Simulation* simulation) { //change
	delete simulation; }

#ifndef DLLFILE
int main()//defined in fileIO.h - need to update this
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