//To-Do
/*
	//- Check CUDA kernels for correctness
	//- Simulation variables, physical constants files (header?)
	!! Needed? Just provide to python - Correct E_B graph data calculator
	//- Maybe need to derive a class
	//- Main body of code (DLL, API) to initialize Simulation instance, run code
	//- Build Simulation::
	//>> serialzeParticleArray
	//>> setElecMagLUT
	//>> calculateFieldsAtTime
	//>> returnResults
	//>> ~Simulation
	//- Modify Python script
	//- One last look for correctness
- Test code for similarity to last sim

- Code to read CSV into array and save into simulation
- New E, B calculator based on new format
*/
#include "include\_simulationvariables.h"
#include "include\api.h"
#include "include\Simulation170925.h"

//double EFieldatZ(double** LUT, double z, double simtime);
double EFieldatZ(double z, double simtime);
double BFieldatZ(double z, double simtime);

DLLEXPORT void resetParticlesEscapedCountWrapper(Simulation170925* simulation)
{
	simulation->resetParticlesEscapedCount();
}

DLLEXPORT double** getPointerToElectricFieldDataWrapper(Simulation170925* simulation)
{
	return simulation->getPointerToElectricFieldData();
}

DLLEXPORT Simulation* createSimulation170925(const char* rootdir)
{
	Simulation* ret = new Simulation170925(2, NUMPARTICLES, 3, DT, rootdir, "ez.out");

	return ret;
}


DLLEXPORT void terminateSimulation170925(Simulation170925* simulation)
{
	delete simulation;
}