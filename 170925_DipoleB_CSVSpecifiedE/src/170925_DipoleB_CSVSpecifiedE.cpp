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
- One last look for correctness
- Modify Python script
- Test code for similarity to last sim

- Code to read CSV into array and save into simulation
- New E, B calculator based on new format
*/
#include "include\_simulationvariables.h"
#include "include\api.h"
#include "include\Simulation170925.h"

double EFieldatZ(double** LUT, double z, double simtime);
double BFieldatZ(double z, double simtime);

DLLEXPORT double** getPointerToElectricFieldDataWrapper(Simulation170925* simulation);
//DLLEXPORT double** getPointerToMagneticFieldDataWrapper(Simulation* simulation);

DLLEXPORT double getEatZ(double** LUT, double z, double simtime)
{
	return EFieldatZ(nullptr, z, simtime);
}

DLLEXPORT double getBatZ(double z, double simtime)
{
	return BFieldatZ(z, simtime);
}

DLLEXPORT double** getPointerToElectricFieldDataWrapper(Simulation170925* simulation)
{
	return simulation->getPointerToElectricFieldData();
}

/*DLLEXPORT double** getPointerToMagneticFieldDataWrapper(Simulation170925* simulation)
{
	return simulation->getPointerToMagneticFieldData();
}*/

DLLEXPORT Simulation* newSimInitialize()
{
	Simulation* ret = new Simulation170925(2, NUMPARTICLES, 3, DT, "ez.out");

	return ret;
}