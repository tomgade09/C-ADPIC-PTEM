#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "SimulationClass\Simulation.h"
#include "FileIO\fileIO.h"

class Simulation170925 : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr };

public:
	Simulation170925(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(numberOfParticleTypes, numberOfParticlesPerType, numberOfAttributesTracked, dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Populate E Field LUT
		std::string LUT{ rootdir_m + "\\in\\" + LUTfilename_m };
		setElecMagLUT(LUT.c_str(), LUTNUMOFENTRS, LUTNUMOFCOLS);

		logFile_m.writeTimeDiffFromNow(0, "End Simulation170925 Constructor");
	}//end constructor

	~Simulation170925() //Destructor
	{
		//Delete double arrays
		for (int iii = 0; iii < 3; iii++)
			delete[] elcFieldLUT_m[iii];
		delete[] elcFieldLUT_m;
		logFile_m.writeTimeDiffFromNow(0, "End Simulation170925 Destructor");
	}

	//One liners
	double*  getPointerToElectricFieldData(int column) { if (elcFieldLUT_m == nullptr)
		{ std::cout << "Array not initialized yet.  Run Simulation::setElecMagLUT.\n"; return nullptr; } return elcFieldLUT_m[column]; }

	//Array tools
	virtual void   setElecMagLUT(const char* filename, int rows, int cols) { elcFieldLUT_m = fileIO::read2DCSV(filename, rows, cols, ' '); }
	virtual double calculateBFieldAtZandTime(double z, double time) { return BFieldatZ(z, time); }
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(elcFieldLUT_m, z, time); }
	virtual void   loadDataFilesIntoParticleArray();

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults();
};

#endif //end header guard