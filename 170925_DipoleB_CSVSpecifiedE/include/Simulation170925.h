#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "SimulationClass\Simulation.h"
#include "FileIO\fileIO.h"

class Simulation170925 : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr }; //use instead of a 2D vector so I don't have to rewrite EFieldAtZ and have two copies (GPU and host)
	//std::vector<std::vector<double>> elcFieldLUT_m;

	virtual void initializeFollowOn() { return; }
	virtual void copyDataToGPUFollowOn() { return; }
	virtual void iterateSimulationFollowOnPreLoop() { return; }
	virtual void iterateSimulationFollowOnInsideLoop() { return; }
	virtual void iterateSimulationFollowOnPostLoop() { return; }
	virtual void copyDataToHostFollowOn() { return; }
	virtual void freeGPUMemoryFollowOn() { return; }

public:
	Simulation170925(double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Populate E Field LUT
		std::string LUT{ rootdir_m + "\\in\\" + LUTfilename_m };
		setElecMagLUT(LUT.c_str(), LUTNUMOFENTRS, LUTNUMOFCOLS);

		logFile_m.writeTimeDiffFromNow(0, "End Simulation170925 Constructor");
	}//end constructor

	~Simulation170925() //Destructor
	{
		//Delete double arrays
		delete[] elcFieldLUT_m[0];
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
	//virtual void   loadDataFilesIntoParticleArray();

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults();
};

#endif //end header guard