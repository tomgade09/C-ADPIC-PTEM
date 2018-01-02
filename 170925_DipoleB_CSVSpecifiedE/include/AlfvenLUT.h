#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "SimulationClass\Simulation.h"
#include "FileIO\fileIO.h"

class AlfvenLUT : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr }; //use instead of a 2D vector so I don't have to rewrite EFieldAtZ and have two copies (GPU and host)

	double omegaE_m{ 20 * PI }; //angular frequency, 10 Hz wave, matches ez.out
	int numOfColsLUT_m{ 3 }; //pass in
	int numOfEntrLUT_m{ 2951 }; //pass in

	virtual void initializeFollowOn();
	virtual void copyDataToGPUFollowOn();
	virtual void iterateSimulationFollowOnPreLoop();
	virtual void iterateSimulationFollowOnInsideLoop();
	virtual void iterateSimulationFollowOnPostLoop();
	virtual void copyDataToHostFollowOn();
	virtual void freeGPUMemoryFollowOn();

public:
	AlfvenLUT(double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Populate E Field LUT
		std::string LUT{ rootdir_m + "\\in\\" + LUTfilename_m };
		setElecMagLUT(LUT.c_str(), numOfEntrLUT_m, numOfColsLUT_m);

		logFile_m.writeTimeDiffFromNow(0, "End Simulation170925 Constructor");
	}//end constructor

	~AlfvenLUT() //Destructor
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
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(elcFieldLUT_m, z, time, omegaE_m, useQSPS_m, useAlfLUT_m); }
};

#endif //end header guard