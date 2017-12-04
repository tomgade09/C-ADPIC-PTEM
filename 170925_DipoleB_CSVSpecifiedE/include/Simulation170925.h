#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "include\Simulation.h"
#include "include\fileIO.h"

double BFieldatZ(double z, double simtime);
double EFieldatZ(double** LUT, double z, double simtime);

class Simulation170925 : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr };
	
	long totalElecEscaped_m{ 0 };
	long totalIonsEscaped_m{ 0 };
	bool mu_m{ 1 };

	std::vector<double> mass_m;

public:
	Simulation170925(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(numberOfParticleTypes, numberOfParticlesPerType, numberOfAttributesTracked, dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Populate mass array
		mass_m.reserve(2);
		mass_m[0] = MASS_ELECTRON;
		mass_m[1] = MASS_PROTON;

		//Populate E Field LUT
		std::string LUT{ rootdir_m + "\\in\\" + LUTfilename_m };
		setElecMagLUT(LUT.c_str(), 2951, 3);

		logFile_m.writeLogTimeDiffFromNow(0, "End Simulation170925 Constructor");
	}//end constructor

	~Simulation170925() //Destructor
	{
		//Save final particle distributions to disk
		std::string fold{ "./particles_final/" };
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}

		//Delete double arrays
		for (int iii = 0; iii < 3; iii++)
			delete[] elcFieldLUT_m[iii];
		delete[] elcFieldLUT_m;
		
		//Delete satellites
		for (int iii = 0; iii < satellites_m.size(); iii++)
			delete satellites_m[iii];
	}

	//One liners
	double*  getPointerToElectricFieldData(int column) { if (elcFieldLUT_m == nullptr)
		{ std::cout << "Array not initialized yet.  Run Simulation::setElecMagLUT.\n"; return nullptr; } return elcFieldLUT_m[column]; }
	virtual void resetParticlesEscapedCount() { totalElecEscaped_m = 0; totalIonsEscaped_m = 0; return; }

	//Array tools
	virtual void setElecMagLUT(const char* filename, int rows, int cols) { elcFieldLUT_m = fileIO::read2DCSV(filename, rows, cols, ' '); }
	virtual double calculateBFieldAtZandTime(double z, double time) { return BFieldatZ(z, time); }
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(elcFieldLUT_m, z, time); }
	virtual void loadDataFilesIntoParticleArray();
	virtual void convertVPerpToMu();
	virtual void convertMuToVPerp();

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults();
};

#endif //end header guard