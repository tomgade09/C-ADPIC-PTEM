#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "include\Simulation.h"
#include "include\fileIO.h"

double BFieldatZ(double z, double simtime);
//double EFieldatZ(double** LUT, double z, double simtime);
double EFieldatZ(double z, double simtime);

class Simulation170925 : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr };
	
	long totalElecEscaped_m{ 0 };
	long totalIonsEscaped_m{ 0 };
	bool mu_m{ 0 };

	std::vector<double> mass_m;

public:
	Simulation170925(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(numberOfParticleTypes, numberOfParticlesPerType, numberOfAttributesTracked, dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Generate particles first
		//2 particle types (electrons and ions) with v_para, mu, and z tracked, generate v_para and v_perp (eventually becomming mu) normally
		double means[4];
		double sigmas[4];
		means[0] = V_DIST_MEAN;
		means[1] = V_DIST_MEAN;
		means[2] = V_DIST_MEAN;
		means[3] = V_DIST_MEAN;
		sigmas[0] = V_SIGMA;
		sigmas[1] = V_SIGMA;
		sigmas[2] = V_SIGMA;
		sigmas[3] = V_SIGMA;
		generateNormallyDistributedValues(2, means, sigmas);

		std::string LUT;
		LUT = rootdir_m + "\\in\\" + LUTfilename_m;
		setElecMagLUT(LUT.c_str(), 2951, 3);

		std::string importdir{ rootdir_m };
		importdir = importdir + "\\in\\data\\";
		std::vector<std::vector<std::string>> files;
		files.resize(2);
		files[0].resize(3);
		files[1].resize(3);
		files[0][0] = "e_vpara.bin";
		files[0][1] = "e_vperp.bin";
		files[0][2] = "e_z.bin";
		files[1][0] = "i_vpara.bin";
		files[1][1] = "i_vperp.bin";
		files[1][2] = "i_z.bin";

		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				loadFileIntoParticleAttribute(iii, jjj, importdir.c_str(), files[iii][jjj].c_str());
		}

		std::string fold;
		fold = "./particles_init/";
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}
	}//end constructor
	~Simulation170925()
	{
		std::string fold;
		fold = "./particles_final/";
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}
	}

	//One liners
	double**  getPointerToElectricFieldData() { if (elcFieldLUT_m == nullptr) { std::cout << "Array not initialized yet.  Run Simulation::setElecMagLUT.\n"; } return elcFieldLUT_m; }
	virtual void resetParticlesEscapedCount() { totalElecEscaped_m = 0; totalIonsEscaped_m = 0; return; }

	//Array tools
	virtual void setElecMagLUT(const char* filename, int rows, int cols);
	virtual double calculateBFieldAtZandTime(double z, double time);
	virtual double calculateEFieldAtZandTime(double z, double time);
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