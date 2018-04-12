#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <iostream>
#include "FileIO\fileIO.h"
#include "utils\loopmacros.h"

class Particle
{
protected:
	std::vector<std::vector<double>> origData_m;
	std::vector<std::vector<double>> currData_m;

	double*  origData1D_d{ nullptr };
	double*  currData1D_d{ nullptr };
	double** origData2D_d{ nullptr };
	double** currData2D_d{ nullptr };

	std::vector<std::string> attributeNames_m;

	long numberOfParticles_m;
	int  numberOfPositionDims_m;
	int  numberOfVelocityDims_m;

	double mass_m;
	double charge_m;
	double normFactor_m;

	std::string name_m;

	virtual void initializeGPU();

	bool initDataLoaded_m{ false };
	bool dataOnGPU_m{ false };

public:
	Particle(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor=1) :
		name_m{ name }, attributeNames_m{ attributeNames }, mass_m{ mass }, charge_m{ charge }, numberOfParticles_m{ numParts }, numberOfPositionDims_m{ posDims },
		numberOfVelocityDims_m{ velDims }, normFactor_m{ normFactor }
	{
		origData_m = std::vector<std::vector<double>>(posDims + velDims, std::vector<double>(numParts));
		currData_m = std::vector<std::vector<double>>(posDims + velDims, std::vector<double>(numParts));

		if (attributeNames.size() < numberOfVelocityDims_m + numberOfPositionDims_m)
		{
			std::cerr << "Particle::Particle: warning: not enough attribute names specified, generic names being generated" << std::endl;
			for (int diff = 0; diff = numberOfVelocityDims_m + numberOfPositionDims_m - (int)(attributeNames.size()); diff++)
				attributeNames_m.push_back(std::to_string(diff) + ".bin");
		}

		//#ifdef SOME_DEFINE_TO_INDICATE_COMPILE_WITH_CUDA
		initializeGPU();
		//copyDataToGPU(); //maybe??  Sometimes data is copied in after the fact...could set it that either a folder is specified in constructor or random particles are generated
		//#endif
	}
	~Particle()
	{ freeGPUMemory(); }

	std::vector<std::vector<double>>& getOrigData() { return origData_m; }
	std::vector<std::vector<double>>& getCurrData() { return currData_m; }
	std::vector<std::string>& getAttrNames() { return attributeNames_m; }
	std::string name() { return name_m; }
	double   mass() { return mass_m; }
	double   charge() { return charge_m; }
	long     getNumberOfParticles() { return numberOfParticles_m; }
	int      getNumberOfAttributes() { return numberOfPositionDims_m + numberOfVelocityDims_m; }
	bool     getInitDataLoaded() { return initDataLoaded_m; }
	double** getOrigDataGPUPtr() { return origData2D_d; }
	double** getCurrDataGPUPtr() { return currData2D_d; }

	virtual int getDimensionIndByName(std::string searchName);
	virtual std::string getDimensionNameByInd(int searchIndx);

	virtual void loadDataFromMem(std::vector<std::vector<double>> data, bool orig = true) { ((orig) ? origData_m = data : currData_m = data); numberOfParticles_m = ((orig) ? (int)origData_m.at(0).size() : (int)currData_m.at(0).size()); }
	virtual void loadDataFromDisk(std::string folder, bool orig = true);
	virtual void saveDataToDisk(std::string folder, bool orig);
	virtual void generateRandomParticles(const std::vector<double>& s, int startInd, int length, double vmean, double kBT_eV, double mass);
	
	virtual void copyDataToGPU();
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void clearGPUMemory();
};

#endif