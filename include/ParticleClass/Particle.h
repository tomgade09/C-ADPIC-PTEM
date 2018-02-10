#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <iostream>
#include "FileIO\fileIO.h"
#include "StandaloneTools\StandaloneTools.h"

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

	long particleCount_m;
	int  numberOfPositionDims_m;
	int  numberOfVelocityDims_m;

	double mass_m;
	double charge_m;
	double normFactor_m;

	std::string name_m;

	bool initDataLoaded_m{ false };
	bool usingGPU_m{ false };
	bool normalized_m{ false };

public:
	Particle(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor=1) :
		name_m{ name }, attributeNames_m{ attributeNames }, mass_m{ mass }, charge_m{ charge }, particleCount_m{ numParts }, numberOfPositionDims_m{ posDims },
		numberOfVelocityDims_m{ velDims }, normFactor_m{ normFactor }
	{
		std::vector<double> tmp;
		tmp.resize(numParts);
		for (int dataind = 0; dataind < velDims + posDims; dataind++)
		{
			origData_m.push_back(tmp);
			currData_m.push_back(tmp);
		}

		if (attributeNames.size() < numberOfVelocityDims_m + numberOfPositionDims_m)
		{
			std::cerr << "Particle::Particle: warning: not enough attribute names specified, generic names being generated" << std::endl;
			for (int diff = 0; diff = numberOfVelocityDims_m + numberOfPositionDims_m - static_cast<int>(attributeNames.size()); diff++)
				attributeNames_m.push_back(std::to_string(diff) + ".bin");
		}
	}
	~Particle()
	{
		if (usingGPU_m)
			freeGPUMemory();
	}

	std::vector<std::vector<double>>& getOrigData() { return origData_m; }
	std::vector<std::vector<double>>& getCurrData() { return currData_m; }
	std::vector<std::string>& getAttrNames() { return attributeNames_m; }
	std::string getName() { return name_m; }
	double   getMass() { return mass_m; }
	double   getCharge() { return charge_m; }
	long     getNumberOfParticles() { return particleCount_m; }
	int      getNumberOfAttributes() { return numberOfPositionDims_m + numberOfVelocityDims_m; }
	bool     getInitDataLoaded() { return initDataLoaded_m; }
	double** getOrigDataGPUPtr() { return origData2D_d; }
	double** getCurrDataGPUPtr() { return currData2D_d; }

	virtual int getDimensionIndByName(std::string searchName);
	virtual std::string getDimensionNameByInd(int searchIndx);

	virtual void loadFilesToArray(std::string folder, bool orig=false);
	virtual void saveArrayToFiles(std::string folder, bool orig);
	virtual void normalizeParticles(bool orig, bool curr, bool inverse=false);

	virtual void initializeGPU();
	virtual void copyDataToGPU();
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
};

#endif