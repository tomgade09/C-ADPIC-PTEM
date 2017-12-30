#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <iostream>
#include "FileIO\fileIO.h"

class Particle
{
protected:
	std::vector<std::vector<double>> origData_m;
	std::vector<std::vector<double>> currData_m;

	long particleCount_m;
	int  numberOfPositionDims_m;
	int  numberOfVelocityDims_m;

	double mass_m;
	double charge_m;
	double normFactor_m;

	std::string name_m;

	bool dataIsNormalized_m;
	bool initDataLoaded_m;

public:
	Particle(std::string name, double mass, double charge, long numParts, int posDims, int velDims, double normFactor=1, bool normalized=false) :
		name_m{ name }, mass_m{ mass }, charge_m{ charge }, particleCount_m{ numParts }, numberOfPositionDims_m{ posDims }, numberOfVelocityDims_m{ velDims },
		normFactor_m{ normFactor }, dataIsNormalized_m{ normalized }
	{
		std::vector<double> tmp;
		tmp.resize(numParts);
		for (int dataind = 0; dataind < velDims + posDims; dataind++)
		{
			origData_m.push_back(tmp);
			currData_m.push_back(tmp);
		}
	}
	~Particle(){}

	std::vector<std::vector<double>>& getOrigData() { return origData_m; }
	std::vector<std::vector<double>>& getCurrData() { return currData_m; }

	virtual void loadFileToArray(std::vector<std::string> filenames, long numOfElements);
	virtual void saveArrayToFile(std::vector<std::string> filenames, long numOfElements);
	virtual void normalizeParticles(bool orig, bool curr, bool divide=true);
};

#endif