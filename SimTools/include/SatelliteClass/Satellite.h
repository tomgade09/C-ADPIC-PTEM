#ifndef SATELLITE_H
#define SATELLITE_H
#include <vector>
#include <iostream>
#include <string>
#include "StandaloneTools\StandaloneTools.h"

class Satellite
{
protected:
	double altitude_m;
	bool upwardFacing_m;
	bool dataReady_m{ false };
	int numberOfAttributes_m;
	long numberOfParticles_m;
	std::vector<std::vector<std::vector<double>>> data_m; //[measurement][attribute][particle]
	std::vector<double**> dblppGPU_m; //double pointers to GPU arrays containing particle data from sim and captured particle data from satellite
	double* satCaptureGPU_m;
	std::string name_m;

	virtual void initializeSatelliteOnGPU();
	virtual void freeGPUMemory();
	virtual void dataAllocateNewMsmtVector() { data_m.push_back(form2DvectorArray(numberOfAttributes_m + 1, numberOfParticles_m)); }

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, long numberOfParticles, double** dataPtrsGPU, std::string name = "Satellite"):
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes }, numberOfParticles_m{ numberOfParticles }, name_m { name }
	{
		dblppGPU_m.resize(2);
		dblppGPU_m.at(0) = dataPtrsGPU;

		initializeSatelliteOnGPU();
	}
	
	virtual ~Satellite() { freeGPUMemory();	}
	
	virtual void iterateDetector(int blockSize, double simtime, double dt); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??

	//Access functions
	std::vector<std::vector<double>> getConsolidatedData(bool removeZeros)
	{
		std::vector<std::vector<double>> tmp2D;

		for (int attrs = 0; attrs < numberOfAttributes_m + 1; attrs++)
			tmp2D.push_back(std::vector<double>());

		LOOP_OVER_3D_ARRAY(data_m.size(), data_m.at(iii).size(), numberOfParticles_m, \
			if (removeZeros) //iii is msmt iterator, jjj is attribute iterator, kk is particle iterator
			{
				size_t tind{ data_m.at(iii).size() - 1 };
				if (data_m.at(iii).at(tind).at(kk) >= 0.0)
					tmp2D.at(jjj).push_back(data_m.at(iii).at(jjj).at(kk));
			}
			else
				tmp2D.at(jjj).push_back(data_m.at(iii).at(jjj).at(kk));
		)

		return tmp2D;
	}

	std::vector<std::vector<std::vector<double>>>& getDataVectorRef() { return data_m; }

	int     getNumberOfAttributes() { return numberOfAttributes_m; }
	long    getNumberOfParticles() { return numberOfParticles_m; }
	double  getAltitude() { return altitude_m; }
	bool	getUpward() { return upwardFacing_m; }
	//void	clearDataReady() { dataReady_m = false; }
	bool	getDataReady() { return dataReady_m; }
	std::string getName() { return name_m; }

	//void    vectorTest(std::vector<double*>& in);
};
#endif