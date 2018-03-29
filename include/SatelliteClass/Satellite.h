#ifndef SATELLITE_H
#define SATELLITE_H
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include "StandaloneTools\StandaloneTools.h"
#include "FileIO\fileIO.h"

class Satellite
{
protected:
	double altitude_m;
	bool upwardFacing_m;
	bool dataReady_m{ false };
	bool dataOnGPU_m{ true };

	int numberOfAttributes_m;
	long numberOfParticles_m;

	std::vector<std::vector<std::vector<double>>> data_m; //[measurement][attribute][particle]
	std::vector<double**> dblppGPU_m; //double pointers to GPU arrays containing particle data from sim and captured particle data from satellite
	double* satCaptureGPU_m;
	
	std::string name_m;

	virtual void initializeSatelliteOnGPU();
	virtual void dataAllocateNewMsmtVector() { data_m.push_back(std::vector<std::vector<double>>(numberOfAttributes_m + 2, std::vector<double>(numberOfParticles_m)));
		/*data_m.push_back(form2DvectorArray(numberOfAttributes_m + 2, numberOfParticles_m));*/ }

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, long numberOfParticles, double** dataPtrsGPU, std::string name = "Satellite"):
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes }, numberOfParticles_m{ numberOfParticles }, name_m { name }
	{
		dblppGPU_m.resize(2);
		dblppGPU_m.at(0) = dataPtrsGPU;

		initializeSatelliteOnGPU();
	}
	
	virtual ~Satellite() { freeGPUMemory(); }
	
	virtual void iterateDetector(double simtime, double dt, int blockSize); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??
	virtual void freeGPUMemory();
	virtual void saveDataToDisk(std::string folder, std::vector<std::string> attrNames, std::function<double(double, double)> BatS, double mass);

	//Access functions
	std::vector<std::vector<double>> getConsolidatedData(bool removeZeros);
	std::vector<std::vector<std::vector<double>>>& getDataVectorRef() { return data_m; }

	int     getNumberOfAttributes() { return numberOfAttributes_m; }
	long    getNumberOfParticles() { return numberOfParticles_m; }
	double  getAltitude() { return altitude_m; }
	bool	getUpward() { return upwardFacing_m; }
	bool	getDataReady() { return dataReady_m; }
	std::string getName() { return name_m; }
};
#endif