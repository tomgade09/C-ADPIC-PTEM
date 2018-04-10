#ifndef SATELLITE_H
#define SATELLITE_H
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include "FileIO\fileIO.h"
#include "utils\loopmacros.h"

class Satellite
{
protected:
	double altitude_m;
	bool upwardFacing_m;
	bool dataReady_m{ false };
	bool dataOnGPU_m{ false };

	int numberOfAttributes_m;
	long numberOfParticles_m;

	std::vector<std::vector<std::vector<double>>> data_m; //[measurement][attribute][particle]
	double*  satCaptrData1D_d{ nullptr }; //flattened satellite capture data
	double** satCaptrData2D_d{ nullptr }; //2D satellite capture data
	double** particleData2D_d{ nullptr };
	
	std::string name_m;

	virtual void initializeGPU();
	virtual void dataAllocateNewMsmtVector() { data_m.push_back(std::vector<std::vector<double>>(numberOfAttributes_m + 2, std::vector<double>(numberOfParticles_m))); }

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, long numberOfParticles, double** partDataGPUPtr, std::string name = "Satellite"):
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes }, numberOfParticles_m{ numberOfParticles }, particleData2D_d{ partDataGPUPtr }, name_m { name }
	{ initializeGPU(); }
	
	virtual ~Satellite() { freeGPUMemory(); }
	
	virtual void iterateDetector(double simtime, double dt, int blockSize); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??
	virtual void freeGPUMemory();
	virtual void saveDataToDisk(std::string folder, std::vector<std::string> attrNames);
	virtual void loadDataFromDisk(std::string folder, std::vector<std::string> attrNames);

	//Access functions
	std::vector<std::vector<double>> getConsolidatedData(bool removeZeros);
	const std::vector<std::vector<std::vector<double>>>& data() { return data_m; }
	double** get2DDataGPUPtr() { return satCaptrData2D_d; }
	double*  get1DDataGPUPtr() { return satCaptrData1D_d; }

	int     getNumberOfAttributes() { return numberOfAttributes_m; }
	long    getNumberOfParticles() { return numberOfParticles_m; }
	double  altitude() { return altitude_m; }
	bool	upward() { return upwardFacing_m; }
	std::string name() { return name_m; }
};
#endif