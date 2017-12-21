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
	bool elecTF_m;
	bool dataReady_m{ false };
	int numberOfAttributes_m;
	int numberOfParticles_m;
	double** data_m{ nullptr };
	std::vector<double**> dblppGPU_m; //double pointers to GPU arrays containing particle data from sim and captured particle data from satellite
	double* satCaptureGPU_m;
	std::string name_m;

	virtual void initializeSatelliteOnGPU();
	virtual void freeGPUMemory();
	
	virtual void allocateData()
	{
		double*** tmp{ form3Darray(1, numberOfAttributes_m + 1, numberOfParticles_m) };
		data_m = tmp[0]; //only need two dimensions - capture pointer to first (only) 2D array within 3D
		delete[] tmp; //discard 3rd dimension
	}

	virtual void deleteData()
	{
		delete[] data_m[0];
		delete[] data_m;
	}

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, int numberOfParticles, double** dataPtrsGPU, bool elecTF, std::string name = "Satellite"):
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes },
		numberOfParticles_m{ numberOfParticles }, elecTF_m{elecTF}, name_m { name }
	{
		allocateData();

		dblppGPU_m.reserve(2);
		dblppGPU_m[0] = dataPtrsGPU;

		initializeSatelliteOnGPU();
	}
	
	virtual ~Satellite()
	{
		deleteData();
		freeGPUMemory();
	}
	
	virtual void iterateDetector(int numberOfBlocks, int blockSize, double simtime); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??

	//Access functions
	double** getDataArrayPointer(bool releaseArray)
	{
		dataReady_m = false;
		if (releaseArray)
		{
			double** ptr{ data_m };
			allocateData();
			return ptr;
		}
		return data_m;
	}

	double  getAltitude() { return altitude_m; }
	bool	getUpward() { return upwardFacing_m; }
	void	clearDataReady() { dataReady_m = false; }
	bool	getDataReady() { return dataReady_m; }
	bool	getElecTF() { return elecTF_m; }
	std::string getName() { return name_m; }

	void    vectorTest(std::vector<double*>& in);
};
#endif