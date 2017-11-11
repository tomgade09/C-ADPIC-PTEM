#ifndef SATELLITE_H
#define SATELLITE_H
#include <vector>
#include <iostream>
#include <string>

class Satellite
{
protected:
	double altitude_m;
	bool upwardFacing_m;
	bool dataReady_m{ false };
	int numberOfAttributes_m;
	int numberOfParticles_m;
	double** data_m;
	std::vector<double*> GPUdata_m;
	std::string name_m;

	virtual void initializeSatelliteOnGPU();
	virtual void freeGPUMemory();

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, int numberOfParticles, std::string name = "Satellite") :
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes }, numberOfParticles_m{ numberOfParticles_m },
		name_m{ name }
	{
		data_m = new double*[numberOfAttributes_m];
		for (int iii = 0; iii < numberOfAttributes_m; iii++)
			data_m[iii] = new double[numberOfParticles_m];
		GPUdata_m.reserve(numberOfAttributes_m);
		
		initializeSatelliteOnGPU();
	}
	
	virtual ~Satellite()
	{
		for (int iii = 0; iii < numberOfAttributes_m; iii++)
			delete[] data_m[iii];
		delete[] data_m;
		
		freeGPUMemory();
	}
	
	virtual void iterateDetector(int numberOfBlocks, int blockSize, double** simData); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??

	//Access functions
	double* getDataArrayPointer(int index) { return data_m[index]; }
	double  getAltitude() { return altitude_m; }
	bool	getUpward() { return upwardFacing_m; }
	void	clearDataReady() { dataReady_m = false; }
	bool	getDataReady() { return dataReady_m; }
};
#endif