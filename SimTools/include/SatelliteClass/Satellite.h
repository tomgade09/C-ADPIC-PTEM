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
	std::vector<double*> origDataGPU_m; //put GPU pointers to particle data here
	std::vector<double*> captureDataGPU_m; //double pointers to GPU arrays containing captured particle attributes - should be numberOfAttributes_m in size
	std::string name_m;

	virtual void initializeSatelliteOnGPU();
	virtual void freeGPUMemory();

public:
	Satellite(double altitude, bool upwardFacing, int numberOfAttributes, int numberOfParticles, double** dataPtrsGPU, std::string name = "Satellite"):
		altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfAttributes_m{ numberOfAttributes },
		numberOfParticles_m{ numberOfParticles }, name_m{ name }
	{
		data_m = new double*[numberOfAttributes_m];
		for (int iii = 0; iii < numberOfAttributes_m; iii++)
		{
			data_m[iii] = new double[numberOfParticles_m];
			for (int jjj = 0; jjj < numberOfParticles_m; jjj++)
				data_m[iii][jjj] = 0.0;
		}
		
		origDataGPU_m.reserve(numberOfAttributes_m);
		captureDataGPU_m.reserve(numberOfAttributes_m + 1);

		for (int iii = 0; iii < numberOfAttributes_m; iii++)
			origDataGPU_m.push_back(dataPtrsGPU[iii]);
		
		initializeSatelliteOnGPU();
	}
	
	virtual ~Satellite()
	{
		for (int iii = 0; iii < numberOfAttributes_m; iii++)
			delete[] data_m[iii];
		delete[] data_m;
		
		freeGPUMemory();
	}
	
	virtual void iterateDetector(int numberOfBlocks, int blockSize, double simtime); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??

	//Access functions
	double* getDataArrayPointer(int index) { dataReady_m = false; return data_m[index]; }
	double  getAltitude() { return altitude_m; }
	bool	getUpward() { return upwardFacing_m; }
	void	clearDataReady() { dataReady_m = false; }
	bool	getDataReady() { return dataReady_m; }

	void    vectorTest(std::vector<double*>& in);
};
#endif