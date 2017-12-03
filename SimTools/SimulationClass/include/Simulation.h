#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include "include\Satellite.h"
#include "include\_simulationvariables.h" //remove this later, replace with passed in variables

struct timeStruct
{
	std::string label;
	std::chrono::steady_clock::time_point tp;
};

class Simulation
{
protected:
	//Simulation Characteristics
	double    simTime_m{ 0 };
	const double dt_m;
	const int numberOfParticleTypes_m;
	const int numberOfParticlesPerType_m;
	const int numberOfAttributesTracked_m;
	std::string rootdir_m;

	//Data Array Pointers
	double*** particles_m{ nullptr };
	bool**    particlesInSim_m{ nullptr };
	int**     particlesEscaped_m{ nullptr };
	std::vector<void*> otherMemoryPointers_m;
	std::vector<Satellite*> satellites_m;
	std::vector<std::vector<std::vector<double*>>> satelliteData_m; //4D satelliteData[recorded measurement number][satellite number][attribute number][particle number]
	std::vector<timeStruct*> timeStructs_m;
	//std::ofstream logFile_m; //more to come
	
	//GPU Memory Pointers
	std::vector<double*> gpuDblMemoryPointers_m { nullptr };
	std::vector<bool*> gpuBoolMemoryPointers_m{ nullptr };
	std::vector<int*> gpuIntMemoryPointers_m { nullptr };
	std::vector<void*> gpuOtherMemoryPointers_m{ nullptr };

	//Calculated Quantities, Flags, Parameters defined in Header
	bool initialized_m{ 0 };
	bool copied_m{ 0 };
	bool resultsPrepared_m{ 0 };
	bool normalizedToRe_m{ (NORMFACTOR > 2.0) ? true : false };
	double simMin_m{ MIN_Z_SIM };
	double simMax_m{ MAX_Z_SIM };
	
public:
	Simulation(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir):
		numberOfParticleTypes_m{ numberOfParticleTypes }, numberOfParticlesPerType_m{ numberOfParticlesPerType },
		numberOfAttributesTracked_m{ numberOfAttributesTracked }, dt_m{ dt }, rootdir_m{ rootdir }
	{
		timeStructs_m.reserve(20);
		timeStructs_m.push_back(createTimeStruct("Begin Simulation Constructor")); //index 0
		
		//Form the particle array, boolean "inSim" array, escaped particle count array, and set "inSim" to true for each particle
		particles_m = new double**[numberOfParticleTypes_m];
		particlesInSim_m = new bool*[numberOfParticleTypes_m];
		particlesEscaped_m = new int*[numberOfParticleTypes_m];
		
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			particles_m[iii] = new double*[numberOfAttributesTracked_m];
			particlesInSim_m[iii] = new bool[numberOfParticlesPerType_m];
			particlesEscaped_m[iii] = new int[numberOfParticlesPerType_m];
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			{
				particles_m[iii][jjj] = new double[numberOfParticlesPerType_m];
				for (int kk = 0; kk < numberOfParticlesPerType_m; kk++)
				{//values to initialize array
					particlesInSim_m[iii][kk] = true;
					particlesEscaped_m[iii][kk] = 0;
					particles_m[iii][jjj][kk] = 0.0;
				}
			}
		}

		//Allocate room in vectors for GPU Memory Pointers
		gpuDblMemoryPointers_m.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m); //holds pointers to GPU memory for particle attributes
		gpuBoolMemoryPointers_m.reserve(numberOfParticleTypes_m); //holds pointers to GPU memory for bool flag of whether particle is in sim or not
		gpuIntMemoryPointers_m.reserve(numberOfParticleTypes_m); //holds pointers to GPU memory for count of escaped particles
	}

	virtual ~Simulation()
	{
		for (int iii = 0; iii < timeStructs_m.size(); iii++)
			delete timeStructs_m[iii];

		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			{
				delete[] particles_m[iii][jjj];
			}
			delete[] particles_m[iii];
			delete[] particlesInSim_m[iii];
			delete[] particlesEscaped_m[iii];
		}
		delete[] particles_m;
		delete[] particlesInSim_m;
		delete[] particlesEscaped_m;

		for (int iii = 0; iii < satelliteData_m.size(); iii++) //number of measurements
		{
			for (int jjj = 0; jjj < satelliteData_m[0].size(); jjj++) //number of satellites
			{
				for (int kk = 0; kk < satelliteData_m[0][0].size(); kk++) //number of attributes
						delete[] satelliteData_m[iii][jjj][kk];
			}
		}

		//for (int iii = 0; iii < timeStructs_m.size(); iii++)
			//delete timeStructs_m[iii];
	};
	//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double getTime() { return simTime_m; }//tested
	double getdt() { return dt_m; };
	void incTime() { simTime_m += dt_m; }
	int getNumberOfParticleTypes() { return numberOfParticleTypes_m; }//tested
	int getNumberOfParticlesPerType() { return numberOfParticlesPerType_m; }//tested
	int getNumberOfAttributesTracked() { return numberOfAttributesTracked_m; }//tested
	bool areResultsPrepared() { return resultsPrepared_m; }

	bool getNormalized(){ return normalizedToRe_m; }

	int getNumberOfSatellites() { return satellites_m.size(); }
	int getNumberOfSatelliteMsmts() { return satelliteData_m.size(); }
	double* getSatelliteDataPointers(int measurementInd, int satelliteInd, int attributeInd) { return satelliteData_m[measurementInd][satelliteInd][attributeInd]; }

	//Pointer one liners
	double*** getPointerTo3DParticleArray() { return particles_m; }//tested
	double**  getPointerToSingleParticleTypeArray(int index) { if (index >= numberOfParticleTypes_m) { return nullptr; } return particles_m[index]; } //eventually print a warning so the user knows
	bool*     getPointerToParticlesInSimArray(int index) { if (index > numberOfParticleTypes_m) {std::cout << "Error: No particle exists for index requested.  ";
		std::cout << "There are only " << numberOfParticleTypes_m << " particle types in the simulation.\n"; return nullptr;} return particlesInSim_m[index]; }
	double*   getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex) { if ((partIndex > numberOfParticleTypes_m) || (attrIndex > numberOfAttributesTracked_m)) { 
		std::cout << numberOfParticleTypes_m << " particle types and " << numberOfAttributesTracked_m << " attributes per particle.  One index is out of bounds.\n"; return nullptr; } return particles_m[partIndex][attrIndex]; }

	///Forward decs for cpp file, or pure virtuals
	//Numerical tools
	virtual void generateNormallyDistributedValues(int numberOfNormalAttributes, double* means, double* sigmas);//tested
	virtual double calculateMeanOfParticleAttribute(int particleIndex, int attributeIndex, bool absValue=false);//tested
	virtual double calculateStdDevOfParticleAttribute(int particleIndex, int attributeIndex);//tested

	//Array tools
	virtual void saveParticleAttributeToDisk(int particleIndex, int attributeIndex, const char* foldername, const char* name);//what is written to disk needs to be tested - make sure right attributes are in the right file
	virtual void loadFileIntoParticleAttribute(int particleIndex, int attributeIndex, const char* foldername, const char* name);
	
	//Field tools
	virtual double calculateBFieldAtZandTime(double z, double time) = 0;
	virtual double calculateEFieldAtZandTime(double z, double time) = 0;

	//Simulation management functions
	virtual void initializeSimulation() = 0;
	virtual void copyDataToGPU() = 0;
	virtual void iterateSimulation(int numberOfIterations) = 0;
	virtual void copyDataToHost() = 0;
	virtual void freeGPUMemory() = 0;
	virtual void prepareResults() = 0;

	//Access functions
	virtual double getSimMin() { return simMin_m; }
	virtual double getSimMax() { return simMax_m; }

	void createSatellite(double altitude, bool upwardFacing, double** GPUdataPointers, std::string name);
	virtual int getSatelliteCount() { return satellites_m.size(); }
	virtual timeStruct* createTimeStruct(std::string label);
	virtual void printTimeNowFromTimeStruct(timeStruct* tS, std::string label);
	virtual void printTimeDiffBtwTwoTimeStructs(timeStruct* startTS, timeStruct* endTS);
};//end class
#endif //end header guard