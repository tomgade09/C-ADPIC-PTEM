#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <random>
#include <cmath>
#include <iostream>
#include <vector>

class Simulation
{
protected:
	//Simulation Characteristics
	double    simTime_m{ 0 };
	double    dt_m;
	int	      numberOfParticleTypes_m;
	int       numberOfParticlesPerType_m;
	int       numberOfAttributesTracked_m;

	//Data Array Pointers
	double*** particles_m{ nullptr };
	double*   particlesSerialized_m{ nullptr };
	bool**    particlesInSim_m{ nullptr };
	int**     particlesEscaped_m{ nullptr };
	double*   fieldsDataForOutput{ nullptr };
	std::vector<void*> otherMemoryPointers;
	
	//GPU Memory Pointers
	double**  gpuDblMemoryPointers_m { nullptr };
	bool**    gpuBoolMemoryPointers_m{ nullptr };
	int**     gpuIntMemoryPointers_m { nullptr };
	void**    gpuOtherMemoryPointers_m{ nullptr };

	//Calculated Quantities
	long totalElecEscaped_m{ 0 };
	long totalIonsEscaped_m{ 0 };
	
public:
	Simulation(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt):
		numberOfParticleTypes_m{ numberOfParticleTypes }, numberOfParticlesPerType_m{ numberOfParticlesPerType },
		numberOfAttributesTracked_m{ numberOfAttributesTracked }, dt_m{ dt }
	{
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
				for (int kkk = 0; kkk < numberOfParticlesPerType_m; kkk++)
				{
					particlesInSim_m[iii][kkk] = true;
					particlesEscaped_m[iii][kkk] = 0;
					particles_m[iii][jjj][kkk] = 0.0;
				}//end loop kkk
			}//end loop jjj
		}//end loop iii

		//Create Arrays for GPU Memory Pointers
		gpuDblMemoryPointers_m  = new double*[numberOfParticleTypes_m * numberOfAttributesTracked_m]; //holds pointers to GPU memory for particle attributes
		gpuBoolMemoryPointers_m = new bool*  [numberOfParticleTypes_m]; //holds pointers to GPU memory for bool flag of whether particle is in sim or not
		gpuIntMemoryPointers_m  = new int*   [numberOfParticleTypes_m]; //holds pointers to GPU memory for count of escaped particles
	}

	virtual ~Simulation() {};//need to write this - or leaked memory will result.
	//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double getTime() { return simTime_m; }//tested
	void incTime() { simTime_m += dt_m; }
	void resetParticlesEscapedCount() { totalElecEscaped_m = 0; totalIonsEscaped_m = 0; return; }
	int getNumberOfParticleTypes() { return numberOfParticleTypes_m; }//tested
	int getNumberOfParticlesPerType() { return numberOfParticlesPerType_m; }//tested
	int getNumberOfAttributesTracked() { return numberOfAttributesTracked_m; }//tested

	//Pointer one liners
	double*** getPointerTo3DParticleArray() { return particles_m; }//tested
	double**  getPointerToSingleParticleTypeArray(int index) { if (index >= numberOfParticleTypes_m) { return nullptr; } return particles_m[index]; }
	double*   getPointerToSerializedParticleArray() { if (particlesSerialized_m == nullptr)//tested
		{ std::cout << "Array not serialized yet.  Run Simulation::serializeParticleArray.\n"; } return particlesSerialized_m; }

	///Forward decs for cpp file, or pure virtuals
	//Numerical tools
	virtual void generateNormallyDistributedValues(int numberOfNormalAttributes, double* means, double* sigmas);//tested
	virtual double calculateMeanOfParticleAttribute(int particleIndex, int attributeIndex, bool absValue=false);//tested
	virtual double calculateStdDevOfParticleAttribute(int particleIndex, int attributeIndex);//tested

	//Array tools
	virtual void saveParticleAttributeToDisk(int particleIndex, int attributeIndex, const char* foldername, const char* name);//what is written to disk needs to be tested - make sure right attributes are in the right file
	virtual void serializeParticleArray();//tested
	virtual double calculateBFieldAtZandTime(double z, double time) = 0;
	virtual double calculateEFieldAtZandTime(double z, double time) = 0;

	//Simulation management functions
	virtual void initializeSimulation() = 0;
	virtual void copyDataToGPU() = 0;
	virtual void iterateSimulation(int numberOfIterations) = 0;
	virtual void copyDataToHost() = 0;
	virtual void terminateSimulation() = 0;
	virtual double* returnResults() = 0;
	
};//end class
#endif //end header guard