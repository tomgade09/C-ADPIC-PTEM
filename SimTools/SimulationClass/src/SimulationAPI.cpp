#include "include\SimulationAPI.h"

///One liner functions
DLLEXPORT double getSimulationTimeWrapper(Simulation* simulation)
{
	return simulation->getTime();
}

DLLEXPORT double getDtWrapper(Simulation* simulation)
{
	return simulation->getdt();
}

DLLEXPORT void incrementSimulationTimeByDtWrapper(Simulation* simulation)
{
	simulation->incTime();
}

DLLEXPORT int getNumberOfParticleTypesWrapper(Simulation* simulation)
{
	return simulation->getNumberOfParticleTypes();
}

DLLEXPORT int getNumberOfParticlesPerTypeWrapper(Simulation* simulation)
{
	return simulation->getNumberOfParticlesPerType();
}

DLLEXPORT int getNumberOfAttributesTrackedWrapper(Simulation* simulation)
{
	return simulation->getNumberOfAttributesTracked();
}

DLLEXPORT bool areResultsPreparedWrapper(Simulation* simulation)
{
	return simulation->areResultsPrepared();
}

DLLEXPORT bool getNormalizedWrapper(Simulation* simulation)
{
	return simulation->getNormalized();
}

DLLEXPORT bool getReplenishWrapper(Simulation* simulation)
{
	return simulation->getReplenish();
}


//Pointer one liners
DLLEXPORT double*** getPointerTo3DParticleArrayWrapper(Simulation* simulation)
{
	return simulation->getPointerTo3DParticleArray();
}

DLLEXPORT double** getPointerToSingleParticleTypeArrayWrapper(Simulation* simulation, int index)
{
	return simulation->getPointerToSingleParticleTypeArray(index);
}

DLLEXPORT double* getPointerToSerializedParticleArrayWrapper(Simulation* simulation)
{
	return simulation->getPointerToSerializedParticleArray();
}

DLLEXPORT bool* getPointerToParticlesInSimArrayWrapper(Simulation* simulation, int index)
{
	return simulation->getPointerToParticlesInSimArray(index);
}

DLLEXPORT double* getPointerToSingleParticleAttributeArrayWrapper(Simulation* simulation, int partIndex, int attrIndex)
{
	return simulation->getPointerToSingleParticleAttributeArray(partIndex, attrIndex);
}

//Numerical tools
DLLEXPORT void generateNormallyDistributedValuesWrapper(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas)
{
	simulation->generateNormallyDistributedValues(numberOfNormalAttributes, means, sigmas);
}

DLLEXPORT double calculateMeanOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex, bool absValue)
{
	return simulation->calculateMeanOfParticleAttribute(particleIndex, attributeIndex, absValue);
}

DLLEXPORT double calculateStdDevOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex)
{
	return simulation->calculateStdDevOfParticleAttribute(particleIndex, attributeIndex);
}

//Array tools
DLLEXPORT void serializeParticleArrayWrapper(Simulation* simulation)
{
	simulation->serializeParticleArray();
}

DLLEXPORT double calculateBFieldAtZandTimeWrapper(Simulation* simulation, double z, double time)
{
	return simulation->calculateBFieldAtZandTime(z, time);
}

DLLEXPORT double calculateEFieldAtZandTimeWrapper(Simulation* simulation, double z, double time)
{
	return simulation->calculateEFieldAtZandTime(z, time);
}

//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationWrapper(Simulation* simulation)
{
	simulation->initializeSimulation();
}

DLLEXPORT void copyDataToGPUWrapper(Simulation* simulation)
{
	simulation->copyDataToGPU();
}

DLLEXPORT void iterateSimulationWrapper(Simulation* simulation, int numberOfIterations)
{
	std::chrono::steady_clock::time_point cudaBegin, cudaEnd;
	cudaBegin = std::chrono::steady_clock::now();
	simulation->iterateSimulation(numberOfIterations);
	cudaEnd = std::chrono::steady_clock::now();
	std::cout << "Parallel Execution Time (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>(cudaEnd - cudaBegin).count() << std::endl;
}

DLLEXPORT void copyDataToHostWrapper(Simulation* simulation)
{
	simulation->copyDataToHost();
}

DLLEXPORT void freeGPUMemoryWrapper(Simulation* simulation)
{
	simulation->freeGPUMemory();
}

DLLEXPORT void prepareResultsWrapper(Simulation* simulation)
{
	simulation->prepareResults();
}

//Satellite functions
DLLEXPORT void createSatelliteWrapper(Simulation* simulation, double altitude, bool upwardFacing, int particleIndex, const char* name)
{
	simulation->createSatellite(altitude, upwardFacing, name);
}

DLLEXPORT int  getNumberOfSatellitesWrapper(Simulation* simulation)
{
	return simulation->getNumberOfSatellites();
}

DLLEXPORT int  getNumberOfSatelliteMsmtsWrapper(Simulation* simulation)
{
	return simulation->getNumberOfSatelliteMsmts();
}

DLLEXPORT double* getSatelliteDataPointersWrapper(Simulation* simulation, int measurementInd, int satelliteInd, int attributeInd)
{
	return simulation->getSatelliteDataPointers(measurementInd, satelliteInd, attributeInd);
}