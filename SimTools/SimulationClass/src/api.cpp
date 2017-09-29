#include "include\api.h"

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
	simulation->iterateSimulation(numberOfIterations);
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
