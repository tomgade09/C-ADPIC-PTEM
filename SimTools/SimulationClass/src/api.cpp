#include "include\api.h"

///One liner functions
DLLEXPORT double getSimulationTimeWrapper(Simulation* simulation)
{
	return simulation->getTime();
}

DLLEXPORT void incrementSimulationTimeByDtWrapper(Simulation* simulation)
{
	simulation->incTime();
}

DLLEXPORT void resetParticlesEscapedCountWrapper(Simulation* simulation)
{
	simulation->resetParticlesEscapedCount();
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

//Numerical tools
DLLEXPORT void generateNormallyDistributedValues(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas)
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
DLLEXPORT void initializeWrapper(Simulation* simulation)
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

DLLEXPORT void terminateSimulationWrapper(Simulation* simulation)
{
	simulation->terminateSimulation();
}

DLLEXPORT double* returnResultsWrapper(Simulation* simulation)
{
	return simulation->returnResults();
}