#ifndef SIMULATIONAPI_H
#define SIMULATIONAPI_H

#include "include\Simulation.h"

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

///One liner functions
DLLEXPORT double getSimulationTimeWrapper(Simulation* simulation);
DLLEXPORT void incrementSimulationTimeByDtWrapper(Simulation* simulation);
DLLEXPORT void resetParticlesEscapedCountWrapper(Simulation* simulation);

//Pointer one liners
DLLEXPORT double*** getPointerTo3DParticleArrayWrapper(Simulation* simulation);
DLLEXPORT double** getPointerToSingleParticleTypeArrayWrapper(Simulation* simulation, int index);
DLLEXPORT double* getPointerToSerializedParticleArrayWrapper(Simulation* simulation);

//Numerical tools
DLLEXPORT void generateNormallyDistributedValues(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas);
DLLEXPORT double calculateMeanOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex, bool absValue);
DLLEXPORT double calculateStdDevOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex);

//Array tools
DLLEXPORT void serializeParticleArrayWrapper(Simulation* simulation);
DLLEXPORT double calculateBFieldAtZandTimeWrapper(Simulation* simulation, double z, double time);
DLLEXPORT double calculateEFieldAtZandTimeWrapper(Simulation* simulation, double z, double time);

//Simulation Management Function Wrappers
DLLEXPORT void initializeWrapper(Simulation* simulation);
DLLEXPORT void copyDataToGPUWrapper(Simulation* simulation);
DLLEXPORT void iterateSimulationWrapper(Simulation* simulation, int numberOfIterations);
DLLEXPORT void copyDataToHostWrapper(Simulation* simulation);
DLLEXPORT void terminateSimulationWrapper(Simulation* simulation);
DLLEXPORT double* returnResultsWrapper(Simulation* simulation);

#endif//end if for header guard