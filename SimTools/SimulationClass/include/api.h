#ifndef SIMULATIONAPI_H
#define SIMULATIONAPI_H

#include <chrono>
#include <iostream>
#include "include\Simulation.h"

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

///One liner functions
DLLEXPORT double getSimulationTimeWrapper(Simulation* simulation);
DLLEXPORT double getDtWrapper(Simulation* simulation);
DLLEXPORT void incrementSimulationTimeByDtWrapper(Simulation* simulation);
DLLEXPORT int getNumberOfParticleTypesWrapper(Simulation* simulation);
DLLEXPORT int getNumberOfParticlesPerTypeWrapper(Simulation* simulation);
DLLEXPORT int getNumberOfAttributesTrackedWrapper(Simulation* simulation);
DLLEXPORT bool areResultsPreparedWrapper(Simulation* simulation);

//Pointer one liners
DLLEXPORT double*** getPointerTo3DParticleArrayWrapper(Simulation* simulation);
DLLEXPORT double** getPointerToSingleParticleTypeArrayWrapper(Simulation* simulation, int index);
DLLEXPORT double* getPointerToSerializedParticleArrayWrapper(Simulation* simulation);
DLLEXPORT bool* getPointerToParticlesInSimArrayWrapper(Simulation* simulation, int index);
DLLEXPORT double* getPointerToSingleParticleAttributeArrayWrapper(Simulation* simulation, int partIndex, int attrIndex);

//Numerical tools
DLLEXPORT void generateNormallyDistributedValuesWrapper(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas);
DLLEXPORT double calculateMeanOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex, bool absValue);
DLLEXPORT double calculateStdDevOfParticleAttributeWrapper(Simulation* simulation, int particleIndex, int attributeIndex);
DLLEXPORT double* getPointerToSingleParticleAttributeArrayWrapper(Simulation* simulation, int partIndex, int attrIndex);

//Array tools
DLLEXPORT void serializeParticleArrayWrapper(Simulation* simulation);
DLLEXPORT double calculateBFieldAtZandTimeWrapper(Simulation* simulation, double z, double time);
DLLEXPORT double calculateEFieldAtZandTimeWrapper(Simulation* simulation, double z, double time);

//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationWrapper(Simulation* simulation);
DLLEXPORT void copyDataToGPUWrapper(Simulation* simulation);
DLLEXPORT void iterateSimulationWrapper(Simulation* simulation, int numberOfIterations);
DLLEXPORT void copyDataToHostWrapper(Simulation* simulation);
DLLEXPORT void freeGPUMemoryWrapper(Simulation* simulation);
DLLEXPORT void prepareResultsWrapper(Simulation* simulation);

#endif//end if for header guard