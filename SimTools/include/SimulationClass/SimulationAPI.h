#ifndef SIMULATIONAPI_H
#define SIMULATIONAPI_H

#include <iostream>
#include "SimulationClass\Simulation.h"

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

///One liner functions
DLLEXPORT double getSimulationTimeAPI(Simulation* simulation);
DLLEXPORT double getDtAPI(Simulation* simulation);
DLLEXPORT void incrementSimulationTimeByDtAPI(Simulation* simulation);
DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation);
DLLEXPORT int getNumberOfParticlesPerTypeAPI(Simulation* simulation);
DLLEXPORT int getNumberOfAttributesTrackedAPI(Simulation* simulation);
DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation);
DLLEXPORT bool getNormalizedAPI(Simulation* simulation);
DLLEXPORT double getSimMinAPI(Simulation* simulation);
DLLEXPORT double getSimMaxAPI(Simulation* simulation);

//Pointer one liners
DLLEXPORT double*** getPointerTo3DParticleArrayAPI(Simulation* simulation);
DLLEXPORT double** getPointerToSingleParticleTypeArrayAPI(Simulation* simulation, int index);
DLLEXPORT double* getPointerToSingleParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex);

//Numerical tools
DLLEXPORT void generateNormallyDistributedValuesAPI(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas);
DLLEXPORT double calculateMeanOfParticleAttributeAPI(Simulation* simulation, int particleIndex, int attributeIndex, bool absValue);
DLLEXPORT double calculateStdDevOfParticleAttributeAPI(Simulation* simulation, int particleIndex, int attributeIndex);
DLLEXPORT double* getPointerToSingleParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex);

//Field tools
DLLEXPORT double calculateBFieldAtZandTimeAPI(Simulation* simulation, double z, double time);
DLLEXPORT double calculateEFieldAtZandTimeAPI(Simulation* simulation, double z, double time);

//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationAPI(Simulation* simulation);
DLLEXPORT void copyDataToGPUAPI(Simulation* simulation);
DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations);
DLLEXPORT void copyDataToHostAPI(Simulation* simulation);
DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation);
DLLEXPORT void prepareResultsAPI(Simulation* simulation);

DLLEXPORT int  getNumberOfSatellitesAPI(Simulation* simulation);
DLLEXPORT int  getNumberOfSatelliteMsmtsAPI(Simulation* simulation);
DLLEXPORT double* getSatelliteDataPointersAPI(Simulation* simulation, int measurementInd, int satelliteInd, int attributeInd);

#endif//end if for header guard