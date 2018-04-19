#ifndef SIMULATIONAPI_H
#define SIMULATIONAPI_H

#include <iostream>
#include "dllexport.h"
#include "SimulationClass\Simulation.h"
#include "utils\write.h"
#include "ErrorHandling\simExceptionMacros.h"

///One liner functions
DLLEXPORT double      getSimulationTimeAPI(Simulation* simulation);
DLLEXPORT double      getDtAPI(Simulation* simulation);
DLLEXPORT double      getSimMinAPI(Simulation* simulation);
DLLEXPORT double      getSimMaxAPI(Simulation* simulation);
DLLEXPORT int         getNumberOfParticleTypesAPI(Simulation* simulation);
DLLEXPORT int         getNumberOfParticlesAPI(Simulation* simulation, int partInd);
DLLEXPORT int         getNumberOfAttributesAPI(Simulation* simulation, int partInd);
DLLEXPORT const char* getParticleNameAPI(Simulation* simulation, int partInd);
DLLEXPORT const char* getSatelliteNameAPI(Simulation* simulation, int satInd);
DLLEXPORT LogFile*    getLogFilePointerAPI(Simulation* simulation);

//Pointer one liners
DLLEXPORT const double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData);

//Field tools
DLLEXPORT double getBFieldAtSAPI(Simulation* simulation, double s, double time);
DLLEXPORT double getEFieldAtSAPI(Simulation* simulation, double s, double time);

//Simulation Management Function Wrappers
DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir);
DLLEXPORT void initializeSimulationAPI(Simulation* simulation);
DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts);
DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation);
DLLEXPORT void terminateSimulationAPI(Simulation* simulation);
DLLEXPORT Simulation* loadCompletedSimDataAPI(const char* fileDir);
DLLEXPORT void setupNormalSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir);
DLLEXPORT void runNormalSimulationAPI(Simulation* sim, int iterations, int printEvery);

//build API functions for resetSimulation and saveDataToDisk

//Fields management
DLLEXPORT void setBFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString); //switch to comma delimited string of variables
DLLEXPORT void addEFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString);

//Particle functions
DLLEXPORT void createParticleTypeAPI(Simulation* simulation, const char* name, double mass, double charge, long numParts, const char* loadFileDir = "");

//Satellite functions
DLLEXPORT void    createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name);
DLLEXPORT int     getNumberOfSatellitesAPI(Simulation* simulation);
DLLEXPORT const double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int msmtInd, int attributeInd);

//CSV functions
DLLEXPORT void writeCommonCSVAPI(Simulation* simulation);

#endif//end if for header guard