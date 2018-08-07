#ifndef SIMULATIONAPI_H
#define SIMULATIONAPI_H

#include "dlldefines.h"
#include "Simulation/Simulation.h"

///One liner functions
DLLEXP_EXTC double      getSimulationTimeAPI(Simulation* simulation);
DLLEXP_EXTC double      getDtAPI(Simulation* simulation);
DLLEXP_EXTC double      getSimMinAPI(Simulation* simulation);
DLLEXP_EXTC double      getSimMaxAPI(Simulation* simulation);
DLLEXP_EXTC int         getNumberOfParticleTypesAPI(Simulation* simulation);
DLLEXP_EXTC int         getNumberOfParticlesAPI(Simulation* simulation, int partInd);
DLLEXP_EXTC int         getNumberOfAttributesAPI(Simulation* simulation, int partInd);
DLLEXP_EXTC const char* getParticleNameAPI(Simulation* simulation, int partInd);
DLLEXP_EXTC const char* getSatelliteNameAPI(Simulation* simulation, int satInd);
DLLEXP_EXTC LogFile*    getLogFilePointerAPI(Simulation* simulation);

//Pointer one liners
DLLEXP_EXTC const double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData);

//Field tools
DLLEXP_EXTC double getBFieldAtSAPI(Simulation* simulation, double s, double time);
DLLEXP_EXTC double getEFieldAtSAPI(Simulation* simulation, double s, double time);

//Simulation Management Function Wrappers
DLLEXP_EXTC Simulation* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir);
DLLEXP_EXTC void initializeSimulationAPI(Simulation* simulation);
DLLEXP_EXTC void __iterateSimCPUAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts);
DLLEXP_EXTC void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts);
DLLEXP_EXTC void freeGPUMemoryAPI(Simulation* simulation);
DLLEXP_EXTC void terminateSimulationAPI(Simulation* simulation);
DLLEXP_EXTC Simulation* loadCompletedSimDataAPI(const char* fileDir);
DLLEXP_EXTC void setupExampleSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir);
DLLEXP_EXTC void runExampleSimulationAPI(Simulation* sim, int iterations, int printEvery);
DLLEXP_EXTC void runSingleElectronAPI(Simulation* sim, double vpara, double vperp, double s, double t_inc, int iterations, int printEvery);

//build API functions for resetSimulation and saveDataToDisk

//Fields management
DLLEXP_EXTC void setBFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString); //switch to comma delimited string of variables
DLLEXP_EXTC void addEFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString);

//Particle functions
DLLEXP_EXTC void createParticleTypeAPI(Simulation* simulation, const char* name, double mass, double charge, long numParts, const char* loadFileDir = "");

//Satellite functions
DLLEXP_EXTC void    createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name);
DLLEXP_EXTC int     getNumberOfSatellitesAPI(Simulation* simulation);
DLLEXP_EXTC const double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int msmtInd, int attributeInd);

//CSV functions
DLLEXP_EXTC void writeCommonCSVAPI(Simulation* simulation);

#endif//end if for header guard