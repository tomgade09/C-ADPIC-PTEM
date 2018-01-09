#include "SimulationClass\SimulationAPI.h"

///One liner functions
DLLEXPORT double getSimulationTimeAPI(Simulation* simulation) {
	return simulation->getTime(); }

DLLEXPORT double getDtAPI(Simulation* simulation) {
	return simulation->getdt(); }

DLLEXPORT double getSimMinAPI(Simulation* simulation) {
	return simulation->getSimMin(); }

DLLEXPORT double getSimMaxAPI(Simulation* simulation) {
	return simulation->getSimMax(); }

DLLEXPORT void incrementSimulationTimeByDtAPI(Simulation* simulation) {
	simulation->incTime(); }

DLLEXPORT void setQSPSAPI(Simulation* simulation, double constE) {
	simulation->setQSPS(constE); }

DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfParticleTypes()); }

DLLEXPORT int getNumberOfParticlesAPI(Simulation* simulation, int partInd) {
	return static_cast<int>(simulation->getNumberOfParticles(partInd)); }

DLLEXPORT int getNumberOfAttributesAPI(Simulation* simulation, int partInd) {
	return static_cast<int>(simulation->getNumberOfAttributes(partInd)); }

DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation) {
	return simulation->areResultsPrepared(); }

DLLEXPORT LogFile* getLogFilePointerAPI(Simulation* simulation) {
	return simulation->getLogFilePointer(); }


//Pointer one liners
DLLEXPORT double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	return simulation->getPointerToParticleAttributeArray(partIndex, attrIndex, originalData); }


//Field tools
DLLEXPORT double calculateBFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateBFieldAtZandTime(z, time); }

DLLEXPORT double calculateEFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateEFieldAtZandTime(z, time); }


//Mu<->VPerp Functions
DLLEXPORT void convertParticleVPerpToMuAPI(Simulation* simulation, int partInd) {
	simulation->convertVPerpToMu(partInd); }

DLLEXPORT void convertParticleMuToVPerpAPI(Simulation* simulation, int partInd) {
	simulation->convertMuToVPerp(partInd); }


//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationAPI(Simulation* simulation) {
	simulation->initializeSimulation(); }

DLLEXPORT void copyDataToGPUAPI(Simulation* simulation) {
	simulation->copyDataToGPU(); }

DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts) {
	simulation->iterateSimulation(numberOfIterations, itersBtwCouts); }

DLLEXPORT void copyDataToHostAPI(Simulation* simulation) {
	simulation->copyDataToHost(); }

DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation) {
	simulation->freeGPUMemory(); }

DLLEXPORT void prepareResultsAPI(Simulation* simulation, bool normalizeToRe) {
	simulation->prepareResults(normalizeToRe); }


//Satellite functions
DLLEXPORT void createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name) {
	simulation->createTempSat(particleInd, altitude, upwardFacing, name); }

DLLEXPORT int  getNumberOfSatellitesAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfSatellites()); }

DLLEXPORT int  getSatNumOfDetectedPartsAPI(Simulation* simulation, int satIndex) {
	return static_cast<int>(simulation->getSatelliteNumberOfDetectedParticles(satIndex)); }

DLLEXPORT double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int attributeInd) {
	return simulation->getSatelliteDataPointers(satelliteInd, attributeInd); }

DLLEXPORT void writeSatelliteDataToCSVAPI(Simulation* simulation) {
	simulation->writeSatelliteDataToCSV(); }

DLLEXPORT void loadCompletedSimDataAPI(Simulation* simulation, const char* fileDir, const char* partNames, const char* attrNames, const char* satNames, int numParts) {
	simulation->loadCompletedSimData(fileDir, constCharToStrVec(partNames), constCharToStrVec(attrNames), constCharToStrVec(satNames), numParts); }

//Particle functions
DLLEXPORT void createParticleTypeAPI(Simulation* simulation, const char* name, const char* attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, const char* loadFileDir) {
	simulation->createParticleType(name, constCharToStrVec(attrNames), mass, charge, numParts, posDims, velDims, normFactor, loadFileDir); }