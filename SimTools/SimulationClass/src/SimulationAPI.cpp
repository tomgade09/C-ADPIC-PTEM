#include "include\SimulationAPI.h"

///One liner functions
DLLEXPORT double getSimulationTimeAPI(Simulation* simulation) {
	return simulation->getTime(); }

DLLEXPORT double getDtAPI(Simulation* simulation) {
	return simulation->getdt(); }

DLLEXPORT void incrementSimulationTimeByDtAPI(Simulation* simulation) {
	simulation->incTime(); }

DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return simulation->getNumberOfParticleTypes(); }

DLLEXPORT int getNumberOfParticlesPerTypeAPI(Simulation* simulation) {
	return simulation->getNumberOfParticlesPerType(); }

DLLEXPORT int getNumberOfAttributesTrackedAPI(Simulation* simulation) {
	return simulation->getNumberOfAttributesTracked(); }

DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation) {
	return simulation->areResultsPrepared(); }

DLLEXPORT bool getNormalizedAPI(Simulation* simulation) {
	return simulation->getNormalized(); }

DLLEXPORT double getSimMinAPI(Simulation* simulation) {
	return simulation->getSimMin(); }

DLLEXPORT double getSimMaxAPI(Simulation* simulation) {
	return simulation->getSimMax(); }

//Pointer one liners
DLLEXPORT double*** getPointerTo3DParticleArrayAPI(Simulation* simulation) {
	return simulation->getPointerTo3DParticleArray(); }

DLLEXPORT double** getPointerToSingleParticleTypeArrayAPI(Simulation* simulation, int index) {
	return simulation->getPointerToSingleParticleTypeArray(index); }

DLLEXPORT double* getPointerToSingleParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex) {
	return simulation->getPointerToSingleParticleAttributeArray(partIndex, attrIndex); }

//Numerical tools
DLLEXPORT void generateNormallyDistributedValuesAPI(Simulation* simulation, int numberOfNormalAttributes, double* means, double* sigmas) {
	simulation->generateNormallyDistributedValues(numberOfNormalAttributes, means, sigmas); }

DLLEXPORT double calculateMeanOfParticleAttributeAPI(Simulation* simulation, int particleIndex, int attributeIndex, bool absValue) {
	return simulation->calculateMeanOfParticleAttribute(particleIndex, attributeIndex, absValue); }

DLLEXPORT double calculateStdDevOfParticleAttributeAPI(Simulation* simulation, int particleIndex, int attributeIndex) {
	return simulation->calculateStdDevOfParticleAttribute(particleIndex, attributeIndex); }

//Field tools
DLLEXPORT double calculateBFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateBFieldAtZandTime(z, time); }

DLLEXPORT double calculateEFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateEFieldAtZandTime(z, time); }

//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationAPI(Simulation* simulation) {
	simulation->initializeSimulation(); }

DLLEXPORT void copyDataToGPUAPI(Simulation* simulation) {
	simulation->copyDataToGPU(); }

DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations) {
	simulation->iterateSimulation(numberOfIterations); }

DLLEXPORT void copyDataToHostAPI(Simulation* simulation) {
	simulation->copyDataToHost(); }

DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation) {
	simulation->freeGPUMemory(); }

DLLEXPORT void prepareResultsAPI(Simulation* simulation) {
	simulation->prepareResults(); }

//Satellite functions
DLLEXPORT int  getNumberOfSatellitesAPI(Simulation* simulation) {
	return simulation->getNumberOfSatellites(); }

DLLEXPORT int  getNumberOfSatelliteMsmtsAPI(Simulation* simulation) {
	return simulation->getNumberOfSatelliteMsmts(); }

DLLEXPORT double* getSatelliteDataPointersAPI(Simulation* simulation, int measurementInd, int satelliteInd, int attributeInd) {
	return simulation->getSatelliteDataPointers(measurementInd, satelliteInd, attributeInd); }