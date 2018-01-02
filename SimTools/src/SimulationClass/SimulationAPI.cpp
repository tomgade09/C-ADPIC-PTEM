#include "SimulationClass\SimulationAPI.h"

///One liner functions
DLLEXPORT double getSimulationTimeAPI(Simulation* simulation) {
	return simulation->getTime(); }

DLLEXPORT double getDtAPI(Simulation* simulation) {
	return simulation->getdt(); }

DLLEXPORT void incrementSimulationTimeByDtAPI(Simulation* simulation) {
	simulation->incTime(); }

DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfParticleTypes()); }

DLLEXPORT int getNumberOfParticlesAPI(Simulation* simulation, int partInd) {
	return static_cast<int>(simulation->getNumberOfParticles(partInd)); }

DLLEXPORT int getNumberOfAttributesAPI(Simulation* simulation, int partInd) {
	return static_cast<int>(simulation->getNumberOfAttributes(partInd)); }

DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation) {
	return simulation->areResultsPrepared(); }

DLLEXPORT bool getNormalizedAPI(Simulation* simulation) {
	return simulation->getNormalized(); }

DLLEXPORT LogFile* getLogFilePointerAPI(Simulation* simulation) {
	return simulation->getLogFilePointer(); }

DLLEXPORT double getSimMinAPI(Simulation* simulation) {
	return simulation->getSimMin(); }

DLLEXPORT double getSimMaxAPI(Simulation* simulation) {
	return simulation->getSimMax(); }

//Pointer one liners
DLLEXPORT double* getPointerToSingleParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	return simulation->getPointerToSingleParticleAttributeArray(partIndex, attrIndex, originalData); }

//Numerical tools
DLLEXPORT double  calculateMeanOfParticleAttributeAPI(double* data, int length, bool absValue) {
	return calculateMeanOfParticleAttribute(data, length, absValue); }

DLLEXPORT double  calculateStdDevOfParticleAttributeAPI(double* data, int length) {
	return calculateStdDevOfParticleAttribute(data, length); }

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
	return static_cast<int>(simulation->getNumberOfSatellites()); }

DLLEXPORT int  getNumberOfSatelliteMsmtsAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfSatelliteMsmts()); }

DLLEXPORT double* getSatelliteDataPointersAPI(Simulation* simulation, int measurementInd, int satelliteInd, int attributeInd) {
	return simulation->getSatelliteDataPointers(measurementInd, satelliteInd, attributeInd); }