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

DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfParticleTypes()); }

DLLEXPORT int getNumberOfParticlesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return static_cast<int>(simulation->getNumberOfParticles(partInd))); }

DLLEXPORT int getNumberOfAttributesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return static_cast<int>(simulation->getNumberOfAttributes(partInd))); }

DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation) { //do I even need this?  maybe change way results are prepared
	return simulation->areResultsPrepared(); }

DLLEXPORT LogFile* getLogFilePointerAPI(Simulation* simulation) {
	return simulation->getLogFilePointer(); }


//Pointer one liners
DLLEXPORT double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	SIM_API_EXCEP_CHECK(return simulation->getPointerToParticleAttributeArray(partIndex, attrIndex, originalData)); }


//Field tools
DLLEXPORT double calculateBFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateBFieldAtZandTime(z, time); }

DLLEXPORT double calculateEFieldAtZandTimeAPI(Simulation* simulation, double z, double time) {
	return simulation->calculateEFieldAtZandTime(z, time); }


//Mu<->VPerp Functions - Should I be exposing this functionality??  Maybe
DLLEXPORT void convertParticleVPerpToMuAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(simulation->convertVPerpToMu(partInd)); }

DLLEXPORT void convertParticleMuToVPerpAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(simulation->convertMuToVPerp(partInd)); }


//Simulation Management Function Wrappers
DLLEXPORT void initializeSimulationAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->initializeSimulation()); }

DLLEXPORT void copyDataToGPUAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->copyDataToGPU()); }

DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(simulation->iterateSimulation(numberOfIterations, itersBtwCouts)); }

DLLEXPORT void copyDataToHostAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->copyDataToHost()); }

DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->freeGPUMemory()); }

DLLEXPORT void prepareResultsAPI(Simulation* simulation, bool normalizeToRe) {
	SIM_API_EXCEP_CHECK(simulation->prepareResults(normalizeToRe)); }


//Satellite functions
DLLEXPORT void createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name) {
	SIM_API_EXCEP_CHECK(simulation->createTempSat(particleInd, altitude, upwardFacing, name)); }

DLLEXPORT int  getNumberOfSatellitesAPI(Simulation* simulation) {
	return static_cast<int>(simulation->getNumberOfSatellites()); }

DLLEXPORT int  getSatNumOfDetectedPartsAPI(Simulation* simulation, int satIndex) {
	SIM_API_EXCEP_CHECK(return static_cast<int>(simulation->getSatelliteNumberOfDetectedParticles(satIndex))); }

DLLEXPORT double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int attributeInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteDataPointers(satelliteInd, attributeInd)); }

DLLEXPORT void writeSatelliteDataToCSVAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->writeSatelliteDataToCSV()); }

DLLEXPORT void loadCompletedSimDataAPI(Simulation* simulation, const char* fileDir, const char* partNames, const char* attrNames, const char* satNames, int numParts) {
	SIM_API_EXCEP_CHECK(simulation->loadCompletedSimData(fileDir, constCharToStrVec(partNames), constCharToStrVec(attrNames), constCharToStrVec(satNames), numParts)); }


//Particle functions
DLLEXPORT void createParticleTypeAPI(Simulation* simulation, const char* name, const char* attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, const char* loadFileDir) {
	SIM_API_EXCEP_CHECK(simulation->createParticleType(name, constCharToStrVec(attrNames), mass, charge, numParts, posDims, velDims, normFactor, loadFileDir)); }


//Simulation creation and deletion
DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, double ionT, double magT, const char* rootdir) {
	SIM_API_EXCEP_CHECK(return new Simulation(dt, simMin, simMax, ionT, magT, rootdir)); }

DLLEXPORT void runNormalSimulationAPI(Simulation* sim, int iterations, int printEvery, const char* loadFileDir)
{
	SIM_API_EXCEP_CHECK( \
	double simMin{ sim->getSimMin() };
	double simMax{ sim->getSimMax() };

	sim->setBFieldModel("DipoleB", { 72.0 });

	sim->createParticleType("elec", { "vpara", "vperp", "s" }, MASS_ELECTRON, -1 * CHARGE_ELEM, 115200, 1, 2, RADIUS_EARTH, loadFileDir);
	sim->createParticleType("ions", { "vpara", "vperp", "s" }, MASS_PROTON,    1 * CHARGE_ELEM, 115200, 1, 2, RADIUS_EARTH, loadFileDir);

	sim->createTempSat(0, simMin * 0.999, true,  "btmElec");
	sim->createTempSat(1, simMin * 0.999, true,  "btmIons");
	sim->createTempSat(0, simMax * 1.001, false, "topElec");
	sim->createTempSat(1, simMax * 1.001, false, "topIons");

	sim->initializeSimulation();
	sim->copyDataToGPU();
	sim->iterateSimulation(iterations, printEvery);
	sim->copyDataToHost();
	sim->freeGPUMemory();
	sim->prepareResults(true);
	); /* SIM_API_EXCEP_CHECK() */
}

DLLEXPORT void terminateSimulationAPI(Simulation* simulation) {
	delete simulation; }

DLLEXPORT void setBFieldModelAPI(Simulation* sim, const char* modelName, double arg1) {
	SIM_API_EXCEP_CHECK(sim->setBFieldModel(modelName, { arg1 })); }

//#undef DLLFILE //uncomment for making an exe file for profiling
#ifndef DLLFILE
int main()
{
	SIM_API_EXCEP_CHECK( \
	Simulation* sim{ createSimulationAPI(0.01, 2030837.49610366, 19881647.2473464, 2.5, 1000.0, "./../") };

	runNormalSimulationAPI(sim, 250, 50, "./../_in/data");

	terminateSimulationAPI(sim);
	); /* SIM_API_EXCEP_CHECK() */

	return 0;
}
#endif