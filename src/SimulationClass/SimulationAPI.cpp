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


//Simulation creation and deletion
DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, double ionT, double magT, const char* rootdir) {
	return new Simulation(dt, simMin, simMax, ionT, magT, rootdir); }

DLLEXPORT void runNormalSimulationAPI(Simulation* sim, int iterations, int printEvery, const char* loadFileDir)
{
	double simMin{ sim->getSimMin() };
	double simMax{ sim->getSimMax() };

	sim->setBFieldModel("DipoleB", { 72.0 });

	sim->createParticleType("elec", { "vpara", "vperp", "s" }, MASS_ELECTRON, -1 * CHARGE_ELEM, 115200, 1, 2, RADIUS_EARTH, loadFileDir);
	sim->createParticleType("ions", { "vpara", "vperp", "s" }, MASS_PROTON,    1 * CHARGE_ELEM, 115200, 1, 2, RADIUS_EARTH, loadFileDir);

	sim->createTempSat(0, simMin * 0.999, true, "bottomElectrons");
	sim->createTempSat(1, simMin * 0.999, true, "bottomIons");
	sim->createTempSat(0, simMax * 1.001, false, "topElectrons");
	sim->createTempSat(1, simMax * 1.001, false, "topIons");

	sim->initializeSimulation();
	sim->copyDataToGPU();
	sim->iterateSimulation(iterations, printEvery);
	sim->copyDataToHost();
	sim->freeGPUMemory();
	sim->prepareResults(true);
}

DLLEXPORT void terminateSimulationAPI(Simulation* simulation) {
	delete simulation; }

DLLEXPORT void setBFieldModelAPI(Simulation* sim, const char* modelName, double arg1) {
	sim->setBFieldModel(modelName, { arg1 }); }

#ifndef DLLFILE
int main()
{
	Simulation* sim{ normalSimulationAPI() };
	delete sim;

	return 0;
}
#endif