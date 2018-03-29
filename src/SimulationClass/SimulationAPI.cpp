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

DLLEXPORT const char* getParticleNameAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return simulation->getParticleName(partInd)); }

DLLEXPORT const char* getSatelliteNameAPI(Simulation* simulation, int satInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteName(satInd)); }

DLLEXPORT bool areResultsPreparedAPI(Simulation* simulation) { //do I even need this?  maybe change way results are prepared
	return simulation->areResultsPrepared(); }

DLLEXPORT LogFile* getLogFilePointerAPI(Simulation* simulation) {
	return simulation->getLogFilePointer(); }


//Pointer one liners
DLLEXPORT double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	SIM_API_EXCEP_CHECK(return simulation->getPointerToParticleAttributeArray(partIndex, attrIndex, originalData)); }


//Field tools
DLLEXPORT double getBFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getBFieldAtS(s, time); }

DLLEXPORT double getEFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getEFieldAtS(s, time); }


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

DLLEXPORT Simulation* loadCompletedSimDataAPI(const char* fileDir, const char* bFieldModel, const char* eFieldElems, const char* partNames, const char* satNames) {
	SIM_API_EXCEP_CHECK( return new PreviousSimulation(fileDir, bFieldModel, constCharToStrVec(eFieldElems), constCharToStrVec(partNames), constCharToStrVec(satNames)) ); }


//Particle functions
DLLEXPORT void createParticleTypeAPI(Simulation* simulation, const char* name, const char* attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, const char* loadFileDir) {
	SIM_API_EXCEP_CHECK(simulation->createParticleType(name, constCharToStrVec(attrNames), mass, charge, numParts, posDims, velDims, normFactor, loadFileDir)); }


//Simulation creation and deletion
DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, double ionT, double magT, const char* rootdir) {
	SIM_API_EXCEP_CHECK(return new Simulation(dt, simMin, simMax, ionT, magT, rootdir)); }

DLLEXPORT void setupNormalSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir)
{
	SIM_API_EXCEP_CHECK(\
	double simMin{ sim->getSimMin() };
	double simMax{ sim->getSimMax() };

	sim->setBFieldModel("DipoleBLUT", { 72.0 });
	//sim->addEFieldModel("QSPS", { 0.0 }, "3185500.0,6185500.0,6556500.0,9556500.0", "0.02,0.04");

	sim->createParticleType("elec", { "vpara", "vperp", "s" }, MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, 1, 2, RADIUS_EARTH, loadFileDir);
	//sim->createParticleType("ions", { "vpara", "vperp", "s" }, MASS_PROTON,    1 * CHARGE_ELEM, numParts, 1, 2, RADIUS_EARTH, loadFileDir);

	sim->createTempSat(0, simMin * 0.999, true, "btmElec");
	//sim->createTempSat(1, simMin * 0.999, true, "btmIons");
	sim->createTempSat(0, simMax * 1.001, false, "topElec");
	//sim->createTempSat(1, simMax * 1.001, false, "topIons");

	sim->createTempSat(0, 1014252.60176003, false, "1e6ElecDown"); //altitude = 1000km, s = what you see to the left
	sim->createTempSat(0, 1014252.60176003, true, "1e6ElecUp");
	sim->createTempSat(0, 3049829.25570638, false, "3e6ElecDown"); //altitude = 3000km
	sim->createTempSat(0, 3049829.25570638, true, "3e6ElecUp");
	sim->createTempSat(0, 4071307.04106411, false, "4e6ElecDown");   //altitude = 4000km
	sim->createTempSat(0, 4071307.04106411, true, "4e6ElecUp");
	); /* SIM_API_EXCEP_CHECK() */
}

DLLEXPORT void runNormalSimulationAPI(Simulation* sim, int iterations, int printEvery)
{
	SIM_API_EXCEP_CHECK( \
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

#ifdef _DEBUG
int main()
{
	SIM_API_EXCEP_CHECK( \
	Simulation* sim{ createSimulationAPI(0.01, 101322.378940846, 19881647.2473464, 2.5, 1000.0, "./out/") };
	setupNormalSimulationAPI(sim, 3456000, "./../in/data");
	runNormalSimulationAPI(sim, 25000, 500);

	terminateSimulationAPI(sim);
	); /* SIM_API_EXCEP_CHECK() */

	return 0;
}
#endif