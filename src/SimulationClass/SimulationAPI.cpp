#include "SimulationClass\SimulationAPI.h"

///One liner functions
DLLEXPORT double getSimulationTimeAPI(Simulation* simulation) {
	return simulation->simtime(); }

DLLEXPORT double getDtAPI(Simulation* simulation) {
	return simulation->dt(); }

DLLEXPORT double getSimMinAPI(Simulation* simulation) {
	return simulation->simMin(); }

DLLEXPORT double getSimMaxAPI(Simulation* simulation) {
	return simulation->simMax(); }

DLLEXPORT int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return (int)(simulation->getNumberOfParticleTypes()); }

DLLEXPORT int getNumberOfParticlesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return (int)(simulation->getNumberOfParticles(partInd))); }

DLLEXPORT int getNumberOfAttributesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return (int)(simulation->getNumberOfAttributes(partInd))); }

DLLEXPORT const char* getParticleNameAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return simulation->getParticleName(partInd).c_str()); }

DLLEXPORT const char* getSatelliteNameAPI(Simulation* simulation, int satInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteName(satInd).c_str()); }

DLLEXPORT LogFile* getLogFilePointerAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(return simulation->logPtr()); }


//Pointer one liners
DLLEXPORT const double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	SIM_API_EXCEP_CHECK(return simulation->getParticleData(partIndex, originalData).at(attrIndex).data()); }


//Field tools
DLLEXPORT double getBFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getBFieldAtS(s, time); }

DLLEXPORT double getEFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getEFieldAtS(s, time); }


//Simulation Management Function Wrappers
DLLEXPORT Simulation* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir) {
	SIM_API_EXCEP_CHECK(return new Simulation(dt, simMin, simMax, rootdir)); }

DLLEXPORT void initializeSimulationAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->initializeSimulation()); }

DLLEXPORT void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(simulation->iterateSimulation(numberOfIterations, itersBtwCouts)); }

DLLEXPORT void freeGPUMemoryAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->freeGPUMemory()); }

DLLEXPORT void saveDataToDiskAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->saveDataToDisk()); }

DLLEXPORT void terminateSimulationAPI(Simulation* simulation) {
	delete simulation; }

DLLEXPORT Simulation* loadCompletedSimDataAPI(const char* fileDir) {
	SIM_API_EXCEP_CHECK( return new Simulation(fileDir) ); }

DLLEXPORT void setupNormalSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir)
{
	SIM_API_EXCEP_CHECK(
	double simMin{ sim->simMin() };
	double simMax{ sim->simMax() };

	sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
	//sim->setBFieldModel("DipoleB", { 72.0 });
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
	SIM_API_EXCEP_CHECK(
	sim->initializeSimulation();
	sim->iterateSimulation(iterations, printEvery);
	); /* SIM_API_EXCEP_CHECK() */
}


//Fields management
DLLEXPORT void setBFieldModelAPI(Simulation* sim, const char* modelName, double arg1) {
	SIM_API_EXCEP_CHECK(sim->setBFieldModel(modelName, { arg1 })); }


//Particle functions
DLLEXPORT void createParticleTypeAPI(Simulation* simulation, const char* name, const char* attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, const char* loadFileDir) {
	SIM_API_EXCEP_CHECK(simulation->createParticleType(name, utils::string::charToStrVec(attrNames), mass, charge, numParts, posDims, velDims, normFactor, loadFileDir)); }


//Satellite functions
DLLEXPORT void createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name) {
	SIM_API_EXCEP_CHECK(simulation->createTempSat(particleInd, altitude, upwardFacing, name)); }

DLLEXPORT int  getNumberOfSatellitesAPI(Simulation* simulation) {
	return (int)(simulation->getNumberOfSatellites()); }

DLLEXPORT const double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int msmtInd, int attributeInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteData(satelliteInd).at(msmtInd).at(attributeInd).data()); }


//CSV functions
DLLEXPORT void writeCommonCSVAPI(Simulation* simulation)
{
	SIM_API_EXCEP_CHECK(
		utils::write::CSV csvtmp("./elecoutput.csv");
		std::vector<std::vector<double>> origData{ simulation->getParticleData(0, true) };
		csvtmp.add(origData, { "vpara orig", "vperp orig", "s orig" });
		csvtmp.addspace();
		
		std::vector<std::vector<double>> btmElecData{ simulation->getSatelliteData(0).at(0) };
		csvtmp.add({ btmElecData.at(3), btmElecData.at(0), btmElecData.at(1), btmElecData.at(2) },
			{ "t_esc btm", "vpara btm", "vperp btm", "s btm" });
		csvtmp.addspace();

		std::vector<std::vector<double>> topElecData{ simulation->getSatelliteData(1).at(0) };
		csvtmp.add({ topElecData.at(3), topElecData.at(0), topElecData.at(1), topElecData.at(2) },
			{ "t_esc top", "vpara top", "vperp top", "s top" });
		csvtmp.addspace();

		std::vector<std::vector<double>> energyPitch(2, std::vector<double>(origData.at(0).size()));
		for (int elem = 0; elem < energyPitch.at(0).size(); elem++)
		{
			energyPitch.at(0).at(elem) = 0.5 * MASS_ELECTRON * (pow(origData.at(0).at(elem), 2) + pow(origData.at(1).at(elem), 2)) / JOULE_PER_EV;
			energyPitch.at(1).at(elem) = atan2(abs(origData.at(1).at(elem)), -origData.at(0).at(elem)) / RADS_PER_DEG;
		}
		csvtmp.add(energyPitch, { "Energy (eV)", "Pitch Angle" });
	);
}

#ifdef _DEBUG
int main()
{
	SIM_API_EXCEP_CHECK(
	Simulation* sim{ createSimulationAPI(0.01, 101322.378940846, 19881647.2473464, 2.5, 1000.0, "./out/") };
	setupNormalSimulationAPI(sim, 3456000, "./../in/data");
	runNormalSimulationAPI(sim, 25000, 500);

	terminateSimulationAPI(sim);
	); /* SIM_API_EXCEP_CHECK() */

	return 0;
}
#endif