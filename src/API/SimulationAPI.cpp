#include "API/SimulationAPI.h"

#include "utils/writeIOclasses.h"
#include "utils/string.h"
#include "ErrorHandling/simExceptionMacros.h"

#include <cmath>

using utils::fileIO::CSV;
using utils::string::strToDblVec;

///One liner functions
DLLEXP_EXTC double getSimulationTimeAPI(Simulation* simulation) {
	return simulation->simtime(); }

DLLEXP_EXTC double getDtAPI(Simulation* simulation) {
	return simulation->dt(); }

DLLEXP_EXTC double getSimMinAPI(Simulation* simulation) {
	return simulation->simMin(); }

DLLEXP_EXTC double getSimMaxAPI(Simulation* simulation) {
	return simulation->simMax(); }

DLLEXP_EXTC int getNumberOfParticleTypesAPI(Simulation* simulation) {
	return (int)(simulation->getNumberOfParticleTypes()); }

DLLEXP_EXTC int getNumberOfParticlesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return (int)(simulation->getNumberOfParticles(partInd))); }

DLLEXP_EXTC int getNumberOfAttributesAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return (int)(simulation->getNumberOfAttributes(partInd))); }

DLLEXP_EXTC const char* getParticleNameAPI(Simulation* simulation, int partInd) {
	SIM_API_EXCEP_CHECK(return simulation->getParticleName(partInd).c_str()); }

DLLEXP_EXTC const char* getSatelliteNameAPI(Simulation* simulation, int satInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteName(satInd).c_str()); }

DLLEXP_EXTC LogFile* getLogFilePointerAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(return simulation->log()); }


//Pointer one liners
DLLEXP_EXTC const double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData) {
	SIM_API_EXCEP_CHECK(return simulation->getParticleData(partIndex, originalData).at(attrIndex).data()); }


//Field tools
DLLEXP_EXTC double getBFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getBFieldAtS(s, time); }

DLLEXP_EXTC double getEFieldAtSAPI(Simulation* simulation, double s, double time) {
	return simulation->getEFieldAtS(s, time); }


//Simulation Management Function Wrappers
DLLEXP_EXTC Simulation* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir) {
	SIM_API_EXCEP_CHECK(return new Simulation(dt, simMin, simMax, rootdir)); }

DLLEXP_EXTC void initializeSimulationAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->initializeSimulation()); }

DLLEXP_EXTC void __iterateSimCPUAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(simulation->__iterateSimCPU(numberOfIterations, itersBtwCouts)); }

DLLEXP_EXTC void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(simulation->iterateSimulation(numberOfIterations, itersBtwCouts)); }

DLLEXP_EXTC void freeGPUMemoryAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->freeGPUMemory()); }

DLLEXP_EXTC void saveDataToDiskAPI(Simulation* simulation) {
	SIM_API_EXCEP_CHECK(simulation->saveDataToDisk()); }

DLLEXP_EXTC void terminateSimulationAPI(Simulation* simulation) {
	delete simulation; }

DLLEXP_EXTC Simulation* loadCompletedSimDataAPI(const char* fileDir) {
	SIM_API_EXCEP_CHECK( return new Simulation(fileDir) ); }

DLLEXP_EXTC void setupNormalSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir)
{
	SIM_API_EXCEP_CHECK(
	double simMin{ sim->simMin() };
	double simMax{ sim->simMax() };

	sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
	//sim->setBFieldModel("DipoleB", { 72.0 });
	//sim->addEFieldModel("QSPS", { 3185500.0, 6185500.0, 0.02, 6556500.0, 9556500.0, 0.04 });

	sim->createParticleType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, loadFileDir);

	sim->createTempSat(0, simMin * 0.999, true, "btmElec");
	sim->createTempSat(0, simMax * 1.001, false, "topElec");

	//sim->createTempSat(0, 1014252.60176003, false, "1e6ElecDown"); //altitude = 1000km, s = what you see to the left
	//sim->createTempSat(0, 1014252.60176003, true, "1e6ElecUp");
	sim->createTempSat(0, 3049829.25570638, false, "3e6ElecDown"); //altitude = 3000km
	sim->createTempSat(0, 3049829.25570638, true, "3e6ElecUp");
	sim->createTempSat(0, 4071307.04106411, false, "4e6ElecDown");   //altitude = 4000km
	sim->createTempSat(0, 4071307.04106411, true, "4e6ElecUp");
	); /* SIM_API_EXCEP_CHECK() */
}

DLLEXP_EXTC void runNormalSimulationAPI(Simulation* sim, int iterations, int printEvery)
{
	SIM_API_EXCEP_CHECK(
	sim->initializeSimulation();
	sim->iterateSimulation(iterations, printEvery);
	); /* SIM_API_EXCEP_CHECK() */
}

DLLEXP_EXTC void runSingleElectronAPI(Simulation* sim, double vpara, double vperp, double s, double t_inc, int iterations, int printEvery)
{
	SIM_API_EXCEP_CHECK(
	sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
	sim->createParticleType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, 1, "", false);
	std::vector<std::vector<double>> attrs{ std::vector<std::vector<double>>({ {vpara}, {vperp}, {s}, {0.0}, {t_inc} }) };
	sim->particle("elec")->loadDataFromMem(attrs, true);
	sim->particle("elec")->loadDataFromMem(attrs, false);
	sim->createTempSat(0, sim->simMin() * 0.999, true, "btmElec");
	sim->initializeSimulation();
	sim->__iterateSimCPU(iterations, printEvery);
	); /* SIM_API_EXCEP_CHECK() */
}


//Fields management
DLLEXP_EXTC void setBFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString) {
	SIM_API_EXCEP_CHECK(sim->setBFieldModel(modelName, strToDblVec(doubleString))); }

DLLEXP_EXTC void addEFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString) {
	SIM_API_EXCEP_CHECK(sim->addEFieldModel(modelName, strToDblVec(doubleString))); }

//Particle functions
DLLEXP_EXTC void createParticleTypeAPI(Simulation* simulation, const char* name, double mass, double charge, long numParts, const char* loadFileDir) {
	SIM_API_EXCEP_CHECK(simulation->createParticleType(name, mass, charge, numParts, loadFileDir)); }


//Satellite functions
DLLEXP_EXTC void createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name) {
	SIM_API_EXCEP_CHECK(simulation->createTempSat(particleInd, altitude, upwardFacing, name)); }

DLLEXP_EXTC int  getNumberOfSatellitesAPI(Simulation* simulation) {
	return (int)(simulation->getNumberOfSatellites()); }

DLLEXP_EXTC const double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int msmtInd, int attributeInd) {
	SIM_API_EXCEP_CHECK(return simulation->getSatelliteData(satelliteInd).at(msmtInd).at(attributeInd).data()); }


//CSV functions
DLLEXP_EXTC void writeCommonCSVAPI(Simulation* simulation)
{
	SIM_API_EXCEP_CHECK(
		CSV csvtmp("./elecoutput.csv");
		std::vector<std::vector<double>> origData{ simulation->getParticleData(0, true) };
		csvtmp.add({ origData.at(0), origData.at(1), origData.at(2) }, { "vpara orig", "vperp orig", "s orig" });
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
	Simulation* sim{ createSimulationAPI(0.01, 101322.378940846, 19881647.2473464, "./out/") };
	setupNormalSimulationAPI(sim, 3456000, "./../_in/data");
	sim->addEFieldModel("QSPS", { 3185500.0, 6185500.0, 0.02, 6556500.0, 9556500.0, 0.04 });
	runNormalSimulationAPI(sim, 500, 50);

	terminateSimulationAPI(sim);
	); /* SIM_API_EXCEP_CHECK() */

	return 0;
}
#endif
