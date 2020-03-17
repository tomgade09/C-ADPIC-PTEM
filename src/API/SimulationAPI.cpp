#include "API/SimulationAPI.h"

#include "utils/writeIOclasses.h"
#include "utils/strings.h"
#include "ErrorHandling/simExceptionMacros.h"

#include <cmath>

using utils::fileIO::CSV;
using utils::strings::strToDblVec;

//Simulation Management Functions
DLLEXP_EXTC Sim* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir)
{
	SIM_API_EXCEP_CHECK(return new Sim(dt, simMin, simMax, rootdir));
	return nullptr; //if above fails
}

DLLEXP_EXTC Sim* loadCompletedSimDataAPI(const char* fileDir)
{
	SIM_API_EXCEP_CHECK(return new Sim(fileDir));
	return nullptr; //if above fails
}

DLLEXP_EXTC void initializeSimulationAPI(Sim* sim) {
	SIM_API_EXCEP_CHECK(sim->initializeSimulation()); }

DLLEXP_EXTC void iterateSimCPUAPI(Sim* sim, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(sim->__iterateSimCPU(numberOfIterations, itersBtwCouts)); }

DLLEXP_EXTC void iterateSimulationAPI(Sim* sim, int numberOfIterations, int itersBtwCouts) {
	SIM_API_EXCEP_CHECK(sim->iterateSimulation(numberOfIterations, itersBtwCouts)); }

DLLEXP_EXTC void freeGPUMemoryAPI(Sim* sim) {
	SIM_API_EXCEP_CHECK(sim->freeGPUMemory()); }

DLLEXP_EXTC void saveDataToDiskAPI(Sim* sim) {
	SIM_API_EXCEP_CHECK(sim->saveDataToDisk()); }

DLLEXP_EXTC void terminateSimulationAPI(Sim* sim) {
	SIM_API_EXCEP_CHECK(delete sim); }

DLLEXP_EXTC void setupExampleSimulationAPI(Sim* sim, int numParts, const char* loadFileDir)
{
	SIM_API_EXCEP_CHECK(
		sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
		//sim->setBFieldModel("DipoleB", { 72.0 });
		//sim->addEFieldModel("QSPS", { 3185500.0, 6185500.0, 0.02, 6556500.0, 9556500.0, 0.04 });

		sim->createParticlesType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, loadFileDir);

		sim->createTempSat(0, sim->simMin(), true, "btmElec");
		sim->createTempSat(0, sim->simMax(), false, "topElec");
		sim->createTempSat(0, 4071307.04106411, false, "4e6ElecUpg");
		sim->createTempSat(0, 4071307.04106411, true, "4e6ElecDng");

		sim->particles(0)->setParticlesSource_s(sim->simMin(), sim->simMax());
	); /* SIM_API_EXCEP_CHECK() */
}

DLLEXP_EXTC void setupSingleElectronAPI(Sim* sim, double vpara, double vperp, double s, double t_inc)
{ //check that satellites and B/E fields have been created here??
	SIM_API_EXCEP_CHECK(
		std::vector<std::vector<double>> attrs{ std::vector<std::vector<double>>({ { vpara },{ vperp },{ s },{ 0.0 },{ t_inc } }) };
		sim->particles("elec")->loadDataFromMem(attrs, true);
		sim->particles("elec")->loadDataFromMem(attrs, false);
	); /* SIM_API_EXCEP_CHECK() */
}


//Field Management Functions
DLLEXP_EXTC void setBFieldModelAPI(Sim* sim, const char* modelName, const char* doubleString) {
	SIM_API_EXCEP_CHECK(sim->setBFieldModel(modelName, strToDblVec(doubleString))); }

DLLEXP_EXTC void addEFieldModelAPI(Sim* sim, const char* modelName, const char* doubleString) {
	SIM_API_EXCEP_CHECK(sim->addEFieldModel(modelName, strToDblVec(doubleString))); }

DLLEXP_EXTC double getBFieldAtSAPI(Sim* sim, double s, double time)
{
	SIM_API_EXCEP_CHECK(return sim->getBFieldAtS(s, time));
	return 0.0; //if above fails
}

DLLEXP_EXTC double getEFieldAtSAPI(Sim* sim, double s, double time)
{
	SIM_API_EXCEP_CHECK(return sim->getEFieldAtS(s, time));
	return 0.0; //if above fails
}


//Particles Management Functions
DLLEXP_EXTC void createParticlesTypeAPI(Sim* sim, const char* name, double mass, double charge, long numParts, const char* loadFileDir) {
	SIM_API_EXCEP_CHECK(sim->createParticlesType(name, mass, charge, numParts, loadFileDir)); }


//Satellite Management Functions
DLLEXP_EXTC void createSatelliteAPI(Sim* sim, int particleInd, double altitude, bool upwardFacing, const char* name) {
	SIM_API_EXCEP_CHECK(sim->createTempSat(particleInd, altitude, upwardFacing, name)); }

DLLEXP_EXTC int  getNumberOfSatellitesAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return (int)(sim->getNumberOfSatellites()));
	return -1; //if above fails
}

DLLEXP_EXTC const double* getSatelliteDataPointersAPI(Sim* sim, int satelliteInd, int attributeInd)
{
	SIM_API_EXCEP_CHECK(return sim->getSatelliteData(satelliteInd).at(attributeInd).data());
	return nullptr; //if above fails
}

DLLEXP_EXTC int getPartIndOfSatAPI(Sim* sim, int satelliteInd)
{
	SIM_API_EXCEP_CHECK(return sim->getParticleIndexOfSat(satelliteInd));
	return -1; //if above fails
}


//Access Functions
DLLEXP_EXTC double getSimTimeAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return sim->simtime());
	return -1.0; //if above fails
}

DLLEXP_EXTC double getDtAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return sim->dt());
	return -1.0; //if above fails
}

DLLEXP_EXTC double getSimMinAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return sim->simMin());
	return -1.0; //if above fails
}

DLLEXP_EXTC double getSimMaxAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return sim->simMax());
	return -1.0; //if above fails
}

DLLEXP_EXTC int getNumberOfParticleTypesAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(return (int)(sim->getNumberOfParticleTypes()));
	return -1; //if above fails
}

DLLEXP_EXTC int getNumberOfParticlesAPI(Sim* sim, int partInd)
{
	SIM_API_EXCEP_CHECK(return (int)(sim->getNumberOfParticles(partInd)));
	return -1; //if above fails
}

DLLEXP_EXTC int getNumberOfAttributesAPI(Sim* sim, int partInd)
{
	SIM_API_EXCEP_CHECK(return (int)(sim->getNumberOfAttributes(partInd)));
	return -1; //if above fails
}

DLLEXP_EXTC const char* getParticlesNameAPI(Sim* sim, int partInd)
{
	SIM_API_EXCEP_CHECK(return sim->getParticlesName(partInd).c_str());
	return nullptr; //if above fails
}

DLLEXP_EXTC const char* getSatelliteNameAPI(Sim* sim, int satInd)
{
	SIM_API_EXCEP_CHECK(return sim->getSatelliteName(satInd).c_str());
	return nullptr; //if above fails
}

DLLEXP_EXTC const double* getPointerToParticlesAttributeArrayAPI(Sim* sim, int partIndex, int attrIndex, bool originalData)
{
	SIM_API_EXCEP_CHECK(return sim->getParticleData(partIndex, originalData).at(attrIndex).data());
	return nullptr; //if above fails
}


//CSV functions
DLLEXP_EXTC void writeCommonCSVAPI(Sim* sim)
{
	SIM_API_EXCEP_CHECK(
		CSV csvtmp("./elecoutput.csv");
		std::vector<std::vector<double>> origData{ sim->getParticleData(0, true) };
		csvtmp.add({ origData.at(0), origData.at(1), origData.at(2) }, { "vpara orig", "vperp orig", "s orig" });
		csvtmp.addspace();
		
		std::vector<std::vector<double>> btmElecData{ sim->getSatelliteData(0).at(0) };
		csvtmp.add({ btmElecData.at(3), btmElecData.at(0), btmElecData.at(1), btmElecData.at(2) },
			{ "t_esc btm", "vpara btm", "vperp btm", "s btm" });
		csvtmp.addspace();

		std::vector<std::vector<double>> topElecData{ sim->getSatelliteData(1).at(0) };
		csvtmp.add({ topElecData.at(3), topElecData.at(0), topElecData.at(1), topElecData.at(2) },
			{ "t_esc top", "vpara top", "vperp top", "s top" });
		csvtmp.addspace();

		std::vector<std::vector<double>> energyPitch(2, std::vector<double>(origData.at(0).size()));
		for (size_t elem = 0; elem < energyPitch.at(0).size(); elem++)
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
	/*SIM_API_EXCEP_CHECK(
		auto sim{ std::make_unique<Simulation>(0.01, 101322.378940846, 19881647.2473464, "./out/") };
		setupExampleSimulationAPI(sim.get(), 3456000, "./../_in/data");
		//sim->addEFieldModel("QSPS", { 3185500.0, 6185500.0, 0.02, 6556500.0, 9556500.0, 0.04 });
		sim->initializeSimulation();
		sim->iterateSimulation(5000, 500);
	);*/ /* SIM_API_EXCEP_CHECK() */

	auto sim{ "../_dataout/200314_11.10.53/" };

	return 0;
}
#endif
