#include "Simulation/Simulation.h"

#include "utils/loopmacros.h"
#include "ErrorHandling/simExceptionMacros.h"

Simulation::Simulation(double dt, double simMin, double simMax, std::string saveRootDir) :
	dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, saveRootDir_m{ saveRootDir + "/" },
	simAttr_m{ std::make_unique<SimAttributes>(saveRootDir + "/Simulation.attr") },
	logFile_m{ std::make_unique<LogFile>(saveRootDir_m + "simulation.log", 20) }
{
	cerrLogOut.open(saveRootDir_m + "errors.log");
	std::cerr.rdbuf(cerrLogOut.rdbuf()); //set cerr output to "errors.log"

	simAttr_m->addData("Simulation", "", {}, {}, { "dt", "simMin", "simMax" }, { dt_m, simMin_m, simMax_m });
}

Simulation::~Simulation()
{
	std::cerr.rdbuf(cerrBufBak); //restore cerr to normal

	if (saveReady_m) { saveDataToDisk(); } //save data if it hasn't been done

	logFile_m->writeTimeDiffFromNow(0, "End Simulation Destructor");
}


void Simulation::printSimAttributes(int numberOfIterations, int itersBtwCouts, std::string GPUName) //protected
{
	//Sim Header (folder) printed from Python - move here eventually
	std::cout << "GPU Name:       " << GPUName << "\n";
	std::cout << "Sim between:    " << simMin_m << "m - " << simMax_m << "m" << std::endl;
	std::cout << "dt:             " << dt_m << "s" << std::endl;
	std::cout << "BField Model:   " << BFieldModel_m->name() << std::endl;
	std::cout << "EField Elems:   " << EFieldModel_m->getEElemsStr() << std::endl;
	std::cout << "Particles:      ";
	for (size_t iii = 0; iii < particles_m.size(); iii++) {
		std::cout << ((iii != 0) ? "                " : "") << particles_m.at(iii)->name() << ": #: " << particles_m.at(iii)->getNumberOfParticles() << ", loaded files?: " << (particles_m.at(iii)->getInitDataLoaded() ? "true" : "false") << std::endl;
	}
	std::cout << "Satellites:     ";
	for (int iii = 0; iii < getNumberOfSatellites(); iii++) {
		std::cout << ((iii != 0) ? "                " : "") << satellite(iii)->name() << ": alt: " << satellite(iii)->altitude() << " m, upward?: " << (satellite(iii)->upward() ? "true" : "false") << std::endl;
	}
	std::cout << "Iterations:     " << numberOfIterations << std::endl;
	std::cout << "Iters Btw Cout: " << itersBtwCouts << std::endl;
	std::cout << "Time to setup:  "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;
	std::cout << "===============================================================" << std::endl;
}


//Access functions
const std::vector<std::vector<double>>& Simulation::getParticleData(int partInd, bool originalData)
{
	if (partInd > (particles_m.size() - 1))
		throw std::out_of_range("Simulation::getParticleData: no particle at the specifed index " + std::to_string(partInd));
	return ((originalData) ? (particles_m.at(partInd)->data(true)) : (particles_m.at(partInd)->data(false)));
}

const std::vector<std::vector<double>>& Simulation::getSatelliteData(int satInd)
{
	if (satInd > (satPartPairs_m.size() - 1))
		throw std::out_of_range("Simulation::getSatelliteData: no satellite at the specifed index " + std::to_string(satInd));
	return satellite(satInd)->data();
}

Particle* Simulation::particle(std::string name) const
{
	for (auto& part : particles_m)
		if (part->name() == name)
			return part.get();
	throw std::invalid_argument("Simulation::particle: no particle of name " + name);
}

Satellite* Simulation::satellite(std::string name) const
{
	for (auto& satPart : satPartPairs_m)
		if (satPart->satellite->name() == name)
			return satPart->satellite.get();
	throw std::invalid_argument("Simulation::satellite: no satellite of name " + name);
}

//Fields
double Simulation::getBFieldAtS(double s, double time) const { return BFieldModel_m->getBFieldAtS(s, time); }
double Simulation::getEFieldAtS(double s, double time) const { return EFieldModel_m->getEFieldAtS(s, time); }

//Class creation functions
void Simulation::createParticleType(std::string name, double mass, double charge, long numParts, std::string loadFilesDir, bool save)
{
	if (simAttr_m == nullptr)
		save = false;

	for (size_t part = 0; part < particles_m.size(); part++)
		if (name == particles_m.at(part).get()->name())
			throw std::invalid_argument("Simulation::createParticleType: particle already exists with the name: " + name);

	logFile_m->writeLogFileEntry("Simulation::createParticleType: Particle Type Created: " + name + ": Mass: " + std::to_string(mass) + ", Charge: " + std::to_string(charge) + ", Number of Parts: " + std::to_string(numParts) + ", Files Loaded?: " + ((loadFilesDir != "") ? "True" : "False"));
	
	std::vector<std::string> attrNames{ "vpara", "vperp", "s", "t_inc", "t_esc" };

	if (save)
	{
		std::vector<std::string> attrLabels;
		for (size_t atr = 0; atr < attrNames.size(); atr++)
			attrLabels.push_back("attrName");
		attrLabels.push_back("loadFilesDir");
		
		std::vector<std::string> namesCopy{ attrNames };
		namesCopy.push_back(loadFilesDir);
		simAttr_m->addData("Particle", name, attrLabels, namesCopy, { "mass", "charge", "numParts" }, { mass, charge, (double)numParts });
	}

	std::shared_ptr<Particle> newPart{ std::make_unique<Particle>(name, attrNames, mass, charge, numParts) };

	if (loadFilesDir != "")
		newPart->loadDataFromDisk(loadFilesDir);

	newPart->__data(true).at(4) = std::vector<double>(newPart->getNumberOfParticles(), -1.0); //sets t_esc to -1.0 - i.e. particles haven't escaped yet
	particles_m.push_back(std::move(newPart));
}

void Simulation::createTempSat(std::string partName, double altitude, bool upwardFacing, std::string name)
{
	for (size_t partInd = 0; partInd < particles_m.size(); partInd++)
	{
		if (particle((int)partInd)->name() == partName)
		{
			createTempSat(partInd, altitude, upwardFacing, name);
			return;
		}
	}
	throw std::invalid_argument("Simulation::createTempSat: no particle of name " + name);
}

void Simulation::createTempSat(int partInd, double altitude, bool upwardFacing, std::string name)
{//"Temp Sats" are necessary to ensure particles are created before their accompanying satellites
	if (initialized_m)
		throw std::runtime_error("Simulation::createTempSat: initializeSimulation has already been called, no satellite will be created of name " + name);
	if (partInd >= particles_m.size())
		throw std::out_of_range("Simulation::createTempSat: no particle at the specifed index " + std::to_string(partInd));

	tempSats_m.push_back(std::make_unique<TempSat>(partInd, altitude, upwardFacing, name));
}

void Simulation::createSatellite(TempSat* tmpsat, bool save) //protected
{
	int partInd{ tmpsat->particleInd };
	double altitude{ tmpsat->altitude };
	bool upwardFacing{ tmpsat->upwardFacing };
	std::string name{ tmpsat->name };

	if (particles_m.size() <= partInd)
		throw std::out_of_range("createSatellite: no particle at the specifed index " + std::to_string(partInd));
	if (particles_m.at(partInd)->getCurrDataGPUPtr() == nullptr)
		throw std::runtime_error("createSatellite: pointer to GPU data is a nullptr of particle " + particles_m.at(partInd)->name() + " - that's just asking for trouble");
	if (simAttr_m == nullptr)
		save = false;

	if (save) { simAttr_m->addData("Satellite", name, {}, {}, { "partInd", "altitude", "upwardFacing" }, { (double)partInd, altitude, (double)upwardFacing }); }

	logFile_m->writeLogFileEntry("Simulation::createSatellite: Created Satellite: " + name + ", Particle tracked: " + particles_m.at(partInd)->name() + ", Altitude: " + std::to_string(altitude) + ", " + ((upwardFacing) ? "Upward" : "Downward") + " Facing Detector");

	std::vector<std::string> attrNames{ "vpara", "vperp", "s", "time", "index" };
	std::shared_ptr<Particle>  part{ particles_m.at(partInd) };
	std::unique_ptr<Satellite> sat{ std::make_unique<Satellite>(name, attrNames, altitude, upwardFacing, part->getNumberOfParticles(), part->getCurrDataGPUPtr()) };
	satPartPairs_m.push_back(std::make_unique<SatandPart>(std::move(sat), std::move(part)));
}

void Simulation::setBFieldModel(std::string name, std::vector<double> args, bool save)
{//add log file messages
	if (BFieldModel_m)
		throw std::invalid_argument("Simulation::setBFieldModel: trying to assign B Field Model when one is already assigned - existing: " + BFieldModel_m->name() + ", attempted: " + name);
	if (args.empty())
		throw std::invalid_argument("Simulation::setBFieldModel: no arguments passed in");
	if (simAttr_m == nullptr)
		save = false;
	
	std::vector<std::string> attrNames;

	if (name == "DipoleB")
	{
		if (args.size() == 1)
		{ //for defaults in constructor of DipoleB
			BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
			args.push_back(((DipoleB*)BFieldModel_m.get())->getErrTol());
			args.push_back(((DipoleB*)BFieldModel_m.get())->getds());
		}
		else if (args.size() == 3)
			BFieldModel_m = std::make_unique<DipoleB>(args.at(0), args.at(1), args.at(2));
		else
			throw std::invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleB: " + std::to_string(args.size()));

		attrNames = { "ILAT", "ds", "errTol" };
	}
	else if (name == "DipoleBLUT")
	{
		if (args.size() == 3)
			BFieldModel_m = std::make_unique<DipoleBLUT>(args.at(0), simMin_m, simMax_m, args.at(1), (int)args.at(2));
		else
			throw std::invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleBLUT: " + std::to_string(args.size()));

		attrNames = { "ILAT", "ds", "numMsmts" };
	}
	else
	{
		std::cout << "Not sure what model is being referenced.  Using DipoleB instead of " << name << std::endl;
		BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
		args.resize(3);
		args.at(1) = ((DipoleB*)BFieldModel_m.get())->getErrTol();
		args.at(1) = ((DipoleB*)BFieldModel_m.get())->getds();
		attrNames = { "ILAT", "ds", "errTol" };
	}

	BFieldModel_d = BFieldModel_m->getPtrGPU();
	if (save) { simAttr_m->addData("BField", name, {}, {}, attrNames, args); }
}

void Simulation::addEFieldModel(std::string name, std::vector<double> args, bool save)
{
	if (EFieldModel_m == nullptr)
		EFieldModel_m = std::make_unique<EField>();
	if (simAttr_m == nullptr)
		save = false;

	std::vector<std::string> attrNames;

	if (name == "QSPS") //need to check to make sure args is formed properly, as well as save to disk
	{
		if (args.size() % 3 != 0) { throw std::invalid_argument("Simulation::addEFieldModel: QSPS: Argument vector is improperly formed.  Proper format is: { altMin, altMax, mag(, altMin, altMax, mag...) }"); }
		
		std::vector<meters> altMin;
		std::vector<meters> altMax;
		std::vector<Vperm> mag;

		for (size_t entry = 0; entry < args.size() / 3; entry++)
		{
			altMin.push_back(args.at(3 * entry));
			altMax.push_back(args.at(3 * entry + 1));
			mag.push_back(args.at(3 * entry + 2));
			attrNames.push_back("altMin");
			attrNames.push_back("altMax");
			attrNames.push_back("magnitude");
		}
		EFieldModel_m->add(std::make_unique<QSPS>(altMin, altMax, mag));
	}
	else if (name == "AlfvenLUT")
	{
		std::cout << "AlfvenLUT not implemented quite yet.  Returning." << std::endl;
		return;
	}
	/*else if (name == "AlfvenCompute")
	{
		std::cout << "AlfvenCompute not implemented quite yet.  Returning." << std::endl;
		return;
	}*/

	if (save) { simAttr_m->addData("EField", name, {}, {}, attrNames, args); }
}


//Other utilities
void Simulation::saveDataToDisk()
{
	if (!initialized_m)
		throw SimFatalException("Simulation::saveDataToDisk: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);
	if (!saveReady_m)
		throw std::runtime_error("Simulation::saveDataToDisk: simulation not iterated and/or copied to host with iterateSmiulation()");

	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->saveDataToDisk(saveRootDir_m + "/bins/particles_init/", true));
	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->saveDataToDisk(saveRootDir_m + "/bins/particles_final/", false));
	LOOP_OVER_1D_ARRAY(getNumberOfSatellites(), satellite(iii)->saveDataToDisk(saveRootDir_m + "/bins/satellites/"));

	simAttr_m = nullptr;
	saveReady_m = false;
}

void Simulation::resetSimulation(bool fields)
{
	while (satPartPairs_m.size() != 0)
		satPartPairs_m.pop_back();
	while (particles_m.size() != 0)
		particles_m.pop_back();

	if (fields)
	{
		BFieldModel_m.reset();
		EFieldModel_m.reset();
	}
}
