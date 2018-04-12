#include "SimulationClass\Simulation.h"

//Access functions
const std::vector<std::vector<double>>& Simulation::getParticleData(int partInd, bool originalData)
{
	if (partInd > (particles_m.size() - 1))
		throw std::out_of_range("Simulation::getParticleData: no particle at the specifed index " + std::to_string(partInd));

	return ((originalData) ? (particles_m.at(partInd)->getOrigData()) : (particles_m.at(partInd)->getCurrData()));
}

const std::vector<std::vector<std::vector<double>>>& Simulation::getSatelliteData(int satInd)
{
	if (satInd > (satellites_m.size() - 1))
		throw std::out_of_range("Simulation::getSatelliteData: no satellite at the specifed index " + std::to_string(satInd));

	return satellites_m.at(satInd)->satellite->data();
}


//Class creation functions
void Simulation::createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir)
{
	//some sort of debug message??  Necessary for daily use?
	logFile_m->writeLogFileEntry("Simulation::createParticleType: Particle Type Created: " + name + ": Mass: " + std::to_string(mass) + ", Charge: " + std::to_string(charge) + ", Number of Parts: " + std::to_string(numParts) + ", Pos Dimensions: " + std::to_string(posDims) + ", Vel Dimensions: " + std::to_string(velDims) + ", Files Loaded?: " + ((loadFilesDir != "") ? "True" : "False"));
	
	fileIO::writeAttrsToFiles({ mass, charge, (double)numParts, (double)posDims, (double)velDims, normFactor },
		{ "mass", "charge", "numParts", "posDims", "velDims", "normFactor", "attrNames", attrNames.at(0), attrNames.at(1), attrNames.at(2), name }, "Particle_" + name, saveRootDir_m + "/_chars/");

	std::shared_ptr<Particle> newPart{ std::make_shared<Particle>(name, attrNames, mass, charge, numParts, posDims, velDims, normFactor) };

	if (loadFilesDir != "")
		newPart->loadDataFromDisk(loadFilesDir);

	particles_m.push_back(std::move(newPart));
}

void Simulation::createTempSat(int partInd, double altitude, bool upwardFacing, std::string name)
{//"Temp Sats" are necessary to ensure particles are created before their accompanying satellites
	if (initialized_m)
		throw std::runtime_error("Simulation::createTempSat: initializeSimulation has already been called, no satellite will be created of name " + name);
	if (partInd >= particles_m.size())
		throw std::out_of_range("Simulation::createTempSat: no particle at the specifed index " + std::to_string(partInd));

	tempSats_m.push_back(std::move(std::make_unique<TempSat>(partInd, altitude, upwardFacing, name)));
}

void Simulation::createSatellite(TempSat* tmpsat) //protected
{
	int partInd{ tmpsat->particleInd };
	double altitude{ tmpsat->altitude };
	bool upwardFacing{ tmpsat->upwardFacing };
	std::string name{ tmpsat->name };

	if (particles_m.size() <= partInd)
		throw std::out_of_range("createSatellite: no particle at the specifed index " + std::to_string(partInd));
	if (particles_m.at(partInd)->getCurrDataGPUPtr() == nullptr)
		throw std::runtime_error("createSatellite: pointer to GPU data is a nullptr of particle " + particles_m.at(partInd)->name() + " - that's just asking for trouble");

	fileIO::writeAttrsToFiles({ (double)partInd, altitude, (double)upwardFacing },
		{ "partInd", "altitude", "upwardFacing", name }, "Satellite_" + name, saveRootDir_m + "/_chars/");

	logFile_m->writeLogFileEntry("Simulation::createSatellite: Created Satellite: " + name + ", Particle tracked: " + particles_m.at(partInd)->name() + ", Altitude: " + std::to_string(altitude) + ", " + ((upwardFacing) ? "Upward" : "Downward") + " Facing Detector");

	std::shared_ptr<Particle> tmpPart{ particles_m.at(partInd) };
	std::unique_ptr<Satellite> newSat{ std::make_unique<Satellite>(altitude, upwardFacing, tmpPart->getNumberOfAttributes(), tmpPart->getNumberOfParticles(), tmpPart->getCurrDataGPUPtr(), name) };
	satellites_m.push_back(std::move(std::make_unique<SatandPart>(std::move(newSat), std::move(tmpPart))));
}

void Simulation::setBFieldModel(std::string name, std::vector<double> args, bool save)
{//add log file messages
	if (BFieldModel_m)
		throw std::invalid_argument("Simulation::setBFieldModel: trying to assign B Field Model when one is already assigned - existing: " + BFieldModel_m->getName() + ", attempted: " + name);
	if (args.empty())
		throw std::invalid_argument("Simulation::setBFieldModel: no arguments passed in");

	std::string attrsDir{ saveRootDir_m + "/_chars/" };
	std::vector<std::string> names;

	if (name == "DipoleB")
	{
		if (args.size() == 1)
		{ //for defaults in constructor of DipoleB
			BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
			args.push_back(BFieldModel_m->getErrTol());
			args.push_back(BFieldModel_m->getds());
		}
		else if (args.size() == 3)
			BFieldModel_m = std::make_unique<DipoleB>(args.at(0), args.at(1), args.at(2));
		else
			throw std::invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleB: " + std::to_string(args.size()));

		names = { "ILAT", "ds", "errTol" };
	}
	else if (name == "DipoleBLUT")
	{
		if (args.size() == 3)
			BFieldModel_m = std::make_unique<DipoleBLUT>(args.at(0), simMin_m, simMax_m, args.at(1), (int)args.at(2));
		else
			throw std::invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleBLUT: " + std::to_string(args.size()));

		names = { "ILAT", "ds", "numMsmts" };
	}
	else if (name == "IGRF")
	{
		//BFieldModel_m = std::make_unique<IGRFB>(args.at(0));
		std::cout << "IGRF not implemented yet!! :D  Using DipoleB" << std::endl;
		BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
		args.resize(3);
		args.at(1) = BFieldModel_m->getErrTol();
		args.at(2) = BFieldModel_m->getds();
		names = { "ILAT", "ds", "errTol" };
	}
	else if (name == "InvRCubedB")
	{
		//BFieldModel_m = std::make_unique<InvRCubedB>(args.at(0));
		std::cout << "InvRCubed not implemented yet!! :D  Using DipoleB" << std::endl;
		BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
		args.resize(3);
		args.at(1) = BFieldModel_m->getErrTol();
		args.at(2) = BFieldModel_m->getds();
		names = { "ILAT", "ds", "errTol" };
	}
	else
	{
		std::cout << "Not sure what model is being referenced.  Using DipoleB instead of " << name << std::endl;
		BFieldModel_m = std::make_unique<DipoleB>(args.at(0));
		args.resize(3);
		args.at(1) = BFieldModel_m->getErrTol();
		args.at(2) = BFieldModel_m->getds();
		names = { "ILAT", "ds", "errTol" };
	}

	BFieldModel_d = BFieldModel_m->getPtrGPU();
	if (save) { fileIO::writeAttrsToFiles(args, names, "BField_" + name, attrsDir); }
}

void Simulation::addEFieldModel(std::string name, std::vector<std::vector<double>> args)
{
	throw std::exception("addEFieldModel: need to code saving parameters");

	if (EFieldModel_m == nullptr)
		EFieldModel_m = std::make_unique<EField>();

	if (name == "QSPS") //need to check to make sure args is formed properly, as well as save to disk
		EFieldModel_m->add(std::make_unique<QSPS>(args.at(0), args.at(1), args.at(2)));
	else if (name == "AlfvenLUT")
	{
		std::cout << "AlfvenLUT not implemented quite yet.  Returning." << std::endl;
		return;
	}
	else if (name == "AlfvenCompute")
	{
		std::cout << "AlfvenCompute not implemented quite yet.  Returning." << std::endl;
		return;
	}
}


//vperp <-> mu
void Simulation::convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, std::vector<double>& t, double mass)
{
	LOOP_OVER_1D_ARRAY(vperp.size(), vperp.at(iii) = 0.5 * mass * vperp.at(iii) * vperp.at(iii) / BFieldModel_m->getBFieldAtS(s.at(iii), t.at(iii)));
}

void Simulation::convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, double mass)
{
	std::vector<double> t((int)vperp.size()); //creates vector of zeroes
	convertVPerpToMu(vperp, s, t, mass);
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, std::vector<double>& t, double mass)
{
	LOOP_OVER_1D_ARRAY(mu.size(), if (mu.at(iii) != 0.0) { mu.at(iii) = sqrt(2 * mu.at(iii) * BFieldModel_m->getBFieldAtS(s.at(iii), t.at(iii)) / mass); });
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, double mass)
{
	std::vector<double> t((int)mu.size()); //creates vector of zeroes
	convertMuToVPerp(mu, s, t, mass);
}


//Other utilities
void Simulation::saveDataToDisk()
{
	if (!initialized_m)
		throw SimFatalException("Simulation::saveDataToDisk: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);
	if (!saveReady_m)
		throw SimFatalException("Simulation::saveDataToDisk: simulation not iterated and/or copied to host with iterateSmiulation()", __FILE__, __LINE__);

	LOOP_OVER_1D_ARRAY(particles_m.size(), particles_m.at(iii)->saveDataToDisk("./bins/particles_init/", true));
	LOOP_OVER_1D_ARRAY(particles_m.size(), particles_m.at(iii)->saveDataToDisk("./bins/particles_final/", false));
	LOOP_OVER_1D_ARRAY(satellites_m.size(), satellites_m.at(iii)->satellite->saveDataToDisk(saveRootDir_m + "/bins/satellites/", { "vpara", "vperp", "s", "time", "index" }));

	saveReady_m = false;
}

void Simulation::resetSimulation(bool fields)
{
	for (int iii = 0; iii < satellites_m.size(); iii++)
		satellites_m.pop_back();
	for (int iii = 0; iii < particles_m.size(); iii++)
		particles_m.pop_back();

	if (fields)
	{
		BFieldModel_m.reset();
		EFieldModel_m.reset();
	}
}

void Simulation::printSimAttributes(int numberOfIterations, int itersBtwCouts) //protected
{
	//Sim Header (folder) printed from Python - move here eventually
	std::cout << "Sim between:    " << simMin_m << "m - " << simMax_m << "m" << std::endl;
	std::cout << "dt:             " << dt_m << "s" << std::endl;
	std::cout << "BField Model:   " << BFieldModel_m->getName() << std::endl;
	std::cout << "EField Elems:   " << ((EFieldModel_m == nullptr) ? ("") : (EFieldModel_m->getEElemsStr())) << std::endl;
	std::cout << "Particles:      "; // << particles_m.at(0)->getName() << ": #: " << particles_m.at(0)->getNumberOfParticles() << ", loaded files?: " << (particles_m.at(0)->getInitDataLoaded() ? "true" : "false") << std::endl;
	for (int iii = 0; iii < particles_m.size(); iii++) {
		std::cout << ((iii != 0) ? "                " : "") << particles_m.at(iii)->name() << ": #: " << particles_m.at(iii)->getNumberOfParticles() << ", loaded files?: " << (particles_m.at(iii)->getInitDataLoaded() ? "true" : "false") << std::endl;
	}
	std::cout << "Satellites:     "; // << satellites_m.at(0)->satellite->getName() << ": alt: " << satellites_m.at(0)->satellite->getAltitude() << " m, upward?: " << (satellites_m.at(0)->satellite->getUpward() ? "true" : "false") << std::endl;
	for (int iii = 0; iii < satellites_m.size(); iii++) {
		std::cout << ((iii != 0) ? "                " : "") << satellites_m.at(iii)->satellite->name() << ": alt: " << satellites_m.at(iii)->satellite->altitude() << " m, upward?: " << (satellites_m.at(iii)->satellite->upward() ? "true" : "false") << std::endl;
	}
	std::cout << "Iterations:     " << numberOfIterations << std::endl;
	std::cout << "Iters Btw Cout: " << itersBtwCouts << std::endl;
	std::cout << "Time to setup:  "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;
	std::cout << "===============================================================" << std::endl;
}