//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

void Simulation::receiveSatelliteData(bool removeZeros)
{
	LOOP_OVER_1D_ARRAY(satellites_m.size(), satellites_m.at(iii)->satellite->copyDataToHost(););
	LOOP_OVER_1D_ARRAY(satellites_m.size(), satelliteData_m.push_back(satellites_m.at(iii)->satellite->getConsolidatedData(removeZeros)););

	//Check particle for index of vperp/mu, iterate over particles
	LOOP_OVER_1D_ARRAY(satellites_m.size(),\
		int vperpInd{ satellites_m.at(iii)->particle->getDimensionIndByName("vperp") };
		int sInd    { satellites_m.at(iii)->particle->getDimensionIndByName("s") };
		int tInd    { satellites_m.at(iii)->satellite->getNumberOfAttributes() };
		convertMuToVPerp(satelliteData_m.at(iii).at(vperpInd), satelliteData_m.at(iii).at(sInd), satelliteData_m.at(iii).at(tInd), satellites_m.at(iii)->particle->getMass());
	);

	LOOP_OVER_2D_ARRAY(satellites_m.size(), satellites_m.at(iii)->satellite->getNumberOfAttributes() + 2,\
		std::string name{ "./bins/satellites/" };
		name += satellites_m.at(iii)->satellite->getName() + "_";
		if (jjj == satellites_m.at(iii)->satellite->getNumberOfAttributes())
			name += "time";
		else if (jjj == satellites_m.at(iii)->satellite->getNumberOfAttributes() + 1)
			name += "index";
		else
			name += satellites_m.at(iii)->particle->getDimensionNameByInd(jjj);
		name += ".bin";
		fileIO::writeDblBin(satelliteData_m.at(iii).at(jjj), name, satellites_m.at(iii)->particle->getNumberOfParticles());
	);
}

//public functions
void Simulation::createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir)
{
	//some sort of debug message??  Necessary for daily use?
	logFile_m.writeLogFileEntry("Simulation::createParticleType: Particle Type Created: " + name + ": Mass: " + std::to_string(mass) + ", Charge: " + std::to_string(charge) + ", Number of Parts: " + std::to_string(numParts) + ", Pos Dimensions: " + std::to_string(posDims) + ", Vel Dimensions: " + std::to_string(velDims) + ", Files Loaded?: " + ((loadFilesDir == "") ? "False" : "True"));

	Particle* newPart = new Particle(name, attrNames, mass, charge, numParts, posDims, velDims, normFactor);
	particleTypes_m.push_back(newPart);

	if (loadFilesDir != "")
		newPart->loadFilesToArray(loadFilesDir);
}

void Simulation::createSatellite(int partInd, double altitude, bool upwardFacing, std::string name)
{//remove elecTF, change to struct
	if (particleTypes_m.size() <= partInd)
	{
		logFile_m.writeErrorEntry("Simulation::createSatellite", "particleTypes.at(" + std::to_string(partInd) + ") does not exist!", { std::to_string(partInd), std::to_string(altitude), std::to_string(upwardFacing), name });
		errorEncountered = true;
	}
	if (particleTypes_m.at(partInd)->getCurrDataGPUPtr() == nullptr)
	{
		logFile_m.writeErrorEntry("Simulation::createSatellite", "Pointer to GPU data is a nullptr.  That's just asking for trouble.", { std::to_string(partInd), std::to_string(altitude), std::to_string(upwardFacing), name });
		errorEncountered = true;
	}

	logFile_m.writeLogFileEntry("Simulation::createSatellite: Created Satellite: " + name + ", Particle tracked: " + particleTypes_m.at(partInd)->getName() + ", Altitude: " + std::to_string(altitude) + ", " + ((upwardFacing) ? "Upward" : "Downward") + " Facing Detector");

	Particle* tmpPart{ particleTypes_m.at(partInd) };
	Satellite* newSat = new Satellite(altitude, upwardFacing, tmpPart->getNumberOfAttributes(), tmpPart->getNumberOfParticles(), tmpPart->getCurrDataGPUPtr(), name);
	SatandPart* newStruct = new SatandPart{ newSat, tmpPart };
	satellites_m.push_back(newStruct);
}

//Vperp <-> Mu conversion tools
void Simulation::convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, double mass)
{
	LOOP_OVER_1D_ARRAY(vperp.size(), vperp.at(iii) = 0.5 * mass * vperp.at(iii) * vperp.at(iii) / BFieldModel->getBFieldAtS(s.at(iii), simTime_m););
}

void Simulation::convertVPerpToMu(Particle* particle)
{
	convertVPerpToMu(particle->getCurrData().at(particle->getDimensionIndByName("vperp")), particle->getCurrData().at(particle->getDimensionIndByName("s")), particle->getMass());
}

void Simulation::convertVPerpToMu(int partInd)
{
	if (partInd > (particleTypes_m.size() - 1))
	{
		logFile_m.writeErrorEntry("Simulation::convertVPerpToMu", "Index out of range.  Returning without change.", { std::to_string(partInd) });
		return;
	}

	convertVPerpToMu(particleTypes_m.at(partInd));
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, double mass)
{
	LOOP_OVER_1D_ARRAY(mu.size(), mu.at(iii) = sqrt(2 * mu.at(iii) * BFieldModel->getBFieldAtS(s.at(iii), simTime_m) / mass););
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, std::vector<double>& t, double mass)
{
	LOOP_OVER_1D_ARRAY(mu.size(), mu.at(iii) = sqrt(2 * mu.at(iii) * BFieldModel->getBFieldAtS(s.at(iii), t.at(iii)) / mass););
}

void Simulation::convertMuToVPerp(Particle* particle)
{
	convertMuToVPerp(particle->getCurrData().at(particle->getDimensionIndByName("vperp")), particle->getCurrData().at(particle->getDimensionIndByName("s")), particle->getMass());
}

void Simulation::convertMuToVPerp(int partInd)
{
	if (partInd > (particleTypes_m.size() - 1))
	{
		logFile_m.writeErrorEntry("Simulation::convertMuToVPerp", "Index out of range.  Returning without change.", { std::to_string(partInd) });
		return;
	}

	convertMuToVPerp(particleTypes_m.at(partInd));
}

void Simulation::writeSatelliteDataToCSV()
{//need to make this more generic
	std::vector<std::string> filename{ "./elecoutput.csv", "./ionsoutput.csv" };

	if (satelliteData_m.size() == 0)
	{
		logFile_m.writeErrorEntry("Simulation::writeSatelliteDataToCSV", "satelliteData size is 0.  This probably means Simulation::receiveSatelliteData hasn't been called yet.  Returning.", {});
		return;
	}
	//need to index the satellite data with the orig data by the saved index number
	for (int hhh = 0; hhh < particleTypes_m.size(); hhh++)
	{
		Particle* tmpPart{ particleTypes_m.at(hhh) };
		int numAttrs{ tmpPart->getNumberOfAttributes() };
		long numParts{ tmpPart->getNumberOfParticles() };

		std::ofstream csv(filename.at(hhh), std::ios::trunc);
		csv << "v_para orig,v_perp orig,z orig,,time escaped top,para top,perp top,z top,,time escaped bottom,para bottom,perp bottom,z bottom,,Energy (eV), Pitch Angle\n";
		csv.close();

		std::vector<std::vector<double>> data;
		std::vector<double> zeros;
		zeros.resize(numParts);

		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(tmpPart->getOrigData().at(iii));); //orig para, perp, s
		data.push_back(zeros); //spacer
		data.push_back(satelliteData_m.at(hhh + particleTypes_m.size()).at(numAttrs)); //time escaped top
		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(satelliteData_m.at(hhh + particleTypes_m.size()).at(iii));); //top para, perp, s
		data.push_back(zeros);
		data.push_back(satelliteData_m.at(hhh).at(numAttrs)); //time escaped bottom
		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(satelliteData_m.at(hhh).at(iii));); //bottom para, perp, s
		data.push_back(zeros);

		int vparaInd{ tmpPart->getDimensionIndByName("vpara") };
		int vperpInd{ tmpPart->getDimensionIndByName("vperp") };
		std::vector<double> tmp;
		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(0.5 * tmpPart->getMass() * ((pow(tmpPart->getOrigData().at(vparaInd).at(iii), 2) + pow(tmpPart->getOrigData().at(vperpInd).at(iii), 2)) * pow(RADIUS_EARTH, 2)) / 1.60218e-19););
		data.push_back(tmp); //Energies in eV
		tmp.clear();

		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(atan2(abs(tmpPart->getOrigData().at(vperpInd).at(iii)), -tmpPart->getOrigData().at(vparaInd).at(iii)) * 180 / PI););
		data.push_back(tmp);

		fileIO::write2DCSV(data, filename.at(hhh), numParts, numAttrs * 3 + 7, ',', false);
	}
}

double* Simulation::getPointerToParticleAttributeArray(int partIndex, int attrIndex, bool originalData)
{
	if (partIndex > (particleTypes_m.size() - 1))
	{
		logFile_m.writeErrorEntry("Simulation::getPointerToParticleAttributeArray", "Particle Index out of bounds.  There aren't that many particle types.  Returning nullptr.", { std::to_string(partIndex), std::to_string(attrIndex), std::to_string(originalData) });
		return nullptr;
	}
	else if (attrIndex > (particleTypes_m.at(partIndex)->getNumberOfAttributes() - 1))
	{
		logFile_m.writeErrorEntry("Simulation::getPointerToParticleAttributeArray", "Attribute Index out of bounds.  There aren't that many attributes for the particle type.  Returning nullptr.", { std::to_string(partIndex), std::to_string(attrIndex), std::to_string(originalData) });
		return nullptr;
	}

	return ((originalData) ? (particleTypes_m.at(partIndex)->getOrigData().at(attrIndex).data()) : (particleTypes_m.at(partIndex)->getCurrData().at(attrIndex).data()));
}

void Simulation::prepareResults(bool normalizeToRe)
{
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), convertMuToVPerp(particleTypes_m.at(iii)););

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_init/", true););
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_final/", false););

	//normalizes m to Re
	if (normalizeToRe)
	{
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->normalizeParticles(true, true););
		LOOP_OVER_2D_ARRAY(satellites_m.size(), satellites_m.at(iii)->satellite->getNumberOfAttributes(), normalizeArray(satelliteData_m.at(iii).at(jjj), RADIUS_EARTH););
	}

	resultsPrepared_m = true;
}

void Simulation::loadCompletedSimData(std::string fileDir, std::vector<std::string> partNames, std::vector<std::string> attrNames, std::vector<std::string> satNames, int numParts)
{
	for (size_t parts = 0; parts < partNames.size(); parts++)
	{
		createParticleType(partNames.at(parts), attrNames, 1, 1, numParts, static_cast<int>(attrNames.size() - 1), 1, 1, fileDir + "particles_final/");
		particleTypes_m.at(parts)->loadFilesToArray(fileDir + "particles_init/", true);
	}

	for (size_t sats = 0; sats < satNames.size(); sats++)
		createSatellite(0, 1, true, satNames.at(sats));

	attrNames.push_back("time");
	attrNames.push_back("index");

	std::vector<std::vector<std::vector<double>>> tmp3D;
	for (size_t sats = 0; sats < satNames.size(); sats++)
	{
		std::vector<std::vector<double>> tmp2D;
		for (size_t attrs = 0; attrs < attrNames.size(); attrs++)
		{
			std::vector<double> tmp;
			tmp.resize(numParts);
			fileIO::readDblBin(tmp, fileDir + "satellites/" + satNames.at(sats) + "_" + attrNames.at(attrs) + ".bin", numParts);
			tmp2D.push_back(tmp);
		}
		tmp3D.push_back(tmp2D);
	}

	satelliteData_m = tmp3D;
}