//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

void Simulation::receiveSatelliteData(bool removeZeros)
{
	LOOP_OVER_1D_ARRAY(satellites_m.size(), satelliteData_m.push_back(satellites_m.at(iii)->satellite->getConsolidatedData(removeZeros));)
	
	//Check particle for index of vperp/mu, iterate over particles
	//sats.at(satInd).at(1).at(partInd) = sqrt(2 * sats.at(satInd).at(1).at(partInd) * BFieldatZ(sats.at(satInd).at(2).at(partInd), simTime_m) / ((satellites_m.at(satInd)->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)));
	LOOP_OVER_1D_ARRAY(satellites_m.size(),\
		int vperpInd{ satellites_m.at(iii)->particle->getDimensionIndByName("vperp") };
		int zInd    { satellites_m.at(iii)->particle->getDimensionIndByName("z") };
		int tInd    { satellites_m.at(iii)->satellite->getNumberOfAttributes() };
		convertMuToVPerp(satelliteData_m.at(iii).at(vperpInd), satelliteData_m.at(iii).at(zInd), satelliteData_m.at(iii).at(tInd), satellites_m.at(iii)->particle->getMass());
	)
}

//public functions
void Simulation::createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir)
{
	//add something to logfile here
	Particle* newPart = new Particle(name, attrNames, mass, charge, numParts, posDims, velDims, normFactor);
	particleTypes_m.push_back(newPart);

	if (loadFilesDir != "")
		newPart->loadFilesToArray(loadFilesDir);
}

void Simulation::createSatellite(int partInd, double altitude, bool upwardFacing, std::string name)
{//remove elecTF, change to struct
	//add something to logfile here
	Particle* tmpPart{ particleTypes_m.at(partInd) };
	Satellite* newSat = new Satellite(altitude, upwardFacing, tmpPart->getNumberOfAttributes(), tmpPart->getNumberOfParticles(), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(partInd)), name);
	SatandPart* newStruct = new SatandPart{ newSat, tmpPart };
	satellites_m.push_back(newStruct);
}

//Vperp <-> Mu conversion tools
void Simulation::convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& z, double mass)
{
	LOOP_OVER_1D_ARRAY(vperp.size(), vperp.at(iii) = 0.5 * mass * vperp.at(iii) * vperp.at(iii) / BFieldatZ(z.at(iii), simTime_m);)
}

void Simulation::convertVPerpToMu(Particle* particle)
{
	convertVPerpToMu(particle->getCurrData().at(particle->getDimensionIndByName("vperp")), particle->getCurrData().at(particle->getDimensionIndByName("z")), particle->getMass());
}

void Simulation::convertVPerpToMu(int partInd)
{
	if (partInd > (particleTypes_m.size() - 1))
	{
		std::cout << "Error: convertVPerpToMu specified index out of range.  Returning without change.\n";
		return;
	}

	convertVPerpToMu(particleTypes_m.at(partInd));
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, double mass)
{
	LOOP_OVER_1D_ARRAY(mu.size(), mu.at(iii) = sqrt(2 * mu.at(iii) * BFieldatZ(z.at(iii), simTime_m) / mass);)
}

void Simulation::convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, std::vector<double>& t, double mass)
{
	LOOP_OVER_1D_ARRAY(mu.size(), mu.at(iii) = sqrt(2 * mu.at(iii) * BFieldatZ(z.at(iii), t.at(iii)) / mass);)
}

void Simulation::convertMuToVPerp(Particle* particle)
{
	convertMuToVPerp(particle->getCurrData().at(particle->getDimensionIndByName("vperp")), particle->getCurrData().at(particle->getDimensionIndByName("z")), particle->getMass());
}

void Simulation::convertMuToVPerp(int partInd)
{
	if (partInd > (particleTypes_m.size() - 1))
	{
		std::cout << "Error: convertVPerpToMu specified index out of range.  Returning without change.\n";
		return;
	}

	convertMuToVPerp(particleTypes_m.at(partInd));
}

void Simulation::writeSatelliteDataToCSV()
{//need to make this more generic
	std::vector<std::string> filename{ "./elecoutput.csv", "./ionsoutput.csv" };

	if (satelliteData_m.size() == 0)
	{
		std::cout << "Error: satelliteData size is 0.  This probably means Simulation::receiveSatelliteData hasn't been called yet.  Returning.\n";
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

		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(tmpPart->getOrigData().at(iii));) //orig para, perp, z
		data.push_back(zeros); //spacer
		data.push_back(satelliteData_m.at(hhh + particleTypes_m.size()).at(numAttrs)); //time escaped top
		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(satelliteData_m.at(hhh + particleTypes_m.size()).at(iii));) //top para, perp, z
		data.push_back(zeros);
		data.push_back(satelliteData_m.at(hhh).at(numAttrs)); //time escaped bottom
		LOOP_OVER_1D_ARRAY(numAttrs, data.push_back(satelliteData_m.at(hhh).at(iii));) //bottom para, perp, z
		data.push_back(zeros);

		int vparaInd{ tmpPart->getDimensionIndByName("vpara") };
		int vperpInd{ tmpPart->getDimensionIndByName("vperp") };
		std::vector<double> tmp;
		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(0.5 * tmpPart->getMass() * ((pow(tmpPart->getOrigData().at(vparaInd).at(iii), 2) + pow(tmpPart->getOrigData().at(vperpInd).at(iii), 2)) * pow(RADIUS_EARTH, 2)) / 1.60218e-19);)
		data.push_back(tmp); //Energies in eV
		tmp.clear();

		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(atan2(abs(tmpPart->getOrigData().at(vperpInd).at(iii)), -tmpPart->getOrigData().at(vparaInd).at(iii)) * 180 / PI);)
		data.push_back(tmp);

		fileIO::write2DCSV(data, filename.at(hhh), numParts, numAttrs * 3 + 7, ',', false);
	}
}

double* Simulation::getPointerToParticleAttributeArray(int partIndex, int attrIndex, bool originalData)
{
	if (partIndex > (particleTypes_m.size() - 1))
	{
		std::cout << "Particle Index out of bounds.  There aren't that many particle types.\n";
		return nullptr;
	}
	else if (attrIndex > (particleTypes_m.at(partIndex)->getNumberOfAttributes() - 1))
	{
		std::cout << "Attribute Index out of bounds.  There aren't that many attributes for the particle type.\n";
		return nullptr;
	}

	return ((originalData) ? (particleTypes_m.at(partIndex)->getOrigData().at(attrIndex).data()) : (particleTypes_m.at(partIndex)->getCurrData().at(attrIndex).data()));
}

void Simulation::prepareResults(bool normalizeToRe)
{
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), convertMuToVPerp(particleTypes_m.at(iii));)

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_init/", true);)
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_final/", false);)

	LOOP_OVER_1D_ARRAY(satelliteData_m.size(), )

	//normalizes m to Re
	if (normalizeToRe)
	{
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->normalizeParticles(true, true);)
		LOOP_OVER_2D_ARRAY(particleTypes_m.size(), satellites_m.at(iii)->satellite->getNumberOfAttributes(), normalizeArray(satelliteData_m.at(iii).at(jjj), RADIUS_EARTH);)
	}

	resultsPrepared_m = true;
}