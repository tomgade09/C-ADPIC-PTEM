//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

//protected functions
void Simulation::receiveSatelliteData()
{
	std::vector<std::vector<std::vector<double>>> sats; //vector of satellites[attributes]

	for (int satInd = 0; satInd < satellites_m.size(); satInd++)
	{
		satellites_m.at(satInd)->copyDataToHost();
		sats.push_back(satellites_m.at(satInd)->getDataArray(true));

		//convert mu to vperp
		bool flag{ 0 };
		for (int partInd = 0; partInd < particleTypes_m.size(); partInd++)
		{
			sats.at(satInd).at(1).at(partInd) = sqrt(2 * sats.at(satInd).at(1).at(partInd) * BFieldatZ(sats.at(satInd).at(2).at(partInd), simTime_m) / ((satellites_m.at(satInd)->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)));
			if ((!flag) && (isnan(sats.at(satInd).at(1).at(partInd))))
			{
				flag = true;
				std::cout << "IsNaN!!  First NaN index: " << partInd << "\nSatellite: " << satellites_m.at(satInd)->getName() << "\nz, BField: " << sats.at(satInd).at(2).at(partInd) << ", " << BFieldatZ(sats.at(satInd).at(2).at(partInd), simTime_m);
				std::cout << "\nElecTF: " << satellites_m.at(satInd)->getElecTF() << ", Mass: " << ((satellites_m.at(satInd)->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)) << "\n";
			}
		}
	}

	satelliteData_m.push_back(sats);
}

//public functions
void Simulation::createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor)
{
	Particle* newPart = new Particle(name, attrNames, mass, charge, numParts, posDims, velDims, normFactor);
	particleTypes_m.push_back(newPart);
}

void Simulation::createSatellite(int partInd, double altitude, bool upwardFacing, double** GPUdataPointers, bool elecTF, std::string name)
{//remove elecTF, change to struct
	Particle* tmp{ particleTypes_m.at(partInd) };
	Satellite* newSat = new Satellite(altitude, upwardFacing, tmp->getNumberOfAttributes(), tmp->getNumberOfParticles(), GPUdataPointers, elecTF, name);
	satellites_m.push_back(newSat);
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
		std::cout << "Warning: satelliteData size is 0.  This probably means Simulation::receiveSatelliteData hasn't been called yet.  Returning.\n";
		return;
	}

	for (int hhh = 0; hhh < particleTypes_m.size(); hhh++)
	{
		Particle* tmpPart{ particleTypes_m.at(hhh) };
		long numParts{ tmpPart->getNumberOfParticles() };

		std::ofstream csv(filename.at(hhh), std::ios::trunc);
		csv << "v_para orig,v_perp orig,z orig,,time escaped top,para top,perp top,z top,,time escaped bottom,para bottom,perp bottom,z bottom,,Energy (eV), Pitch Angle\n";
		csv.close();

		std::vector<std::vector<double>> data;
		std::vector<double> zeros;
		zeros.resize(numParts);

		LOOP_OVER_1D_ARRAY(3, data.push_back(tmpPart->getOrigData().at(iii));) //orig para, perp, z
		data.push_back(zeros); //spacer
		data.push_back(satelliteData_m.at(0).at(hhh + 2).at(3)); //time escaped top
		LOOP_OVER_1D_ARRAY(3, data.push_back(satelliteData_m.at(0).at(hhh + 2).at(iii));) //top para, perp, z
		data.push_back(zeros);
		data.push_back(satelliteData_m.at(0).at(hhh).at(3)); //time escaped bottom
		LOOP_OVER_1D_ARRAY(3, data.push_back(satelliteData_m.at(0).at(hhh).at(iii));) //bottom para, perp, z
		data.push_back(zeros);

		std::vector<double> tmp;
		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(0.5 * tmpPart->getMass() * ((pow(tmpPart->getOrigData().at(0).at(iii), 2) + pow(tmpPart->getOrigData().at(1).at(iii), 2)) * pow(RADIUS_EARTH, 2)) / 1.60218e-19);)
		data.push_back(tmp); //Energies in eV
		tmp.clear();

		LOOP_OVER_1D_ARRAY(numParts, tmp.push_back(atan2(abs(tmpPart->getOrigData().at(1).at(iii)), -tmpPart->getOrigData().at(0).at(iii)) * 180 / PI);)
		data.push_back(tmp);

		fileIO::write2DCSV(data, filename.at(hhh), numParts, 16, ',', false);
	}
}

double* Simulation::getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex, bool originalData)
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

void Simulation::prepareResults()
{
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), convertMuToVPerp(particleTypes_m.at(iii));)

		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_init/", true);)
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->saveArrayToFiles("./bins/particles_final/", false);)

		//normalizes m to Re
		if (normalizedToRe_m)
		{
			LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->normalizeParticles(true, true);)

				for (int lll = 0; lll < satelliteData_m.size(); lll++)
				{//loop over number of measurements
					LOOP_OVER_1D_ARRAY(satellites_m.size(), convertMuToVPerp(satelliteData_m.at(lll).at(iii).at(1), satelliteData_m.at(lll).at(iii).at(2), particleTypes_m.at(iii % 2)->getMass());)
						LOOP_OVER_3D_ARRAY(satellites_m.size(), satellites_m.at(iii)->getNumOfAttr(), satellites_m.at(iii)->getNumOfParts(), satelliteData_m.at(lll).at(iii).at(jjj).at(kk) /= RADIUS_EARTH;)
				}
		}

	resultsPrepared_m = true;
}