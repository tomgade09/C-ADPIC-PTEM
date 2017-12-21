//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

//protected functions
void Simulation::receiveSatelliteData()
{
	std::vector<std::vector<double*>> tmpcont; //vector of satellites[attributes]
	tmpcont.reserve(satellites_m.size());
	for (int iii = 0; iii < satellites_m.size(); iii++)
	{
		satellites_m[iii]->copyDataToHost();
		std::vector<double*> tmp; //vector of attributes[individual particles (through double*)]

		double** satDat{ satellites_m[iii]->getDataArrayPointer(true) };
		for (int jjj = 0; jjj < numberOfAttributesTracked_m + 1; jjj++)
			tmp.push_back(satDat[jjj]);

		//convert mu to vperp
		bool flag{ 0 };
		for (int jjj = 1; jjj < numberOfParticlesPerType_m; jjj++)
		{
			tmp[1][jjj] = sqrt(2 * tmp[1][jjj] * BFieldatZ(tmp[2][jjj], simTime_m) / ((satellites_m[iii]->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)));
			if ((!flag) && (isnan(tmp[1][jjj])))
			{
				flag = true;
				std::cout << "IsNaN!!  First NaN index: " << jjj << "\nSatellite: " << satellites_m[iii]->getName() << "\nz, BField: " << tmp[2][jjj] << ", " << BFieldatZ(tmp[2][jjj], simTime_m);
				std::cout << "\nElecTF: " << satellites_m[iii]->getElecTF() << ", Mass: " << ((satellites_m[iii]->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)) << "\n";
			}
		}
		tmpcont.push_back(tmp);
	}
	satelliteData_m.push_back(tmpcont);

	//legacy (read: less efficient) way of doing this
	/*std::vector<std::vector<double*>> tmpcont; //vector of satellites[attributes]
	tmpcont.reserve(satellites_m.size());
	for (int iii = 0; iii < satellites_m.size(); iii++)
	{
		satellites_m[iii]->copyDataToHost();
		std::vector<double*> tmp; //vector of attributes[individual particles (through double*)]
		
		double* satDat{ satellites_m[iii]->getDataArrayPointer() };
		double* dbltmp = new double[LENGTHSATDATA];
		
		std::copy(&satDat[0], &satDat[LENGTHSATDATA - 1], &dbltmp[0]);
		
		for (int jjj = 0; jjj < numberOfAttributesTracked_m + 1; jjj++)
			tmp.push_back(&dbltmp[jjj * numberOfParticlesPerType_m]);

		//convert mu to vperp
		bool flag{ 0 };
		for (int jjj = 1; jjj < numberOfParticlesPerType_m; jjj++)
		{
			tmp[1][jjj] = sqrt(2 * tmp[1][jjj] * BFieldatZ(tmp[2][jjj], simTime_m) / ((satellites_m[iii]->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)));
			if ((!flag) && (isnan(tmp[1][jjj])))
			{
				flag = true;
				std::cout << "IsNaN!!  First NaN index: " << jjj << "\nSatellite: " << satellites_m[iii]->getName() << "\nz, BField: " << tmp[2][jjj] << ", " << BFieldatZ(tmp[2][jjj], simTime_m);
				std::cout << "\nElecTF: " << satellites_m[iii]->getElecTF() << ", Mass: " << ((satellites_m[iii]->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)) << "\n";
			}
		}
		tmpcont.push_back(tmp);
	}
	satelliteData_m.push_back(tmpcont);*/
}

//public functions
void Simulation::createSatellite(double altitude, bool upwardFacing, double** GPUdataPointers, bool elecTF, std::string name)
{
	Satellite* newSat = new Satellite(altitude, upwardFacing, numberOfAttributesTracked_m, numberOfParticlesPerType_m, GPUdataPointers, elecTF, name);
	satellites_m.push_back(newSat);
}

void Simulation::convertVPerpToMu(int vind, int zind)
{
	if (mu_m)
	{
		std::cout << "v_perp has already been converted to mu.  Returning with no changes.\n";
		return;
	}

	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfParticlesPerType_m, particles_m[iii][vind][jjj] = 0.5 * mass_m[iii] * particles_m[iii][vind][jjj] * particles_m[iii][vind][jjj] / BFieldatZ(particles_m[iii][zind][jjj], simTime_m);)
	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfParticlesPerType_m, particlesorig_m[iii][vind][jjj] = 0.5 * mass_m[iii] * particlesorig_m[iii][vind][jjj] * particlesorig_m[iii][vind][jjj] / BFieldatZ(particlesorig_m[iii][zind][jjj], simTime_m);)

	mu_m = true;
}

void Simulation::convertMuToVPerp(int vind, int zind)
{
	if (!mu_m)
	{
		std::cout << "Quantity is v_perp (not mu).  Returning with no changes.\n";
		return;
	}

	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfParticlesPerType_m, particles_m[iii][vind][jjj] = sqrt(2 * particles_m[iii][vind][jjj] * BFieldatZ(particles_m[iii][zind][jjj], simTime_m) / mass_m[iii]);)

	mu_m = false;
}