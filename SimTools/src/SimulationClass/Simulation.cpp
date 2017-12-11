//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

//protected functions
void Simulation::receiveSatelliteData()
{
	std::cout << "Not removing zeros.\n";
	std::vector<std::vector<double*>> tmpcont; //vector of satellites[attributes]
	tmpcont.reserve(satellites_m.size());
	for (int iii = 0; iii < satellites_m.size(); iii++)
	{
		satellites_m[iii]->copyDataToHost(false);
		std::vector<double*> tmp; //vector of attributes[individual particles (through double*)]
		for (int jjj = 0; jjj < numberOfAttributesTracked_m + 1; jjj++)
		{
			double* satDat{ satellites_m[iii]->getDataArrayPointer(jjj) };
			//double* dbltmp = new double[static_cast<int>(satDat[0]) + 1];
			//std::copy(&satDat[0], &satDat[static_cast<int>(satDat[0])], &dbltmp[0]);
			double* dbltmp = new double[numberOfParticlesPerType_m];
			std::copy(&satDat[0], &satDat[numberOfParticlesPerType_m], &dbltmp[0]);
			tmp.push_back(dbltmp);
		}

		//convert mu to vperp
		//for (int jjj = 1; jjj < static_cast<int>(tmp[0][0]) + 1; jjj++)
		for (int jjj = 1; jjj < numberOfParticlesPerType_m; jjj++)
			tmp[1][jjj] = sqrt(2 * tmp[1][jjj] * BFieldatZ(tmp[2][jjj], simTime_m) / ((satellites_m[iii]->getElecTF()) ? (MASS_ELECTRON) : (MASS_PROTON)));

		tmpcont.push_back(tmp);
	}
	std::cout << "Don't forget to change receiveSatelliteData to reflect if whole array is used or just non-zero.\n";
	satelliteData_m.push_back(tmpcont);
}

double*** Simulation::form3Darray()
{
	double*** array3D;
	array3D = new double**[numberOfParticleTypes_m];
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		array3D[iii] = new double*[numberOfAttributesTracked_m];
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
		{
			array3D[iii][jjj] = new double[numberOfParticlesPerType_m];
			for (int kk = 0; kk < numberOfParticlesPerType_m; kk++)//values to initialize array
				array3D[iii][jjj][kk] = 0.0;
		}
	}
	return array3D;
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
	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfParticlesPerType_m, particlesorig_m[iii][vind][jjj] = sqrt(2 * particlesorig_m[iii][vind][jjj] * BFieldatZ(particlesorig_m[iii][zind][jjj], simTime_m) / mass_m[iii]);)
	
	mu_m = false;
}