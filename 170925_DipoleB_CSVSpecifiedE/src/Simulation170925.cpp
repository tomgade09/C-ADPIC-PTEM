#include "include\Simulation170925.h"

void Simulation170925::setElecMagLUT(const char* filename, int rows, int cols)
{
	 elcFieldLUT_m = fileIO::read2DCSV(filename, rows, cols, ' ');
}

double Simulation170925::calculateBFieldAtZandTime(double z, double time)
{
	return BFieldatZ(z, time);
}

double Simulation170925::calculateEFieldAtZandTime(double z, double time)
{
	//return EFieldatZ(getPointerToElectricFieldData(), z, time);
	return EFieldatZ(z, time);
}

void Simulation170925::convertVPerpToMu()
{
	if (mu_m)
	{
		std::cout << "v_perp has already been converted to mu.  Returning with no changes.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			particles_m[iii][1][jjj] = 0.5 * mass_m[iii] * particles_m[iii][1][jjj] * particles_m[iii][1][jjj] / BFieldatZ(particles_m[iii][2][jjj], simTime_m);
	}
	mu_m = true;
}

void Simulation170925::convertMuToVPerp()
{
	if (!mu_m)
	{
		std::cout << "Quantity is v_perp (not mu).  Returning with no changes.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			particles_m[iii][1][jjj] = sqrt(2 * particles_m[iii][1][jjj] * BFieldatZ(particles_m[iii][2][jjj], simTime_m) / mass_m[iii]);
	}
	
	mu_m = false;
}

void Simulation170925::prepareResults()
{
	convertMuToVPerp();
	serializeParticleArray();
	resultsPrepared_m = true;
}