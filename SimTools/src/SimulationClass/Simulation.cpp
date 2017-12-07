//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "SimulationClass\Simulation.h"

void Simulation::createSatellite(double altitude, bool upwardFacing, double** GPUdataPointers, std::string name)
{
	Satellite* newSat = new Satellite(altitude, upwardFacing, numberOfAttributesTracked_m, numberOfParticlesPerType_m, GPUdataPointers, name);
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