// 170919_GenerateNormalDistributionFile.cpp
// Generates a normalized distribution of particles, and writes them to a file

#include "include\_simulationvariables.h"
#include "include\numericaltools.h"
#include "include\fileIO.h"

int main()
{
	double** particles;
	particles = normalDistribution_v_z(NUMPARTICLES, V_DIST_MEAN, V_SIGMA, Z_DIST_MEAN, Z_SIGMA);
	
	double* partArray1D = new double[NUMPARTICLES * 3];

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{
		partArray1D[iii + 0 * NUMPARTICLES] = particles[0][iii];
		partArray1D[iii + 1 * NUMPARTICLES] = particles[1][iii];
		partArray1D[iii + 2 * NUMPARTICLES] = particles[2][iii];
	}

	delete[] particles[0];
	delete[] particles[1];
	delete[] particles[2];
	delete[] particles[3];
	delete[] particles;

	dblBinIO::writeDblBin("normDist.bin", partArray1D, NUMPARTICLES * 3);

	delete[] partArray1D;

	return 0;
}

