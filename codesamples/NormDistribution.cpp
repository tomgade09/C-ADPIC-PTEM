#include <random>
#include <cmath>
#include <iostream>

#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

constexpr double C_PI = 3.141592653589793238463;

DLLEXPORT double* normalDistribution_v_z(int numOfParticles, double vmean, double vsigma, double zmean, double zsigma)
{
	std::random_device randDev;
	std::mt19937 mtgen(randDev());

	std::normal_distribution<> vpara_nd(vmean, vsigma);
	std::normal_distribution<> vperp_nd(vmean, vsigma);
	std::normal_distribution<> z_nd(zmean, zsigma);

	double* vpara_vperp_z_pitch = new double[numOfParticles * 4];

	for (int iii = 0; iii < numOfParticles; iii++)
	{
		vpara_vperp_z_pitch[iii * 4] = vpara_nd(mtgen);
		vpara_vperp_z_pitch[iii * 4 + 1] = vperp_nd(mtgen);
		vpara_vperp_z_pitch[iii * 4 + 2] = z_nd(mtgen);
		vpara_vperp_z_pitch[iii * 4 + 3] = atan2(vpara_vperp_z_pitch[iii * 4 + 1], vpara_vperp_z_pitch[iii * 4]) * 180 / C_PI;
		//std::cout << vpara_vperp_z_pitch[iii * 4] << " " << vpara_vperp_z_pitch[iii * 4 + 1] << " " << vpara_vperp_z_pitch[iii * 4 + 2] << " " << vpara_vperp_z_pitch[iii * 4 + 3] << "\n";
	}

	//std::cout << "Done with Cpp.\n";

	return vpara_vperp_z_pitch;
}