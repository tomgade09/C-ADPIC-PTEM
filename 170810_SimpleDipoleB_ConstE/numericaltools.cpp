#include <random>
#include <cmath>

#include "include\numericaltools.h"
#include "include\_simulationvariables.h"

double fourthOrderRungeKutta1D(FncPnt1d_t funcPointer, double* funcArg, int arrayLen, double h)
{	// funcArg requirements: [dt, y_0, ...] where dt = {0, h/2, h}, initial dt should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = funcPointer(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = h / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = funcPointer(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = funcPointer(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = h;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = funcPointer(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns units of y, not dy / dt
}

double** normalDistribution_v_z(int numOfParticles, double vmean, double vsigma, double zmean, double zsigma)
{
	std::random_device randDev;
	std::mt19937 mtgen(randDev());

	std::normal_distribution<> vpara_nd(vmean, VPARACONST * vsigma);
	std::normal_distribution<> vperp_nd(vmean, vsigma);

	double* vpara = new double[numOfParticles];
	double* vperp = new double[numOfParticles];
	double* z = new double[numOfParticles];
	double** vpara_vperp_z = new double*[4];

	for (int iii = 0; iii < numOfParticles; iii++)
	{
		vpara[iii] = vpara_nd(mtgen);
		vperp[iii] = vperp_nd(mtgen);
		if (iii % 2 == 0)
		{ //alternates placing the electron, ion at the Ionosphere / Magnetosphere and setting vpara in the direction away from the boundary
			z[iii] = IONSPH_MIN_Z + 0.1;
			vpara[iii] = abs(vpara[iii]);
		}
		else
		{
			z[iii] = MAGSPH_MAX_Z - 0.1;
			vpara[iii] = -abs(vpara[iii]);
		}
	}

	vpara_vperp_z[0] = vpara;
	vpara_vperp_z[1] = vperp;
	vpara_vperp_z[2] = z;
	vpara_vperp_z[3] = nullptr;

	return vpara_vperp_z;
}