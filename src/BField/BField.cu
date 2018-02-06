#include "BField\BField.h"

//constexpr double B0ATTHETA{ BFIELD_EARTH *  1.9102530 };//sqrt(1 + 3 * pow(cos(20.0 * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))
/*__host__ __device__ double BFieldatZ(double z, double simtime)
{//for now, a simple dipole field
	if (z == 0)
		return 0.0; //add an error here if this case is true, at some point

	double norm{ RADIUS_EARTH };

	if ((z < RADIUS_EARTH) && (z > 0))
		norm = 1.0;

	return B0ATTHETA / pow(z / norm, 3); //Bz = B0 at theta * (1/rz(in Re))^3
}*/

//F_mir = -args[2] * B0ATTHETA * (-3 / (pow(ztmp / RADIUS_EARTH, 4))) / RADIUS_EARTH; //mu in [kg.m^2 / s^2.T] = [N.m / T]

/*__host__ __device__ double qspsEatZ(double s, double simtime, double constE)
{
	//if ((s > E_RNG_CENTER + E_RNG_DELTA) || (s < E_RNG_CENTER - E_RNG_DELTA))
		//return 0.0;
	return constE;
}*/

/*__host__ __device__ double EFieldatZ(double** LUT, double s, double simtime, double omegaE, double constE, bool qsps, bool alfven)
{
	bool alfLUT{ false };
	bool alfCalc{ false };

	if (LUT == nullptr && alfven)
		alfCalc = true;
	else if (LUT != nullptr && alfven)
		alfLUT = true;

	return 0.0; //(qsps ? (qspsEatZ(s, simtime, constE)) : (0.0)) + (alfLUT ? (alfvenWaveEbyLUT(LUT, s, simtime, omegaE)) : (0.0)) + (alfCalc ? (alfvenWaveEbyCompute(s, simtime)) : (0.0));
}*/