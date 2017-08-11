#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "include\_simulationvariables.h"

double*** mainCUDA(double** electrons, double** ions)
{
	//allocate arrays on device
}

__global__ double accel1dCUDA(double* args, int len) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [dt, vz, mu, q, m, pz_0]
	if (len != 6)
	{
		std::cout << "Array is not the right length.  Proper array format is: [dt, vz, mu, q, m, pz_0, E].  Returning zero.\n";
		return 0.0;
	}

	double F_lor, F_mir;
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * CONSTEFIELD / RADIUS_EARTH; //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * (-3 / pow(args[5], 4)) * DIPOLECONST; //have function for gradB based on dipole B field - will need to change later

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__global__ double foRungeKuttaCUDA(double* funcArg, int arrayLen, double h)
{
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = accel1dCUDA(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = h / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = h;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns units of y, not dy / dt
}