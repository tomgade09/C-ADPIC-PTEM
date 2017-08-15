#include <iostream>
#include <cmath>
#include <iostream>

#include "include\_simulationvariables.h"

//this file will get bigger once the simulations become more complex

double accel1DCB(double* args, int len) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [dt, vz, mu, q, m, pz_0]
	if (len != 6)
	{
		std::cout << "Array is not the right length.  Proper array format is: [dt, vz, mu, q, m, pz_0, E].  Returning zero.\n";
		return 0.0;
	}

	double F_lor, F_mir;
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * CONSTEFIELD / NORMFACTOR; //will need to replace E with a function to calculate in more complex models
	
	//Mirror force
	F_mir = -args[2] * (-3 / pow(args[5],4)) * DIPOLECONST; //have function for gradB based on dipole B field - will need to change later

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field