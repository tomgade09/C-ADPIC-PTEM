#ifndef DIPOLEB_BFIELD_H
#define DIPOLEB_BFIELD_H

#include <cmath>
#include <iostream>
#include "BField\BField.h"
#include "physicalconstants.h"

__host__ __device__ double getSAtLambda_DipoleB(double* consts, int arrayLength, double s);
__host__ __device__ double getLambdaAtS_DipoleB(double* consts, int arrayLength, double s);
__host__ __device__ double BFieldAtS_DipoleB(double* consts, int arrayLength, double s, double simtime);

class DipoleB : public BField
{
protected:
	//Field simulation constants
	double B0_m{ 3.12e-5 }; //B0 specified for the dipole model of the geomagnetic field
	double L_m{ 0.0 }; //this will be populated with the L constant in the constructor
	double L_norm_m{ 0.0 }; //same
	double s_max_m{ 0.0 }; //same

	//specified variables
	double ILATDegrees_m{ 0.0 };
	double ds_m{ 0.0 }; //For Grad B, seems to me to be a reasonable ds, but needs a test
	double errorTolerance_m{ 0.0 };

	//protected functions
	virtual void setupEnvironment();

public:
	DipoleB(double ILATDegrees, double errorTolerance = 1e-4, double normFactor = RADIUS_EARTH, double ds = RADIUS_EARTH / 1000) :
		BField(), ILATDegrees_m{ ILATDegrees }, ds_m{ ds }, errorTolerance_m{ errorTolerance }
	{
		L_m = RADIUS_EARTH / pow(cos(ILATDegrees * PI / 180.0), 2);
		L_norm_m = L_m / RADIUS_EARTH;
		
		fieldConstArray_m = { B0_m, ILATDegrees_m, L_m, L_norm_m, 0.0, ds_m, errorTolerance_m };
		s_max_m = getSAtLambda_DipoleB(fieldConstArray_m.data(), fieldConstArray_m.size(), ILATDegrees_m);
		fieldConstArray_m.at(4) = s_max_m;
		
		setupEnvironment();
	}
	~DipoleB() {}

	virtual double getBFieldAtS(double s, double t);
	virtual double getGradBAtS(double s, double t);

	//do you need access functions here?  Return the various constants?  Prob not
};

#endif