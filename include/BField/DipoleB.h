#ifndef DIPOLEB_BFIELD_H
#define DIPOLEB_BFIELD_H

#include <cmath>
#include <iostream>
#include "BField\BField.h"
#include "physicalconstants.h"

constexpr double B0{ 3.12e-5 }; //won't change from sim to sim

class DipoleB : public BField
{
protected:
	//Field simulation constants
	double L_m{ 0.0 }; //this will be populated with the L constant in the constructor
	double L_norm_m{ 0.0 }; //same
	double s_max_m{ 0.0 }; //same

	//specified variables
	double ILATDegrees_m{ 0.0 };
	double ds_m{ 0.0 };
	double errorTolerance_m{ 0.0 };

	//protected functions
	__host__ virtual    void   setupEnvironment();
	__host__ virtual    void   deleteEnvironment();
	__host__ __device__ double getSAtLambda(const double lambdaDegrees);
	__host__ __device__ double getLambdaAtS(const double s);

public:
	__host__ __device__ DipoleB(double ILATDegrees, double errorTolerance = 1e-4, double ds = RADIUS_EARTH / 1000.0) :
		BField(), ILATDegrees_m{ ILATDegrees }, ds_m{ ds }, errorTolerance_m{ errorTolerance }
	{
		L_m = RADIUS_EARTH / pow(cos(ILATDegrees * PI / 180.0), 2);
		L_norm_m = L_m / RADIUS_EARTH;
		s_max_m = getSAtLambda(ILATDegrees_m);
		modelName_m = "DipoleB";

		#ifndef __CUDA_ARCH__ //host code
		setupEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ virtual ~DipoleB()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ virtual void setds(double ds) { ds_m = ds; }

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t);
	__host__ __device__ virtual double getGradBAtS (const double s, const double t);

	__host__ virtual double getErrTol() { return errorTolerance_m; }
	__host__ virtual double getds() { return ds_m; }
};

#endif