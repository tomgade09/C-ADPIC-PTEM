#ifndef DIPOLEBLUT_BFIELD_H
#define DIPOLEBLUT_BFIELD_H

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include "BField\DipoleB.h"
#include "physicalconstants.h"

class DipoleBLUT : public BField
{
protected:
	//specified variables
	double ILATDegrees_m{ 0.0 };
	double ds_msmt_m{ 0.0 };
	double ds_gradB_m{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	std::vector<double> altitude_m;
	std::vector<double> magnitude_m;
	#endif /* !__CUDA_ARCH__ */

	//on device variables
	double* altitude_d{ nullptr };
	double* magnitude_d{ nullptr };
	
	double  simMin_m{ 0.0 };
	double  simMax_m{ 0.0 };
	int     numMsmts_m{ 0 };

	//protected functions
	__host__ virtual void setupEnvironment() {} //need to get the compiler to stop whining about setupEnvironment not being overwritten
	__host__ virtual void setupEnvironment(double ds_dipoleB);
	__host__ virtual void deleteEnvironment();

public:
	__host__ __device__ DipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_dipoleB, int numberOfMeasurements) :
		BField(), ILATDegrees_m{ ILATDegrees }, simMin_m{ simMin }, simMax_m{ simMax }, ds_gradB_m{ ds_dipoleB }, numMsmts_m{ numberOfMeasurements }
	{
		#ifndef __CUDA_ARCH__ //host code
		modelName_m = "DipoleBLUT";
		setupEnvironment(ds_dipoleB);
		#endif /* !__CUDA_ARCH__ */
		
		ds_msmt_m = (simMax_m * 1.001 - simMin_m * 0.999) / (numMsmts_m - 1);
	}

	__host__ __device__ virtual ~DipoleBLUT()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t);
	__host__ __device__ virtual double getGradBAtS(const double s, const double t);

	//do you need access functions here?  Return the various constants?  Prob not
};

#endif /* !DIPOLEBLUT_BFIELD_H */