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
	__host__ virtual void setupEnvironment();
	__host__ virtual void deleteEnvironment();

public:
	__host__ __device__ DipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_gradB, int numberOfMeasurements) :
		BField(), ILATDegrees_m{ ILATDegrees }, simMin_m{ simMin }, simMax_m{ simMax }, ds_gradB_m{ ds_gradB }, numMsmts_m{ numberOfMeasurements }
	{
		ds_msmt_m = (simMax_m - simMin_m) / (numMsmts_m - 1);
		modelName_m = "DipoleBLUT";
		
		#ifndef __CUDA_ARCH__ //host code
		setupEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ virtual ~DipoleBLUT()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t);
	__host__ __device__ virtual double getGradBAtS(const double s, const double t);

	__device__ void setAltArray(double* altArray) { altitude_d = altArray; }
	__device__ void setMagArray(double* magArray) { magnitude_d = magArray; }

	__host__ virtual double getErrTol() { return 1.0e-10; }
	__host__ virtual double getds() { return ds_gradB_m; }
};

#endif /* !DIPOLEBLUT_BFIELD_H */