#ifndef DIPOLEBLUT_BFIELD_H
#define DIPOLEBLUT_BFIELD_H

#include <vector>
#include "BField/BField.h"
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
	__host__ void setupEnvironment() override;
	__host__ void deleteEnvironment() override;

public:
	__host__ __device__ DipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_gradB, int numberOfMeasurements);
	__host__ __device__ ~DipoleBLUT();
	__host__ __device__ DipoleBLUT(const DipoleBLUT&) = delete;
	__host__ __device__ DipoleBLUT& operator=(const DipoleBLUT&) = delete;

	//for testing
	double ILAT() { return ILATDegrees_m; }
	double ds_msmt() { return ds_msmt_m; }
	double ds_gradB() { return ds_gradB_m; }

	__host__ __device__ double getBFieldAtS(const double s, const double t) const override;
	__host__ __device__ double getGradBAtS(const double s, const double t) const override;

	__device__ void setAltArray(double* altArray) { altitude_d = altArray; }
	__device__ void setMagArray(double* magArray) { magnitude_d = magArray; }

	__host__ double getErrTol() const { return 1.0e-10; }
	__host__ double getds()     const { return ds_gradB_m; }
};

#endif /* !DIPOLEBLUT_BFIELD_H */