#ifndef QSPS_EFIELD_H
#define QSPS_EFIELD_H

#include "EField\EField.h"

//class QSPS;
//__global__ void deleteEnvironmentGPU_QSPS(EElem** qsps);

class QSPS : public EElem
{
protected:
	#ifndef __CUDA_ARCH__ //host code
	std::vector<double> altMin_m;
	std::vector<double> altMax_m;
	std::vector<double> magnitude_m;
	#endif /* !__CUDA_ARCH__ */
	
	int numRegions_m{ 0 };

	double* altMin_d; //on host this stores the pointer to the data on GPU, on GPU ditto
	double* altMax_d;
	double* magnitude_d;

	__host__ void setupEnvironment();
	__host__ void deleteEnvironment();

public:
	__host__ QSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude);

	__device__ QSPS(double* altMin, double* altMax, double* magnitude, int numRegions) :
		EElem(), altMin_d{ altMin }, altMax_d{ altMax }, magnitude_d{ magnitude }, numRegions_m{ numRegions } {}

	__host__ __device__ virtual ~QSPS();

	__host__ __device__ double getEFieldAtS(const double s, const double t);

	__host__ const std::vector<double>& altMin();
	__host__ const std::vector<double>& altMax();
	__host__ const std::vector<double>& magnitude();
};

#endif /* !QSPS_EFIELD_H */