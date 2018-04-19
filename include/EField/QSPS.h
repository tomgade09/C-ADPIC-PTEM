#ifndef QSPS_EFIELD_H
#define QSPS_EFIELD_H

#include "EField\EField.h"
#include "utils\string.h"

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
	__host__ QSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude) :
		EElem(), numRegions_m{ (int)magnitude.size() }
	{
		if (magnitude.size() != altMin.size() || magnitude.size() != altMax.size())
			throw std::invalid_argument("QSPS::QSPS: invalid parameters passed in magnitude, altMin, altMax: resolved vector lengths are not equal");
		
		#ifndef __CUDA_ARCH__ //host code
		altMin_m = altMin;       //unfortunately this wrapping is necessary
		altMax_m = altMax;       //as the vectors above also have to be wrapped
		magnitude_m = magnitude; //in an ifndef/endif block so this will compile
		modelName_m = "QSPS";
		#endif /* !__CUDA_ARCH__ */

		setupEnvironment();
	}

	__device__ QSPS(double* altMin, double* altMax, double* magnitude, int numRegions) :
		EElem(), altMin_d{ altMin }, altMax_d{ altMax }, magnitude_d{ magnitude }, numRegions_m{ numRegions } {}

	__host__ __device__ virtual ~QSPS()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ double getEFieldAtS(const double s, const double t);

	__host__ const std::vector<double>& altMin();
	__host__ const std::vector<double>& altMax();
	__host__ const std::vector<double>& magnitude();
};

#endif /* !QSPS_EFIELD_H */