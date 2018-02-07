#ifndef QSPS_EFIELD_H
#define QSPS_EFIELD_H

#include "EField\EField.h"
#include "StandaloneTools\binaryfiletools.h"

//class QSPS;
//__global__ void deleteEnvironmentGPU_QSPS(EElem** qsps);

class QSPS : public EElem
{
protected:
	#ifndef __CUDA_ARCH__ //host code
	std::vector<double> altMinMax_m;
	std::vector<double> magnitude_m;
	#endif /* !__CUDA_ARCH__ */
	
	int numRegions_m{ 0 };

	double* altMinMax_d; //on host this stores the pointer to the data on GPU, on GPU ditto
	double* magnitude_d;

	__host__ void setupEnvironment();
	__host__ void deleteEnvironment();

public:
	#ifndef __CUDA_ARCH__ //host code
	__host__ QSPS(std::string altMinMaxStr, std::string magStr): EElem()
	{
		modelName_m = "QSPS";
		altMinMax_m = constCharToDblVec(altMinMaxStr.c_str());
		magnitude_m = constCharToDblVec(magStr.c_str());

		if (magnitude_m.size() != altMinMax_m.size() / 2)
			throw std::invalid_argument ("QSPS::QSPS: invalid parameters passed in magStr, altMinMaxStr: resolved vector lengths are not equal");

		setupEnvironment();
	}

	__host__ QSPS(std::vector<double> altMinMax, std::vector<double> magnitude) : EElem(), altMinMax_m{ altMinMax }, magnitude_m{ magnitude }
	{
		modelName_m = "QSPS";
		if (magnitude_m.size() != altMinMax_m.size() / 2)
			throw std::invalid_argument("QSPS::QSPS: invalid parameters passed in magStr, altMinMaxStr: resolved vector lengths are not equal");

		setupEnvironment();
	}
	#endif /* !__CUDA_ARCH__ */

	__device__ QSPS(double* altMinMax, double* magnitude, int numRegions) :
		EElem(), altMinMax_d{ altMinMax }, magnitude_d{ magnitude }, numRegions_m{ numRegions } {}

	__host__ __device__ virtual ~QSPS()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__ __device__ double getEFieldAtS(const double s, const double t);
};

#endif /* !QSPS_EFIELD_H */