#ifndef QSPS_EFIELD_H
#define QSPS_EFIELD_H

#include "EField\EField.h"

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

	__host__ void setupEnvironment()  override;
	__host__ void deleteEnvironment() override;

public:
	__host__ QSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude);

	__device__ QSPS(double* altMin, double* altMax, double* magnitude, int numRegions) :
		EElem(), altMin_d{ altMin }, altMax_d{ altMax }, magnitude_d{ magnitude }, numRegions_m{ numRegions } {}

	__host__ __device__ ~QSPS();

	__host__ __device__ double getEFieldAtS(const double s, const double t) const override;

	__host__ const std::vector<double>& altMin() const;
	__host__ const std::vector<double>& altMax() const;
	__host__ const std::vector<double>& magnitude() const;
};

#endif /* !QSPS_EFIELD_H */