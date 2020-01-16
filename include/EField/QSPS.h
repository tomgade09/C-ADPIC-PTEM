#ifndef QSPS_EFIELD_H
#define QSPS_EFIELD_H

#include "EField/EField.h"

using std::vector;

class QSPS : public EElem
{
protected:
	#ifndef __CUDA_ARCH__ //host code
	vector<meters> altMin_m;
	vector<meters> altMax_m;
	vector<double> magnitude_m;
	#endif /* !__CUDA_ARCH__ */
	
	int numRegions_m{ 0 };

	meters* altMin_d; //on host this stores the pointer to the data on GPU, on GPU ditto
	meters* altMax_d;
	double* magnitude_d;

	bool useGPU_m{ true };

	__host__ void setupEnvironment()  override;
	__host__ void deleteEnvironment() override;
	__host__ void deserialize(string serialFolder, int nameIndex) override;

public:
	__host__ QSPS(vector<meters> altMin, vector<meters> altMax, vector<double> magnitude);
	__device__ QSPS(meters* altMin, meters* altMax, double* magnitude, int numRegions);
	__host__ __device__ ~QSPS();
	__host__ __device__ QSPS(const QSPS&) = delete;
	__host__ __device__ QSPS& operator=(const QSPS&) = delete;

	__host__ __device__ Vperm getEFieldAtS(const meters s, const seconds t) const override;

	__host__ const vector<meters>& altMin() const;
	__host__ const vector<meters>& altMax() const;
	__host__ const vector<double>& magnitude() const;

	__host__ void serialize(string serialFolder) const override;
};

#endif /* !QSPS_EFIELD_H */