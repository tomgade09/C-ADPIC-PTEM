#ifndef BFIELD_H
#define BFIELD_H

#include <string>

//CUDA includes
#include "host_defines.h"

//Project includes
#include "dlldefines.h"

class BField
{
protected:
	BField** this_d{ nullptr };
	
	#ifndef __CUDA_ARCH__ //host code
	std::string modelName_m;
	#else //device code
	const char* modelName_m; //placeholder, not used
	#endif /* !__CUDA_ARCH__ */

	__host__ virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__ virtual void deleteEnvironment() = 0;

	__host__ __device__ BField() {}

public:
	__host__ __device__ virtual ~BField() {};

	__host__ __device__ BField(const BField&) = delete;
	__host__ __device__ BField& operator=(const BField&) = delete;

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t) const = 0;
	__host__ __device__ virtual double getGradBAtS (const double s, const double t) const = 0;

	__host__ virtual std::string name()   const { return modelName_m; }
	__host__ virtual BField** getPtrGPU() const { return this_d; } //once returned, have to cast it to the appropriate type
};

#endif