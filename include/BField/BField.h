#ifndef BFIELD_H
#define BFIELD_H

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"

class BField
{
protected:
	BField** this_d{ nullptr };
	
	#ifndef __CUDA_ARCH__ //host code
	std::string modelName_m;
	#else //device code
	const char* modelName_m;
	#endif /* !__CUDA_ARCH__ */

	__host__ virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d

public:
	__host__ __device__ BField() {}
	__host__ __device__ virtual ~BField() {}

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t)=0;
	__host__ __device__ virtual double getGradBAtS (const double s, const double t)=0;

	__host__ virtual std::string getName() { return modelName_m; }
	__host__ virtual BField** getPtrGPU()  { return this_d; } //once returned, have to cast it to the appropriate type
};

#endif