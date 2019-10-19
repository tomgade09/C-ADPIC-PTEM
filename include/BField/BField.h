#ifndef BFIELD_H
#define BFIELD_H

#include <string>

//CUDA includes
#include "cuda_runtime.h"

class BField
{
protected:
	BField** this_d{ nullptr };
	
	const char* modelName_m;

	__host__ virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__ virtual void deleteEnvironment() = 0;

	__host__ __device__ BField(const char* modelName) : modelName_m{ modelName } {}

public:
	__host__ __device__ virtual ~BField() {}
	__host__ __device__ BField(const BField&) = delete;
	__host__ __device__ BField& operator=(const BField&) = delete;

	__host__ __device__ virtual double getBFieldAtS(const double s, const double t) const = 0;
	__host__ __device__ virtual double getGradBAtS (const double s, const double t) const = 0;
	__host__ __device__ virtual double getSAtAlt(const double alt_fromRe) const = 0;

	virtual double ILAT() const = 0;

	__host__ virtual std::string name()   const { return modelName_m; }
	__host__ virtual BField** getPtrGPU() const { return this_d; } //once returned, have to cast it to the appropriate type
};

#endif