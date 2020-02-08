#ifndef BFIELD_H
#define BFIELD_H

#include <string>

//CUDA includes
#include "cuda_runtime.h"

#include "utils/unitsTypedefs.h"

using std::string;

class BField
{
protected:
	BField** this_d{ nullptr }; //pointer to device-side instance
	
	const char* name_m{ nullptr };

	__host__            virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__            virtual void deleteEnvironment() = 0;
	__host__            virtual void deserialize(string serialFolder) = 0;

	__host__ __device__ BField(const char* modelName);

public:
	__host__ __device__ virtual ~BField();
	__host__ __device__ BField(const BField&) = delete;
	__host__ __device__ BField& operator=(const BField&) = delete;

	__host__ __device__ virtual tesla  getBFieldAtS(const meters s, const seconds t) const = 0;
	__host__ __device__ virtual double getGradBAtS (const meters s, const seconds t) const = 0;
	__host__ __device__ virtual meters getSAtAlt(const meters alt_fromRe) const = 0;

	__host__            virtual meters ILAT() const = 0;

	__host__            virtual string name()   const;
	__host__            virtual BField** getPtrGPU() const; //once returned, have to cast it to the appropriate type

	__host__            virtual void serialize(string serialFolder) const = 0;
};

#endif