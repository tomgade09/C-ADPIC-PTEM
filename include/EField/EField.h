#ifndef EFIELD_H
#define EFIELD_H

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

typedef double(*callbackFcn)(const double, const double);

//on GPU global variables
extern __device__ callbackFcn EFieldFcnPtr_GPU;

class EField
{
protected:
	void setupArrayOfEFieldElements(int numOfElements);

public:
	EField() {}
	~EField() {}

};

#endif /* EFIELD_H */