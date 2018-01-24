#include "BField\BField.h"
#include <cstdlib> //for EXIT_FAILURE
#include <iostream>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)

__host__ __device__ double getBFieldAtS(double(*fcnPtr)(double*, int), double* args, int count)
{
	return fcnPtr(args, count);
}

void BField::setupCallbacksonGPU()
{
	CUDA_CALL(cudaMalloc((void **)&BFieldFcnPtr_d, sizeof(callbackFcn)));
	CUDA_CALL(cudaMalloc((void **)&gradBFcnPtr_d,  sizeof(callbackFcn)));

	callSetupCallbacksKernel();
}