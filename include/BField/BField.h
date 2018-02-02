#ifndef BFIELD_H
#define BFIELD_H

#include <vector>
#include <iostream>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"

typedef double(*callbackFcn)(const double, const double);

//on GPU global variables
extern __device__   callbackFcn BFieldFcnPtr_GPU;
extern __device__   callbackFcn gradBFcnPtr_GPU;

class BField
{
protected:
	//constant array
	std::vector<double> fieldConstArray_m;
	
	//GPU pointer to field constants array
	//double* fieldConstants_d{ nullptr };

	std::string modelName_m{ "" };
	//callback function GPU pointer
	/*
	  I would imagine the process is like so:
	  - Create pointer on GPU with cudaMalloc
	  - Run CUDA kernel to assign pointer to the function to the memory location
	  - Pass allocated pointer to CUDA kernel responsible for running the callback
	  - Run as normal within the CUDA kernel

	  Above is an option, but I did it differently:
	  - Declare global variables on GPU
	  - Call CUDA kernel to set globals to the function desired (from derived classes)

	  Can also use function templates with CUDA!!
	*/

	virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d

public:
	BField() {}
	~BField() {}

	virtual double      getBFieldAtS(double s, double t)=0;
	virtual double      getGradBAtS(double s, double t)=0;
	virtual std::string getName() { return modelName_m; }
};

#endif