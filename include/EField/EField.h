#ifndef EFIELD_H
#define EFIELD_H

#include <memory>
#include <vector>
#include <string>

//CUDA includes
#include "host_defines.h"
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include "cuda_profiler_api.h"

//#include "ErrorHandling\cudaErrorCheck.h"
//#include "ErrorHandling\cudaDeviceMacros.h"

class EElem //inherit from this class
{
protected:
	EElem** this_d{ nullptr }; //not really used on device
	
	#ifndef __CUDA_ARCH__ //host code
	std::string modelName_m;
	#else //device code
	const char* modelName_m; //placeholder so code will compile...not actually used
	#endif /* !__CUDA_ARCH__ */

	__host__ virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__ virtual void deleteEnvironment() = 0;

	__host__ __device__ EElem() {}

public:
	__host__ __device__ virtual ~EElem() {}

	__host__ __device__ virtual double getEFieldAtS(const double s, const double t) = 0;

	__host__ virtual std::string name() { return modelName_m; }
	__host__ virtual EElem** getPtrGPU() { return this_d; }
};

class EField
{
private: //not meant to be inherited from
	EField** this_d{ nullptr }; //pointer to EField instance on GPU
	
	#ifndef __CUDA_ARCH__ //host code - need this since there is no CUDA device version of vectors
	std::vector<std::unique_ptr<EElem>> Eelems_m; //Eelems pointers on host
	std::vector<std::string> modelNames_m;        //holds names of elements - or could iterate over elements...
	#endif /* !__CUDA_ARCH__ */

	//GPU container of EElems variables
	EElem*** Eelems_d; //host: holds ptr to on GPU array, used to increase size, device: holds ptr to on GPU array, used to access elements
	int capacity_d{ 5 }; //denotes size and capacity of E element array on device
	int size_d{ 0 };
	//End container variables

	__host__ void setupEnvironment();
	__host__ void deleteEnvironment();

public:
	__host__ __device__ EField()
	{
		#ifndef __CUDA_ARCH__ //host code
		setupEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}
	
	__host__ __device__ ~EField()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#endif /* !__CUDA_ARCH__ */
	}

	__host__   void add(std::unique_ptr<EElem> elem);
	__device__ void add(EElem** elem);

	__device__ int      capacity()                 { return capacity_d; }
	__device__ void     capacity(int cap)          { capacity_d = cap; }
	__device__ int      size()                     { return size_d; }
	__device__ EElem*** elemArray()                { return Eelems_d; }
	__device__ void     elemArray(EElem*** eelems) { Eelems_d = eelems; }
	
	__host__ __device__ double   getEFieldAtS(const double s, const double t);
	__host__            EField** getPtrGPU() { return this_d; }
	
	__host__ EElem*      element(int ind);
	#ifndef __CUDA_ARCH__ //host code
	__host__ std::string getEElemsStr();
	#endif /* !__CUDA_ARCH__ */
};

#endif /* EFIELD_H */