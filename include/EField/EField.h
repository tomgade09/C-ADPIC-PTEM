#ifndef EFIELD_H
#define EFIELD_H

#include <memory>
#include <vector>
#include <string>
#include <sstream>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"
#include "ErrorHandling\cudaDeviceMacros.h"

class EElem //inherit from this class
{
protected:
	EElem** this_d{ nullptr }; //not really used on device
	
	#ifndef __CUDA_ARCH__ //host code
	std::string modelName_m;
	#else //device code
	const char* modelName_m;
	#endif /* !__CUDA_ARCH__ */

	__host__ virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__ virtual void deleteEnvironment() = 0;

public:
	__host__ __device__ EElem() {}
	__host__ __device__ virtual ~EElem() {}

	__host__ __device__ virtual double getEFieldAtS(const double s, const double t)=0;

	__host__ virtual std::string getName() { return modelName_m; }
	__host__ EElem** getPtrGPU() { return this_d; }
};

class EField //not meant to be inherited from
{
private:
	EField** this_d{ nullptr };
	
	#ifndef __CUDA_ARCH__ //host code
	std::vector<std::unique_ptr<EElem>> Eelems_m;
	std::vector<std::string> modelNames_m; //could be useful, definitely not essential
	#endif /* !__CUDA_ARCH__ */

	EElem*** Eelems_d; //host: holds ptr to on GPU array, used to increase size, device: holds ptr to on GPU array, used to access elements
	int capacity_d{ 5 }; //denotes size and capacity of E element array on device
	int size_d{ 0 };

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
	
	__device__ int      elemCapacity() { return capacity_d; }
	__device__ int      elemSize() { return size_d; }
	__device__ EElem*** getElemArray() { return Eelems_d; }
	__device__ void     setCap(int capacity) { capacity_d = capacity; }
	__device__ void     setElemArray(EElem*** eelems) { Eelems_d = eelems; }
	
	#ifndef __CUDA_ARCH__ //host code
	#else  //device code
	//__device__ void setElemArray(EElem*** eelems) { Eelems_d = eelems; }
	#endif /* !__CUDA_ARCH__ */
	
	__host__ __device__ double getEFieldAtS(const double s, const double t);
	__host__ EField** getPtrGPU() { return this_d; }
	
	__host__ std::string getEElemsStr()
	#ifndef __CUDA_ARCH__ //host code
	{ std::stringstream out; for (int elem = 0; elem < Eelems_m.size(); elem++) { out << Eelems_m.at(elem)->getName() << ", "; }; return out.str(); }
	#else  //device code
	;
	#endif /* !__CUDA_ARCH__ */
};

#endif /* EFIELD_H */