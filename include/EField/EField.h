#ifndef EFIELD_H
#define EFIELD_H

#include <memory>
#include <vector>
#include <string>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"

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
	std::vector<EElem**> Eelems_d;
	std::vector<std::unique_ptr<EElem>> Eelems_m;
	std::vector<std::string> modelNames_m;
	#else //device code
	EElem*** Eelems_m; //array of E Field elements
	#endif /* !__CUDA_ARCH__ */

	int capacity_m{ 0 }; //need these on device and it's seriously a pain to enclose everything in #ifndef's
	int size_m{ 0 };  //so these will be on both host and device
	
	void setupEnvironment();
	void deleteEnvironment();

public:
	__host__ __device__ EField(int reserveNumElems)
	{
		#ifndef __CUDA_ARCH__ //host code
		capacity_m = reserveNumElems;
		Eelems_m.reserve(capacity_m);
		setupEnvironment();
		#else  //device code
		Eelems_m = new EElem**[reserveNumElems];
		capacity_m = reserveNumElems;
		#endif /* !__CUDA_ARCH__ */
	}
	
	__host__ __device__ ~EField()
	{
		#ifndef __CUDA_ARCH__ //host code
		deleteEnvironment();
		#else  //device code
		delete[] Eelems_m;
		#endif /* !__CUDA_ARCH__ */
	}

	#ifndef __CUDA_ARCH__ //host code
	__host__            void   add(std::unique_ptr<EElem> elem);
	#else //device code
	__device__          void   add(EElem** elem);
	#endif /* !__CUDA_ARCH__ */
	
	__host__ __device__ double getEFieldAtS(const double s, const double t);
	__host__ EField** getPtrGPU() { return this_d; }
};

#endif /* EFIELD_H */