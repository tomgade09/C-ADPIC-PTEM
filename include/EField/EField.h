#ifndef EFIELD_H
#define EFIELD_H

#include <memory>
#include <vector>
#include <string>

//CUDA includes
#include "host_defines.h"

//Project includes
#include "dlldefines.h"

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
	__host__ __device__ EElem(const EElem&) = delete;
	__host__ __device__ EElem& operator=(const EElem&) = delete;

	__host__ __device__ virtual double getEFieldAtS(const double s, const double t) const = 0;

	__host__ virtual std::string name()  const { return modelName_m; }
	__host__ virtual EElem** getPtrGPU() const { return this_d; }
};

class EField final //not meant to be inherited from
{
private:
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
	__host__ __device__ EField();
	__host__ __device__ ~EField();
	__host__ __device__ EField(const EField&) = delete;
	__host__ __device__ EField& operator=(const EField&) = delete;

	__host__   void add(std::unique_ptr<EElem> elem);
	__device__ void add(EElem** elem);

	__host__ __device__ int      capacity()  const { return capacity_d; }
	__host__ __device__ int      size()      const { return size_d; }
	__device__ void     capacity(int cap)          { capacity_d = cap; }
	__device__ EElem*** elemArray()          const { return Eelems_d; }
	__device__ void     elemArray(EElem*** eelems) { Eelems_d = eelems; }
	
	__host__ __device__ double   getEFieldAtS(const double s, const double t) const;
	__host__            EField** getPtrGPU() const { return this_d; }
	
	__host__ EElem*      element(int ind) const;
	__host__ std::string getEElemsStr() const;
};

#endif /* EFIELD_H */