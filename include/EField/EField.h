#ifndef EFIELD_H
#define EFIELD_H

#include <memory>
#include <vector>
#include <string>

//CUDA includes
#include "cuda_runtime.h"

//Project includes
#include "dlldefines.h"
#include "utils/unitsTypedefs.h"

using std::string;
using std::vector;
using std::unique_ptr;

class EElem //inherit from this class
{
protected:
	EElem** this_d{ nullptr }; //not really used on device
	
	const char* name_m;

	__host__            virtual void setupEnvironment() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d
	__host__            virtual void deleteEnvironment() = 0;
	__host__            virtual void deserialize(string serialFolder, int nameIndex) = 0;
	
	__host__ __device__ EElem(const char* modelName);

public:
	__host__ __device__ virtual ~EElem();
	__host__ __device__ EElem(const EElem&) = delete;
	__host__ __device__ EElem& operator=(const EElem&) = delete;

	__host__ __device__ virtual Vperm getEFieldAtS(const meters s, const seconds t) const = 0;

	__host__            virtual string name() const;
	__host__            virtual EElem** getPtrGPU() const;

	__host__            virtual void serialize(string serialFolder) const = 0;
};

class EField final //not meant to be inherited from
{
private:
	EField** this_d{ nullptr }; //pointer to EField instance on GPU
	
	#ifndef __CUDA_ARCH__ //host code - need this since there is no CUDA device version of vectors
	vector<unique_ptr<EElem>> Eelems_m; //Eelems pointers on host
	vector<string> modelNames_m;        //holds names of elements - or could iterate over elements...
	#endif /* !__CUDA_ARCH__ */  //the above vector references need these dev exclusion blocks or there is a cudaIllegalMemoryAccess error

	//GPU container of EElems variables
	EElem*** Eelems_d; //host: holds ptr to on GPU array, used to increase size, device: holds ptr to on GPU array, used to access elements
	int capacity_d{ 5 }; //denotes size and capacity of E element array on device
	int size_d{ 0 };
	//End container variables

	bool useGPU_m{ true };

	__host__ void setupEnvironment();
	__host__ void deleteEnvironment();

	__host__ void deserialize(string serialFolder);

public:
	__host__ __device__ EField(bool useGPU = true);
	__host__            EField(string serialFolder);
	__host__ __device__ ~EField();
	__host__ __device__ EField(const EField&) = delete;
	__host__ __device__ EField& operator=(const EField&) = delete;

	__host__            void add(unique_ptr<EElem> elem);
	__device__          void add(EElem** elem);

	__host__ __device__ int capacity() const;
	__host__ __device__ int size() const;
	__device__          void capacity(int cap);
	__device__          EElem*** elemArray() const;
	__device__          void elemArray(EElem*** eelems);
	
	__host__ __device__ Vperm    getEFieldAtS(const meters s, const seconds t) const;
	__host__            EField** getPtrGPU() const;
	
	__host__            EElem* element(int ind) const;
	__host__            string getEElemsStr() const;

	__host__            void serialize(string serialFolder) const;
};

#endif /* EFIELD_H */
