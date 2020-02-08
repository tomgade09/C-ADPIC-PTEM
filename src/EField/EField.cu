#include "EField/EField.h"
#include "EField/QSPS.h"
//#include "EField/AlfvenLUT.h"

#include <sstream>

//CUDA includes
#include "device_launch_parameters.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/cudaDeviceMacros.h"

using std::string;
using std::to_string;
using std::stringstream;

__host__ __device__ EElem::EElem(const char* modelName) : name_m{ modelName }
{

}

__host__ __device__ EElem::~EElem()
{

}

__host__ string EElem::name() const
{
	return name_m;
}

__host__ EElem** EElem::getPtrGPU() const
{
	return this_d;
}


//EField ctor, dtor
__host__ __device__ EField::EField(bool useGPU) : useGPU_m{ useGPU }
{
	#ifndef __CUDA_ARCH__ //host code
	if(useGPU_m) setupEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ EField::EField(string serialFolder)
{
	deserialize(serialFolder);
	if(useGPU_m) setupEnvironment();
}

__host__ __device__ EField::~EField()
{
	#ifndef __CUDA_ARCH__ //host code
	if(useGPU_m)deleteEnvironment();
	#endif /* !__CUDA_ARCH__ */
}


//device global kernels
namespace EField_d
{
	__global__ void setupEnvironmentGPU(EField** efield, EElem*** eelems)
	{
		ZEROTH_THREAD_ONLY(
			(*efield) = new EField();
			(*efield)->elemArray(eelems);
		);
	}

	__global__ void deleteEnvironmentGPU(EField** efield)
	{
		ZEROTH_THREAD_ONLY(delete (*efield));
	}

	__global__ void addGPU(EField** efield, EElem** elem)
	{
		ZEROTH_THREAD_ONLY((*efield)->add(elem));
	}

	__global__ void increaseCapacity(EField** efield, EElem*** newArray, int capacity)
	{
		ZEROTH_THREAD_ONLY(
			EElem*** oldArray{ (*efield)->elemArray() };

			for (int elem = 0; elem < (*efield)->size(); elem++)
				newArray[elem] = oldArray[elem];

			(*efield)->capacity(capacity);
			(*efield)->elemArray(newArray); //still retaining the pointer to this memory on host, so no big deal if it's lost here
		);
	}
}

//EField functions
__host__ string EField::getEElemsStr() const
{
	stringstream out;
	for (int elem = 0; elem < size_d; elem++) { out << element(elem)->name() << ", "; }
	return out.str();
}

void EField::setupEnvironment()
{
	CUDA_API_ERRCHK(cudaMalloc((void**)&this_d, sizeof(EField**)));               //allocate memory for EField**
	CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //allocate memory for EElem** array
	CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));       //clear memory
	
	EField_d::setupEnvironmentGPU <<< 1, 1 >>> (this_d, Eelems_d);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void EField::deleteEnvironment()
{
	EField_d::deleteEnvironmentGPU <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(Eelems_d));
}

#ifndef __CUDA_ARCH__ //host code
__host__ EElem* EField::element(int ind) const
{
	return Eelems_m.at(ind).get();
}

__host__ void EField::add(unique_ptr<EElem> eelem)
{
	if (capacity_d == size_d)
	{
		EElem*** oldArray{ Eelems_d }; //retain so we can cudaFree at the end
		capacity_d += 5;
		
		CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //create new array that is 5 larger in capacity than the previous
		CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));

		EField_d::increaseCapacity <<< 1, 1 >>> (this_d, Eelems_d, capacity_d);
		CUDA_KERNEL_ERRCHK();

		CUDA_API_ERRCHK(cudaFree(oldArray));
	}

	//add elem to dev
	EField_d::addGPU <<< 1, 1 >>> (this_d, eelem->getPtrGPU());
	CUDA_KERNEL_ERRCHK_WSYNC();
	
	//add elem to host
	Eelems_m.push_back(std::move(eelem));
	size_d++;
}
#endif /* !__CUDA_ARCH__ */

__device__ void EField::add(EElem** newElem)
{
	Eelems_d[size_d] = newElem;
	size_d++;
}

__host__ __device__ int EField::capacity() const
{
	return capacity_d;
}

__host__ __device__ int EField::size() const
{
	return size_d;
}

__device__ void EField::capacity(int cap)
{
	capacity_d = cap;
}

__device__ EElem*** EField::elemArray() const
{
	return Eelems_d;
}

__device__ void EField::elemArray(EElem*** eelems)
{
	Eelems_d = eelems;
}
	
__host__ EField** EField::getPtrGPU() const
{
	return this_d;
}

__host__ __device__ Vperm EField::getEFieldAtS(const meters s, const seconds t) const
{
	tesla ret{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	for (auto& elem : Eelems_m) //vector of unique_ptr<EElem>'s
		ret += elem->getEFieldAtS(s, t);
	#else //device code
	for (int elem = 0; elem < size_d; elem++) //c-style array of EElem*'s
		ret += (*(Eelems_d[elem]))->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}