#include "EField/EField.h"

#include <sstream>

//CUDA includes
#include "device_launch_parameters.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/cudaDeviceMacros.h"


__host__ __device__ EField::EField()
{
	#ifndef __CUDA_ARCH__ //host code
	setupEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ __device__ EField::~EField()
{
	#ifndef __CUDA_ARCH__ //host code
	deleteEnvironment();
	#endif /* !__CUDA_ARCH__ */
}


//device global kernels
__global__ void setupEnvironmentGPU_EField(EField** efield, EElem*** eelems)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_EField",
		(*efield) = new EField();
		(*efield)->elemArray(eelems);
	);
}

__global__ void deleteEnvironmentGPU_EField(EField** efield)
{
	ZEROTH_THREAD_ONLY("deleteEnvironmentGPU_EField", delete (*efield));
}

__global__ void addGPU_EField(EField** efield, EElem** elem)
{
	ZEROTH_THREAD_ONLY("addGPU_EField", (*efield)->add(elem));
}

__global__ void increaseCapacity_EField(EField** efield, EElem*** newArray, int capacity)
{
	ZEROTH_THREAD_ONLY("increaseCapacity_EField",
		EElem*** oldArray{ (*efield)->elemArray() };

		for (int elem = 0; elem < (*efield)->size(); elem++)
			newArray[elem] = oldArray[elem];

		(*efield)->capacity(capacity);
		(*efield)->elemArray(newArray); //still retaining the pointer to this memory on host, so no big deal if it's lost here
	);
}


//EField functions
__host__ std::string EField::getEElemsStr() const
{
	std::stringstream out;
	for (int elem = 0; elem < size_d; elem++) { out << element(elem)->name() << ", "; }
	return out.str();
}

void EField::setupEnvironment()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(EField**)));              //allocate memory for EField**
	CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //allocate memory for EElem** array
	CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));       //clear memory

	setupEnvironmentGPU_EField <<< 1, 1 >>> (this_d, Eelems_d);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void EField::deleteEnvironment()
{
	deleteEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(Eelems_d));
}

#ifndef __CUDA_ARCH__ //host code
__host__ EElem* EField::element(int ind) const
{
	return Eelems_m.at(ind).get();
}

__host__ void EField::add(std::unique_ptr<EElem> eelem)
{
	if (capacity_d == size_d)
	{
		EElem*** oldArray{ Eelems_d }; //retain so we can cudaFree at the end
		capacity_d += 5;
		
		CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //create new array that is 5 larger in capacity than the previous
		CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));

		increaseCapacity_EField <<< 1, 1 >>> (this_d, Eelems_d, capacity_d);
		CUDA_KERNEL_ERRCHK();

		CUDA_API_ERRCHK(cudaFree(oldArray));
	}

	//add elem to dev
	addGPU_EField <<< 1, 1 >>> (this_d, eelem->getPtrGPU());
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

__host__ __device__ double EField::getEFieldAtS(const double s, const double t) const
{
	double ret{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	for (auto& elem : Eelems_m) //std::vector of std::unique_ptr<EElem>'s
		ret += elem->getEFieldAtS(s, t);
	#else //device code
	for (int elem = 0; elem < size_d; elem++) //c-style array of EElem*'s
		ret += (*(Eelems_d[elem]))->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}