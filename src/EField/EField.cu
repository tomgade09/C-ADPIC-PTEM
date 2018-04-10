#include "EField\EField.h"

//device global kernels
__global__ void setupEnvironmentGPU_EField(EField** efield)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_EField", (*efield) = new EField());
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
	EElem*** oldArray{ (*efield)->getElemArray() };
	int size{ (*efield)->elemSize() };

	for (int elem = 0; elem < size; elem++)
		newArray[elem] = oldArray[elem];

	(*efield)->setCap(capacity);
	(*efield)->setElemArray(newArray); //still retaining the pointer to this memory on host, so no big deal if it's lost here
	);
}


//EField functions
void EField::setupEnvironment()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(EField**)));
	setupEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * 5));
	CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * 5));
}

void EField::deleteEnvironment()
{
	deleteEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(Eelems_d));
}

#ifndef __CUDA_ARCH__ //host code
__host__ void EField::add(std::unique_ptr<EElem> eelem)
{
	if (capacity_d == size_d)
	{
		EElem*** oldArray{ Eelems_d }; //retain so we can cudaFree at the end
		CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * (capacity_d + 5))); //create new array that is 5 larger in capacity than the previous
		CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * (capacity_d + 5)));
		capacity_d += 5;

		increaseCapacity_EField <<< 1, 1 >>> (this_d, Eelems_d, capacity_d);
		CUDA_KERNEL_ERRCHK();

		CUDA_API_ERRCHK(cudaFree(oldArray));
	}

	//add elem to host
	addGPU_EField <<< 1, 1 >>> (this_d, eelem->getPtrGPU());
	CUDA_KERNEL_ERRCHK_WSYNC();
	Eelems_m.push_back(std::move(eelem));
	size_d++;
}
#endif /* !__CUDA_ARCH__ */

__device__ void EField::add(EElem** newElem) //needs this because capacity and size are in blocks
{
	Eelems_d[size_d] = newElem;
	size_d++;
}

__host__ __device__ double EField::getEFieldAtS(const double s, const double t)
{
	double ret{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	for (int elem = 0; elem < Eelems_m.size(); elem++) //std::vector of std::unique_ptr<EElem>'s
		ret += Eelems_m.at(elem)->getEFieldAtS(s, t);
	#else //device code
	for (int elem = 0; elem < size_d; elem++) //c-style array of EElem**'s
		ret += (*Eelems_d[elem])->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}