#include "EField\EField.h"

//device global setup kernels
__global__ void setupEnvironmentGPU_EField(EField** efield, int reserveNumElems)
{
	(*efield) = new EField(reserveNumElems);
}

__global__ void deleteEnvironmentGPU_EField(EField** efield)
{
	delete (*efield);
}

//EField functions
void EField::setupEnvironment()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(EField**)));
	setupEnvironmentGPU_EField <<< 1, 1 >>> (this_d, capacity_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void EField::deleteEnvironment()
{
	deleteEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
}

#ifndef __CUDA_ARCH__ //host code
__global__ void addGPU_EField(EField** efield, EElem** eelem);

__host__ void EField::add(std::unique_ptr<EElem> eelem)
{
	addGPU_EField <<< 1, 1 >>> (this_d, eelem->getPtrGPU());
	CUDA_KERNEL_ERRCHK_WSYNC();
	Eelems_m.push_back(std::move(eelem));
}
#else //device code
__global__ void addGPU_EField(EField** efield, EElem** elem)
{
	(*efield)->add(elem); //has to be enclosed in ifndef/else macro because on host add() takes a unique_ptr
}

__device__ void EField::add(EElem** elem)
{
	if (capacity_m == size_m)
	{
		printf("CUDA: EField::add: array of EField pointers is full, resizing array: size: %f, capacity: %f\n", (float)size_m, (float)capacity_m);
		EElem*** tmp = new EElem**[capacity_m + 5];
		
		for (int elem = 0; elem < size_m; elem++)
			tmp[elem] = Eelems_m[elem];

		delete[] Eelems_m;
		Eelems_m = tmp;
		capacity_m += 5;
	}

	Eelems_m[size_m] = elem;
	size_m++;
}
#endif /* !__CUDA_ARCH__ */

__host__ __device__ double EField::getEFieldAtS(const double s, const double t)
{
	double ret{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	for (int elem = 0; elem < Eelems_m.size(); elem++) //std::vector of std::unique_ptr<EElem>'s
		ret += Eelems_m.at(elem)->getEFieldAtS(s, t);
	#else //device code
	for (int elem = 0; elem < size_m; elem++) //c-style array of EElem**'s
		ret += (*Eelems_m[elem])->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}