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
	setupEnvironmentGPU_EField <<< 1, 1 >>> (this_d, numElems_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void EField::deleteEnvironment()
{
	deleteEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

#ifndef __CUDA_ARCH__
__global__ void addGPU_EField(EField** efield, EElem** eelem);

__host__ void EField::add(std::unique_ptr<EElem> elem)
{
	Eelems_m.push_back(std::move(elem));
	addGPU_EField <<< 1, 1 >>> (this_d, elem->getPtrGPU());
	CUDA_KERNEL_ERRCHK_WSYNC();
}
#else
__global__ void addGPU_EField(EField** efield, EElem** elem)
{
	(*efield)->add(elem);
}

__device__ void EField::add(EElem** elem)
{
	if (elemsSize_m == numElems_m)
	{
		printf("CUDA: EFieldContainer::add: array of EField pointers is full, resizing array\n");
		EElem*** tmp = new EElem**[elemsSize_m + 5];
		
		for (int elem = 0; elem < numElems_m; elem++)
			tmp[elem] = Eelems_m[elem];

		delete[] Eelems_m;
		Eelems_m = tmp;
		elemsSize_m += 5;
	}

	Eelems_m[numElems_m] = elem;
	numElems_m++;
}
#endif /* !__CUDA_ARCH__ */

__host__ __device__ double EField::getEFieldAtS(const double s, const double t)
{
	double ret{ 0.0 };

	#ifndef __CUDA_ARCH__
	for (int elem = 0; elem < Eelems_m.size(); elem++)
		ret += Eelems_m.at(elem)->getEFieldAtS(s, t);
	#else
	for (int elem = 0; elem < numElems_m; elem++)
		ret += (*Eelems_m[elem])->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}