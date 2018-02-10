#include "EField\QSPS.h"

__global__ void setupEnvironmentGPU_QSPS(EElem** qsps, double* altMinMax, double* magnitude, int numRegions)
{
	if (threadIdx.x == 0 && blockIdx.x == 0)
		(*qsps) = new QSPS(altMinMax, magnitude, numRegions); //this overloaded constructor is only compiled in the case where __CUDA_ARCH__ is defined
}

__global__ void deleteEnvironmentGPU_QSPS(EElem** qsps)
{
	delete (*qsps);
}

__host__ void QSPS::setupEnvironment()
{
	#ifndef __CUDA_ARCH__ //host code
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(QSPS**))); //malloc for ptr to ptr to GPU QSPS Obj
	CUDA_API_ERRCHK(cudaMalloc((void **)&altMinMax_d, altMinMax_m.size() * sizeof(double))); //array of altitude min bounds and max bounds
	CUDA_API_ERRCHK(cudaMalloc((void **)&magnitude_d, magnitude_m.size() * sizeof(double))); //array of E magnitude between above min/max
	CUDA_API_ERRCHK(cudaMemcpy(altMinMax_d, altMinMax_m.data(), altMinMax_m.size() * sizeof(double), cudaMemcpyHostToDevice));
	CUDA_API_ERRCHK(cudaMemcpy(magnitude_d, magnitude_m.data(), magnitude_m.size() * sizeof(double), cudaMemcpyHostToDevice));

	setupEnvironmentGPU_QSPS <<< 1, 1 >>> (this_d, altMinMax_d, magnitude_d, static_cast<int>(magnitude_m.size()));
	CUDA_KERNEL_ERRCHK_WSYNC(); //creates GPU instance of QSPS
	#endif /* !__CUDA_ARCH__ */
}

__host__ void QSPS::deleteEnvironment()
{
	deleteEnvironmentGPU_QSPS <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(altMinMax_d)); //On device
	CUDA_API_ERRCHK(cudaFree(magnitude_d));
}

__host__ __device__ double QSPS::getEFieldAtS(const double s, const double t)
{
	#ifndef __CUDA_ARCH__ //host code
	for (int ind = 0; ind < magnitude_m.size(); ind++)
	{
		if (s > altMinMax_m.at(2 * ind) && s < altMinMax_m.at(2 * ind + 1))
			return magnitude_m.at(ind);
	}

	return 0.0;
	#else //device code
	for (int ind = 0; ind < numRegions_m; ind++)
	{
		if (s > altMinMax_d[2 * ind] && s < altMinMax_d[2 * ind + 1])
			return magnitude_d[ind];
	}

	return 0.0;
	#endif /* !__CUDA_ARCH__ */
}