#include "EField\QSPS.h"

#include "device_launch_parameters.h"
#include "ErrorHandling\cudaErrorCheck.h"
#include "ErrorHandling\cudaDeviceMacros.h"

__global__ void setupEnvironmentGPU_QSPS(EElem** qsps, double* altMin, double* altMax, double* magnitude, int numRegions)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_QSPS", (*qsps) = new QSPS(altMin, altMax, magnitude, numRegions)); //this overloaded constructor is only compiled in the case where __CUDA_ARCH__ is defined
}

__global__ void deleteEnvironmentGPU_QSPS(EElem** qsps)
{
	ZEROTH_THREAD_ONLY("deleteEnvironmentGPU_QSPS", delete (*qsps));
}

#ifndef __CUDA_ARCH__ //host code
__host__ const std::vector<double>& QSPS::altMin() const
{
	return altMin_m;
}

__host__ const std::vector<double>& QSPS::altMax() const 
{
	return altMax_m;
}

__host__ const std::vector<double>& QSPS::magnitude() const 
{
	return magnitude_m;
}
#endif

__host__ QSPS::QSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude) :
	EElem(), numRegions_m{ (int)magnitude.size() }
{
	if (magnitude.size() != altMin.size() || magnitude.size() != altMax.size())
		throw std::invalid_argument("QSPS::QSPS: invalid parameters passed in magnitude, altMin, altMax: resolved vector lengths are not equal");

	#ifndef __CUDA_ARCH__ //host code
	altMin_m = altMin;       //unfortunately this wrapping is necessary
	altMax_m = altMax;       //as the vectors above also have to be wrapped
	magnitude_m = magnitude; //in an ifndef/endif block so this will compile
	modelName_m = "QSPS";
	#endif /* !__CUDA_ARCH__ */

	setupEnvironment();
}

__host__ __device__ QSPS::~QSPS()
{
	#ifndef __CUDA_ARCH__ //host code
	deleteEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ void QSPS::setupEnvironment()
{
	#ifndef __CUDA_ARCH__ //host code
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(QSPS**))); //malloc for ptr to ptr to GPU QSPS Obj
	CUDA_API_ERRCHK(cudaMalloc((void **)&altMin_d, altMin_m.size() * sizeof(double))); //array of altitude min bounds
	CUDA_API_ERRCHK(cudaMalloc((void **)&altMax_d, altMax_m.size() * sizeof(double)));
	CUDA_API_ERRCHK(cudaMalloc((void **)&magnitude_d, magnitude_m.size() * sizeof(double))); //array of E magnitude between above min/max
	CUDA_API_ERRCHK(cudaMemcpy(altMin_d, altMin_m.data(), altMin_m.size() * sizeof(double), cudaMemcpyHostToDevice));
	CUDA_API_ERRCHK(cudaMemcpy(altMax_d, altMax_m.data(), altMax_m.size() * sizeof(double), cudaMemcpyHostToDevice));
	CUDA_API_ERRCHK(cudaMemcpy(magnitude_d, magnitude_m.data(), magnitude_m.size() * sizeof(double), cudaMemcpyHostToDevice));

	setupEnvironmentGPU_QSPS <<< 1, 1 >>> (this_d, altMin_d, altMax_d, magnitude_d, (int)(magnitude_m.size()));
	CUDA_KERNEL_ERRCHK_WSYNC(); //creates GPU instance of QSPS
	#endif /* !__CUDA_ARCH__ */
}

__host__ void QSPS::deleteEnvironment()
{
	deleteEnvironmentGPU_QSPS <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(altMin_d)); //On device
	CUDA_API_ERRCHK(cudaFree(altMax_d));
	CUDA_API_ERRCHK(cudaFree(magnitude_d));
}

__host__ __device__ double QSPS::getEFieldAtS(const double s, const double t) const
{
	#ifndef __CUDA_ARCH__ //host code
	for (int ind = 0; ind < magnitude_m.size(); ind++)
	{
		if (s >= altMin_m.at(ind) && s <= altMax_m.at(ind))
			return magnitude_m.at(ind);
	}
	#else //device code
	for (int ind = 0; ind < numRegions_m; ind++)
	{
		if (s >= altMin_d[ind] && s <= altMax_d[ind])
			return magnitude_d[ind];
	}
	#endif /* !__CUDA_ARCH__ */

	return 0.0;
}