#include "include\AlfvenLUT.h"

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)
__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries);

extern const int SIMCHARSIZE;

__host__ __device__ double alfvenWaveEbyLUT(double** LUT, double z, double simtime, double omegaE)
{//E Field in the direction of B (radially outward)
	if (LUT == nullptr)
		return 0.0;
	if (z > RADIUS_EARTH) //in case z is passed in as m, not Re, convert to Re
		z = z / RADIUS_EARTH;
	if (z < LUT[0][0] || z > LUT[0][2950])
		return 0.0;

	double offset{ LUT[0][0] };
	int stepsFromZeroInd{ static_cast<int>(floor((z - offset) / (LUT[0][1] - LUT[0][0]))) }; //only works for constant bin size - if the binsize changes throughout LUT, need to iterate which will take longer

	//y = mx + b
	double linearInterpReal{ ((LUT[1][stepsFromZeroInd + 1] - LUT[1][stepsFromZeroInd]) / (LUT[0][stepsFromZeroInd + 1] - LUT[0][stepsFromZeroInd])) *
		(z - LUT[0][stepsFromZeroInd]) + LUT[1][stepsFromZeroInd] };
	double linearInterpImag{ ((LUT[2][stepsFromZeroInd + 1] - LUT[2][stepsFromZeroInd]) / (LUT[0][stepsFromZeroInd + 1] - LUT[0][stepsFromZeroInd])) *
		(z - LUT[0][stepsFromZeroInd]) + LUT[2][stepsFromZeroInd] };

	//E-par = (column 2)*cos(omega*t) + (column 3)*sin(omega*t), omega - angular frequency of wave
	return (linearInterpReal * cos(omegaE * simtime) + linearInterpImag * sin(omegaE * simtime)) / 1000; //LUT E is in mV / m
}

void AlfvenLUT::initializeFollowOn()
{
	useAlfLUT_m = true;
	CUDA_CALL(cudaMalloc((void **)&gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), numOfColsLUT_m * numOfEntrLUT_m * sizeof(double)));
	CUDA_CALL(cudaMalloc((void **)&gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1), numOfColsLUT_m * sizeof(double*)));
}

void AlfvenLUT::copyDataToGPUFollowOn()
{
	CUDA_CALL(cudaMemcpy(gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), elcFieldLUT_m[0], numOfColsLUT_m * numOfEntrLUT_m * sizeof(double), cudaMemcpyHostToDevice));
	setup2DArray <<< 1, 1 >>> (gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1)), numOfColsLUT_m, numOfEntrLUT_m);
	CUDA_CALL(cudaMemcpy(gpuDblMemoryPointers_m.at(2 * particleTypes_m.size()) + (SIMCHARSIZE/sizeof(double) - 1), &omegaE_m, sizeof(double), cudaMemcpyHostToDevice));
}

void AlfvenLUT::iterateSimulationFollowOnPreLoop()
{
	if (useAlfLUT_m && gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1) == nullptr)
		std::cout << "Warning: LUT pointer is a nullptr.  Alfven wave function will return 0.0.  Continuing.\n";
}

void AlfvenLUT::iterateSimulationFollowOnInsideLoop()
{
	return;
}

void AlfvenLUT::iterateSimulationFollowOnPostLoop()
{
	return;
}

void AlfvenLUT::copyDataToHostFollowOn()
{
	return;
}

void AlfvenLUT::freeGPUMemoryFollowOn()
{
	return;
}