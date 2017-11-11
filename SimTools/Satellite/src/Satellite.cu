//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h" //didn't add to this vs project - each project this class is attached to will have its own variables header
#include "include\Satellite.h"



#include <stdio.h>
#include <stdlib.h>
#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
    printf("Error at %s:%d  Error Number: %d\n",__FILE__,__LINE__, EXIT_FAILURE);}} while(0)



__global__ void satelliteDetector(double* v_d, double* mu_d, double* z_d, double* detected_v_d, double* detected_mu_d, double* detected_z_d, double altitude, bool upward)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	double dz{ v_d[thdInd] * DT };
	double z_minus_vdt{ z_d[thdInd] - dz };
	
	bool detected{
		detected_z_d[thdInd] < 1 && ( //no detected particle is in the data array at the thread's index already AND
		//particle is traveling upward, detector is facing down, and particle crosses altitude
		(dz > 0 && !upward && z_d[thdInd] > altitude && z_minus_vdt < altitude)
		|| //OR
		//particle is traveling downward, detector is facing up, and particle crosses altitude
		(dz < 0 && upward && z_d[thdInd] < altitude && z_minus_vdt > altitude) ) };
	
	if (detected)
	{
		detected_v_d[thdInd] = v_d[thdInd];
		detected_mu_d[thdInd] = mu_d[thdInd];
		detected_z_d[thdInd] = z_d[thdInd];
	}//particle not removed from sim
}

void Satellite::initializeSatelliteOnGPU()
{
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
	{
		CUDA_CALL(cudaMalloc((void **)&GPUdata_m[iii], sizeof(double) * numberOfParticles_m)); //makes room for data of detected particles
		CUDA_CALL(cudaMemset(GPUdata_m[iii], 0, sizeof(double) * numberOfParticles_m)); //sets values to 0
	}
}

void Satellite::iterateDetector(int numberOfBlocks, int blockSize, double** simData)
{
	satelliteDetector <<< numberOfBlocks, blockSize >>> (simData[0], simData[1], simData[2], GPUdata_m[0], GPUdata_m[1], GPUdata_m[2], altitude_m, upwardFacing_m);
}

void Satellite::copyDataToHost()
{
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
	{
		CUDA_CALL(cudaMemcpy(data_m[iii], GPUdata_m[iii], sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemset(GPUdata_m[iii], 0, sizeof(double) * numberOfParticles_m)); //sets values to 0
		dataReady_m = true;
	}
}

void Satellite::freeGPUMemory()
{
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
		cudaFree(GPUdata_m[iii]);
}