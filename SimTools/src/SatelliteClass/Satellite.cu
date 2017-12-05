//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h" //didn't add to this vs project - each project this class is attached to will have its own variables header
#include "include\Satellite.h"

__global__ void satelliteDetector(double* v_d, double* mu_d, double* z_d, double* detected_v_d, double* detected_mu_d, double* detected_z_d, double altitude, bool upward)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	double z_minus_vdt{ z_d[thdInd] - v_d[thdInd] * DT };//dz };
	
	bool detected{
		detected_z_d[thdInd] < 1 && ( //no detected particle is in the data array at the thread's index already AND
		//detector is facing down and particle crosses altitude in dt
		(!upward && z_d[thdInd] > altitude && z_minus_vdt < altitude)
		|| //OR
		//detector is facing up and particle crosses altitude in dt
		(upward && z_d[thdInd] < altitude && z_minus_vdt > altitude) ) };
	
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
		cudaMalloc((void **)&captureDataGPU_m[iii], sizeof(double) * numberOfParticles_m); //makes room for data of detected particles
		cudaMemset(captureDataGPU_m[iii], 0, sizeof(double) * numberOfParticles_m); //sets values to 0
	}
}

void Satellite::iterateDetector(int numberOfBlocks, int blockSize) {
	satelliteDetector <<< numberOfBlocks, blockSize >>>
		(origDataGPU_m[0], origDataGPU_m[1], origDataGPU_m[2], captureDataGPU_m[0],	captureDataGPU_m[1], captureDataGPU_m[2], altitude_m, upwardFacing_m); }

void Satellite::copyDataToHost()
{
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
	{
		cudaMemcpy(data_m[iii], captureDataGPU_m[iii], sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost);
		cudaMemset(captureDataGPU_m[iii], 0, sizeof(double) * numberOfParticles_m); //sets values to 0
		dataReady_m = true;
	}
	//std::cout << name_m << " Copied.\n" << data_m[0][0] << " " << data_m[0][100] << " " << data_m[1][0] << " " << data_m[1][100] << "\n\n";
}

void Satellite::freeGPUMemory()
{
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
		cudaFree(captureDataGPU_m[iii]);
}

void Satellite::vectorTest(std::vector<double*>& in)
{
	int wrong{ 0 };
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
		for (int jjj = 0; jjj < numberOfParticles_m; jjj++)
			if (in[iii][jjj] != data_m[iii][jjj]) { wrong++; }

	std::cout << "Wrong: " << wrong << "\n";
}