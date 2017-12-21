//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h" //didn't add to this vs project - each project this class is attached to will have its own variables header
#include "SatelliteClass\Satellite.h"

__global__ void setupKernel(double* array1D, double** array2D, int cols, int entrs)
{
	if (blockIdx.x * blockDim.x + threadIdx.x != 0)
		return;

	for (int iii = 0; iii < cols; iii++)
		array2D[iii] = &array1D[iii * entrs];
}

__global__ void satelliteDetector(double** data_d, double** capture_d, double simtime, double altitude, bool upward)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;

	double* v_d; double* mu_d; double* z_d; double* simtime_d;
	double* detected_v_d; double* detected_mu_d; double* detected_z_d;
	v_d = data_d[0]; mu_d = data_d[1]; z_d = data_d[2]; simtime_d = capture_d[3];
	detected_v_d = capture_d[0]; detected_mu_d = capture_d[1]; detected_z_d = capture_d[2];

	double z_minus_vdt{ z_d[thdInd] - v_d[thdInd] * DT };
	
	if (simtime == 0) //not sure I fully like this, but it works
		simtime_d[thdInd] = -1.0;

	if (
		(detected_z_d[thdInd] < 1) && ( //no detected particle is in the data array at the thread's index already AND
		//detector is facing down and particle crosses altitude in dt
		((!upward) && (z_d[thdInd] > altitude) && (z_minus_vdt < altitude))
		|| //OR
		//detector is facing up and particle crosses altitude in dt
		((upward) && (z_d[thdInd] < altitude) && (z_minus_vdt > altitude)) ) 
		)
	{
		detected_v_d[thdInd] = v_d[thdInd];
		detected_mu_d[thdInd] = mu_d[thdInd];
		detected_z_d[thdInd] = z_d[thdInd];
		simtime_d[thdInd] = simtime;
	}//particle not removed from sim
}

void Satellite::initializeSatelliteOnGPU()
{
	cudaMalloc((void **)&satCaptureGPU_m, sizeof(double) * (numberOfAttributes_m + 1) * numberOfParticles_m); //makes room for data of detected particles
	cudaMemset((void **)&satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 1) * numberOfParticles_m); //sets values to 0
	cudaMalloc((void **)&dblppGPU_m[1], sizeof(double*) * numberOfAttributes_m);

	setupKernel <<< 1, 1 >>> (satCaptureGPU_m, dblppGPU_m[1], numberOfAttributes_m + 1, numberOfParticles_m);
}

void Satellite::iterateDetector(int numberOfBlocks, int blockSize, double simtime) {
	satelliteDetector <<< numberOfBlocks, blockSize >>>	(dblppGPU_m[0], dblppGPU_m[1], simtime, altitude_m, upwardFacing_m); }

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, z, time][particle number]
	cudaMemcpy(data_m[0], satCaptureGPU_m, sizeof(double) * (numberOfAttributes_m + 1) * numberOfParticles_m, cudaMemcpyDeviceToHost);
	cudaMemset(satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 1) * numberOfParticles_m); //sets values to 0

	dataReady_m = true;
}

void Satellite::freeGPUMemory()
{
	cudaFree(satCaptureGPU_m);
	cudaFree(dblppGPU_m[1]); //DO NOT FREE dblppGPU_m[0] - this is the 2D data array that the sim uses (not the satellite)
}

void Satellite::vectorTest(std::vector<double*>& in)
{
	int wrong{ 0 };
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
		for (int jjj = 0; jjj < numberOfParticles_m; jjj++)
			if (in[iii][jjj] != data_m[iii][jjj]) { wrong++; }

	std::cout << "Wrong: " << wrong << "\n";
}