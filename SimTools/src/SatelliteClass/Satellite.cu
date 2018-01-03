//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h"
#include "SatelliteClass\Satellite.h"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)

__global__ void setupKernel(double* array1D, double** array2D, int cols, int entrs)
{
	if (blockIdx.x * blockDim.x + threadIdx.x != 0)
		return;

	for (int iii = 0; iii < cols; iii++)
		array2D[iii] = &array1D[iii * entrs];
}

__global__ void satelliteDetector(double** data_d, double** capture_d, double simtime, double dt, double altitude, bool upward)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	
	double* v_d; double* mu_d; double* z_d; double* simtime_d; double* index_d;
	double* detected_v_d; double* detected_mu_d; double* detected_z_d;
	v_d = data_d[0]; mu_d = data_d[1]; z_d = data_d[2]; simtime_d = capture_d[3]; index_d = capture_d[4];
	detected_v_d = capture_d[0]; detected_mu_d = capture_d[1]; detected_z_d = capture_d[2];

	double z_minus_vdt{ z_d[thdInd] - v_d[thdInd] * dt };
	
	if (simtime == 0) //not sure I fully like this, but it works
		simtime_d[thdInd] = -1.0;

	if (
		(detected_z_d[thdInd] < 1) &&
			( //no detected particle is in the data array at the thread's index already AND
				//detector is facing down and particle crosses altitude in dt
				((!upward) && (z_d[thdInd] > altitude) && (z_minus_vdt < altitude))
					|| //OR
				//detector is facing up and particle crosses altitude in dt
				((upward) && (z_d[thdInd] < altitude) && (z_minus_vdt > altitude))
			)
		)
	{
		detected_v_d[thdInd] = v_d[thdInd];
		detected_mu_d[thdInd] = mu_d[thdInd];
		detected_z_d[thdInd] = z_d[thdInd];
		simtime_d[thdInd] = simtime;
		index_d[thdInd] = thdInd;
	}//particle not removed from sim
}

void Satellite::initializeSatelliteOnGPU()
{
	dataAllocateNewMsmtVector(); //make room for the first measurement data set

	CUDA_CALL(cudaMalloc((void **)&satCaptureGPU_m, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //makes room for data of detected particles
	CUDA_CALL(cudaMemset(satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //sets values to 0
	CUDA_CALL(cudaMalloc((void **)&dblppGPU_m[1], sizeof(double*) * numberOfAttributes_m));

	setupKernel <<< 1, 1 >>> (satCaptureGPU_m, dblppGPU_m[1], numberOfAttributes_m + 2, numberOfParticles_m);
}

void Satellite::iterateDetector(int blockSize, double simtime, double dt)
{
	if (numberOfParticles_m % blockSize != 0)
		std::cout << "Warning: " << name_m << ": Satellite::iterateDetector: numberOfParticles is not a whole multiple of blocksize.  Best case: some particles aren't checked.  Worst case: undefined.\n";
	
	satelliteDetector <<< numberOfParticles_m / blockSize, blockSize >>> (dblppGPU_m.at(0), dblppGPU_m.at(1), simtime, dt, altitude_m, upwardFacing_m);
}

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, z, time][particle number]
	std::vector<std::vector<double>>& mostRecent{ data_m.at(data_m.size() - 1) };

	for (int satattr = 0; satattr < numberOfAttributes_m + 1; satattr++)
		CUDA_CALL(cudaMemcpy(mostRecent.at(satattr).data(), satCaptureGPU_m + satattr * numberOfParticles_m, sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost));
	
	CUDA_CALL(cudaMemset(satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 1) * numberOfParticles_m)); //sets values to 0

	dataAllocateNewMsmtVector();
}

void Satellite::freeGPUMemory()
{
	CUDA_CALL(cudaFree(satCaptureGPU_m));
	CUDA_CALL(cudaFree(dblppGPU_m.at(1))); //DO NOT FREE dblppGPU_m[0] - this is the 2D data array that the sim uses (not the satellite)
}

/*void Satellite::vectorTest(std::vector<double*>& in)
{
	int wrong{ 0 };
	for (int iii = 0; iii < numberOfAttributes_m; iii++)
		for (int jjj = 0; jjj < numberOfParticles_m; jjj++)
			if (in.at(iii)[jjj] != data_m.at(iii)[jjj]) { wrong++; }

	std::cout << "Wrong: " << wrong << "\n";
}*/