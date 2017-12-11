//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h" //didn't add to this vs project - each project this class is attached to will have its own variables header
#include "SatelliteClass\Satellite.h"

__global__ void satelliteDetector(double* v_d, double* mu_d, double* z_d, double* detected_v_d, double* detected_mu_d, double* detected_z_d, double* simtime_d, double simtime, double altitude, bool upward)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	double z_minus_vdt{ z_d[thdInd] - v_d[thdInd] * DT };
	
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
		simtime_d[thdInd] = simtime;
	}//particle not removed from sim
}

void Satellite::initializeSatelliteOnGPU()
{
	for (int iii = 0; iii < numberOfAttributes_m + 1; iii++)
	{//[0] = v_para, [1] = v_perp, [2] = z, [3] = simtime
		cudaMalloc((void **)&captureDataGPU_m[iii], sizeof(double) * numberOfParticles_m); //makes room for data of detected particles
		cudaMemset(captureDataGPU_m[iii], 0, sizeof(double) * numberOfParticles_m); //sets values to 0
	}
}

void Satellite::iterateDetector(int numberOfBlocks, int blockSize, double simtime) {
	satelliteDetector <<< numberOfBlocks, blockSize >>>
		(simDataPtrsGPU_m[0], simDataPtrsGPU_m[1], simDataPtrsGPU_m[2], //Pointers to simulation data on GPU
			captureDataGPU_m[0], captureDataGPU_m[1], captureDataGPU_m[2], captureDataGPU_m[3], //Pointers on GPU to arrays that capture data if the criteria is met
				simtime, altitude_m, upwardFacing_m); } //other variables

void Satellite::copyDataToHost(bool removeZeros)
{// data_m array: [v_para, mu, z, time][particle number]
	if (removeZeros)
	{
		if (data_m != nullptr)
			deleteData();
		allocateData();
	}
	
	for (int iii = 0; iii < numberOfAttributes_m + 1; iii++)
	{
		cudaMemcpy(data_m[iii], captureDataGPU_m[iii], sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost);
		cudaMemset(captureDataGPU_m[iii], 0, sizeof(double) * numberOfParticles_m); //sets values to 0
	}

	if (removeZeros)
	{
		//record indicies of particles captured - throw out the rest
		std::vector<int> ind;
		ind.reserve(numberOfParticles_m);

		for (int iii = 0; iii < numberOfParticles_m; iii++)
			if (data_m[2][iii] > 1)
				ind.push_back(iii);

		//remove zeroes from array
		double** tmp2D = new double*[numberOfAttributes_m + 1];
		for (int iii = 0; iii < numberOfAttributes_m + 1; iii++)
		{
			tmp2D[iii] = new double[ind.size() + 1];
			tmp2D[iii][0] = static_cast<double>(ind.size());
			for (int jjj = 0; jjj < ind.size(); jjj++)
				tmp2D[iii][jjj + 1] = data_m[iii][ind[jjj]]; //index is jjj + 1 because element 0 is the length of the array
		}

		deleteData();
		data_m = tmp2D;
	}
	dataReady_m = true;
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