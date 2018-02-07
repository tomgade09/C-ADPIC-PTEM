//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

//Project specific includes
#include "SatelliteClass\Satellite.h"
#include "ErrorHandling\cudaErrorCheck.h"

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries);

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
		index_d[thdInd] = static_cast<double>(thdInd);
	}//particle not removed from sim
}

void Satellite::initializeSatelliteOnGPU()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&satCaptureGPU_m, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //makes room for data of detected particles
	CUDA_API_ERRCHK(cudaMemset(satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //sets values to 0
	CUDA_API_ERRCHK(cudaMalloc((void **)&dblppGPU_m.at(1), sizeof(double*) * (numberOfAttributes_m + 2)));

	setup2DArray <<< 1, 1 >>> (satCaptureGPU_m, dblppGPU_m.at(1), numberOfAttributes_m + 2, numberOfParticles_m, true);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void Satellite::iterateDetector(double simtime, double dt, int blockSize)
{
	if (numberOfParticles_m % blockSize != 0)
		throw std::invalid_argument ("Satellite::iterateDetector: numberOfParticles is not a whole multiple of blocksize, some particles will not be checked");
	
	satelliteDetector <<< numberOfParticles_m / blockSize, blockSize >>> (dblppGPU_m.at(0), dblppGPU_m.at(1), simtime, dt, altitude_m, upwardFacing_m);
}

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, z, time, partindex][particle number]
	dataAllocateNewMsmtVector();
	std::vector<std::vector<double>>& mostRecent{ data_m.at(data_m.size() - 1) };

	for (int satattr = 0; satattr < numberOfAttributes_m + 2; satattr++)
		CUDA_API_ERRCHK(cudaMemcpy(mostRecent.at(satattr).data(), satCaptureGPU_m + satattr * numberOfParticles_m, sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost));
	
	CUDA_API_ERRCHK(cudaMemset(satCaptureGPU_m, 0, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //sets values to 0

	dataReady_m = true; //sets to true the first time called
}

void Satellite::freeGPUMemory()
{
	if (!dataOnGPU_m)
		return;

	dataOnGPU_m = false;
	CUDA_API_ERRCHK(cudaFree(satCaptureGPU_m));
	CUDA_API_ERRCHK(cudaFree(dblppGPU_m.at(1))); //DO NOT FREE dblppGPU_m.at(0) - this is the 2D data array that the sim uses (not the satellite)
}

std::vector<std::vector<double>> Satellite::getConsolidatedData(bool removeZeros)
{
	if (!dataReady_m)
		copyDataToHost();

	std::vector<std::vector<double>> tmp2D;

	for (int attrs = 0; attrs < numberOfAttributes_m + 2; attrs++)
		tmp2D.push_back(std::vector<double>());

	LOOP_OVER_3D_ARRAY(data_m.size(), data_m.at(iii).size(), numberOfParticles_m, \
		if (removeZeros) //iii is msmt iterator, jjj is attribute iterator, kk is particle iterator
		{
			size_t tind{ data_m.at(iii).size() - 1 };
			if (data_m.at(iii).at(tind).at(kk) >= 0.0)
				tmp2D.at(jjj).push_back(data_m.at(iii).at(jjj).at(kk));
		}
		else
			tmp2D.at(jjj).push_back(data_m.at(iii).at(jjj).at(kk));
	)

		return tmp2D;
}