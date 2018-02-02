//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"

#include "ParticleClass\Particle.h"

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries);

void Particle::loadFilesToArray(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
		fileIO::readDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", particleCount_m);

	initDataLoaded_m = true;
}

void Particle::saveArrayToFiles(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
		fileIO::writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", particleCount_m);
}

void Particle::normalizeParticles(bool orig, bool curr, bool inverse)
{
	if (!orig && !curr)
		return;
	
	if (normalized_m == true)
		std::cerr << "Particle::normalizeParticles: warning: at least one of the data sets is already normalized" << std::endl;
	if (normFactor_m == 1.0)
		return;

	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
	{
		for (int parts = 0; parts < particleCount_m; parts++)
		{//normalize -> divide by normalization factor
			if (orig) { origData_m.at(attrs).at(parts) /= (inverse ? (1 / normFactor_m) : (normFactor_m)); }
			if (curr) { currData_m.at(attrs).at(parts) /= (inverse ? (1 / normFactor_m) : (normFactor_m)); }
		}
	}

	normalized_m = true;
}

int Particle::getDimensionIndByName(std::string searchName)
{
	for (int name = 0; name < attributeNames_m.size(); name++)
	{
		if (searchName == attributeNames_m.at(name))
			return name;
	}
	//throw exception?
	return -1;
}

std::string Particle::getDimensionNameByInd(int searchIndx)
{
	if (!(searchIndx <= (attributeNames_m.size() - 1) && (searchIndx >= 0)))
		return std::to_string(searchIndx);

	return attributeNames_m.at(searchIndx);
}

//CUDA Functions here
void Particle::initializeGPU()
{
	size_t memSize{ particleCount_m * (getNumberOfAttributes()) * sizeof(double) };
	CUDA_API_ERRCHK(cudaMalloc((void **)&origData1D_d, memSize));
	CUDA_API_ERRCHK(cudaMalloc((void **)&currData1D_d, memSize));
	CUDA_API_ERRCHK(cudaMalloc((void **)&origData2D_d, getNumberOfAttributes() * sizeof(double)));
	CUDA_API_ERRCHK(cudaMalloc((void **)&currData2D_d, getNumberOfAttributes() * sizeof(double)));
	
	CUDA_API_ERRCHK(cudaMemset(origData1D_d, 0, memSize));
	CUDA_API_ERRCHK(cudaMemset(currData1D_d, 0, memSize));

	setup2DArray <<< 1, 1 >>> (origData1D_d, origData2D_d, getNumberOfAttributes(), particleCount_m);
	setup2DArray <<< 1, 1 >>> (currData1D_d, currData2D_d, getNumberOfAttributes(), particleCount_m);

	usingGPU = true;
}

void Particle::copyDataToGPU()
{
	if (!usingGPU)
		throw std::logic_error ("copyDataToGPU: GPU memory has not been initialized yet for particle " + name_m);
	if (!initDataLoaded_m)
		return;

	size_t memSize{ particleCount_m * sizeof(double) };
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_API_ERRCHK(cudaMemcpy(currData1D_d + particleCount_m * iii, currData_m.at(iii).data(), memSize, cudaMemcpyHostToDevice)));
}

void Particle::copyDataToHost()
{
	if (!usingGPU)
		throw std::logic_error ("copyDataToHost: GPU memory has not been initialized yet for particle " + name_m);
	
	size_t memSize{ particleCount_m * sizeof(double) };
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_API_ERRCHK(cudaMemcpy(origData_m.at(iii).data(), origData1D_d + particleCount_m * iii, memSize, cudaMemcpyDeviceToHost)));
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_API_ERRCHK(cudaMemcpy(currData_m.at(iii).data(), currData1D_d + particleCount_m * iii, memSize, cudaMemcpyDeviceToHost)));
}

void Particle::freeGPUMemory()
{
	usingGPU = false;

	CUDA_API_ERRCHK(cudaFree(origData1D_d));
	CUDA_API_ERRCHK(cudaFree(currData1D_d));
	CUDA_API_ERRCHK(cudaFree(origData2D_d));
	CUDA_API_ERRCHK(cudaFree(currData2D_d));
}