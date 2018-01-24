//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ParticleClass\Particle.h"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)
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
{//eventually work in the duplicate function and remove the excess code
	if (!orig && !curr)
		return;
	
	if (normalized_m == true)
	{//need a way to write to log file
		std::cout << "Warning: At least one of the data sets is already normalized.  Make sure neither data set is being normalized twice!" << std::endl;
	}

	if (normFactor_m == 1)
	{
		std::cout << "Warning: Norm factor is 1.  Normalizing will have no effect.  Returning." << std::endl;
		return;
	}

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
	CUDA_CALL(cudaMalloc((void **)&origData1D_d, memSize));
	CUDA_CALL(cudaMalloc((void **)&currData1D_d, memSize));
	CUDA_CALL(cudaMalloc((void **)&origData2D_d, getNumberOfAttributes() * sizeof(double)));
	CUDA_CALL(cudaMalloc((void **)&currData2D_d, getNumberOfAttributes() * sizeof(double)));
	
	CUDA_CALL(cudaMemset(origData1D_d, 0, memSize));
	CUDA_CALL(cudaMemset(currData1D_d, 0, memSize));

	setup2DArray <<< 1, 1 >>> (origData1D_d, origData2D_d, getNumberOfAttributes(), particleCount_m);
	setup2DArray <<< 1, 1 >>> (currData1D_d, currData2D_d, getNumberOfAttributes(), particleCount_m);

	usedGPU = true;
}

void Particle::copyDataToGPU()
{
	if (!usedGPU)
	{
		std::cout << "Error: Particle::copyDataToGPU: Memory has not been initialized yet.  Call Particle::initializeGPU before calling this function.  Returning with no changes." << std::endl;
		return;
	}
	if (!initDataLoaded_m)
		return;

	size_t memSize{ particleCount_m * sizeof(double) };
	/* SHOULD I ALSO COPY ORIG DATA? */
	//LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(origData1D_d + particleCount_m * iii, origData_m.at(iii).data(), memSize, cudaMemcpyHostToDevice)););
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(currData1D_d + particleCount_m * iii, currData_m.at(iii).data(), memSize, cudaMemcpyHostToDevice)););
}

void Particle::copyDataToHost()
{
	if (!usedGPU)
	{
		std::cout << "Error: Particle::copyDataToGPU: Memory has not been initialized yet.  Call Particle::initializeGPU before calling this function.  Returning with no changes." << std::endl;
		return;
	}
	
	size_t memSize{ particleCount_m * sizeof(double) };
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(origData_m.at(iii).data(), origData1D_d + particleCount_m * iii, memSize, cudaMemcpyDeviceToHost)););
	LOOP_OVER_1D_ARRAY(getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(currData_m.at(iii).data(), currData1D_d + particleCount_m * iii, memSize, cudaMemcpyDeviceToHost)););
}

void Particle::freeGPUMemory()
{
	CUDA_CALL(cudaFree(origData1D_d));
	CUDA_CALL(cudaFree(currData1D_d));
	CUDA_CALL(cudaFree(origData2D_d));
	CUDA_CALL(cudaFree(currData2D_d));
}