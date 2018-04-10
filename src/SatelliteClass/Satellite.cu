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
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
	
	double* v_d; double* mu_d; double* s_d; double* simtime_d; double* index_d;
	double* detected_v_d; double* detected_mu_d; double* detected_s_d;
	v_d = data_d[0]; mu_d = data_d[1]; s_d = data_d[2]; simtime_d = capture_d[3]; index_d = capture_d[4];
	detected_v_d = capture_d[0]; detected_mu_d = capture_d[1]; detected_s_d = capture_d[2];

	double s_minus_vdt{ s_d[thdInd] - v_d[thdInd] * dt };
	
	if (simtime == 0.0) //not sure I fully like this, but it works
	{
		simtime_d[thdInd] = -1.0;
		index_d[thdInd] = -1.0;
	}

	if (//if a particle is detected, the slot in array[thdInd] is filled and no further detection happens for that index
		(simtime_d[thdInd] < -0.1) &&
			( //no detected particle is in the data array at the thread's index already AND
				//detector is facing down and particle crosses altitude in dt
				((!upward) && (s_d[thdInd] >= altitude) && (s_minus_vdt < altitude))
					|| //OR
				//detector is facing up and particle crosses altitude in dt
				((upward) && (s_d[thdInd] <= altitude) && (s_minus_vdt > altitude))
			)
		)
	{
		detected_v_d[thdInd] = v_d[thdInd];
		detected_mu_d[thdInd] = mu_d[thdInd];
		detected_s_d[thdInd] = s_d[thdInd];
		simtime_d[thdInd] = simtime;
		index_d[thdInd] = (double)(thdInd);
	}//particle not removed from sim
}

void Satellite::initializeGPU()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&satCaptrData1D_d, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //makes room for data of detected particles
	CUDA_API_ERRCHK(cudaMalloc((void **)&satCaptrData2D_d, sizeof(double*) * (numberOfAttributes_m + 2)));
	CUDA_API_ERRCHK(cudaMemset(satCaptrData1D_d, 0, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //sets values to 0

	setup2DArray <<< 1, 1 >>> (satCaptrData1D_d, satCaptrData2D_d, numberOfAttributes_m + 2, numberOfParticles_m);
	CUDA_KERNEL_ERRCHK_WSYNC();

	dataOnGPU_m = true;
}

void Satellite::iterateDetector(double simtime, double dt, int blockSize)
{
	if (numberOfParticles_m % blockSize != 0)
		throw std::invalid_argument("Satellite::iterateDetector: numberOfParticles is not a whole multiple of blocksize, some particles will not be checked");
	
	satelliteDetector <<< numberOfParticles_m / blockSize, blockSize >>> (particleData2D_d, satCaptrData2D_d, simtime, dt, altitude_m, upwardFacing_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, z, time, partindex][particle number]
	dataAllocateNewMsmtVector();
	std::vector<std::vector<double>>& mostRecent{ data_m.at(data_m.size() - 1) };

	for (int satattr = 0; satattr < numberOfAttributes_m + 2; satattr++)
		CUDA_API_ERRCHK(cudaMemcpy(mostRecent.at(satattr).data(), satCaptrData1D_d + satattr * numberOfParticles_m, sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost));
	
	CUDA_API_ERRCHK(cudaMemset(satCaptrData1D_d, 0, sizeof(double) * (numberOfAttributes_m + 2) * numberOfParticles_m)); //sets values to 0

	dataReady_m = true; //sets to true the first time called
}

void Satellite::freeGPUMemory()
{
	if (!dataOnGPU_m) { return; }

	CUDA_API_ERRCHK(cudaFree(satCaptrData1D_d));
	CUDA_API_ERRCHK(cudaFree(satCaptrData2D_d));
	//DO NOT FREE particleData2D_d - this is the 2D data array that the sim uses (not the satellite)

	satCaptrData1D_d = nullptr;
	satCaptrData2D_d = nullptr;
	particleData2D_d = nullptr;

	dataOnGPU_m = false;
}

std::vector<std::vector<double>> Satellite::getConsolidatedData(bool removeZeros)
{//redo this/look into at some point
	if (!dataReady_m)
		copyDataToHost();

	std::vector<std::vector<std::vector<double>>> data_copy{ data_m };

	if (removeZeros)
	{
		for (int iii = 0; iii < data_copy.size(); iii++) //msmts
		{
			std::vector<double> timeTmp{ data_copy.at(iii).at(3) }; //make a copy, because array will iterate over index 3 (removing -1.0) before iterating over index 4 (not removing anything)

			for (int jjj = 0; jjj < data_copy.at(iii).size(); jjj++)
			{
				std::vector<double>& tmp{ data_copy.at(iii).at(jjj) };
				tmp.erase(std::remove_if(tmp.begin(), tmp.end(), //below searches time vector for -1.0 (abs within 1e-6) and is true if so (i.e. removes the element)
					[&timeTmp, &tmp](double& x) { return ((abs(timeTmp.at(&x - &*tmp.begin()) + 1.0)) <= 1e-6); }), tmp.end());				
			}
		}
	}

	int sumLen{ 0 };
	for (auto msmt = data_copy.begin(); msmt < data_copy.end(); msmt++) //msmts
		sumLen += (int)(*msmt).at(0).size();

	std::vector<std::vector<double>> ret(numberOfAttributes_m + 2, std::vector<double>(sumLen));

	//append elements
	int offset{ 0 };
	for (int iii = 0; iii < data_copy.size(); iii++) //msmts
	{
		for (int jjj = 0; jjj < data_copy.at(iii).size(); jjj++) //attrs
			std::copy(data_copy.at(iii).at(jjj).begin(), data_copy.at(iii).at(jjj).end(), ret.at(jjj).begin() + offset);
		offset += data_copy.at(iii).at(0).size();
	}

	return ret;
}

void Satellite::saveDataToDisk(std::string folder, std::vector<std::string> attrNames) //move B and mass to getConsolidatedData and have it convert back (or in gpu?)
{
	std::vector<std::vector<double>> results{ getConsolidatedData(false) }; //try removing zeroes, fix loadDataFromDisk to "reexpand"

	for (int attr = 0; attr < results.size(); attr++)
		fileIO::writeDblBin(results.at(attr), folder + name_m + "_" + attrNames.at(attr) + ".bin", (int)results.at(attr).size());
}

void Satellite::loadDataFromDisk(std::string folder, std::vector<std::string> attrNames)
{
	std::vector<std::vector<double>> load(attrNames.size());

	for (int attr = 0; attr < attrNames.size(); attr++)
		fileIO::readDblBin(load.at(attr), folder + name_m + "_" + attrNames.at(attr) + ".bin");

	int offset{ 0 };
	for (int iii = 1; iii < load.at(4).size(); iii++) //split condensed data into measurements
	{
		if ((int)load.at(4).at(iii) <= (int)load.at(4).at(iii - 1) || iii == load.at(4).size() - 1)
		{
			std::vector<std::vector<double>> tmp2D;
			std::vector<double> tmp(iii - offset + ((iii == load.at(4).size() - 1) ? 1 : 0)); //appropriately sized array of zeroes
			for (auto attr = load.begin(); attr < load.end(); attr++)
			{
				std::copy((*attr).begin() + offset, ((iii == load.at(4).size() - 1) ? (*attr).end() : (*attr).begin() + iii), tmp.begin());
				tmp2D.push_back(tmp);
			}
			data_m.push_back(tmp2D);
			offset = iii;
		}
	}

	for (auto msmt = data_m.begin(); msmt < data_m.end(); msmt++) //add zeroes to the array where they are missing
	{//measurements
		for (auto attr = (*msmt).begin(); attr < (*msmt).end(); attr++) { (*attr).resize(numberOfParticles_m); }
		for (int part = numberOfParticles_m - 1; part >= 0; part--) //particles
		{
			int ind{ (int)(*msmt).at(4).at(part) }; //at(4) is "particle index" attribute

			if (ind == 0 && (int)(*msmt).at(0).at(part) == 0 && (int)(*msmt).at(1).at(part) == 0)
			{
				(*msmt).at(3).at(part) = -1;
				(*msmt).at(4).at(part) = -1;
			}
			else if ((ind != part) && (ind != -1))
			{
				if ((int)((*msmt).at(0).at(ind) != 0.0)) { throw std::runtime_error("FUNCTIONNAME HERE: data is being overwritten in reconstructing array - something is wrong " + std::to_string(ind) + std::to_string((int)(*msmt).at(0).at(ind))); }
				for (auto attr = (*msmt).begin(); attr < (*msmt).end(); attr++)
				{
					(*attr).at(ind) = (*attr).at(part);
					(*attr).at(part) = (attr - (*msmt).begin() < 3) ? 0.0 : -1.0;
				}
			}
		}
	}
}