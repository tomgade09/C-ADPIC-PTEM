//STL includes
#include <iterator> //for back_inserter
#include <algorithm>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

//Project specific includes
#include "Satellite/Satellite.h"
#include "utils/fileIO.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/simExceptionMacros.h"

using utils::fileIO::writeDblBin;
using utils::fileIO::readDblBin;

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries);

__global__ void satelliteDetector(double** data_d, double** capture_d, double simtime, double dt, double altitude, bool upward)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
	
	double* detected_t_d{ capture_d[3] };    //do this first before creating a bunch of pointers
	
	if (simtime == 0.0) //not sure I fully like this, but it works
	{
		double* detected_ind_d{ capture_d[4] };
		detected_t_d[thdInd] = -1.0;
		detected_ind_d[thdInd] = -1.0;
	}
	
	//guard to prevent unnecessary variable creation, if condition checks
	if (detected_t_d[thdInd] > -0.1) return; //if the slot in detected_t[thdInd] is filled (gt or equal to 0), return

	const double* v_d{ data_d[0] };
	const double* s_d{ data_d[2] };
	const double* s0_d{ data_d[5] };

	if (//no detected particle is in the data array at the thread's index already AND
		((!upward) && (s_d[thdInd] >= altitude) && (s0_d[thdInd] < altitude)) //detector is facing down and particle crosses altitude in dt
		  || //OR
		(( upward) && (s_d[thdInd] <= altitude) && (s0_d[thdInd] > altitude)) //detector is facing up and particle crosses altitude in dt
	   )
	{
		const double* mu_d{ data_d[1] };

		double* detected_v_d{ capture_d[0] }; double* detected_mu_d{ capture_d[1] };
		double* detected_s_d{ capture_d[2] }; double* detected_ind_d{ capture_d[4] };

		detected_v_d[thdInd] = v_d[thdInd];
		detected_mu_d[thdInd] = mu_d[thdInd];
		detected_s_d[thdInd] = s_d[thdInd];
		detected_t_d[thdInd] = simtime;
		detected_ind_d[thdInd] = (double)(thdInd);
	}//particle not removed from sim
}

void Satellite::initializeGPU()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&satCaptrData1D_d, sizeof(double) * getNumberOfAttributes() * numberOfParticles_m)); //makes room for data of detected particles
	CUDA_API_ERRCHK(cudaMalloc((void **)&satCaptrData2D_d, sizeof(double*) * getNumberOfAttributes()));
	CUDA_API_ERRCHK(cudaMemset(satCaptrData1D_d, 0, sizeof(double) * getNumberOfAttributes() * numberOfParticles_m)); //sets values to 0

	setup2DArray <<< 1, 1 >>> (satCaptrData1D_d, satCaptrData2D_d, (int)attrNames_m.size(), numberOfParticles_m);
	CUDA_KERNEL_ERRCHK_WSYNC();

	dataOnGPU_m = true;
}

void Satellite::iterateDetectorCPU(const std::vector<std::vector<double>>& partdata, double simtime, double dt) //eventually take a reference to partdata into Satellite??
{
	if (simtime == 0.0 && data_m.size() != 0) data_m.clear(); //if data exists, clear it - temporary work around...
	if (simtime == 0.0 && data_m.size() == 0) dataAllocateNewMsmtVector(); //temporary, if called after 0.0s, this still needs to happen
	if (simtime == 0.0 && !dataReady_m) dataReady_m = true;
	SIM_API_EXCEP_CHECK(satelliteDetectorCPU(partdata, simtime, dt));
}


void Satellite::iterateDetector(double simtime, double dt, int blockSize)
{
	if (numberOfParticles_m % blockSize != 0)
		throw std::invalid_argument("Satellite::iterateDetector: numberOfParticles is not a whole multiple of blocksize, some particles will not be checked");
	
	if (dataReady_m) { dataReady_m = false; }

	satelliteDetector <<< numberOfParticles_m / blockSize, blockSize >>> (particleData2D_d, satCaptrData2D_d, simtime, dt, altitude_m, upwardFacing_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, s, time, partindex][particle number]
	dataAllocateNewMsmtVector();
	
	for (int satattr = 0; satattr < getNumberOfAttributes(); satattr++)
		CUDA_API_ERRCHK(cudaMemcpy(data_m.back().at(satattr).data(), satCaptrData1D_d + satattr * numberOfParticles_m, sizeof(double) * numberOfParticles_m, cudaMemcpyDeviceToHost));
	
	CUDA_API_ERRCHK(cudaMemset(satCaptrData1D_d, 0, sizeof(double) * getNumberOfAttributes() * numberOfParticles_m)); //sets values to 0

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
{//GOING TO HAVE TO REMOVE TO SOMEWHERE - IMPLEMENTATION DEFINED, NOT GENERIC
	if (!dataReady_m)
		copyDataToHost();

	std::vector<std::vector<std::vector<double>>> data_copy{ data_m }; //don't want to do this to the live data so create a copy

	if (removeZeros)
	{
		for (auto& msmt : data_copy) //msmts (vector<vector<double>>) reference
		{
			std::vector<double> timeTmp{ msmt.at(3) }; //make a copy, because array will iterate over index 3 (removing -1.0) before iterating over index 4 (not removing anything)

			for (auto& attr : msmt) //attrs (vector<double>) reference
				attr.erase(std::remove_if(attr.begin(), attr.end(), //below searches time vector for -1.0 (abs within 1e-6) and is true if so (i.e. removes the element)
					[&](double& x) { return ((abs(timeTmp.at(&x - &(*attr.begin())) + 1.0)) <= 1e-6); }), attr.end());
		}
	}

	while (data_copy.size() > 1)
	{
		std::vector<std::vector<double>>& bak{ data_copy.back() };
		for (int attr = 0; attr < bak.size(); attr++) //vector<double>
			data_copy.front().at(attr).insert(data_copy.front().at(attr).end(), bak.at(attr).begin(), bak.at(attr).end());
		data_copy.pop_back();
	}

	return data_copy.at(0);
}

void Satellite::saveDataToDisk(std::string folder) //move B and mass to getConsolidatedData and have it convert back (or in gpu?)
{
	std::vector<std::vector<double>> results{ getConsolidatedData(true) }; //try removing zeroes, fix loadDataFromDisk to "reexpand"

	for (int attr = 0; attr < results.size(); attr++)
		writeDblBin(results.at(attr), folder + name_m + "_" + attrNames_m.at(attr) + ".bin", (int)results.at(attr).size());
}

void Satellite::loadDataFromDisk(std::string folder)
{ //WILL HAVE TO REMOVE EXPANDING PORTION - IMPLEMENTATION DEFINED, NOT GENERIC
	std::vector<std::vector<double>> load(attrNames_m.size());

	for (int attr = 0; attr < attrNames_m.size(); attr++)
		readDblBin(load.at(attr), folder + name_m + "_" + attrNames_m.at(attr) + ".bin");

	bool expand;
	if (load.at(0).size() == numberOfParticles_m)
		expand = false;
	else if (load.at(0).size() < numberOfParticles_m)
		expand = true;
	else
		throw std::logic_error("Satellite::loadDataFromDisk: number of particles loaded from disk is greater than specified numberOfParticles_m.  "
			+ std::string("That means that the wrong data was loaded or wrong number of particles was specified.  Not loading data."));

	if (!expand)
	{
		data_m.push_back(load);
	}
	else
	{
		double prevInd{ -1.0 };
		for (auto index = load.at(4).begin(); index != load.at(4).end(); index++)
		{
			if (*index <= prevInd)
			{
				int len{ (int)(index - load.at(4).begin()) };
				std::vector<std::vector<double>> tmp2D(getNumberOfAttributes());
				for (int attr = 0; attr < tmp2D.size(); attr++) //not an iterator - a vector<double> reference!!
				{
					std::vector<double>& tmp{ tmp2D.at(attr) };
					std::vector<double>& ld{ load.at(attr) };
					tmp.insert(tmp.end(), ld.begin(), ld.begin() + len);
				}
				data_m.push_back(tmp2D); //resets iterator since the preceding data was erased
			}
			else if (index == load.at(4).end() - 1)
			{
				std::vector<std::vector<double>> tmp2D(getNumberOfAttributes());
				for (int attr = 0; attr < tmp2D.size(); attr++)
					tmp2D.at(attr).insert(tmp2D.at(attr).end(), load.at(attr).begin(), load.at(attr).end());
				data_m.push_back(tmp2D);
			}
			prevInd = (*index);
		}

		for (auto& msmt : data_m) //add zeroes to the array where they are missing
		{//measurements
			for (auto attr = msmt.begin(); attr < msmt.end(); attr++) { (*attr).resize(numberOfParticles_m); }
			for (int part = numberOfParticles_m - 1; part >= 0; part--) //particles, iterating backwards
			{
				int ind{ (int)msmt.at(4).at(part) }; //at(4) is "particle index" attribute

				if (ind == 0 && (int)msmt.at(0).at(part) == 0 && (int)msmt.at(1).at(part) == 0)
				{
					msmt.at(3).at(part) = -1;
					msmt.at(4).at(part) = -1;
				}
				else if ((ind != part) && (ind != -1))
				{
					if ((int)msmt.at(0).at(ind) != 0.0) { throw std::runtime_error("FUNCTIONNAME HERE: data is being overwritten in reconstructing array - something is wrong " + std::to_string(ind) + std::to_string((int)msmt.at(0).at(ind))); }
					for (auto& attr : msmt)
					{
						attr.at(ind) = attr.at(part); //move attr at the current location in iteration - part - to the index where it should be - ind
						attr.at(part) = (&attr - &*msmt.begin() < 3) ? 0.0 : -1.0; //overwrite the old data with 0s and -1s
					}
				}
			}
		}
	}
}
