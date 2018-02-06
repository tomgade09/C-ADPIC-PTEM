#include "EField\AlfvenLUT.h"

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "ErrorHandling\cudaErrorCheck.h"

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries);

__device__ double   omegaE_AlfvenLUT;
__device__ double** alfLUT_AlfvenLUT;

__host__ __device__ double EAlfven_AlfvenLUT(const double s, const double simtime)
{//E Field in the direction of B (radially outward)
	if (alfLUT_AlfvenLUT == nullptr)
		return 0.0;
	if (s > RADIUS_EARTH) //in case z is passed in as m, not Re, convert to Re
		s = s / RADIUS_EARTH;
	if (s < alfLUT_AlfvenLUT[0][0] || s > alfLUT_AlfvenLUT[0][2950])
		return 0.0;

	double offset{ alfLUT_AlfvenLUT[0][0] };
	int stepsFrZero{ static_cast<int>(floor((s - offset) / (alfLUT_AlfvenLUT[0][1] - alfLUT_AlfvenLUT[0][0]))) }; //only works for constant bin size - if the binsize changes throughout LUT, need to iterate which will take longer

	//y = mx + b
	double linearInterpReal{ ((alfLUT_AlfvenLUT[1][stepsFrZero + 1] - alfLUT_AlfvenLUT[1][stepsFrZero]) / (alfLUT_AlfvenLUT[0][stepsFrZero + 1] - alfLUT_AlfvenLUT[0][stepsFrZero])) *
		(s - alfLUT_AlfvenLUT[0][stepsFrZero]) + alfLUT_AlfvenLUT[1][stepsFrZero] };
	double linearInterpImag{ ((alfLUT_AlfvenLUT[2][stepsFrZero + 1] - alfLUT_AlfvenLUT[2][stepsFrZero]) / (alfLUT_AlfvenLUT[0][stepsFrZero + 1] - alfLUT_AlfvenLUT[0][stepsFrZero])) *
		(s - alfLUT_AlfvenLUT[0][stepsFrZero]) + alfLUT_AlfvenLUT[2][stepsFrZero] };

	//E-par = (column 2)*cos(omega*t) + (column 3)*sin(omega*t), omega - angular frequency of wave
	return (linearInterpReal * cos(omegaE_AlfvenLUT * simtime) + linearInterpImag * sin(omegaE_AlfvenLUT * simtime)) / 1000.0; //LUT E is in mV / m
}

__global__ void setupEnvironmentGPU_AlfvenLUT(double** LUT, int ind)
{
	if (ind >= MAXEFIELDELEMS)
	{
		fprintf(stderr, "CUDA: __global__ void assignElemToIndex: cannot assign pointer to index specified");
		return;
	}

	EFieldFcnPtrs_EField[ind] = EAlfven_AlfvenLUT;
}



void AlfvenLUT::setupEnvironment()
{
	names_m.push_back("AlfvenLUT");

	int num{ 0 };
	cudaMemcpyFromSymbol(&num, numEFldElems_EField, sizeof(int));

	if (num >= MAXEFIELDELEMS)
		throw SimException("addEFieldElem: too many field elements, cannot add more", __FILE__, __LINE__);

	setupEnvironmentGPU_AlfvenLUT <<< 1, 1 >>> (EFieldLUT2D_d, num);

	num++;
	cudaMemcpyToSymbol(numEFldElems_EField, &num, sizeof(int));
	cudaMemcpyToSymbol(omegaE_AlfvenLUT, &omegaE_m, sizeof(double));
	cudaMemcpyToSymbol(alfLUT_AlfvenLUT, )
}

//GPU kernel setup - 2D array, etc
//Host call GPU kernel setup
//getEFieldAtS(), device, host

/*void AlfvenLUT::initializeFollowOn()
{
	CUDA_CALL(cudaMalloc((void **)&elcFieldLUT1D_d, numOfColsLUT_m * numOfEntrLUT_m * sizeof(double)));
	CUDA_CALL(cudaMalloc((void **)&elcFieldLUT_d,   numOfColsLUT_m * sizeof(double*)));
	CUDA_CALL(cudaMalloc((void **)&omegaE_d, sizeof(double)));
}

void AlfvenLUT::copyDataToGPUFollowOn()
{
	CUDA_CALL(cudaMemcpy(elcFieldLUT1D_d, elcFieldLUT_m[0], numOfColsLUT_m * numOfEntrLUT_m * sizeof(double), cudaMemcpyHostToDevice));
	setup2DArray <<< 1, 1 >>> (elcFieldLUT1D_d, elcFieldLUT_d, numOfColsLUT_m, numOfEntrLUT_m);
	CUDA_CALL(cudaMemcpy(omegaE_d, &omegaE_m, sizeof(double), cudaMemcpyHostToDevice));
}

void AlfvenLUT::iterateSimulationFollowOnPreLoop()
{
	if (useAlfLUT_m && elcFieldLUT_d == nullptr)
		logFile_m.writeLogFileEntry("AlfvenLUT::iterateSimulationFollowOnPreLoop:  Warning: LUT pointer is a nullptr.  Alfven wave function will return 0.0.  Continuing.");
}

void AlfvenLUT::iterateSimulationFollowOnInsideLoop()
{
	return;
}

void AlfvenLUT::iterateSimulationFollowOnPostLoop()
{
	return;
}

void AlfvenLUT::copyDataToHostFollowOn()
{
	return;
}

void AlfvenLUT::freeGPUMemoryFollowOn()
{
	return;
}*/