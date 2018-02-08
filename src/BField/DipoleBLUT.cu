#include "BField\DipoleBLUT.h"

__host__ __device__ double DipoleBLUT::getBFieldAtS(const double s, const double simtime)
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	if (s > simMin_m * 0.999)
		return 0.0;

	int startInd{ (int)((s - (simMin_m * 0.999)) / ds_msmt_m) }; //c-style cast to int basically == floor()
	double altBin{ (simMin_m * 0.999 + ds_msmt_m * startInd) };  //or altitude_d[startInd] - not sure which will be faster

	//while ((simMin_m + startInd * ds_msmt_m) > s) // this shouldn't have to be executed - just in case, remove later
		//startInd--; //if location at startInd is over the actual s, subtract one to get the index below
	
	return (s - altBin) * (magnitude_d[startInd + 1] - magnitude_d[startInd]) / (ds_msmt_m); //B = ms + b(0)
}

__host__ __device__ double DipoleBLUT::getGradBAtS(const double s, const double simtime)
{
	return (getBFieldAtS(s + ds_gradB_m, simtime) - getBFieldAtS(s - ds_gradB_m, simtime)) / (2 * ds_gradB_m);
}

//setup CUDA kernels
__global__ void setupEnvironmentGPU_DipoleBLUT(BField** this_d, double ILATDeg, double simMin, double simMax, double ds_gradB, int numMsmts)
{
	if (threadIdx.x == 0 && blockIdx.x == 0)
		(*this_d) = new DipoleBLUT(ILATDeg, simMin, simMax, ds_gradB, numMsmts);
}

__global__ void deleteEnvironmentGPU_DipoleBLUT(BField** this_d)
{
	delete (*this_d);
}

__global__ void calcBarray_DipoleBLUT(BField** dipole, double* altitude, double* magnitude, double simMin, double ds)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
	double s{ simMin +  ds * thdInd };

	altitude[thdInd] = s;
	magnitude[thdInd] = (*dipole)->getBFieldAtS(s, 0.0);
}

//DipoleB class member functions
void DipoleBLUT::setupEnvironment(double ds_dipoleB)
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	std::unique_ptr<DipoleB> dip = std::make_unique<DipoleB>(ILATDegrees_m, 1e-6, ds_dipoleB);

	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(BField**)));
	CUDA_API_ERRCHK(cudaMalloc((void **)&altitude_d, sizeof(double) * numMsmts_m));
	CUDA_API_ERRCHK(cudaMalloc((void **)&magnitude_d, sizeof(double) * numMsmts_m));

	int blocksize{ 0 };
	if (numMsmts_m % 256 == 0)
		blocksize = 256;
	else if (numMsmts_m % 128 == 0)
		blocksize = 128;
	else if (numMsmts_m % 64 == 0)
		blocksize = 64;
	else if (numMsmts_m % 16 == 0)
		blocksize = 16;
	else if (numMsmts_m % 4 == 0)
		blocksize = 4;
	else if (numMsmts_m % 2 == 0)
		blocksize = 2;
	else
		blocksize = 1;

	calcBarray_DipoleBLUT <<< numMsmts_m / blocksize, blocksize >>>(dip->getPtrGPU(), altitude_d, magnitude_d, simMin_m, ds_msmt_m);
	CUDA_KERNEL_ERRCHK_WSYNC();

	setupEnvironmentGPU_DipoleBLUT <<< 1, 1 >>> (this_d, ILATDegrees_m, simMin_m, simMax_m, ds_gradB_m, numMsmts_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void DipoleBLUT::deleteEnvironment()
{
	deleteEnvironmentGPU_DipoleBLUT <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(altitude_d));
	CUDA_API_ERRCHK(cudaFree(magnitude_d));
}