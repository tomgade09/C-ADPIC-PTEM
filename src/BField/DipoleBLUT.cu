#include "BField/DipoleBLUT.h"
#include "BField/DipoleB.h"

#include <memory>

#include "device_launch_parameters.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/cudaDeviceMacros.h"

//setup CUDA kernels
__global__ void setupEnvironmentGPU_DipoleBLUT(BField** this_d, double ILATDeg, double simMin, double simMax, double ds_gradB, int numMsmts, double* altArray, double* magArray)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_DipoleBLUT",
		DipoleBLUT* tmp_d = new DipoleBLUT(ILATDeg, simMin, simMax, ds_gradB, numMsmts);
		tmp_d->setAltArray(altArray);
		tmp_d->setMagArray(magArray);
		(*this_d) = tmp_d;
		);
}

__global__ void deleteEnvironmentGPU_DipoleBLUT(BField** this_d)
{
	ZEROTH_THREAD_ONLY("deleteEnvironmentGPU_DipoleBLUT", delete ((DipoleBLUT*)(*this_d)));
}

__global__ void calcBarray_DipoleBLUT(BField** dipole, double* altitude, double* magnitude, double simMin, double ds)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
	double s{ simMin + ds * thdInd };

	altitude[thdInd] = s;
	magnitude[thdInd] = (*dipole)->getBFieldAtS(s, 0.0);
}
//end

__host__ __device__ DipoleBLUT::DipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_gradB, int numberOfMeasurements) :
	BField("DipoleBLUT"), ILATDegrees_m{ ILATDegrees }, simMin_m{ simMin }, simMax_m{ simMax }, ds_gradB_m{ ds_gradB }, numMsmts_m{ numberOfMeasurements }
{
	ds_msmt_m = (simMax_m - simMin_m) / (numMsmts_m - 1);

	#ifndef __CUDA_ARCH__ //host code
	setupEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ __device__ DipoleBLUT::~DipoleBLUT()
{
	#ifndef __CUDA_ARCH__ //host code
	deleteEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ __device__ double DipoleBLUT::getBFieldAtS(const double s, const double simtime) const
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	int startInd{ 0 };
	if (s <= simMin_m)
		startInd = 0;
	else if (s >= simMax_m)
		startInd = numMsmts_m - 2; //if s is above simMax, we interpolate based on the highest indicies ([numMsmts - 2] to [numMsmts - 1])
	else
		startInd = (int)((s - simMin_m) / ds_msmt_m); //c-style cast to int basically == floor()

	// deltaB_bin / deltas_bin^3 * (s'(dist up from altBin)) + B@altBin
	#ifndef __CUDA_ARCH__ //host code
	return (s - altitude_m.at(startInd)) * (magnitude_m.at(startInd + 1) - magnitude_m.at(startInd)) / ds_msmt_m + magnitude_m.at(startInd); //B = ms + b(0)
	#else
	return (s - altitude_d[startInd]) * (magnitude_d[startInd + 1] - magnitude_d[startInd]) / ds_msmt_m + magnitude_d[startInd]; //B = ms + b(0)
	#endif /* !__CUDA_ARCH__ */
}

__host__ __device__ double DipoleBLUT::getGradBAtS(const double s, const double simtime) const
{
	return (getBFieldAtS(s + ds_gradB_m, simtime) - getBFieldAtS(s - ds_gradB_m, simtime)) / (2 * ds_gradB_m);
}

__host__ __device__ double DipoleBLUT::getSAtAlt(const double alt_fromRe) const
{
	//admittedly, this is a pretty inefficient way of doing this...but...
	return DipoleB(ILATDegrees_m, 1.0e-4, RADIUS_EARTH / 1000.0, false).getSAtAlt(alt_fromRe);
}


//DipoleB class member functions
void DipoleBLUT::setupEnvironment()
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	std::unique_ptr<DipoleB> dip = std::make_unique<DipoleB>(ILATDegrees_m, 1e-10, ds_gradB_m); //destroyed at end of function

	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(BField**)));
	CUDA_API_ERRCHK(cudaMalloc((void **)&altitude_d, sizeof(double) * numMsmts_m));
	CUDA_API_ERRCHK(cudaMalloc((void **)&magnitude_d, sizeof(double) * numMsmts_m));

	int blocksize{ 0 };
	if (numMsmts_m % 256 == 0) //maybe add log entry at this point
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

	calcBarray_DipoleBLUT <<< numMsmts_m / blocksize, blocksize >>> (dip->getPtrGPU(), altitude_d, magnitude_d, simMin_m, ds_msmt_m);
	CUDA_KERNEL_ERRCHK_WSYNC();

	setupEnvironmentGPU_DipoleBLUT <<< 1, 1 >>> (this_d, ILATDegrees_m, simMin_m, simMax_m, ds_gradB_m, numMsmts_m, altitude_d, magnitude_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	//copy generated data back to host into _m arrays
	#ifndef __CUDA_ARCH__ //host code
	altitude_m.resize(numMsmts_m);
	magnitude_m.resize(numMsmts_m);
	CUDA_API_ERRCHK(cudaMemcpy(altitude_m.data(), altitude_d, sizeof(double) * numMsmts_m, cudaMemcpyDeviceToHost));
	CUDA_API_ERRCHK(cudaMemcpy(magnitude_m.data(), magnitude_d, sizeof(double) * numMsmts_m, cudaMemcpyDeviceToHost));
	#endif /* !__CUDA_ARCH__ */
}

void DipoleBLUT::deleteEnvironment()
{
	deleteEnvironmentGPU_DipoleBLUT <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(altitude_d));
	CUDA_API_ERRCHK(cudaFree(magnitude_d));
}