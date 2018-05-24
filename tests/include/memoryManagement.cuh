#ifndef TESTS_MEMORYMANAGEMENT_H
#define TESTS_MEMORYMANAGEMENT_H

#include "device_launch_parameters.h"
#include "ErrorHandling/cudaErrorCheck.h" //includes cuda runtime header
#include "ErrorHandling/cudaDeviceMacros.h"

namespace test
{
	template <class BASE>
	class fieldshell;

	template <typename BASE>
	__global__ void makeFSGPU(BASE** this_d, const char* name, double val)
	{
		ZEROTH_THREAD_ONLY("makeFSGPU", (*this_d) = new fieldshell<BASE>(name, val));
	}

	template <typename BASE>
	__global__ void freeFSGPU(BASE** this_d)
	{
		ZEROTH_THREAD_ONLY("freeFSGPU", delete (*this_d));
	}

	template <class BASE>
	class fieldshell : public BASE
	{ //barebones class for testing BField and EElem independent of derived class implementations
	protected:
		__host__ void setupEnvironment() override
		{
			CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(BASE**)));
			makeFSGPU<BASE><<< 1, 1 >>>(this_d, name_m, val_m); //won't compile...need to fix eventually
		}// for now, just won't instantiate a fieldshell class
		__host__ void deleteEnvironment() override
		{
			freeFSGPU<BASE><<< 1, 1 >>>(this_d);
			CUDA_API_ERRCHK(cudaFree(this_d));
		}
		double val_m;
		const char* name_m;

	public:
		__host__ __device__ fieldshell<BASE>(const char* name, double val) : BASE(name), name_m{ name }
		{
			#ifndef __CUDA_ARCH__ //host code
			setupEnvironment();
			#endif
		}

		__host__ __device__ ~fieldshell<BASE>()
		{
			#ifndef __CUDA_ARCH__ //host code
			deleteEnvironment();
			#endif
		}

		__host__ __device__ double getBFieldAtS(const double s, const double t) const { return val_m; }
		__host__ __device__ double getEFieldAtS(const double s, const double t) const { return val_m; }
		__host__ __device__ double getGradBAtS (const double s, const double t) const { return val_m / 10.0; }
	};
}

#endif /* !TESTS_MEMORYMANAGEMENT_H */