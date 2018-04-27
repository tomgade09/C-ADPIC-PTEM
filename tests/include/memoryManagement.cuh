#ifndef TESTS_MEMORYMANAGEMENT_H
#define TESTS_MEMORYMANAGEMENT_H

#include "device_launch_parameters.h"
#include "ErrorHandling\cudaErrorCheck.h" //includes cuda runtime header

namespace test
{
	template <class T>
	class fieldshell;

	template <typename T>
	__global__ void makeFSGPU(T** this_d, double val)
	{
		(*this_d) = new fieldshell<T>(val);
	}

	template <typename T>
	__global__ void freeFSGPU(T** this_d)
	{
		delete (*this_d);
	}

	template <class T>
	class fieldshell : public T
	{ //barebones class for testing BField and EElem independent of derived class implementations
	protected:
		__host__ void setupEnvironment() override
		{
			CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(fieldshell<T>**)));
			makeFSGPU<T><<< 1, 1 >>>(this_d, val_m); //won't compile...need to fix eventually
		}// for now, just won't instantiate a fieldshell class
		__host__ void deleteEnvironment() override
		{
			freeFSGPU<T><<< 1, 1 >>>(this_d);
			CUDA_API_ERRCHK(cudaFree(this_d));
		}
		double val_m;

	public:
		__host__ __device__ fieldshell<T>(double val)
		{
			#ifndef __CUDA_ARCH__ //host code
			setupEnvironment();
			#endif
		}

		__host__ __device__ ~fieldshell<T>()
		{
			#ifndef __CUDA_ARCH__ //host code
			deleteEnvironment();
			#endif
		}

		__host__ __device__ double getBFieldAtS(const double s, const double t) const { return val_m; }
		__host__ __device__ double getEFieldAtS(const double s, const double t) const { return val_m; }
		__host__ __device__ double getGradBAtS (const double s, const double t) const { return val_m / 10.0; }

		#ifndef __CUDA_ARCH__ //host code
		//operator overload
		bool operator==(const T& E) const { return true; } //these don't really test anything
		bool operator==(const T* E) const { return true; } //they are just placeholders so the code will compile
		#endif /* !__CUDA_ARCH__ */
	};
}

#endif /* !TESTS_MEMORYMANAGEMENT_H */