#ifndef CUDAERRORCHECK_H
#define CUDAERRORCHECK_H

/*
	Code originally posted at: https://codeyarns.com/2011/03/02/how-to-do-error-checking-in-cuda/
	Modified slightly
*/

#include <iostream>

//CUDA includes
#include "cuda_runtime.h"

#define CUDA_API_ERRCHK( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CUDA_KERNEL_ERRCHK() __cudaCheckError( __FILE__, __LINE__ )
#define CUDA_KERNEL_ERRCHK_WSYNC() __cudaCheckError( __FILE__, __LINE__, true )

inline void __cudaSafeCall(cudaError err, const char* file, const int line)
{
	if (cudaSuccess != err)
		std::cout << file << ":" << line << " : " << "CUDA API error: " << cudaGetErrorString(err) << std::endl;

	return;
}

inline void __cudaCheckError(const char* file, const int line, bool sync=false)
{
	if (sync)
	{
		cudaError err = cudaDeviceSynchronize();
		if (cudaSuccess != err)
			std::cout << file << ":" << line << " : " << "CUDA Kernel error: " << cudaGetErrorString(err) << std::endl;
	}
	else
	{
		cudaError err = cudaGetLastError();
		if (cudaSuccess != err)
			std::cout << "Kernel error: " << file << ":" << line << " : " << cudaGetErrorString(err) << std::endl;
	}

	return;
}

#endif /* CUDAERRORCHECK_H */