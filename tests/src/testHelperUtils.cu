#include <iostream>

#include "memoryManagement.cuh"

namespace test
{
	void checkGPUMemory(size_t& free, size_t& total) {
		CUDA_API_ERRCHK(cudaMemGetInfo(&free, &total));	}
}