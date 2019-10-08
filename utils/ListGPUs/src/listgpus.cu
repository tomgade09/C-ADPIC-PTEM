#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>

#define NL {cout << "\n";}

using std::cout;

int main()
{
	int devCount{ 0 };

	if (cudaGetDeviceCount(&devCount) != cudaSuccess)
	{
		cout << "Get Device Count Error: " << cudaGetErrorName(cudaGetLastError()) << "  Exiting.\n";
		return 1;
	}

	for (int devIdx = 0; devIdx < devCount; devIdx++)
	{
		cudaDeviceProp dev;

		if (cudaGetDeviceProperties(&dev, devIdx) != cudaSuccess)
		{
			cout << "Get Device Properties Error: Device Index: " << devIdx << "  Error: " << cudaGetErrorName(cudaGetLastError()) << "  Exiting.\n";
			return 1;
		}

		cout << "================ Device: " << devIdx << " ================"; NL; 
		cout << "\tDevice Name:              " << dev.name; NL;
		cout << "\tTotal Global Mem:         " << dev.totalGlobalMem / 1024 / 1024 / 1024 << " GB"; NL;
		cout << "\tMax Threads Per Block:    " << dev.maxThreadsPerBlock; NL;
		cout << "\tWarp Size:                " << dev.warpSize; NL;
		cout << "\tClock Rate:               " << dev.clockRate / 1024 << " MHz"; NL;
		cout << "\tMemory Clock Rate:        " << dev.memoryClockRate / 1024 << " MHz"; NL;
		cout << "\tMemory Bus Width:         " << dev.memoryBusWidth << " bit"; NL;
		cout << "\tCompute Capability:       " << dev.major << "." << dev.minor; NL;
		NL; NL;
	}

    return 0;
}