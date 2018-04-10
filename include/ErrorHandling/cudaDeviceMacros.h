#ifndef CUDADEVICEMACROS_H
#define CUDADEVICEMACROS_H

//guard macros meant to execute on device
#define ZEROTH_THREAD_ONLY(i, j) if (threadIdx.x == 0 && blockIdx.x == 0) { j; } \
	else { printf("Trying to run %s on multiple threads.  These kernels are not meant for this.  Use <<< 1, 1 >>> only!", i); }

#endif