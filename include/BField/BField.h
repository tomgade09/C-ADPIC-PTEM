#ifndef BFIELD_H
#define BFIELD_H

typedef double(*callbackFcn)(double*, int);

class BField
{
protected:
	//callback function pointer
	callbackFcn BFieldFcnPtr_m;
	callbackFcn gradBFcnPtr_m;
	//assign with: BFieldFcnPtr = getBatSt; call with: result = BFieldFcnPtr(args, count); etc

	//callback function GPU pointer
	/*I would imagine the process is like so:
	  - Create pointer on GPU with cudaMalloc
	  - Run CUDA kernel to assign pointer to the function to the memory location
	  - Pass allocated pointer to CUDA kernel responsible for running the callback
	  - Run as normal within the CUDA kernel
	*/
	callbackFcn BFieldFcnPtr_d;
	callbackFcn gradBFcnPtr_d;

	virtual void setupCallbacksonGPU();
	virtual void callSetupCallbacksKernel() = 0; //define this function in derived classes to assign a pointer to that function's B Field code to the location indicated by BFieldFcnPtr_d and gradBFcnPtr_d

public:
	BField() {}
	~BField() {}
};

#endif