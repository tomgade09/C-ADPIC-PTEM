//Standard Library includes
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "_simulationvariables.h"
#include "SimulationClass\Simulation.h"

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)
__device__ double getBFieldAtS(double s, double t);
__device__ double getGradBAtS(double s, double t);
//__device__ double getEFieldAtS(double s, double t);

//CUDA Variables - if you change these, don't forget to change the associated curand code/blocks/etc
// For Geforce 960M (author's computer) - maximum 1024 threads per block - try this to see if it results in faster code execution sometime
constexpr int  BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance
constexpr int  NUMRNGSTATES{ 64 * BLOCKSIZE };

//Commonly used values
extern const int SIMCHARSIZE{ 6 * sizeof(double) };

__device__ double accel1dCUDA(double* args, int len) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RKiter, vs, mu, q, m, ps_0, simtime, dt, omega E, const E]
	double F_lor, F_mir, stmp;
	stmp = args[5] + args[1] * args[0]; //ps_0 + vs * t_RK
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * 0.0; //getEFieldAtS(some arguments, etc)

	//Mirror force
	F_mir = -args[2] * getGradBAtS(stmp, args[0] + args[6]); //-mu * gradB

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(double* funcArg, int arrayLen)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	//args array: [t_RKiter, vs, mu, q, m, ps_0, simtime, dt, omega E, const E]
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = accel1dCUDA(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt
	
	funcArg[0] = funcArg[7] / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = funcArg[7];
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * funcArg[7] / 6; //returns delta y, not dy / dt, not total y
}

__global__ void initCurand(curandStateMRG32k3a* state, long long seed)
{
	long long id = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(seed, id, 0, &state[id]);
}

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries)
{//run once on only one thread
	if (blockIdx.x * blockDim.x + threadIdx.x != 0)
		return;
	
	for (int iii = 0; iii < cols; iii++)
		array2D[iii] = &array1D[iii * entries];
}

__device__ void ionosphereGenerator(double* v_part, double* mu_part, double* s_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{//takes pointers to single particle location in attribute arrays (ex: particle 100: ptr to v[100], ptr to mu[100], ptr to s[100], elecTF, ptr to crndStateA[rnd index of thread]
	double2 v_norm; //v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState); //more efficient to generate two doubles in the one function than run curand_normal_double twice according to CUDA docs
	v_norm.x = v_norm.x * sqrt(simConsts[3] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[3] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	///change to Maxwellian E and pitch angle
	*s_part = simConsts[1];
	*v_part = abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void magnetosphereGenerator(double* v_part, double* mu_part, double* s_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{
	double2 v_norm; //two normal dist values returned to v_norm.x and v_norm.y; v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState);
	v_norm.x = v_norm.x * sqrt(simConsts[4] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[4] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	///change to Maxwellian E and pitch angle
	*s_part = simConsts[2];
	*v_part = -abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void ionosphereScattering(double* v_part, double* mu_part, double* s_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{
	//
	//
	//
	//Physics needs to be improved
	ionosphereGenerator(v_part, mu_part, s_part, simConsts, mass, rndState);
}

__global__ void computeKernel(double** currData_d, double** origData_d, double* simConsts, curandStateMRG32k3a* crndStateA,	double simtime, double mass, double charge, long numParts)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };

	double* v_d; double* mu_d; double* s_d;
	double* v_orig; double* vperp_orig; double* s_orig;
	v_d = currData_d[0]; mu_d = currData_d[1]; s_d = currData_d[2];
	v_orig = origData_d[0]; vperp_orig = origData_d[1]; s_orig = origData_d[2];

	if (s_d[thdInd] < 0.001) //if s is zero (or pretty near zero to account for FP error), generate particles - every other starting at bottom/top of sim
	{//previous way to index curandStates: (blockIdx.x * 2) + (threadIdx.x % 2) - this leads to each block accessing two curand states - 128 threads call the same state simultaneously and end up with the same values
		if (thdInd < numParts / 2) //need perhaps a better way to determine distribution of ionosphere/magnetosphere particles
			ionosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &s_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		else
			magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &s_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		
		v_orig[thdInd] = v_d[thdInd];
		vperp_orig[thdInd] = mu_d[thdInd];
		s_orig[thdInd] = s_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / getBFieldAtS(s_d[thdInd], simtime);
	}
	else if (simtime == 0) //copies data to arrays that track the initial distribution - if data is loaded in, the above block won't be called
	{
		v_orig[thdInd] = v_d[thdInd];
		vperp_orig[thdInd] = mu_d[thdInd];
		s_orig[thdInd] = s_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / getBFieldAtS(s_d[thdInd], simtime);
	}
	else if (s_d[thdInd] < simConsts[1] * 0.999) //out of sim to the bottom, particle has 50% chance of reflecting, 50% chance of new particle
		//ionosphereScattering(&v_d[thdInd], &mu_d[thdInd], &s_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	else if (s_d[thdInd] > simConsts[2] * 1.001) //out of sim to the top, particle is lost, new one generated in its place
		//magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &s_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	
	//args array: [t_RKiter, vs, mu, q, m, ps_0, simtime, dt, omega E, const E]
	//simConsts: dt, sim min, sim max, t ion, t mag, v mean, omega E Alfven, QSPS const E
	double args[10];
	args[0] = 0.0;
	args[1] = v_d[thdInd]; //velocity along s
	args[2] = mu_d[thdInd]; //mu
	args[3] = charge;
	args[4] = mass;
	args[5] = s_d[thdInd]; //position along s
	args[6] = simtime;
	args[7] = simConsts[0]; //dt

	v_d[thdInd] += foRungeKuttaCUDA(args, 8);
	s_d[thdInd] += v_d[thdInd] * simConsts[0];
}

void Simulation::initializeSimulation()
{	
	logFile_m.createTimeStruct("Start Sim Init"); //index 1
	logFile_m.writeTimeDiff(0, 1);

	//
	//Check for user error
	if (particleTypes_m.size() == 0)
	{
		logFile_m.writeErrorEntry("Simulation::initializeSimulation", "No particles in sim.  You need to add particles before calling this function.  Returning.", {});
		errorEncountered = true;
		return;
	}

	if (errorEncountered)
		return;
	//Check for user error complete
	//

	//Allocate room in vectors for GPU Memory Pointers
	satelliteData_m.reserve(100); //not resize...Don't know the exact size here so need to use push_back

	//Allocate memory on GPU for elec/ions variables
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->initializeGPU(););

	if (tempSats_m.size() > 0)
	{
		LOOP_OVER_1D_ARRAY(tempSats_m.size(), createSatellite(tempSats_m.at(iii)->particleInd, tempSats_m.at(iii)->altitude, tempSats_m.at(iii)->upwardFacing, tempSats_m.at(iii)->name););
	}
	else
		logFile_m.writeLogFileEntry("Warning: Simulation::initializeSimulation: No satellites created.  That's odd.");

	//Array of sim characteristics - dt, sim min, sim max, t ion, t mag, v mean, omega E Alfven, QSPS const E
	CUDA_CALL(cudaMalloc((void **)&simConstants_d, SIMCHARSIZE)); //for right now, fine

	//Array of random number generator states
	CUDA_CALL(cudaMalloc((void **)&curandRNGStates_d, NUMRNGSTATES * sizeof(curandStateMRG32k3a))); //sizeof(curandStateMRG32k3a) is 72 bytes

	//For derived classes to add code
	initializeFollowOn(); //need to remove these

	initialized_m = true;
	logFile_m.createTimeStruct("End Sim Init");
	logFile_m.writeTimeDiff(1, 2);
}

void Simulation::copyDataToGPU()
{//copies particle distribution and associated data to GPU in preparation of iterative calculations over the data
	logFile_m.writeLogFileEntry("Simulation::copyDataToGPU: Start copy to GPU");

	//
	//Check for user error
	if (!initialized_m)
	{
		logFile_m.writeErrorEntry("Simulation::copyDataToGPU", "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.", {});
		errorEncountered = true;
		return;
	}

	if (errorEncountered)
		return;
	//Check for user error complete
	//

	//Copies initial data of particles to GPU, if loaded
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->copyDataToGPU(););
	
	
	//Copies array of sim characteristics to GPU - dt, sim min, sim max, t ion, t mag, v mean
	double data[]{ dt_m, simMin_m, simMax_m, ionT_m, magT_m, vmean_m, };
	CUDA_CALL(cudaMemcpy(simConstants_d, data, SIMCHARSIZE, cudaMemcpyHostToDevice));
	
	//Prepare curand states for random number generation
	long long seed = time(NULL);
	initCurand <<< NUMRNGSTATES / 256, 256 >>> (static_cast<curandStateMRG32k3a*>(curandRNGStates_d), seed);
	
	//For derived classes to add code
	copyDataToGPUFollowOn(); //need to get rid of later

	copied_m = true;
	
	logFile_m.writeLogFileEntry("Simulation::copyDataToGPU: End copy to GPU");
}

void Simulation::iterateSimulation(int numberOfIterations, int itersBtwCouts)
{//conducts iterative calculations of data previously copied to GPU - runs the data through the computeKernel
	logFile_m.createTimeStruct("Start Iterate " + std::to_string(numberOfIterations)); //index 3
	logFile_m.writeLogFileEntry("Simulation::iterateSimulation: Start Iteration of Sim:  " + std::to_string(numberOfIterations));
	
	//
	//Check for user error
	if (!initialized_m)
	{
		logFile_m.writeErrorEntry("Simulation::iterateSimulation", "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.  You also need to copy data to the GPU with Simulation::copyDataToGPU.", { std::to_string(numberOfIterations) });
		errorEncountered = true;
		return;
	}
	if (!copied_m)
	{
		logFile_m.writeErrorEntry("Simulation::iterateSimulation", "You haven't copied any data to the GPU with Simulation::copyDataToGPU.  Do that first or the GPU has no numbers to work on.", { std::to_string(numberOfIterations) });
		errorEncountered = true;
		return;
	}
	if (errorEncountered)
		return;
	//Check for user error complete
	//

	//For derived classes to add code
	iterateSimulationFollowOnPreLoop(); //remove this later

	size_t numParts{ particleTypes_m.size() };

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{	
		for (int parts = 0; parts < particleTypes_m.size(); parts++)
		{
			Particle* tmpPart{ particleTypes_m.at(parts) };

			computeKernel <<< tmpPart->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> (tmpPart->getCurrDataGPUPtr(), tmpPart->getOrigDataGPUPtr(), simConstants_d,
				static_cast<curandStateMRG32k3a*>(curandRNGStates_d),	simTime_m, tmpPart->getMass(), tmpPart->getCharge(), tmpPart->getNumberOfParticles());

			cudaDeviceSynchronize();
		}

		for (int sats = 0; sats < satellites_m.size(); sats++)
			satellites_m.at(sats)->satellite->iterateDetector(BLOCKSIZE, simTime_m, dt_m);
		
		cudaloopind++;
		incTime();

		if (cudaloopind % itersBtwCouts == 0)
		{
			std::stringstream out;
			out << std::setw(std::to_string(numberOfIterations).size()) << cudaloopind; //not sure if I like the setw(std::to_string(blah)) solution...
			std::cout << out.str() << " / " << numberOfIterations << "  |  Sim Time (s): ";
			out.str(""); out.clear();
			out << std::setw(std::to_string(static_cast<double>(numberOfIterations) * dt_m).size()) << std::fixed << simTime_m; //not sure if I like the setw(std::to_string(blah)) solution...
			std::cout << out.str() << "  |  Real Time Elapsed (s): ";
			logFile_m.printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << "\n";
		}

		//if (cudaloopind % 2500 == 0)//need better conditional
			//receiveSatelliteData();

		//For derived classes to add code
		iterateSimulationFollowOnInsideLoop(); //remove this eventually
	}
	receiveSatelliteData(false);
	std::cout << "\nReceive sat data outside main loop.  Remove after.\n";

	//For derived classes to add code
	iterateSimulationFollowOnPostLoop(); //remove later

	logFile_m.createTimeStruct("End Iterate " + std::to_string(numberOfIterations)); //index 4
	logFile_m.writeTimeDiffFromNow(3, "End Iterate " + std::to_string(numberOfIterations));
	logFile_m.writeLogFileEntry("Simulation::iterateSimulation: End Iteration of Sim:  " + std::to_string(numberOfIterations));
}

void Simulation::copyDataToHost()
{//copies data back to host from GPU
	logFile_m.writeLogFileEntry("Simulation::copyDataToHost: Copy simulation data from GPU back to host");
	
	//
	//Check for user error
	if (!initialized_m)
	{
		logFile_m.writeErrorEntry("Simulation::copyDataToHost", "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.", {});
		errorEncountered = true;
		return;
	}
	
	if (errorEncountered)
		return;
	//Check for user error complete
	//

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->copyDataToHost(););

	//For derived classes to add code
	copyDataToHostFollowOn(); //remove later

	logFile_m.writeLogFileEntry("Simulation::copyDataToHost: Done with copying.");
}

void Simulation::freeGPUMemory()
{//used to free the memory on the GPU that's no longer needed
	logFile_m.writeLogFileEntry("Simulation::freeGPUMemory: Start free GPU Memory.");

	//
	//Check for user error
	if (!initialized_m)
	{
		logFile_m.writeErrorEntry("Simulation::freeGPUMemory", "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.", {});
		return;
	}
	//Check for user error complete
	//

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->freeGPUMemory(););

	//For derived classes to add code
	freeGPUMemoryFollowOn(); //remove later

	freedGPUMem_m = true;

	logFile_m.writeLogFileEntry("Simulation::freeGPUMemory: End free GPU Memory.");

	CUDA_CALL(cudaProfilerStop()); //For profiling the profiler in the CUDA bundle}
}