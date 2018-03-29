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
#include "physicalconstants.h"
#include "SimulationClass\Simulation.h"
#include "ErrorHandling\cudaErrorCheck.h"
#include "ErrorHandling\SimFatalException.h"

//CUDA Variables - if you change these, don't forget to change the associated curand code/blocks/etc
// For Geforce 960M (author's computer) - maximum 1024 threads per block - try this to see if it results in faster code execution sometime
constexpr int  BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance
constexpr int  NUMRNGSTATES{ 64 * BLOCKSIZE };

//Commonly used values
extern const int SIMCHARSIZE{ 6 * sizeof(double) };

__global__ void initCurand(curandStateMRG32k3a* state, long long seed)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	//if (id > 16380) { printf("Index too high: %i, %i\n", id, sizeof(state[id - 1000])); }
	curand_init(seed, id, 0, &state[id]);
}

__global__ void setup2DArray(double* array1D, double** array2D, int cols, int entries)
{//run once on only one thread
	if (blockIdx.x * blockDim.x + threadIdx.x != 0)
		return;

	for (int iii = 0; iii < cols; iii++)
		array2D[iii] = &array1D[iii * entries];
}

__device__ double accel1dCUDA(const double vs_RK, const double t_RK, const double* args, BField** bfield, EField** efield) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [ps_0, mu, q, m, simtime]
	double F_lor, F_mir, stmp;
	stmp = args[0] + vs_RK * t_RK; //ps_0 + vs_RK * t_RK
	
	//Mirror force
	F_mir = -args[1] * (*bfield)->getGradBAtS(stmp, t_RK + args[4]); //-mu * gradB(pos, runge-kutta time + simtime)
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[2] * (efield == nullptr) ? (0.0) : ((*efield)->getEFieldAtS(stmp, t_RK + args[4])); //q * EFieldatS
	
	return (F_lor + F_mir) / args[3];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(const double y_0, const double h, const double* funcArg, BField** bfield, EField** efield)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	// args array: [t_RKiter, vs, mu, q, m, ps_0, simtime, dt, omega E, const E]
	double k1, k2, k3, k4; double y{ y_0 }; double t_RK{ 0.0 };

	k1 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k1 = f(t_n, y_n), units of dy / dt
	
	t_RK = h / 2;
	y = y_0 + k1 * t_RK;
	k2 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k2 = f(t_n + h/2, y_n + h/2 k1)

	y = y_0 + k2 * t_RK;
	k3 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k3 = f(t_n + h/2, y_n + h/2 k2)

	t_RK = h;
	y = y_0 + k3 * t_RK;
	k4 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns delta y, not dy / dt, not total y
}

__device__ void ionosphereGenerator(double& v_part, double& mu_part, double& s_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{//takes pointers to single particle location in attribute arrays (ex: particle 100: ptr to v[100], ptr to mu[100], ptr to s[100], elecTF, ptr to crndStateA[rnd index of thread]
	double2 v_norm; //v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState); //more efficient to generate two doubles in the one function than run curand_normal_double twice according to CUDA docs
	v_norm.x = v_norm.x * sqrt(simConsts[3] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[3] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	///change to Maxwellian E and pitch angle
	s_part = simConsts[1];
	v_part = abs(v_norm.x);
	mu_part = v_norm.y;
}

__device__ void magnetosphereGenerator(double& v_part, double& mu_part, double& s_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{
	double2 v_norm; //two normal dist values returned to v_norm.x and v_norm.y; v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState);
	v_norm.x = v_norm.x * sqrt(simConsts[4] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[4] * JOULE_PER_EV / mass) + simConsts[5]; //normal dist -> maxwellian
	///change to Maxwellian E and pitch angle
	s_part = simConsts[2];
	v_part = -abs(v_norm.x);
	mu_part = v_norm.y;
}

__global__ void computeKernel(double** currData_d, double** origData_d, double* simConsts, curandStateMRG32k3a* crndStateA, BField** bfield, EField** efield, const double simtime, const double mass, const double charge, const long numParts)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };

	//double v_d = currData_d[0][thdInd]; //not sure why, but switching to kernel-local variables causes a bunch of "unspecified launch failure" errors
	//double mu_d = currData_d[1][thdInd];
	//double s_d = currData_d[2][thdInd];

	double* v_d; double* mu_d; double* s_d;
	v_d = currData_d[0]; mu_d = currData_d[1]; s_d = currData_d[2];
	double* v_orig; double* vperp_orig; double* s_orig;
	v_orig = origData_d[0]; vperp_orig = origData_d[1]; s_orig = origData_d[2];

	if (s_d[thdInd] < 0.001 && simtime == 0.0) //if s is zero (or pretty near zero to account for FP error), generate particles - every other starting at bottom/top of sim
	{
		if (thdInd < numParts / 2) //need perhaps a better way to determine distribution of ionosphere/magnetosphere particles
			ionosphereGenerator(v_d[thdInd], mu_d[thdInd], s_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		else
			magnetosphereGenerator(v_d[thdInd], mu_d[thdInd], s_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		
		v_orig[thdInd] = v_d[thdInd];
		vperp_orig[thdInd] = mu_d[thdInd];
		s_orig[thdInd] = s_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / (*bfield)->getBFieldAtS(s_d[thdInd], simtime);
	}
	else if (simtime == 0.0) //copies data to arrays that track the initial distribution - if data is loaded in, the above block won't be called
	{
		v_orig[thdInd] = v_d[thdInd];
		vperp_orig[thdInd] = mu_d[thdInd];
		s_orig[thdInd] = s_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / (*bfield)->getBFieldAtS(s_d[thdInd], simtime);
	}
	else if (s_d[thdInd] < simConsts[1] * 0.999) //out of sim to the bottom
		return;
	else if (s_d[thdInd] > simConsts[2] * 1.001) //out of sim to the top
		return;

	//args array: [ps_0, mu, q, m, simtime]
	const double args[]{ s_d[thdInd], mu_d[thdInd], charge, mass, simtime, (double)thdInd };

	v_d[thdInd] += foRungeKuttaCUDA(v_d[thdInd], simConsts[0], args, bfield, efield);
	s_d[thdInd] += v_d[thdInd] * simConsts[0];

	//currData_d[0][thdInd] = v_d;
	//currData_d[2][thdInd] = s_d;
}

void Simulation::initializeSimulation()
{
	logFile_m->createTimeStruct("Start Sim Init"); //index 1
	logFile_m->writeTimeDiff(0, 1);

	//size_t free, total;
	//CUDA_API_ERRCHK(cudaMemGetInfo(&free, &total));
	//std::cout << "Pre-Initialize cudaMemGetInfo: free: " << free << ", total: " << total << std::endl;

	//
	//Check for user error
	if (BFieldModel_m == nullptr)
		throw SimFatalException ("initializeSimulation: no Magnetic Field model specified", __FILE__, __LINE__);
	if (particleTypes_m.size() == 0)
		throw SimFatalException ("initializeSimulation: no particles in simulation, sim cannot be initialized without particles", __FILE__, __LINE__);
	//Check for user error complete
	//

	//Allocate room in vectors for GPU Memory Pointers
	satelliteData_m.reserve(100); //not resize...Don't know the exact size here so need to use push_back

	BFieldModel_d = BFieldModel_m->getPtrGPU(); //set pointers that will be passed into CUDA computeKernel
	EFieldModel_d = (EFieldModel_m == nullptr) ? (nullptr) : (EFieldModel_m->getPtrGPU());

	//Allocate memory on GPU for elec/ions variables
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->initializeGPU());

	if (tempSats_m.size() > 0)
	{
		LOOP_OVER_1D_ARRAY(tempSats_m.size(), createSatellite(tempSats_m.at(iii)->particleInd, tempSats_m.at(iii)->altitude, tempSats_m.at(iii)->upwardFacing, tempSats_m.at(iii)->name));
	}
	else
		std::cerr << "Simulation::initializeSimulation: warning: no satellites created" << std::endl;

	CUDA_API_ERRCHK(cudaMalloc((void **)&simConstants_d, SIMCHARSIZE)); //Array of sim characteristics - dt, sim min, sim max, t ion, t mag, v mean
	CUDA_API_ERRCHK(cudaMalloc((void **)&curandRNGStates_d, NUMRNGSTATES * sizeof(curandStateMRG32k3a))); //Array of random number generator states - sizeof(curandStateMRG32k3a) is 72 bytes

	initialized_m = true;
	logFile_m->createTimeStruct("End Sim Init");
	logFile_m->writeTimeDiff(1, 2);

	//some cout statements about the sim here
	//Sim Header printed from Python - move here eventually
	std::cout << "Sim between:    " << simMin_m << "m - " << simMax_m << "m" << std::endl;
	std::cout << "dt:             " << dt_m << "s" << std::endl;
	std::cout << "BField Model:   " << BFieldModel_m->getName() << std::endl;
	std::cout << "EField Elems:   " << ((EFieldModel_m == nullptr) ? ("") : (EFieldModel_m->getEElemsStr())) << std::endl;
	std::cout << "Particles:      " << particleTypes_m.at(0)->getName() << ": #: " << particleTypes_m.at(0)->getNumberOfParticles() << ", loaded files?: " << (particleTypes_m.at(0)->getInitDataLoaded() ? "true" : "false") << std::endl;
	for (int iii = 1; iii < particleTypes_m.size(); iii++) {
		std::cout << "                " << particleTypes_m.at(iii)->getName() << ": #: " << particleTypes_m.at(iii)->getNumberOfParticles() << ", loaded files?: " << (particleTypes_m.at(iii)->getInitDataLoaded() ? "true" : "false") << std::endl; }
	std::cout << "Satellites:     " << satellites_m.at(0)->satellite->getName() << ": alt: " << satellites_m.at(0)->satellite->getAltitude() << " m, upward?: " << (satellites_m.at(0)->satellite->getUpward() ? "true" : "false") << std::endl;
	for (int iii = 1; iii < satellites_m.size(); iii++) {
		std::cout << "                " << satellites_m.at(iii)->satellite->getName() << ": alt: " << satellites_m.at(iii)->satellite->getAltitude() << " m, upward?: " << (satellites_m.at(iii)->satellite->getUpward() ? "true" : "false") << std::endl; }
}

void Simulation::copyDataToGPU()
{//copies particle distribution and associated data to GPU in preparation of iterative calculations over the data
	logFile_m->writeLogFileEntry("Simulation::copyDataToGPU: Start copy to GPU");

	//
	//Check for user error
	if (!initialized_m)
		throw SimFatalException ("copyDataToGPU: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);
	//Check for user error complete
	//

	//Copies initial data of particles to GPU, if loaded
	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->copyDataToGPU());
	
	//Copies array of sim characteristics to GPU - dt, sim min, sim max, t ion, t mag, v mean
	double data[]{ dt_m, simMin_m, simMax_m, ionT_m, magT_m, vmean_m, };
	CUDA_API_ERRCHK(cudaMemcpy(simConstants_d, data, SIMCHARSIZE, cudaMemcpyHostToDevice));
	
	//Prepare curand states for random number generation
	long long seed = time(NULL);
	//std::cout << "seed: " << seed << ", num RNG states: " << NUMRNGSTATES << std::endl;
	initCurand <<< NUMRNGSTATES / 256, 256 >>> (static_cast<curandStateMRG32k3a*>(curandRNGStates_d), seed);
	CUDA_KERNEL_ERRCHK();

	copied_m = true;
	
	logFile_m->writeLogFileEntry("Simulation::copyDataToGPU: End copy to GPU");
}

void Simulation::iterateSimulation(int numberOfIterations, int itersBtwCouts)
{//conducts iterative calculations of data previously copied to GPU - runs the data through the computeKernel
	logFile_m->createTimeStruct("Start Iterate " + std::to_string(numberOfIterations)); //index 3
	logFile_m->writeLogFileEntry("Simulation::iterateSimulation: Start Iteration of Sim:  " + std::to_string(numberOfIterations));

	//
	//Check for user error
	if (!initialized_m)
		throw SimFatalException ("iterateSimulation: sim not initialized with initializeSimulation(), and data not copied to GPU with copyDataToGPU()", __FILE__, __LINE__);
	if (!copied_m)
		throw SimFatalException ("iterateSimulation: data not copied to the GPU with copyDataToGPU()", __FILE__, __LINE__);
	//Check for user error complete
	//

	std::cout << "Iterations:     " << numberOfIterations << std::endl;
	std::cout << "Iters Btw Cout: " << itersBtwCouts << std::endl;
	std::cout << "Time to setup:  "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;
	std::cout << "===============================================================" << std::endl;

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{	
		for (int parts = 0; parts < particleTypes_m.size(); parts++)
		{
			Particle* tmpPart{ particleTypes_m.at(parts).get() };

			computeKernel <<< tmpPart->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> (tmpPart->getCurrDataGPUPtr(), tmpPart->getOrigDataGPUPtr(), simConstants_d,
				static_cast<curandStateMRG32k3a*>(curandRNGStates_d), BFieldModel_d, EFieldModel_d, simTime_m, tmpPart->getMass(), tmpPart->getCharge(), tmpPart->getNumberOfParticles());
		}

		CUDA_KERNEL_ERRCHK_WSYNC_WABORT(); //side effect: cudaDeviceSynchronize() needed for computeKernel to function properly, which this macro provides

		for (int sats = 0; sats < satellites_m.size(); sats++)
			satellites_m.at(sats)->satellite->iterateDetector(simTime_m, dt_m, BLOCKSIZE);

		CUDA_KERNEL_ERRCHK_WSYNC();
		
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
			logFile_m->printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << "\n";
		}

		//if (cudaloopind % 2500 == 0)//need better conditional
			//receiveSatelliteData();
	}
	std::cout << "Total sim time: "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;

	receiveSatelliteData(false);
	std::cout << "\nReceive sat data outside main loop.  Remove after.\n";

	logFile_m->createTimeStruct("End Iterate " + std::to_string(numberOfIterations)); //index 4
	logFile_m->writeTimeDiffFromNow(3, "End Iterate " + std::to_string(numberOfIterations));
	logFile_m->writeLogFileEntry("Simulation::iterateSimulation: End Iteration of Sim:  " + std::to_string(numberOfIterations));
}

void Simulation::copyDataToHost()
{//copies data back to host from GPU
	logFile_m->writeLogFileEntry("Simulation::copyDataToHost: Copy simulation data from GPU back to host");

	//
	//Check for user error
	if (!initialized_m)
		throw SimFatalException("copyDataToGPU: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);
	//Check for user error complete
	//

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->copyDataToHost());
	//maybe add transferring data from satellite over here

	logFile_m->writeLogFileEntry("Simulation::copyDataToHost: Done with copying.");
}

void Simulation::freeGPUMemory()
{//used to free the memory on the GPU that's no longer needed
	logFile_m->writeLogFileEntry("Simulation::freeGPUMemory: Start free GPU Memory.");

	//
	//Check for user error
	if (!initialized_m)
		throw SimFatalException ("freeGPUMemory: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);
	//Check for user error complete
	//

	LOOP_OVER_1D_ARRAY(particleTypes_m.size(), particleTypes_m.at(iii)->freeGPUMemory());
	LOOP_OVER_1D_ARRAY(satellites_m.size(), satellites_m.at(iii)->satellite->freeGPUMemory());

	CUDA_API_ERRCHK(cudaFree(simConstants_d));
	CUDA_API_ERRCHK(cudaFree(static_cast<curandStateMRG32k3a*>(curandRNGStates_d)));

	//size_t free, total;
	//CUDA_API_ERRCHK(cudaMemGetInfo(&free, &total));
	//std::cout << "Post-Free      cudaMemGetInfo: free: " << free << ", total: " << total << std::endl;

	freedGPUMem_m = true;
	logFile_m->writeLogFileEntry("Simulation::freeGPUMemory: End free GPU Memory.");

	CUDA_API_ERRCHK(cudaProfilerStop()); //For profiling with the CUDA bundle
}