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
#include "include\_simulationvariables.h"
#include "include\Simulation170925.h"

//Array Size Variables
//constexpr int DBLARRAY_BYTES { NUMPARTICLES * sizeof(double) };
//constexpr int BOOLARRAY_BYTES{ NUMPARTICLES * sizeof(bool) };
//constexpr int INTARRAY_BYTES { NUMPARTICLES * sizeof(int) };

//Commonly used values
constexpr double EOMEGA{ 20 * PI }; //10 Hz wave, matches ez.out
constexpr double B_AT_MIN_Z{ B0ATTHETA / (MIN_Z_NORM * MIN_Z_NORM * MIN_Z_NORM) };
constexpr double B_AT_MAX_Z{ B0ATTHETA / (MAX_Z_NORM * MAX_Z_NORM * MAX_Z_NORM) };

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { printf("Error %d at %s:%d\n",EXIT_FAILURE,__FILE__,__LINE__);}} while(0)

__host__ __device__ double qspsEatZ(double z, double simtime)
{
	//if ((z > E_RNG_CENTER + E_RNG_DELTA) || (z < E_RNG_CENTER - E_RNG_DELTA))
		//return 0.0;
	return CONSTEFIELD;
}

//merge these two below...
__host__ __device__ double alfvenWaveEbyLUT(double** LUT, double z, double simtime)
//__host__ double EFieldatZ(double** LUT, double z, double simtime)
{//E Field in the direction of B (radially outward)
	//E-par = (column 2)*cos(omega*t) + (column 3)*sin(omega*t), omega = 20 PI
	if (z > RADIUS_EARTH) //in case z is passed in as m, not Re
		z = z / RADIUS_EARTH; //convert to Re
	
	if (z < LUT[0][0] || z > LUT[0][2950])
		return 0.0;
	double offset{ LUT[0][0] };
	int stepsFromZeroInd{ static_cast<int>(floor((z - offset) / (LUT[0][1] - LUT[0][0]))) }; //only works for constant bin size - if the binsize changes throughout LUT, need to iterate which will take longer
	
	//y = mx + b
	double linearInterpReal{ ((LUT[1][stepsFromZeroInd + 1] - LUT[1][stepsFromZeroInd]) / (LUT[0][stepsFromZeroInd + 1] - LUT[0][stepsFromZeroInd])) *
		(z - LUT[0][stepsFromZeroInd]) + LUT[1][stepsFromZeroInd] };
	double linearInterpImag{ ((LUT[2][stepsFromZeroInd + 1] - LUT[2][stepsFromZeroInd]) / (LUT[0][stepsFromZeroInd + 1] - LUT[0][stepsFromZeroInd])) *
		(z - LUT[0][stepsFromZeroInd]) + LUT[2][stepsFromZeroInd] };
	
	return (linearInterpReal * cos(EOMEGA * simtime) + linearInterpImag * sin(EOMEGA * simtime)) / 1000; //LUT E is in mV / m
}

__host__ __device__ double EFieldatZ(double** LUT, double z, double simtime)
{
	return ((false) ? (qspsEatZ(z, simtime)) : (0.0)) + ((false) ? (alfvenWaveEbyLUT(LUT, z, simtime)) : (0.0));
}

__host__ __device__ double BFieldatZ(double z, double simtime)
{//for now, a simple dipole field
	if (z == 0)
		return 0.0; //add an error here if this case is true, at some point

	double norm{ RADIUS_EARTH };

	if ((z < MAX_Z_NORM * 1.001 ) && (z > MIN_Z_NORM * 0.999))
		norm = 1.0;

	return B0ATTHETA / pow(z / norm, 3); //Bz = B0 at theta * (1/rz(in Re))^3
}

__device__ double accel1dCUDA(double* args, int len, double** LUT) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	double F_lor, F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * EFieldatZ(LUT, ztmp, args[6] + args[0]); //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * B0ATTHETA * (-3 / (pow(ztmp / RADIUS_EARTH, 4))) / RADIUS_EARTH; //mu in [kg.m^2 / s^2.T] = [N.m / T]

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(double* funcArg, int arrayLen, double** LUT)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = accel1dCUDA(funcArg, arrayLen, LUT); //k1 = f(t_n, y_n), units of dy / dt
	
	funcArg[0] = DT / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen, LUT); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen, LUT); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = DT;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen, LUT); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * DT / 6; //returns delta y, not dy / dt, not total y
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

__device__ void ionosphereGenerator(double* v_part, double* mu_part, double* z_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{//takes pointers to single particle location in attribute arrays (ex: particle 100: ptr to v[100], ptr to mu[100], ptr to z[100], elecTF, ptr to crndStateA[rnd index of thread]
	double2 v_norm; //v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState); //more efficient to generate two doubles in the one function than run curand_normal_double twice according to CUDA docs
	v_norm.x = v_norm.x * sqrt(simConsts[3] * JOULE_PER_EV / mass) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[3] * JOULE_PER_EV / mass) + V_DIST_MEAN; //normal dist -> maxwellian
	
	*z_part = MIN_Z_SIM;
	*v_part = abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void magnetosphereGenerator(double* v_part, double* mu_part, double* z_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{
	double2 v_norm; //two normal dist values returned to v_norm.x and v_norm.y; v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState);
	v_norm.x = v_norm.x * sqrt(simConsts[4] * JOULE_PER_EV / mass) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = v_norm.y * sqrt(simConsts[4] * JOULE_PER_EV / mass) + V_DIST_MEAN; //normal dist -> maxwellian

	*z_part = MAX_Z_SIM;
	*v_part = -abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void ionosphereScattering(double* v_part, double* mu_part, double* z_part, double* simConsts, double mass, curandStateMRG32k3a* rndState)
{	
	//some array + 1
	
	ionosphereGenerator(v_part, mu_part, z_part, simConsts, mass, rndState);
	//scattering distribution of some sort here - right now, reflects half of the time, generates a new particle the other half
	/*if (curand_normal_double(rndState) > 0)
	{
		*v_part *= -1; //sufficient to reflect?
		*z_part = MIN_Z_SIM; //necessary?
	}
	else
		ionosphereGenerator(v_part, mu_part, z_part, elecTF, rndState);*/
}

__global__ void computeKernel(double** partData_d, double** origData_d, double** LUT, double* simConsts, curandStateMRG32k3a* crndStateA, double simtime, double mass, double charge, long numParts)
{
	unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
	//if (thdInd > 20000 * simtime / DT) //add 20000 particles per timestep - need something other than "magic" 20000
		//return;

	double* v_d; double* mu_d; double* z_d;
	double* v_orig; double* mu_orig; double* z_orig;
	v_d = partData_d[0]; mu_d = partData_d[1]; z_d = partData_d[2];
	v_orig = origData_d[0]; mu_orig = origData_d[1]; z_orig = origData_d[2];

	if (z_d[thdInd] < 0.001) //if z is zero (or pretty near zero to account for FP error), generate particles - every other starting at bottom/top of sim
	{//previous way to index curandStates: (blockIdx.x * 2) + (threadIdx.x % 2) - this leads to each block accessing two curand states - 128 threads call the same state simultaneously and end up with the same values
		if (thdInd < numParts / 2) //need perhaps a better way to determine distribution of ionosphere/magnetosphere particles
			ionosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		else
			magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], simConsts, mass, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		
		v_orig[thdInd] = v_d[thdInd];
		mu_orig[thdInd] = mu_d[thdInd];
		z_orig[thdInd] = z_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / ((thdInd < numParts / 2) ? (B_AT_MIN_Z) : (B_AT_MAX_Z));
	}
	else if (simtime == 0) //copies data to arrays that track the initial distribution - if data is loaded in, the above block won't be called
	{
		v_orig[thdInd] = v_d[thdInd];
		mu_orig[thdInd] = mu_d[thdInd];
		z_orig[thdInd] = z_d[thdInd];
		mu_d[thdInd] = 0.5 * mass * mu_d[thdInd] * mu_d[thdInd] / ((thdInd < numParts / 2) ? (B_AT_MIN_Z) : (B_AT_MAX_Z));
	}
	else if (z_d[thdInd] < simConsts[1] * 0.999) //out of sim to the bottom, particle has 50% chance of reflecting, 50% chance of new particle
		//ionosphereScattering(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	else if (z_d[thdInd] > simConsts[2] * 1.001) //out of sim to the top, particle is lost, new one generated in its place
		//magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	
	//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime]
	double args[7];
	args[0] = 0.0;
	args[1] = v_d[thdInd];
	args[2] = mu_d[thdInd];
	args[3] = charge;
	args[4] = mass;
	args[5] = z_d[thdInd];
	args[6] = simtime;

	v_d[thdInd] += foRungeKuttaCUDA(args, 7, LUT);
	z_d[thdInd] += v_d[thdInd] * DT;
}

void Simulation170925::initializeSimulation()
{	
	logFile_m.createTimeStruct("Start Sim Init"); //index 1
	logFile_m.writeTimeDiff(0, 1);

	if (particleTypes_m.size() == 0)
	{
		std::cout << "Error: No particles in sim.  You need to add particles before calling this function.  Returning.\n";
		return;
	}

	//Allocate room in vectors for GPU Memory Pointers
	gpuDblMemoryPointers_m.resize(2 * particleTypes_m.size() + 2); //part 0 curr data, part 1 curr data... part 0 orig data, part 1 orig data... simconsts, LUT
	gpuOtherMemoryPointers_m.resize(2 * particleTypes_m.size() + 2); //part 0 curr 2D, part 1 curr 2D... part 0 orig 2D, part 1 orig 2D... curand, LUT 2D
	satelliteData_m.reserve(100); //not resize...Don't know the exact size here so need to use push_back

	//Allocate memory on GPU for elec/ions variables
	for (int parts = 0; parts < 2 * particleTypes_m.size(); parts++) //[0] = e data, [1] = i data, [2] = e orig data, [3] = i orig data
	{
		Particle* partTmp{ particleTypes_m.at(parts % particleTypes_m.size()) };
		size_t memSize{ partTmp->getNumberOfParticles() * partTmp->getNumberOfAttributes() * sizeof(double) };
		
		CUDA_CALL(cudaMalloc((void **)&gpuDblMemoryPointers_m.at(parts), memSize));
		CUDA_CALL(cudaMemset(gpuDblMemoryPointers_m.at(parts), 0, memSize));
		CUDA_CALL(cudaMalloc((void **)&gpuOtherMemoryPointers_m.at(parts), partTmp->getNumberOfAttributes() * sizeof(double*))); //2D array
	}

	//Array of sim characteristics - dt, sim min, sim max, t ion, t mag, v mean
	CUDA_CALL(cudaMalloc((void **)&gpuDblMemoryPointers_m.at(2 * particleTypes_m.size()), 6 * sizeof(double)));

	//Code to prepare random number generator to produce pseudo-random numbers (for normal dist)
	long long seed = time(NULL);
	CUDA_CALL(cudaMalloc((void **)&gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size()), NUMRNGSTATES * sizeof(curandStateMRG32k3a))); //sizeof(curandStateMRG32k3a) is 72 bytes
	initCurand <<< NUMRNGSTATES / 256, 256 >>> (reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size())), seed);

	//Location of LUTarray
	CUDA_CALL(cudaMalloc((void **)&gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), LUTNUMOFENTRS * LUTNUMOFCOLS * sizeof(double)));
	CUDA_CALL(cudaMalloc((void **)&gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1), LUTNUMOFCOLS * sizeof(double*))); ///changed from double** to double* - if you have weird errors check here first

	//Create Satellites for observation - export to python
	createSatellite(particleTypes_m.at(0), MIN_Z_SIM * 0.999, true, reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(0)), true, "bottomElectrons");
	createSatellite(particleTypes_m.at(1), MIN_Z_SIM * 0.999, true, reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(1)), false, "bottomIons");
	createSatellite(particleTypes_m.at(0), MAX_Z_SIM * 1.001, false, reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(0)), true, "topElectrons");
	createSatellite(particleTypes_m.at(1), MAX_Z_SIM * 1.001, false, reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(1)), false, "topIons");

	//For derived classes to add code
	initializeFollowOn();

	initialized_m = true;
	logFile_m.createTimeStruct("End Sim Init"); //index 2
	logFile_m.writeTimeDiff(1, 2);
}

void Simulation170925::copyDataToGPU()
{//copies particle distribution and associated data to GPU in preparation of iterative calculations over the data
	logFile_m.writeLogFileEntry("copyDataToGPU", "Start copy to GPU");
	
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	//Copies initial data of particles to GPU, if loaded
	for (int parts = 0; parts < particleTypes_m.size(); parts++)
	{
		if (particleTypes_m.at(parts)->getInitDataLoaded())
		{
			Particle* tmpPart{ particleTypes_m.at(parts) };
			size_t memSize{ tmpPart->getNumberOfParticles() * sizeof(double) };
			LOOP_OVER_1D_ARRAY(tmpPart->getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(gpuDblMemoryPointers_m.at(parts) + tmpPart->getNumberOfParticles() * iii, tmpPart->getOrigData().at(iii).data(), memSize, cudaMemcpyHostToDevice));)
		}//look over one more time later
	}

	//Copies array of sim characteristics to GPU - dt, sim min, sim max, t ion, t mag, v mean
	double data[]{ dt_m, simMin_m, simMax_m, tIon_m, tMag_m, vmean_m };
	CUDA_CALL(cudaMemcpy(gpuDblMemoryPointers_m.at(2 * particleTypes_m.size()), data, 6 * sizeof(double), cudaMemcpyHostToDevice));
	
	for (int iii = 0; iii < 2 * particleTypes_m.size(); iii++)
		setup2DArray <<< 1, 1 >>> (gpuDblMemoryPointers_m.at(iii), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(iii)), particleTypes_m.at(iii % particleTypes_m.size())->getNumberOfAttributes(), particleTypes_m.at(iii % particleTypes_m.size())->getNumberOfParticles());
	
	//copies E field LUT to the GPU
	CUDA_CALL(cudaMemcpy(gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), elcFieldLUT_m[0], LUTNUMOFCOLS * LUTNUMOFENTRS * sizeof(double), cudaMemcpyHostToDevice));
	setup2DArray <<< 1, 1 >>> (gpuDblMemoryPointers_m.at(2 * particleTypes_m.size() + 1), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1)), LUTNUMOFCOLS, LUTNUMOFENTRS);

	//For derived classes to add code
	copyDataToGPUFollowOn();

	copied_m = true;
	
	logFile_m.writeLogFileEntry("copyDataToGPU", "End copy to GPU");
}

void Simulation170925::iterateSimulation(int numberOfIterations)
{//conducts iterative calculations of data previously copied to GPU - runs the data through the computeKernel
	logFile_m.createTimeStruct("Start Iterate " + std::to_string(numberOfIterations)); //index 3
	logFile_m.writeLogFileEntry("iterateSimulation", "Start Iteration of Sim:  " + std::to_string(numberOfIterations));
	
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.  You also need to copy data to the GPU with Simulation::copyDataToGPU.\n";
		return;
	}

	if (!copied_m)
	{
		std::cout << "You haven't copied any data to the GPU with Simulation::copyDataToGPU.  Do that first or the GPU has no numbers to work on.\n";
		return;
	}

	//For derived classes to add code
	iterateSimulationFollowOnPreLoop();

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{
		//computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(0)), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(2)),
			//reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size())), reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1)), simTime_m, true);
		//replace with for loop over particleTypes_m.size()//
		//computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(1)), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(3)),
			//reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(4)), reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m.at(5)), simTime_m, false);
		
		for (int parts = 0; parts < particleTypes_m.size(); parts++)
			computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(parts)), reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(parts + particleTypes_m.size())),
				reinterpret_cast<double**>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size() + 1)), gpuDblMemoryPointers_m.at(2 * particleTypes_m.size()),
				reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m.at(2 * particleTypes_m.size())), simTime_m,	particleTypes_m.at(parts)->getMass(),
				particleTypes_m.at(parts)->getCharge(), particleTypes_m.at(parts)->getNumberOfParticles());

		for (int sats = 0; sats < satellites_m.size(); sats++)
			satellites_m.at(sats)->iterateDetector(NUMBLOCKS, BLOCKSIZE, simTime_m);
		
		cudaloopind++;
		incTime();

		if (cudaloopind % LOOPS_BTW_PROGRESS_COUT == 0)
		{
			//std::string indstr{ std::to_string(cudaloopind) };
			//std::string timestr{ std::to_string(simTime_m) };
			//stringPadder(indstr, 5);
			//stringPadder(timestr, 3);
			std::cout << cudaloopind << " / " << numberOfIterations << "  |  Sim Time (s): " << simTime_m << "  |  Real Time Elapsed (s): ";
			logFile_m.printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << "\n";
		}

		//if (cudaloopind % 2500 == 0)//need better conditional
			//receiveSatelliteData();

		//For derived classes to add code
		iterateSimulationFollowOnInsideLoop();
	}
	receiveSatelliteData(); //not super generic - fix later
	std::cout << "Receive sat data outside main loop.  Remove after.\n";

	//For derived classes to add code
	iterateSimulationFollowOnPostLoop();

	logFile_m.createTimeStruct("End Iterate " + std::to_string(numberOfIterations)); //index 4
	logFile_m.writeTimeDiffFromNow(3, "End Iterate " + std::to_string(numberOfIterations));
	logFile_m.writeLogFileEntry("iterateSimulation", "End Iteration of Sim:  " + std::to_string(numberOfIterations));
}

void Simulation170925::copyDataToHost()
{//copies data back to host from GPU
	logFile_m.writeLogFileEntry("copyDataToHost", "Copy simulation data from GPU back to host");
	
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	for (int parts = 0; parts < particleTypes_m.size(); parts++)
	{
		Particle* tmpPart{ particleTypes_m.at(parts) };
		size_t memSize{ tmpPart->getNumberOfParticles() * sizeof(double) };
		long numParts{ tmpPart->getNumberOfParticles() };
		LOOP_OVER_1D_ARRAY(tmpPart->getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(tmpPart->getCurrData().at(iii).data(), gpuDblMemoryPointers_m.at(parts) + numParts * iii, memSize, cudaMemcpyDeviceToHost));)
		LOOP_OVER_1D_ARRAY(tmpPart->getNumberOfAttributes(), CUDA_CALL(cudaMemcpy(tmpPart->getOrigData().at(iii).data(), gpuDblMemoryPointers_m.at(parts + particleTypes_m.size()) + numParts * iii, memSize, cudaMemcpyDeviceToHost));)
	}

	//For derived classes to add code
	copyDataToHostFollowOn();

	logFile_m.writeLogFileEntry("copyDataToHost", "Done with copying.");
}

void Simulation170925::freeGPUMemory()
{//used to free the memory on the GPU that's no longer needed
	logFile_m.writeLogFileEntry("freeGPUMemory", "Start free GPU Memory.");

	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	LOOP_OVER_1D_ARRAY(gpuDblMemoryPointers_m.size(), CUDA_CALL(cudaFree(gpuDblMemoryPointers_m.at(iii)));)
	LOOP_OVER_1D_ARRAY(gpuOtherMemoryPointers_m.size(), CUDA_CALL(cudaFree(gpuOtherMemoryPointers_m.at(iii)));)

	//For derived classes to add code
	freeGPUMemoryFollowOn();

	logFile_m.writeLogFileEntry("freeGPUMemory", "End free GPU Memory.");

	CUDA_CALL(cudaProfilerStop()); //For profiling the profiler in the CUDA bundle}
}