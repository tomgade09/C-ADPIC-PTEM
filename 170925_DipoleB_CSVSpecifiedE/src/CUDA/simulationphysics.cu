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
constexpr int DBLARRAY_BYTES { NUMPARTICLES * sizeof(double) };
//constexpr int BOOLARRAY_BYTES{ NUMPARTICLES * sizeof(bool) };
//constexpr int INTARRAY_BYTES { NUMPARTICLES * sizeof(int) };

//Commonly used values
constexpr double EOMEGA{ 20 * PI }; //10 Hz wave
constexpr double B_AT_MIN_Z{ B0ATTHETA / (MIN_Z_NORM * MIN_Z_NORM * MIN_Z_NORM) };
constexpr double B_AT_MAX_Z{ B0ATTHETA / (MAX_Z_NORM * MAX_Z_NORM * MAX_Z_NORM) };

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

__device__ void ionosphereGenerator(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{//takes pointers to single particle location in attribute arrays (ex: particle 100: ptr to v[100], ptr to mu[100], ptr to z[100], elecTF, ptr to crndStateA[rnd index of thread]
	double2 v_norm; //v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState); //more efficient to generate two doubles in the one function than run curand_normal_double twice according to CUDA docs
	v_norm.x = v_norm.x * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = v_norm.y * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) + V_DIST_MEAN; //normal dist -> maxwellian
	
	*z_part = MIN_Z_SIM;
	*v_part = abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void magnetosphereGenerator(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{
	double2 v_norm; //two normal dist values returned to v_norm.x and v_norm.y; v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(rndState);
	v_norm.x = (v_norm.x * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) * sqrt(T_RATIO)) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = (v_norm.y * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) * sqrt(T_RATIO)) + V_DIST_MEAN; //normal dist -> maxwellian

	*z_part = MAX_Z_SIM;
	*v_part = -abs(v_norm.x);
	*mu_part = v_norm.y;
}

__device__ void ionosphereScattering(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{	
	//some array + 1
	
	ionosphereGenerator(v_part, mu_part, z_part, elecTF, rndState);
	//scattering distribution of some sort here - right now, reflects half of the time, generates a new particle the other half
	/*if (curand_normal_double(rndState) > 0)
	{
		*v_part *= -1; //sufficient to reflect?
		*z_part = MIN_Z_SIM; //necessary?
	}
	else
		ionosphereGenerator(v_part, mu_part, z_part, elecTF, rndState);*/
}

__global__ void computeKernel(double** partData_d, double** origData_d, double** LUT, curandStateMRG32k3a* crndStateA, double simtime, bool elecTF)
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
		if (thdInd % 2 == 0) //need perhaps a better way to determine distribution of ionosphere/magnetosphere particles
			ionosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		else
			magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x % (NUMRNGSTATES / BLOCKSIZE)) * blockDim.x + (threadIdx.x)]);
		
		v_orig[thdInd] = v_d[thdInd];
		mu_orig[thdInd] = mu_d[thdInd];
		z_orig[thdInd] = z_d[thdInd];
		mu_d[thdInd] = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * mu_d[thdInd] * mu_d[thdInd] / ((thdInd % 2 == 0) ? (B_AT_MIN_Z) : (B_AT_MAX_Z));
	}
	else if (simtime == 0) //copies data to arrays that track the initial distribution - if data is loaded in, the above block won't be called
	{
		v_orig[thdInd] = v_d[thdInd];
		mu_orig[thdInd] = mu_d[thdInd];
		z_orig[thdInd] = z_d[thdInd];
		mu_d[thdInd] = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * mu_d[thdInd] * mu_d[thdInd] / ((thdInd % 2 == 0) ? (B_AT_MIN_Z) : (B_AT_MAX_Z));
	}
	else if (z_d[thdInd] < MIN_Z_SIM * 0.999) //out of sim to the bottom, particle has 50% chance of reflecting, 50% chance of new particle
		//ionosphereScattering(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	else if (z_d[thdInd] > MAX_Z_SIM * 1.001) //out of sim to the top, particle is lost, new one generated in its place
		//magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		return;
	
	//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime]
	double args[7];
	args[0] = 0.0;
	args[1] = v_d[thdInd];
	args[2] = mu_d[thdInd];
	args[3] = CHARGE_ELEM * ((elecTF) ? (-1.0) : (1.0));
	args[4] = (elecTF) ? (MASS_ELECTRON) : (MASS_PROTON);
	args[5] = z_d[thdInd];
	args[6] = simtime;

	v_d[thdInd] += foRungeKuttaCUDA(args, 7, LUT);
	z_d[thdInd] += v_d[thdInd] * DT;
}

void Simulation170925::initializeSimulation()
{	
	logFile_m.createTimeStruct("Start Sim Init"); //index 1
	logFile_m.writeTimeDiff(0, 1);
	
	//Allocate memory on GPU for elec/ions variables
	gpuDblMemoryPointers_m.reserve(2 * numberOfParticleTypes_m * numberOfAttributesTracked_m + 1);
	for (int iii = 0; iii < 2 * numberOfParticleTypes_m; iii++) //[0] = e data, [1] = i data, [2] = e orig data, [3] = i orig data
	{
		cudaMalloc((void **)&gpuDblMemoryPointers_m[iii], numberOfAttributesTracked_m * DBLARRAY_BYTES);
		cudaMemset((void **)&gpuDblMemoryPointers_m[iii], 0, numberOfAttributesTracked_m * DBLARRAY_BYTES);
		cudaMalloc((void **)&gpuOtherMemoryPointers_m[iii], numberOfAttributesTracked_m * sizeof(double*)); //2D array
	}

	//Location of LUTarray
	cudaMalloc((void **)&gpuDblMemoryPointers_m[4], LUTNUMOFENTRS * LUTNUMOFCOLS * sizeof(double));
	cudaMalloc((void **)&gpuOtherMemoryPointers_m[4], LUTNUMOFCOLS * sizeof(double**));

	//Code to prepare random number generator to produce pseudo-random numbers (for normal dist)
	gpuOtherMemoryPointers_m.reserve(1);
	curandStateMRG32k3a* mrgStates_dev;
	long long seed = time(NULL);
	cudaMalloc((void **)&mrgStates_dev, NUMRNGSTATES * sizeof(curandStateMRG32k3a)); //sizeof(curandStateMRG32k3a) is 72 bytes
	initCurand <<< NUMRNGSTATES / 256, 256 >>> (mrgStates_dev, seed);
	gpuOtherMemoryPointers_m[5] = mrgStates_dev;

	//Create Satellites for observation
	//createSatellite(2 * RADIUS_EARTH, true, elec, true, "downwardElectrons");
	//createSatellite(2 * RADIUS_EARTH, true, ions, false, "downwardIons");
	//createSatellite(2 * RADIUS_EARTH, false, elec, true, "upwardElectrons");
	//createSatellite(2 * RADIUS_EARTH, false, ions, false, "upwardIons");
	createSatellite(MIN_Z_SIM * 0.99999, true, reinterpret_cast<double**>(gpuOtherMemoryPointers_m[0]), true, "bottomElectrons");
	createSatellite(MIN_Z_SIM * 0.99999, true, reinterpret_cast<double**>(gpuOtherMemoryPointers_m[1]), false, "bottomIons");
	createSatellite(MAX_Z_SIM * 1.00001, false, reinterpret_cast<double**>(gpuOtherMemoryPointers_m[0]), true, "topElectrons");
	createSatellite(MAX_Z_SIM * 1.00001, false, reinterpret_cast<double**>(gpuOtherMemoryPointers_m[1]), false, "topIons");

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

	if (particles_m[0][2][0] != 0) //replace this with a flag specified in Simulation - although this works, not completely error proof
	{
		LOOP_OVER_1D_ARRAY(numberOfParticleTypes_m, cudaMemcpy(gpuDblMemoryPointers_m[iii], particles_m[iii][0], numberOfAttributesTracked_m * DBLARRAY_BYTES, cudaMemcpyHostToDevice);)
	}

	//copies E field LUT to the GPU
	cudaMemcpy(gpuDblMemoryPointers_m[4], elcFieldLUT_m[0], LUTNUMOFCOLS * LUTNUMOFENTRS * sizeof(double), cudaMemcpyHostToDevice);
	
	for (int iii = 0; iii < 2 * numberOfParticleTypes_m; iii++)
		setup2DArray <<< 1, 1 >>> (gpuDblMemoryPointers_m[iii], reinterpret_cast<double**>(gpuOtherMemoryPointers_m[iii]), numberOfAttributesTracked_m, numberOfParticlesPerType_m);
	setup2DArray <<< 1, 1 >>> (gpuDblMemoryPointers_m[4], reinterpret_cast<double**>(gpuOtherMemoryPointers_m[4]), LUTNUMOFCOLS, LUTNUMOFENTRS);

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

	//Make room for 100 measurements
	satelliteData_m.reserve(100);

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{
		computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (reinterpret_cast<double**>(gpuOtherMemoryPointers_m[0]), reinterpret_cast<double**>(gpuOtherMemoryPointers_m[2]),
			reinterpret_cast<double**>(gpuOtherMemoryPointers_m[4]), reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[5]), simTime_m, true);
		
		computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (reinterpret_cast<double**>(gpuOtherMemoryPointers_m[1]), reinterpret_cast<double**>(gpuOtherMemoryPointers_m[3]),
			reinterpret_cast<double**>(gpuOtherMemoryPointers_m[4]), reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[5]), simTime_m, false);
		
		for (int iii = 0; iii < satellites_m.size(); iii++)
			satellites_m[iii]->iterateDetector(NUMBLOCKS, BLOCKSIZE, simTime_m);
		
		//cudaDeviceSynchronize();
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
	}
	receiveSatelliteData();
	std::cout << "Receive sat data outside main loop.  Remove after.\n";
	
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
	
	mu_m = true;

	LOOP_OVER_1D_ARRAY(numberOfParticleTypes_m, cudaMemcpy(particles_m[iii][0], gpuDblMemoryPointers_m[iii], DBLARRAY_BYTES * numberOfAttributesTracked_m, cudaMemcpyDeviceToHost);)
	LOOP_OVER_1D_ARRAY(numberOfParticleTypes_m, cudaMemcpy(particlesorig_m[iii][0], gpuDblMemoryPointers_m[iii + 2], DBLARRAY_BYTES * numberOfAttributesTracked_m, cudaMemcpyDeviceToHost);)
	//LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, cudaMemcpy(particles_m[iii][jjj], &gpuDblMemoryPointers_m[iii][jjj * numberOfParticlesPerType_m], DBLARRAY_BYTES, cudaMemcpyDeviceToHost);)
	//LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, cudaMemcpy(particlesorig_m[iii][jjj], &gpuDblMemoryPointers_m[iii+2][jjj * numberOfParticlesPerType_m], DBLARRAY_BYTES, cudaMemcpyDeviceToHost);)
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

	//Destroy previously created rn generator
	cudaFree(gpuOtherMemoryPointers_m[0]);

	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m + 1; iii++)
		cudaFree(gpuDblMemoryPointers_m[iii]);

	logFile_m.writeLogFileEntry("freeGPUMemory", "End free GPU Memory.");

	cudaProfilerStop(); //For profiling the profiler in the CUDA bundle
}