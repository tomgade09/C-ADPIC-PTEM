//Standard Library includes
#include <string>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <time.h>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include "curand_kernel.h"

//Project specific includes
#include "include\_simulationvariables.h"
#include "include\Simulation170925.h"

//Array Size Variables
const int DBLARRAY_BYTES { NUMPARTICLES * sizeof(double) }; //Global vars - not my favorite solution, but I suppose it works for now
const int BOOLARRAY_BYTES{ NUMPARTICLES * sizeof(bool) };
const int INTARRAY_BYTES { NUMPARTICLES * sizeof(int) };

__host__ __device__ double EFieldatZ(double z, double simtime)
{//E Field in the direction of B (radially outward)
	if ((z > E_RNG_CENTER + E_RNG_DELTA) || (z < E_RNG_CENTER - E_RNG_DELTA))
		return 0.0;
	return CONSTEFIELD;
}

__host__ __device__ double BFieldatZ(double z, double simtime) //this will change in future iterations
{//for now, a simple dipole field
	return DIPOLECONST / pow(z / (RADIUS_EARTH / NORMFACTOR), 3);
}

__global__ void initCurand(curandStateMRG32k3a* state, long long seed)
{
	long long id = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(seed, id, 0, &state[id]);
}

__device__ double normalGeneratorCUDA(curandStateMRG32k3a* state, long long id, double mean, double sigma)
{
	curandStateMRG32k3a localState = state[id];
	
	double res = sigma * curand_normal_double(&localState) + mean;
	state[id] = localState;

	return res;
}

__device__ double accel1dCUDA(double* args, int len) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	double F_lor, F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * EFieldatZ(ztmp, args[6]) / NORMFACTOR; //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * (-3 / (ztmp * pow(ztmp / (RADIUS_EARTH / NORMFACTOR), 3))) * DIPOLECONST; //mu in [kg.m^2 / s^2.T] = [N.m / T]

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(double* funcArg, int arrayLen)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = accel1dCUDA(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = DT / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = DT;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)
	
	return (k1 + 2 * k2 + 2 * k3 + k4) * DT / 6; //returns units of y, not dy / dt
}

__global__ void computeKernel(double* v_d, double* mu_d, double* z_d, bool* inSimBool, int* numEscaped, bool elecTF, curandStateMRG32k3a* crndStateA, double simtime)
{
	int iii = blockIdx.x * blockDim.x + threadIdx.x;
	int nrmGenIdx = (blockIdx.x * 2) + (threadIdx.x % 2);
	double mass;
	double q;

	if (elecTF)//would have to be changed in the event of multiple ion species
	{
		mass = MASS_ELECTRON;
		q = -1.0;
	}
	else
	{
		mass = MASS_PROTON;
		q = 1.0;
	}

#ifdef CUDANORMAL_TEST
	v_d[iii] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ));
	mu_d[iii] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ)) * 1e-21;

	if (iii % 2 == 0)
	{
		z_d[iii] = IONSPH_MIN_Z + 0.1;
		v_d[iii] = abs(v_d[iii]);
	}
	else
	{
		z_d[iii] = MAGSPH_MAX_Z - 0.1;
		v_d[iii] = -abs(v_d[iii]);
	}
	inSimBool[iii] = true;
	return;
#endif

	inSimBool[iii] = ((z_d[iii] < MAGSPH_MAX_Z) && (z_d[iii] > IONSPH_MIN_Z)); //Makes sure particles are within bounds

	double args[7];

	if (REPLENISH_E_I)
	{
		if (!inSimBool[iii])
		{
			inSimBool[iii] = true;
			v_d[iii] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ) * VPARACONST);
			numEscaped[iii] += 1;
			if (z_d[iii] < IONSPH_MIN_Z)
			{
				z_d[iii] = IONSPH_MIN_Z + 0.01;
				v_d[iii] = abs(v_d[iii]);
			}
			else
			{
				z_d[iii] = MAGSPH_MAX_Z - 0.01;
				v_d[iii] = -abs(v_d[iii]);
			}
			mu_d[iii] = pow(normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ)), 2) * 0.5 * mass / BFieldatZ(z_d[iii], simtime);
		}
	}

	if (inSimBool[iii])
	{//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime]
		args[0] = 0.0;
		args[1] = v_d[iii];
		args[2] = mu_d[iii];
		args[3] = CHARGE_ELEM * q;
		args[4] = mass;
		args[5] = z_d[iii];
		args[6] = simtime;

		v_d[iii] += foRungeKuttaCUDA(args, 7);
		z_d[iii] += v_d[iii] * DT;
	}
}

void Simulation170925::initializeSimulation()
{
	//Generate z values and convert v_perp to mu here
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
		{
			if (jjj % 2 == 0) //Generate z, every other particle starts at top/bottom of sim respectively
			{
				particles_m[iii][2][jjj] = IONSPH_MIN_Z + 0.01;
				particles_m[iii][0][jjj] = abs(particles_m[iii][0][jjj]);
			}
			else
			{
				particles_m[iii][2][jjj] = MAGSPH_MAX_Z - 0.01;
				particles_m[iii][0][jjj] = -abs(particles_m[iii][0][jjj]);
			}
		}//end for jjj
	}//end for iii
	
	mass_m.reserve(2);
	mass_m[0] = MASS_ELECTRON;
	mass_m[1] = MASS_PROTON;

	convertVPerpToMu();
	
	//Allocate memory on GPU for elec/ions variables
	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m; iii++)
	{//[0] = v_e_para, [1] = mu_e_para, [2] = z_e, [3-5] = same attributes for ions
		cudaMalloc((void **)&gpuDblMemoryPointers_m[iii], DBLARRAY_BYTES);
		if (iii < numberOfParticleTypes_m)
		{
			cudaMalloc((void **)&gpuBoolMemoryPointers_m[iii], BOOLARRAY_BYTES); //for inSim bool per particle
			cudaMalloc((void **)&gpuIntMemoryPointers_m[iii], INTARRAY_BYTES); //for escaped particle count
			cudaMemset(gpuIntMemoryPointers_m[iii], 0, INTARRAY_BYTES); //setting escaped particle count to 0
		}
	}

	//Code to prepare random number generator to produce pseudo-random numbers (for normal dist)
	gpuOtherMemoryPointers_m.reserve(1);
	if (REPLENISH_E_I)
	{
		curandStateMRG32k3a* mrgStates_dev;
		long long seed = time(NULL);
		cudaMalloc((void **)&mrgStates_dev, 392 * 2 * sizeof(curandStateMRG32k3a));
		initCurand <<< 49, 16 >>> (mrgStates_dev, seed); //2 per block, 128 threads per random generator
		gpuOtherMemoryPointers_m[0] = mrgStates_dev;
	}
	else
		gpuOtherMemoryPointers_m[0] = nullptr;

#ifdef NO_NORMALIZE_M
	std::string unitstring{ " m" };
#else
	std::string unitstring{ " Re" };
#endif

	//Print things related to simulation characteristics
	std::cout << "Sim between:      " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << unitstring << "\n";
	std::cout << "E Field between:  " << (E_RNG_CENTER - E_RNG_DELTA) << " - " << (E_RNG_CENTER + E_RNG_DELTA) << unitstring << "\n";
	std::cout << "Const E:          " << CONSTEFIELD << " V/m\n\n";
	std::cout << "Particle Number:  " << numberOfParticlesPerType_m << "\n";
	std::cout << "Iteration Number: " << NUMITERATIONS << "\n";
	std::cout << "Replenish lost p: "; (REPLENISH_E_I) ? (std::cout << "True\n\n") : (std::cout << "False\n\n");

	initialized_m = true;
}

void Simulation170925::copyDataToGPU()
{
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
		{
			cudaMemcpy(gpuDblMemoryPointers_m[iii * numberOfAttributesTracked_m + jjj], particles_m[iii][jjj], DBLARRAY_BYTES, cudaMemcpyHostToDevice);
			//std::cout << "Dbls cudaMemcpyHostToDevice: " << iii * numberOfAttributesTracked_m + jjj << "  ";
		}
		//std::cout << "\n";
		cudaMemcpy(gpuBoolMemoryPointers_m[iii], particlesInSim_m[iii], BOOLARRAY_BYTES, cudaMemcpyHostToDevice);
	}
}

void Simulation170925::iterateSimulation(int numberOfIterations)
{
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	long cudaloopind{ 0 };
	//Loop code
	while (cudaloopind < numberOfIterations)
	{//__global__ void computeKernel(double* v_d, double* mu_d, double* z_d, bool* inSimBool, int* numEscaped, bool elecTF, curandStateMRG32k3a* crndStateA, double simtime)
		computeKernel <<< numberOfParticlesPerType_m / BLOCKSIZE, BLOCKSIZE >>>(gpuDblMemoryPointers_m[0], gpuDblMemoryPointers_m[1], gpuDblMemoryPointers_m[2], 
			gpuBoolMemoryPointers_m[0], gpuIntMemoryPointers_m[0], 1, reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m);
		computeKernel <<< numberOfParticlesPerType_m / BLOCKSIZE, BLOCKSIZE >>> (gpuDblMemoryPointers_m[3], gpuDblMemoryPointers_m[4], gpuDblMemoryPointers_m[5], 
			gpuBoolMemoryPointers_m[1], gpuIntMemoryPointers_m[1], 0, reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m);
		cudaloopind++;
		incTime();

		if (cudaloopind % 1000 == 0)
			std::cout << cudaloopind << " / " << numberOfIterations << "\n";
	}
}

void Simulation170925::copyDataToHost()
{
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
		{
			cudaMemcpy(particles_m[iii][jjj], gpuDblMemoryPointers_m[iii * numberOfAttributesTracked_m + jjj], DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
			//std::cout << "Dbls cudaMemcpyDeviceToHost: " << iii * numberOfAttributesTracked_m + jjj << "  ";
		}
		//std::cout << "\n";
		cudaMemcpy(particlesInSim_m[iii], gpuBoolMemoryPointers_m[iii], BOOLARRAY_BYTES, cudaMemcpyDeviceToHost);
		cudaMemcpy(particlesEscaped_m[iii], gpuIntMemoryPointers_m[iii], INTARRAY_BYTES, cudaMemcpyDeviceToHost);
	}

	//not generic but oh well
	for (int iii = 0; iii < numberOfParticlesPerType_m; iii++)
	{
		totalElecEscaped_m += particlesEscaped_m[0][iii];
		totalIonsEscaped_m += particlesEscaped_m[1][iii];
	}

	std::cout << "Electrons escaped: " << totalElecEscaped_m << "\n";
	std::cout << "Ions escaped:      " << totalIonsEscaped_m << "\n";
}

void Simulation170925::freeGPUMemory()
{	
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	//Destroy previously created rn generator
	if (REPLENISH_E_I)
		cudaFree(gpuOtherMemoryPointers_m[0]);

	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m; iii++)
	{
		cudaFree(gpuDblMemoryPointers_m[iii]);
		if (iii < numberOfParticleTypes_m)
		{
			cudaFree(gpuBoolMemoryPointers_m[iii]);
			cudaFree(gpuIntMemoryPointers_m[iii]);
		}
	}

	int e_in_sim{ 0 };
	int i_in_sim{ 0 };
	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{
		if (particlesInSim_m[0][iii])
			e_in_sim++;
		if (particlesInSim_m[1][iii])
			i_in_sim++;
	}

	std::cout << "C++: " << e_in_sim << " " << i_in_sim << " " << ((e_in_sim + i_in_sim) * 3) + 4 << "\n";
	//std::cout << "Sim Characteristics post compute - simTime, dt, # Types, # Particles/Type, # attr tracked: \n" << simTime_m << ", " << dt_m << ", " << numberOfParticleTypes_m;
	//std::cout << ", " << numberOfParticlesPerType_m << ", " << numberOfAttributesTracked_m << "\n";
	cudaProfilerStop(); //For profiling
}