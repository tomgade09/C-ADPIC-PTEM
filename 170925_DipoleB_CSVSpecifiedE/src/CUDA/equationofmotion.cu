//Standard Library includes
#include <string>
#include <iostream>
#include <stdio.h>
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
const int DBLARRAY_BYTES { NUMPARTICLES * sizeof(double) }; //Global vars - not my favorite solution, but I suppose it works for now
const int BOOLARRAY_BYTES{ NUMPARTICLES * sizeof(bool) };
const int INTARRAY_BYTES { NUMPARTICLES * sizeof(int) };

//Angular Freq of Efield
constexpr double OMEGA{ 2 * PI / 10 };

__host__ double EFieldatZ(double** LUT, double z, double simtime)
{//E Field in the direction of B (radially outward)
	//E-par = (column 2)*cos(omega*t) + (column 3)*sin(omega*t), omega = 2 PI / 10
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
	
	return linearInterpReal * cos(OMEGA * simtime) + linearInterpImag * sin(OMEGA * simtime);
}

__device__ double EFieldatZ(double* LUT, double z, double simtime)//biggest concern here
{//E Field in the direction of B (radially outward)
	if (z > RADIUS_EARTH) //in case z is passed in as m, not Re
		z = z / RADIUS_EARTH; //convert to Re
	
	if (z < LUT[0] || z > LUT[2950])
		return 0.0;
	double offset{ LUT[0] };
	int stepsFromZeroInd{ static_cast<int>(floor((z - offset) / (LUT[1] - LUT[0]))) }; //only works for constant bin size - if the binsize changes throughout LUT, need to iterate which will take longer
	
	//y = mx + b
	double linearInterpReal{ ((LUT[2951 + stepsFromZeroInd + 1] - LUT[2951 + stepsFromZeroInd]) / (LUT[stepsFromZeroInd + 1] - LUT[stepsFromZeroInd])) * 
		(z - LUT[stepsFromZeroInd]) + LUT[2951 + stepsFromZeroInd] };
	double linearInterpImag{ ((LUT[2 * 2951 + stepsFromZeroInd + 1] - LUT[2 * 2951 + stepsFromZeroInd]) / (LUT[stepsFromZeroInd + 1] - LUT[stepsFromZeroInd])) * 
		(z - LUT[stepsFromZeroInd]) + LUT[2 * 2951 + stepsFromZeroInd] };

	return linearInterpReal * cos(OMEGA * simtime) + linearInterpImag * sin(OMEGA * simtime);
}

__host__ __device__ double BFieldatZ(double z, double simtime) //this will change in future iterations
{//for now, a simple dipole field
	return DIPOLECONST * pow((RADIUS_EARTH / NORMFACTOR) / z, 3); //Bz = B0 * (r0/rz)^3, r0=DIPOLECONST - B at 1 Re
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

__device__ double accel1dCUDA(double* args, int len, double* LUT) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	double F_lor, F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * EFieldatZ(LUT, ztmp, args[6] + args[0]) / NORMFACTOR; //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * (-3 / (ztmp * pow(ztmp / (RADIUS_EARTH / NORMFACTOR), 3))) * DIPOLECONST; //mu in [kg.m^2 / s^2.T] = [N.m / T]

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(double* funcArg, int arrayLen, double* LUT)
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

__global__ void computeKernel(double* v_d, double* mu_d, double* z_d, bool* inSimBool, int* numEscaped, bool elecTF, curandStateMRG32k3a* crndStateA, double simtime, double* LUT)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	int nrmGenIdx = (blockIdx.x * 2) + (threadIdx.x % 2);//256 threads per block, 2 random generators per block, 128 threads per RG
	double mass;
	double q;
	double v_sigma;

	if (elecTF)//would have to be changed in the event of multiple ion species
	{
		mass = MASS_ELECTRON;
		q = -1.0;
		v_sigma = sqrt(V_SIGMA_SQ_ELEC);
	}
	else
	{
		mass = MASS_PROTON;
		q = 1.0;
		v_sigma = sqrt(V_SIGMA_SQ_IONS);
	}

#ifdef CUDANORMAL_TEST
	v_d[thdInd] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ));
	mu_d[thdInd] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ)) * 1e-21;

	if (thdInd % 2 == 0)
	{
		z_d[thdInd] = IONSPH_MIN_Z + 0.1;
		v_d[thdInd] = abs(v_d[thdInd]);
	}
	else
	{
		z_d[thdInd] = MAGSPH_MAX_Z - 0.1;
		v_d[thdInd] = -abs(v_d[thdInd]) * sqrt(T_RATIO);
	}
	inSimBool[thdInd] = true;
	return;
#endif

	inSimBool[thdInd] = ((z_d[thdInd] < MAX_Z_SIM) && (z_d[thdInd] > MIN_Z_SIM)); //Makes sure particles are within bounds

	double args[7];

	if (REPLENISH_E_I)
	{
		if (!inSimBool[thdInd])
		{
			inSimBool[thdInd] = true;
			v_d[thdInd] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, v_sigma * VPARACONST);
			numEscaped[thdInd] += 1;
			if (z_d[thdInd] < MIN_Z_SIM)
			{
				z_d[thdInd] = MIN_Z_SIM + 0.01 * (RADIUS_EARTH / NORMFACTOR);
				v_d[thdInd] = abs(v_d[thdInd]);
			}
			else
			{
				z_d[thdInd] = MAX_Z_SIM - 0.01 * (RADIUS_EARTH / NORMFACTOR);
				mu_d[thdInd] *= T_RATIO;
				v_d[thdInd] = -abs(v_d[thdInd]) * sqrt(T_RATIO);
			}
			mu_d[thdInd] = pow(normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, v_sigma), 2) * 0.5 * mass / BFieldatZ(z_d[thdInd], simtime);
		}
	}

	if (inSimBool[thdInd])
	{//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime]
		args[0] = 0.0;
		args[1] = v_d[thdInd];
		args[2] = mu_d[thdInd];
		args[3] = CHARGE_ELEM * q;
		args[4] = mass;
		args[5] = z_d[thdInd];
		args[6] = simtime;

		v_d[thdInd] += foRungeKuttaCUDA(args, 7, LUT);
		z_d[thdInd] += v_d[thdInd] * DT;
	}
}

void Simulation170925::initializeSimulation()
{
	//Allocate memory on GPU for elec/ions variables
	gpuDblMemoryPointers_m.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m + 1);
	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m + 1; iii++)
	{//[0] = v_e_para, [1] = mu_e_para, [2] = z_e, [3-5] = same attributes for ions, [6] = E Field LUT
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

	initialized_m = true;
}

void Simulation170925::copyDataToGPU()
{//copies particle distribution and associated data to GPU in preparation of iterative calculations over the data
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	//copies double arrays associated with particle distribution
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			cudaMemcpy(gpuDblMemoryPointers_m[iii * numberOfAttributesTracked_m + jjj], particles_m[iii][jjj], DBLARRAY_BYTES, cudaMemcpyHostToDevice);

		cudaMemcpy(gpuBoolMemoryPointers_m[iii], particlesInSim_m[iii], BOOLARRAY_BYTES, cudaMemcpyHostToDevice);
	}

	//copies E field LUT to the GPU
	double LUTtmp[3 * 2951];
	for (int iii = 0; iii < 3; iii++)
	{
		for (int jjj = 0; jjj < 2951; jjj++)
			LUTtmp[iii * 2951 + jjj] = elcFieldLUT_m[iii][jjj];
	}

	cudaMemcpy(gpuDblMemoryPointers_m[6], LUTtmp, 3 * 2951 * sizeof(double), cudaMemcpyHostToDevice);
	
	copied_m = true;
	
#ifdef LUTCOPY_TEST
	int serializedLUTerr{ 0 };
	int LUTtblerr{ 0 };
	double LUTfromGPU[3 * 2951];
	int tablerow{ -1 };
	cudaMemcpy(LUTfromGPU, gpuDblMemoryPointers_m[6], 3 * 2951 * sizeof(double), cudaMemcpyDeviceToHost);
	for (int iii = 0; iii < 3 * 2951; iii++)
	{
		if (iii % 2951 == 0)
			tablerow++;
		int tablecol{ iii % 2951 };
		if (LUTfromGPU[iii] != LUTtmp[iii])
			serializedLUTerr++;
		if (LUTfromGPU[iii] != elcFieldLUT_m[tablerow][tablecol])
			LUTtblerr++;
	}
	std::cout << "\n\nSerialized Errors: " << serializedLUTerr << "\nTable Errors: " << LUTtblerr << "\n\n";
#endif
}

void Simulation170925::iterateSimulation(int numberOfIterations)
{//conducts iterative calculations of data previously copied to GPU - runs the data through the computeKernel
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

	double** satelliteGPUDataPtrs[2];
	for (int iii = 0; iii < satellites_m.size(); iii++)
		satelliteGPUDataPtrs[iii] = new double*[3];
	satelliteGPUDataPtrs[0][0] = gpuDblMemoryPointers_m[0];//electrons
	satelliteGPUDataPtrs[0][1] = gpuDblMemoryPointers_m[1];
	satelliteGPUDataPtrs[0][2] = gpuDblMemoryPointers_m[2];
	satelliteGPUDataPtrs[1][0] = gpuDblMemoryPointers_m[3];//ions
	satelliteGPUDataPtrs[1][1] = gpuDblMemoryPointers_m[4];
	satelliteGPUDataPtrs[1][2] = gpuDblMemoryPointers_m[5];
	
	//Make room for 100 measurements
	satelliteData_m.reserve(100);

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{
		computeKernel <<< numberOfParticlesPerType_m / BLOCKSIZE, BLOCKSIZE >>> (gpuDblMemoryPointers_m[0], gpuDblMemoryPointers_m[1], gpuDblMemoryPointers_m[2], 
			gpuBoolMemoryPointers_m[0], gpuIntMemoryPointers_m[0], 1, reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m, gpuDblMemoryPointers_m[6]);
		computeKernel <<< numberOfParticlesPerType_m / BLOCKSIZE, BLOCKSIZE >>> (gpuDblMemoryPointers_m[3], gpuDblMemoryPointers_m[4], gpuDblMemoryPointers_m[5], 
			gpuBoolMemoryPointers_m[1], gpuIntMemoryPointers_m[1], 0, reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m, gpuDblMemoryPointers_m[6]);
		for (int iii = 0; iii < satellites_m.size(); iii++)
			satellites_m[iii]->iterateDetector(numberOfParticlesPerType_m / BLOCKSIZE, BLOCKSIZE, satelliteGPUDataPtrs[iii % 2]);
		
		cudaloopind++;
		incTime();

		if (cudaloopind % 1000 == 0)
			std::cout << cudaloopind << " / " << numberOfIterations << "  Sim Time: " << simTime_m << "\n";

		if (cudaloopind % 3000 == 0)
		{
			std::vector<std::vector<double*>> tmpcont; //vector of satellites[attributes]
			tmpcont.reserve(satellites_m.size());
			for (int iii = 0; iii < satellites_m.size(); iii++)
			{
				satellites_m[iii]->copyDataToHost();

				std::vector<double*> tmp; //vector of attributes[individual particles (through double*)]
				for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				{
					double* dbltmp = new double[NUMPARTICLES];
					double* satDat{ satellites_m[iii]->getDataArrayPointer(jjj) };
					std::copy(&satDat[0], &satDat[NUMPARTICLES - 1], &dbltmp[0]);
					tmp.push_back(dbltmp);
				}

				for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
					tmp[1][jjj] = sqrt(2 * tmp[1][jjj] * BFieldatZ(tmp[2][jjj], simTime_m) / mass_m[iii % 2]);
				tmpcont.push_back(tmp);
			}
			satelliteData_m.push_back(tmpcont);
		}//end if (cudaloop...
	}//end while (cudaloop...
}

void Simulation170925::copyDataToHost()
{//copies data back to host from GPU
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			cudaMemcpy(particles_m[iii][jjj], gpuDblMemoryPointers_m[iii * numberOfAttributesTracked_m + jjj], DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
		
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
{//used to free the memory on the GPU that's no longer needed
	if (!initialized_m)
	{
		std::cout << "You haven't initialized the simulation yet with Simulation::initializeSimulation.  Do that first.\n";
		return;
	}

	//Destroy previously created rn generator
	if (REPLENISH_E_I)
		cudaFree(gpuOtherMemoryPointers_m[0]);

	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m + 1; iii++)
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
	cudaProfilerStop(); //For profiling the profiler in the CUDA bundle
}