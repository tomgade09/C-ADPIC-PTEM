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
constexpr double EOMEGA{ 2 * PI / 10 };
constexpr double B_AT_MIN_Z{ B0ATTHETA / (MIN_Z_NORM * MIN_Z_NORM * MIN_Z_NORM) };
constexpr double B_AT_MAX_Z{ B0ATTHETA / (MAX_Z_NORM * MAX_Z_NORM * MAX_Z_NORM) };

//merge these two below...
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
	
	return linearInterpReal * cos(EOMEGA * simtime) + linearInterpImag * sin(EOMEGA * simtime);
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

	return linearInterpReal * cos(EOMEGA * simtime) + linearInterpImag * sin(EOMEGA * simtime);
}

__host__ __device__ double BFieldatZ(double z, double simtime)
{//for now, a simple dipole field
	return B0ATTHETA / pow(z / RADIUS_EARTH, 3); //Bz = B0 at theta * (1/rz(in Re))^3
}

__device__ double accel1dCUDA(double* args, int len, double* LUT) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	double F_lor, F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * EFieldatZ(LUT, ztmp, args[6] + args[0]); //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * B0ATTHETA * (-3 / (pow(ztmp / RADIUS_EARTH, 4))); //mu in [kg.m^2 / s^2.T] = [N.m / T]

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

__global__ void initCurand(curandStateMRG32k3a* state, long long seed)
{
	long long id = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(seed, id, 0, &state[id]);
}

__device__ void ionosphereGenerator(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{//takes pointers to single particle location in attribute arrays (ex: particle 100: ptr to v[100], ptr to mu[100], ptr to z[100], elecTF, ptr to crndStateA[rnd index of thread]
	curandStateMRG32k3a localState = *rndState;

	double2 v_norm; //v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(&localState); //more efficient to generate two doubles in the one function than run curand_normal_double twice according to CUDA docs
	v_norm.x = v_norm.x * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = v_norm.y * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) + V_DIST_MEAN; //normal dist -> maxwellian
	
	*z_part = MIN_Z_SIM;
	*v_part = abs(v_norm.x);
	*mu_part = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * v_norm.y * v_norm.y / B_AT_MIN_Z;

	*rndState = localState;
}

__device__ void magnetosphereGenerator(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{
	curandStateMRG32k3a localState = *rndState; //code I want to check

	double2 v_norm; //two normal dist values returned to v_norm.x and v_norm.y; v_norm.x = v_para; v_norm.y = v_perp
	v_norm = curand_normal2_double(&localState);
	v_norm.x = (v_norm.x * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) * sqrt(T_RATIO)) + V_DIST_MEAN; //normal dist -> maxwellian
	v_norm.y = (v_norm.x * (elecTF ? sqrt(V_SIGMA_SQ_ELEC) : sqrt(V_SIGMA_SQ_IONS)) * sqrt(T_RATIO)) + V_DIST_MEAN; //normal dist -> maxwellian

	*z_part = MAX_Z_SIM;
	*v_part = -abs(v_norm.x);
	*mu_part = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * v_norm.y * v_norm.y / B_AT_MAX_Z;

	*rndState = localState; //code I want to check
}

__device__ void ionosphereScattering(double* v_part, double* mu_part, double* z_part, bool elecTF, curandStateMRG32k3a* rndState)
{	
	//scattering distribution of some sort here - right now, reflects half of the time, generates a new particle the other half
	if (curand_normal_double(rndState) > 0)
	{
		*v_part *= -1; //sufficient to reflect?
		*z_part = MIN_Z_SIM; //necessary?
	}
	else
		ionosphereGenerator(v_part, mu_part, z_part, elecTF, rndState);
}

__global__ void computeKernel(double* v_d, double* mu_d, double* z_d, bool elecTF, curandStateMRG32k3a* crndStateA, double simtime, double* LUT)
{
	int thdInd = blockIdx.x * blockDim.x + threadIdx.x;
	if (thdInd > 20000 * simtime / DT) //add 20000 particles per timestep - need something other than "magic" 20000
		return;
	
	if (z_d[thdInd] < 0.001) //if z is zero (or pretty near zero to account for FP error), generate particles - every other starting at bottom/top of sim
	{
		if (thdInd % 2 == 0) //need perhaps a better way to determine distribution of ionosphere/magnetosphere particles
			ionosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
		else
			magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
	}
	else if (z_d[thdInd] < MIN_Z_SIM * 0.999) //out of sim to the bottom, particle has 50% chance of reflecting, 50% chance of new particle
		ionosphereScattering(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);
	else if (z_d[thdInd] > MAX_Z_SIM * 1.001) //out of sim to the top, particle is lost, new one generated in its place
		magnetosphereGenerator(&v_d[thdInd], &mu_d[thdInd], &z_d[thdInd], elecTF, &crndStateA[(blockIdx.x * 2) + (threadIdx.x % 2)]);

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
	gpuDblMemoryPointers_m.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m + 1);
	for (int iii = 0; iii < numberOfParticleTypes_m * numberOfAttributesTracked_m + 1; iii++) //[0] = v_e_para, [1] = mu_e_para, [2] = z_e, [3-5] = same attributes for ions, [6] = E Field LUT
		cudaMalloc((void **)&gpuDblMemoryPointers_m[iii], DBLARRAY_BYTES);

	//Code to prepare random number generator to produce pseudo-random numbers (for normal dist)
	gpuOtherMemoryPointers_m.reserve(1);
	curandStateMRG32k3a* mrgStates_dev;
	long long seed = time(NULL);
	cudaMalloc((void **)&mrgStates_dev, 392 * 2 * sizeof(curandStateMRG32k3a));
	initCurand <<< 49, 16 >>> (mrgStates_dev, seed); //2 per block, 128 threads per random generator
	gpuOtherMemoryPointers_m[0] = mrgStates_dev;

	double* elec[3];
	double* ions[3];
	
	for (int iii = 0; iii < 3; iii++)
	{
		elec[iii] = gpuDblMemoryPointers_m[iii];
		ions[iii] = gpuDblMemoryPointers_m[iii + 3];
	}

	//Create Satellites for observation
	createSatellite(2 * RADIUS_EARTH, true, elec, "downwardElectrons"); //Later code will take advantage of the interleaved particle order of the satellites
	createSatellite(2 * RADIUS_EARTH, true, ions, "downwardIons");	    //Look for [SOMEINDEX % 2]
	createSatellite(2 * RADIUS_EARTH, false, elec, "upwardElectrons");
	createSatellite(2 * RADIUS_EARTH, false, ions, "upwardIons");

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

	//copies E field LUT to the GPU
	double LUTtmp[3 * 2951];
	LOOP_OVER_2D_ARRAY(3, 2951, LUTtmp[iii * 2951 + jjj] = elcFieldLUT_m[iii][jjj];)

	cudaMemcpy(gpuDblMemoryPointers_m[6], LUTtmp, 3 * 2951 * sizeof(double), cudaMemcpyHostToDevice);
	
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
		computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (gpuDblMemoryPointers_m[0], gpuDblMemoryPointers_m[1], gpuDblMemoryPointers_m[2], 1,
			reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m, gpuDblMemoryPointers_m[6]);
		
		computeKernel <<< NUMBLOCKS, BLOCKSIZE >>> (gpuDblMemoryPointers_m[3], gpuDblMemoryPointers_m[4], gpuDblMemoryPointers_m[5], 0,
			reinterpret_cast<curandStateMRG32k3a*>(gpuOtherMemoryPointers_m[0]), simTime_m, gpuDblMemoryPointers_m[6]);
		
		for (int iii = 0; iii < satellites_m.size(); iii++)
			satellites_m[iii]->iterateDetector(NUMBLOCKS, BLOCKSIZE, simTime_m);
		
		cudaloopind++;
		incTime();

		if (cudaloopind % LOOPS_BTW_PROGRESS_COUT == 0)
		{
			std::cout << cudaloopind << " / " << numberOfIterations << "  |  Sim Time (s): " << simTime_m << "  |  Real Time Elapsed (s): ";
			logFile_m.printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << "\n";
		}

		if (cudaloopind % 2500 == 0)//need better conditional
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

				//convert mu to vperp
				for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
					if (tmp[2][jjj] != 0)
						tmp[1][jjj] = sqrt(2 * tmp[1][jjj] * BFieldatZ(tmp[2][jjj], simTime_m) / mass_m[iii % 2]);
				
				tmpcont.push_back(tmp);
			}
			satelliteData_m.push_back(tmpcont);
		}//end if (cudaloop...


	}//end while (cudaloop...
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

	LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, cudaMemcpy(particles_m[iii][jjj], gpuDblMemoryPointers_m[iii * numberOfAttributesTracked_m + jjj], DBLARRAY_BYTES, cudaMemcpyDeviceToHost);)

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