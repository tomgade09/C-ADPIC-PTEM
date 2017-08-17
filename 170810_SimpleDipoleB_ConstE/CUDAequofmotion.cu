//Standard Library includes
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

__host__ __device__ double EFieldatZ(double z)
{
	if ((z > E_RNG_CENTER + E_RNG_DELTA) || (z < E_RNG_CENTER - E_RNG_DELTA))
		return 0.0;
	return CONSTEFIELD;
}

__host__ __device__ double BFieldatZ(double z) //this will change in future iterations
{//for now, a simple dipole field
	return DIPOLECONST / pow(z, 3);
}

__global__ void initKernel(curandStateMRG32k3a* state, long long seed)
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
{//args array: [dt, vz, mu, q, m, pz_0]
	double F_lor, F_mir;
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * EFieldatZ(args[5]) / NORMFACTOR; //will need to replace E with a function to calculate in more complex models

	//Mirror force
	F_mir = -args[2] * (-3 / pow(args[5], 4)) * DIPOLECONST; //have function for gradB based on dipole B field - will need to change later

	return (F_lor + F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

__device__ double foRungeKuttaCUDA(double* funcArg, int arrayLen, double h)
{
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = accel1dCUDA(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = h / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = h;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)
	
	return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns units of y, not dy / dt
}

__global__ void computeKernel(double* v_d, double* mu_d, double* z_d, bool* inSimBool, bool elecTF, curandStateMRG32k3a* crndStateA)
{
	int iii = blockIdx.x * blockDim.x + threadIdx.x;
	int nrmGenIdx = (blockIdx.x * 2) + (threadIdx.x % 2);
	double mass;
	double q;

	if (elecTF)
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

	double args[6];

	if (REPLENISH_E_I)
	{
		if (!inSimBool[iii])
		{
			inSimBool[iii] = true;
			v_d[iii] = normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ) * VPARACONST);
			if (z_d[iii] < IONSPH_MIN_Z)
			{
				z_d[iii] = IONSPH_MIN_Z + 0.1;
				v_d[iii] = abs(v_d[iii]);
			}
			else
			{
				z_d[iii] = MAGSPH_MAX_Z - 0.1;
				v_d[iii] = -abs(v_d[iii]);
			}
			mu_d[iii] = pow(normalGeneratorCUDA(crndStateA, nrmGenIdx, V_DIST_MEAN, sqrt(V_SIGMA_SQ)), 2) * 0.5 * mass / BFieldatZ(z_d[iii]);
		}
	}

	if (inSimBool[iii])
	{
		args[0] = 0.0;
		args[1] = v_d[iii];
		args[2] = mu_d[iii];
		args[3] = CHARGE_ELEM * q;
		args[4] = mass;
		args[5] = z_d[iii];

		v_d[iii] += foRungeKuttaCUDA(args, 6, DT);
		z_d[iii] += v_d[iii] * DT;
	}
}

void mainCUDA(double** electrons, double** ions, bool* elec_in_sim_host, bool* ions_in_sim_host)
{
	long cudaloopind{ 0 };
	double* v_e_para_host; double* mu_e_host; double* z_e_host; double* v_i_para_host; double* mu_i_host; double* z_i_host;
	double* v_e_para_dev; double* mu_e_dev; double* z_e_dev; double* v_i_para_dev; double* mu_i_dev; double* z_i_dev;
	bool* elec_in_sim_dev; bool* ions_in_sim_dev;
	v_e_para_host = electrons[0];
	mu_e_host = electrons[1];
	z_e_host = electrons[2];
	v_i_para_host = ions[0];
	mu_i_host = ions[1];
	z_i_host = ions[2];

	const int DBLARRAY_BYTES{ NUMPARTICLES * sizeof(double) };
	const int BOOLARRAY_BYTES{ NUMPARTICLES * sizeof(bool) };

	//allocate memory on GPU
	cudaMalloc((void **) &v_e_para_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &mu_e_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &z_e_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &v_i_para_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &mu_i_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &z_i_dev, DBLARRAY_BYTES);
	cudaMalloc((void **) &elec_in_sim_dev, BOOLARRAY_BYTES);
	cudaMalloc((void **) &ions_in_sim_dev, BOOLARRAY_BYTES);

	//copy memory to device
	cudaMemcpy(v_e_para_dev, v_e_para_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(mu_e_dev, mu_e_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(z_e_dev, z_e_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(v_i_para_dev, v_i_para_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(mu_i_dev, mu_i_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(z_i_dev, z_i_host, DBLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(elec_in_sim_dev, elec_in_sim_host, BOOLARRAY_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy(ions_in_sim_dev, ions_in_sim_host, BOOLARRAY_BYTES, cudaMemcpyHostToDevice);
	
	//Code to prepare random number generator to produce pseudo-random numbers (for normal dist)
	curandStateMRG32k3a* mrgStates_dev;
	
	if (REPLENISH_E_I)
	{
		long long seed = time(NULL);
		cudaMalloc((void **) &mrgStates_dev, 392 * 2 * sizeof(curandStateMRG32k3a));
		initKernel<<< 49, 16 >>>(mrgStates_dev, seed); //2 per block, 128 threads per random generator
	}
	
	//Loop code
	while (cudaloopind < NUMITERATIONS)
	{
		computeKernel<<< NUMPARTICLES / BLOCKSIZE, BLOCKSIZE >>>(v_e_para_dev, mu_e_dev, z_e_dev, elec_in_sim_dev, 1, mrgStates_dev);
		computeKernel<<< NUMPARTICLES / BLOCKSIZE, BLOCKSIZE >>>(v_i_para_dev, mu_i_dev, z_i_dev, ions_in_sim_dev, 0, mrgStates_dev);
		cudaloopind++;
		cudaDeviceSynchronize();
		if (cudaloopind % 1000 == 0)
			std::cout << cudaloopind << " / " << NUMITERATIONS << "\n";
	}

	//Destroy previously created rn generator
	if (REPLENISH_E_I)
		cudaFree(mrgStates_dev);

	//Copy data back out
	cudaMemcpy(v_e_para_host, v_e_para_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(mu_e_host, mu_e_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(z_e_host, z_e_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(v_i_para_host, v_i_para_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(mu_i_host, mu_i_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(z_i_host, z_i_dev, DBLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(elec_in_sim_host, elec_in_sim_dev, BOOLARRAY_BYTES, cudaMemcpyDeviceToHost);
	cudaMemcpy(ions_in_sim_host, ions_in_sim_dev, BOOLARRAY_BYTES, cudaMemcpyDeviceToHost);

	//Free memory
	cudaFree(v_e_para_dev);
	cudaFree(mu_e_dev);
	cudaFree(z_e_dev);
	cudaFree(v_i_para_dev);
	cudaFree(mu_i_dev);
	cudaFree(z_i_dev);
	cudaFree(elec_in_sim_dev);
	cudaFree(ions_in_sim_dev);

	cudaProfilerStop(); //For profiling
}