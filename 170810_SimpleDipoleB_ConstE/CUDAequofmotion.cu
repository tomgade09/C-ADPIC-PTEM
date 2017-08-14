#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

#include "include\_simulationvariables.h"

__device__ double accel1dCUDA(double* args, int len) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [dt, vz, mu, q, m, pz_0]
	double F_lor, F_mir;
	//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
	F_lor = args[3] * CONSTEFIELD / NORMFACTOR; //will need to replace E with a function to calculate in more complex models

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

__global__ void computeKernel(double* v_e_d, double* mu_e_d, double* z_e_d, double* v_i_d, double* mu_i_d, double* z_i_d, bool* elecBool, bool* ionsBool)
							//double** elecDbls, double** ionsDbls, bool* elecBool, bool* ionsBool) //each thread does one electron and one ion
{//dblArray: v_e_para, mu_e, z_e, v_i_para, mu_i, z_i
	int iii = blockIdx.x * blockDim.x + threadIdx.x;
	//make sure particle is still in bounds
	elecBool[iii] = ((z_e_d[iii] < MAGSPH_MAX_Z) && (z_e_d[iii] > IONSPH_MIN_Z)); //Makes sure particles are within bounds
	ionsBool[iii] = ((z_i_d[iii] < MAGSPH_MAX_Z) && (z_i_d[iii] > IONSPH_MIN_Z));

	double elecargs[6];
	double ionargs[6];

	if (elecBool[iii] == true)
	{
		elecargs[0] = 0.0;
		elecargs[1] = v_e_d[iii];
		elecargs[2] = mu_e_d[iii];
		elecargs[3] = CHARGE_ELEM * -1;
		elecargs[4] = MASS_ELECTRON;
		elecargs[5] = z_e_d[iii];

		v_e_d[iii] += foRungeKuttaCUDA(elecargs, 6, DT);
		z_e_d[iii] += v_e_d[iii] * DT;
	}
	
	if (ionsBool[iii] == true)
	{
		ionargs[0] = 0.0;
		ionargs[1] = v_i_d[iii];
		ionargs[2] = mu_i_d[iii];
		ionargs[3] = CHARGE_ELEM;
		ionargs[4] = MASS_PROTON;
		ionargs[5] = z_i_d[iii];

		v_i_d[iii] += foRungeKuttaCUDA(ionargs, 6, DT);
		z_i_d[iii] += v_i_d[iii] * DT;
	}
}

void mainCUDA(double** electrons, double** ions, bool* elec_in_sim_host, bool* ions_in_sim_host)
{
	unsigned long int cudaloopind{ 0 };
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

	while (cudaloopind < NUMITERATIONS)
	{									//computeKernel(v_e_d,			mu_e_d, z_e_d,		v_i_d,		mu_i_d,   z_i_d,	elecBool,		  ionsBool)
		computeKernel<<< NUMPARTICLES / 1000, 1000 >>> (v_e_para_dev, mu_e_dev, z_e_dev, v_i_para_dev, mu_i_dev, z_i_dev, elec_in_sim_dev, ions_in_sim_dev);
		cudaloopind++;
		if (cudaloopind % 1000 == 0)
			std::cout << cudaloopind << " / " << NUMITERATIONS << "\n";
	}

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

	//return values
	//double** eDataPtrArray = new double*[3];
	//double** iDataPtrArray = new double*[3];
	//double*** ret = new double**[2];
	return; //modify original array, so no need to pass anything back
}