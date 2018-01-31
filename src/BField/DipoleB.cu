#include "BField\DipoleB.h"

//B Field related kernels
__host__ __device__ double getSAtLambda_DipoleB(double* consts, int arrayLength, double lambdaDegrees)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{
	double x{ asinh(sqrt(3.0) * sinpi(lambdaDegrees / 180.0)) };

	return (0.5 * consts[2] / sqrt(3.0)) * (x + sinh(x) * cosh(x)); /* L */
}

__host__ __device__ double getLambdaAtS_DipoleB(double* consts, int arrayLength, double s)
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_tmp{ (-consts[1] / consts[4]) * s + consts[1] }; //-ILAT / s_max * s + ILAT
	double s_tmp{ consts[4] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp) };
	double dlambda{ 1.0 };
	bool   over{ 0 };

	while (abs((s_tmp - s) / s) > consts[6]) //errorTolerance
	{
		while (1)
		{
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = consts[4] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp);
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = consts[4] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp);
				if (s_tmp >= s)
					break;
			}
		}
		if (dlambda < consts[6] / 100.0) //errorTolerance
			break;
		dlambda /= 5.0; //through trial and error, this reduces the number of calculations usually (compared with 2, 2.5, 3, 4, 10)
	}

	return lambda_tmp;
}

__host__ __device__ double BFieldAtS_DipoleB(double* consts, int arrayLength, double s, double simtime)
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_deg{ getLambdaAtS_DipoleB(consts, arrayLength, s) };
	double rnorm{ consts[3] * cospi(lambda_deg / 180.0) * cospi(lambda_deg / 180.0) };

	return -consts[0] / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sinpi(lambda_deg / 180.0) * sinpi(lambda_deg / 180.0));
}

__host__ __device__ double gradBAtS_DipoleB(double* consts, int arrayLength, double s, double simtime)
{
	return (BFieldAtS_DipoleB(consts, arrayLength, s + consts[5], simtime) - BFieldAtS_DipoleB(consts, arrayLength, s - consts[5], simtime)) / (2 * consts[5]);
}

//setup CUDA kernels
__global__ void setupEnvironmentGPU_DipoleB(double* constArrayPtr)
{
	BFieldFcnPtr_GPU = BFieldAtS_DipoleB;
	gradBFcnPtr_GPU = gradBAtS_DipoleB;

	arraySize_GPU = 7;
	fieldConstArray_GPU = constArrayPtr;
}

//DipoleB class member functions
void DipoleB::setupEnvironment()
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	CUDA_CALL(cudaMalloc((void **)&fieldConstants_d, fieldConstArray_m.size() * sizeof(double)));
	CUDA_CALL(cudaMemcpy(fieldConstants_d, fieldConstArray_m.data(), fieldConstArray_m.size() * sizeof(double), cudaMemcpyHostToDevice));

	setupEnvironmentGPU_DipoleB <<< 1, 1 >>> (fieldConstants_d);
	cudaDeviceSynchronize();
}

double DipoleB::getBFieldAtS(double s, double t)
{
	return BFieldAtS_DipoleB(fieldConstArray_m.data(), fieldConstArray_m.size(), s, t);
}

double DipoleB::getGradBAtS(double s, double t)
{
	return gradBAtS_DipoleB(fieldConstArray_m.data(), fieldConstArray_m.size(), s, t);
}