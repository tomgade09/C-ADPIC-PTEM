#include "BField\DipoleB.h"

constexpr int GPUARRAYSIZE{ 6 };
__constant__ double constarr_DipoleB[GPUARRAYSIZE];

//B Field related kernels
__host__ __device__ double getSAtLambda_DipoleB(const double* consts, const int arrayLength, const double lambdaDegrees)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{
	//double x{ asinh(sqrt(3.0) * sinpi(lambdaDegrees / 180.0)) };
	double sinh_x{ sqrt(3.0) * sinpi(lambdaDegrees / 180.0) };
	double x{ log(sinh_x + sqrt(sinh_x * sinh_x + 1)) }; //trig identity for asinh - a bit faster - asinh(x) == ln(x + sqrt(x*x + 1))

	return (0.5 * consts[1] / sqrt(3.0)) * (x + 0.25 * (exp(2.0*x)-exp(-2.0*x))); /* L */ //0.25 * (exp(2*x)-exp(-2*x)) == sinh(x) * cosh(x) and is faster
}

__host__ __device__ double getLambdaAtS_DipoleB(const double* consts, const int arrayLength, const double s)
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_tmp{ (-consts[0] / consts[3]) * s + consts[0] }; //-ILAT / s_max * s + ILAT
	double s_tmp{ consts[3] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp) };
	double dlambda{ 1.0 };
	bool   over{ 0 };

	while (abs((s_tmp - s) / s) > consts[5]) //errorTolerance
	{
		while (1)
		{
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = consts[3] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp);
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = consts[3] - getSAtLambda_DipoleB(consts, arrayLength, lambda_tmp);
				if (s_tmp >= s)
					break;
			}
		}
		if (dlambda < consts[5] / 100.0) //errorTolerance
			break;
		dlambda /= 5.0; //through trial and error, this reduces the number of calculations usually (compared with 2, 2.5, 3, 4, 10)
	}

	return lambda_tmp;
}

__device__ double BFieldAtS_DipoleB(const double s, const double simtime)
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_deg{ getLambdaAtS_DipoleB(constarr_DipoleB, 6, s) };
	double rnorm{ constarr_DipoleB[2] * cospi(lambda_deg / 180.0) * cospi(lambda_deg / 180.0) };

	return -B0 / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sinpi(lambda_deg / 180.0) * sinpi(lambda_deg / 180.0));
}

__device__ double gradBAtS_DipoleB(const double s, const double simtime)
{
	return (BFieldAtS_DipoleB(s + constarr_DipoleB[4], simtime) - BFieldAtS_DipoleB(s - constarr_DipoleB[4], simtime)) / (2 * constarr_DipoleB[4]);
}

//setup CUDA kernels
__global__ void setupEnvironmentGPU_DipoleB()
{
	BFieldFcnPtr_GPU = BFieldAtS_DipoleB;
	gradBFcnPtr_GPU = gradBAtS_DipoleB;
}

//DipoleB class member functions
void DipoleB::setupEnvironment()
{// consts: [ ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	CUDA_API_ERRCHK(cudaMemcpyToSymbol(constarr_DipoleB, fieldConstArray_m.data(), GPUARRAYSIZE * sizeof(double)));
	setupEnvironmentGPU_DipoleB <<< 1, 1 >>> ();
	CUDA_KERNEL_ERRCHK_WSYNC();
}

double DipoleB::getBFieldAtS(double s, double t)
{
	double lambda_deg{ getLambdaAtS_DipoleB(fieldConstArray_m.data(), 6, s) };
	double rnorm{ L_norm_m * cospi(lambda_deg / 180.0) * cospi(lambda_deg / 180.0) };

	return -B0 / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sinpi(lambda_deg / 180.0) * sinpi(lambda_deg / 180.0));
}

double DipoleB::getGradBAtS(double s, double t)
{
	return (getBFieldAtS(s + ds_m, t) - getBFieldAtS(s - ds_m, t)) / (2 * ds_m);
}