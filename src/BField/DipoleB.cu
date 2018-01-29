#include "BField\DipoleB.h"

//B Field related kernels
__host__ __device__ double getSAtLambda_DipoleB(double* consts, int arrayLength, double lambdaDegrees)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{//returns s in units of L
	//double xtmp{ sqrt(3.0) * sin(lambdaDegrees * PI / 180) };
	//double x{ log(xtmp + sqrt(xtmp * xtmp + 1)) };
	double x{ asinh(sqrt(3.0) * sin(lambdaDegrees * PI / 180.0)) };

	return (0.5 * /* L */consts[2] / sqrt(3.0)) * (x + sinh(x) * cosh(x));
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
	double lambda_rad{ lambda_deg * PI / 180.0 };
	double rnorm{ consts[3] * pow(cos(lambda_rad), 2) };

	return -consts[0] / pow(rnorm, 3) * sqrt(1.0 + 3 * pow(sin(lambda_rad), 2));
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
	getSAtLambdaPtr_GPU = getSAtLambda_DipoleB;

	arraySize_GPU = 7;
	fieldConstArray_GPU = constArrayPtr;

	printf("%d, %f, %f, %f, %f, %f, %f, %f\n", arraySize_GPU, fieldConstArray_GPU[0], fieldConstArray_GPU[1], fieldConstArray_GPU[2], fieldConstArray_GPU[3], fieldConstArray_GPU[4], fieldConstArray_GPU[5], fieldConstArray_GPU[6]);
	printf("%.6e, %.6e\n", BFieldFcnPtr_GPU(fieldConstArray_GPU, arraySize_GPU, 0.0, 0.0), gradBFcnPtr_GPU(fieldConstArray_GPU, arraySize_GPU, 6.3712e6, 0.0));
	printf("%f, %f\n", BFieldAtS_DipoleB(fieldConstArray_GPU, arraySize_GPU, 0.0, 0.0), gradBAtS_DipoleB(fieldConstArray_GPU, arraySize_GPU, 6.3712e6, 0.0));
	printf("%f\n", 85670894.1 - getSAtLambdaPtr_GPU(fieldConstArray_GPU, arraySize_GPU, 70.29323259));
}

//DipoleB class member functions
void DipoleB::setupEnvironment()
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	std::cout << "Size: " << fieldConstArray_m.size() << ": " << fieldConstArray_m.at(0) << ", " << fieldConstArray_m.at(1) << ", " << fieldConstArray_m.at(2) << ", " << fieldConstArray_m.at(3) << ", " << fieldConstArray_m.at(4) << ", " << fieldConstArray_m.at(5) << std::endl;
	CUDA_CALL(cudaMalloc((void **)&fieldConstants_d, fieldConstArray_m.size() * sizeof(double)));
	CUDA_CALL(cudaMemcpy(fieldConstants_d, fieldConstArray_m.data(), fieldConstArray_m.size() * sizeof(double), cudaMemcpyHostToDevice));

	std::cout << "================= GPU" << std::endl;
	setupEnvironmentGPU_DipoleB <<< 1, 1 >>> (fieldConstants_d);
	cudaDeviceSynchronize();
	std::cout << "================= HOST" << std::endl;
	std::cout << BFieldAtS_DipoleB(fieldConstArray_m.data(), 7, 0.0, 0.0) << ", " << gradBAtS_DipoleB(fieldConstArray_m.data(), 7, 6.3712e6, 0.0) << std::endl;
	std::cout << 85670894.1 - getSAtLambda_DipoleB(fieldConstArray_m.data(), 7, 70.29323259) << std::endl;
}

double DipoleB::getBFieldAtS(double s, double t)
{
	return BFieldAtS_DipoleB(fieldConstArray_m.data(), fieldConstArray_m.size(), s, t);
}

double DipoleB::getGradBAtS(double s, double t)
{
	return gradBAtS_DipoleB(fieldConstArray_m.data(), fieldConstArray_m.size(), s, t);
}