#include "BField/DipoleB.h"

#include "device_launch_parameters.h"
#include "utils/serializationHandlers.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/cudaDeviceMacros.h"

using std::string;
using namespace utils::fileIO::serialize;

constexpr double B0{ 3.12e-5 }; //B_0 for Earth dipole B model

//setup CUDA kernels
__global__ void setupEnvironmentGPU_DipoleB(BField** this_d, double ILATDeg, double errTol, double ds)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_DipoleB", (*this_d) = new DipoleB(ILATDeg, errTol, ds));
}

__global__ void deleteEnvironmentGPU_DipoleB(BField** dipoleb)
{
	ZEROTH_THREAD_ONLY("deleteEnvironmentGPU_DipoleB", delete ((DipoleB*)(*dipoleb)));
}

__host__ __device__ DipoleB::DipoleB(degrees ILAT, double lambdaErrorTolerance, meters ds, bool useGPU) :
	BField("DipoleB"), ILAT_m{ ILAT }, ds_m{ ds }, lambdaErrorTolerance_m{ lambdaErrorTolerance }, useGPU_m{ useGPU }
{
	L_m = RADIUS_EARTH / pow(cos(ILAT_m * RADS_PER_DEG), 2);
	L_norm_m = L_m / RADIUS_EARTH;
	s_max_m = getSAtLambda(ILAT_m);

	#ifndef __CUDA_ARCH__ //host code
	if (useGPU_m) setupEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ DipoleB::DipoleB(string serialFolder)
{
	deserialize(serialFolder);
	if (useGPU_m) setupEnvironment();
}

__host__ DipoleB::~DipoleB()
{
	if (useGPU_m) deleteEnvironment();
}

__host__ degrees DipoleB::ILAT() const override
{
	return ILAT_m;
}


//B Field related kernels
__host__ __device__ meters DipoleB::getSAtLambda(const degrees lambda) const
{
	//double x{ asinh(sqrt(3.0) * sinpi(lambdaDegrees / 180.0)) }; //asinh triggers an odd cuda 8.x bug that is resolved in 9.x+
	double sinh_x{ sqrt(3.0) * sinpi(lambda / 180.0) };
	double x{ log(sinh_x + sqrt(sinh_x * sinh_x + 1)) }; //trig identity for asinh - a bit faster - asinh(x) == ln(x + sqrt(x*x + 1))

	return (0.5 * L_m / sqrt(3.0)) * (x + 0.25 * (exp(2.0*x)-exp(-2.0*x))); /* L */ //0.25 * (exp(2*x)-exp(-2*x)) == sinh(x) * cosh(x) and is faster
}

__host__ __device__ degrees DipoleB::getLambdaAtS(const meters s) const
{
	degrees lambda_tmp{ (-ILATDegrees_m / s_max_m) * s + ILATDegrees_m }; //-ILAT / s_max * s + ILAT
	meters  s_tmp{ s_max_m - getSAtLambda(lambda_tmp) };
	degrees dlambda{ 1.0 };
	bool   over{ 0 };

	while (abs((s_tmp - s) / s) > lambdaErrorTolerance_m)
	{
		while (1)
		{
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = s_max_m - getSAtLambda(lambda_tmp);
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = s_max_m - getSAtLambda(lambda_tmp);
				if (s_tmp >= s)
					break;
			}
		}
		if (dlambda < lambdaErrorTolerance_m / 100.0)
			break;
		dlambda /= 5.0; //through trial and error, this reduces the number of calculations usually (compared with 2, 2.5, 3, 4, 10)
	}

	return lambda_tmp;
}

__host__ __device__ tesla DipoleB::getBFieldAtS(const meters s, const seconds simtime) const
{// consts: [ ILATDeg, L, L_norm, s_max, ds, lambdaErrorTolerance ]
	double lambda_deg{ getLambdaAtS(s) };
	double rnorm{ L_norm_m * cospi(lambda_deg / 180.0) * cospi(lambda_deg / 180.0) };

	return -B0 / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sinpi(lambda_deg / 180.0) * sinpi(lambda_deg / 180.0));
}

__host__ __device__ double DipoleB::getGradBAtS(const meters s, const seconds simtime) const
{
	return (getBFieldAtS(s + ds_m, simtime) - getBFieldAtS(s - ds_m, simtime)) / (2 * ds_m);
}

__host__ __device__ meters DipoleB::getSAtAlt(const meters alt_fromRe) const
{
	double lambda{ acos(sqrt((alt_fromRe + RADIUS_EARTH) / L_m)) / RADS_PER_DEG };
	return s_max_m - getSAtLambda(lambda);
}

__host__ double DipoleB::getErrTol() const
{
	return lambdaErrorTolerance_m;
}

__host__ meters DipoleB::getds() const
{
	return ds_m;
}

__host__ void DipoleB::serialize(string serialFolder) const override
{
	string filename{ serialFolder + string("BField_DipoleB.ser") };

	if (std::filesystem::exists(filename))
		cerr << "DipoleB::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.";
	
	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("DipoleB::serialize: unable to create file: " + filename);
	
	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	out.write(reinterpret_cast<const char*>(this), sizeof(DipoleB));
	writeStrBuf(serializeString(string(name_m)));

	out.close();
}

__host__ void DipoleB::deserialize(string serialFolder) override
{
	string filename{ serialFolder + string("/BField_DipoleB.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("DipoleB::deserialize: unable to open file: " + filename);

	DipoleB* dipb{ nullptr };
	vector<char> dipbchar(sizeof(DipoleB));

	in.read(dipbchar.data(), sizeof(DipoleB));
	dipb = reinterpret_cast<DipoleB*>(dipbchar.data());
	
	this_d = nullptr;
	name_m = deserializeStr(in);
	
	L_m = dipb->L_m;
	L_norm_m = dipb->L_norm_m;
	s_max_m = dipb->s_max_m;
	ILAT_m = dipb->ILAT_m;
	ds_m = dipb->ds_m;
	lambdaErrorTolerance_m = dipb->lambdaErrorTolerance;

	useGPU_m = dipb->useGPU_m;
}

//DipoleB protected member functions
void DipoleB::setupEnvironment()
{// consts: [ ILATDeg, L, L_norm, s_max, ds, lambdaErrorTolerance ]
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(BField**)));
	setupEnvironmentGPU_DipoleB <<< 1, 1 >>> (this_d, ILAT_m, lambdaErrorTolerance_m, ds_m);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void DipoleB::deleteEnvironment()
{
	deleteEnvironmentGPU_DipoleB <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
}