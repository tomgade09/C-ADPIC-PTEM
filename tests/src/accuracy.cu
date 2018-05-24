#include "test.h"
#include "testHelperUtils.h"
#include "utils/fileIO.h"
#include "BField/DipoleB.h"
#include "BField/DipoleBLUT.h"
#include "EField/QSPS.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "device_launch_parameters.h"

#include <vector>
#include <memory>

using utils::fileIO::readDblBin;
using utils::fileIO::writeDblBin;

#define TEST_RESULTS_WITH_PASS(name, p) { if (p) { TEST_RESULTS(name, true); } else { pass = false; TEST_RESULTS(name, false) } }

namespace test
{
	__global__ void runB(BField** b, double smin, double ds, double* ret)
	{
		unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };
		ret[thdInd] = (*b)->getBFieldAtS(ds * thdInd + smin, 0.0);
	}

	//outside accessible tests
	bool dipoleAccuracy(bool save)
	{
		//Setup test
		constexpr int    NUMITERS{ 1024 };
		constexpr double SMIN{ 100000 };
		constexpr double SMAX{ 20000000 };
		constexpr double DS{ (SMAX - SMIN) / NUMITERS };

		std::vector<double> dipChkVals(NUMITERS);
		std::vector<double> lutChkVals(NUMITERS);
		TEST_EXCEP_CHECK(readDblBin(dipChkVals, "./../tests/data/dipole.bin", NUMITERS));
		TEST_EXCEP_CHECK(readDblBin(lutChkVals, "./../tests/data/dipolelut.bin", NUMITERS));

		std::vector<double> dipResults(NUMITERS);
		std::vector<double> lutResults(NUMITERS);
		std::vector<double> dipGPUResults(NUMITERS);
		std::vector<double> lutGPUResults(NUMITERS);

		std::unique_ptr<DipoleB> dip{ std::make_unique<DipoleB>(72.0) };
		std::unique_ptr<DipoleBLUT> lut{ std::make_unique<DipoleBLUT>(72.0, SMIN, SMAX, RADIUS_EARTH / 100.0, 1024) };

		auto runCUDAB = [](BField** b, double smin, double ds, int numiters, std::vector<double>& resOut)
		{
			if ((int)resOut.size() != numiters) resOut.resize(numiters);

			size_t cudamemsize{ numiters * sizeof(double) };
			double* tmp;

			CUDA_API_ERRCHK(cudaMalloc((void**)&tmp, cudamemsize)); //allocate memory and set to 0
			CUDA_API_ERRCHK(cudaMemset(tmp, 0, cudamemsize));

			runB<<< numiters / 256, 256 >>>(b, smin, ds, tmp); //iterate/get values
			CUDA_KERNEL_ERRCHK_WSYNC();

			CUDA_API_ERRCHK(cudaMemcpy(resOut.data(), tmp, cudamemsize, cudaMemcpyDeviceToHost));
			CUDA_API_ERRCHK(cudaFree(tmp));
		};

		//Iterate
		for (int iii = 0; iii < NUMITERS; iii++)
		{
			dipResults.at(iii) = dip->getBFieldAtS(DS * iii + SMIN, 0.0);
			lutResults.at(iii) = lut->getBFieldAtS(DS * iii + SMIN, 0.0);
		}

		runCUDAB(dip->getPtrGPU(), SMIN, DS, NUMITERS, dipGPUResults);
		runCUDAB(lut->getPtrGPU(), SMIN, DS, NUMITERS, lutGPUResults);

		//Check equality
		bool pass{ true };
		TEST_RESULTS_WITH_PASS("DipoleB CPU", fuzzyEq(dipChkVals, dipResults));
		TEST_RESULTS_WITH_PASS("DipoleB GPU", fuzzyEq(dipChkVals, dipGPUResults));
		TEST_RESULTS_WITH_PASS("DipoleBLUT CPU", fuzzyEq(lutChkVals, lutResults));
		TEST_RESULTS_WITH_PASS("DipoleBLUT GPU", fuzzyEq(lutChkVals, lutGPUResults));

		//Save
		if (save)
		{
			TEST_EXCEP_CHECK(writeDblBin(dipResults, "./../tests/data/dipole.bin", NUMITERS));
			TEST_EXCEP_CHECK(writeDblBin(lutResults, "./../tests/data/dipolelut.bin", NUMITERS));
		}

		return pass;
	}

	__global__ void runE(EField** e, double s, double* ret)
	{
		*ret = (*e)->getEFieldAtS(s, 0.0);
	}

	bool QSPSAccuracy()
	{
		std::unique_ptr<EField> e{ std::make_unique<EField>() };
		e->add(std::make_unique<QSPS>(std::vector<double>{4.0e6}, std::vector<double>{4.4e6}, std::vector<double>{1234.0}));
		e->add(std::make_unique<QSPS>(std::vector<double>{4.2e6}, std::vector<double>{4.6e6}, std::vector<double>{0.5678}));

		//hard-coded result values
		std::vector<double> alts{ 3.9e6, 4.0e6, 4.2e6, 4.5e6, 4.7e6 };
		std::vector<double> results{ 0.0, 1234.0, 1234.5678, 0.5678, 0.0 };

		std::vector<double> resultsCPU(5);
		std::vector<double> resultsGPU(5);

		double* results_d;
		CUDA_API_ERRCHK(cudaMalloc((void **)&results_d, sizeof(double) * 5));

		for (int alt = 0; alt < 5; alt++)
		{
			resultsCPU.at(alt) = e->getEFieldAtS(alts.at(alt), 0.0);
			runE<<< 1, 1 >>>(e->getPtrGPU(), alts.at(alt), results_d + alt);
			CUDA_KERNEL_ERRCHK_WSYNC();
		}

		CUDA_API_ERRCHK(cudaMemcpy(resultsGPU.data(), results_d, sizeof(double) * 5, cudaMemcpyDeviceToHost));
		CUDA_API_ERRCHK(cudaFree(results_d));

		bool pass{ true };
		TEST_RESULTS_WITH_PASS("EField (QSPS) CPU", fuzzyEq(results, resultsCPU));
		TEST_RESULTS_WITH_PASS("EField (QSPS) GPU", fuzzyEq(results, resultsGPU));
		return pass;
	}
}