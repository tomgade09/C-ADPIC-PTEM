#include <memory>
#include "utils\postprocess.h"
#include "CDFFileClass.h"
#include "ErrorHandling/simExceptionMacros.h"
#include "utils\numerical.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>

constexpr int CDFNEBINS       { 48 };
constexpr int CDFNANGLEBINS   { 18 };
constexpr double NFLUXIONRATIO{ 90.0 / 90.0 }; //numerator (90) is normal range (0-90), denominator is the range of the distribution
constexpr double NFLUXMAGRATIO{ 16.0 / 90.0 }; //so here, the distribution is from 16 degrees to zero (not 90 to 0), meaning the dist has more particles per angle
//const std::string SIMDATADIR  { "./../_dataout/180716_14.09.22/" };
const std::string PARTNAME    { "elec" };
const std::string BTMSATNM    { "btmElec" };
const std::string UPGSATNM    { "4e6ElecDown" }; //"Down" is downward facing detector (so upgoing), not downgoing particles
const std::string DNGSATNM    { "4e6ElecUp" }; //same


using postprocess::PPData;
using postprocess::ParticleData;
using postprocess::Maxwellian;
using utils::numerical::generateSpacedValues;

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "\nBinToCDF Usage:\n" << "\tBinToCDF.exe datadir\n" << "\n\tdatadir:\tThe directory where the data to process resides.\nExiting.\n";
		exit(1);
	}
	std::string simdatadir{ argv[1] };
	simdatadir += "\\";
	
	auto printVec = [](const std::vector<double>& x, int start = 0, int end = 0, int intvl = 1)
	{ //lambda which prints values from a vector
		if (end == 0) end = (int)x.size();
		for (int iii = start; iii < end; iii += intvl)
			std::cout << x.at(iii) << (((iii == (int)x.size() - 1) || iii + intvl == end) ? "\n" : ",");
	};

	auto err = [](double x, double y) { if (x == 0.0 && y == 0.0) return 0.0; return abs((x - y) / x); }; //check on error between two values

	Maxwellian maxwellian(4.0 / 95.0); //dlogE of distribution - 4.0 / 95.0, dlogE of bins - 4.0 / 47.0
	maxwellian.push_back_ion(10.0,  7.00e7, 5000);
	maxwellian.push_back_mag(10.0,  2.00e7, 5000);
	maxwellian.push_back_mag(5.0e3, 1.00e8, 5000);
	maxwellian.magModFactor = NFLUXMAGRATIO; // * magcm2Ratio //pitch angle space density difference from ionosphere, pitch range is from 0-16, not 0-90

	PPData ppdata{ maxwellian, generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true), generateSpacedValues(5.0, 175.0, CDFNANGLEBINS, false, true),
		simdatadir, PARTNAME, BTMSATNM, UPGSATNM, DNGSATNM };

	std::vector<std::vector<double>> fluxData;
	SIM_API_EXCEP_CHECK(fluxData = postprocess::steadyFlux(ppdata));


	/* Prep Data for CDF - For some reason, array of vector.data() doesn't work */
	double cntArray2D[CDFNANGLEBINS][CDFNEBINS];
	for (int ang = 0; ang < CDFNANGLEBINS; ang++)
		for (int eng = 0; eng < CDFNEBINS; eng++)
			cntArray2D[ang][eng] = fluxData.at(ang).at(eng);

	/* Create CDF file and setup with appropriate variables, write data */
	std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

	cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS }, (void*)ppdata.energyBins.data());
	cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, (void*)ppdata.pitchBins.data());
	cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS }, cntArray2D);

	return 0;
}