#include <memory>
#include "utils\ionosphere.h"
#include "utils\ionosphereUtils.h"
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
constexpr int DSTNEBINS       { 96 };
constexpr int DSTNANGLEBINS   { 36000 };
constexpr double NFLUXIONRATIO{ 90.0 / 90.0 }; //numerator (90) is normal range (0-90), denominator is the range of the distribution
constexpr double NFLUXMAGRATIO{ 16.0 / 90.0 }; //so here, the distribution is from 16 degrees to zero (not 90 to 0), meaning the dist has more particles per angle
const std::string PARTNAME    { "elec" };
const std::string BTMSATNM    { "btmElec" };
const std::string UPGSATNM    { "4e6ElecUpg" };
const std::string DNGSATNM    { "4e6ElecDng" };

using std::string;
using ionosphere::Bins;
using ionosphere::EOMSimData;
using ionosphere::ParticleData;
using ionosphere::MaxwellianSpecs;
using ionosphere::IonosphereSpecs;
using utils::numerical::generateSpacedValues;

#define DUALREGIONLOGLINEFIT(bound, logm_upper, logb_upper, logm_lower, logb_lower) [](double s){ if (s > bound) return pow(10.0, logm_upper * s + logb_upper); else return pow(10.0, logm_lower * s + logb_lower); }

struct funcargs
{
	string simdatadir;
	bool cdf{ true };
};

funcargs parseMainArgs(int argc, char* argv[])
{
	funcargs args;

	if (argc < 2)
	{
		std::cout 
			<< "\nBinToCDF Usage:\n" << "\tBinToCDF.exe datadir [options]\n"
			<< "\n\tdatadir:\tThe directory where the data to process resides.\n\n"
			<< "\t[options]:\n"
			<< "\t--no-cdf:\tDo not output a cdf file after sim run.\n\n"
			<< "Exiting.\n";
		exit(1);
	}

	args.simdatadir = argv[1];
	args.simdatadir += "\\";

	for (int arg = 2; arg < argc; arg++)
	{
		if (string(argv[arg]) == string("--no-cdf"))
			args.cdf = false;
		else
		{
			std::cout << "Unrecognized flag: " << string(argv[arg]) << "  Exiting.\n";
			exit(1);
		}
	}

	return args;
}

int main(int argc, char* argv[])
{
	funcargs args{ parseMainArgs(argc, argv) };

	auto printVec = [](const std::vector<double>& x, int start = 0, int end = 0, int intvl = 1)
	{ //lambda which prints values from a vector
		if (end == 0) end = (int)x.size();
		for (int iii = start; iii < end; iii += intvl)
			std::cout << x.at(iii) << (((iii == (int)x.size() - 1) || iii + intvl == end) ? "\n" : ",");
	};

	auto err = [](double x, double y) { if (x == 0.0 && y == 0.0) return 0.0; return abs((x - y) / x); }; //check on error between two values

	// Form Maxwellian
	MaxwellianSpecs maxwellian(4.0 / 95.0); //dlogE of distribution - 4.0 / 95.0, dlogE of bins - 4.0 / 47.0
	maxwellian.push_back_ion(2.5,   6.00e7, 13500); //12000
	maxwellian.push_back_mag(2.5,   6.00e7, 18500 * 0.3); //4250
	maxwellian.push_back_mag(2.5e3, 1.20e8, 2000); //4250
	maxwellian.magModFactor = NFLUXMAGRATIO; //pitch angle space density difference from ionosphere, pitch range is from 0-16, not 0-90

	// Form Postprocessing Data
	Bins distbins(generateSpacedValues(0.5, 4.5, DSTNEBINS, true, true), generateSpacedValues(179.9975, 0.0025, DSTNANGLEBINS, false, true));
	Bins satbins (generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true), generateSpacedValues(5.0, 175.0, CDFNANGLEBINS, false, true));


	auto O = [](double s) {
		if (s > 135000.0)
			return pow(10.0, (4.0 - log10(3.0e10)) / (800000.0 - 135000.0) * s + 11.79203);
		else if (s > 98000.0)
			return pow(10.0, (log10(3.0e10) - log10(6.0e11)) / (135000.0 - 98000.0) * s + 15.22412);
		else
			return pow(10.0, (log10(6.0e11) - 8.0) / (98000.0 - 76000.0) * s - 5.0518);
	};

	/* Physical ionosphere model */
	IonosphereSpecs ionsph(58, 620000.0, 50000.0);
	//Ionosphere ionsph(232, 620000.0, 50000.0);
	ionsph.addSpecies("N2", 14.0, DUALREGIONLOGLINEFIT(130000.0, (4.0 - 11.0) / (540000.0 - 130000.0), 13.21951, (11.0 - 19.0) / (130000 - 8000), 19.52459));
	ionsph.addSpecies("He", 3.0,  DUALREGIONLOGLINEFIT(120000.0, (6.0 - log10(5.0e7)) / (1000000.0 - 120000.0), 7.930648, (log10(5.0e7) - 14) / (120000.0), 14.0));
	ionsph.addSpecies("O2", 16.0, DUALREGIONLOGLINEFIT(130000.0, (10.0 - 4.0) / (130000.0 - 450000.0), 12.4375, (18.0 - 10.0) / (20000.0 - 130000.0), 19.45455));
	ionsph.addSpecies("O",  8.0,  O);

	/* Test ionosphere where everything scatters */
	//Ionosphere ionsph(2, 620000.0, 619999.9999);
	//ionsph.addSpecies("ScatterAll", 1.0e6, [](double s) { return 1.0e30; });
	/* End test ionosphere */


	EOMSimData eomdata{ ionsph, maxwellian, distbins, satbins,
		args.simdatadir, PARTNAME, BTMSATNM, UPGSATNM, DNGSATNM };

	//printVec(ppdata.maxWeights, 0, 96);
	//std::cout << "\n";
	//printVec(ppdata.maxWeights, 1728000, 1728096);
	//exit(1);

	// Run Post Process Code
	std::vector<std::vector<double>> fluxData;
	SIM_API_EXCEP_CHECK(fluxData = steadyFlux(eomdata));


	if (args.cdf)
	{
		/* Prep Data for CDF - For some reason, array of vector.data() doesn't work */
		double cntArray2D[CDFNANGLEBINS][CDFNEBINS];
		for (int ang = 0; ang < CDFNANGLEBINS; ang++)
			for (int eng = 0; eng < CDFNEBINS; eng++)
				cntArray2D[ang][eng] = fluxData.at(ang).at(eng);

		/* Create CDF file and setup with appropriate variables, write data */
		std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

		cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS }, (void*)eomdata.satbins.E.data());
		cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, (void*)satbins.PA.data());
		cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS }, cntArray2D);

		//cdf is written to disk on class destruction
	}

	return 0;
}