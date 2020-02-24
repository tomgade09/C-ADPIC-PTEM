#include <memory>
#include "ionosphere/ionosphere.h"
#include "ionosphere/ionosphereUtils.h"
#include "CDFFileClass.h"
#include "ErrorHandling/simExceptionMacros.h"
#include "utils/numerical.h"
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
constexpr double NFLUXMAGRATIO{ 90.0 / 90.0 }; //so here, the distribution is from 16 degrees to zero (not 90 to 0), meaning the dist has more particles per angle
const std::string PARTNAME    { "elec" };
const std::string BTMSATNM    { "btmElec" };
const std::string UPGSATNM    { "4e6ElecUpg" };
const std::string DNGSATNM    { "4e6ElecDng" };

using std::string;
using std::logic_error;
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

namespace debug
{
	void setIdealSatDists(EOMSimData& eomdata)
	{
		ParticleData btm(eomdata.bottom.energy.size());
		ParticleData upw(eomdata.upward.energy.size());
		ParticleData dnw(eomdata.dnward.energy.size());

		for (int part = 0; part < eomdata.initial.energy.size(); part++)
		{
			if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)
			{
				double satPA{ ionosphere::dEFlux::newPA(eomdata.initial.pitch.at(part), eomdata.B_ion, eomdata.B_sat) };

				if (satPA > 0)
				{
					upw.pitch.at(part) = satPA;
					upw.energy.at(part) = eomdata.initial.energy.at(part);
				}
			}
			else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)
			{
				double satPA{ ionosphere::dEFlux::newPA(eomdata.initial.pitch.at(part), eomdata.B_mag, eomdata.B_sat) };
				double btmPA{ ionosphere::dEFlux::newPA(eomdata.initial.pitch.at(part), eomdata.B_mag, eomdata.B_ion/*-0.0000587126*/) };

				if (satPA > 0)
				{
					dnw.pitch.at(part) = satPA;
					dnw.energy.at(part) = eomdata.initial.energy.at(part);
					if (btmPA < 0)
					{
						upw.pitch.at(part) = 180.0 - satPA;
						upw.energy.at(part) = eomdata.initial.energy.at(part);
					}
				}
				if (btmPA > 0)
				{
					btm.pitch.at(part) = btmPA;
					btm.energy.at(part) = eomdata.initial.energy.at(part);
				}
			}
			else
			{
				throw logic_error("debug::generateIdealSatDists : particle is not ionospheric or magnetospheric source");
			}
		}

		eomdata.bottom = btm;
		eomdata.upward = upw;
		eomdata.dnward = dnw;
	}

	void setRealMaxwellians(EOMSimData& eomdata)
	{
		vector<double> realMagMaxwellian = {
			511.7803363,
			455.366866,
			404.1659709,
			357.6960119,
			315.5198528,
			277.2407488,
			242.498614,
			210.9666344,
			182.3481933,
			156.3740811,
			132.7999634,
			111.4040819,
			91.7746781,
			73.55799355,
			57.02452115,
			42.01873297,
			34.36415945,
			28.31990953,
			22.99218389,
			18.56252397,
			14.73461663,
			11.99012355,
			9.875617326,
			8.350030687,
			7.817520858,
			7.352563096,
			7.057355075,
			6.732956547,
			5.444957563,
			4.258616122,
			3.148063762,
			2.667671894,
			2.31874741,
			2.413216721,
			2.510346733,
			2.16919973,
			1.822258224,
			1.538644415,
			1.454668358,
			1.448302276,
			1.421289335,
			1.392400846,
			1.356434703,
			1.3223143,
			1.340996193,
			1.24936111,
			1.082697097,
			1.027704468,
			1.022389203,
			0.954603057,
			0.853591162,
			0.787414014,
			0.712480444,
			0.618899297,
			0.613108903,
			0.676524001,
			0.741544197,
			0.809125578,
			0.826634801,
			0.844583081,
			0.909628356,
			0.96415381,
			1.006445782,
			1.033757784,
			1.034400169 * 1.480727,
			1.042600148 * 1.480727,
			1.055748627 * 1.428545,
			1.07152601  * 1.428545,
			1.078133553 * 1.470468,
			1.058734323 * 1.470468,
			1.039323547 * 1.438624,
			1.012305997 * 1.438624,
			0.994205528 * 1.442204,
			0.985836018 * 1.442204,
			0.910594196 * 1.312604,
			0.838235396 * 1.312604,
			0.759978758 * 1.234802,
			0.688727757 * 1.234802,
			0.582535504 * 1.223116,
			0.484989966 * 1.223116,
			0.393631204 * 1.319446,
			0.330308124 * 1.319446,
			0.27732655  * 1.287453,
			0.225898359 * 1.287453,
			0.178965883, //84
			0.142250867,
			0.110109027,
			0.082409802,
			0.060637842,
			0.042555514,
			0.027887484,
			0.0150541,
			0.008638964,
			0.004889727,
			0.002865042,
			0.010868959
		};

		vector<double> realIonMaxwellian = {
			628.8979407 / 2.4786570,
			565.8638197 / 2.4786570,
			508.6540186 / 2.1904623,
			456.7303733 / 2.1904623,
			409.6044458 / 1.9406344,
			366.8329293 / 1.9406344,
			328.0134787 / 1.7238276,
			292.780925  / 1.7238276,
			260.8038409 / 1.5356057,
			231.7814228 / 1.5356057,
			205.4406609 / 1.3725660,
			181.5337716 / 1.3725660,
			159.4832008 / 1.2364910,
			138.7981913 / 1.2364910,
			120.0244659 / 1.1444663,
			102.9854228 / 1.1444663,
			89.92642941 / 1.0706112,
			78.43829219 / 1.0706112,
			67.92556403,
			58.16316132,
			49.00339734,
			39.55476535,
			32.43565065,
			27.49226718,
			23.65374146,
			20.06519194,
			16.08473353,
			12.55060252,
			10.72488715,
			9.128684634,
			7.798530255,
			6.437362853,
			5.176560025,
			4.356319064,
			3.620845049,
			3.394283882,
			3.146182049,
			2.535331134,
			2.166546899,
			1.906599314,
			1.725538745,
			1.572503401,
			1.42385702,
			1.291013895,
			1.33179657,
			1.213164987,
			0.98581798,
			0.891060771,
			0.856747842,
			0.813801222,
			0.767419376,
			0.762158834,
			0.727575643,
			0.647540963,
			0.613472219,
			0.616017952,
			0.584045464,
			0.515672013,
			0.529429898,
			0.548721841,
			0.510828937,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0
		};
		
		//vector<double> realIonMaxwellian(96, 0);

		for (size_t part = 0; part < eomdata.initial.pitch.size(); part++)
		{
			if (eomdata.initial.pitch.at(part) < 90.0)
			{
				eomdata.maxwellian.at(part) = realMagMaxwellian.at(part % 96)/* * NFLUXMAGRATIO*/;
			}
			else
			{
				eomdata.maxwellian.at(part) = realIonMaxwellian.at(part % 96) / NFLUXMAGRATIO /* * NFLUXIONRATIO*/;
			}
		}
	}
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
	maxwellian.push_back_mag(2.5,   6.00e7, 18500 * 3 / 10); //4250
	maxwellian.push_back_mag(2.5e3, 1.20e8, 2000); //4250
	maxwellian.magModFactor = NFLUXMAGRATIO; //pitch angle space density difference from ionosphere, pitch range is from 0-16, not 0-90

	
	// Form Postprocessing Data
	Bins distbins(generateSpacedValues(0.5, 4.5, DSTNEBINS, true, true), generateSpacedValues(179.9975, 0.0025, DSTNANGLEBINS, false, true));
	Bins satbins (generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true), generateSpacedValues(/*5.0, 175.0*/180.0 / (CDFNANGLEBINS) / 2,180.0 - 180.0 / (CDFNANGLEBINS) / 2, CDFNANGLEBINS, false, true));
	
	//
	//
	// Remove ~3.1 to 5ish eV bins (to more closely approximate FAST)
	vector<double> newEsatbin;
	for (size_t bin = 3; bin < satbins.E.size(); bin++)
		newEsatbin.push_back(satbins.E.at(bin));

	satbins.E = newEsatbin;
	//
	//
	//

	/* Physical ionosphere model */
	auto O = [](double s) {
		if (s > 135000.0)
			return pow(10.0, (4.0 - log10(3.0e10)) / (800000.0 - 135000.0) * s + 11.79203);
		else if (s > 98000.0)
			return pow(10.0, (log10(3.0e10) - log10(6.0e11)) / (135000.0 - 98000.0) * s + 15.22412);
		else
			return pow(10.0, (log10(6.0e11) - 8.0) / (98000.0 - 76000.0) * s - 5.0518);
	};

	IonosphereSpecs ionsph(58, 620000.0, 50000.0); //this is in alt - converted to s (dist along field line) once it's inserted into eomdata below (constructor handles this)
	//Ionosphere ionsph(232, 620000.0, 50000.0);
	//These functions return number density in cm^-3
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

	//debug::setIdealSatDists(eomdata); //set non-time dependent equation - calculated distribution
	debug::setRealMaxwellians(eomdata);


	// Run Post Process Code
	std::vector<std::vector<double>> fluxData;
	SIM_API_EXCEP_CHECK(fluxData = steadyFlux(eomdata));


	if (args.cdf)
	{
		/* Prep Data for CDF - For some reason, array of vector.data() doesn't work */
		double cntArray2D[CDFNANGLEBINS][CDFNEBINS - 3];

		for (int ang = 0; ang < CDFNANGLEBINS; ang++)
			for (int eng = 0; eng < CDFNEBINS - 3; eng++)
				cntArray2D[ang][eng] = fluxData.at(ang).at(eng);

		/* Create CDF file and setup with appropriate variables, write data */
		std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

		cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS - 3 }, (void*)eomdata.satbins.E.data());
		cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, (void*)satbins.PA.data());
		cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS - 3 }, cntArray2D);

		//cdf is written to disk on class destruction
	}

	return 0;
}