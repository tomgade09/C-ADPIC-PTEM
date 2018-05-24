#include <memory>
#include "Simulation\Simulation.h"
#include "ErrorHandling\simExceptionMacros.h"
#include "utils\postprocess.h"
#include "CDFFileClass.h"
#include "utils\numerical.h"
#include "API\SimulationAPI.h"
#include "utils\silenceStreamMacros.h"
#include <sstream>
#include <algorithm>
#include <iterator>

constexpr int CDFNEBINS{ 48 };
constexpr int CDFNANGLEBINS{ 18 };

using postprocess::PPData;
using postprocess::ParticleData;
using postprocess::MaxwellianData;
using utils::numerical::generateSpacedValues;

int main()
{
	/*std::vector<double> pitches(180);
	std::vector<double> energies(96);
	for (int iii = 0; iii < 96; iii++)
		energies.at(iii) = pow(10, iii * 4.0 / 95.0 + 0.5);
	for (int iii = 0; iii < 180; iii++)
		pitches.at(iii) = (double)iii;

	double dangle{ 10.0 };
	double dlogE{ 4.0 / 47.0 };

	for (int iii = 0; iii < 96; iii++)
	{
		int engbin{ (int)((log10(energies.at(iii)) - (0.5 - 0.5 * dlogE)) / dlogE) }; //ditto

		std::cout << engbin << "  ";
	}

	std::cout << "\n";

	for (int iii = 0; iii < 180; iii++)
	{
		int angbin{ (int)(pitches.at(iii) / dangle) }; //this should give the bin index

		std::cout << angbin << "  ";
	}

	std::cout << "\n";
	

	exit(1);*/

	auto printVec = [](const std::vector<double>& x, int start = 0, int end = 0, int intvl = 1)
	{
		if (end == 0) end = (int)x.size();
		for (int iii = start; iii < end; iii += intvl)
			std::cout << x.at(iii) << ((iii != (int)x.size() - 1) ? "," : "\n");
	};

	std::unique_ptr<Simulation> sim;
	
	SILENCE_COUT(SIM_API_EXCEP_CHECK(sim = std::make_unique<Simulation>("./../_dataout/180423_09.34.25/")));

	std::vector<double> energyBins;
	std::vector<double> pitchBins;

	SIM_API_EXCEP_CHECK(energyBins = generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true));
	SIM_API_EXCEP_CHECK(pitchBins  = generateSpacedValues(5.0, 175.0, CDFNANGLEBINS, false, true));

	std::cout << "test:        " << sim->particle("elec")->data(true).at(0).at(0) << " " << sim->particle("elec")->data(true).at(1).at(0) << " " << sim->particle("elec")->data(true).at(2).at(0) << " " << sim->particle("elec")->mass() << std::endl;
	std::cout << "btmelec:     " << sim->satellite("btmElec")->data().at(0).at(0).at(3455000) << " " << sim->satellite("btmElec")->data().at(0).at(1).at(3455000) << " " << sim->particle("elec")->mass() << std::endl;
	std::cout << "4e6ElecDown: " << sim->satellite("4e6ElecDown")->data().at(0).at(0).at(0) << " " << sim->satellite("4e6ElecDown")->data().at(0).at(1).at(0) << " " << sim->particle("elec")->mass() << std::endl;
	std::cout << "4e6ElecUp:   " << sim->satellite("4e6ElecUp")->data().at(0).at(0).at(3455000) << " " << sim->satellite("4e6ElecUp")->data().at(0).at(1).at(3455000) << " " << sim->particle("elec")->mass() << std::endl;


	ParticleData init{ sim->particle("elec")->data(true).at(0), sim->particle("elec")->data(true).at(1), sim->particle("elec")->mass() };
	init.s_pos = sim->particle("elec")->data(true).at(2);
	ParticleData btm { sim->satellite("btmElec")->data().at(0).at(0), sim->satellite("btmElec")->data().at(0).at(1), sim->particle("elec")->mass() }; //for now, don't need t_esc, but maybe later
	ParticleData up  { sim->satellite("4e6ElecDown")->data().at(0).at(0), sim->satellite("4e6ElecDown")->data().at(0).at(1), sim->particle("elec")->mass() };
	ParticleData dn  { sim->satellite("4e6ElecUp")->data().at(0).at(0), sim->satellite("4e6ElecUp")->data().at(0).at(1), sim->particle("elec")->mass() };

	//std::cout << "init: " << init.energy.at(0) << " " << init.pitch.at(0) << "\n";
	std::cout << "init: " << init.energy.at(3450000) << " " << init.pitch.at(3450000) << "\n";
	std::cout << "btm : " << btm.energy.at(3450000) << " " << btm.pitch.at(3450000) << "\n";
	//std::cout << "up  : " << up.energy.at(0) << " " << up.pitch.at(0) << "\n";
	std::cout << "dn  : " << dn.energy.at(3450000) << " " << dn.pitch.at(3450000) << "\n";

	exit(1);
	auto err = [](double x, double y) { if (x == 0.0 && y == 0.0) return 0.0; return abs((x - y) / x); };

	for (int iii = 0; iii < 96 * 36000; iii++)
		if (err(init.energy.at(iii), dn.energy.at(iii)) > 1.0e-4)
			std::cout << init.energy.at(iii) << " : " << dn.energy.at(iii) << " : " << err(init.energy.at(iii), dn.energy.at(iii)) << std::endl;

	exit(1);

	double alt{ sim->satellite("4e6ElecDown")->altitude() };
	double ioncm2Ratio{ sqrt(sim->Bmodel()->getBFieldAtS(alt, 0.0) / sim->Bmodel()->getBFieldAtS(sim->simMin(), 0.0)) };
	double magcm2Ratio{ sqrt(sim->Bmodel()->getBFieldAtS(alt, 0.0) / sim->Bmodel()->getBFieldAtS(sim->simMax(), 0.0)) };

	MaxwellianData maxData(sim->simMin(), sim->simMax(), 4.0 / 95.0);
	maxData.push_back_ion(10.0,  7.00e7 /* * ioncm2Ratio */);
	maxData.push_back_mag(10.0,  2.00e7 /* * magcm2Ratio *// 5.625); //pitch angle space density difference from ionosphere
	maxData.push_back_mag(5.0e3, 1.00e8 /* * magcm2Ratio *// 5.625); //pitch range is from 0-16, not 0-90

	PPData ppdata{ sim->simMin(), sim->simMax(), sim->getBFieldAtS(sim->simMin(), 0.0), sim->getBFieldAtS(alt, 0.0), sim->getBFieldAtS(sim->simMax(), 0.0),
		energyBins, pitchBins, maxData, init, btm, up, dn };

	/*printVec(init.energy, 0, 96, 1);
	std::cout << std::endl << std::endl;
	printVec(init.pitch, 0, 96*96, 96);
	std::cout << std::endl << std::endl;
	printVec(btm.energy, 0, 96, 1);

	exit(1);*/
	//printVec(ppdata.maxCounts, 0, 96);
	//std::cout << std::endl << std::endl;
	//printVec(ppdata.maxCounts, (int)ppdata.initial.energy.size() / 2, (int)ppdata.initial.energy.size() / 2 + 96, 1);
	//std::cout << std::endl << std::endl;
	//printVec(ppdata.maxCounts, 0, 1000, 96);
	//exit(1);

	std::vector<std::vector<double>> fluxData;
	SIM_API_EXCEP_CHECK(fluxData = postprocess::steadyFlux(ppdata));



	/* Prep Data for CDF - For some reason, array of vector.data() doesn't work */
	double cntArray2D[CDFNANGLEBINS][CDFNEBINS];
	for (int ang = 0; ang < CDFNANGLEBINS; ang++)
		for (int eng = 0; eng < CDFNEBINS; eng++)
			cntArray2D[ang][eng] = fluxData.at(ang).at(eng);

	/* Create CDF file and setup with appropriate variables, write data */
	std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

	cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS }, energyBins.data());
	cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, pitchBins.data());
	cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS }, cntArray2D);

	return 0;
}