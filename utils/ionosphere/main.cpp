#include <memory>
#include "Simulation\Simulation.h"
#include "utils\postprocess.h"
#include "src\CDFFileClass.h"
#include "utils\numerical.h"
#include "API\SimulationAPI.h"

constexpr int CDFNEBINS{ 48 };
constexpr int CDFNANGLEBINS{ 18 };

int main()
{
	//std::unique_ptr<Simulation> sim{ std::make_unique<Simulation>("./../../../../../_dataout/180419_15.13.26/") };
	Simulation* sim{ loadCompletedSimDataAPI("./../../../../../_dataout/180419_15.13.26/") };

	std::vector<double> energyBins{ utils::numerical::generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true) };
	std::vector<double> pitchBins{ utils::numerical::generateSpacedValues(5.0, 175.0, CDFNANGLEBINS, false, true) };

	postprocess::ParticleData init{ sim->particle(0)->data(true).at(0), sim->particle(0)->data(true).at(1), sim->particle(0)->mass() };
	init.s_pos = sim->particle(0)->data(true).at(2);
	postprocess::ParticleData btm { sim->satellite("btmElec")->data().at(0).at(0), sim->satellite("btmElec")->data().at(0).at(1), sim->particle("elec")->mass() }; //for now, don't need t_esc, but maybe later
	postprocess::ParticleData up  { sim->satellite("4e6ElecDown")->data().at(0).at(0), sim->satellite("4e6ElecDown")->data().at(0).at(1), sim->particle("elec")->mass() };
	postprocess::ParticleData dn  { sim->satellite("4e6ElecUp")->data().at(0).at(0), sim->satellite("4e6ElecUp")->data().at(0).at(1), sim->particle("elec")->mass() };

	double alt{ sim->satellite("4e6ElecDown")->altitude() };
	double ioncm2Ratio{ sqrt(sim->Bmodel()->getBFieldAtS(alt, 0.0) / sim->Bmodel()->getBFieldAtS(sim->simMin(), 0.0)) };
	double magcm2Ratio{ sqrt(sim->Bmodel()->getBFieldAtS(alt, 0.0) / sim->Bmodel()->getBFieldAtS(sim->simMax(), 0.0)) };

	postprocess::MaxwellianData maxData(sim->simMin(), sim->simMax(), 4.0 / 95.0);
	maxData.push_back_ion(10.0,  7.00e7 * ioncm2Ratio);
	maxData.push_back_mag(10.0,  2.00e7 * magcm2Ratio / 5.625); //pitch angle space density difference from ionosphere
	maxData.push_back_mag(5.0e3, 1.00e8 * magcm2Ratio / 5.625); //pitch range is from 0-16, not 0-90

	postprocess::PPData ppdata{ sim->simMin(), sim->simMax(), sim->getBFieldAtS(sim->simMin(), 0.0),
		sim->getBFieldAtS(alt, 0.0), sim->getBFieldAtS(sim->simMax(), 0.0), energyBins, pitchBins, maxData, init, btm, up, dn };
	std::vector<std::vector<double>> fluxData;
	fluxData = postprocess::steadyFlux(ppdata);

	terminateSimulationAPI(sim);

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