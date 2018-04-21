#include <memory>
#include "SimulationClass\Simulation.h"
#include "postprocess.h"
#include "CDFFileClass.h"
#include "utils\numerical.h"

constexpr int CDFNEBINS{ 48 };
constexpr int CDFNANGLEBINS{ 18 };

int main()
{
	std::unique_ptr<Simulation> sim{ std::make_unique<Simulation>("./../../../../../_dataout/180419_15.13.26/", true) };

	std::vector<double> energyBins{ utils::numerical::generateSpacedValues(0.5, 4.5, CDFNEBINS, true, true) };
	std::vector<double> pitchBins{ utils::numerical::generateSpacedValues(5.0, 175.0, CDFNANGLEBINS, false, true) };

	postprocess::ParticleData init{ sim->particle(0)->getOrigData().at(0), sim->particle(0)->getOrigData().at(1), sim->particle(0)->mass() };
	init.s_pos = sim->particle(0)->getOrigData().at(2);
	postprocess::ParticleData btm { sim->satellite("btmElec")->data().at(0).at(0), sim->satellite("btmElec")->data().at(0).at(1), sim->particle("elec")->mass }; //for now, don't need t_esc, but maybe later
	postprocess::ParticleData up  { sim->satellite("4e6ElecDown")->data().at(0).at(0), sim->satellite("4e6ElecDown")->data().at(0).at(1), sim->particle("elec")->mass };
	postprocess::ParticleData dn  { sim->satellite("4e6ElecUp")->data().at(0).at(0), sim->satellite("4e6ElecUp")->data().at(0).at(1), sim->particle("elec")->mass };

	postprocess::MaxwellianData maxData;
	maxData.push_back_ion(10.0,  1.25e5);
	maxData.push_back_mag(10.0,  1.50e5 / 56);
	maxData.push_back_mag(5.0e3, 7.00e5 / 56);

	std::vector<double> max{ postprocess::maxwellian::formCountsVector(init, maxData, sim->simMin, sim->simMax, 4.0 / 95.0) };

	postprocess::PPData ppdata{ sim->simMin, sim->simMax, energyBins, pitchBins, max, init, btm, up, dn };
	std::vector<std::vector<double>> fluxData;
	fluxData = postprocess::steadyFlux(ppdata);

	exit(1);

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