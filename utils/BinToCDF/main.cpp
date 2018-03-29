#include <memory>
#include "postprocess.h"
#include "CDFFileClass.h"

constexpr int CDFNEBINS{ 48 };
constexpr int CDFNANGLEBINS{ 18 };
constexpr int SIMNEBINS{ 96 };
constexpr int SIMNANGLEBINS{ 3600 };
constexpr int PARTICLECOUNT{ SIMNEBINS * SIMNANGLEBINS };

int main()
{
	std::vector<double> binAngles;
	std::vector<double> binEnergies;
	vecDbl2D fluxData;
	TRYCATCHSTDEXP(fluxData = postprocess::steadyFlux("./../../../../../_dataout/180328_16.42.53/", { 0.0, 180.0 }, CDFNANGLEBINS, { 0.5, 4.5 }, CDFNEBINS, { 10.0, 1.25e5 }, { 5.0e3, 7.0e5, 10.0, 1.5e5 }, 9.10938356e-31, PARTICLECOUNT, binAngles, binEnergies));

	/* Prep Data for CDF */
	double cntArray2D[CDFNANGLEBINS][CDFNEBINS];
	for (int ang = 0; ang < CDFNANGLEBINS; ang++)
		for (int eng = 0; eng < CDFNEBINS; eng++)
			cntArray2D[ang][eng] = fluxData.at(ang).at(eng);

	/* Create CDF file and setup with appropriate variables, write data */
	std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

	cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS }, binEnergies.data());
	cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, binAngles.data());
	cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS }, cntArray2D);

	return 0;
}