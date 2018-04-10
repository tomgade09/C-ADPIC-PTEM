#include "utils\utilsAPI.h"

DLLEXPORT utils::write::ParticleDistribution* createParticleDistributionAPI(const char* saveFolder, const char* attrNames, const char* particleName, double mass)
{
	SIM_API_EXCEP_CHECK(
	std::string save{ saveFolder };
	std::vector<std::string> attrs{ utils::string::charToStrVec(attrNames) };
	std::string part{ particleName };

	utils::write::ParticleDistribution* ret;
	if (save == "" && attrs.size() == 0 && part == "" && mass == 0.0)
		ret = new utils::write::ParticleDistribution();
	else
		ret = new utils::write::ParticleDistribution(save, attrs, part, mass);

	return ret;
	);
}

DLLEXPORT void addPDEnergyRangeAPI(utils::write::ParticleDistribution* pd, int energyBins, double Emin, double Emax, bool logE) {
	SIM_API_EXCEP_CHECK(pd->addEnergyRange(energyBins, Emin, Emax, logE)); }

DLLEXPORT void addPDPitchRangeAPI(utils::write::ParticleDistribution* pd, int pitchBins, double PAmin, double PAmax, bool midBin) {
	SIM_API_EXCEP_CHECK(pd->addPitchRange(pitchBins, PAmin, PAmax, midBin)); }

DLLEXPORT void generatePDAPI(utils::write::ParticleDistribution* pd, double s_ion, double s_mag) {
	SIM_API_EXCEP_CHECK(pd->generate(s_ion, s_mag)); }

DLLEXPORT void writePDAPI(utils::write::ParticleDistribution* pd) {
	delete pd; }