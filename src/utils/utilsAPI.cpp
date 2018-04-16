#include "utils\utilsAPI.h"

//ParticleDistribution functions
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

DLLEXPORT void padExtraPDAttrsAPI(utils::write::ParticleDistribution* pd)
{
	SIM_API_EXCEP_CHECK(
		int attrssize{ (int)pd->data().size() };
		int partsize{ (int)pd->data().at(0).size() };
		std::vector<double> zeroes(partsize);
		for (int iii = 0; iii < attrssize - 3; iii++)
			pd->setattr(zeroes, 3 + iii);
	);
}

DLLEXPORT void writePDAPI(utils::write::ParticleDistribution* pd) {
	delete pd; }


//DistributionFromDisk functions
DLLEXPORT utils::load::DistributionFromDisk* loadDistributionFromDiskAPI(const char* name, const char* loadFolder, const char* attrNames, const char* particleName) {
	SIM_API_EXCEP_CHECK(return new utils::load::DistributionFromDisk(name, loadFolder, particleName, utils::string::charToStrVec(attrNames))); }

DLLEXPORT void DistFromDiskPrintAPI(utils::load::DistributionFromDisk* dfd, int at) {
	SIM_API_EXCEP_CHECK(dfd->print(at)); }

DLLEXPORT void DistFromDiskPrintDiffAPI(utils::load::DistributionFromDisk* dfd_this, utils::load::DistributionFromDisk* dfd_other, int at) {
	SIM_API_EXCEP_CHECK(dfd_this->printdiff(*dfd_other, at)); }

DLLEXPORT void DistFromDiskZeroesAPI(utils::load::DistributionFromDisk* dfd) {
	SIM_API_EXCEP_CHECK(dfd->zeroes()); }

DLLEXPORT void DistFromDiskCompareAPI(utils::load::DistributionFromDisk* dfd_this, utils::load::DistributionFromDisk* dfd_other) {
	SIM_API_EXCEP_CHECK(dfd_this->compare(*dfd_other)); }

DLLEXPORT void deleteDistFromDiskAPI(utils::load::DistributionFromDisk* dfd) {
	delete dfd; }