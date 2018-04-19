#ifndef UTILS_API_H
#define UTILS_API_H

#include "dllexport.h"
#include "utils\load.h"
#include "utils\write.h"
#include "utils\string.h"
#include "ErrorHandling\simExceptionMacros.h"

//ParticleDistribution functions
DLLEXPORT utils::write::ParticleDistribution* createParticleDistributionAPI(const char* saveFolder, const char* attrNames, const char* particleName, double mass);
DLLEXPORT void addPDEnergyRangeAPI(utils::write::ParticleDistribution* pd, int energyBins, double Emin, double Emax, bool logE = true);
DLLEXPORT void addPDPitchRangeAPI(utils::write::ParticleDistribution* pd, int pitchBins, double PAmin, double PAmax, bool midBin = true);
DLLEXPORT void generatePDAPI(utils::write::ParticleDistribution* pd, double s_ion, double s_mag);
DLLEXPORT void padExtraPDAttrsAPI(utils::write::ParticleDistribution* pd);
DLLEXPORT void writePDAPI(utils::write::ParticleDistribution* pd);

//DistributionFromDisk functions
DLLEXPORT utils::load::DistributionFromDisk* loadDistributionFromDiskAPI(const char* name, const char* loadFolder, const char* attrNames, const char* particleName);
DLLEXPORT const double* DistFromDiskDataAPI(utils::load::DistributionFromDisk* dfd, int attrInd);
DLLEXPORT void DistFromDiskPrintAPI(utils::load::DistributionFromDisk* dfd, int at);
DLLEXPORT void DistFromDiskPrintDiffAPI(utils::load::DistributionFromDisk* dfd_this, utils::load::DistributionFromDisk* dfd_other, int at);
DLLEXPORT void DistFromDiskZeroesAPI(utils::load::DistributionFromDisk* dfd);
DLLEXPORT void DistFromDiskCompareAPI(utils::load::DistributionFromDisk* dfd_this, utils::load::DistributionFromDisk* dfd_other);
DLLEXPORT void deleteDistFromDiskAPI(utils::load::DistributionFromDisk* dfd);

#endif /* UTILS_API_H */