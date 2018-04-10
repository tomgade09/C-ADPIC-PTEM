#ifndef UTILS_API_H
#define UTILS_API_H

#include "dllexport.h"
#include "utils\write.h"
#include "utils\string.h"
#include "ErrorHandling\simExceptionMacros.h"

DLLEXPORT utils::write::ParticleDistribution* createParticleDistributionAPI(const char* saveFolder, const char* attrNames, const char* particleName, double mass);
DLLEXPORT void addPDEnergyRangeAPI(utils::write::ParticleDistribution* pd, int energyBins, double Emin, double Emax, bool logE = true);
DLLEXPORT void addPDPitchRangeAPI(utils::write::ParticleDistribution* pd, int pitchBins, double PAmin, double PAmax, bool midBin = true);
DLLEXPORT void generatePDAPI(utils::write::ParticleDistribution* pd, double s_ion, double s_mag);
DLLEXPORT void writePDAPI(utils::write::ParticleDistribution* pd);

#endif /* UTILS_API_H */