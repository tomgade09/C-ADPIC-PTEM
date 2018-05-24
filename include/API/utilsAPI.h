#ifndef UTILS_API_H
#define UTILS_API_H

#include "dlldefines.h"
#include "utils/fileIO.h"

typedef utils::fileIO::ParticleDistribution PD;
typedef utils::fileIO::DistributionFromDisk DFD;

//ParticleDistribution functions
DLLEXP_EXTC PD* createParticleDistributionAPI(const char* saveFolder, const char* attrNames, const char* particleName, double mass);
DLLEXP_EXTC void addPDEnergyRangeAPI(PD* pd, int energyBins, double Emin, double Emax, bool logE = true);
DLLEXP_EXTC void addPDPitchRangeAPI(PD* pd, int pitchBins, double PAmin, double PAmax, bool midBin = true);
DLLEXP_EXTC void generatePDAPI(PD* pd, double s_ion, double s_mag);
DLLEXP_EXTC void padExtraPDAttrsAPI(PD* pd);
DLLEXP_EXTC void writePDAPI(PD* pd);

//DistributionFromDisk functions
DLLEXP_EXTC DFD* loadDistributionFromDiskAPI(const char* name, const char* loadFolder, const char* attrNames, const char* particleName);
DLLEXP_EXTC const double* DistFromDiskDataAPI(DFD* dfd, int attrInd);
DLLEXP_EXTC void DistFromDiskPrintAPI(DFD* dfd, int at);
DLLEXP_EXTC void DistFromDiskPrintDiffAPI(DFD* dfd_this, DFD* dfd_other, int at);
DLLEXP_EXTC void DistFromDiskZeroesAPI(DFD* dfd);
DLLEXP_EXTC void DistFromDiskCompareAPI(DFD* dfd_this, DFD* dfd_other);
DLLEXP_EXTC void deleteDistFromDiskAPI(DFD* dfd);

#endif /* UTILS_API_H */