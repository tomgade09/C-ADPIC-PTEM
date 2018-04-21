#ifndef UTILS_NUMERICAL_H
#define UTILS_NUMERICAL_H

#include <string>
#include <vector>
#include <iostream>
#include "dllexport.h"
#include "physicalconstants.h"

namespace utils
{
	namespace numerical
	{
		DLLEXP_NOEXTC void v2DtoEPitch(const std::vector<double>& vpara, const std::vector<double>& vperp, double mass, std::vector<double>& energies_eV, std::vector<double>& pitches_deg);
		DLLEXP_NOEXTC void EPitchTov2D(const std::vector<double>& energies_eV, const std::vector<double>& pitches_deg, double mass, std::vector<double>& vpara, std::vector<double>& vperp);
		DLLEXP_NOEXTC std::vector<double> generateSpacedValues(double min, double max, int number, bool logSpaced, bool endInclusive);
		DLLEXP_NOEXTC void normalize(std::vector<double>& normalizeMe, double normFactor, bool inverse = false);
		DLLEXP_NOEXTC double calcMean(const std::vector<double>& calcMyMean, bool absValue = false);
		DLLEXP_NOEXTC double calcStdDev(const std::vector<double>& calcMyStdDev);
		DLLEXP_NOEXTC void coutMinMaxErr(const std::vector<double>& basevals, const std::vector<double>& testvals, std::string label="", bool skipzeroes = true);
		DLLEXP_NOEXTC void coutNumAboveErrEps(const std::vector<double>& basevals, const std::vector<double>& testvals, double errEps, std::string label="", bool skipzeroes = true);
	}
}

#endif /* !UTILS_NUMERICAL_H */