#ifndef UTILS_NUMERICAL_H
#define UTILS_NUMERICAL_H

#include <string>
#include <vector>
#include <iostream>
#include "dlldefines.h"
#include "physicalconstants.h"

namespace utils
{
	namespace numerical
	{
		DLLEXP void v2DtoEPitch(const std::vector<double>& vpara, const std::vector<double>& vperp, double mass, std::vector<double>& energies_eV, std::vector<double>& pitches_deg);
		DLLEXP void EPitchTov2D(const std::vector<double>& energies_eV, const std::vector<double>& pitches_deg, double mass, std::vector<double>& vpara, std::vector<double>& vperp);
		DLLEXP std::vector<double> generateSpacedValues(double min, double max, int number, bool logSpaced, bool endInclusive);
		DLLEXP void normalize(std::vector<double>& normalizeMe, double normFactor, bool inverse = false);
		DLLEXP double calcMean(const std::vector<double>& calcMyMean, bool absValue = false);
		DLLEXP double calcStdDev(const std::vector<double>& calcMyStdDev);
		DLLEXP void coutMinMaxErr(const std::vector<double>& basevals, const std::vector<double>& testvals, std::string label="", bool skipzeroes = true);
		DLLEXP void coutNumAboveErrEps(const std::vector<double>& basevals, const std::vector<double>& testvals, double errEps, std::string label="", bool skipzeroes = true);
	}
}

#endif /* !UTILS_NUMERICAL_H */