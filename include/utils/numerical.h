#ifndef UTILS_NUMERICAL_H
#define UTILS_NUMERICAL_H

#include <string>
#include <vector>
#include <iostream>
#include "dllexport.h"

namespace utils
{
	namespace numerical
	{
		DLLEXP_NOEXTC void normalize(std::vector<double>& normalizeMe, double normFactor, bool inverse=false);
		DLLEXP_NOEXTC double calcMean(const std::vector<double>& calcMyMean, bool absValue);
		DLLEXP_NOEXTC double calcStdDev(std::vector<double>& calcMyStdDev);
		DLLEXP_NOEXTC void coutMinMaxErr(const std::vector<double>& basevals, const std::vector<double>& testvals, std::string label="", bool skipzeroes=true);
		DLLEXP_NOEXTC void coutNumAboveErrEps(const std::vector<double>& basevals, const std::vector<double>& testvals, double errEps, std::string label="", bool skipzeroes=true);
	}
}

#endif /* !UTILS_NUMERICAL_H */