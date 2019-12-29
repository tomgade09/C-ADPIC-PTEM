#ifndef UTILS_STRING_H
#define UTILS_STRING_H

#include <string>
#include <vector>
#include "dlldefines.h"

namespace utils
{
	namespace strings
	{
		DLLEXP size_t findAttrInd(std::string attr, std::vector<std::string> allAttrs);
		DLLEXP std::vector<std::string> strToStrVec(std::string str, const char delim = ',');
		DLLEXP std::string strVecToStr(std::vector<std::string> strVec, const char delim = ',');
		DLLEXP std::vector<double> strToDblVec(std::string str, const char delim = ',');
		DLLEXP void stringPadder(std::string& in, size_t totalStrLen, int indEraseFrom = 0);
	}
}

#endif /* !UTILS_STRING_H */