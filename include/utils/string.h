#ifndef UTILS_STRING_H
#define UTILS_STRING_H

#include <string>
#include <vector>
#include "dllexport.h"

namespace utils
{
	namespace string
	{
		DLLEXP_NOEXTC int findAttrInd(std::string attr, std::vector<std::string> allAttrs);
		DLLEXP_NOEXTC std::vector<std::string> charToStrVec(std::string str, const char delim = ',');
		DLLEXP_NOEXTC std::vector<double> charToDblVec(std::string str, const char delim = ',');
		DLLEXP_NOEXTC void stringPadder(std::string& in, int totalStrLen, int indEraseFrom = 0);
	}
}

#endif /* !UTILS_STRING_H */