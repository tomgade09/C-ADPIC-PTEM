#include "utils\string.h"

namespace utils
{
	namespace string
	{
		DLLEXP_NOEXTC int findAttrInd(std::string attr, std::vector<std::string> allAttrs)
		{
			for (int ind = 0; ind < allAttrs.size(); ind++)
			{
				if (allAttrs.at(ind) == attr)
					return ind;
			}

			std::string allAttrsStr;
			for (int attr = 0; attr < allAttrs.size(); attr++)
				allAttrsStr += allAttrs.at(attr);

			throw std::invalid_argument("utils::string::findAttrInd: cannot find attribute " + attr + " in string " + allAttrsStr);
		}

		DLLEXP_NOEXTC std::vector<std::string> charToStrVec(std::string str, const char delim)
		{
			std::vector<std::string> charVec;

			if (str == "")
				return charVec;

			size_t loc{ 0 };
			while (loc != std::string::npos)
			{
				loc = str.find(delim);
				charVec.push_back(str.substr(0, loc));
				str.erase(0, loc + 1);
				while (str.at(0) == ' ')
					str.erase(0, 1);
			}

			return charVec;
		}

		DLLEXP_NOEXTC std::vector<double> charToDblVec(std::string str, const char delim)
		{
			std::vector<std::string> strVec{ charToStrVec(str, delim) };
			std::vector<double> ret;

			if (strVec.size() == 0)
				return ret;

			for (int str = 0; str < strVec.size(); str++)
				ret.push_back(atof(strVec.at(str).c_str()));

			return ret;
		}

		DLLEXP_NOEXTC void stringPadder(std::string& in, int totalStrLen, int indEraseFrom)
		{
			if (totalStrLen <= 0 || indEraseFrom < 0)
				return;

			size_t txtlen = in.length();

			if ((totalStrLen - txtlen) > 0)
			{
				for (int iii = 0; iii < (totalStrLen - txtlen); iii++)
					in += ' ';
			}
			else
				in.erase(indEraseFrom, txtlen - totalStrLen);
		}
	}
}