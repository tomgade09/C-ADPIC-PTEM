#include "utils\string.h"

namespace utils
{
	namespace string
	{
		DLLEXP int findAttrInd(std::string attr, std::vector<std::string> allAttrs)
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

		DLLEXP std::vector<std::string> strToStrVec(std::string str, const char delim)
		{
			std::vector<std::string> strVec;

			if (str == "")
				return strVec;

			size_t loc{ 0 };
			while (loc != std::string::npos)
			{
				loc = str.find(delim);
				strVec.push_back(str.substr(0, loc));
				str.erase(0, loc + 1);
				while (str.at(0) == ' ')
					str.erase(0, 1);
			}

			return strVec;
		}

		DLLEXP std::string strVecToStr(std::vector<std::string> strVec, const char delim)
		{
			std::string ret;
			for (auto& str : strVec)
				ret += str + ((str != strVec.back()) ? std::string{ delim } : "");
			
			return ret;
		}

		DLLEXP std::vector<double> strToDblVec(std::string str, const char delim)
		{
			std::vector<std::string> strVec{ strToStrVec(str, delim) };
			std::vector<double> ret;

			if (strVec.size() == 0)
				return ret;

			for (int str = 0; str < strVec.size(); str++)
				ret.push_back(atof(strVec.at(str).c_str()));

			return ret;
		}

		DLLEXP void stringPadder(std::string& in, int totalStrLen, int indEraseFrom)
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