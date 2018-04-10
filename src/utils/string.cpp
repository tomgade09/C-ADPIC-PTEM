#include "utils\string.h"

namespace utils
{
	namespace string
	{
		DLLEXP_NOEXTC std::string discoverBFieldType(std::string attrDir)
		{
			std::string ret;

			for (int iii = 0; iii < BELEMNAMES.size(); iii++)
			{
				try
				{
					std::string tmp;
					fileIO::readTxtFile(tmp, "BField_" + BELEMNAMES.at(iii) + ".bin"); //will throw exception if file doesn't exist and next line won't execute
					ret = BELEMNAMES.at(iii); //file exists - therefore this model is used
					break; //only one model possible per sim, so we can break here
				}
				catch (std::invalid_argument&)
				{//indicates file doesn't exist - therefore this isn't the B Field model used
					continue;
				}
				catch (...)
				{//unknown error - pass up to caller
					throw;
				}
			}

			if (ret.empty())
				throw std::runtime_error("utils::string::discoverBFieldType: No BField attribute files are found in the specified folder: " + attrDir);

			return ret;
		}

		DLLEXP_NOEXTC std::vector<std::string> discoverEFieldTypes(std::string attrDir)
		{
			std::vector<std::string> ret;

			for (int iii = 0; iii < EELEMNAMES.size(); iii++)
			{
				try
				{
					std::string tmp;
					fileIO::readTxtFile(tmp, "EField_" + EELEMNAMES.at(iii) + ".bin"); //will throw exception if file doesn't exist and next line won't execute
					ret.push_back(EELEMNAMES.at(iii)); //file exists - therefore this model is used
				}
				catch (std::invalid_argument&)
				{//indicates file doesn't exist - therefore this isn't the B Field model used
					continue;
				}
				catch (...)
				{//unknown error - pass up to caller
					throw;
				}
			}

			return ret;
		}

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

		DLLEXP_NOEXTC std::vector<std::string> charToStrVec(const char* str, const char delim)
		{
			std::string tmp{ str };
			std::vector<std::string> charVec;

			if (tmp == "")
				return charVec;

			size_t loc{ 0 };
			while (loc != std::string::npos)
			{
				loc = tmp.find(delim);
				charVec.push_back(tmp.substr(0, loc));
				tmp.erase(0, loc + 1);
				while (tmp.at(0) == ' ')
					tmp.erase(0, 1);
			}

			return charVec;
		}

		DLLEXP_NOEXTC std::vector<double> charToDblVec(const char* str, const char delim)
		{
			std::vector<std::string> strVec{ charToStrVec(str, delim) };
			std::vector<double> ret;

			if (strVec.size() == 0)
				return ret;

			for (int str = 0; str < strVec.size(); str++)
				ret.push_back(atof(strVec.at(str).c_str()));

			return ret;
		}

		DLLEXP_NOEXTC int sizeofStrVecFromFile(std::string fileName)
		{
			std::string str;
			fileIO::readTxtFile(str, fileName);
			return (int)std::vector<std::string>{charToStrVec(str.c_str())}.size();
		}

		DLLEXP_NOEXTC void stringPadder(std::string& in, int totalStrLen, int indEraseFrom)
		{
			if (totalStrLen <= 0 || indEraseFrom <= 0)
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