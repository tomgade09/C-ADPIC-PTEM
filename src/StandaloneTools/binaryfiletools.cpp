#include "StandaloneTools\binaryfiletools.h"

std::vector<std::string> constCharToStrVec(const char* str, const char delim)
{
	std::string tmp{ str };
	std::vector<std::string> charVec;

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

void stringPadder(std::string& in, int totalStrLen, int indEraseFrom)
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

std::vector<std::vector<std::vector<double>>> form3DvectorArray(int numPartTypes, int numAttr, int numParts)
{
	std::vector<std::vector<std::vector<double>>> ret;

	std::vector<double> partVec;
	partVec.resize(numParts);

	for (int type = 0; type < numPartTypes; type++)
	{
		std::vector<std::vector<double>> attrVec;
		for (int attrs = 0; attrs < numAttr; attrs++)
			attrVec.push_back(partVec);
		ret.push_back(attrVec);
	}

	return ret;
}

std::vector<std::vector<double>> form2DvectorArray(int numAttr, int numParts)
{
	std::vector<std::vector<double>> ret;

	std::vector<double> partVec;
	partVec.resize(numParts);

	for (int attrs = 0; attrs < numAttr; attrs++)
		ret.push_back(partVec);

	return ret;
}