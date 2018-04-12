#include "CDFFileClass.h"

void CDFFileClass::writeNewZVar(std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD)
{//can take up to two dimensions
	if (dimSizes.size() > 2)
		throw std::invalid_argument("writezVar1Rec: function is not able to handle more than 2 dimensions at this time");

	long zVarIndex{ static_cast<long>(zVarIndicies_m.size()) }; //size before pushback is equivalent to size - 1 after pushback, since vector is zero indexed, we need -1
	zVarNameStrs_m.push_back(varName);
	zVarIndicies_m.push_back(0);

	std::vector<int> intervals(dimSizes.size(), 1);
	std::vector<int> indicies(dimSizes.size(), 0);
	int dimVar[]{ VARY, VARY };

	exitStatus_m = CDFcreatezVar(cdfid_m, varName.c_str(), cdftype, 1, (long)dimSizes.size(), dimSizes.data(), VARY, dimVar, &zVarIndicies_m.at(zVarIndex)); CDFCHKERR();
	exitStatus_m = CDFhyperPutzVarData(cdfid_m, zVarIndicies_m.at(zVarIndex), 0, 1, 1, indicies.data(), dimSizes.data(), intervals.data(), arrayXD); CDFCHKERR();
}

void CDFFileClass::writeNewGlobalAttr(std::string attrName, long cdftype, long datalen, void* data)
{
	long attrIndex{ static_cast<long>(attrIndicies_m.size()) };
	attrNameStrs_m.push_back(attrName);
	attrIndicies_m.push_back(0);

	exitStatus_m = CDFcreateAttr(cdfid_m, attrName.c_str(), GLOBAL_SCOPE, &attrIndicies_m.at(attrIndex)); CDFCHKERR();
	exitStatus_m = CDFputAttrgEntry(cdfid_m, attrIndicies_m.at(attrIndex), 0, cdftype, datalen, data); CDFCHKERR();
}

void CDFFileClass::writeNewZVarAttr(long varNum, std::string attrName, long cdftype, long datalen, void* data)
{
	long attrIndex{ static_cast<long>(attrIndicies_m.size()) };
	attrNameStrs_m.push_back(attrName);
	attrIndicies_m.push_back(0);

	exitStatus_m = CDFcreateAttr(cdfid_m, attrName.c_str(), VARIABLE_SCOPE, &attrIndicies_m.at(attrIndex)); CDFCHKERR();
	exitStatus_m = CDFputAttrzEntry(cdfid_m, attrIndicies_m.at(attrIndex), varNum, cdftype, datalen, data); CDFCHKERR();
}

void CDFFileClass::writeNewZVarAttr(std::string varName, std::string attrName, long cdftype, long datalen, void* data)
{
	long varNum{ findzVarCDFIndexByName(varName) };
	writeNewZVarAttr(varNum, attrName, cdftype, datalen, data);
}

long CDFFileClass::findzVarCDFIndexByName(std::string varName)
{
	for (int zVar = 0; zVar < zVarNameStrs_m.size(); zVar++)
	{
		if (varName == zVarNameStrs_m.at(zVar))
			return zVarIndicies_m.at(zVar);
	}

	throw std::invalid_argument ("CDFFileClass::findAttrCDFIndexByName: invalid argument - cannot find zVariable of name " + varName);
}