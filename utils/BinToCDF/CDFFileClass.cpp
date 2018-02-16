#include "CDFFileClass.h"

void CDFFileClass::writeNewZVar(std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD)
{//can take up to two dimensions
	if (dimSizes.size() > 2)
		throw std::invalid_argument("writezVar1Rec: function is not able to handle more than 2 dimensions at this time");

	long zVarIndex{ static_cast<long>(zVarIndicies.size()) }; //size before pushback is equivalent to size - 1 after pushback, since vector is zero indexed, we need -1
	zVarNameStrs.push_back(varName);
	zVarIndicies.push_back(0);

	std::vector<int> intervals(dimSizes.size(), 1);
	std::vector<int> indicies(dimSizes.size(), 0);
	int dimVar[]{ VARY, VARY };

	//CDFstatus CDFcreatezVar(CDFid id, char* varName, long dataType, long numElements, long numDims, long dimSizes[], long recVariance, long dimVariances[], long* varNum)
	//CDFstatus CDFhyperPutzVarData(CDFid id, long varNum, long recStart, long recCount, long recInterval, long indicies[], long counts[], long intervals[], void* buffer)
	exitStatus_m = CDFcreatezVar(cdfid_m, varName.c_str(), cdftype, 1, (long)dimSizes.size(), dimSizes.data(), VARY, dimVar, &zVarIndicies.at(zVarIndex)); CDFCHKERR();
	exitStatus_m = CDFhyperPutzVarData(cdfid_m, zVarIndicies.at(zVarIndex), 0, 1, 1, indicies.data(), dimSizes.data(), intervals.data(), arrayXD); CDFCHKERR();
}

void CDFFileClass::createZVarAttr(long varNum, std::string attrName, long cdftype, void* data)
{
	long attrIndex{ static_cast<long>(attrIndicies.size()) };
	attrNameStrs.push_back(attrName);
	attrIndicies.push_back(0);

	exitStatus_m = CDFcreateAttr(cdfid_m, attrName.c_str(), VARIABLE_SCOPE, &attrIndicies.at(attrIndex)); CDFCHKERR();
	exitStatus_m = CDFputAttrzEntry(cdfid_m, attrIndicies.at(attrIndex), varNum, cdftype, 1, data); CDFCHKERR();
}

void CDFFileClass::createZVarAttr(std::string varName, std::string attrName, long cdftype, void* data)
{
	long varNum{ findzVarCDFIndexByName(varName) };
	createZVarAttr(varNum, attrName, cdftype, data);
}

long CDFFileClass::findzVarCDFIndexByName(std::string varName)
{
	for (int zVar = 0; zVar < zVarNameStrs.size(); zVar++)
	{
		if (varName == zVarNameStrs.at(zVar))
			return zVarIndicies.at(zVar);
	}

	throw std::invalid_argument ("CDFFileClass::findAttrCDFIndexByName: invalid argument - cannot find zVariable of name " + varName);
}