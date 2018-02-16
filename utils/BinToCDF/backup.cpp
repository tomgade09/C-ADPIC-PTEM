#include <iostream>
#include <vector>
#define WIN32
#include "cdf.h"

CDFstatus exitStatus;
char errTxt[CDF_STATUSTEXT_LEN + 1];
#define CDFCHKERR() { if (exitStatus != CDF_OK) { CDFgetStatusText(exitStatus, errTxt); std::cout << "Something went wrong: " << errTxt << " : " << __FILE__ << ":" << __LINE__ << std::endl; } }

CDFid setupCDF(char* fileName, char* dataName, int energyBins, int angleBins)
{
	CDFid     cdfid;

	/* Create File */
	exitStatus = CDFcreateCDF(fileName, &cdfid); CDFCHKERR();

	/* Write global attributes */
	//data
	//energy
	//theta
	//geom
	//denergy
	//dtheta

	int nameAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "DATA_NAME", GLOBAL_SCOPE, &nameAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nameAttrNum, 0, CDF_CHAR, strlen(dataName), dataName); CDFCHKERR();

	int validAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "VALID", GLOBAL_SCOPE, &validAttrNum); CDFCHKERR();
	int valid{ 1 };
	exitStatus = CDFputAttrgEntry(cdfid, validAttrNum, 0, CDF_INT4, 1, &valid); CDFCHKERR();

	int projNameNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "PROJECT_NAME", GLOBAL_SCOPE, &projNameNum); CDFCHKERR();
	char* proj{ "FAST SIM" };
	exitStatus = CDFputAttrgEntry(cdfid, projNameNum, 0, CDF_CHAR, strlen(proj), proj); CDFCHKERR();

	int unitsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "UNITS_NAME", GLOBAL_SCOPE, &unitsAttrNum); CDFCHKERR();
	char* units{ "Counts" };
	exitStatus = CDFputAttrgEntry(cdfid, unitsAttrNum, 0, CDF_CHAR, strlen(units), units); CDFCHKERR();

	int timeAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "TIME", GLOBAL_SCOPE, &timeAttrNum); CDFCHKERR();
	double t0{ 0.0 };
	exitStatus = CDFputAttrgEntry(cdfid, timeAttrNum, 0, CDF_DOUBLE, 1, &t0); CDFCHKERR();

	int endTimeAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "END_TIME", GLOBAL_SCOPE, &endTimeAttrNum); CDFCHKERR();
	double tfinal{ 250.0 };
	exitStatus = CDFputAttrgEntry(cdfid, endTimeAttrNum, 0, CDF_DOUBLE, 1, &tfinal); CDFCHKERR();

	int nEBinsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "NEBINS", GLOBAL_SCOPE, &nEBinsAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nEBinsAttrNum, 0, CDF_INT4, 1, &energyBins); CDFCHKERR();

	int nAngBinsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "NANGBINS", GLOBAL_SCOPE, &nAngBinsAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nAngBinsAttrNum, 0, CDF_INT4, 1, &angleBins); CDFCHKERR();

	int massAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "MASS", GLOBAL_SCOPE, &massAttrNum); CDFCHKERR();
	double mass{ 9.10938356e-31 };
	exitStatus = CDFputAttrgEntry(cdfid, massAttrNum, 0, CDF_DOUBLE, 1, &mass); CDFCHKERR();

	return cdfid;
}

void createWritezVar1Rec(CDFid id, std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD, long& variableInd)
{//can take up to two dimensions
	if (dimSizes.size() > 2)
		throw std::invalid_argument("writezVar1Rec: function is not able to handle more than 2 dimensions at this time");

	std::vector<int> intervals(dimSizes.size(), 1);
	std::vector<int> indicies(dimSizes.size(), 0);
	int dimVar[]{ VARY, VARY };

	//CDFstatus CDFcreatezVar(CDFid id, char* varName, long dataType, long numElements, long numDims, long dimSizes[], long recVariance, long dimVariances[], long* varNum)
	//CDFstatus CDFhyperPutzVarData(CDFid id, long varNum, long recStart, long recCount, long recInterval, long indicies[], long counts[], long intervals[], void* buffer)
	exitStatus = CDFcreatezVar(id, varName.c_str(), cdftype, 1, (long)dimSizes.size(), dimSizes.data(), VARY, dimVar, &variableInd); CDFCHKERR();
	exitStatus = CDFhyperPutzVarData(id, variableInd, 0, 1, 1, indicies.data(), dimSizes.data(), intervals.data(), arrayXD); CDFCHKERR();
}

void createzVarAttr(CDFid id, long varNum, std::string attrName, long cdftype, void* data, long& attrNum)
{
	exitStatus = CDFcreateAttr(id, attrName.c_str(), VARIABLE_SCOPE, &attrNum); CDFCHKERR();
	exitStatus = CDFputAttrzEntry(id, attrNum, varNum, cdftype, 1, data); CDFCHKERR();
}