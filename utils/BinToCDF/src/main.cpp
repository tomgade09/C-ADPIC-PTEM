#include "FileIO\fileIO.h"

#define WIN32
#include "cdf.h"

#define CDFCHKERR() { if (exitStatus != CDF_OK) { CDFgetStatusText(exitStatus, errTxt); std::cout << "Something went wrong: " << errTxt << " : " << __FILE__ << ":" << __LINE__ << std::endl; } }

CDFid setupCDF(char* fileName, char* dataName, long energyBins, long angleBins)
{
	std::cout << "In setupCDF" << std::endl;

	CDFid     cdfid;
	CDFstatus exitStatus;
	char      errTxt[CDF_STATUSTEXT_LEN + 1];

	/* Create File */
	exitStatus = CDFcreateCDF(fileName, &cdfid);
	if (exitStatus != CDF_OK) { CDFgetStatusText(exitStatus, errTxt); std::cout << "Something went wrong with createCDF: " << errTxt << std::endl; }

	/* Write global attributes */
	//data
	//energy
	//theta
	//geom
	//denergy
	//dtheta

	long nameAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "DATA_NAME", GLOBAL_SCOPE, &nameAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nameAttrNum, 0, CDF_CHAR, strlen(dataName), dataName); CDFCHKERR();

	long validAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "VALID", GLOBAL_SCOPE, &validAttrNum); CDFCHKERR();
	int valid{ 1 };
	exitStatus = CDFputAttrgEntry(cdfid, validAttrNum, 0, CDF_INT4, 1, &valid); CDFCHKERR();

	long projNameNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "PROJECT_NAME", GLOBAL_SCOPE, &projNameNum); CDFCHKERR();
	char* proj{ "FAST SIM" };
	exitStatus = CDFputAttrgEntry(cdfid, projNameNum, 0, CDF_CHAR, strlen(proj), proj); CDFCHKERR();

	long unitsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "UNITS_NAME", GLOBAL_SCOPE, &unitsAttrNum); CDFCHKERR();
	char* units{ "Counts" };
	exitStatus = CDFputAttrgEntry(cdfid, unitsAttrNum, 0, CDF_CHAR, strlen(units), units); CDFCHKERR();

	long timeAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "TIME", GLOBAL_SCOPE, &timeAttrNum); CDFCHKERR();
	double t0{ 0.0 };
	exitStatus = CDFputAttrgEntry(cdfid, timeAttrNum, 0, CDF_DOUBLE, 1, &t0); CDFCHKERR();
	
	long endTimeAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "END_TIME", GLOBAL_SCOPE, &endTimeAttrNum); CDFCHKERR();
	double tfinal{ 250.0 };
	exitStatus = CDFputAttrgEntry(cdfid, endTimeAttrNum, 0, CDF_DOUBLE, 1, &tfinal); CDFCHKERR();

	long nEBinsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "NEBINS", GLOBAL_SCOPE, &nEBinsAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nEBinsAttrNum, 0, CDF_INT4, 1, &energyBins); CDFCHKERR();
	
	long nAngBinsAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "NANGBINS", GLOBAL_SCOPE, &nAngBinsAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrgEntry(cdfid, nAngBinsAttrNum, 0, CDF_INT4, 1, &angleBins); CDFCHKERR();

	long massAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "MASS", GLOBAL_SCOPE, &massAttrNum); CDFCHKERR();
	double mass{ 9.10938356e-31 };
	exitStatus = CDFputAttrgEntry(cdfid, massAttrNum, 0, CDF_DOUBLE, 1, &mass); CDFCHKERR();

	std::cout << nameAttrNum << validAttrNum << projNameNum << unitsAttrNum << timeAttrNum << endTimeAttrNum << nEBinsAttrNum << nAngBinsAttrNum << massAttrNum << std::endl;

	return cdfid;
}

int main()
{
	constexpr int NEBINS{ 96 };
	constexpr int NANGLEBINS{ 1800 };
	constexpr int PARTICLECOUNT{ NEBINS * NANGLEBINS };

	std::string savefolder{ "./../../../../../_dataout/180209_17.45.00/" };

	std::vector<double> btmElec_time (PARTICLECOUNT);
	std::vector<double> btmElec_vpara(PARTICLECOUNT);
	std::vector<double> btmElec_vperp(PARTICLECOUNT);

	std::vector<double> topElec_time (PARTICLECOUNT);
	std::vector<double> topElec_vpara(PARTICLECOUNT);
	std::vector<double> topElec_vperp(PARTICLECOUNT);

	fileIO::readDblBin(btmElec_time,  savefolder + "bins/satellites/btmElec_time.bin",  PARTICLECOUNT);
	fileIO::readDblBin(btmElec_vpara, savefolder + "bins/satellites/btmElec_vpara.bin", PARTICLECOUNT);
	fileIO::readDblBin(btmElec_vperp, savefolder + "bins/satellites/btmElec_vperp.bin", PARTICLECOUNT);

	fileIO::readDblBin(topElec_time,  savefolder + "bins/satellites/topElec_time.bin",  PARTICLECOUNT);
	fileIO::readDblBin(topElec_vpara, savefolder + "bins/satellites/topElec_vpara.bin", PARTICLECOUNT);
	fileIO::readDblBin(topElec_vperp, savefolder + "bins/satellites/topElec_vperp.bin", PARTICLECOUNT);

	CDFid     cdfid;
	CDFstatus exitStatus;
	char      errTxt[CDF_STATUSTEXT_LEN + 1];

	cdfid = setupCDF("mirrorFelectrons", "Binned E, Pitch - Mirror F only", 96, 1800);
	
	//process data from files
	//get energies and pitch, separate times, consolidate times
	//count in bins, mult by maxwellian factor

	//write data here - maybe make zVar or do in setupCDF

	CDFcloseCDF(cdfid);

	return 0;
}