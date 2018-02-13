#include "FileIO\fileIO.h"

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
		throw std::invalid_argument ("writezVar1Rec: function is not able to handle more than 2 dimensions at this time");

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

std::vector<std::vector<long long>> countBins(std::vector<double> energies, std::vector<double> pitches, int particlecount, int nebins, int nanglebins, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::vector<long long>                ecount(nebins);
	std::vector<std::vector<long long>>   eAngCount(nanglebins);   //vec[angbin][ebin]
	std::vector<double>             energiesLogMidBin(nebins);
	std::vector<double>             anglesMidBin(nanglebins);

	for (int angbin = 0; angbin < nanglebins; angbin++)
		eAngCount.at(angbin) = ecount;

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / nanglebins};

	double logeMn{ 0.49 };
	double logeMx{ 4.5 };
	double eBinSz{ (logeMx - logeMn) / nebins };

	long long partCountTotal{ 0 };

	for (int part = 0; part < particlecount; part++)
	{
		double emsmt{ energies.at(part) };
		double angmsmt{ pitches.at(part) };

		bool found{ false };

		for (int angbin = 0; angbin < nanglebins; angbin++)
		{
			double angmin{ angbin * aBinSz };
			double angmax{ (angbin + 1) * aBinSz };

			anglesMidBin.at(angbin) = angmin + 0.5 * aBinSz; //writes angle in degrees

			for (int ebin = 0; ebin < nebins; ebin++)
			{
				double emin{ /*((ebin==0) ? (-100.0) : */(ebin * eBinSz + logeMn)/*)*/ };
				double emax{ /*((ebin==NEBINS-1) ? (100.0) : */((ebin + 1) * eBinSz + logeMn)/*)*/ };

				if (angbin == 0) { energiesLogMidBin.at(ebin) = pow(10, emin + 0.5 * eBinSz); } //writes energies in eV
				if (angbin == nanglebins - 1) { angmax += 0.001; } //I want the top angle bin to include 180 and right now the conditional is ... < 180, not ...<= 180
				if (angmsmt >= angmin && angmsmt < angmax && emsmt >= emin && emsmt < emax)
				{
					partCountTotal++;
					eAngCount.at(angbin).at(ebin) += 1;
					found = true;
				}
			}
		}
		//if (!found) { std::cout << "Not found: Ind: " << part << ": E: " << energies.at(part) << ": Pitch: " << pitches.at(part) << std::endl; }
	}

	binEnergies = energiesLogMidBin;
	binAngles = anglesMidBin;

	return eAngCount;
}

std::vector<std::vector<long long>> weightMaxwellian(std::vector<std::vector<long long>> counts, std::vector<double> eBinMids, double logeMin, double dloge, double angleMin, double dangle, double sigmaT_eV)
{
	double c{ 1 / (sqrt(2 * 3.14159265358979323846) * sigmaT_eV) };
	double eBinSize{ pow(10, logeMin + dloge * counts.at(0).size()) - pow(10, logeMin + dloge * (counts.at(0).size() - 1)) };
	double expNormFactor{ eBinSize * c * exp(- pow(10, logeMin + dloge * (counts.at(0).size() - 0.5)) / sigmaT_eV) }; //logeMin + dloge * ebinsTot gets you to the top - loge = 4.5; we want the last bin which is a half step back from the end point, hence -0.5
	
	for (int angs = 0; angs < counts.size(); angs++)
	{
		for (int engs = 0; engs < counts.at(0).size(); engs++)
		{
			double eBinMin{ pow(10, logeMin + dloge * engs) };
			double eBinMax{ pow(10, logeMin + dloge * (engs + 1)) };
			long long multFact{ (long long)round((eBinMax - eBinMin) * c * exp(-eBinMids.at(engs) / sigmaT_eV) / expNormFactor) }; //also have to mult by bin width
			counts.at(angs).at(engs) *= multFact;
			//if (angs == 0) { std::cout << multFact << "  "; }
			//if (angs == counts.size() - 1) { std::cout << counts.at(angs).at(engs) << "  "; }
		}
		//if (angs == 0 || angs == counts.size() - 1) { std::cout << std::endl << std::endl; }
	}

	std::cout << std::endl << "Maxwellian count out:" << std::endl;
	for (int angbin = 0; angbin < counts.size(); angbin++)
	{
		for (int ebin = 0; ebin < counts.at(0).size(); ebin++)
		{
			std::cout << counts.at(angbin).at(ebin) << " ";
		}
		std::cout << std::endl;
	}

	std::vector<std::vector<long long>> ret{ counts };
	return ret;
}

int main()
{
	/* Load data from completed sim */
	constexpr int NEBINS{ 48 };
	constexpr int NANGLEBINS{ 18 };
	constexpr int PARTICLECOUNT{ 96 * 1800 };

	std::string savefolder{ "./../../../../../_dataout/180212_12.18.59/" };

	std::vector<double> btmElec_time(PARTICLECOUNT);
	std::vector<double> btmElec_vpara(PARTICLECOUNT);
	std::vector<double> btmElec_vperp(PARTICLECOUNT);

	std::vector<double> topElec_time(PARTICLECOUNT);
	std::vector<double> topElec_vpara(PARTICLECOUNT);
	std::vector<double> topElec_vperp(PARTICLECOUNT);

	std::vector<double> final_vpara(PARTICLECOUNT);
	std::vector<double> final_vperp(PARTICLECOUNT);

	try {
		fileIO::readDblBin(btmElec_time,  savefolder + "bins/satellites/btmElec_time.bin",  PARTICLECOUNT);
		fileIO::readDblBin(btmElec_vpara, savefolder + "bins/satellites/btmElec_vpara.bin", PARTICLECOUNT);
		fileIO::readDblBin(btmElec_vperp, savefolder + "bins/satellites/btmElec_vperp.bin", PARTICLECOUNT);

		fileIO::readDblBin(topElec_time,  savefolder + "bins/satellites/topElec_time.bin",  PARTICLECOUNT);
		fileIO::readDblBin(topElec_vpara, savefolder + "bins/satellites/topElec_vpara.bin", PARTICLECOUNT);
		fileIO::readDblBin(topElec_vperp, savefolder + "bins/satellites/topElec_vperp.bin", PARTICLECOUNT);

		fileIO::readDblBin(final_vpara, savefolder + "bins/particles_final/elec_vpara.bin", PARTICLECOUNT);
		fileIO::readDblBin(final_vperp, savefolder + "bins/particles_final/elec_vperp.bin", PARTICLECOUNT);
	}
	catch (std::exception& exp)
	{
		std::cout << exp.what() << std::endl;
		exit(1);
	}


	/* Process data */
	std::vector<double> energies(PARTICLECOUNT); //log(eV)
	std::vector<double> pitches(PARTICLECOUNT);  //deg

	for (int part = 0; part < PARTICLECOUNT; part++)
	{
		energies.at(part) = log10(0.5 * 9.10938356e-31 * (final_vpara.at(part) * final_vpara.at(part) + final_vperp.at(part) * final_vperp.at(part)) / 1.60218e-19);
		pitches.at(part) = atan2(abs(final_vperp.at(part)), -final_vpara.at(part)) * 180.0 / 3.14159265358979323846;
	}

	std::vector<double> binEnergies;
	std::vector<double> binAngles;
	std::vector<std::vector<long long>> ang_eCounts{ countBins(energies, pitches, PARTICLECOUNT, NEBINS, NANGLEBINS, binEnergies, binAngles) };

	std::cout << std::endl << "Count out:" << std::endl;
	for (int angbin = 0; angbin < ang_eCounts.size(); angbin++)
	{
		for (int ebin = 0; ebin < ang_eCounts.at(0).size(); ebin++)
		{
			std::cout << ang_eCounts.at(angbin).at(ebin) << " ";
		}
		std::cout << std::endl;
	}

	//count in bins, mult by maxwellian factor
	std::vector<std::vector<long long>> maxWeighted{ weightMaxwellian(ang_eCounts, binEnergies, 0.49, (4.5 - 0.49) / 48, 0, 10, 1000.0) };
	
	long long maxArray2D[NANGLEBINS][NEBINS];
	long long cntArray2D[NANGLEBINS][NEBINS];
	for (int ang = 0; ang < NANGLEBINS; ang++)
		for (int eng = 0; eng < NEBINS; eng++)
		{
			cntArray2D[ang][eng] = ang_eCounts.at(ang).at(eng);
			maxArray2D[ang][eng] = maxWeighted.at(ang).at(eng);
		}

	/* Create CDF file and setup with appropriate variables, write data */
	CDFid     cdfid;
	std::vector<long> zVarInds(10);
	std::vector<long> attrInds(10);

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / NANGLEBINS };

	double logeMn{ 0.49 };
	double logeMx{ 4.5 };
	double eBinSz{ (logeMx - logeMn) / NEBINS };

	cdfid = setupCDF("test", "Electrons Escape Energy/Pitch Angle Count", NEBINS, NANGLEBINS);
	//void createWritezVar1Rec(CDFid id, std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD, int& variableInd)
	//void createzVarAttr(CDFid id, long varNum, std::string attrName, long cdftype, void* data, int& attrNum)
	createWritezVar1Rec(cdfid, "Mid-Bin Energies (eV)",                                          CDF_DOUBLE, { NEBINS }, binEnergies.data(), zVarInds.at(0));
	createWritezVar1Rec(cdfid, "Mid-Bin Angles (Degrees)",                                       CDF_DOUBLE, { NANGLEBINS }, binAngles.data(), zVarInds.at(1));
	createWritezVar1Rec(cdfid, "Electrons Escape Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_INT8, { NANGLEBINS, NEBINS }, maxArray2D, zVarInds.at(2));
	createWritezVar1Rec(cdfid, "Electrons Escape Energy/Pitch Angle Count",                      CDF_INT8, { NANGLEBINS, NEBINS }, cntArray2D, zVarInds.at(3));

	createzVarAttr(cdfid, zVarInds.at(0), "LOG E BIN SIZE (LOG EV)", CDF_DOUBLE, &eBinSz, attrInds.at(0));
	createzVarAttr(cdfid, zVarInds.at(1), "ANGLE BIN SIZE (DEGREES)", CDF_DOUBLE, &aBinSz, attrInds.at(1));

	CDFcloseCDF(cdfid);

	return 0;
}