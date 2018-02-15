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

std::vector<std::vector<double>> countBins(std::vector<double> energies, std::vector<double> pitches, int particlecount, int nebins, int nanglebins, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::vector<double>                ecount(nebins);
	std::vector<std::vector<double>>   eAngCount(nanglebins);   //vec[angbin][ebin]
	std::vector<double>                energiesLogMidBin(nebins);
	std::vector<double>                anglesMidBin(nanglebins);

	for (int angbin = 0; angbin < nanglebins; angbin++) //zeros vector array
		eAngCount.at(angbin) = ecount;

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / nanglebins };

	double logeMidBinMn{ 0.5 };
	double logeMidBinMx{ 4.5 };
	double eBinSz{ (logeMidBinMx - logeMidBinMn) / (nebins - 1) };
	
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
				double emin{ /*((ebin==0) ? (-100.0) : */((ebin - 0.5) * eBinSz + logeMidBinMn)/*)*/ };
				double emax{ /*((ebin==NEBINS-1) ? (100.0) : */((ebin + 0.5) * eBinSz + logeMidBinMn)/*)*/ };

				if (part == 0 && angbin == 0) { energiesLogMidBin.at(ebin) = pow(10, ebin * eBinSz + logeMidBinMn); /*std::cout << energiesLogMidBin.at(ebin); ((ebin != nebins - 1) ? (std::cout << " ") : (std::cout << std::endl));*/ } //writes energies in eV
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

std::vector<std::vector<double>> weightMaxwellian(std::vector<std::vector<double>> counts, std::vector<double> eBinMids, double dEflux_kT, double sigma_kT, double logeMin, double dloge, double angleMin, double dangle)
{
	double logTeV{ log10(sigma_kT) };
	double binWidth{ pow(10, logTeV + dloge * 0.5) - pow(10, logTeV - dloge * 0.5) };
	double multFactor{ dEflux_kT / (exp(-1.0) * binWidth) };
	
	for (int angs = 0; angs < counts.size(); angs++)
	{
		for (int engs = 0; engs < counts.at(0).size(); engs++)
		{
			double eBinWid{ pow(10, logeMin + dloge * (engs + 0.5)) - pow(10, logeMin + dloge * (engs - 0.5)) };
			double multFact{ eBinWid * exp(-eBinMids.at(engs) / sigma_kT) * multFactor };
			counts.at(angs).at(engs) *= multFact;
			//if (angs == 0) { std::cout << multFact << "  "; }
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

	std::vector<std::vector<double>> ret{ counts };
	return ret;
}

void loadVparaPerp(std::string savefolder, std::string filePrefix, int count, std::vector<double>& vpara, std::vector<double>& vperp)
{
	try
	{
		fileIO::readDblBin(vpara, savefolder + filePrefix + "vpara.bin", count);
		fileIO::readDblBin(vperp, savefolder + filePrefix + "vperp.bin", count);
	}
	catch (std::exception& exp)
	{
		std::cout << exp.what() << std::endl;
		exit(1);
	}
}

int main()
{
	/* Load data from completed sim */
	constexpr int NEBINS{ 48 };
	constexpr int NANGLEBINS{ 18 };
	constexpr int PARTICLECOUNT{ 96 * 1800 };

	std::string savefolder{ "./../../../../../_dataout/180213_15.58.50/" };

	std::vector<double> satdown_vpara(PARTICLECOUNT);
	std::vector<double> satdown_vperp(PARTICLECOUNT);
	std::vector<double> satup_vpara(PARTICLECOUNT);
	std::vector<double> satup_vperp(PARTICLECOUNT);

	loadVparaPerp(savefolder, "bins/satellites/4e6ElecDown_", PARTICLECOUNT, satdown_vpara, satdown_vperp);
	loadVparaPerp(savefolder, "bins/satellites/4e6ElecUp_", PARTICLECOUNT, satup_vpara, satup_vperp);

	std::vector<double> final_vpara_ion(PARTICLECOUNT/2);
	std::vector<double> final_vperp_ion(PARTICLECOUNT/2);
	std::vector<double> final_vpara_mag(PARTICLECOUNT/2);
	std::vector<double> final_vperp_mag(PARTICLECOUNT/2);

	std::vector<int> doubles_ion;
	std::vector<int> zeroes_ion;
	std::vector<int> doubles_mag;
	std::vector<int> zeroes_mag;

	for (int part = 0; part < PARTICLECOUNT; part++)
	{
		if (satdown_vpara.at(part) != 0.0)
		{
			(part < PARTICLECOUNT / 2) ? (final_vpara_ion.at(part) = satdown_vpara.at(part)) : (final_vpara_mag.at(part - PARTICLECOUNT / 2) = satdown_vpara.at(part - PARTICLECOUNT / 2));
			(part < PARTICLECOUNT / 2) ? (final_vperp_ion.at(part) = satdown_vperp.at(part)) : (final_vperp_mag.at(part - PARTICLECOUNT / 2) = satdown_vperp.at(part - PARTICLECOUNT / 2));
		}
		else if (satup_vpara.at(part) != 0.0)
		{
			(part < PARTICLECOUNT / 2) ? (final_vpara_ion.at(part) = satup_vpara.at(part)) : (final_vpara_mag.at(part - PARTICLECOUNT / 2) = satup_vpara.at(part - PARTICLECOUNT / 2));
			(part < PARTICLECOUNT / 2) ? (final_vperp_ion.at(part) = satup_vperp.at(part)) : (final_vperp_mag.at(part - PARTICLECOUNT / 2) = satup_vperp.at(part - PARTICLECOUNT / 2));
		}
		else
			//std::cout << "Both vpara elements 0 at index: " << part << std::endl;
			(part < PARTICLECOUNT / 2) ? (zeroes_ion.push_back(part)) : (zeroes_mag.push_back(part - PARTICLECOUNT / 2));
		
		if ((satdown_vpara.at(part) != 0.0) && (satup_vpara.at(part) != 0.0))
			(part < PARTICLECOUNT / 2) ? (doubles_ion.push_back(part)) : (doubles_mag.push_back(part - PARTICLECOUNT / 2));
	}
	std::cout << "final values: " << final_vpara_ion.at(0) << " " << final_vpara_ion.at(100) << " " << final_vpara_mag.at(0) << final_vpara_mag.at(100) << std::endl;
	std::cout << "zeroes: " << zeroes_ion.size() << "  " << zeroes_mag.size() << "doubles: " << doubles_ion.size() << "  " << doubles_mag.size() << std::endl;
	std::cout << "zeroes indicies: " << zeroes_mag.at(0) << "  " << zeroes_mag.at(1) << std::endl;
	if (doubles_ion.size() > zeroes_ion.size() || doubles_mag.size() > zeroes_mag.size())
		std::cout << "less zero values than doubles" << std::endl;
	for (int dbl = 0; dbl < doubles_ion.size(); dbl++)
	{
		final_vpara_ion.at(zeroes_ion.at(dbl)) = satup_vpara.at(doubles_ion.at(dbl));
		final_vperp_ion.at(zeroes_ion.at(dbl)) = satup_vperp.at(doubles_ion.at(dbl));
	}
	for (int dbl = 0; dbl < doubles_mag.size(); dbl++)
	{
		final_vpara_mag.at(zeroes_mag.at(dbl)) = satup_vpara.at(doubles_mag.at(dbl));
		final_vperp_mag.at(zeroes_mag.at(dbl)) = satup_vperp.at(doubles_mag.at(dbl));
	}

	std::cout << "biggest not sensed: " << zeroes_mag.at(zeroes_mag.size() - 1) << "  " << zeroes_mag.at(zeroes_mag.size() - 2) << std::endl << std::endl;


	/* Process data */
	std::vector<double> energies_ion(PARTICLECOUNT/2); //log(eV)
	std::vector<double> pitches_ion(PARTICLECOUNT/2);  //deg
	std::vector<double> energies_mag(PARTICLECOUNT/2);
	std::vector<double> pitches_mag(PARTICLECOUNT/2);

	//separate magsph and ionsph here

	for (int part = 0; part < PARTICLECOUNT; part++)
	{
		if (part < PARTICLECOUNT / 2)
		{
			if (final_vpara_ion.at(part) != 0.0 && final_vperp_ion.at(part) != 0.0)
			{
				energies_ion.at(part) = log10(0.5 * 9.10938356e-31 * (final_vpara_ion.at(part) * final_vpara_ion.at(part) + final_vperp_ion.at(part) * final_vperp_ion.at(part)) / 1.60218e-19);
				pitches_ion.at(part) = atan2(abs(final_vperp_ion.at(part)), -final_vpara_ion.at(part)) * 180.0 / 3.14159265358979323846;
			}
		}
		else
		{
			if (final_vpara_mag.at(part - PARTICLECOUNT / 2) != 0.0 && final_vperp_mag.at(part - PARTICLECOUNT / 2) != 0.0)
			{
				energies_mag.at(part - PARTICLECOUNT/2) = log10(0.5 * 9.10938356e-31 * (final_vpara_mag.at(part - PARTICLECOUNT/2) * final_vpara_mag.at(part - PARTICLECOUNT/2) + final_vperp_mag.at(part - PARTICLECOUNT/2) * final_vperp_mag.at(part - PARTICLECOUNT/2)) / 1.60218e-19);
				pitches_mag.at(part - PARTICLECOUNT/2) = atan2(abs(final_vperp_mag.at(part - PARTICLECOUNT/2)), -final_vpara_mag.at(part - PARTICLECOUNT/2)) * 180.0 / 3.14159265358979323846;
			}
		}
	}

	std::vector<double> binEnergies;
	std::vector<double> binAngles;
	std::vector<std::vector<double>> ang_eCounts_ion{ countBins(energies_ion, pitches_ion, PARTICLECOUNT/2, NEBINS, NANGLEBINS, binEnergies, binAngles) };
	std::vector<std::vector<double>> ang_eCounts_mag{ countBins(energies_mag, pitches_mag, PARTICLECOUNT/2, NEBINS, NANGLEBINS, binEnergies, binAngles) };

	std::cout << std::endl << "Count out ion:" << std::endl;
	for (int angbin = 0; angbin < ang_eCounts_ion.size(); angbin++)
	{
		for (int ebin = 0; ebin < ang_eCounts_ion.at(0).size(); ebin++)
		{
			std::cout << ang_eCounts_ion.at(angbin).at(ebin) << " ";
		}
		std::cout << std::endl;
	}

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / NANGLEBINS };

	double logeMidBinMn{ 0.5 };
	double logeMidBinMx{ 4.5 };
	double eBinSz{ (logeMidBinMx - logeMidBinMn) / (NEBINS - 1) };

	//count in bins, mult by maxwellian factor
	std::vector<std::vector<double>> maxWeightedIon{ weightMaxwellian(ang_eCounts_ion, binEnergies, 1.0e7, 1.0, logeMidBinMn, eBinSz, angleMin, aBinSz) };
	std::vector<std::vector<double>> maxWeightedMag{ weightMaxwellian(ang_eCounts_mag, binEnergies, 1.0e8, 1000.0, logeMidBinMn, eBinSz, angleMin, aBinSz) };
	
	std::vector<std::vector<double>> ang_eCounts;
	std::vector<std::vector<double>> maxWeighted;

	for (int iii = 0; iii < maxWeightedIon.size(); iii++)
	{
		std::vector<double> angTmp;
		std::vector<double> maxTmp;
		for (int jjj = 0; jjj < maxWeightedIon.at(0).size(); jjj++)
		{
			angTmp.push_back(ang_eCounts_ion.at(iii).at(jjj) + ang_eCounts_mag.at(iii).at(jjj));
			maxTmp.push_back(maxWeightedIon.at(iii).at(jjj) + maxWeightedMag.at(iii).at(jjj));
		}
		ang_eCounts.push_back(angTmp);
		maxWeighted.push_back(maxTmp);
	}

	double maxArray2D[NANGLEBINS][NEBINS];
	double cntArray2D[NANGLEBINS][NEBINS];
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

	cdfid = setupCDF("4e6Altitude", "Electrons Energy/Pitch Angle at 4000km Altitude", NEBINS, NANGLEBINS);
	//void createWritezVar1Rec(CDFid id, std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD, int& variableInd)
	//void createzVarAttr(CDFid id, long varNum, std::string attrName, long cdftype, void* data, int& attrNum)
	createWritezVar1Rec(cdfid, "Mid-Bin Energies (eV)",                                          CDF_DOUBLE, { NEBINS }, binEnergies.data(), zVarInds.at(0));
	createWritezVar1Rec(cdfid, "Mid-Bin Angles (Degrees)",                                       CDF_DOUBLE, { NANGLEBINS }, binAngles.data(), zVarInds.at(1));
	createWritezVar1Rec(cdfid, "Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { NANGLEBINS, NEBINS }, maxArray2D, zVarInds.at(2));
	createWritezVar1Rec(cdfid, "Electrons Energy/Pitch Angle Count",                      CDF_DOUBLE, { NANGLEBINS, NEBINS }, cntArray2D, zVarInds.at(3));

	createzVarAttr(cdfid, zVarInds.at(0), "LOG E BIN SIZE (LOG EV)", CDF_DOUBLE, &eBinSz, attrInds.at(0));
	createzVarAttr(cdfid, zVarInds.at(1), "ANGLE BIN SIZE (DEGREES)", CDF_DOUBLE, &aBinSz, attrInds.at(1));

	CDFcloseCDF(cdfid);

	return 0;
}