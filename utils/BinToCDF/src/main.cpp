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

std::vector<std::vector<double>> countBins(std::vector<double> energies, std::vector<double> pitches, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::vector<double>                ecount(binEnergies.size());
	std::vector<std::vector<double>>   eAngCount(binAngles.size());   //vec[angbin][ebin]

	for (int angbin = 0; angbin < eAngCount.size(); angbin++) //zeros vector array
		eAngCount.at(angbin) = ecount;

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / binAngles.size() }; //sets bin min to 0 and bin max to 180 with midbins at 5 + 10*ind

	double logeMidBinMin{ 0.5 };
	double logeMidBinMax{ 4.5 };
	double eBinSz{ (logeMidBinMax - logeMidBinMin) / (binEnergies.size() - 1) }; //sets midbins from 0.5 - 4.5 log eV

	for (int part = 0; part < energies.size(); part++)
	{
		for (int angbin = 0; angbin < binAngles.size(); angbin++)
		{
			double angmin{ aBinSz * angbin + angleMin };
			double angmax{ aBinSz * (angbin + 1) + angleMin };

			if (part == 0) { binAngles.at(angbin) = (angbin + 0.5) * aBinSz + angleMin; } //writes angle in degrees

			for (int ebin = 0; ebin < binEnergies.size(); ebin++)
			{
				double emin{ ((ebin - 0.5) * eBinSz + logeMidBinMin) };
				double emax{ ((ebin + 0.5) * eBinSz + logeMidBinMin) };

				if (part == 0 && angbin == 0) { binEnergies.at(ebin) = pow(10, ebin * eBinSz + logeMidBinMin); /*std::cout << energiesLogMidBin.at(ebin); ((ebin != nebins - 1) ? (std::cout << " ") : (std::cout << std::endl));*/ } //writes energies in eV
				if (angbin == binAngles.size() - 1) { angmax += 0.001; } //I want the top angle bin to include 180 and right now the conditional is ... < 180, not ...<= 180
				if (pitches.at(part) >= angmin && pitches.at(part) < angmax && energies.at(part) >= emin && energies.at(part) < emax)
					eAngCount.at(angbin).at(ebin) += 1;
			}
		}
	}

	std::cout << std::endl << "Count out:" << std::endl;
	for (int angbin = 0; angbin < eAngCount.size(); angbin++)
	{
		for (int ebin = 0; ebin < eAngCount.at(0).size(); ebin++)
		{
			std::cout << eAngCount.at(angbin).at(ebin) << " ";
		}
		std::cout << std::endl;
	}

	return eAngCount;
}

std::vector<std::vector<double>> weightMaxwellian(std::vector<std::vector<double>> counts, std::vector<double> binEnergies, double dEflux_kT, double sigma_kT)
{
	double logeMin{ log10(binEnergies.at(0)) };
	double dloge{ log10(binEnergies.at(1)) - logeMin };

	double logTeV{ log10(sigma_kT) };
	double binWidth{ pow(10, logTeV + dloge * 0.5) - pow(10, logTeV - dloge * 0.5) };
	double multFactor{ dEflux_kT / (exp(-1.0) * binWidth) };
	
	for (int angs = 0; angs < counts.size(); angs++)
	{
		for (int engs = 0; engs < counts.at(0).size(); engs++)
		{
			double eBinWid{ pow(10, dloge * (engs + 0.5) + logeMin) - pow(10, dloge * (engs - 0.5) + logeMin) };
			double multFact{ eBinWid * exp(-binEnergies.at(engs) / sigma_kT) * multFactor };
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

void processUpDown(std::vector<double> down_vpara, std::vector<double> down_vperp, std::vector<double> up_vpara, std::vector<double> up_vperp, std::vector<double>& out_vpara, std::vector<double>& out_vperp)
{
	int doubles{ 0 };
	int zeroes{ 0 };
	std::vector<int> zeroesInds;
	std::vector<int> doublesInds;

	for (int part = 0; part < down_vpara.size(); part++)
	{
		bool downNonZero{ down_vpara.at(part) != 0.0 || down_vperp.at(part) != 0.0 };
		bool upNonZero{ up_vpara.at(part) != 0.0 || up_vperp.at(part) != 0.0 };

		if (downNonZero)
		{
			out_vpara.at(part) = down_vpara.at(part);
			out_vperp.at(part) = down_vperp.at(part);
		}
		else if (upNonZero)
		{
			out_vpara.at(part) = up_vpara.at(part);
			out_vperp.at(part) = up_vperp.at(part);
		}
		else
		{
			zeroes++;
			zeroesInds.push_back(part);
		}

		if (downNonZero && upNonZero)
		{
			//out_vpara.push_back(up_vpara.at(part));
			//out_vperp.push_back(up_vperp.at(part));
			doublesInds.push_back(part);
			doubles++;
		}
	}

	for (int iii = 0; iii < doublesInds.size(); iii++)
	{
		out_vpara.at(zeroesInds.at(iii)) = up_vpara.at(doublesInds.at(iii));
		out_vperp.at(zeroesInds.at(iii)) = up_vperp.at(zeroesInds.at(iii));
	}

	std::cout << "Summary: doubles: " << doubles << " zeroes: " << zeroes << " size of vpara/perp: " << out_vpara.size() << std::endl;
}

void getParticlePitchLogE(std::vector<double> vpara, std::vector<double> vperp, double mass, std::vector<double>& particleEnergies, std::vector<double>& particleAngles)
{
	for (int part = 0; part < vpara.size(); part++)
	{
		bool nonZero{ vpara.at(part) != 0.0 || vperp.at(part) != 0.0 };

		if (nonZero)
		{
			particleEnergies.at(part) = log10(0.5 * mass * (vpara.at(part) * vpara.at(part) + vperp.at(part) * vperp.at(part)) / 1.60218e-19);
			particleAngles.at(part) = atan2(abs(vperp.at(part)), -vpara.at(part)) * 180.0 / 3.14159265358979323846;
		}

	}
}

void splitVectorInHalf(std::vector<double> in, std::vector<double>& out1, std::vector<double>& out2)
{
	for (int iii = 0; iii < in.size(); iii++)
		(iii < in.size() / 2) ? (out1.push_back(in.at(iii))) : (out2.push_back(in.at(iii)));
}

int main()
{
	/* Load data from completed sim */
	constexpr int NEBINS{ 48 };
	constexpr int NANGLEBINS{ 18 };
	constexpr int PARTICLECOUNT{ 96 * 1800 };

	std::string savefolder{ "./../../../../../_dataout/180214_15.37.09/" };

	std::vector<double> satdown_vpara(PARTICLECOUNT);
	std::vector<double> satdown_vperp(PARTICLECOUNT);
	std::vector<double> satup_vpara(PARTICLECOUNT);
	std::vector<double> satup_vperp(PARTICLECOUNT);

	//loadVparaPerp(savefolder, "bins/particles_init/elec_", PARTICLECOUNT, satdown_vpara, satdown_vperp);
	loadVparaPerp(savefolder, "bins/satellites/4e6ElecDown_", PARTICLECOUNT, satdown_vpara, satdown_vperp);
	loadVparaPerp(savefolder, "bins/satellites/4e6ElecUp_", PARTICLECOUNT, satup_vpara, satup_vperp);

	int downnegzerocount{ 0 };
	int upppnegzerocount{ 0 };
	for (int iii = 0; iii < PARTICLECOUNT; iii++)
	{
		if (satdown_vperp.at(iii) == -0.0)
		{
			downnegzerocount++;
			satdown_vperp.at(iii) = 1.0;
			satdown_vperp.at(iii) = 0.0;
		}
		if (satup_vperp.at(iii) = -0.0)
		{
			upppnegzerocount++;
			satup_vperp.at(iii) = 1.0;
			satup_vperp.at(iii) = 0.0;
		}
	}
	std::cout << "negative zeros: " << downnegzerocount << "  " << upppnegzerocount << std::endl;

	std::vector<double> satdown_vpara_ion;
	std::vector<double> satdown_vperp_ion;
	std::vector<double> satdown_vpara_mag;
	std::vector<double> satdown_vperp_mag;
	std::vector<double> satup_vpara_ion;
	std::vector<double> satup_vperp_ion;
	std::vector<double> satup_vpara_mag;
	std::vector<double> satup_vperp_mag;

	splitVectorInHalf(satdown_vpara, satdown_vpara_ion, satdown_vpara_mag);
	splitVectorInHalf(satdown_vperp, satdown_vperp_ion, satdown_vperp_mag);
	splitVectorInHalf(satup_vpara, satup_vpara_ion, satup_vpara_mag);
	splitVectorInHalf(satup_vperp, satup_vperp_ion, satup_vperp_mag);

	for (int iii = 0; iii < PARTICLECOUNT/2; iii++)
	{
		if (satdown_vpara.at(iii) != satdown_vpara_ion.at(iii))
			std::cout << "Error down vpara - ind: " << iii << std::endl;
		if (satdown_vpara.at(iii + PARTICLECOUNT / 2) != satdown_vpara_mag.at(iii))
			std::cout << "Error down vpara - ind: " << iii + PARTICLECOUNT / 2 << std::endl;
		if (satdown_vperp.at(iii) != satdown_vperp_ion.at(iii))
			std::cout << "Error down vperp - ind: " << iii << std::endl;
		if (satdown_vperp.at(iii + PARTICLECOUNT / 2) != satdown_vperp_mag.at(iii))
			std::cout << "Error down vperp - ind: " << iii + PARTICLECOUNT / 2 << std::endl;

		if (satup_vpara.at(iii) != satup_vpara_ion.at(iii))
			std::cout << "Error up vpara - ind: " << iii << std::endl;
		if (satup_vpara.at(iii + PARTICLECOUNT / 2) != satup_vpara_mag.at(iii))
			std::cout << "Error up vpara - ind: " << iii + PARTICLECOUNT / 2 << std::endl;
		if (satup_vperp.at(iii) != satup_vperp_ion.at(iii))
			std::cout << "Error up vperp - ind: " << iii << std::endl;
		if (satup_vperp.at(iii + PARTICLECOUNT / 2) != satup_vperp_mag.at(iii))
			std::cout << "Error up vperp - ind: " << iii + PARTICLECOUNT / 2 << std::endl;
	}

	//std::cout << "Split check:" << std::endl;
	//std::cout << "down_para: " << satdown_vpara.at(0) << " " << satdown_vpara_ion.at(0) << " " << satdown_vpara.at(PARTICLECOUNT / 2 + 1) << " " << satdown_vpara_mag.at(1) << std::endl;
	//std::cout << "down_perp: " << satdown_vperp.at(0) << " " << satdown_vperp_ion.at(0) << " " << satdown_vperp.at(PARTICLECOUNT / 2 + 1) << " " << satdown_vperp_mag.at(1) << std::endl;
	//std::cout << "up_para: " << satup_vpara.at(PARTICLECOUNT - 1) << " " << satup_vpara_mag.at(PARTICLECOUNT / 2 - 1) << std::endl;
	//std::cout << "up_perp: " << satup_vperp.at(PARTICLECOUNT - 200) << " " << satup_vperp_mag.at(PARTICLECOUNT / 2 - 200) << std::endl;

	std::vector<double> final_vpara_ion(PARTICLECOUNT/2);
	std::vector<double> final_vperp_ion(PARTICLECOUNT/2);
	std::vector<double> final_vpara_mag(PARTICLECOUNT/2);
	std::vector<double> final_vperp_mag(PARTICLECOUNT/2);
	processUpDown(satdown_vpara_ion, satdown_vperp_ion, satup_vpara_ion, satup_vperp_ion, final_vpara_ion, final_vperp_ion);
	processUpDown(satdown_vpara_mag, satdown_vperp_mag, satup_vpara_mag, satup_vperp_mag, final_vpara_mag, final_vperp_mag);

	for (int iii = 0; iii < PARTICLECOUNT / 2; iii++)
	{
		if (final_vpara_ion.at(iii) != satdown_vpara_ion.at(iii) && final_vpara_ion.at(iii) != satup_vpara_ion.at(iii))
			std::cout << "Mismatch final_vpara_ion: " << final_vpara_ion.at(iii) << " " << satdown_vpara_ion.at(iii) << " " << satup_vpara_ion.at(iii) << std::endl;
		if (final_vperp_ion.at(iii) != satdown_vperp_ion.at(iii) && final_vperp_ion.at(iii) != satup_vperp_ion.at(iii))
			std::cout << "Mismatch final_vperp_ion: " << final_vperp_ion.at(iii) << " " << satdown_vperp_ion.at(iii) << " " << satup_vperp_ion.at(iii) << std::endl;

		if (final_vpara_mag.at(iii) != satdown_vpara_mag.at(iii) && final_vpara_mag.at(iii) != satup_vpara_mag.at(iii) && (satdown_vpara_mag.at(iii) != 0.0 || satup_vpara_mag.at(iii) != 0.0))
			std::cout << "Mismatch final_vpara_mag: " << final_vpara_mag.at(iii) << " " << satdown_vpara_mag.at(iii) << " " << satup_vpara_mag.at(iii) << std::endl;
			//count++;
		if (final_vperp_mag.at(iii) != satdown_vperp_mag.at(iii) && final_vperp_mag.at(iii) != satup_vperp_mag.at(iii) && (satdown_vperp_mag.at(iii) != 0.0 || satup_vperp_mag.at(iii) != 0.0))
			std::cout << "Mismatch final_vperp_mag: " << final_vperp_mag.at(iii) << " " << satdown_vperp_mag.at(iii) << " " << satup_vperp_mag.at(iii) << std::endl;
			//count++;
	}
	
	for (int iii = 0; iii < 10 * 98; iii += 98)
	{
		std::cout << "ion index: " << iii << " vpara: " << final_vpara_ion.at(iii) << " vperp: " << final_vperp_ion.at(iii) << std::endl;
		std::cout << "mag index: " << iii << " vpara: " << final_vperp_ion.at(iii) << " vperp: " << final_vperp_mag.at(iii) << std::endl;
		printf("%.4e\n", final_vperp_mag.at(iii));
	}

	std::cout << std::endl;

	std::vector<double> pitches_ion(final_vpara_ion.size());
	std::vector<double> energies_ion(final_vpara_ion.size());
	getParticlePitchLogE(final_vpara_ion, final_vperp_ion, 9.10938356e-31, energies_ion, pitches_ion);

	std::vector<double> pitches_mag(final_vpara_mag.size());
	std::vector<double> energies_mag(final_vpara_mag.size());
	getParticlePitchLogE(final_vpara_mag, final_vperp_mag, 9.10938356e-31, energies_mag, pitches_mag);

	std::cout << std::endl;

	for (int iii = 0; iii < 10 * 98; iii += 98)
	{
		std::cout << "ion index: " << iii << " pitch: " << pitches_ion.at(iii) << " energy: " << energies_ion.at(iii) << std::endl;
		std::cout << "mag index: " << iii << " pitch: " << pitches_mag.at(iii) << " energy: " << energies_mag.at(iii) << std::endl;
	}

	std::cout << "Doneski" << std::endl;
	//return 0;

	std::vector<double> binEnergies(NEBINS);
	std::vector<double> binAngles(NANGLEBINS);
	std::vector<std::vector<double>> ang_eCounts_ion{ countBins(energies_ion, pitches_ion, binEnergies, binAngles) };
	std::vector<std::vector<double>> ang_eCounts_mag{ countBins(energies_mag, pitches_mag, binEnergies, binAngles) };

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / NANGLEBINS };

	double logeMidBinMin{ 0.5 };
	double logeMidBinMax{ 4.5 };
	double eBinSz{ (logeMidBinMax - logeMidBinMin) / (NEBINS - 1) };

	//count in bins, mult by maxwellian factor
	std::vector<std::vector<double>> maxWeightedIon{ weightMaxwellian(ang_eCounts_ion, binEnergies, 1.0e7, 1.0) };
	std::vector<std::vector<double>> maxWeightedMag{ weightMaxwellian(ang_eCounts_mag, binEnergies, 1.0e8, 1000.0) };
	
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

	cdfid = setupCDF("4e6Altitude", "Electrons Energy/Pitch Angle", NEBINS, NANGLEBINS);
	//void createWritezVar1Rec(CDFid id, std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD, int& variableInd)
	//void createzVarAttr(CDFid id, long varNum, std::string attrName, long cdftype, void* data, int& attrNum)
	createWritezVar1Rec(cdfid, "Mid-Bin Energies (eV)",                                   CDF_DOUBLE, { NEBINS }, binEnergies.data(), zVarInds.at(0));
	createWritezVar1Rec(cdfid, "Mid-Bin Angles (Degrees)",                                CDF_DOUBLE, { NANGLEBINS }, binAngles.data(), zVarInds.at(1));
	createWritezVar1Rec(cdfid, "Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { NANGLEBINS, NEBINS }, maxArray2D, zVarInds.at(2));
	createWritezVar1Rec(cdfid, "Electrons Energy/Pitch Angle Count",                      CDF_DOUBLE, { NANGLEBINS, NEBINS }, cntArray2D, zVarInds.at(3));

	createzVarAttr(cdfid, zVarInds.at(0), "LOG E BIN SIZE (LOG EV)", CDF_DOUBLE, &eBinSz, attrInds.at(0));
	createzVarAttr(cdfid, zVarInds.at(1), "ANGLE BIN SIZE (DEGREES)", CDF_DOUBLE, &aBinSz, attrInds.at(1));

	CDFcloseCDF(cdfid);

	return 0;
}