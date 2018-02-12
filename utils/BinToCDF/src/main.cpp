#include "FileIO\fileIO.h"

#define WIN32
#include "cdf.h"

CDFstatus exitStatus;
char errTxt[CDF_STATUSTEXT_LEN + 1];
#define CDFCHKERR() { if (exitStatus != CDF_OK) { CDFgetStatusText(exitStatus, errTxt); std::cout << "Something went wrong: " << errTxt << " : " << __FILE__ << ":" << __LINE__ << std::endl; } }

CDFid setupCDF(char* fileName, char* dataName, int energyBins, int angleBins, double eBinSize, double aBinSize, std::vector<int>& indsOut)
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

	int dimSizes[]{ 18, 48 };
	int dimVar[]{ VARY, VARY };
	//exitStatus = CDFcreatezVar(cdfid, "Ionospheric-Source Electrons Escape Energy/Pitch Angle Count", CDF_INT4, 1, 2, dimSizes, VARY, dimVar, ionIndOut); CDFCHKERR();
	//exitStatus = CDFcreatezVar(cdfid, "Magnetospheric-Source Electrons Escape Energy/Pitch Angle Count", CDF_INT4, 1, 2, dimSizes, VARY, dimVar, magIndOut); CDFCHKERR();
	
	//CDFstatus CDFcreatezVar(CDFid id, char* varName, int dataType, int numElements, int numDims, int dimSizes[], int recVariance, int dimVariances[], int* varNum)
	exitStatus = CDFcreatezVar(cdfid, "Electrons Escape Energy/Pitch Angle Count", CDF_INT4, 1, 2, dimSizes, VARY, dimVar, &indsOut.at(0)); CDFCHKERR();
	exitStatus = CDFcreatezVar(cdfid, "Mid-Bin Energies", CDF_DOUBLE, 1, 1, &dimSizes[1], VARY, &dimVar[0], &indsOut.at(1)); CDFCHKERR();
	exitStatus = CDFcreatezVar(cdfid, "Mid-Bin Angles", CDF_DOUBLE, 1, 1, &dimSizes[0], VARY, &dimVar[1], &indsOut.at(2)); CDFCHKERR();

	int eBinzAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "LOG E BIN SIZE (LOG EV)", VARIABLE_SCOPE, &eBinzAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrzEntry(cdfid, eBinzAttrNum, indsOut.at(1), CDF_DOUBLE, 1, &eBinSize); CDFCHKERR();

	int aBinzAttrNum{ 0 };
	exitStatus = CDFcreateAttr(cdfid, "ANG BIN SIZE (LOG EV)", VARIABLE_SCOPE, &aBinzAttrNum); CDFCHKERR();
	exitStatus = CDFputAttrzEntry(cdfid, aBinzAttrNum, indsOut.at(2), CDF_DOUBLE, 1, &aBinSize); CDFCHKERR();

	return cdfid;
}

void writezVar1Rec(CDFid id, int variableInd, std::vector<int> dimSizes, void* arrayXD)
{
	std::vector<int> intervals(dimSizes.size(), 1);
	std::vector<int> indicies(dimSizes.size(), 0);

	//CDFstatus CDFhyperPutzVarData(CDFid id, long varNum, long recStart, long recCount, long recInterval, long indicies[], long counts[], long intervals[], void* buffer)
	exitStatus = CDFhyperPutzVarData(id, variableInd, 0, 1, 1, indicies.data(), dimSizes.data(), intervals.data(), arrayXD); CDFCHKERR();
}

std::vector<std::vector<int>> countBins(std::vector<double> energies, std::vector<double> pitches, int particlecount, int nebins, int nanglebins, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::cout << sizeof(int) << "  " << sizeof(int) << std::endl << std::endl;

	std::vector<int>                ecount(nebins);
	std::vector<std::vector<int>>   eAngCount(nanglebins);   //vec[angbin][ebin]
	std::vector<double>             energiesLogMidBin(nebins);
	std::vector<double>             anglesMidBin(nanglebins);

	for (int angbin = 0; angbin < nanglebins; angbin++)
		eAngCount.at(angbin) = ecount;

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / nanglebins};

	double logeMn{ 0.5 };
	double logeMx{ 4.5 };
	double eBinSz{ (logeMx - logeMn) / nebins };

	std::cout << aBinSz << "  " << eBinSz << std::endl;

	int partCountTotal{ 0 };

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

	std::cout << partCountTotal << "  total" << std::endl;
	std::cout << "Data count out:" << std::endl;
	for (int angbin = 0; angbin < nanglebins; angbin++)
	{
		for (int ebin = 0; ebin < nebins; ebin++)
		{
			std::cout << eAngCount.at(angbin).at(ebin) << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl << "Angle Bins:" << std::endl;
	for (int angbin = 0; angbin < nanglebins; angbin++)
		std::cout << anglesMidBin.at(angbin) << "  ";
	std::cout << std::endl << std::endl << "Energy Bins:" << std::endl;
	for (int ebin = 0; ebin < nebins; ebin++)
		std::cout << std::setprecision(6) << energiesLogMidBin.at(ebin) << "  ";
	std::cout << std::endl;

	return eAngCount;
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
	std::vector<std::vector<int>> ang_eCounts{ countBins(energies, pitches, PARTICLECOUNT, NEBINS, NANGLEBINS, binEnergies, binAngles) };

	//count in bins, mult by maxwellian factor
	
	int array2D[NANGLEBINS][NEBINS];
	for (int ang = 0; ang < NANGLEBINS; ang++)
		for (int eng = 0; eng < NEBINS; eng++)
			array2D[ang][eng] = ang_eCounts.at(ang).at(eng);


	/* Create CDF file and setup with appropriate variables, write data */
	CDFid     cdfid;
	std::vector<int> zVarIndicies(10); //[all particles, mid-bin energies, mid-bin angles]

	double angleMin{ 0.0 };
	double angleMax{ 180.0 };
	double aBinSz{ (angleMax - angleMin) / NANGLEBINS };

	double logeMn{ 0.5 };
	double logeMx{ 4.5 };
	double eBinSz{ (logeMx - logeMn) / NEBINS };

	cdfid = setupCDF("electrons", "Electrons Escape Energy/Pitch Angle Count", NEBINS, NANGLEBINS, eBinSz, aBinSz, zVarIndicies);
	writezVar1Rec(cdfid, zVarIndicies.at(0), { NANGLEBINS, NEBINS }, array2D); //Counts data
	writezVar1Rec(cdfid, zVarIndicies.at(1), { NEBINS }, binEnergies.data());  //mid-bin energies
	writezVar1Rec(cdfid, zVarIndicies.at(2), { NANGLEBINS }, binAngles.data());//mid-bin angles

	CDFcloseCDF(cdfid);

	return 0;
}