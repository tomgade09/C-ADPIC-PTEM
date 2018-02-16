#include <memory>

#include "FileIO\fileIO.h"
#include "CDFFileClass.h"

std::vector<std::vector<double>> countBins(std::vector<long> indicies, std::vector<double> logenergies, std::vector<double> pitches, std::vector<double> ionweights, std::vector<double> magweights, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::vector<double>                ecount(binEnergies.size());
	std::vector<std::vector<double>>   eAngCount(binAngles.size());   //vec[angbin][ebin]

	for (int angbin = 0; angbin < eAngCount.size(); angbin++) //zeros vector array
		eAngCount.at(angbin) = ecount;

	double aBinSize{ 180.0 / 18.0 }; //sets bin min to 0 and bin max to 180 with midbins at 5 + 10*ind
	double logeBinSize{ (4.5 - 0.5) / 47 }; //sets midbins from 0.5 - 4.5 log eV

	std::cout << "logenergies sample: " << logenergies.at(0) << " " << logenergies.at(200) << " pitch sample: " << pitches.at(0) << " " << pitches.at(200) << std::endl;
	std::cout << "logenergies size: " << logenergies.size() << std::endl;

	for (int part = 0; part < logenergies.size(); part++)
	{
		bool found{ false };
		for (int angbin = 0; angbin < binAngles.size(); angbin++)
		{
			double angmin{ angbin * aBinSize };
			double angmax{ (angbin + 1) * aBinSize };

			if (part == 0) { binAngles.at(angbin) = (angbin + 0.5) * aBinSize; } //writes angle in degrees
			//if (part == 0) { std::cout << "angle bins: " << angmin << " " << binAngles.at(angbin) << " " << angmax << std::endl; }

			for (int ebin = 0; ebin < binEnergies.size(); ebin++)
			{
				double logemin{ ((ebin - 0.5) * logeBinSize + 0.5) };
				double logemax{ ((ebin + 0.5) * logeBinSize + 0.5) };

				if (part == 0 && angbin == 0) { binEnergies.at(ebin) = pow(10, ebin * logeBinSize + 0.5); /*std::cout << energiesLogMidBin.at(ebin); ((ebin != nebins - 1) ? (std::cout << " ") : (std::cout << std::endl));*/ } //writes energies in eV
				//if (part == 0 && angbin == 0) { std::cout << "e bins: " << pow(10, logemin) << " " << binEnergies.at(ebin) << " " << pow(10, logemax) << std::endl; }
				if (angbin == binAngles.size() - 1 && ebin == 0) { angmax += 0.001; } //I want the top angle bin to include 180 and right now the conditional is ... < 180, not ...<= 180
				//if (part <= 191 && ebin < 1 && angbin < 1) { std::cout << part % magweights.size() << " " << magweights.at(part % magweights.size()) << std::endl; }
				if ((pitches.at(part) >= angmin) &&	(pitches.at(part) < angmax) &&
					(logenergies.at(part) >= logemin) && (logenergies.at(part) < logemax) && !found)
				{
					(part < logenergies.size() / 2) ?
						(eAngCount.at(angbin).at(ebin) += ionweights.at(indicies.at(part) % ionweights.size())) :
						(eAngCount.at(angbin).at(ebin) += magweights.at(indicies.at(part) % magweights.size()));
					found = true; //this shouldn't be necessary but just in case
					//eAngCount.at(angbin).at(ebin) += 1;
				}
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

std::vector<double> maxwellianWeights(std::vector<double> initLogEnergies, double dEflux_kT, double sigma_kT)
{
	double dloge{ (initLogEnergies.at(1) - initLogEnergies.at(0)) };

	double logTeV{ log10(sigma_kT) };
	double binWidth_kT{ pow(10, logTeV + dloge * 0.5) - pow(10, logTeV - dloge * 0.5) };
	double multFactor{ dEflux_kT / (exp(-1.0) * binWidth_kT) };
	
	std::vector<double> ret;
	for (int engInd = 0; engInd < initLogEnergies.size(); engInd++)
	{
		double eBinWid{ pow(10, dloge * (engInd + 0.5) + 0.5) - pow(10, dloge * (engInd - 0.5) + 0.5) };
		ret.push_back( eBinWid * exp(-pow(10, initLogEnergies.at(engInd)) / sigma_kT) * multFactor );
	}

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

void processUpDown(std::vector<double> down_vpara, std::vector<double> down_vperp, std::vector<double> up_vpara, std::vector<double> up_vperp, std::vector<long>& indicies, std::vector<double>& out_vpara, std::vector<double>& out_vperp)
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
			doubles++;
			doublesInds.push_back(part);
		}
	}

	for (int iii = 0; iii < doublesInds.size(); iii++)
	{
		indicies.at(zeroesInds.at(iii)) = indicies.at(doublesInds.at(iii));
		out_vpara.at(zeroesInds.at(iii)) = up_vpara.at(doublesInds.at(iii));
		out_vperp.at(zeroesInds.at(iii)) = up_vperp.at(doublesInds.at(iii));
	}

	std::cout << "Summary: doubles: " << doubles << " zeroes: " << zeroes << " size of vpara/perp: " << out_vpara.size() << std::endl;
}

void getParticlePitchLogE(std::vector<double> vpara, std::vector<double> vperp, double mass, std::vector<double>& particleLogEnergies, std::vector<double>& particleAngles)
{
	for (int part = 0; part < vpara.size(); part++)
	{
		bool nonZero{ vpara.at(part) != 0.0 || vperp.at(part) != 0.0 };

		if (nonZero)
		{
			particleLogEnergies.at(part) = log10(0.5 * mass * (vpara.at(part) * vpara.at(part) + vperp.at(part) * vperp.at(part)) / 1.60218e-19);
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
	constexpr int CDFNEBINS{ 48 };
	constexpr int CDFNANGLEBINS{ 18 };
	constexpr int SIMNEBINS{ 192 };
	constexpr int SIMNANGLEBINS{ 1800 };
	constexpr int PARTICLECOUNT{ SIMNEBINS * SIMNANGLEBINS };

	std::string savefolder{ "./../../../../../_dataout/180216_13.35.56/" };

	std::vector<double> satdown_vpara(PARTICLECOUNT);
	std::vector<double> satdown_vperp(PARTICLECOUNT);
	std::vector<double> satup_vpara(PARTICLECOUNT);
	std::vector<double> satup_vperp(PARTICLECOUNT);
	std::vector<double> init_vpara(PARTICLECOUNT);
	std::vector<double> init_vperp(PARTICLECOUNT);
	std::vector<long> indicies(PARTICLECOUNT);

	for (int iii = 0; iii < PARTICLECOUNT; iii++)
		indicies.at(iii) = iii;

	//loadVparaPerp(savefolder, "bins/particles_init/elec_", PARTICLECOUNT, satdown_vpara, satdown_vperp);
	loadVparaPerp(savefolder, "bins/satellites/4e6ElecDown_", PARTICLECOUNT, satdown_vpara, satdown_vperp);
	loadVparaPerp(savefolder, "bins/satellites/4e6ElecUp_", PARTICLECOUNT, satup_vpara, satup_vperp);
	loadVparaPerp(savefolder, "bins/particles_init/elec_", PARTICLECOUNT, init_vpara, init_vperp);

	std::vector<double> final_vpara(PARTICLECOUNT);
	std::vector<double> final_vperp(PARTICLECOUNT);
	processUpDown(satdown_vpara, satdown_vperp, satup_vpara, satup_vperp, indicies, final_vpara, final_vperp);

	for (int iii = 0; iii < PARTICLECOUNT / 2; iii++)
	{
		if (final_vpara.at(iii) != satdown_vpara.at(iii) && final_vpara.at(iii) != satup_vpara.at(iii) && (satdown_vpara.at(iii) != 0.0 || satup_vpara.at(iii) != 0.0))
			std::cout << "Mismatch final_vpara: " << final_vpara.at(iii) << " " << satdown_vpara.at(iii) << " " << satup_vpara.at(iii) << std::endl;
		if (final_vperp.at(iii) != satdown_vperp.at(iii) && final_vperp.at(iii) != satup_vperp.at(iii) && (satdown_vperp.at(iii) != 0.0 || satup_vperp.at(iii) != 0.0))
			std::cout << "Mismatch final_vperp: " << final_vperp.at(iii) << " " << satdown_vperp.at(iii) << " " << satup_vperp.at(iii) << std::endl;
	}
	
	for (int iii = 0; iii < 10 * 960; iii += 960)
	{
		std::cout << "index: " << iii << " vpara: " << final_vpara.at(iii) << " vperp: " << final_vperp.at(iii) << std::endl;
	}

	std::cout << std::endl;

	std::vector<double> pitches(final_vpara.size());
	std::vector<double> logenergies(final_vpara.size());
	getParticlePitchLogE(final_vpara, final_vperp, 9.10938356e-31, logenergies, pitches);

	std::vector<double> initpitches(PARTICLECOUNT);
	std::vector<double> initlogenergies(PARTICLECOUNT);
	getParticlePitchLogE(init_vpara, init_vperp, 9.10938356e-31, initlogenergies, initpitches);

	std::cout << "energies done" << std::endl;

	for (int iii = 0; iii < 10 * 96; iii += 96)
	{
		std::cout << "index: " << iii << " pitch: " << pitches.at(iii) << " energy: " << logenergies.at(iii) << std::endl;
	}

	//count in bins, mult by maxwellian factor
	std::vector<double> initEbins;
	for (int iii = 0; iii < 192; iii++)
		initEbins.push_back(initlogenergies.at(iii));

	std::vector<double> maxWeightsIon{ maxwellianWeights(initEbins, 1.0e7, 1.0) };
	std::vector<double> maxWeightsMag{ maxwellianWeights(initEbins, 1.0e8, 1000.0) };

	std::cout << "Sample max weights in main func: " << maxWeightsIon.at(0) << " " << maxWeightsIon.at(10) << " " << maxWeightsMag.at(0) << " " << maxWeightsMag.at(10) << std::endl;
	
	std::vector<double> binEnergies(CDFNEBINS);
	std::vector<double> binAngles(CDFNANGLEBINS);
	std::vector<std::vector<double>> ang_eCounts{ countBins(indicies, logenergies, pitches, maxWeightsIon, maxWeightsMag, binEnergies, binAngles) };

	double cntArray2D[CDFNANGLEBINS][CDFNEBINS];
	for (int ang = 0; ang < CDFNANGLEBINS; ang++)
		for (int eng = 0; eng < CDFNEBINS; eng++)
			cntArray2D[ang][eng] = ang_eCounts.at(ang).at(eng);

	/* Create CDF file and setup with appropriate variables, write data */
	std::unique_ptr<CDFFileClass> cdf = std::make_unique<CDFFileClass>("4e6Altitude");

	cdf->writeNewZVar("Mid-Bin Energies (eV)", CDF_DOUBLE, { CDFNEBINS }, binEnergies.data());
	cdf->writeNewZVar("Mid-Bin Angles (Degrees)", CDF_DOUBLE, { CDFNANGLEBINS }, binAngles.data());
	cdf->writeNewZVar("Electrons Energy/Pitch Angle Count, Maxwellian-Weighted", CDF_DOUBLE, { CDFNANGLEBINS, CDFNEBINS }, cntArray2D);

	return 0;
}