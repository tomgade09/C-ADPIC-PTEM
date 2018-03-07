#include <memory>

#include "FileIO\fileIO.h"
#include "CDFFileClass.h"

constexpr int CDFNEBINS{ 48 };
constexpr int CDFNANGLEBINS{ 18 };
constexpr int SIMNEBINS{ 96 };
constexpr int SIMNANGLEBINS{ 3600 };
constexpr int PARTICLECOUNT{ SIMNEBINS * SIMNANGLEBINS };

double integralF_bs(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb)
{
	if (upper > incidentE) { throw std::invalid_argument ("integralF_bs: upper limit of integration is higher than incident energy - function is zero above incident energy"); }
	double integral_sec{ (pow(upper, scnd_fact + 1) - pow(lower, scnd_fact + 1)) * pow(10, scnd_logb)     / (scnd_fact + 1) };
	double integral_prm{ (pow(upper, prim_fact + 1) - pow(lower, prim_fact + 1)) * pow(10, prim_logb + 4) / ((prim_fact + 1) * pow(incidentE, prim_fact + 1)) };
	return integral_sec + integral_prm;
}

std::vector<std::vector<double>> backscatter(std::vector<double> ion_init_energies, std::vector<double> ion_init_pitches, std::vector<double> sat_energies,
	std::vector<double> sat_pitches, std::vector<double> escape_energies, std::vector<double> maxWeights_mag, std::vector<double> binEnergies, std::vector<double> binAngles)
{
	/*
		Based on equations:
		BS_scnd(E) = 10 ^ (scnd_factor * logE + scnd_logb)
			=== d * x ^ (scnd_factor), where d = 10 ^ (scnd_logb)
		BS_prim(E) = 10 ^ (prim_factor * log(E / E_incident) + prim_logb) * (10000 / E_incident)
			=== d * x ^ (prim_factor), where d = 10 ^ (prim_logb + 4) / E_incident ^ (prim_factor + 1)
	*/

	double secondary_factor{ -2.1 }; //obtained these from observing excel plotting, what seemed the closest to 'Evans, 1974' backscatter graphs
	double secondary_logb  {  0.3 };  //f = 10 ^ (scnd_factor * log(E_back) + scnd_logb)
	double primary_factor  {  1.5 };  //f = 10 ^ (prim_factor * log(E_back/E_incident) + prim_logb)
	double primary_logb    { -4.0 };

	double dloge{ log10(binEnergies.at(1)) - log10(binEnergies.at(0)) };
	double dangle{ binAngles.at(1) - binAngles.at(0) };

	/* Count number of escaped particles in bin */
	std::vector<double> counts(binEnergies.size());
	for (int part = 0; part < escape_energies.size(); part++)
	{
		for (int ener = 0; ener < binEnergies.size(); ener++)
		{
			double engmin{ pow(10, log10(binEnergies.at(ener)) - 0.5 * dloge) };
			double engmax{ pow(10, log10(binEnergies.at(ener)) + 0.5 * dloge) };

			if (escape_energies.at(part) >= engmin && escape_energies.at(part) < engmax)
				counts.at(ener) += maxWeights_mag.at(part % maxWeights_mag.size());
		}
	}
	
	std::cout << std::endl << "BS counts out:" << std::endl;
	for (int angbin = 0; angbin < counts.size(); angbin++)
		std::cout << counts.at(angbin) << " ";
	std::cout << std::endl;

	/* Calculate BS E non-normalized probability/amplitude, divide isotropically across pitch angles (divide by number of pitch angles) */
	std::vector<double> bsEflux(counts.size());
	for (int ampbin = 0; ampbin < bsEflux.size(); ampbin++)
	{
		double engmin{ pow(10, log10(binEnergies.at(ampbin)) - 0.5 * dloge) };
		double engmax{ pow(10, log10(binEnergies.at(ampbin)) + 0.5 * dloge) };
		
		for (int cntbin = ampbin; cntbin < bsEflux.size(); cntbin++)
		{
			double incidentE{ pow(10, log10(binEnergies.at(cntbin)) + 0.5 * dloge) }; //incident E is upper limit of bin - particles in here can be up to this energy, so why not?
			double intF{ integralF_bs(engmin, engmax, incidentE, primary_factor, primary_logb, secondary_factor, secondary_logb) };
			//double integralF_bs(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb)
			if (ampbin == 0 && cntbin == 0) { std::cout << "integralF_bs args: " << engmin << ", " << engmax << ", " << incidentE << ", " << primary_factor << ", " << primary_logb << ", " << secondary_factor << ", " << secondary_logb << std::endl;
				std::cout << "integralF_bs: " << intF << std::endl; }
			bsEflux.at(ampbin) += intF * counts.at(cntbin) / (engmax - engmin);
		}

		bsEflux.at(ampbin) /= ((double)binAngles.size() / 2); //divide by number of pitch angle bins
	}

	std::cout << std::endl << "bsEflux out:" << std::endl;
	for (int iii = 0; iii < bsEflux.size(); iii++)
		std::cout << bsEflux.at(iii) << "   ";
	std::cout << std::endl;

	/* Convert ion source particle e/pitch to sat e/pitch */
	int ebinfact{ SIMNEBINS / CDFNEBINS };
	int angbinfact{ SIMNANGLEBINS / CDFNANGLEBINS };
	std::vector<std::vector<double>> ret(binAngles.size());
	for (int iii = 0; iii < ret.size(); iii++)
		ret.at(iii) = std::vector<double>(binEnergies.size());
	
	for (int retAbin = 0; retAbin < ret.size(); retAbin++)
	{
		for (int retEbin = 0; retEbin < ret.at(0).size(); retEbin++)
		{
			int fiveOffset{ (int)(5.0 / (180.0/((double)SIMNANGLEBINS))) * SIMNEBINS };
			int binPartDataInd{ retAbin * angbinfact * SIMNEBINS + retEbin * ebinfact + fiveOffset }; //where (index number) is the bin particle in the simulation input data? PARTICLES START AT 180 DEGREES for now
			if (retEbin >= ret.at(0).size() / 2) { binPartDataInd++; }
			double binPartEng{ sat_energies.at(binPartDataInd) }; //gets the above particle at satellite
			double binPartAng{ sat_pitches.at(binPartDataInd) };

			/* 
				Particles above are indexed like so:

				Angles ------->>>>>
				     5...     155       165       175
				-----------------------------------------   Energies
				|   ...   |    7    |    4    |    1    | 3.1  |
				|         |         |         |         |      |
				-----------------------------------------      |
				|         |    8    |    5    |    2    | 3.4 \ /
				|         |         |         |         | ...  `
				-----------------------------------------
				|         |    9    |    6    |    3    | 13k
				|         |         |         |         |
				-----------------------------------------
			*/
			//std::cout << "Found particle: " << retAbin << "," << retEbin << ": init: " << ion_init_energies.at(binPartDataInd) << ", " << ion_init_pitches.at(binPartDataInd) << " bin: " << binEnergies.at(retEbin) << ", " << binAngles.at(ret.size() - 1 - retAbin) << " at sat: " << binPartEng << ", " << binPartAng << std::endl;
			//binAngles starts at 5 degrees (mag source)... sim particles data array starts at 180.0 - PITA to deal with
			for (int angleInd = 0; angleInd < binAngles.size(); angleInd++) //iterate through angles and energies to find out what bin the particle ends up in at satellite
			{
				double angmin{ binAngles.at(angleInd) - 0.5 * dangle };
				double angmax{ binAngles.at(angleInd) + 0.5 * dangle };

				for (int engInd = 0; engInd < binEnergies.size(); engInd++)
				{
					double engmin{ pow(10, log10(binEnergies.at(engInd)) - 0.5 * dloge) };
					double engmax{ pow(10, log10(binEnergies.at(engInd)) + 0.5 * dloge) };

					if (binPartAng > 90.0 && binPartEng >= engmin && binPartEng < engmax && binPartAng >= angmin && binPartAng < angmax)
						ret.at(angleInd).at(engInd) += bsEflux.at(retEbin); //think this is ok...
				}
			}
		}
	}
	
	for (int iii = 0; iii < ret.size(); iii++)
		for (int jjj = 0; jjj < ret.at(0).size(); jjj++)
			ret.at(iii).at(jjj) /= abs(cos(3.14159265 / 180.0 * binAngles.at(iii)));

	std::cout << std::endl << "BS out:" << std::endl;
	for (int angbin = 0; angbin < ret.size(); angbin++)
	{
		for (int ebin = 0; ebin < ret.at(0).size(); ebin++)
			std::cout << ret.at(angbin).at(ebin) << " ";
		std::cout << std::endl;
	}

	return ret;
}

std::vector<std::vector<double>> countBins(std::vector<long> indicies, std::vector<double> energies, std::vector<double> pitches, std::vector<double> ionweights, std::vector<double> magweights, std::vector<double>& binEnergies, std::vector<double>& binAngles)
{
	std::vector<std::vector<double>>   eAngCount(binAngles.size());   //vec[angbin][ebin]

	for (int angbin = 0; angbin < eAngCount.size(); angbin++) //zeros vector array
		eAngCount.at(angbin) = std::vector<double>(binEnergies.size());

	double aBinSize{ 180.0 / 18.0 }; //sets bin min to 0 and bin max to 180 with midbins at 5 + 10*ind
	double logeBinSize{ (4.5 - 0.5) / 47 }; //sets midbins from 0.5 - 4.5 log eV

	std::cout << "energies sample: " << energies.at(0) << " " << energies.at(200) << " pitch sample: " << pitches.at(0) << " " << pitches.at(200) << std::endl;
	std::cout << "energies size: " << energies.size() << std::endl;

	for (int part = 0; part < energies.size(); part++)
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
					(energies.at(part) >= pow(10, logemin)) && (energies.at(part) < pow(10, logemax)) && !found)
				{
					(part < energies.size() / 2) ?
						(eAngCount.at(angbin).at(ebin) += ionweights.at(indicies.at(part) % ionweights.size())) : // / abs(cos(pitches.at(part) * 3.14159265 / 180.0))) : //* (1 - cos(pitches.at(part) * 3.1415926 / 180.0)) / pow(sin(pitches.at(part) * 3.1415926 / 180.0), 2) :
						(eAngCount.at(angbin).at(ebin) += magweights.at(indicies.at(part) % magweights.size()));  // / abs(cos(pitches.at(part) * 3.14159265 / 180.0)));  //* (1 - cos(pitches.at(part) * 3.1415926 / 180.0)) / pow(sin(pitches.at(part) * 3.1415926 / 180.0), 2);
					found = true;
					//eAngCount.at(angbin).at(ebin) += 1;
				}
			}
		}
	}

	for (int iii = 0; iii < eAngCount.size(); iii++)
		for (int jjj = 0; jjj < eAngCount.at(0).size(); jjj++)
			eAngCount.at(iii).at(jjj) /= abs(cos(3.14159265 / 180.0 * binAngles.at(iii)));

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

std::vector<double> maxwellianWeights(std::vector<double> initEnergies, double dEflux_kT, double sigma_kT)
{
	double dloge{ log10(initEnergies.at(1)) - log10(initEnergies.at(0)) };

	double logTeV{ log10(sigma_kT) };
	double binWidth_kT{ pow(10, logTeV + dloge * 0.5) - pow(10, logTeV - dloge * 0.5) };
	double multFactor{ dEflux_kT / (exp(-1.0) * binWidth_kT) };
	
	std::vector<double> ret;
	for (int engInd = 0; engInd < initEnergies.size(); engInd++)
	{
		double eBinWid{ pow(10, dloge * (engInd + 0.5) + 0.5) - pow(10, dloge * (engInd - 0.5) + 0.5) };
		ret.push_back( eBinWid * exp(-initEnergies.at(engInd) / sigma_kT) * multFactor );
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

void getParticlePitchE(std::vector<double> vpara, std::vector<double> vperp, double mass, std::vector<double>& particleEnergies, std::vector<double>& particleAngles)
{
	for (int part = 0; part < vpara.size(); part++)
	{
		bool nonZero{ vpara.at(part) != 0.0 || vperp.at(part) != 0.0 };

		if (nonZero)
		{
			particleEnergies.at(part) = (0.5 * mass * (vpara.at(part) * vpara.at(part) + vperp.at(part) * vperp.at(part)) / 1.60218e-19);
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
	std::string savefolder{ "./../../../../../_dataout/180220_15.55.11/" };

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
	std::vector<double> energies(final_vpara.size());
	getParticlePitchE(final_vpara, final_vperp, 9.10938356e-31, energies, pitches);

	std::vector<double> initpitches(PARTICLECOUNT);
	std::vector<double> initenergies(PARTICLECOUNT);
	getParticlePitchE(init_vpara, init_vperp, 9.10938356e-31, initenergies, initpitches);

	std::cout << "energies done" << std::endl;

	for (int iii = 0; iii < 10 * 96; iii += 96)
	{
		std::cout << "index: " << iii << " pitch: " << pitches.at(iii) << " energy: " << energies.at(iii) << std::endl;
	}

	//count in bins, mult by maxwellian factor
	std::vector<double> initEbins;
	for (int iii = 0; iii < SIMNEBINS; iii++)
		initEbins.push_back(initenergies.at(iii));

	std::vector<double> maxWeightsIon{ maxwellianWeights(initEbins, 1.25e7, 5.0) };
	std::vector<double> maxWeightsMag{ maxwellianWeights(initEbins, 1.00e8, 6000.0) };
	std::vector<double> maxWtsSecndry{ maxwellianWeights(initEbins, 5.00e7, 5.0) };
	std::vector<double> maxWtsMagPlSc;

	for (int iii = 0; iii < maxWeightsMag.size(); iii++)
		maxWtsMagPlSc.push_back(maxWeightsMag.at(iii) + maxWtsSecndry.at(iii));

	std::cout << "Sample max weights in main func: " << maxWeightsIon.at(0) << " " << maxWeightsIon.at(10) << " " << maxWeightsMag.at(0) << " " << maxWeightsMag.at(10) << std::endl;
	
	std::vector<double> binEnergies(CDFNEBINS);
	std::vector<double> binAngles(CDFNANGLEBINS);
	std::vector<std::vector<double>> ang_eCounts{ countBins(indicies, energies, pitches, maxWeightsIon, maxWtsMagPlSc, binEnergies, binAngles) };


	/* Backscatter addition */
	std::vector<double> bottom_vpara(PARTICLECOUNT);
	std::vector<double> bottom_vperp(PARTICLECOUNT);
	loadVparaPerp(savefolder, "bins/satellites/btmElec_", PARTICLECOUNT, bottom_vpara, bottom_vperp);
	std::vector<double> bottom_energies(PARTICLECOUNT);
	std::vector<double> bottom_pitches(PARTICLECOUNT);
	getParticlePitchE(bottom_vpara, bottom_vperp, 9.10938356e-31, bottom_energies, bottom_pitches);
	//std::vector<std::vector<double>> backscatter(std::vector<double> ion_init_energies, std::vector<double> ion_init_pitches, std::vector<double> sat_energies,
		//std::vector<double> sat_pitches, std::vector<double> escape_energies, std::vector<double> binEnergies, std::vector<double> binAngles)
	std::vector<std::vector<double>> bs{ backscatter(initenergies, initpitches, energies, pitches, bottom_energies, maxWeightsMag, binEnergies, binAngles) };

	for (int iii = 0; iii < bs.size(); iii++)
		for (int jjj = 0; jjj < bs.at(0).size(); jjj++)
			ang_eCounts.at(iii).at(jjj) += bs.at(iii).at(jjj);

	/* Prep Data for CDF */
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