#include "postprocess.h"

constexpr double PI_C{ 3.14159265358979323846 };
constexpr double JPEREV_C{ 1.6021766209e-19 };
constexpr double DEGTORAD_C{ PI_C / 180.0 };

namespace postprocess
{
	vecDbl2D steadyFlux(std::string datarootfolder, std::vector<double> pitchBin_Min_Max, int pitchBinNum, std::vector<double> logEBin_Min_Max, int logEBinNum,
		std::vector<double> ionMaxwellian_kT_dEflux, std::vector<double> magMaxwellian_kT_dEflux, double mass, int particlecount, std::vector<double>& binAnglesOut, std::vector<double>& binEnergiesOut)
	{
		/* Load data from files */
		vecDbl2D satdn(5); //vpara, vperp, pitches, energies, time
		vecDbl2D satup(5);
		vecDbl2D bottom(5);
		vecDbl2D init(5);  //vpara, vperp, pitches, energies, s

		for (int iii = 0; iii < satdn.size(); iii++)
		{
			satdn.at(iii)  = std::vector<double>(particlecount);
			satup.at(iii)  = std::vector<double>(particlecount);
			init.at(iii)   = std::vector<double>(particlecount);
			bottom.at(iii) = std::vector<double>(particlecount);
		}
		
		//need to deal with a difference in labeling: (SATNAME)Down == downward facing detector - i.e. upward flux
		//I will change the label over when dealing with flux
		utils::loadvFromDisk(datarootfolder, "bins/satellites/4e6ElecDown_", particlecount, satdn.at(0), satdn.at(1));
		numerical::vToEPitch(satdn.at(0), satdn.at(1), mass, satdn.at(2), satdn.at(3));
		TRYCATCHSTDEXP(fileIO::readDblBin(satdn.at(4), datarootfolder + "/bins/satellites/4e6ElecDown_time.bin", particlecount));

		utils::loadvFromDisk(datarootfolder, "bins/satellites/4e6ElecUp_", particlecount, satup.at(0), satup.at(1));
		numerical::vToEPitch(satup.at(0), satup.at(1), mass, satup.at(2), satup.at(3));
		TRYCATCHSTDEXP(fileIO::readDblBin(satup.at(4), datarootfolder + "/bins/satellites/4e6ElecUp_time.bin", particlecount));

		utils::loadvFromDisk(datarootfolder, "bins/satellites/btmElec_", particlecount, bottom.at(0), bottom.at(1));
		numerical::vToEPitch(bottom.at(0), bottom.at(1), mass, bottom.at(2), bottom.at(3));
		TRYCATCHSTDEXP(fileIO::readDblBin(bottom.at(4), datarootfolder + "/bins/satellites/btmElec_time.bin", particlecount));

		utils::loadvFromDisk(datarootfolder, "bins/particles_init/elec_", particlecount, init.at(0), init.at(1));
		numerical::vToEPitch(init.at(0), init.at(1), mass, init.at(2), init.at(3));
		TRYCATCHSTDEXP(fileIO::readDblBin(init.at(4), datarootfolder + "/bins/particles_init/elec_s.bin", particlecount));

		#ifdef PPCOUTDBG
		std::cout << std::setprecision(6);
		std::cout << "1. Load previous sim data: " << datarootfolder << std::endl;
		std::cout << "Sample loaded vpara, vperp data:" << std::endl;
		std::cout << "SatDown: " << satdn.at(0).at(0) << ", " << satdn.at(1).at(0) << ", SatUp:  " << satup.at(0).at(332650)  << ", " << satup.at(1).at(332650)  << std::endl;
		std::cout << "Init:    " << init.at(0).at(0)  << ", " << init.at(1).at(0)  << ", Bottom: " << bottom.at(0).at(332650) << ", " << bottom.at(1).at(332650) << " Init s: " << init.at(4).at(0) << std::endl;
		std::cout << "1a. Convert to Energy, Pitch:" << std::endl;
		std::cout << "SatDown: " << satdn.at(2).at(0) << ", " << satdn.at(3).at(0) << ", SatUp:  " << satup.at(2).at(332650)  << ", " << satup.at(3).at(332650)  << std::endl;
		std::cout << "Init:    " << init.at(2).at(0)  << ", " << init.at(3).at(0)  << ", Bottom: " << bottom.at(2).at(332650) << ", " << bottom.at(3).at(332650) << std::endl;
		std::cout << std::endl;
		#endif /* PPCOUTDBG */

		/* Generate mid-bin Pitches and Energies */
		std::vector<double> binAngles  { numerical::generatePitchBins(pitchBin_Min_Max.at(0), pitchBin_Min_Max.at(1), pitchBinNum) };
		std::vector<double> binEnergies{ numerical::generateLogSpacedEnergyBins(logEBin_Min_Max.at(0), logEBin_Min_Max.at(1), logEBinNum) };

		#ifdef PPCOUTDBG
		std::cout << "2. Generate bins:" << std::endl;
		std::cout << "Mid-Bin Angles: ";
		for (int ang = 0; ang < binAngles.size(); ang++) { std::cout << binAngles.at(ang) << " "; }
		std::cout << std::endl << "Mid-Bin Energies: ";
		for (int eng = 0; eng < binEnergies.size(); eng++) { std::cout << binEnergies.at(eng) << " "; }
		std::cout << std::endl << std::endl;
		#endif /* PPCOUTDBG */

		/* Test energies */
		/*for (int iii = 0; iii < init.at(3).size(); iii++)
		{
			double TMPEPS{ 1.5e-1 };
			std::cout << std::setprecision(10);
			if (abs(init.at(3).at(iii) - satup.at(3).at(iii)) / init.at(3).at(iii) > TMPEPS && satup.at(3).at(iii) != 0.0)
				std::cout << "satup value off!! ind: " << iii << " val: " << satup.at(3).at(iii) << " should be: " << init.at(3).at(iii) << " err: " << abs(init.at(3).at(iii) - satup.at(3).at(iii)) / init.at(3).at(iii) << std::endl;
			if (abs(init.at(3).at(iii) - satdn.at(3).at(iii)) / init.at(3).at(iii) > TMPEPS && satdn.at(3).at(iii) != 0.0)
				std::cout << "satdn value off!! ind: " << iii << " val: " << satdn.at(3).at(iii) << " should be: " << init.at(3).at(iii) << " err: " << abs(init.at(3).at(iii) - satdn.at(3).at(iii)) / init.at(3).at(iii) << std::endl;
			if (abs(init.at(3).at(iii) - bottom.at(3).at(iii)) / init.at(3).at(iii) > TMPEPS && bottom.at(3).at(iii) != 0.0)
				std::cout << "bottom value off!! ind: " << iii << " val: " << bottom.at(3).at(iii) << " should be: " << init.at(3).at(iii) << " err: " << abs(init.at(3).at(iii) - bottom.at(3).at(iii)) / init.at(3).at(iii) << std::endl;
		}*/
		
		//utils::numerical::coutMinMaxErr(init.at(3), satup.at(3), "satup");
		//utils::numerical::coutMinMaxErr(init.at(3), satdn.at(3), "satdn");
		//utils::numerical::coutMinMaxErr(init.at(3), bottom.at(3), "bottm");
		//exit(1);

		/* Generate Maxwellian counts for particles */
		std::vector<double> ionsphE;
		std::vector<double> magsphE;
		TRYCATCHSTDEXP(numerical::splitIonMagEngs(init.at(4), init.at(3), ionsphE, magsphE)); //split energies into magnetospheric source and ionospheric source

		std::vector<double> maxwCounts(ionsphE.size());    //ionosphere
		std::vector<double> maxwCountsMag(magsphE.size()); //magnetosphere
		for (int iii = 0; iii < ionMaxwellian_kT_dEflux.size(); iii += 2)
		{
			std::vector<double> wtTmp{ numerical::maxwellianCounts(ionsphE, ionMaxwellian_kT_dEflux.at(iii), ionMaxwellian_kT_dEflux.at(iii + 1)) };
			
			for (int jjj = 0; jjj < ionsphE.size(); jjj++)
				maxwCounts.at(jjj) += wtTmp.at(jjj); //sum together all specified maxwellians in source
		}
		for (int iii = 0; iii < magMaxwellian_kT_dEflux.size(); iii += 2)
		{
			std::vector<double> wtTmp{ numerical::maxwellianCounts(magsphE, magMaxwellian_kT_dEflux.at(iii), magMaxwellian_kT_dEflux.at(iii + 1)) };

			for (int jjj = 0; jjj < magsphE.size(); jjj++)
				maxwCountsMag.at(jjj) += wtTmp.at(jjj); //sum together all specified maxwellians in source
		}
		//put arrays together
		maxwCounts.insert(maxwCounts.end(), maxwCountsMag.begin(), maxwCountsMag.end()); //assumes ionosphere particles are first...need to have a way to check because this may not always be the case - maybe pass out min and max from splitIonMagEngs
		
		#ifdef PPCOUTDBG
		std::cout << "3. Generate Maxwellian Weights for particles:" << std::endl;
		std::cout << "First 96 ionsphere weights:" << std::endl;
		for (int energy = 0; energy < 96; energy++) { std::cout << maxwCounts.at(energy) << "  "; if (energy % 8 == 7) { std::cout << std::endl; } }
		std::cout << "First 96 magnetosphere weights:" << std::endl;
		for (int energy = 0; energy < 96; energy++) { std::cout << maxwCounts.at(energy + maxwCounts.size() / 2) << "  "; if (energy % 8 == 7) { std::cout << std::endl; } }
		std::cout << std::endl << std::endl;
		#endif /* PPCOUTDBG */

		for (int iii = 0; iii < init.at(2).size(); iii++)
		{
			if (iii < init.at(2).size() / 2) //ionospheric source
				maxwCounts.at(iii) *= -cos(init.at(2).at(iii) * DEGTORAD_C);
			else //magnetospheric source
				maxwCounts.at(iii) *= 1.0 / cos(satup.at(2).at(iii) * DEGTORAD_C);
		}

		//double dAngSrc{ abs(init.at(2).at(0) - init.at(2).at(100)) };
		//double Bsrc{ 5.73447893584459e-5 }; //B at ionosphere - 100km altitude, 101.322 s
		//for (int iii = 0; iii < init.at(2).size(); iii++)
		//{
			//if (iii == init.at(2).size() / 2) { Bsrc = 8.2358677163765e-7; } //B at magnetosphere - 3Re altitude, 3.12Re s
			//maxwCounts.at(iii) /* * */ /= Bsrc * (cos((init.at(2).at(iii) - 0.5 * dAngSrc) * DEGTORAD_C) - cos((init.at(2).at(iii) + 0.5 * dAngSrc) * DEGTORAD_C));
		//}

		vecDbl2D simfluxupward{ steady::simEnergyFlux(satdn, binAngles, binEnergies, maxwCounts, 9.10938356e-31, -1.6021766209e-19, -1.3597036e-5) }; //calculate binned flux for downward facing detector data, upward flux
		vecDbl2D simfluxdnward{ steady::simEnergyFlux(satup, binAngles, binEnergies, maxwCounts, 9.10938356e-31, -1.6021766209e-19, -1.3597036e-5) }; //calculate binned flux for upward facing detector data
		//vecDbl2D backscatflux { steady::bsEnergyFlux (init, satdn, bottom, binAngles, binEnergies, maxwCounts, 9.10938356e-31, -1.6021766209e-19, -1.3597036e-5) }; //calculate backscatter flux
		
		#ifdef PPCOUTDBG
		std::cout << "4. Calculate upward, downward, and backscatter fluxes per bin:" << std::endl;
		std::cout << "Upward flux by bin:" << std::endl;
		for (int ang = 0; ang < simfluxupward.size(); ang++) { for (int eng = 0; eng < simfluxupward.at(0).size(); eng++) { std::cout << simfluxupward.at(ang).at(eng) << "  "; } std::cout << std::endl; }
		std::cout << "Downward flux by bin:" << std::endl;
		for (int ang = 0; ang < simfluxdnward.size(); ang++) { for (int eng = 0; eng < simfluxdnward.at(0).size(); eng++) { std::cout << simfluxdnward.at(ang).at(eng) << "  "; } std::cout << std::endl; }
		std::cout << "Backscatter flux by bin:" << std::endl;
		//for (int ang = 0; ang < backscatflux.size(); ang++)  { for (int eng = 0; eng < backscatflux.at(0).size(); eng++)  { std::cout << backscatflux.at(ang).at(eng) << "  "; }  std::cout << std::endl; }
		#endif /* PPCOUTDBG */
		
		for (int iii = 0; iii < simfluxupward.size(); iii++)
			for (int jjj = 0; jjj < simfluxupward.at(0).size(); jjj++)
				simfluxupward.at(iii).at(jjj) += simfluxdnward.at(iii).at(jjj) /*+ backscatflux.at(iii).at(jjj)*/;

		#ifdef PPCOUTDBG
		std::cout << "Total flux by bin:" << std::endl;
		for (int ang = 0; ang < simfluxupward.size(); ang++) { for (int eng = 0; eng < simfluxupward.at(0).size(); eng++) { std::cout << simfluxupward.at(ang).at(eng) << "  "; } std::cout << std::endl; }
		std::cout << std::endl;
		#endif /* PPCOUTDBG */

		binAnglesOut = binAngles;
		binEnergiesOut = binEnergies;

		return simfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}

	vecDbl2D timedepFlux()
	{
		throw std::logic_error ("postprocess::timedepFlux: error: time dependent case has not been coded yet.");
		//
		//need to code timedep::satEFlux() and timedep::backscatterFlux() before removing the above throw statement
		//

		//load data
		//convert to energy/pitch
		//generate maxwellian weights
		//find energy fluxes

		vecDbl2D sat{ timedep::simEnergyFlux() };
		vecDbl2D bs { timedep::bsEnergyFlux() };

		for (int iii = 0; iii < sat.size(); iii++)
			for (int jjj = 0; jjj < sat.at(0).size(); jjj++)
				sat.at(iii).at(jjj) += bs.at(iii).at(jjj);

		return sat;
	}

	namespace steady
	{
		vecDbl2D simEnergyFlux(const vecDbl2D& particleData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts,
			double mass, double charge, double BatXSection)
		{
			const std::vector<double>& particlePitches { particleData.at(2) };
			const std::vector<double>& particleEnergies{ particleData.at(3) };
			
			vecDbl2D ret{ numerical::countInBinsWeighted(particlePitches, particleEnergies, binAngles, binEnergies, maxwCounts) };

			//convert to dEflux
			//double Bsat{ 1.35970360638581e-5 }; //B at 4e6 altitude, 4.0713e6 s (along field line from Re)
			//double dAngBin{ binAngles.at(1) - binAngles.at(0) };
			for (int ang = 0; ang < ret.size(); ang++)
			{
				//double binSolidAngle{ cos((binAngles.at(ang) - 0.5 * dAngBin) * DEGTORAD_C) - cos((binAngles.at(ang) + 0.5 * dAngBin) * DEGTORAD_C) };
				for (int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng);// * /* / */ (Bsat * binSolidAngle);
			}

			return ret;
		}

		vecDbl2D bsEnergyFlux(const vecDbl2D& initialData, const vecDbl2D& satData, const vecDbl2D& escapeData,	const std::vector<double>& binAngles,
			const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge, double BatXSection)
		{
			vecDbl2D escCntBins{ numerical::countInBinsWeighted(escapeData.at(2), escapeData.at(3), binAngles, binEnergies, maxwCounts) };
			double Bion{ 5.73447893584459e-5 }; //B at ionosphere - 100km altitude, 101.322 s
			std::vector<double> escCntBinsSum(escCntBins.at(0).size()); //each bin in units of counts -> # particles
			for (int iii = 0; iii < escCntBinsSum.size(); iii++)
			{
				for (int jjj = 0; jjj < escCntBins.size(); jjj++)
					escCntBinsSum.at(iii) += escCntBins.at(jjj).at(iii);
				escCntBinsSum.at(iii) /= (Bion * escCntBins.size()); //isotropically space across pitch angle bins
			}

			/*for (int iii = 0; iii < escCntBinsSum.size(); iii++)
				std::cout << escCntBinsSum.at(iii) << "  ";
			std::cout << std::endl;*/

			std::vector<double> numFluxByBin{ backscat::sumIntegralsOfNumFluxFcnsPerBin(escCntBinsSum, binEnergies, 1.5, -4.0, -2.1, 0.3) }; //obtained by log linefitting Evans, 1974 - these seem closest

			double dAngBin{ binAngles.at(1) - binAngles.at(0) };
			vecDbl2D numFlux(binAngles.size(), std::vector<double>(binEnergies.size()));
			for (int ang = 0; ang < binAngles.size(); ang++)
				for (int eng = 0; eng < binEnergies.size(); eng++)
					numFlux.at(ang).at(eng) = numFluxByBin.at(eng) * Bion * (cos((binAngles.at(ang) - 0.5 * dAngBin) * DEGTORAD_C) - cos((binAngles.at(ang) + 0.5 * dAngBin) * DEGTORAD_C));

			vecDbl2D ret{ backscat::matchIonBSToSatAndCount(numFlux, initialData, satData, binAngles, binEnergies) };
			
			//convert to dEflux
			double Bsat{ 1.35970360638581e-5 }; //B at 4e6 altitude, 4.0713e6 s (along field line from Re)
			for (int ang = 0; ang < ret.size(); ang++)
			{
				double binSolidAngle{ cos((binAngles.at(ang) - 0.5 * dAngBin) * DEGTORAD_C) - cos((binAngles.at(ang) + 0.5 * dAngBin) * DEGTORAD_C) };
				for (int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng) / (Bsat * binSolidAngle);
			}

			return ret;
		}
	}

	namespace timedep
	{
		//to be written
		//for time dependent simulations, such as including Alfven waves and other time dep field phenomena
		vecDbl2D simEnergyFlux()
		{

		}

		vecDbl2D bsEnergyFlux()
		{

		}
	}

	namespace numerical
	{
		void vToEPitch(const std::vector<double>& vpara, const std::vector<double>& vperp, double mass, std::vector<double>& particlePitches, std::vector<double>& particleEnergies)
		{
			for (int part = 0; part < vpara.size(); part++)
			{
				bool nonZero{ vpara.at(part) != 0.0 || vperp.at(part) != 0.0 };

				if (nonZero) //check this or else the function can produce "NaN" in some indicies (I think atan2 is responsible) -> if false, the data at that index will be left 0
				{
					particleEnergies.at(part) = (0.5 * mass * (vpara.at(part) * vpara.at(part) + vperp.at(part) * vperp.at(part)) / JPEREV_C);
					particlePitches.at(part) = atan2(abs(vperp.at(part)), -vpara.at(part)) / DEGTORAD_C;
				}
			}
		}

		std::vector<double> generatePitchBins(double pitchMin, double pitchMax, int numBins)
		{
			/*
			   Returns vector of angles mid-bin
			    |-----|
				 dangle
				-------------------------
				|     |     |     |     |
				-------------------------
				^pitchMin               ^
				                pitchMax^
			*/

			std::vector<double> ret;
			double dangle{ (pitchMax - pitchMin) / (double)numBins };

			for (int ang = 0; ang < numBins; ang++)
				ret.push_back((ang + 0.5) * dangle + pitchMin);

			return ret;
		}

		std::vector<double> generateLogSpacedEnergyBins(double logEmidBinMin, double logEmidBinMax, int numBins)
		{
			/*
			   Returns vector of energy at mid-bin
			   Calculated in log Energy, vector returned in Energy
			    |-----|
				 dlogE
				-------------------------
				|     |     |     |     |
				-------------------------
				   ^logEmidBinMin    ^
				        logEmidBinMax^
			*/

			std::vector<double> ret;
			double dlogE{ (logEmidBinMax - logEmidBinMin) / (numBins - 1) };

			for (int eng = 0; eng < numBins; eng++)
				ret.push_back(pow(10, eng * dlogE + logEmidBinMin));

			return ret;
		}

		void splitIonMagEngs(const std::vector<double>& s_init, const std::vector<double>& E_init, std::vector<double>& ionsphE, std::vector<double>& magsphE)
		{
			double s_min{ 1.0e10 };
			double s_max{ 0.0 };

			for (int part = 0; part < s_init.size(); part++)
			{
				if (s_init.at(part) < s_min)
					s_min = s_init.at(part);
				else if (s_init.at(part) > s_max)
					s_max = s_init.at(part);
			}

			for (int part = 0; part < s_init.size(); part++)
			{
				if (s_init.at(part) == s_min)
					ionsphE.push_back(E_init.at(part));
				else if (s_init.at(part) == s_max)
					magsphE.push_back(E_init.at(part));
			}

			if (ionsphE.size() + magsphE.size() != E_init.size()) //later, this condition may not apply - if particles are generated mid-sim, for example
				throw std::runtime_error("postprocess::numerical::splitIonMagEngs: sum number of entries in ion, mag energy arrays is not equal to the number of entries in initial energy array - this means that some particles were not detected -> "
					+ std::to_string(ionsphE.size()) + " + " + std::to_string(magsphE.size()) + " != " + std::to_string(E_init.size()));
		}

		std::vector<double> maxwellianCounts(const std::vector<double>& initEnergies, double sigma_kT, double dEflux_kT) //or maybe maxwellian::getWeights
		{
			double logEBinMin{ log10(initEnergies.at(0)) };
			double dlogE{ log10(initEnergies.at(1)) - logEBinMin };

			double logTeV{ log10(sigma_kT) };
			double binWidth_kT{ pow(10, logTeV + dlogE * 0.5) - pow(10, logTeV - dlogE * 0.5) };
			double multFactor{ dEflux_kT / (exp(-1.0) * binWidth_kT) };

			std::vector<double> ret(initEnergies.size());
			for (int engInd = 0; engInd < initEnergies.size(); engInd++)
			{
				if (initEnergies.at(engInd) == 0.0) { continue; }
				double eBinWid{ pow(10, log10(initEnergies.at(engInd)) + 0.5 * dlogE) - pow(10, log10(initEnergies.at(engInd)) - 0.5 * dlogE) };
				ret.at(engInd) = eBinWid * exp(-initEnergies.at(engInd) / sigma_kT) * multFactor / initEnergies.at(engInd);
			}

			return ret;
		}

		vecDbl2D countInBinsWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts)
		{
			vecDbl2D ret(binAngles.size());
			//vecDbl2D tmp(binAngles.size());

			double logEMinBinMid{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEMinBinMid };
			double dangle{ binAngles.at(1) - binAngles.at(0) };
			double allEmin{ pow(10, logEMinBinMid - 0.5 * dlogE) };
			double allEmax{ pow(10, log10(binEnergies.at(binEnergies.size() - 1)) + 0.5 * dlogE) };
			
			for (int part = 0; part < particleEnergies.size(); part++)
			{
				double partEnerg{ particleEnergies.at(part) };
				double partPitch{ particlePitches.at(part) };
				
				if (part != 0 && (partEnerg < allEmin || partEnerg > allEmax)) { continue; } //should speed up in case of zero, below min, or above max

				for (int angbin = 0; angbin < binAngles.size(); angbin++)
				{
					if (part == 0) { ret.at(angbin) = std::vector<double>(binEnergies.size()); }
					//if (part == 0) { tmp.at(angbin) = std::vector<double>(binEnergies.size()); }

					double angmin{ angbin *       dangle }; //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
					double angmax{ (angbin + 1) * dangle };
					if (angbin == (binAngles.size() - 1)) { angmax += 0.001; } //I want the top angle bin to include 180 and right now the conditional is ... < 180, not ...<= 180

					for (int ebin = 0; ebin < binEnergies.size(); ebin++)
					{
						double emin{ pow(10, (ebin - 0.5) * dlogE + logEMinBinMid) }; //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
						double emax{ pow(10, (ebin + 0.5) * dlogE + logEMinBinMid) };

						if ((partPitch >= angmin) && (partPitch < angmax) &&
							(partEnerg >= emin)   && (partEnerg < emax))
						{
							ret.at(angbin).at(ebin) += maxwCounts.at(part); //maxwCounts needs to be as long as particleEnergies/particlePitches
							//tmp.at(angbin).at(ebin) += 1.0;
							angbin = binAngles.size(); //no need to iterate over other bins once particle bin is found
							ebin = binEnergies.size();
						}
					}
				}
			}

			for (int iii = 0; iii < ret.size(); iii++)
			{
				for (int jjj = 0; jjj < ret.at(0).size(); jjj++)
					std::cout << ret.at(iii).at(jjj) << "  ";
				std::cout << std::endl;
			}
			std::cout << std::endl;

			//std::cout << "returning tmp (not ret)" << std::endl;
			return ret;
		}

		void countsToEFlux(vecDbl2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection)
		{
			double dlogE{ log10(binEnergies.at(1)) - log10(binEnergies.at(0)) };

			for (int ang = 0; ang < binAngles.size(); ang++)
			{
				for (int eng = 0; eng < binEnergies.size(); eng++)
				{
					double v_perp{ sqrt(2 * binEnergies.at(eng) * JPEREV_C / mass) };
					double A_gyro{ PI_C * pow((mass * v_perp / (charge * BatXSection)), 2) }; //don't need to absolute value because it's squared
					double binWidth{ pow(10, log10(binEnergies.at(eng)) + 0.5 * dlogE) - pow(10, log10(binEnergies.at(eng)) - 0.5 * dlogE) };
					//double solidAngle{ write this }; //involves pitch angle
					double solidAngle{ 1.0 };

					energyData.at(ang).at(eng) *= abs(cos(binAngles.at(ang) * DEGTORAD_C)) / (A_gyro * solidAngle * binWidth);
				}
			}
		}

		void divBinsByCosPitch(vecDbl2D& data, std::vector<double> binAnglesDegrees)
		{
			for (int iii = 0; iii < data.size(); iii++)
				for (int jjj = 0; jjj < data.at(0).size(); jjj++)
					data.at(iii).at(jjj) /= abs(cos(DEGTORAD_C * binAnglesDegrees.at(iii)));
		}
	} //end namespace postprocess::numerical

	namespace backscat //postprocess::backscat
	{
		/*
			Based on equations:
			BS_scnd(x) = 10 ^ (scnd_logm * logE + scnd_logb)
				=== d_sec * x ^ (scnd_logm), where d_sec = 10 ^ (scnd_logb)
			BS_prim(x) = 10 ^ (prim_logm * log(E / E_incident) + prim_logb) * (10000 / E_incident)
				=== d_pri * x ^ (prim_logm), where d_pri = 10 ^ (prim_logb + 4) / E_incident ^ (prim_logm + 1)

			Integral:
			BS'_scnd(x) = d_sec / (scnd_logm + 1) * x ^ (scnd_logm + 1)
			BS'_prim(x) = d_pri / (prim_logm + 1) * x ^ (prim_logm + 1), if x > E_incident, BS_scnd||prim = 0
		*/

		double F_flux(double evalE, double incidentE, double incidentCnt, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb) { //describes a log linefit of the backscatter curves detailed in Evans, 1974
			return (pow(10, scnd_logb) * pow(evalE, scnd_logm) + pow(10, prim_logb + 4) / pow(incidentE, prim_logm + 1) * pow(evalE, prim_logm)) * incidentCnt; }

		double integralF_flux(double lower, double upper, double incidentE, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{
			if (upper > incidentE * (1 + FLT_EPSILON)) { throw std::invalid_argument("integralF_bs: upper limit of integration is higher than incident energy - function is zero above incident energy - upper limit: " + std::to_string(upper) + ", incidentE: " + std::to_string(incidentE)); }
			double integral_sec{ (pow(upper, scnd_logm + 1) - pow(lower, scnd_logm + 1)) * pow(10, scnd_logb) / (scnd_logm + 1) };
			double integral_prm{ (pow(upper, prim_logm + 1) - pow(lower, prim_logm + 1)) * pow(10, prim_logb + 4) / ((prim_logm + 1) * pow(incidentE, prim_logm + 1)) };
			return integral_sec + integral_prm;
		}

		std::vector<double> sumIntegralsOfNumFluxFcnsPerBin(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb)
		{
			double logEBinMin{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEBinMin };

			std::vector<double> bsNflux(binEnergies.size());
			for (int ampbin = 0; ampbin < bsNflux.size(); ampbin++)
			{
				double engmin{ pow(10, (ampbin - 0.5) * dlogE + logEBinMin) }; //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
				double engmax{ pow(10, (ampbin + 0.5) * dlogE + logEBinMin) };

				for (int cntbin = ampbin; cntbin < bsNflux.size(); cntbin++)
				{
					double incidentE{ pow(10, log10(binEnergies.at(cntbin)) + 0.5 * dlogE) }; //incident E is upper limit of bin - particles in here can be up to this energy, so why not use this instead of mid bin?
					double intF{ integralF_flux(engmin, engmax, incidentE, primary_logm, primary_logb, secondary_logm, secondary_logb) };
					
					bsNflux.at(ampbin) += intF * binCounts.at(cntbin) / (engmax - engmin);
				}
			}
			
			return bsNflux;
		}

		vecDbl2D matchIonBSToSatAndCount(const vecDbl2D& bsNumFluxBins, const vecDbl2D& initialData, const vecDbl2D& satDownData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies)
		{
			//double ang_epsilon{ 0.026 }; //these can be simulation specific
			//double eng_epsilon{ 0.05 };  //if more bins are used, these may not yield the closest result

			vecDbl2D ret(binAngles.size());
			for (int iii = 0; iii < ret.size(); iii++)
				ret.at(iii) = std::vector<double>(binEnergies.size());

			double logEMinBinMid{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEMinBinMid };
			double dangle{ binAngles.at(1) - binAngles.at(0) };

			#ifdef PPCOUTDBG_VERBOSE
			std::cout << "4inner. Match bin particle to closest sim particle:" << std::endl;
			#endif /* PPCOUTDBG_VERBOSE */

			for (int ionAngBin = binAngles.size() / 2; ionAngBin < binAngles.size(); ionAngBin++)
			{ //ionAngBin starts at 95 degrees
				for (int ionEngBin = 0; ionEngBin < binEnergies.size(); ionEngBin++)
				{
					bool found{ false };
					int partInd{ -1 }; //index of the sim particle that is closest to the bin particle
					double angDiff{ 1e10 };
					double engDiff{ 1e10 };

					double tmpBinAng{ binAngles.at(ionAngBin) };
					double tmpBinEng{ binEnergies.at(ionEngBin) };

					//find initial particle that is closest to bin particle
					for (int part = 0; part < initialData.at(2).size(); part++)
					{ //don't know if I'm wild about iterating over every particle for every bin particle, but it's the most general and should work regardless of sim particle arrangement
						double tmpAngDiff{ abs(tmpBinAng - initialData.at(2).at(part)) };
						double tmpEngDiff{ abs(tmpBinEng - initialData.at(3).at(part)) };

						if (tmpAngDiff * (1 - FLT_EPSILON) <= angDiff && tmpEngDiff * (1 - FLT_EPSILON) <= engDiff) { angDiff = tmpAngDiff; engDiff = tmpEngDiff; partInd = part;/*std::cout << partInd << ":" << engDiff << "," << angDiff << ":" << tmpBinEng << "," << tmpBinAng << std::endl;*/  }
						if (tmpAngDiff < FLT_EPSILON && tmpEngDiff < FLT_EPSILON) { part = initialData.at(2).size(); std::cout << "Sys EPS" << std::endl; } //if both differences are less than FLT_EPSILON - ~1e-7 on Windows 10, it's close enough - there's no bin 1e-7eV wide, nor is there likely to ever be
						//if (tmpAngDiff < ang_epsilon && tmpEngDiff / tmpBinEng < eng_epsilon) { part = initialData.at(2).size(); /*std::cout << "My EPS" << std::endl;*/ } //ends iterations if sim particle representing bin particle is close enough as defined above
					}

					double foundAng{ satDownData.at(2).at(partInd) }; //angle of the matching particle at satellite
					double foundEng{ satDownData.at(3).at(partInd) }; //energy of the matching particle at satellite

					#ifdef PPCOUTDBG_VERBOSE
					std::cout << "Bin: " << binEnergies.at(ionEngBin) << "," << binAngles.at(ionAngBin) << " Init: " << initialData.at(3).at(partInd) << "," << initialData.at(2).at(partInd) << " Sim: " << foundEng << "," << foundAng << std::endl;
					#endif /* PPCOUTDBG_VERBOSE */

					//assign flux at ionosphere in a bin to the appropriate bin at satellite
					for (int satAngBin = 0; satAngBin < binAngles.size(); satAngBin++)
					{
						double tmpAngMin{ (satAngBin)     * dangle }; //assumes starting at 0 degrees, pretty reasonable assumption
						double tmpAngMax{ (satAngBin + 1) * dangle };

						if (satAngBin == binAngles.size() - 1) { tmpAngMax += 0.001; } //includes 180 degrees explicitly

						for (int satEngBin = 0; satEngBin < binEnergies.size(); satEngBin++)
						{
							double tmpEngMin{ pow(10, (satEngBin - 0.5) * dlogE + logEMinBinMid) };
							double tmpEngMax{ pow(10, (satEngBin + 0.5) * dlogE + logEMinBinMid) };

							if (foundAng >= tmpAngMin && foundAng < tmpAngMax &&
								foundEng >= tmpEngMin && foundEng < tmpEngMax)// && !found)
							{
								ret.at(satAngBin).at(satEngBin) += bsNumFluxBins.at(ionAngBin).at(ionEngBin);
								found = true;
								satAngBin = binAngles.size();  //prevents unnecessary iteration
								satEngBin = binEnergies.size();//and also double counting
							}
						}
					}

					//check to make sure each bin particle was found in a bin at satellite, if not, throw
					if (!found)
						throw std::runtime_error("postprocess::backscat::matchIonBSToSatAndCount: sim particle was not added to bin: ionosphere Eng,Ang: "
							+ std::to_string(ionAngBin) + "," + std::to_string(ionEngBin) + ":: matching sim particle Eng,Ang: " + std::to_string(foundEng) + "," + std::to_string(foundAng));
				}
			}

			return ret;
		}
	} //end namespace postprocess::backscat

	namespace utils //postprocess::utils
	{
		void loadvFromDisk(std::string savefolder, std::string filePrefix, int count, std::vector<double>& vpara, std::vector<double>& vperp)
		{
			try
			{
				fileIO::readDblBin(vpara, savefolder + filePrefix + "vpara.bin", count);
				fileIO::readDblBin(vperp, savefolder + filePrefix + "vperp.bin", count);
			}
			catch (std::exception& exp)
			{
				std::cout << exp.what() << " Exiting." << std::endl;
				exit(1);
			}
		}
	} //end namespace postprocess::utils
}