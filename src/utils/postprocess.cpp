#include "utils/postprocess.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

using std::pow;

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(const PPData& ppdata)
	{
		// 1. Adjust Maxwellian by Cos, SQRT(B) Factors
		double Aratio_ion_sat{ std::sqrt(ppdata.B_ion / ppdata.B_alt) };
		double Aratio_mag_sat{ std::sqrt(ppdata.B_alt / ppdata.B_mag) };
		double Aratio_ion_bs { std::sqrt(ppdata.B_ion / ppdata.B_ion) * std::sqrt(ppdata.B_alt / ppdata.B_ion) };
		double Aratio_mag_bs { std::sqrt(ppdata.B_ion / ppdata.B_mag) * std::sqrt(ppdata.B_alt / ppdata.B_ion) };
		double bsScale{ 0.1 }; //for now, arbitrary factor to get in the ballpark

		std::vector<double> satWeights{ ppdata.maxCounts };
		std::vector<double> bsWeights { ppdata.maxCounts };

		for (unsigned int iii = 0; iii < ppdata.initial.s_pos.size(); iii++) //isotropize counts -> 3D
		{
			if (ppdata.initial.s_pos.at(iii) < ppdata.s_ion * 1.001)     //ionospheric source
			{
				satWeights.at(iii) *= -cos(ppdata.initial.pitch.at(iii) * RADS_PER_DEG) * Aratio_ion_sat;
				bsWeights.at(iii)  *= -cos(ppdata.initial.pitch.at(iii) * RADS_PER_DEG) * Aratio_ion_bs * bsScale; //without QSPS, there shouldn't be any ionospheric-source particles influencing the backscatter
			}
			else if (ppdata.initial.s_pos.at(iii) > ppdata.s_mag * 0.999)//magnetospheric source
			{
				satWeights.at(iii) *= 1.0 / cos(ppdata.dnward.pitch.at(iii) * RADS_PER_DEG) * Aratio_mag_sat;
				bsWeights.at(iii)  *= 1.0 / cos(ppdata.bottom.pitch.at(iii) * RADS_PER_DEG) * Aratio_mag_bs * bsScale;
			}
			else
				throw std::logic_error("postprocess::steadyFlux : particle is not ionospheric or magnetospheric source");
		}

		// 2. Calculate dEfluxes
		dblVec2D distfluxdnward{ EFlux::satdEFlux(ppdata.dnward, ppdata.pitchBins, ppdata.energyBins, satWeights) };
		dblVec2D distfluxupward{ EFlux::satdEFlux(ppdata.upward, ppdata.pitchBins, ppdata.energyBins, satWeights) };
		dblVec2D backfluxupward{ EFlux::bksdEFlux(ppdata.initial, ppdata.upward, ppdata.bottom, ppdata.pitchBins, ppdata.energyBins, bsWeights) };
		
		// 3. Sum dEfluxes
		for (unsigned int iii = 0; iii < distfluxupward.size(); iii++)
			for (unsigned int jjj = 0; jjj < distfluxupward.at(iii).size(); jjj++)
				distfluxupward.at(iii).at(jjj) += distfluxdnward.at(iii).at(jjj) + backfluxupward.at(iii).at(jjj);

		return distfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satUpwardData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies)
		{
			const double ANGMAXERR{ 0.1 }; //max error in degrees
			const double ENYMAXERR{ 0.1 }; //max error in eV
			auto err = [](double base, double diff) { return std::abs((base - diff) / base); };

			dblVec2D ret(binAngles.size());
			for (unsigned int iii = 0; iii < ret.size(); iii++)
				ret.at(iii) = std::vector<double>(binEnergies.size());

			const double logEMinBinMid{ log10(binEnergies.at(0)) };
			const double dlogE{ log10(binEnergies.at(1)) - logEMinBinMid };
			const double dangle{ binAngles.at(1) - binAngles.at(0) };

			std::vector<double> pitches(binAngles.size() * binEnergies.size());
			std::vector<double> energies(binAngles.size() * binEnergies.size());
			std::vector<double> weights(binAngles.size() * binEnergies.size());

			for (unsigned int ionAngBin = binAngles.size() / 2; ionAngBin < binAngles.size(); ionAngBin++) //some other time, just get the repeating energies and angles, fit between those, calculate index
			{ //ionAngBin starts at 95 degrees
				for (unsigned int ionEnyBin = 0; ionEnyBin < binEnergies.size(); ionEnyBin++)
				{
					int partInd{ -1 }; //index of the sim particle that is closest to the bin particle
					double angBestErr{ 1e10 };
					double enyBestErr{ 1e10 };

					double tmpBinAng{ binAngles.at(ionAngBin) };
					double tmpBinEny{ binEnergies.at(ionEnyBin) };

					//find initial particle that is closest to bin particle
					for (unsigned int part = 0; part < initialData.pitch.size(); part++)
					{ //don't know if I'm wild about iterating over every particle for every bin particle, but it's the most general and should work regardless of sim particle arrangement
						double tmpAngErr{ err(tmpBinAng, initialData.pitch.at(part)) };
						double tmpEnyErr{ err(tmpBinEny, initialData.energy.at(part)) };

						//int printEvery{ 96 * 1000 };
						//if (part % printEvery == 0)
						//{
							//std::cout << initialData.pitch.at(part) << ", " << initialData.energy.at(part) << "  :  " << angBestErr << ", " << enyBestErr
								//<< " : " << tmpAngErr << ", " << tmpEnyErr << " : " << partInd << ", "
								//<< "  :  " << ((tmpAngErr * (1 - FLT_EPSILON) <= angBestErr) && (tmpEnyErr * (1 - FLT_EPSILON) <= enyBestErr))
								//<< ", " << (tmpAngErr * (1 - FLT_EPSILON) <= angBestErr) << ", " << (tmpEnyErr * (1 - FLT_EPSILON) <= enyBestErr) << "\n";
						//}

						if ((tmpAngErr * (1 - FLT_EPSILON) <= angBestErr) && (tmpEnyErr * (1 - FLT_EPSILON) <= enyBestErr))
						{
							angBestErr = tmpAngErr;
							enyBestErr = tmpEnyErr;
							partInd = (int)part;
						}

						//if (part == 10000 * printEvery) exit(1);
					}

					pitches.at (ionAngBin * binEnergies.size() + ionEnyBin) = satUpwardData.pitch.at(partInd);
					energies.at(ionAngBin * binEnergies.size() + ionEnyBin) = satUpwardData.energy.at(partInd);
					weights.at (ionAngBin * binEnergies.size() + ionEnyBin) = bsNumFluxBins.at(ionAngBin).at(ionEnyBin);
					//std::cout << "BinE, FoundE: " << tmpBinEny << ", " << initialData.energy.at(partInd) << "  :  BinA, FoundA: " << tmpBinAng << ", " << initialData.pitch.at(partInd) << "  :  PartInd: " << partInd << "\n";
					//std::cout << "Err: A, E: " << angBestErr << ", " << enyBestErr << "\n";

					//throw if over max err
					if (angBestErr > ANGMAXERR)
						throw std::runtime_error("postprocess::steady::bsSrcToSat - Minimum angle error is more than max error allowed - something is wrong: Best Err, Max Err " + std::to_string(angBestErr) + ", " + std::to_string(ANGMAXERR));
					if (enyBestErr > ENYMAXERR)
						throw std::runtime_error("postprocess::steady::bsSrcToSat - Minimum energy error is more than max error allowed - something is wrong: Best Err, Max Err " + std::to_string(enyBestErr) + ", " + std::to_string(ENYMAXERR));
				}
			}

			//exit(1);
			return binning::binWeighted(pitches, energies, binAngles, binEnergies, weights);
		}
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satdEFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWeight)
		{
			// Section 1 - Get dNflux at Satellite
			// 1. Bin Particles by PA, Energy
			dblVec2D ret{ binning::binWeighted(sat.pitch, sat.energy, binAngles, binEnergies, numWeight) };
			// Section 1 End

			//diagnostic - to print the vector and verify accuracy
			auto printVec2D = [](const dblVec2D& prt, std::string name)
			{ //lambda function to print the results
				std::cout << name << ":\n";
				for (auto& dblVec : prt)
				{
					for (auto& elem : dblVec)
						std::cout << elem << ",";
					std::cout << "\n";
				}
				std::cout << "\n";
			};
			
			printVec2D(ret, "nFlux");
			//end diagnostic - can probably get rid of once everything is verified

			// 2. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng);

			return ret;
		}

		DLLEXP dblVec2D bksdEFlux(const ParticleData& initialData, const ParticleData& satData, const ParticleData& escapeData,
			const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts)
		{
			//diagnostic
			auto printVec2D = [](const dblVec2D& prt, std::string name)
			{ //lambda function to print the results
				std::cout << name << ":\n";
				for (auto& dblVec : prt)
				{
					for (auto& elem : dblVec)
						std::cout << elem << ",";
					std::cout << "\n";
				}
				std::cout << "\n";
			};
			//end diagnostic

			// Section 1 - Get dNflux at Satellite
			// 1. Bin Escaped Particles by PA, Energy
			dblVec2D escapeEPABins{ binning::binWeighted(escapeData.pitch, escapeData.energy, binAngles, binEnergies, maxwCounts) };
			//printVec2D(escapeEPABins, "Escaped Particles (dNflux) per E, PA bin");
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle

			// 2. Sum dNflux over PA Bins, Per E Bin and Average
			std::vector<double> escapeEBins(escapeEPABins.at(0).size()); //Sum of escaped particles at each energy, units of dNflux
			for (unsigned int egy = 0; egy < escapeEBins.size(); egy++)  //iterate over energies
			{
				for (auto& ang : escapeEPABins)                          //sum over angles, add to sum vector at energyBin
					escapeEBins.at(egy) += ang.at(egy);
				escapeEBins.at(egy) /= (double)escapeEPABins.size() / 2; //isotropically space across pitch angle bins - divide by # ion PA bins, this value later put in each ionospheric angle bin
			}
			//printVec2D({ escapeEBins }, "Escaped Particles (dNflux) per E bin");
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins

			// 3. Calculate BS dNflux from dNflux Incident to Ionosphere
			std::vector<double> bsdNFluxByEBin{ backscat::dNflux_bs(escapeEBins, binEnergies, 1.5, -4.0, -2.1, 0.3) }; //obtained by log linefitting Evans, 1974 - these seem closest
			//printVec2D({ bsdNFluxByEBin }, "BS NFlux per E bin at Ionosphere");
			// output: 1D vector of backscatter dNFlux(?) by energy bin at ionosphere

			// 4. Distribute BS dNflux Equally Over Pitch Bins
			dblVec2D numFlux(binAngles.size());
			for (unsigned int ang = 0; ang < numFlux.size(); ang++)
			{
				if (ang < numFlux.size() / 2)
					numFlux.at(ang) = std::vector<double>(binEnergies.size());
				else
				{
					numFlux.at(ang) = bsdNFluxByEBin;
					for (unsigned int eny = 0; eny < numFlux.at(ang).size(); eny++)
						numFlux.at(ang).at(eny) *= -cos(binAngles.at(ang) * RADS_PER_DEG);
				}

			}
			//printVec2D(numFlux, "BS NFlux at Ionosphere");
			// output: 2D vector of bs dNFlux at ionosphere per pitch, energy bin - should only be upward (90-180)

			// 5. Translate BS dNflux at Source to dNflux at Satellite
			dblVec2D ret{ steady::bsSrcToSat(numFlux, initialData, satData, binAngles, binEnergies) };
			//printVec2D(ret, "satNFlux");
			// output: 2D vector of bs dNFlux at satellite per pitch, energy bin - should only be upward (90-180)
			// Section 1 End

			// 6. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
			{
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng);
			}
			printVec2D(ret, "bs dEflux");
			// output: 2D vector of bs dEflux

			return ret;
		}
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& counts)
		{
			dblVec2D ret(binAngles.size(), std::vector<double>(binEnergies.size()));

			double logMinEBinMid{ log10(binEnergies.at(0)) };
			double dlogE_bin{ log10(binEnergies.at(1)) - logMinEBinMid };
			double dangle_bin{ binAngles.at(1) - binAngles.at(0) };

			int outsideLimits{ 0 };

			for (unsigned int part = 0; part < particleEnergies.size(); part++) //iterate over particles
			{
				double partEnerg{ particleEnergies.at(part) };
				double partPitch{ particlePitches.at(part) };

				if (partEnerg == 0.0 && partPitch == 0.0) continue; //guards - if particle E, PA is zero, it wasnt detected - just skip it
				if (partEnerg > pow(10, log10(binEnergies.back() ) + (0.5 * dlogE_bin)) || //if particle is outside E, PA measured limits - skip it
					partEnerg < pow(10, log10(binEnergies.front()) - (0.5 * dlogE_bin)) || //PA shouldn't be an issue, there just in case
					partPitch > binAngles.back()  + (0.5 * dangle_bin) ||                  //also, code counts how many are outside limits for diagnostics
					partPitch < binAngles.front() - (0.5 * dangle_bin))
				{ outsideLimits++; continue; }

				int angbin{ (int)(partPitch / dangle_bin) }; //this should give the bin index
				int engbin{ (int)((log10(partEnerg) - (logMinEBinMid - 0.5 * dlogE_bin)) / dlogE_bin) }; //ditto

				if (partEnerg >= pow(10, (engbin - 0.5) * dlogE_bin + logMinEBinMid) &&
					partEnerg < pow(10, (engbin + 0.5) * dlogE_bin + logMinEBinMid) &&
					partPitch >= angbin * dangle_bin &&
					partPitch < (angbin + 1) * dangle_bin)
				{
					ret.at(angbin).at(engbin) += counts.at(part);
					continue;
				}
				else //this shouldn't ever execute, guards should prevent zero and out of limits values
					throw std::logic_error(std::to_string(part) + " >> " + std::to_string(particleEnergies.at(part)) + ":" + std::to_string(engbin) + ",  " + std::to_string(particlePitches.at(part)) + ":" + std::to_string(angbin));
			}

			if (outsideLimits > 0) std::cout << "postprocess::binning::binWeighted : Particles out of limits: " << outsideLimits << "\n";

			return ret;
		}

		DLLEXP void symmetricBins0To360(dblVec2D& data, std::vector<double>& binAngles) //takes bins from 0-180 and extends to 360 symmetrically reflected across 180 degrees
		{//data[angle][energy]
			binAngles.resize(2 * binAngles.size()); //double size of angles vector (goes from 0-180 now)
			data.resize(binAngles.size()); //set data (vector of double vectors) to size of binAngles (now double sized)
			
			for (unsigned int iii = binAngles.size() / 2; iii < binAngles.size(); iii++)
			{
				binAngles.at(iii) = 2 * binAngles.at(iii - 1) - binAngles.at(iii - 2);
				data.at(iii) = data.at(binAngles.size() - iii - 1);
			}
		} //needs a good testing
	} //end namespace postprocess::numerical

	namespace backscat //postprocess::backscat
	{
		/*
			Based on equations:
			BS_scnd(E_eval) = 10 ^ (scnd_logm * logE_eval + scnd_logb)
				=== d_sec * E_eval ^ (scnd_logm), where d_sec = 10 ^ (scnd_logb)

			BS_prim(E_eval) = 10 ^ (prim_logm * log(E_eval / E_incident) + prim_logb) * (10000 / E_incident)
				=== d_pri * E_eval ^ (prim_logm), where d_pri = 10 ^ (prim_logb + 4) / E_incident ^ (prim_logm + 1)

			Integral:
			BS'_scnd(x) = d_sec / (scnd_logm + 1) * x ^ (scnd_logm + 1)
			BS'_prim(x) = d_pri / (prim_logm + 1) * x ^ (prim_logm + 1), if x > E_incident, BS_scnd||prim = 0
		*/

		/*double john_flux(double E_eval, double E_incident, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{
			double scnd{ 3.62 * pow(E_eval, -1.947) };
			double prim{ 0.000183 * pow(E_eval / E_incident, 6.015) * 10000.0 / E_incident };
			return prim + scnd;
		}*/

		DLLEXP double evans_flux(double E_eval, double E_incident, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{// Describes a log linefit of the backscatter curves detailed in Evans, 1974
		 // E_incident defines the primary BS curve
		 // E_eval sets where we want to know the value of the function
			return (pow(10.0, scnd_logb) * pow(E_eval, scnd_logm) + //secondary BS +
				    pow(10.0, prim_logb + 4.0) / pow(E_incident, prim_logm + 1.0) * pow(E_eval, prim_logm)); //primary BS
		}

		DLLEXP double integralEvans_flux(double lower, double upper, double E_incident, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{
			//if (upper > E_incident * (1 + FLT_EPSILON))
				//throw std::invalid_argument("integralF_bs: upper limit of integration is higher than incident energy - function is zero above incident energy - upper limit: " + std::to_string(upper) + ", incidentE: " + std::to_string(E_incident));
			
			double integral_sec{ (pow(upper, scnd_logm + 1.0) - pow(lower, scnd_logm + 1.0)) * pow(10.0, scnd_logb) / (scnd_logm + 1.0) };
			double integral_prm{ (pow(upper, prim_logm + 1.0) - pow(lower, prim_logm + 1.0)) * pow(10.0, prim_logb + 4.0) / ((prim_logm + 1.0) * pow(E_incident, prim_logm + 1.0)) };
			return integral_sec + integral_prm;
		}

		DLLEXP std::vector<double> dNflux_bs(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb)
		{ //converts incident dNflux to bs Nflux
			double logEBinMin{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEBinMin };                           //depends on a logarithmically spaced E, won't work otherwise

			std::vector<double> binCntxE{ binCounts };
			for (unsigned int ebin = 0; ebin < binCounts.size(); ebin++)
				binCntxE.at(ebin) *= binEnergies.at(ebin);

			std::vector<double> bsdNflux(binEnergies.size());
			for (unsigned int dNFluxBin = 0; dNFluxBin < bsdNflux.size(); dNFluxBin++)       //bins that contain the number flux of the backscatter in the energy bin of the same index
			{
				double engmin{ pow(10, (dNFluxBin - 0.5) * dlogE + logEBinMin) };            //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
				double engmax{ pow(10, (dNFluxBin + 0.5) * dlogE + logEBinMin) };

				for (unsigned int incEbin = dNFluxBin; incEbin < bsdNflux.size(); incEbin++) //nFlux bins are contributed to by all incident particles with E higher than the E associated with the nFlux bin
				{//change to mid bin
					double incidentE{ binEnergies.at(incEbin) };                             //incident E is upper limit of bin
					double intF{ integralEvans_flux(engmin, engmax, incidentE, primary_logm, primary_logb, secondary_logm, secondary_logb) };

					bsdNflux.at(dNFluxBin) += intF * binCntxE.at(incEbin) / (engmax - engmin);
					//bsdNflux.at(dNFluxBin) += john_flux(binEnergies.at(dNFluxBin), binEnergies.at(incEbin), primary_logm, primary_logb, secondary_logm, secondary_logb) * binCntxE.at(incEbin);
				}
			}

			return bsdNflux;
		}
	} //end namespace postprocess::backscat
}
