#include "utils/postprocess.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(PPData ppdata)
	{
		dblVec2D distfluxdnward{ EFlux::satEFlux(ppdata.dnward, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts) };
		dblVec2D distfluxupward{ EFlux::satEFlux(ppdata.upward, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts) };
		dblVec2D backfluxupward{ EFlux::backEFlux(ppdata.initial, ppdata.upward, ppdata.bottom, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts, ppdata.mass, ppdata.charge) }; //calculate backscatter flux
		
		for (unsigned int iii = 0; iii < distfluxupward.size(); iii++)
			for (unsigned int jjj = 0; jjj < distfluxupward.at(iii).size(); jjj++)
				distfluxupward.at(iii).at(jjj) += distfluxdnward.at(iii).at(jjj) + backfluxupward.at(iii).at(jjj);

		return distfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satUpwardData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies)
		{
			dblVec2D ret(binAngles.size());
			for (unsigned int iii = 0; iii < ret.size(); iii++)
				ret.at(iii) = std::vector<double>(binEnergies.size());

			double logEMinBinMid{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEMinBinMid };
			double dangle{ binAngles.at(1) - binAngles.at(0) };

			std::vector<double> pitches(binAngles.size() * binEnergies.size());
			std::vector<double> energies(binAngles.size() * binEnergies.size());
			std::vector<double> weights(binAngles.size() * binEnergies.size());

			for (unsigned int ionAngBin = binAngles.size() / 2; ionAngBin < binAngles.size(); ionAngBin++) //some other time, just get the repeating energies and angles, fit between those, calculate index
			{ //ionAngBin starts at 95 degrees
				for (unsigned int ionEngBin = 0; ionEngBin < binEnergies.size(); ionEngBin++)
				{
					bool found{ false };
					int partInd{ -1 }; //index of the sim particle that is closest to the bin particle
					double angDiff{ 1e10 };
					double engDiff{ 1e10 };

					double tmpBinAng{ binAngles.at(ionAngBin) };
					double tmpBinEng{ binEnergies.at(ionEngBin) };

					//find initial particle that is closest to bin particle
					for (unsigned int part = 0; part < initialData.pitch.size(); part++)
					{ //don't know if I'm wild about iterating over every particle for every bin particle, but it's the most general and should work regardless of sim particle arrangement
						double tmpAngDiff{ std::abs(tmpBinAng - initialData.pitch.at(part)) };
						double tmpEngDiff{ std::abs(tmpBinEng - initialData.energy.at(part)) };

						if (tmpAngDiff <= angDiff && tmpEngDiff <= engDiff)
						{
							angDiff = tmpAngDiff;
							engDiff = tmpEngDiff;
							partInd = (int)part;
						}
						if (tmpAngDiff < FLT_EPSILON && tmpEngDiff < FLT_EPSILON)
						{
							part = initialData.pitch.size();
							std::cout << "Sys EPS" << std::endl;
						} //if both differences are less than FLT_EPSILON - ~1e-7 on Windows 10, it's close enough - there's no bin 1e-7eV wide, nor is there likely to ever be
					}

					pitches.at (ionAngBin * binEnergies.size() + ionEngBin) = satUpwardData.pitch.at(partInd);
					energies.at(ionAngBin * binEnergies.size() + ionEngBin) = satUpwardData.energy.at(partInd);
					weights.at (ionAngBin * binEnergies.size() + ionEngBin) = bsNumFluxBins.at(ionAngBin).at(ionEngBin);
				}
			}

			return binning::binWeighted(pitches, energies, binAngles, binEnergies, weights);
		}
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satEFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWt)
		{
			dblVec2D ret{ binning::binWeighted(sat.pitch, sat.energy, binAngles, binEnergies, numWt) };
			
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

			//convert to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng);

			return ret;
		}

		DLLEXP dblVec2D backEFlux(const ParticleData& initialData, const ParticleData& satData, const ParticleData& escapeData,
			const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge)
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

			dblVec2D escapeEPABins{ binning::binWeighted(escapeData.pitch, escapeData.energy, binAngles, binEnergies, maxwCounts) };
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle

			printVec2D(escapeEPABins, "Escaped Particles (dNflux) per E, PA bin");

			std::vector<double> escapeEBins(escapeEPABins.at(0).size()); //Sum of escaped particles at each energy, units of dNflux
			for (unsigned int egy = 0; egy < escapeEBins.size(); egy++)  //iterate over energies
			{
				for (auto& ang : escapeEPABins)                          //sum over angles, add to sum vector at energyBin
					escapeEBins.at(egy) += ang.at(egy);
				escapeEBins.at(egy) /= (double)escapeEPABins.size() / 2; //isotropically space across pitch angle bins - divide by # ion PA bins, this value later put in each ionospheric angle bin
			}
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins

			printVec2D({ escapeEBins }, "Escaped Particles (dNflux) per E bin");

			std::vector<double> bsdNFluxByEBin{ backscat::Nflux_bs(escapeEBins, binEnergies, 1.5, -4.0, -2.1, 0.3) }; //obtained by log linefitting Evans, 1974 - these seem closest
			// output: 1D vector of backscatter dNFlux(?) by energy bin at ionosphere
			
			printVec2D({ bsdNFluxByEBin }, "BS NFlux per E bin at Ionosphere");

			double dAngBin{ binAngles.at(1) - binAngles.at(0) };
			dblVec2D numFlux(binAngles.size(), std::vector<double>(binEnergies.size()));
			for (unsigned int ang = binAngles.size() / 2; ang < binAngles.size(); ang++)
			{
				double binSolidAngle{ cos((binAngles.at(ang) - 0.5 * dAngBin) * RADS_PER_DEG) - cos((binAngles.at(ang) + 0.5 * dAngBin) * RADS_PER_DEG) };
				for (unsigned int eng = 0; eng < binEnergies.size(); eng++)
					numFlux.at(ang).at(eng) = bsdNFluxByEBin.at(eng) * binSolidAngle;
			}
			// output: 2D vector of bs NFlux(?) at ionosphere per pitch, energy bin - should only be upward (90-180)

			printVec2D(numFlux, "BS NFlux at Ionosphere");

			dblVec2D ret{ steady::bsSrcToSat(numFlux, initialData, satData, binAngles, binEnergies) };
			// output: 2D vector of bs NFlux(?) at satellite per pitch, energy bin - should only be upward (90-180)

			printVec2D(ret, "satNFlux");

			//convert to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
			{
				double binSolidAngle{ cos((binAngles.at(ang) - 0.5 * dAngBin) * RADS_PER_DEG) - cos((binAngles.at(ang) + 0.5 * dAngBin) * RADS_PER_DEG) };
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng) / (binSolidAngle);
			}
			// output: 2D vector of bs dEflux(?)

			printVec2D(ret, "dEflux");

			exit(1);

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

		DLLEXP void countsToEFlux(dblVec2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection)
		{
			double dlogE{ log10(binEnergies.at(1)) - log10(binEnergies.at(0)) };

			for (unsigned int ang = 0; ang < binAngles.size(); ang++)
			{
				for (unsigned int eng = 0; eng < binEnergies.size(); eng++)
				{
					double v_perp{ sqrt(2 * binEnergies.at(eng) * JOULE_PER_EV / mass) };
					double A_gyro{ PI * pow((mass * v_perp / (charge * BatXSection)), 2) }; //don't need to absolute value because it's squared
					double binWidth{ pow(10, log10(binEnergies.at(eng)) + 0.5 * dlogE) - pow(10, log10(binEnergies.at(eng)) - 0.5 * dlogE) };
					//double solidAngle{ write this }; //involves pitch angle
					double solidAngle{ 1.0 };

					energyData.at(ang).at(eng) *= std::abs(cos(binAngles.at(ang) * RADS_PER_DEG)) / (A_gyro * solidAngle * binWidth);
				}
			}
		}

		DLLEXP void divBinsByCosPitch(dblVec2D& data, std::vector<double> binAnglesDegrees)
		{
			for (unsigned int iii = 0; iii < data.size(); iii++)
				for (unsigned int jjj = 0; jjj < data.at(0).size(); jjj++)
					data.at(iii).at(jjj) /= std::abs(cos(RADS_PER_DEG * binAnglesDegrees.at(iii)));
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
			BS_scnd(x) = 10 ^ (scnd_logm * logE + scnd_logb)
				=== d_sec * x ^ (scnd_logm), where d_sec = 10 ^ (scnd_logb)
			BS_prim(x) = 10 ^ (prim_logm * log(E / E_incident) + prim_logb) * (10000 / E_incident)
				=== d_pri * x ^ (prim_logm), where d_pri = 10 ^ (prim_logb + 4) / E_incident ^ (prim_logm + 1)

			Integral:
			BS'_scnd(x) = d_sec / (scnd_logm + 1) * x ^ (scnd_logm + 1)
			BS'_prim(x) = d_pri / (prim_logm + 1) * x ^ (prim_logm + 1), if x > E_incident, BS_scnd||prim = 0
		*/

		DLLEXP double F_flux(double evalE, double incidentE, double incidentCnt, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{ //describes a log linefit of the backscatter curves detailed in Evans, 1974
			return (pow(10, scnd_logb) * pow(evalE, scnd_logm) + //secondary BS +
				    pow(10, prim_logb + 4) / pow(incidentE, prim_logm + 1) * pow(evalE, prim_logm)) * //primary BS *
				    incidentCnt; //incident e count
		}

		DLLEXP double integralF_flux(double lower, double upper, double incidentE, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb)
		{
			if (upper > incidentE * (1 + FLT_EPSILON))
				throw std::invalid_argument("integralF_bs: upper limit of integration is higher than incident energy - function is zero above incident energy - upper limit: " + std::to_string(upper) + ", incidentE: " + std::to_string(incidentE));
			
			double integral_sec{ (pow(upper, scnd_logm + 1) - pow(lower, scnd_logm + 1)) * pow(10, scnd_logb) / (scnd_logm + 1) };
			double integral_prm{ (pow(upper, prim_logm + 1) - pow(lower, prim_logm + 1)) * pow(10, prim_logb + 4) / ((prim_logm + 1) * pow(incidentE, prim_logm + 1)) };
			return integral_sec + integral_prm;
		}

		DLLEXP std::vector<double> Nflux_bs(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb)
		{ //converts incident dNflux to bs Nflux
			double logEBinMin{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEBinMin };                             //depends on a logarithmically spaced E, won't work otherwise

			std::vector<double> bsNflux(binEnergies.size());
			for (unsigned int nFluxBin = 0; nFluxBin < bsNflux.size(); nFluxBin++)             //bins that contain the number flux of the backscatter in the energy bin of the same index
			{
				double engmin{ pow(10, (nFluxBin - 0.5) * dlogE + logEBinMin) };               //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
				double engmax{ pow(10, (nFluxBin + 0.5) * dlogE + logEBinMin) };

				for (unsigned int incEbin = nFluxBin; incEbin < bsNflux.size(); incEbin++)     //nFlux bins are contributed to by all incident particles with E higher than the E associated with the nFlux bin
				{
					double incidentE{ pow(10, log10(binEnergies.at(incEbin)) + 0.5 * dlogE) }; //incident E is upper limit of bin - particles in here can be up to this energy, so why not use this instead of mid bin?
					double intF{ integralF_flux(engmin, engmax, incidentE, primary_logm, primary_logb, secondary_logm, secondary_logb) };

					bsNflux.at(nFluxBin) += intF * binCounts.at(incEbin) / (engmax - engmin);
				}
			}

			return bsNflux;
		}
	} //end namespace postprocess::backscat
}
