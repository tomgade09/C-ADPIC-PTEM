#include "postprocess.h"
#include <algorithm>
#include <iterator>
#include "utils\numerical.h"

namespace postprocess
{
	//struct defines
	ParticleData::ParticleData(const std::vector<double>& v_para, const std::vector<double>& v_perp, double mass) //auto calculate E, Pitch
		{ utils::numerical::v2DtoEPitch(v_para, v_perp, mass, energy, pitch); } //vpara, vperp not even assigned - assumed they aren't needed outside of calculating E, pitch

	void ParticleData::free()
	{
		vpara.clear(); //clear data in vectors
		vperp.clear();
		energy.clear();
		pitch.clear();
		t_esc.clear();
		s_pos.clear();

		vpara.shrink_to_fit(); //reduce capacity to 0, freeing memory
		vperp.shrink_to_fit();
		energy.shrink_to_fit();
		pitch.shrink_to_fit();
		t_esc.shrink_to_fit();
		s_pos.shrink_to_fit();
	}
	
	MaxwellianData::MaxwellianData(double sion, double smag, double dlogE) : s_ion{ sion }, s_mag{ smag }, dlogE_dist{ dlogE } {}

	void MaxwellianData::push_back_ion(double peak, double magnitude) { ionEPeak.push_back(peak); iondEMag.push_back(magnitude); }
	void MaxwellianData::push_back_mag(double peak, double magnitude) { magEPeak.push_back(peak); magdEMag.push_back(magnitude); }
	
	PPData::PPData(double sion, double smag, std::vector<double> EBins, std::vector<double> PABins, MaxwellianData maxData, ParticleData init, ParticleData btm, ParticleData up, ParticleData dn) :
		s_ion{ sion }, s_mag{ smag }, energyBins{ EBins }, pitchBins{ PABins }, initial{ init }, bottom{ btm }, upward{ up }, dnward{ dn }
	{
		maxCounts = maxwellian::formCountsVector(initial, dnward, maxData);
	}
	//end struct defines


	dblVec2D steadyFlux(PPData ppdata)
	{
		dblVec2D simfluxupward{ steady::simEnergyFlux(ppdata.dnward, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts, -1.3597036e-5) }; //calculate binned flux for downward facing detector data, upward flux
		dblVec2D simfluxdnward{ steady::simEnergyFlux(ppdata.upward, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts, -1.3597036e-5) }; //calculate binned flux for upward facing detector data
		//dblVec2D backscatflux { steady::bsEnergyFlux (ppdata.initial, ppdata.dnward, ppdata.bottom, ppdata.pitchBins, ppdata.energyBins, ppdata.maxCounts, 9.10938356e-31, -1.6021766209e-19, -1.3597036e-5) }; //calculate backscatter flux
		
		for (int iii = 0; iii < simfluxupward.size(); iii++)
			for (int jjj = 0; jjj < simfluxupward.at(0).size(); jjj++)
				simfluxupward.at(iii).at(jjj) += simfluxdnward.at(iii).at(jjj) /*+ backscatflux.at(iii).at(jjj)*/;

		return simfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}

	namespace steady
	{
		dblVec2D simEnergyFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWt,
			double BatXSection)
		{
			dblVec2D ret{ maxwellian::countInBinsWeighted(sat.pitch, sat.energy, binAngles, binEnergies, numWt) };

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

		dblVec2D bsEnergyFlux(const dblVec2D& initialData, const dblVec2D& satData, const dblVec2D& escapeData, const std::vector<double>& binAngles,
			const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge, double BatXSection)
		{
			dblVec2D escCntBins{ maxwellian::countInBinsWeighted(escapeData.at(2), escapeData.at(3), binAngles, binEnergies, maxwCounts) };
			double Bion{ 5.73447893584459e-5 }; //B at ionosphere - 100km altitude, 101.322 s
			std::vector<double> escCntBinsSum(escCntBins.at(0).size()); //each bin in units of counts -> # particles
			for (int iii = 0; iii < escCntBinsSum.size(); iii++)
			{
				for (int jjj = 0; jjj < escCntBins.size(); jjj++)
					escCntBinsSum.at(iii) += escCntBins.at(jjj).at(iii);
				escCntBinsSum.at(iii) /= (Bion * escCntBins.size()); //isotropically space across pitch angle bins
			}

			std::vector<double> numFluxByBin{ backscat::sumIntegralsOfNumFluxFcnsPerBin(escCntBinsSum, binEnergies, 1.5, -4.0, -2.1, 0.3) }; //obtained by log linefitting Evans, 1974 - these seem closest

			double dAngBin{ binAngles.at(1) - binAngles.at(0) };
			dblVec2D numFlux(binAngles.size(), std::vector<double>(binEnergies.size()));
			for (int ang = 0; ang < binAngles.size(); ang++)
				for (int eng = 0; eng < binEnergies.size(); eng++)
					numFlux.at(ang).at(eng) = numFluxByBin.at(eng) * Bion * (cos((binAngles.at(ang) - 0.5 * dAngBin) * RADS_PER_DEG) - cos((binAngles.at(ang) + 0.5 * dAngBin) * RADS_PER_DEG));

			dblVec2D ret{ backscat::matchIonBSToSatAndCount(numFlux, initialData, satData, binAngles, binEnergies) };
			
			//convert to dEflux
			double Bsat{ 1.35970360638581e-5 }; //B at 4e6 altitude, 4.0713e6 s (along field line from Re)
			for (int ang = 0; ang < ret.size(); ang++)
			{
				double binSolidAngle{ cos((binAngles.at(ang) - 0.5 * dAngBin) * RADS_PER_DEG) - cos((binAngles.at(ang) + 0.5 * dAngBin) * RADS_PER_DEG) };
				for (int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= binEnergies.at(eng) / (Bsat * binSolidAngle);
			}

			return ret;
		}
	}

	namespace maxwellian
	{
		double counts(double E, double binWidth_E, double sigma_kT, double dEflux_kT, double binWidth_kT) //or maybe maxwellian::getWeights
		{
			if (E == 0.0) { return 0.0; }
			return (binWidth_E * exp(-E / sigma_kT)) / E * (dEflux_kT / (exp(-1.0) * binWidth_kT));
		}

		std::vector<double> formCountsVector(const ParticleData& init, const ParticleData& dnward, const MaxwellianData& maxData)
		{
			std::vector<double> max(init.energy.size());
			
			double s_ion{ maxData.s_ion };
			double s_mag{ maxData.s_mag };
			double dlogE_dist{ maxData.dlogE_dist };

			auto maxwrapper = [&dlogE_dist, &maxData](double E, double Epeak, double PkMag, bool zero)
			{
				if (zero)
					return 0.0;
				else
					return maxwellian::counts(E,
						pow(10, log10(E) + 0.5 * dlogE_dist) - pow(10, log10(E) - 0.5 * dlogE_dist), //bin width at E, according to specified dlogE
						Epeak, PkMag,
						pow(10, log10(E) + 0.5 * 4.0 / 47.0) - pow(10, log10(E) - 0.5 * 4.0 / 47.0)); }; //bin width at kT, according to satellite bins

			std::vector<double> maxtmp;
			for (int entr = 0; entr < maxData.ionEPeak.size(); entr++) //iterate over ionospheric maxwellian specifications
			{
				auto ioncnt = [&](double E, double s) { return maxwrapper(E, maxData.ionEPeak.at(entr), maxData.iondEMag.at(entr), (s > s_ion * 1.001)); };

				std::transform(init.energy.begin(), init.energy.end(), init.s_pos.begin(), std::back_inserter(maxtmp), ioncnt); //generate maxwellian count if ionospheric particle
				std::transform(maxtmp.begin(), maxtmp.end(), max.begin(), max.begin(), [](double x, double y) { return x + y; }); //add final vector and the tmp vector together
				maxtmp.clear();
			}

			for (int entr = 0; entr < maxData.magEPeak.size(); entr++) //iterate over magnetospheric maxwellian specifications
			{
				auto maxcnt = [&](double E, double s) { return maxwrapper(E, maxData.magEPeak.at(entr), maxData.magdEMag.at(entr), (s < s_mag * 0.999)); };

				std::transform(init.energy.begin(), init.energy.end(), init.s_pos.begin(), std::back_inserter(maxtmp), maxcnt);
				std::transform(maxtmp.begin(), maxtmp.end(), max.begin(), max.begin(), [](double x, double y) { return x + y; });
				maxtmp.clear();
			}

			for (int iii = 0; iii < init.pitch.size(); iii++) //isotropize counts -> 3D
			{
				if (init.s_pos.at(iii) < maxData.s_ion) //ionospheric source
					max.at(iii) *= -cos(init.pitch.at(iii) * RADS_PER_DEG);
				else //magnetospheric source
					max.at(iii) *= 1.0 / cos(dnward.pitch.at(iii) * RADS_PER_DEG); //satup...is this right?
			}

			return max;
		}

		dblVec2D countInBinsWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts)
		{
			dblVec2D ret(binAngles.size(), std::vector<double>(binEnergies.size()));

			double logMinEBinMid{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logMinEBinMid };
			double dangle{ binAngles.at(1) - binAngles.at(0) };
			
			int endedEarly{ 0 };
			int endedLater{ 0 };

			for (int part = 0; part < particleEnergies.size(); part++)
			{
				double partEnerg{ particleEnergies.at(part) };
				double partPitch{ particlePitches.at(part) };
				
				{
					int angbin{ (int)(partPitch / dangle) }; //this should give the bin index
					int engbin{ (int)((log10(partEnerg) - logMinEBinMid) / dlogE) }; //ditto

					if (partEnerg >= pow(10, (engbin - 0.5) * dlogE + logMinEBinMid) &&
						partEnerg <  pow(10, (engbin + 0.5) * dlogE + logMinEBinMid) &&
						partPitch >= engbin * dangle &&
						partPitch <  (engbin + 1) * dangle)
					{
						endedEarly++;
						ret.at(angbin).at(engbin) += maxwCounts.at(part);
						continue;
					}
					endedLater++;
				}
				
				for (int angbin = 0; angbin < binAngles.size(); angbin++)
				{
					double angmin{ angbin *       dangle }; //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
					double angmax{ (angbin + 1) * dangle };
					if (angbin == (binAngles.size() - 1)) { angmax += 0.001; } //I want the top angle bin to include 180 and right now the conditional is ... < 180, not ...<= 180

					for (int ebin = 0; ebin < binEnergies.size(); ebin++)
					{
						double emin{ pow(10, (ebin - 0.5) * dlogE + logMinEBinMid) }; //assumes evenly spaced angle bins - a reasonable assumption and probably commonly the case
						double emax{ pow(10, (ebin + 0.5) * dlogE + logMinEBinMid) };

						if ((partPitch >= angmin) && (partPitch < angmax) &&
							(partEnerg >= emin)   && (partEnerg < emax))
						{
							ret.at(angbin).at(ebin) += maxwCounts.at(part); //maxwCounts needs to be as long as particleEnergies/particlePitches
							angbin = binAngles.size(); //no need to iterate over other bins once particle bin is found
							ebin = binEnergies.size();
						}
					}
				}
			}

			std::cout << "Count in bins weighted: ended early: " << endedEarly << ", ended later: " << endedLater << std::endl;

			return ret;
		}

		void countsToEFlux(dblVec2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection)
		{
			double dlogE{ log10(binEnergies.at(1)) - log10(binEnergies.at(0)) };

			for (int ang = 0; ang < binAngles.size(); ang++)
			{
				for (int eng = 0; eng < binEnergies.size(); eng++)
				{
					double v_perp{ sqrt(2 * binEnergies.at(eng) * JOULE_PER_EV / mass) };
					double A_gyro{ PI * pow((mass * v_perp / (charge * BatXSection)), 2) }; //don't need to absolute value because it's squared
					double binWidth{ pow(10, log10(binEnergies.at(eng)) + 0.5 * dlogE) - pow(10, log10(binEnergies.at(eng)) - 0.5 * dlogE) };
					//double solidAngle{ write this }; //involves pitch angle
					double solidAngle{ 1.0 };

					energyData.at(ang).at(eng) *= abs(cos(binAngles.at(ang) * RADS_PER_DEG)) / (A_gyro * solidAngle * binWidth);
				}
			}
		}

		void divBinsByCosPitch(dblVec2D& data, std::vector<double> binAnglesDegrees)
		{
			for (int iii = 0; iii < data.size(); iii++)
				for (int jjj = 0; jjj < data.at(0).size(); jjj++)
					data.at(iii).at(jjj) /= abs(cos(RADS_PER_DEG * binAnglesDegrees.at(iii)));
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

		dblVec2D matchIonBSToSatAndCount(const dblVec2D& bsNumFluxBins, const dblVec2D& initialData, const dblVec2D& satDownData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies)
		{
			dblVec2D ret(binAngles.size());
			for (int iii = 0; iii < ret.size(); iii++)
				ret.at(iii) = std::vector<double>(binEnergies.size());

			double logEMinBinMid{ log10(binEnergies.at(0)) };
			double dlogE{ log10(binEnergies.at(1)) - logEMinBinMid };
			double dangle{ binAngles.at(1) - binAngles.at(0) };

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
}