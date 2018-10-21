#include "utils/postprocess.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

using std::pow;
using std::string;
using utils::numerical::generateSpacedValues;

constexpr double EVANS_PRIM_LOGM{ 1.5 };  //obtained by log linefitting Evans, 1974 - these seem closest
constexpr double EVANS_PRIM_LOGB{ -4.0 };
constexpr double EVANS_SECD_LOGM{ -2.1 };
constexpr double EVANS_SECD_LOGB{ 0.3 };

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(const PPData& ppdata)
	{
		auto printVec2D = [](const dblVec2D& prt, string name)
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

		//Add some guards to make sure ppdata is properly formed

		// 1. Adjust Maxwellian by Cos, SQRT(B) Factors
		double Aratio_ion_sat{ std::sqrt(ppdata.B_ion / ppdata.B_sat) };
		double Aratio_mag_sat{ std::sqrt(ppdata.B_sat / ppdata.B_mag) };
		double Aratio_ion_bs { std::sqrt(ppdata.B_ion / ppdata.B_ion) * std::sqrt(ppdata.B_sat / ppdata.B_ion) };
		double Aratio_mag_bs { std::sqrt(ppdata.B_ion / ppdata.B_mag) * std::sqrt(ppdata.B_sat / ppdata.B_ion) };
		double bsScale{ 0.1 }; //for now, arbitrary factor to get in the ballpark

		dblVec weights_atSat{ ppdata.maxWeights }; //scaled by decrease in cross-sectional area A
		dblVec weights_atIon{ ppdata.maxWeights }; //as the particle moves down the B field line

		for (unsigned int iii = 0; iii < ppdata.initial.s_pos.size(); iii++) //isotropize counts -> 3D
		{
			if (ppdata.initial.s_pos.at(iii) < ppdata.s_ion * 1.001)     //ionospheric source
			{
				weights_atSat.at(iii) *= -cos(ppdata.initial.pitch.at(iii) * RADS_PER_DEG) * Aratio_ion_sat;
				weights_atIon.at(iii) *= -cos(ppdata.initial.pitch.at(iii) * RADS_PER_DEG) * Aratio_ion_bs * bsScale; //without QSPS, there shouldn't be any ionospheric-source particles influencing the backscatter
			}
			else if (ppdata.initial.s_pos.at(iii) > ppdata.s_mag * 0.999)//magnetospheric source
			{
				weights_atSat.at(iii) *= 1.0 / cos(ppdata.dnward.pitch.at(iii) * RADS_PER_DEG) * Aratio_mag_sat;
				weights_atIon.at(iii) *= 1.0 / cos(ppdata.bottom.pitch.at(iii) * RADS_PER_DEG) * Aratio_mag_bs * bsScale;
			}
			else
				throw std::logic_error("postprocess::steadyFlux : particle is not ionospheric or magnetospheric source");
		}

		
		// 2. Calculate dEfluxes
		dblVec2D distfluxdnward{ EFlux::satdEFlux(ppdata.dnward, ppdata.satbins, weights_atSat) };
		dblVec2D distfluxupward{ EFlux::satdEFlux(ppdata.upward, ppdata.satbins, weights_atSat) };
		dblVec2D backfluxupward{ EFlux::bksdEFlux(ppdata, weights_atIon) };

		printVec2D(distfluxdnward, "Dnward Flux at Satellite");
		printVec2D(distfluxupward, "Upward Flux at Satellite");
		printVec2D(backfluxupward, "Upward Backscatter Flux at Satellite");

		// 3. Sum dEfluxes
		for (unsigned int iii = 0; iii < distfluxupward.size(); iii++)
			for (unsigned int jjj = 0; jjj < distfluxupward.at(iii).size(); jjj++)
				distfluxupward.at(iii).at(jjj) += distfluxdnward.at(iii).at(jjj) + backfluxupward.at(iii).at(jjj);

		return distfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const Bins& dist, const Bins& sat, const dblVec2D& bsCounts, const ParticleData& initialData, const ParticleData& satUpwardData)
		{
			const double ANGMAXERR{ 0.1 }; //max error in degrees
			const double ENYMAXERR{ 0.1 }; //max error in eV
			auto err = [](double base, double diff) { return std::abs((base - diff) / base); };

			const double logEMinBinMid{ log10(dist.E.at(0)) };
			const double dlogE{ log10(dist.E.at(1)) - logEMinBinMid };
			const double dangle{ dist.PA.at(1) - dist.PA.at(0) };

			ParticleData particles;
			dblVec       counts;

			for (unsigned int ang = 0; ang < dist.PA.size(); ang++)
			{
				for (unsigned int eny = 0; eny < dist.E.size(); eny++) //this works because ionospheric bins are same as distribution
				{
					particles.pitch.push_back (satUpwardData.pitch.at (ang * dist.E.size() + eny));
					particles.energy.push_back(satUpwardData.energy.at(ang * dist.E.size() + eny));
					counts.push_back(bsCounts.at(ang).at(eny));
				}
			}

			return binning::binWeighted(particles, sat, counts);
		}
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satdEFlux(const ParticleData& sat, const Bins& satBins, const dblVec& weights_atSat)
		{
			// Section 1 - Get dNflux at Satellite
			// 1. Bin Particles by Satellite Detector PA, Energy Bins
			dblVec2D ret{ binning::binWeighted(sat, satBins, weights_atSat) };
			// Section 1 End

			// 2. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= satBins.E.at(eng);


			return ret;
		}

		DLLEXP dblVec2D bksdEFlux(const PPData& ppdata, const dblVec& weights_atIon)
		{
			// Section 1 - Get dNflux at Satellite
			// 1.1. Bin Escaped Particles by EOM simulation (distribution) PA, Energy Bins (for high res backscatter)
			dblVec2D escapeCountBinned{ binning::binWeighted(ppdata.bottom, ppdata.distbins, weights_atIon) };
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle


			// 1.2. Calculate BS dNflux from dNflux Incident to Ionosphere
			dblVec2D dNflux_BS{ backscat::dNflux_bs_ion(ppdata.distbins, escapeCountBinned) };
			// output: 2D vector of backscatter dNFlux by dist bins at ionosphere (upgoing)


			// 1.3. Translate BS dNflux at Ionosphere (dist binned) to dNflux at Satellite (sat binned)
			dblVec2D ret{ steady::bsSrcToSat(ppdata.distbins, ppdata.satbins, dNflux_BS, ppdata.initial, ppdata.upward) };
			// output: 2D vector of bs dNFlux at satellite per PA, E (sat binned) - should only be upward (90-180)
			// Section 1 End


			// 2.1. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
			{
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= ppdata.satbins.E.at(eng);
			}
			// output: 2D vector of bs dEflux


			return ret;
		}
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const ParticleData& particles, const Bins& bins, const dblVec& counts)
		{
			dblVec2D ret(bins.PA.size(), dblVec(bins.E.size()));

			bool Eascending{ (bins.E.back() > bins.E.front()) };   //determines whether or not bin E's ascend from less E to more E as ind increases
			bool Aascending{ (bins.PA.back() > bins.PA.front()) }; //same for angle

			double Emax{ Eascending ? bins.E.back() : bins.E.front() };
			double Emin{ !Eascending ? bins.E.back() : bins.E.front() };
			double Amax{ Aascending ? bins.PA.back() : bins.PA.front() };
			double Amin{ !Aascending ? bins.PA.back() : bins.PA.front() };

			double logMinEBinMid{ log10(Emin) };
			double dlogE_bin{ std::abs(log10(bins.E.at(1)) - log10(bins.E.at(0))) };
			double dangle_bin{ std::abs(bins.PA.at(1) - bins.PA.at(0)) };

			int outsideLimits{ 0 };

			for (unsigned int part = 0; part < particles.energy.size(); part++) //iterate over particles
			{
				double partEnerg{ particles.energy.at(part) };
				double partPitch{ particles.pitch.at(part) };

				if ((partEnerg == 0.0 && partPitch == 0.0) || counts.at(part) == 0.0) continue; //guards - if particle E, PA is zero, it wasnt detected - just skip it
				if (partEnerg > pow(10, log10(Emax) + (0.5 * dlogE_bin)) ||                     //if particle is outside E, PA measured limits - skip it
					partEnerg < pow(10, log10(Emin) - (0.5 * dlogE_bin)) ||                     //PA shouldn't be an issue, there just in case
					partPitch > Amax + (0.5 * dangle_bin) ||                                    //also, code counts how many are outside limits for diagnostics
					partPitch < Amin - (0.5 * dangle_bin))
				{ outsideLimits++; continue; }

				//calculate bin index for E, PA of particle
				int angbin{ (int)(partPitch / dangle_bin) }; //this should give the bin index
				int engbin{ (int)((log10(partEnerg) - (logMinEBinMid - 0.5 * dlogE_bin)) / dlogE_bin) }; //ditto
				if (!Eascending) engbin = bins.E.size()  - 1 - engbin; //reverses the bin index if E is not ascending
				if (!Aascending) angbin = bins.PA.size() - 1 - angbin; //ditto

				if (partEnerg >= pow(10, log10(bins.E.at(engbin)) - 0.5 * dlogE_bin) &&
					partEnerg <  pow(10, log10(bins.E.at(engbin)) + 0.5 * dlogE_bin) &&
					partPitch >= bins.PA.at(angbin) - 0.5 * dangle_bin &&
					partPitch <  bins.PA.at(angbin) + 0.5 * dangle_bin)
				{
					ret.at(angbin).at(engbin) += counts.at(part);
					continue;
				}
				else //this shouldn't ever execute, guards should prevent zero and out of limits values
					throw std::logic_error(std::to_string(part) + " >> " + std::to_string(particles.energy.at(part)) + ":" + std::to_string(engbin) + ",  " + std::to_string(particles.pitch.at(part)) + ":" + std::to_string(angbin));
			}

			if (outsideLimits > 0) std::cout << "postprocess::binning::binWeighted : Particles out of limits: " << outsideLimits << "\n";

			return ret;
		}

		DLLEXP void symmetricBins0To360(dblVec2D& data, dblVec& binAngles) //takes bins from 0-180 and extends to 360 symmetrically reflected across 180 degrees
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

		DLLEXP double evans_flux(double E_eval, double E_incident)
		{// Describes a log linefit of the backscatter curves detailed in Evans, 1974
		 // E_incident defines the primary BS curve
		 // E_eval sets where we want to know the value of the function
			return (pow(10.0, EVANS_SECD_LOGB) * pow(E_eval, EVANS_SECD_LOGM) + //secondary BS +
				    pow(10.0, EVANS_PRIM_LOGB + 4.0) / pow(E_incident, EVANS_PRIM_LOGM + 1.0) * pow(E_eval, EVANS_PRIM_LOGM)); //primary BS
		}

		DLLEXP double integralEvans_flux(double lower, double upper, double E_incident)
		{
			//if (upper > E_incident * (1 + FLT_EPSILON))
				//throw std::invalid_argument("integralF_bs: upper limit of integration is higher than incident energy - function is zero above incident energy - upper limit: " + std::to_string(upper) + ", incidentE: " + std::to_string(E_incident));
			
			double integral_sec{ (pow(upper, EVANS_SECD_LOGM + 1.0) - pow(lower, EVANS_SECD_LOGM + 1.0)) * pow(10.0, EVANS_SECD_LOGB) / (EVANS_SECD_LOGM + 1.0) };
			double integral_prm{ (pow(upper, EVANS_PRIM_LOGM + 1.0) - pow(lower, EVANS_PRIM_LOGM + 1.0)) * pow(10.0, EVANS_PRIM_LOGB + 4.0) / ((EVANS_PRIM_LOGM + 1.0) * pow(E_incident, EVANS_PRIM_LOGM + 1.0)) };
			return integral_sec + integral_prm;
		}

		DLLEXP dblVec2D dNflux_bs_ion(const Bins& dist, const dblVec2D& escapeCountBinned)
		{ //converts downward dNflux at ionosphere (dist binned) to bs (upward) dNflux (also dist binned)
			// 1. Sum dNflux over PA Bins (dist bins), Per E Bin and Average
			dblVec escapeCountPerE(dist.E.size());              //Sum of escaped particles at each energy, units of dNflux
			for (unsigned int egy = 0; egy < escapeCountPerE.size(); egy++) //iterate over energies
			{
				for (auto& ang : escapeCountBinned)                                //sum over angles, add to sum vector at energyBin
					escapeCountPerE.at(egy) += ang.at(egy);
				escapeCountPerE.at(egy) /= (double)escapeCountBinned.size() / 2;   //isotropically space across pitch angle bins - divide by # ion PA bins, this value later put in each ionospheric angle bin
			}
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins


			// 2. Calculate upward dNflux (backscatter) per E bin
			double logEBinMin{ log10(dist.E.at(0)) };         //depends on an array where E is minimum at index 0, max at last index
			double dlogE{ log10(dist.E.at(1)) - logEBinMin }; //depends on a logarithmically spaced E, won't work otherwise

			dblVec upwardCountPerE{ escapeCountPerE };
			for (unsigned int ebin = 0; ebin < escapeCountPerE.size(); ebin++) //why is this the case!!!
				upwardCountPerE.at(ebin) *= dist.E.at(ebin);  //convert to dEflux escaping into the layer

			dblVec dNfluxPerE_bs(dist.E.size());
			for (unsigned int dNFluxBin = 0; dNFluxBin < dNfluxPerE_bs.size(); dNFluxBin++)       //bins that contain the number flux of the backscatter in the energy bin of the same index
			{
				double engmin{ pow(10, (dNFluxBin - 0.5) * dlogE + logEBinMin) };                 //assumes evenly log spaced angle bins
				double engmax{ pow(10, (dNFluxBin + 0.5) * dlogE + logEBinMin) };

				for (unsigned int incEbin = dNFluxBin; incEbin < dNfluxPerE_bs.size(); incEbin++) //nFlux bins are contributed to by all incident particles with E higher than the E associated with the nFlux bin
				{//change to mid bin
					double incidentE{ dist.E.at(incEbin) };                                       //incident E is upper limit of bin
					double intF{ integralEvans_flux(engmin, engmax, incidentE) };

					dNfluxPerE_bs.at(dNFluxBin) += intF * upwardCountPerE.at(incEbin) / (engmax - engmin);
					//dNfluxPerE_bs.at(dNFluxBin) += john_flux(binEnergies.at(dNFluxBin), binEnergies.at(incEbin), primary_logm, primary_logb, secondary_logm, secondary_logb) * binCntxE.at(incEbin);
				}
			}
			// output: 1D vector of the upgoing (backscatter) dNflux per E


			// 3. Distribute BS dNflux Isotropically Over Pitch Bins
			dblVec2D dNfluxPerEPA_bs(dist.PA.size());
			for (unsigned int ang = 0; ang < dNfluxPerEPA_bs.size(); ang++)
			{
				if (dist.PA.at(ang) <= 90.0)
					dNfluxPerEPA_bs.at(ang) = dblVec(dist.E.size()); //empty vector of the right size
				else
				{
					dNfluxPerEPA_bs.at(ang) = dNfluxPerE_bs;
					for (unsigned int eny = 0; eny < dNfluxPerEPA_bs.at(ang).size(); eny++)
					{
						dNfluxPerEPA_bs.at(ang).at(eny) *= -cos(dist.PA.at(ang) * RADS_PER_DEG);
					}
				}
			}
			// output: 2D vector of bs dNFlux at ionosphere per pitch, energy bin - should only be upward (90-180)


			return dNfluxPerEPA_bs;
		}
	} //end namespace postprocess::backscat

	/*namespace multLevelBS
	{
		//DLLEXP dblVec2D scatterMain(const dblVec& bins_PA, const dblVec& bins_E, const dblVec2D& escapeBins, const dblVec& alt_levels, const dblVec& Z_levels, const dblVec& p_levels, double B_ionTop)
		DLLEXP dblVec2D scatterMain(const Ionosphere& ionsph, const Bins& escape)
		{
			//if (alt_levels.size() != Z_levels.size()) throw std::logic_error("postprocess::multLevelBS::scatterMain: altitude and Z vectors are not the same size.  Each altitude needs to have an associated Z and vice versa.");

			dblVec2D sumCollidePercent(escape.PA.size(), dblVec(escape.E.size()));
			dblVec2D alt_refl{ alt_reflect(escape, B, ionsph.B.at(0), t) };
			dblVec2D bsSumAtIonsphTop(escape.PA.size(), dblVec(escape.E.size()));

			//loop over levels of ionosphere
			for (unsigned int level; level < ionsph.s.size(); level++)
			{
				// >> level calculate
				dblVec2D bs_level{ bsAtLevel(ionsph, escape, alt_refl, sumCollidePercent, level) };

				// >> >> adjust PA and add to 620km bins

			}

		}

		DLLEXP dblVec2D bsAtLevel(const Ionosphere& ionsph, const Bins& ppbins, const dblVec2D& alt_refl, dblVec2D& sumCollideAbove, unsigned int level)
		{
			ParticleData downgoing;
			ParticleData upgoing;
			dblVec downcount;
			dblVec upcount;
			vector<double*> sumCollideAbove1D;

			// >> >> adjust PA, check if reflected
			for (unsigned int ang = 0; ang < ppbins.PA.size(); ang++)
			{
				if (level == 0) continue;

				double sinpa{ ionsph.B.at(level) / ionsph.B.at(0) * cos(ppbins.PA.at(ang) * RADS_PER_DEG) };
				
				if (sinpa > 1) continue; //particle reflects before this level

				if (ionsph.B.at(level + 1) / ionsph.B.at(0) * cos(ppbins.PA.at(ang) * RADS_PER_DEG) > 1)
				{ //particle reflects before next level
					for (unsigned int eny = 0; eny < ppbins.E.size(); eny++)
					{
						upgoing.energy.push_back(ppbins.E.at(eny));
						upgoing.pitch.push_back(180.0 - ppbins.PA.at(ang));
						upcount.push_back(ppbins.counts.at(ang).at(eny));
					}
				}

				for (unsigned int eny = 0; eny < ppbins.E.size(); eny++)
				{
					downgoing.energy.push_back(ppbins.E.at(eny));
					downgoing.pitch.push_back(asin(sinpa) / RADS_PER_DEG);
					downcount.push_back(ppbins.counts.at(ang).at(eny));
					sumCollideAbove1D.push_back(&sumCollideAbove.at(ang).at(eny));
				}
			}

			// >> >> calculate scattered
			for (unsigned int part = 0; part < downgoing.energy.size(); part++)
			{
				double sct{ scatterPct(*sumCollideAbove1D.at(part), ionsph.Z.at(level), ionsph.p.at(level), ionsph.h.at(level), levelbins.PA.at(ang), levelbins.E.at(eny)) };
				*sumCollideAbove1D.at(part) += sct;
				downcount.at(part) *= sct;
			}

			// >> >> bin particles at level
			dblVec2D downBinned{ binning::binWeighted(downgoing, ppbins, downcount) };
			dblVec2D upBinned{ binning::binWeighted(upgoing, ppbins, upcount) };

			// >> >> calculate backscatter
			Bins bsbins(ppbins.PA, ppbins.E, downBinned);
			dblVec2D ret{ backscat::dNflux_bs_ion(bsbins) };
			
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eny = 0; eny < ret.at(0).size(); eny++)
					ret.at(ang).at(eny) += upBinned.at(ang).at(eny);

			return ret;
		}

		DLLEXP dblVec2D alt_reflect(const Bins& bins, BField* B, double B_ion, double t)
		{
			auto sAtB_bracket = [&](double B_target, double err = 1.0e-10)
			{
				//setup variables for execution
				double s_guess{ 0.5 * RADIUS_EARTH };          //some initial guess - 0.5 RE will be in the sim
				double B_guess{ B->getBFieldAtS(s_guess, t) }; //B at s_guess
				double delta_s{ 1.0e6 };                       //initial value that is subtracted/added to s_guess
				unsigned int count{ 0 };                       //count number of loop iterations to prevent inf loops

				bool over{ 0 };

				while (abs((B_target - B_guess) / B_target) > err) //loop while B_guess is not close enough to B_target
				{
					//first: check where B_guess is in relation to B_target
					over = (B_guess > B_target);

					//second: subtract (if over) or add (if not) delta_s to s_guess until B_guess is less than/greater than (respectively) B_target
					while (1) //loop until break
					{
						if (over)
						{
							s_guess -= delta_s;
							if (s_guess < 0) { s_guess += delta_s; break; } //guard in case s_guess - delta_s is less than 0, which is unphysical
							B_guess = B->getBFieldAtS(s_guess, t);

							if (B_guess <= B_target) break; //if over, subtract delta_s until under
						}
						else
						{
							s_guess += delta_s;
							B_guess = B->getBFieldAtS(s_guess, t);

							if (B_guess >= B_target) break; //if under, add delta_s until over
						}
						count++;
						if (count > 1000000) throw std::logic_error("postprocess::multiLevelBS::alt_reflect: loop count is over a million - that shouldn't happen: delta_s, s_guess, B_guess: " +
							std::to_string(delta_s) + ", " + std::to_string(s_guess) + ", " + std::to_string(B_guess));
					}
					
					//third: shrink delta_s, reset loop count
					delta_s /= 5.0;
					count = 0;

					//guard in case delta_s is ever less than err - small changes in s will result in very small B changes
					//if delta_s is as low as err, B is really not changing much - trying to prevent infinite loops
					if (delta_s < err) throw std::logic_error("postprocess::multiLevelBS::alt_reflect: delta_s is less than error in bracketing lambda: B_target, delta_s, err"
						+ std::to_string(B_target) + ", " + std::to_string(delta_s) + ", " + std::to_string(err));
				}

				//return the result
				return s_guess;
			};

			dblVec2D s_ref(bins.PA.size(), dblVec(bins.E.size()));

			for (auto& ang : bins.PA)
			{
				double s_ang{ sAtB_bracket(B_ion / sin(ang * RADS_PER_DEG)) };

				for (auto& eny : bins.E)
				{
					s_ref.at(ang).at(eny) = s_ang;
				}
			}

			return s_ref;
		}

		DLLEXP double scatterPct(double sumCollideAbove, double Z, double p, double h, double E, double PA)
		{
			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	}*/
}
