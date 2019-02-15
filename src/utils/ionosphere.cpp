#include "utils/ionosphere.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

#include "utils/fileIO.h"

using std::pow;
using std::cout;
using std::string;
using std::to_string;
using std::logic_error;
using utils::numerical::generateSpacedValues;

constexpr double EVANS_PRIM_LOGM{ 1.5 };  //obtained by log linefitting Evans, 1974 - these seem closest
constexpr double EVANS_PRIM_LOGB{ -4.0 };
constexpr double EVANS_SECD_LOGM{ -2.1 };
constexpr double EVANS_SECD_LOGB{ 0.3 };

inline void printIonosphere(const ionosphere::IonosphereSpecs& ionsph)
{
	cout << "==================== Backscatter Simulation ====================" << "\n";
	cout << "Min, Max s (m): " << *(ionsph.s.end() - 1) << ", " << ionsph.s.front() << "\n";
	cout << "Min, Max B (T): " << *(ionsph.B.end() - 1) << ", " << ionsph.B.front() << "\n";
	cout << "Num of Layers : " << ionsph.s.size() - 1 << "\n";
	cout << "Atomic Species: " << "Fix later.\n";
	cout << "================================================================" << "\n";
}

inline void printLayer(const ionosphere::IonosphereSpecs& ionsph, unsigned int layer)
{
	cout << "Layer: " << layer << " / " << ionsph.s.size() - 2 << ", s: " << ionsph.s.at(layer) << ", B: " << ionsph.B.at(layer) << "\n";
}

namespace ionosphere
{
	DLLEXP dEflux_v2D steadyFlux(const EOMSimData& eomdata)
	{
		printIonosphere(eomdata.ionsph);

		// 1. Adjust Maxwellian by Cos, SQRT(B ratio) Factors
		double Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph to satellite
		double Aratio_mag_sat{ std::sqrt(eomdata.B_sat / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to satellite
		double Aratio_ion_ion{ std::sqrt(eomdata.B_ion / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph (up, reflect, down) to ionsph
		double Aratio_mag_ion{ std::sqrt(eomdata.B_ion / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to ionsph

		dNflux_v1D maxwellian_sat{ eomdata.maxwellian }; //scaled by decrease in gyroradius cross-sectional area A below
		dNflux_v1D maxwellian_ion{ eomdata.maxwellian };

		for (unsigned int part = 0; part < eomdata.initial.s_pos.size(); part++) //isotropize counts -> 3D
		{
			if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)     //ionospheric source, upgoing
			{
				maxwellian_sat.at(part) *= -cos(eomdata.initial.pitch.at(part) * RADS_PER_DEG) * Aratio_ion_sat;
				maxwellian_ion.at(part) *= Aratio_ion_ion; //should just be equal to 1, but gives it symmetry
			}
			else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)//magnetospheric source, downgoing
			{
				maxwellian_sat.at(part) *= 1.0 / cos(eomdata.dnward.pitch.at(part) * RADS_PER_DEG) * Aratio_mag_sat;
				maxwellian_ion.at(part) *= Aratio_mag_ion;
			}
			else
				throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source");
		}

		TESTVEC_NOTNEGWHOLEVEC(maxwellian_sat, "steadyFlux::maxwellian_sat");

		// 2. Calculate dEfluxes
		dEflux_v2D distfluxdnward{ dEFlux::satellite(eomdata.dnward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D distfluxupward{ dEFlux::satellite(eomdata.upward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D backfluxupward{ dEFlux::backscatr(eomdata, maxwellian_ion) };

		printVec2D(distfluxdnward, "Dnward Flux at Satellite");
		printVec2D(distfluxupward, "Upward Flux at Satellite");
		printVec2D(backfluxupward, "Flux Due to Backscatter at Satellite");


		// 3. Sum dEfluxes
		for (unsigned int ang = 0; ang < distfluxupward.size(); ang++)
			for (unsigned int eny = 0; eny < distfluxupward.at(ang).size(); eny++)
				distfluxupward.at(ang).at(eny) += distfluxdnward.at(ang).at(eny) + backfluxupward.at(ang).at(eny);

		return distfluxupward; //really, instead of just upward data, this is the total (see the nested loop above)
	}


	namespace dEFlux
	{
		DLLEXP dEflux_v2D satellite(const ParticleData& particles, const Bins& satBins, const dNflux_v1D& dNatSat)
		{
			// Section 1 - Get dNflux at Satellite
			// 1. Bin Particles by Satellite Detector PA, Energy Bins
			dNflux_v2D ret{ binning::binParticles(particles, satBins, dNatSat) };
			// Section 1 End

			// 2. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= satBins.E.at(eng);


			return ret;
		}

		DLLEXP dEflux_v2D backscatr(const EOMSimData& eomdata, const dNflux_v1D& dNatIonsph)
		{
			// Section 1 - Get dNflux at Satellite
			// 1.1. Bin Escaped Particles by EOM simulation (distribution) PA, Energy Bins (for high res backscatter)
			dNflux_v2D downwardAtIonsph{ binning::binParticles(eomdata.bottom, eomdata.distbins, dNatIonsph) };
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle

			TESTVEC_ISZEROFRSTHALF(downwardAtIonsph, "backscatr::downwardAtIonsph");

			// 1.2. Calculate BS dNflux from dNflux Incident to Ionosphere
			dNflux_v2D BSatIonsph{ multiLevelBS::scatterMain(eomdata, downwardAtIonsph) }; //new multi-level hotness
			// output: 2D vector of backscatter dNFlux by dist bins at ionosphere (upgoing)

			TESTVEC_ISZEROLASTHALF(BSatIonsph, "backscatr::BSatIonsph");

			// 1.3. Translate BS dNflux at Ionosphere (dist binned) to dNflux at Satellite (sat binned)
			dNflux_v2D ret{ backscat::ionsphToSatellite(eomdata, BSatIonsph) };
			// output: 2D vector of bs dNFlux at satellite per PA, E (sat binned) - should only be upward (90-180)
			// Section 1 End

			// 2.1. Convert from dNflux to dEflux
			for (unsigned int ang = 0; ang < ret.size(); ang++)
			{
				for (unsigned int eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= eomdata.satbins.E.at(eng);
			}
			// output: 2D vector of bs dEflux


			return ret;
		}
	}

	namespace binning
	{
		DLLEXP dNflux_v2D binParticles(const ParticleData& particles, const Bins& bins, const dNflux_v1D& countPerParticle)
		{
			dNflux_v2D ret(bins.PA.size(), dNflux_v1D(bins.E.size()));

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
			int epsilonNeeded{ 0 };

			for (unsigned int part = 0; part < particles.energy.size(); part++) //iterate over particles
			{
				double partEnerg{ particles.energy.at(part) };
				double partPitch{ particles.pitch.at(part) };

				if ((partEnerg == 0.0 && partPitch == 0.0) || countPerParticle.at(part) == 0.0)
					continue;                                               //guards - if particle E, PA is zero, it wasnt detected - just skip it
				if (partEnerg > pow(10, log10(Emax) + (0.5 * dlogE_bin)) || //if particle is outside E, PA measured limits - skip it
					partEnerg < pow(10, log10(Emin) - (0.5 * dlogE_bin)) || //PA shouldn't be an issue, there just in case
					partPitch > Amax + (0.5 * dangle_bin) ||                //also, code counts how many are outside limits for diagnostics
					partPitch < Amin - (0.5 * dangle_bin))
				{
					outsideLimits++;
					continue;
				}

				//calculate bin index for E, PA of particle
				int angbin{ (int)(partPitch / dangle_bin) }; //this should give the bin index
				int engbin{ (int)((log10(partEnerg) - (logMinEBinMid - 0.5 * dlogE_bin)) / dlogE_bin) }; //ditto
				if (!Eascending) engbin = (int)bins.E.size()  - 1 - engbin; //reverses the bin index if E is descending
				if (!Aascending) angbin = (int)bins.PA.size() - 1 - angbin; //ditto for PA

				double Ebin_min{ pow(10, log10(bins.E.at(engbin)) - 0.5 * dlogE_bin) };
				double Ebin_max{ pow(10, log10(bins.E.at(engbin)) + 0.5 * dlogE_bin) };
				double Abin_min{ bins.PA.at(angbin) - 0.5 * dangle_bin };
				double Abin_max{ bins.PA.at(angbin) + 0.5 * dangle_bin };
				double bin_epsilon{ 10000.0 * DBL_EPSILON };

				if (partEnerg >= Ebin_min &&
					partEnerg <  Ebin_max &&
					partPitch >= Abin_min &&
					partPitch <  Abin_max)
				{
					ret.at(angbin).at(engbin) += countPerParticle.at(part);
				}
				else if (partEnerg >= Ebin_min / (1.0 + bin_epsilon) &&
						 partEnerg <  Ebin_max * (1.0 + bin_epsilon) &&
					  	 partPitch >= Abin_min / (1.0 + bin_epsilon) &&
						 partPitch <  Abin_max * (1.0 + bin_epsilon))
				{ //expands bin boundaries by bin_epsilon (something like 2.2e-10 % -> 1 + 2.2e-12) to mitigate floating point error
					ret.at(angbin).at(engbin) += countPerParticle.at(part);
					epsilonNeeded++;
				}
				else //this shouldn't ever execute, guards should prevent zero and out of limits values
					throw logic_error("ionosphere::binning::binWeighted: Particle does not belong in bin identified for it. "
						+ string("\nInd: ") + to_string(part)
						+ "\ndangle, dlogE: " + to_string(dangle_bin) + ", " + to_string(dlogE_bin)
						+ "\nPart Energy: " + to_string(particles.energy.at(part))
						+ "\nEnergy Bin: " + to_string(engbin)
						+ "\nBin Min, Max: " + to_string(pow(10, log10(bins.E.at(engbin)) - 0.5 * dlogE_bin))
						     + ", " + to_string(pow(10, log10(bins.E.at(engbin)) + 0.5 * dlogE_bin))
						+ "\nConditions - GTE, LT: " + to_string(partEnerg >= pow(10, log10(bins.E.at(engbin)) - 0.5 * dlogE_bin)) + ", "
						     + to_string(partEnerg < pow(10, log10(bins.E.at(engbin)) + 0.5 * dlogE_bin))
						+ "\nPart Pitch: " + to_string(particles.pitch.at(part))
						+ "\nPitch Bin: " + to_string(angbin)
					    + "\nBin Min, Max: " + to_string(bins.PA.at(angbin) - 0.5 * dangle_bin)
						     + ", " + to_string(bins.PA.at(angbin) + 0.5 * dangle_bin)
						+ "\nConditions - GTE, LT: " + to_string(partPitch >= (bins.PA.at(angbin) - 0.5 * dangle_bin))
						     + ", " + to_string(partPitch < bins.PA.at(angbin) + 0.5 * dangle_bin));
			}

			if (outsideLimits > 0) cout << "ionosphere::binning::binWeighted : Particles out of limits: " << outsideLimits << "\n";
			if (epsilonNeeded > 0) cout << "ionosphere::binning::binWeighted : Particles needing bin boundaries expanded by ~2.2e-10 %: " << epsilonNeeded << "\n";

			return ret;
		}

		DLLEXP void symmetricBins0To360(dEflux_v2D& data, double_v1D& binAngles) //takes bins from 0-180 and extends to 360 symmetrically reflected across 180 degrees
		{//data[angle][energy]
			binAngles.resize(2 * binAngles.size()); //double size of angles vector (goes from 0-180 now)
			data.resize(binAngles.size()); //set data (vector of double vectors) to size of binAngles (now double sized)
			
			for (unsigned int ang = (unsigned int)binAngles.size() / 2; ang < (unsigned int)binAngles.size(); ang++)
			{
				binAngles.at(ang) = 2 * binAngles.at(ang - 1) - binAngles.at(ang - 2);
				data.at(ang) = data.at(binAngles.size() - ang - 1);
			}
		} //needs a good testing
	} //end namespace ionosphere::numerical

	namespace backscat //ionosphere::backscat
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

		constexpr double JOHND_PRIM_LOGM{ 3.5 }; //was 3.5 (from John) //0.8, 0.4 matches high E bump
		constexpr double JOHND_PRIM_LOGB{ -5.8 }; //was -5.8 (from John) //-5.5 matches high E bump
		constexpr double JOHND_SECD_LOGM_LT10{ -1.0 }; //was -1.0 (from John)
		constexpr double JOHND_SECD_LOGB_LT10{ -2.5 }; //was -2.5 (from John)
		constexpr double JOHND_SECD_LOGM_GT10{ -2.25 }; //was -2.25 (from John) //was at -4.0
		constexpr double JOHND_SECD_LOGB_GT10{ -1.7 }; //was -1.7 (from John) //was at -1.0

		DLLEXP dNflux johnd_flux(double E_eval, double E_incident, dEflux dE_incident)
		{
			if (E_eval > E_incident * 1.001) throw logic_error("johnd_flux: E_eval is higher than E_incident.  Not physical.  Eval, Incident: " + to_string(E_eval) + " , " + to_string(E_incident));

			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			return dE_incident * pow(10.0, secd_logm * log10(E_eval) + secd_logb) + //secondary BS
				   dE_incident * (10000.0 / E_incident) * pow(10.0, JOHND_PRIM_LOGM * log10(E_eval / E_incident) + JOHND_PRIM_LOGB); //primary BS
		}

		DLLEXP dEflux integralJohnd_flux(double lower, double upper, double E_incident)
		{
			throw std::exception("integralJohnd_flux used");//make sure this isnt used for now
			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			double integral_sec{ (pow(upper, secd_logm + 1.0) - pow(lower, secd_logm + 1.0)) * pow(10.0, secd_logb) / (secd_logm + 1.0) };
			double integral_prm{ (pow(upper, JOHND_PRIM_LOGM + 1.0) - pow(lower, JOHND_PRIM_LOGM + 1.0)) * pow(10.0, JOHND_PRIM_LOGB + 4.0) / ((JOHND_PRIM_LOGM + 1.0) * pow(E_incident, JOHND_PRIM_LOGM + 1.0)) };
			return integral_sec + integral_prm;
		}

		DLLEXP dNflux_v2D downwardToBackscatter(const Bins& dist, const dNflux_v2D& downwardAtIonsph)
		{ //converts downward dNflux at ionosphere (dist binned) to bs (upward) dNflux (also dist binned)
			// 1. Sum dNflux over PA Bins (dist bins), Per E Bin and Average
			dNflux_v1D dNsumPerE(dist.E.size());                   //Sum of escaped particles at each energy, units of dNflux
			for (unsigned int egy = 0; egy < dist.E.size(); egy++) //iterate over energies
			{
				for (auto& angvec : downwardAtIonsph)              //sum over angles, add to sum vector at energyBin
					dNsumPerE.at(egy) += angvec.at(egy);

				dNsumPerE.at(egy) /= ((double)dist.PA.size() / 2); //average over bins
			}
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins

			// 2. Calculate upward dNflux (backscatter) per E bin
			dNflux_v1D dNbsPerE(dist.E.size());
			for (unsigned int incEbin = 0; incEbin < dist.E.size(); incEbin++)
			{
				double E_incident{ dist.E.at(incEbin) };
				double dEflux_incBin{ dNsumPerE.at(incEbin) * E_incident };

				for (unsigned int evalEbin = 0; evalEbin <= incEbin; evalEbin++)
				{
					double E_eval{ dist.E.at(evalEbin) };
					
					dNbsPerE.at(evalEbin) += johnd_flux(E_eval, E_incident, dEflux_incBin);
				}
			}
			// output: 1D vector of the upgoing (backscatter) dNflux per E
			
			// 3. Distribute BS dNflux Isotropically Over Pitch Bins
			dNflux_v2D BS(dist.PA.size());
			for (unsigned int ang = 0; ang < BS.size(); ang++)
			{
				if (dist.PA.at(ang) <= 90.0)
					BS.at(ang) = dNflux_v1D(dist.E.size()); //empty vector of the right size
				else
				{
					BS.at(ang) = dNbsPerE;

					for (unsigned int eny = 0; eny < BS.at(ang).size(); eny++)
						BS.at(ang).at(eny) *= -cos(dist.PA.at(ang) * RADS_PER_DEG);

				}
			}
			// output: 2D vector of bs dEFlux at ionosphere per pitch, energy bin - should only be upward (90-180)

			return BS;
		}

		DLLEXP dNflux_v2D ionsphToSatellite(const EOMSimData& eomdata, const dNflux_v2D& bsCounts)
		{
			double Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) };

			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleData particles;

			for (unsigned int ang = 0; ang < eomdata.distbins.PA.size(); ang++)
			{
				for (unsigned int eny = 0; eny < eomdata.distbins.E.size(); eny++) //this works because ionospheric bins are same as distribution
				{ //ionospheric source particles detected moving upward
					particles.pitch.push_back(eomdata.upward.pitch.at(ang * eomdata.distbins.E.size() + eny));
					particles.energy.push_back(eomdata.upward.energy.at(ang * eomdata.distbins.E.size() + eny));
					particles.count.push_back(bsCounts.at(ang).at(eny) * Aratio_ion_sat);

					if (eomdata.dnward.energy.at(ang * eomdata.distbins.E.size() + eny) > 0.0)
					{ //ionospheric source particles detected moving downward (QSPS, Alfven)
						particles.pitch.push_back(eomdata.dnward.pitch.at(ang * eomdata.distbins.E.size() + eny));
						particles.energy.push_back(eomdata.dnward.energy.at(ang * eomdata.distbins.E.size() + eny));
						particles.count.push_back(bsCounts.at(ang).at(eny) * Aratio_ion_sat);
					}
				}
			}

			return binning::binParticles(particles, eomdata.satbins, particles.count);
		}
	} //end namespace ionosphere::backscat

	namespace multiLevelBS
	{
		DLLEXP dNflux_v2D scatterMain(const EOMSimData& eom, const dNflux_v2D& ionsphTopLvl)
		{
			double_v2D pctScattered(eom.distbins.PA.size(), dNflux_v1D(eom.distbins.E.size())); //% scattered per bin

			ParticleData particles;

			auto newPA = [](degrees PA_init, tesla B_init, tesla B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle has reflects before B_final

				degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations

				return ret;
			};

			// >> level calculate
			for (unsigned int level = 0; level < eom.ionsph.s.size() - 1; level++)
			{//for now, depends on adding one extra level (to see if particles reflect somewhere within the last layer)
				printLayer(eom.ionsph, level);

				dNflux_v2D bs_level{ bsAtLevel(eom, ionsphTopLvl, pctScattered, level) };

				TESTVEC_ISZEROLASTHALF(bs_level, "scatterMain::bs_level");

				double Aratio_bslevel_ion{ std::sqrt(eom.B_ion / eom.ionsph.B.at(level)) };

				// >> >> adjust PA
				for (unsigned int ang = 0; ang < eom.distbins.PA.size(); ang++)
				{
					double newpa{ newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(level), eom.ionsph.B.at(0)) };
					if (newpa < 0.0) throw logic_error(string("ionosphere::multiLevelBS::scatterMain: return value of newPA is -1 ")
						+ "indicating reflection before B_final.  This shouldn't happen as B_final is less in magnitude than B_initial. "
						+ " The particle can't physically reflect in this case, but simply be pushed closer to 180.  Something is wrong.");

					for (unsigned int eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (bs_level.at(ang).at(eny) != 0.0)
						{
							//#define count s_pos <- defined above in bsSrcToSat
							//allows me to use an existing vector in ParticleData with a name that makes sense
							//so, particles.count uses particles.s_pos, but renames it to how I want to use it
							particles.pitch.push_back(newpa);
							particles.energy.push_back(eom.distbins.E.at(eny));
							particles.count.push_back(bs_level.at(ang).at(eny) * Aratio_bslevel_ion); //from level to top ionsph
							//these vectors need to get cleaned out as the sim runs - otherwise all layers
							//have their entries in the vectors leading to excessive memory use
							//maybe multithread it?  Multiple sets of vectors?
						}
					}
				}
			}
			
			return binning::binParticles(particles, eom.distbins, particles.count);
		}

		DLLEXP dNflux_v2D bsAtLevel(const EOMSimData& eom, const dNflux_v2D& ionsphTopLvl, double_v2D& sumCollideAbove, unsigned int level)
		{
			ParticleData dngoing;
			ParticleData upgoing;

			double Aratio_ion_bslevel{ std::sqrt(eom.ionsph.B.at(level) / eom.B_ion) };

			//lambda that generates new PA from an initial PA based on B field strength change
			auto newPA = [](degrees PA_init, tesla B_init, tesla B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

				degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

				return ret;
			};

			//lambda that generates scatter/reflect percentage and updates total scatter % of defined particle
			auto sctReflPct = [&](double& sumSctAbove, eV E, degrees pitch)
			{
				if (pitch > 90.0) //upgoing - whatever hasn't scattered so far, reflects
					return (1.0 - sumSctAbove); //so return 1.0 - "sum % scattered in layers above"

				//dngoing - percent that scatters in this layer
				double sct{ 0.0 };

				for (size_t species = 0; species < eom.ionsph.p.size(); species++)
				{
					sct += scatterPct(sumSctAbove, eom.ionsph.Z.at(species), eom.ionsph.p.at(species).at(level),
						eom.ionsph.h.at(level), E, pitch);
				}

				if (sct + sumSctAbove > 1.0)
				{ //no more than 100% of particles can have scattered
					sct = 1.0 - sumSctAbove; //sct% for this layer is equal to however much is left until the total scattered is 100%
					sumSctAbove = 1.0;
				}
				else
				{ //sum is less than 100%, so add to the total as normal
					sumSctAbove += sct;
				}

				if (sct < 0.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is < 0.0 - sct%, %collideAbove: "
					+ to_string(sct) + ", " + to_string(sumSctAbove));
				if (sct > 1.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is > 1.0 - sct%, %collideAbove: "
					+ to_string(sct) + ", " + to_string(sumSctAbove));

				return sct;
			};

			// >> >> adjust PA, check if reflected
			for (unsigned int ang = 0; ang < eom.distbins.PA.size(); ang++)
			{
				if (eom.distbins.PA.at(ang) > 90.0) continue; //there should be no upgoing in ionsphTopLvl

				degrees pa_level{ newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level)) };
				degrees pa_nextlevel{ newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level + 1)) };

				if (pa_level > 90.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: downgoing particle ended up with pitch > 90.0 - PA_bin, PA_level, level: "
					+ to_string(eom.distbins.PA.at(ang)) + ", " + to_string(pa_level) + ", " + to_string(level));
				
				if (pa_level < 0.0)
				{ //particle reflects before this level
					continue;
				}
				else if (pa_nextlevel < 0.0)
				{ //particle reflects before next level, add all particles of this pitch moving in the opposite direction
					for (unsigned int eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (ionsphTopLvl.at(ang).at(eny) != 0.0 && sumCollideAbove.at(ang).at(eny) < 1.0)
						{ //these particles are upgoing
							degrees pa_sat{ newPA(180.0 - pa_level, eom.ionsph.B.at(level), eom.B_sat) }; //pitch at satellite - eventually may need to run through sourceToSatellite
							dNflux layerEPAcount{ ionsphTopLvl.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* sctReflPct(sumCollideAbove.at(ang).at(eny), eom.distbins.E.at(eny), 180.0 - pa_level) //percent reflected
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ -cos(pa_sat * RADS_PER_DEG) }; //cos of pitch at satellite (where detected)
							
							//#define count s_pos <- defined above in sourceToSatellite
							//allows me to use an existing vector in ParticleData with a name that makes sense
							//so, particles.count uses particles.s_pos, but renames it to how I want to use it
							upgoing.energy.push_back(eom.distbins.E.at(eny));
							upgoing.pitch.push_back(180.0 - pa_level);
							upgoing.count.push_back(layerEPAcount);
						}
					}
				}
				else
				{ //particle makes it to the next level - invoke scattering
					for (unsigned int eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (ionsphTopLvl.at(ang).at(eny) != 0.0 && sumCollideAbove.at(ang).at(eny) < 1.0)
						{ //these particles are downgoing
							dNflux layerEPAcount{ ionsphTopLvl.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* sctReflPct(sumCollideAbove.at(ang).at(eny), eom.distbins.E.at(eny), pa_level) //percent scattered
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ cos(pa_level * RADS_PER_DEG) }; //pitch at lowest level (level of scattering)

							dngoing.energy.push_back(eom.distbins.E.at(eny));
							dngoing.pitch.push_back(pa_level);
							dngoing.count.push_back(layerEPAcount);
						}
					}
				}
			}	

			TESTVEC_NOTNEGWHOLEVEC(upgoing.count, "bsAtLevel:: " + to_string(level) + " :: upgoing.count");
			TESTVEC_NOTNEGWHOLEVEC(dngoing.count, "bsAtLevel:: " + to_string(level) + " :: dngoing.count");

			// >> >> bin particles at level
			dNflux_v2D dnBinned{ binning::binParticles(dngoing, eom.distbins, dngoing.count) };
			dNflux_v2D upBinned{ binning::binParticles(upgoing, eom.distbins, upgoing.count) };
			
			TESTVEC_ISZEROFRSTHALF(dnBinned, "bsAtLevel:: " + to_string(level) + " ::dnBinned");
			TESTVEC_ISZEROLASTHALF(upBinned, "bsAtLevel:: " + to_string(level) + " ::upBinned");
			
			// >> >> calculate backscatter
			dNflux_v2D ret{ backscat::downwardToBackscatter(eom.distbins, dnBinned) };

			TESTVEC_ISZEROLASTHALF(ret, "bsAtLevel:: " + to_string(level) + " ::backscatter_level");
			
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eny = 0; eny < ret.at(0).size(); eny++)
					ret.at(ang).at(eny) += upBinned.at(ang).at(eny);

			return ret;
		}

		DLLEXP double_v2D alt_reflect(const Bins& distbins, BField* B, double B_ion, double t)
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
						if (count > 1000000) throw logic_error("ionosphere::multiLevelBS::alt_reflect: loop count is over a million - that shouldn't happen: delta_s, s_guess, B_guess: " +
							to_string(delta_s) + ", " + to_string(s_guess) + ", " + to_string(B_guess));
					}
					
					//third: shrink delta_s, reset loop count
					delta_s /= 5.0;
					count = 0;

					//guard in case delta_s is ever less than err - small changes in s will result in very small B changes
					//if delta_s is as low as err, B is really not changing much - trying to prevent infinite loops
					if (delta_s < err) throw logic_error("ionosphere::multiLevelBS::alt_reflect: delta_s is less than error in bracketing lambda: B_target, delta_s, err"
						+ to_string(B_target) + ", " + to_string(delta_s) + ", " + to_string(err));
				}

				//return the result
				return s_guess;
			};

			double_v2D s_ref(distbins.PA.size(), double_v1D(distbins.E.size()));

			for (size_t ang = 0; ang < distbins.PA.size(); ang++)
			{
				double s_ang{ sAtB_bracket(B_ion / sin(distbins.PA.at(ang) * RADS_PER_DEG)) };

				for (size_t eny = 0; eny < distbins.E.size(); eny++)
				{
					s_ref.at(ang).at(eny) = s_ang;
				}
			}

			return s_ref;
		}

		DLLEXP double scatterPct(double sumCollideAbove, double Z, double p, double h, eV E, degrees PA)
		{
			if (sumCollideAbove >= 1.0) throw logic_error("ionosphere::multiLevelBS::scatterPct: sumCollideAbove is greater than / equal to 1.0.  "
				+ string("100% has scattered already.  Conditions in bsAtLevel should have prevented this from happening."));

			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	} //end namespace multiLevelBS

} //end namespace ionosphere
