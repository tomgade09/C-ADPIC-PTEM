#include "ionosphere/ionosphere.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

#include "utils/fileIO.h"
#include "utils/writeIOclasses.h"
#include "ErrorHandling/simExceptionMacros.h"

using std::pow;
using std::cout;
using std::string;
using std::to_string;
using std::exception;
using std::logic_error;
using std::out_of_range;
using utils::fileIO::readDblBin;
using utils::numerical::generateSpacedValues;
using ionosphere::Bins;
using ionosphere::ParticlesBinned;
using ionosphere::IonosphereSpecs;

//diagnostic things
//remove later
vector<ParticlesBinned<dNflux>> UPWARDDATA_G;
vector<eV> EBINS_DIAG{ generateSpacedValues(0.5, 4.5, 5, true, true) };
//vector<degrees> PABINS_DIAG{ generateSpacedValues(2.5, 177.5, 36, false, true) };
//vector<degrees> PABINS_DIAG{ generateSpacedValues(1.25, 178.75, 72, false, true) };
//vector<degrees> PABINS_DIAG{ generateSpacedValues(0.625, 179.375, 144, false, true) };
//vector<degrees> PABINS_DIAG{ generateSpacedValues(90.0 / 216.0, 180.0 - 90.0 / 216.0, 216, false, true) };
vector<degrees> PABINS_DIAG{ generateSpacedValues(0.3125, 179.6875, 288, false, true) };
//
//

inline void printIonosphere(const IonosphereSpecs& ionsph)
{
	auto stringConcat = [&]()
	{
		string ret{ "" };

		for (size_t spec = 0; spec < ionsph.names.size(); spec++)
		{
			ret += ionsph.names.at(spec);
			ret += (spec == ionsph.names.size() - 1) ? ("") : (", ");
		}

		return ret;
	};

	cout << "==================== Ionosphere Simulation ====================" << "\n";
	cout << "Min, Max s (m): " << *(ionsph.s.end() - 2) << ", " << ionsph.s.front() << "\n";
	cout << "Min, Max B (T): " << *(ionsph.B.end() - 2) << ", " << ionsph.B.front() << "\n";
	cout << "Num of Layers : " << ionsph.s.size() - 1 << "\n";
	cout << "Atomic Species: " << stringConcat() << "\n";
	cout << "===============================================================" << "\n";
}

inline void printLayer(const IonosphereSpecs& ionsph, size_t layer)
{
	cout << "Layer: " << layer << " / " << ionsph.s.size() - 2 << ", s: ";
	cout << ionsph.s.at(layer) << ", B: " << ionsph.B.at(layer) << "\n";
}

namespace ionosphere
{
	namespace debug
	{
		void eomError(const EOMSimData& ideal, const EOMSimData& eomsim)
		{
			cout << "Need to finish.  Returning.\n";
			return;

			auto err = [](double base, double dev) { return abs((base - dev) / dev); };

			vector<double> maxErr(8); //topPA, btmPA, upwPA, dnwPA, topE, btmE, upwE, dnwE
			vector<int> extra(4); //top, btm, upw, dnw
			
			//load top data
			ParticleList eomtop;
			readDblBin(eomtop.vpara, eomsim.datadir + "bins\\satellites\\topElec_vpara.bin");
			readDblBin(eomtop.vperp, eomsim.datadir + "bins\\satellites\\topElec_vperp.bin");
			readDblBin(eomtop.s_pos, eomsim.datadir + "bins\\satellites\\topElec_s.bin");
			utils::numerical::v2DtoEPitch(eomtop.vpara, eomtop.vperp, eomsim.mass, eomtop.energy, eomtop.pitch);

			for (int part = 0; part < eomsim.initial.energy.size(); part++)
			{
				//does ideal dist have more detected particles?
				//top condition here - will require newPA
				if ((ideal.ionosph.pitch.at(part) != 0.0) && (eomsim.ionosph.pitch.at(part) == 0.0))
					extra.at(1)++;
				if ((ideal.upgoing.pitch.at(part) != 0.0) && (eomsim.upgoing.pitch.at(part) == 0.0))
					extra.at(2)++;
				if ((ideal.dngoing.pitch.at(part) != 0.0) && (eomsim.dngoing.pitch.at(part) == 0.0))
					extra.at(3)++;

				//compare pitch, E of top, btm, upw, dnw
				if ((ideal.ionosph.pitch.at(part) != 0.0) && (eomsim.ionosph.pitch.at(part) != 0.0) && (err(ideal.ionosph.pitch.at(part), eomsim.ionosph.pitch.at(part)) > maxErr.at(1)))
					maxErr.at(1) = err(ideal.ionosph.pitch.at(part), eomsim.ionosph.pitch.at(part));
				if ((ideal.upgoing.pitch.at(part) != 0.0) && (eomsim.upgoing.pitch.at(part) != 0.0) && (err(ideal.upgoing.pitch.at(part), eomsim.upgoing.pitch.at(part)) > maxErr.at(2)))
					maxErr.at(2) = err(ideal.upgoing.pitch.at(part), eomsim.upgoing.pitch.at(part));
				if ((ideal.dngoing.pitch.at(part) != 0.0) && (eomsim.dngoing.pitch.at(part) != 0.0) && (err(ideal.dngoing.pitch.at(part), eomsim.dngoing.pitch.at(part)) > maxErr.at(3)))
					maxErr.at(3) = err(ideal.dngoing.pitch.at(part), eomsim.dngoing.pitch.at(part));
			}
		}

		void outputVectorsToCSV(string filename, const vector<vector<double>>& data, const vector<string>& labels, int writeEvery = 1, size_t startIdx = 0, size_t stopIdx = 0)
		{
			utils::fileIO::CSV file{ filename };

			if (writeEvery == 1 && startIdx == 0 && stopIdx == 0)
			{
				file.add(data, labels);
			}
			else
			{
				vector<vector<double>> dataSmall;
				vector<string> labelsSmall;

				bool labelsError{ false };

				if (stopIdx == 0) stopIdx = data.size();
				for (size_t idx = startIdx; idx < stopIdx; idx += writeEvery)
				{
					dataSmall.push_back(data.at(idx));
					try
					{
						labelsSmall.push_back(labels.at(idx));
					}
					catch (out_of_range)
					{
						if (!labelsError)
							cout << "Warning: File name " + filename + " does not have as many labels as data indicies.  Things may be labeled incorrectly.\n";
						
						labelsSmall.push_back("");
						labelsError = true;
					}
				}

				file.add(dataSmall, labelsSmall);
			}
		}

		void outputVectorsToCSV(string filename, const vector<vector<double>>& data, const vector<double>& labels, int writeEvery = 1, size_t startIdx = 0, size_t stopIdx = 0)
		{
			vector<string> strlabels;

			for (size_t idx = 0; idx < labels.size(); idx++)
				strlabels.push_back(to_string(labels.at(idx)));

			outputVectorsToCSV(filename, data, strlabels, writeEvery, startIdx, stopIdx);
		}

		void outputVectorsToCSV(string filename, const function<double(double, double)> datafunc, const vector<double>& outeridx, const vector<double>& inneridx, const vector<string> labels, int writeEvery = 1, size_t startIdx = 0, size_t stopIdx = 0)
		{
			vector<vector<double>> data(outeridx.size(), vector<double>(inneridx.size()));

			for (size_t outer = 0; outer < outeridx.size(); outer++)
			{
				for (size_t inner = 0; inner < inneridx.size(); inner++)
				{
					data.at(outer).at(inner) = datafunc(outeridx.at(outer), inneridx.at(inner));
				}
			}

			outputVectorsToCSV(filename, data, labels, writeEvery, startIdx, stopIdx);
		}

		void outputVectorsToCSV(string filename, const function<double(double, double)> datafunc, const vector<double>& outeridx, const vector<double>& inneridx, const vector<double> labels, int writeEvery = 1, size_t startIdx = 0, size_t stopIdx = 0)
		{
			vector<string> strlabels;

			for (size_t idx = 0; idx < labels.size(); idx++)
				strlabels.push_back(to_string(labels.at(idx)));

			outputVectorsToCSV(filename, datafunc, outeridx, inneridx, strlabels, writeEvery, startIdx, stopIdx);
		}
	}

	DLLEXP vector<vector<dEflux>> steadyFlux(const EOMSimData& eomdata)
	{
		//
		//
		// Spit out csv files (for diagnostics)
		{
			auto gendNflux = [](double E_eval, double E_incident)
			{
				return backscat::johnd_flux(E_eval, E_incident, 1.0);
			};

			debug::outputVectorsToCSV("debug/backscatter/dNflux_1electron.csv", gendNflux, eomdata.distbins.E, eomdata.distbins.E, eomdata.distbins.E);
			debug::outputVectorsToCSV("debug/ionospheric densities/densities.csv", eomdata.ionsph.p, eomdata.ionsph.names);
		}
		// End spit out csv files
		//
		//

		printIonosphere(eomdata.ionsph);
		
		TESTVEC_NOTNEGWHOLEVEC(eomdata.maxwellian, "Maxwellian distribution");
		
		cout << "\n\nUsing nasty global.  Make sure to remove after diagnostics.\n\n\n";

		// 2. Calculate dEfluxes
		ParticlesBinned<dEflux> ptemflux{ dEFlux::satellite(eomdata) }; // dE flux due to PTEM simulation
		ParticlesBinned<dEflux> isphflux{ dEFlux::backscatr(eomdata) }; // dE flux due to reflection and backscatter in ionosphere

		printVec2D(ptemflux.binnedData, "dE flux at Satellite");
		printVec2D(isphflux.binnedData, "Flux Due to Backscatter at Satellite");

		//debug::outputVectorsToCSV("debug/final dE/downward.csv", distfluxdnward, eomdata.satbins.PA);
		//debug::outputVectorsToCSV("debug/final dE/upward.csv", distfluxupward, eomdata.satbins.PA);
		//debug::outputVectorsToCSV("debug/final dE/backscatter.csv", bkscfluxupward, eomdata.satbins.PA);

		#define totalflux ptemflux //allows me to rename the data to what I'm using it for
		// 3. Sum dEfluxes
		totalflux += isphflux;

		TESTVEC_NOTNEGWHOLEVEC(totalflux.binnedData, "Total Flux at End");

		return totalflux.binnedData;
		#undef totalflux //undefine so the rest of the code isn't affected
	}


	namespace dEFlux
	{
		DLLEXP ParticlesBinned<dEflux> satellite(const EOMSimData& eomdata)
		{
			// Section 1 - Adjust dNflux for altitude and cos factor
			ratio Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph to satellite
			ratio Aratio_mag_sat{ std::sqrt(eomdata.B_sat / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to satellite
			if (isnan(Aratio_ion_sat) || isnan(Aratio_mag_sat)) throw logic_error("dEFlux::satellite: Aratio is NaN - negative value inside square root");

			dNflux_v1D distribution_sat{ eomdata.maxwellian }; //becomes dE flux when multiplied by energy in below loop

			for (size_t part = 0; part < eomdata.initial.s_pos.size(); part++) //isotropize counts -> 3D
			{
				if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)    //ionospheric source, upgoing
				{
					if (eomdata.initial.pitch.at(part) < 90.0) throw logic_error("dEFlux::satellite: ionospheric source particle has downgoing pitch");
					if (eomdata.upgoing.pitch.at(part) < 90.0) throw logic_error("dEFlux::satellite: ionospheric source particle is downgoing at upgoing detector");
					distribution_sat.at(part) *= -cos(eomdata.initial.pitch.at(part) * RADS_PER_DEG) * Aratio_ion_sat
						* eomdata.upgoing.energy.at(part);
				}
				else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)//magnetospheric source, downgoing
				{
					if (eomdata.initial.pitch.at(part) > 90.0) throw logic_error("dEFlux::satellite: magnetospheric source particle has downgoing pitch");
					if (eomdata.dngoing.pitch.at(part) > 90.0) throw logic_error("dEFlux::satellite: magnetospheric source particle is upgoing at downgoing detector");
					distribution_sat.at(part) *= 1.0 / cos(eomdata.dngoing.pitch.at(part) * RADS_PER_DEG) * Aratio_mag_sat
						* eomdata.dngoing.energy.at(part);
				}
				else
					throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source - it's somewhere between the ionosphere and magnetosphere");

				if (eomdata.ionosph.pitch.at(part) >= 90.0)
					throw logic_error("ionosphere::steadyFlux : pitch of particle at ionosphere is >= 90 deg");
			}

			TESTVEC_NOTNEGWHOLEVEC(distribution_sat, "satellite::maxwellian_sat");

			// Section 2 - Get dEflux at Satellite - upward and downward
			// 1. Bin Particles by Satellite Detector PA, Energy Bins
			ParticlesBinned<dEflux> upgoing{ eomdata.satbins.binParticleList(eomdata.upgoing, distribution_sat) };
			ParticlesBinned<dEflux> dngoing{ eomdata.satbins.binParticleList(eomdata.dngoing, distribution_sat) };
			
			TESTVEC_ISZEROFRSTHALF(upgoing.binnedData, "satellite::upgoing");
			TESTVEC_ISZEROLASTHALF(dngoing.binnedData, "satellite::dngoing");

			// 2. Sum upward and downward; convert from dNflux to dEflux
			#define total upgoing //using upward as the total / return array - just renaming to make that more obvious
			total += dngoing;

			return total;
			#undef total //don't want this affecting code below
		}

		DLLEXP ParticlesBinned<dEflux> backscatr(const EOMSimData& eomdata)
		{
			// Section 1 - Adjust dNflux for altitude
			ratio Aratio_ion_ion{ std::sqrt(eomdata.B_ion / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph (up, reflect, down) to ionsph
			ratio Aratio_mag_ion{ std::sqrt(eomdata.B_ion / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to ionsph
			if (isnan(Aratio_ion_ion) || isnan(Aratio_mag_ion)) throw logic_error("dEFlux::backscatr: Aratio is NaN - negative value inside square root");

			dNflux_v1D distribution_ion{ eomdata.maxwellian };

			for (size_t part = 0; part < eomdata.initial.s_pos.size(); part++) //isotropize counts -> 3D
			{
				if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)     //ionospheric source, upgoing
					distribution_ion.at(part) *= Aratio_ion_ion; //should just be equal to 1, but gives it symmetry
				else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)//magnetospheric source, downgoing
					distribution_ion.at(part) *= Aratio_mag_ion; //keep in mind: no cos factor yet - will be determined in reflect/scatter code
				else
					throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source - it's somewhere between the ionosphere and magnetosphere");
			}

			TESTVEC_NOTNEGWHOLEVEC(distribution_ion, "satellite::distribution_ion");

			// Section 2 - Get escaping dNflux at ionosphere top
			// 2.1. Bin Escaped Particles by EOM simulation (distribution) PA, Energy Bins (for high res backscatter)
			ParticlesBinned<dNflux> escapeIntoIonsph{ eomdata.distbins.binParticleList(eomdata.ionosph, distribution_ion) };
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle

			TESTVEC_ISZEROFRSTHALF(escapeIntoIonsph.binnedData, "backscatr::escapeIntoIonsph"); //should be no upgoing escaping particles

			// 2.2. Calculate upgoing dNflux due to reflection, scattering from incident dNflux to ionosphere
			ParticlesBinned<dNflux> upgoingIonosphere{ multiLevelBS::scatterMain(eomdata, escapeIntoIonsph) }; //new multi-level hotness
			// output: 2D vector of reflection + backscatter dNFlux by dist bins at ionosphere (upgoing)

			TESTVEC_ISZEROLASTHALF(upgoingIonosphere.binnedData, "backscatr::upgoingIonosphere"); //should be no downgoing backscatter

			// 2.3. Translate BS dNflux at ionosphere (dist binned) to dEflux at satellite (sat binned)
			return backscat::ionsphToSatellite(eomdata, upgoingIonosphere);
			// output: 2D vector of bs dNFlux at satellite per PA, E (sat binned) - should only be upward (90-180)
		}

		DLLEXP double newPA(const degrees PA_init, const tesla B_init, const tesla B_final)
		{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
			if ((PA_init < 90.0 && B_init < B_final) ||
				(PA_init > 90.0 && B_init > B_final))
			{
				cout << "PA_init: " << PA_init << "\n";
				cout << "B_init:  " << B_init << "\n";
				cout << "B_final: " << B_final << "\n";
				throw logic_error(string("ionosphere::dEFlux::newPA: down[up]going particle with ")
					+ "decreasing[increasing] magnitude B Field.  This can't happen in this sim without "
					+ "time-independent fields which violates the conditions required for the function to be valid.");
			}

			double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

			if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

			degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

			if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

			return ret;
		}
	}


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

		constexpr double JOHND_PRIM_LOGM{ 0.505 };
		constexpr double JOHND_PRIM_LOGB{ -5.16 };
		constexpr double JOHND_SECD_LOGM_LT25{ -0.975 };
		constexpr double JOHND_SECD_LOGB_LT25{ -1.47 };
		constexpr double JOHND_SECD_LOGM_GT25{ -1.95 };
		constexpr double JOHND_SECD_LOGB_GT25{ -0.11 };

		DLLEXP dNflux johnd_flux(eV E_eval, eV E_incident, dNflux dN_incident)
		{
			if (E_eval > E_incident * (1.0 + FLT_EPSILON)) return 0.0;

			double secd_logm{ (E_incident <= 25.0) ? JOHND_SECD_LOGM_LT25 : JOHND_SECD_LOGM_GT25 };
			double secd_logb{ (E_incident <= 25.0) ? JOHND_SECD_LOGB_LT25 : JOHND_SECD_LOGB_GT25 };

			//
			//
			if (E_incident <= 25.0)
				return E_incident * dN_incident * pow(10.0, secd_logm * log10(25.0) + secd_logb) + //secondary BS
					   E_incident * dN_incident * (10000.0 / E_incident) * pow(10.0, JOHND_PRIM_LOGM * log10(E_eval / E_incident) + JOHND_PRIM_LOGB); //primary BS;
			//
			//

			return E_incident * dN_incident * pow(10.0, secd_logm * log10(E_eval) + secd_logb) + //secondary BS
				   E_incident * dN_incident * (10000.0 / E_incident) * pow(10.0, JOHND_PRIM_LOGM * log10(E_eval / E_incident) + JOHND_PRIM_LOGB); //primary BS
		}

		DLLEXP ParticlesBinned<dNflux> downwardToBackscatter(const Bins& dist, const ParticlesBinned<dNflux>& dNpointofScatter)
		{ //converts downward dNflux at point of scatter (dist binned) to bs (upward) dNflux (also dist binned)
			// 1. Sum dNflux over PA Bins (dist bins), Per E Bin and Average
			dNflux_v1D dNsumPerE(dist.E.size());             //Sum of escaped particles at each energy, units of dNflux

			for (size_t egy = 0; egy < dist.E.size(); egy++) //iterate over energies
			{
				for (size_t ang = 0; ang < dist.PA.size(); ang++)
				{
					if (dist.PA.at(ang) > 90.0) continue; //there's no downward flux in the upward direction

					degrees dangle{ (ang == 0) ?
						(abs(dist.PA.at(ang + 1) - dist.PA.at(ang))) : //if we are at first index, ang - 1 doesn't exist
						(abs(dist.PA.at(ang) - dist.PA.at(ang - 1))) };//if we are at last index, ang + 1 doesn't exist

					dNsumPerE.at(egy) += dNpointofScatter.binnedData.at(ang).at(egy) / (dist.PA.size() / 2) * 7500 * 4 * PI *
						abs(sin(dist.PA.at(ang) * RADS_PER_DEG) * sin(dangle / 2.0 * RADS_PER_DEG) *
							cos(dist.PA.at(ang) * RADS_PER_DEG) * cos(dangle / 2.0 * RADS_PER_DEG));
				}
			}
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins

			// 2. Calculate upward dNflux (backscatter) per E bin
			dNflux_v1D dNbsPerE(dist.E.size());
			for (size_t incEbin = 0; incEbin < dist.E.size(); incEbin++)
			{
				double E_incident{ dist.E.at(incEbin) };
				double dNflux_incBin{ dNsumPerE.at(incEbin) };

				for (size_t evalEbin = 0; evalEbin <= incEbin; evalEbin++)
				{
					double E_eval{ dist.E.at(evalEbin) };

					dNbsPerE.at(evalEbin) += johnd_flux(E_eval, E_incident, dNflux_incBin);
				}
			}
			// output: 1D vector of the upgoing (backscatter) dNflux per E

			// 3. Distribute BS dNflux Isotropically Over Pitch Bins
			ParticlesBinned<dNflux> BS;
			BS.binnedData = vector<vector<dNflux>>(dist.PA.size(), vector<dNflux>(dist.E.size()));

			for (size_t ang = 0; ang < dist.PA.size(); ang++)
			{
				if (dist.PA.at(ang) <= 90.0)
					continue; //empty vector of the right size
				else
				{
					for (size_t eny = 0; eny < dist.E.size(); eny++)
					{
						BS.binnedData.at(ang).at(eny) = dNbsPerE.at(eny) * -cos(dist.PA.at(ang) * RADS_PER_DEG) / (4 * PI);
					}
				}
			}
			// output: 2D vector of bs dEFlux at ionosphere per pitch, energy bin - should only be upward (90-180)

			return BS;
		}

		DLLEXP ParticlesBinned<dEflux> ionsphToSatellite(const EOMSimData& eomdata, const ParticlesBinned<dNflux>& bsdN)
		{
			ratio Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) };

			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleList particles;

			for (size_t ang = 0; ang < eomdata.distbins.PA.size(); ang++)
			{
				if (eomdata.distbins.PA.at(ang) <= 90.0) continue;

				for (size_t eny = 0; eny < eomdata.distbins.E.size(); eny++) //this works because ionospheric bins are same as distribution
				{ //ionospheric source particles detected moving upward
					size_t idx1D{ eomdata.initial1DToDistbins2DMapping.binnedData.at(ang).at(eny) };

					particles.pitch.push_back(eomdata.upgoing.pitch.at(idx1D));
					particles.energy.push_back(eomdata.upgoing.energy.at(idx1D));
					particles.count.push_back(bsdN.binnedData.at(ang).at(eny) * Aratio_ion_sat * eomdata.upgoing.energy.at(idx1D));

					if (eomdata.dngoing.energy.at(idx1D) > 0.0)
					{ //ionospheric source particles detected moving downward (QSPS)
						particles.pitch.push_back(eomdata.dngoing.pitch.at(idx1D));
						particles.energy.push_back(eomdata.dngoing.energy.at(idx1D));
						particles.count.push_back(bsdN.binnedData.at(ang).at(eny) * Aratio_ion_sat * eomdata.dngoing.energy.at(idx1D));
					}
					else if (eomdata.dngoing.energy.at(idx1D) < 0.0)
						throw logic_error("ionosphere::backscat::ionosphToSatellite: E is < 0.0, idx: " + to_string(idx1D));
				}
			}

			return eomdata.satbins.binParticleList(particles, particles.count);
			#undef count
		}
	} //end namespace ionosphere::backscat

	namespace multiLevelBS
	{
		DLLEXP ParticlesBinned<dNflux> scatterMain(const EOMSimData& eom, const ParticlesBinned<dNflux>& dNionsphTop)
		{
			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			
			//
			//
			// Spit out csv files (for diagnostics)
			debug::outputVectorsToCSV("debug/angle, pitch bins/satbins.csv", { eom.satbins.E, eom.satbins.PA }, vector<string>{ "E", "PA" });
			debug::outputVectorsToCSV("debug/angle, pitch bins/distbins.csv", { eom.distbins.E, eom.distbins.PA }, vector<string>{ "E", "PA" });
			debug::outputVectorsToCSV("debug/dNdwnwdIonosphereTop.csv", dNionsphTop.binnedData, eom.distbins.PA, 600, 18000, 36000);
			// End spit out csv files
			//
			//

			double_v2D pctScattered(eom.distbins.PA.size(), double_v1D(eom.distbins.E.size())); //% scattered so far per bin

			ParticleList particles;

			// >> level calculate
			for (size_t level = 0; level < eom.ionsph.s.size() - 1; level++)
			{//for now, depends on adding one extra level (to see if particles reflect somewhere within the last layer)
				printLayer(eom.ionsph, level);

				ParticlesBinned<dNflux> upg_level{ upgAtLevel(eom, dNionsphTop, pctScattered, level) };
				UPWARDDATA_G.push_back(upg_level);
				TESTVEC_ISZEROLASTHALF(upg_level.binnedData, "scatterMain::bs_level");

				ratio Aratio_bslevel_ion{ std::sqrt(eom.B_ion / eom.ionsph.B.at(level)) };

				// >> >> adjust PA
				for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
				{
					if (eom.distbins.PA.at(ang) <= 90.0) continue; //all particles are upgoing

					degrees PA_ionsphTop{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(level), eom.B_ion) };
					if (PA_ionsphTop < 0.0) throw logic_error(string("ionosphere::multiLevelBS::scatterMain: return value of newPA is -1 ")
						+ "indicating reflection before B_final.  This shouldn't happen as B_final is less in magnitude than B_initial. "
						+ " The particle can't physically reflect in this case, but simply be pushed closer to 180.  Something is wrong.");

					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						//doesn't gain anything - backscatter is deposited in every upgoing bin - consider removing the if condition
						//if (upg_level.at(ang).at(eny) != 0.0)
						//{
							particles.pitch.push_back(PA_ionsphTop);
							particles.energy.push_back(eom.distbins.E.at(eny));
							particles.count.push_back(upg_level.binnedData.at(ang).at(eny) * Aratio_bslevel_ion); //from level to top ionsph
							//these vectors need to get cleaned out as the sim runs - otherwise all layers
							//have their entries in the vectors leading to excessive memory use
							//maybe multithread it?  Multiple sets of vectors?
						//}
					}
				}
			}

			//
			// diagnostics
			// mirror / sum all below layers, output to csv
			{
				Bins diagbins(EBINS_DIAG, PABINS_DIAG);

				ParticlesBinned<dNflux> lastLayer;
				lastLayer.binnedData = vector<vector<dNflux>>(eom.distbins.PA.size(), vector<dNflux>(eom.distbins.E.size()));

				size_t toLevel = eom.ionsph.s.size() - 1;
				do
				{
					toLevel--;
					cout << "level: " << toLevel << "\n";
					
					ParticleList belowParticles;
					const ratio Aratio_from_to{ std::sqrt(eom.ionsph.B.at(toLevel) / eom.ionsph.B.at(toLevel + 1)) };

					for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
					{
						for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
						{
							if (eom.distbins.PA.at(ang) < 90.0)
							{
								if (UPWARDDATA_G.at(toLevel).binnedData.at(ang).at(eny) != 0.0 ||
									lastLayer.binnedData.at(ang).at(eny) != 0.0)
									throw logic_error("multiLevelBS::diags: downward distribution is non-zero");
								
								continue;
							}

							if (lastLayer.binnedData.at(ang).at(eny) == 0.0) continue;

							const degrees toPA{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(toLevel + 1), eom.ionsph.B.at(toLevel)) };

							belowParticles.pitch.push_back(toPA);
							belowParticles.energy.push_back(eom.distbins.E.at(eny));
							belowParticles.count.push_back(lastLayer.binnedData.at(ang).at(eny) * Aratio_from_to);
						}
					}

					lastLayer = eom.distbins.binParticleList(belowParticles, belowParticles.count);
					//lastLayer += UPWARDDATA_G.at(toLevel);
					for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
						for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
							lastLayer.binnedData.at(ang).at(eny) += UPWARDDATA_G.at(toLevel).binnedData.at(ang).at(eny); // add current layer dN flux

					ParticleList output;
					for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
					{
						if (eom.distbins.PA.at(ang) < 90.0) continue;
						for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
						{
							output.pitch.push_back(eom.distbins.PA.at(ang));
							output.energy.push_back(eom.distbins.E.at(eny));
							output.count.push_back(lastLayer.binnedData.at(ang).at(eny) * eom.distbins.E.at(eny));
						}
					}

					ParticlesBinned<dEflux> lastLayerDiagBinned{ diagbins.binParticleList(output, output.count) };

					debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/upSum_" + to_string((int)eom.ionsph.s.at(toLevel)) + ".csv",
						lastLayerDiagBinned.binnedData, diagbins.PA);
				} while (toLevel != 0);
			}
			//
			//
			//
			
			return eom.distbins.binParticleList(particles, particles.count);
			#undef count
		}

		DLLEXP ParticlesBinned<dNflux> upgAtLevel(const EOMSimData& eom, const ParticlesBinned<dNflux>& dNionsphTop, double_v2D& pctScatteredAbove, size_t level)
		{
			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			//
			//
			// Spit out csv files (for diagnostics)
			{
				auto sumSctPct = [&](double PA, double E)
				{
					double prob{ 0.0 };

					for (size_t species = 0; species < eom.ionsph.p.size(); species++)
					{
						prob += scatterPct(0.0, eom.ionsph.Z.at(species), eom.ionsph.p.at(species).at(level), eom.ionsph.h.at(species),
							E, PA);
					}

					return prob;
				};
				
				debug::outputVectorsToCSV("debug/level scatter probs/scatter_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
					sumSctPct, eom.satbins.PA, eom.satbins.E, eom.satbins.PA);
			}
			// End spit out csv files
			//
			//

			//
			// for diagnostics
			//
			double_v2D diagcopyPctScatteredAbove{ pctScatteredAbove };
			//
			//
			//

			ParticleList scattered;
			ParticleList reflected;

			ratio Aratio_ion_bslevel{ std::sqrt(eom.ionsph.B.at(level) / eom.B_ion) };

			//lambdas that generate scatter/reflect percentage and updates total scatter % of defined particle
			auto scatterPercent = [&](percent& sumScatteredAbove, const eV E, const degrees pitch)
			{
				if (pitch > 90.0) throw logic_error("bsAtLevel::scatterPercent: scattering pitch is > 90.0 degrees");

				percent scatter{ 0.0 };

				for (size_t species = 0; species < eom.ionsph.p.size(); species++)
				{
					scatter += scatterPct(sumScatteredAbove, eom.ionsph.Z.at(species),
						eom.ionsph.p.at(species).at(level),	eom.ionsph.h.at(level), E, pitch);
				}

				if (scatter + sumScatteredAbove > 1.0)
				{
					scatter = 1.0 - sumScatteredAbove;
					sumScatteredAbove = 1.0;
				}
				else
					sumScatteredAbove += scatter;

				if (scatter < 0.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is < 0.0 - sct%, %collideAbove: "
					+ to_string(scatter) + ", " + to_string(sumScatteredAbove));
				if (scatter > 1.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is > 1.0 - sct%, %collideAbove: "
					+ to_string(scatter) + ", " + to_string(sumScatteredAbove));

				return scatter;
			};

			auto reflectPercent = [&](percent& sumScatteredAbove, const degrees pitch)
			{
				if (pitch < 90.0) throw logic_error("bsAtLevel::reflectPercent: reflected particle pitch is downgoing");

				percent reflected{ 1.0 - sumScatteredAbove };
				sumScatteredAbove = 1.0;
				return reflected;
			};

			// >> >> adjust PA, check if reflected
			for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
			{
				if (eom.distbins.PA.at(ang) > 90.0) continue; //there should be no upgoing in dNionsphTop

				degrees pa_level{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level)) };
				degrees pa_nextlevel{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level + 1)) };

				if (pa_level > 90.0 || pa_nextlevel > 90.0)
					throw logic_error("ionosphere::multiLevelBS::bsAtLevel: downgoing particle ended up with pitch > 90.0 - PA_bin, PA_level, level: "
					+ to_string(eom.distbins.PA.at(ang)) + ", " + to_string(pa_level) + ", " + to_string(level));
				
				if (pa_level < 0.0)
				{ //particle reflects before this level
					continue;
				}
				else if (pa_nextlevel < 0.0)
				{ //particle reflects before next level, add all remaining particles of this pitch moving in the opposite direction
					//
					//
					// this code has been validated: treating the bottom of the ionosphere as a hard boundary,
					// the graph matches the graph produced with scattering turned off
					//
					//
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (dNionsphTop.binnedData.at(ang).at(eny) != 0.0 && pctScatteredAbove.at(ang).at(eny) < 1.0)
						{ //these particles are reflected and end up upgoing before the bottom of the layer
							if (pa_level >= 90.0) //makes sure the original PA is not upgoing - which would indicate a particle that has reflected in lower levels
								throw logic_error("ionosphere::multiLevelBS::bsAtLevel: pa_level > 90.0, logic error somewhere");

							degrees pa_sat{ dEFlux::newPA(180.0 - pa_level, eom.ionsph.B.at(level), eom.B_sat) }; //pitch at satellite - eventually may need to run through sourceToSatellite

							dNflux reflect{ dNionsphTop.binnedData.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* reflectPercent(pctScatteredAbove.at(ang).at(eny), 180.0 - pa_level)
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ -cos(pa_sat * RADS_PER_DEG) }; //cos of pitch at satellite (where detected)

							if (reflect < 0.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: dN calculated is negative");

							//#define count s_pos <- defined above in sourceToSatellite
							//allows me to use an existing vector in ParticleData with a name that makes sense
							//so, particles.count uses particles.s_pos, but renames it to how I want to use it
							reflected.energy.push_back(eom.distbins.E.at(eny));
							reflected.pitch.push_back(180.0 - pa_level);
							reflected.count.push_back(reflect);
						}
					}
				}
				else
				{ //particle makes it to the next level - invoke scattering
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (dNionsphTop.binnedData.at(ang).at(eny) != 0.0 && pctScatteredAbove.at(ang).at(eny) < 1.0)
						{ //these particles are scattered
							dNflux scatter{ dNionsphTop.binnedData.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* scatterPercent(pctScatteredAbove.at(ang).at(eny), eom.distbins.E.at(eny), pa_level) //percent scattered
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ cos(pa_level * RADS_PER_DEG) }; //pitch at lowest level (level of scattering)

							if (scatter < 0.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: dN calculated is negative");

							scattered.energy.push_back(eom.distbins.E.at(eny));
							scattered.pitch.push_back(pa_level);
							scattered.count.push_back(scatter);
						}
					}
				}
			}	

			TESTVEC_NOTNEGWHOLEVEC(reflected.count, "bsAtLevel:: " + to_string(level) + " :: reflected.count");
			TESTVEC_NOTNEGWHOLEVEC(scattered.count, "bsAtLevel:: " + to_string(level) + " :: scattered.count");

			// >> >> bin particles at level
			ParticlesBinned<dNflux> scatBinned{ eom.distbins.binParticleList(scattered, scattered.count) }; //downgoing
			ParticlesBinned<dNflux> reflBinned{ eom.distbins.binParticleList(reflected, reflected.count) }; //upgoing

			TESTVEC_ISZEROFRSTHALF(scatBinned.binnedData, "bsAtLevel:: " + to_string(level) + " ::scatBinned");
			TESTVEC_ISZEROLASTHALF(reflBinned.binnedData, "bsAtLevel:: " + to_string(level) + " ::reflBinned");
			
			// >> >> calculate backscatter
			ParticlesBinned<dNflux> ret{ backscat::downwardToBackscatter(eom.distbins, scatBinned) };

			//
			//
			// Spit out csv files (for diagnostics)
			debug::outputVectorsToCSV("debug/level sctpct sum/sumpct_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				pctScatteredAbove, eom.distbins.PA, 1200);
			debug::outputVectorsToCSV("debug/level dN/bksc_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				ret.binnedData, eom.distbins.PA, 600, 0, 18000);
			debug::outputVectorsToCSV("debug/level dN/scat_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				scatBinned.binnedData, eom.distbins.PA, 600, 18000, 36000);
			debug::outputVectorsToCSV("debug/level dN/refl_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				reflBinned.binnedData, eom.distbins.PA, 1, 17500, 18000); //should be very close to 90 degrees
			
			if (level == 0) cout << "\n\n================\n\nNeed constructor to force ParticlesBinned to take a Bins instance\n\n================\n\n";

			{
				ParticleList refldiag;
				ParticleList scatdiag;
				ParticleList upwddiag;
				ParticleList downdiag;
				ParticleList initial_ionsph_top;
				
				for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
				{
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						refldiag.energy.push_back(reflBinned.bins.E.at(eny));
						refldiag.pitch.push_back(reflBinned.bins.PA.at(ang));
						refldiag.count.push_back(reflBinned.binnedData.at(ang).at(eny) * reflBinned.bins.E.at(eny));

						scatdiag.energy.push_back(scatBinned.bins.E.at(eny));
						scatdiag.pitch.push_back(scatBinned.bins.PA.at(ang));
						scatdiag.count.push_back(scatBinned.binnedData.at(ang).at(eny) * scatBinned.bins.E.at(eny));
						
						upwddiag.energy.push_back(eom.distbins.E.at(eny));
						upwddiag.pitch.push_back(eom.distbins.PA.at(ang));
						upwddiag.count.push_back((reflBinned.binnedData.at(ang).at(eny) + ret.binnedData.at(ang).at(eny))
							* eom.distbins.E.at(eny));// / abs(cos(eom.distbins.PA.at(ang) * RADS_PER_DEG)));

						if (eom.distbins.PA.at(ang) < 90.0)
						{
							if (reflBinned.binnedData.at(ang).at(eny) != 0.0) throw logic_error("refldiag");
							if (ret.binnedData.at(ang).at(eny) != 0.0) throw logic_error("upwddiag");
						}
						else if (eom.distbins.PA.at(ang) > 90.0)
						{
							if (scatBinned.binnedData.at(ang).at(eny) != 0.0) throw logic_error("scatdiag");
							if (dNionsphTop.binnedData.at(ang).at(eny) != 0.0) throw logic_error("dNionsphTop");
						}

						if (eom.distbins.PA.at(ang) > 90.0) continue; //no downward particles above 90
						degrees downwardNewPA{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level)) };

						if (downwardNewPA < 0.0) continue;
						if (diagcopyPctScatteredAbove.at(ang).at(eny) > (1.0 - FLT_EPSILON)) continue;
						
						downdiag.energy.push_back(dNionsphTop.bins.E.at(eny));
						downdiag.pitch.push_back(downwardNewPA);
						downdiag.count.push_back((1.0 - diagcopyPctScatteredAbove.at(ang).at(eny))
							* dNionsphTop.binnedData.at(ang).at(eny)* dNionsphTop.bins.E.at(eny));// / abs(cos(eom.distbins.PA.at(ang) * RADS_PER_DEG)));
						
						if (level == 0)
						{
							initial_ionsph_top.energy.push_back(eom.distbins.E.at(eny));
							initial_ionsph_top.pitch.push_back(eom.distbins.PA.at(ang));
							initial_ionsph_top.count.push_back(dNionsphTop.binnedData.at(ang).at(eny) * eom.distbins.E.at(eny));
						}
					}
				}

				Bins diagbins(EBINS_DIAG, PABINS_DIAG);

				ParticlesBinned<dNflux> reflbin{ diagbins.binParticleList(refldiag, refldiag.count) };
				ParticlesBinned<dNflux> scatbin{ diagbins.binParticleList(scatdiag, scatdiag.count) };
				ParticlesBinned<dNflux> downbin{ diagbins.binParticleList(downdiag, downdiag.count) };
				ParticlesBinned<dNflux> upwdbin{ diagbins.binParticleList(upwddiag, upwddiag.count) };

				if (level == 0)
				{
					cout << "_ionsph_top.csv\n";
					ParticlesBinned<dNflux> top_ionsph_downward{ diagbins.binParticleList(initial_ionsph_top, initial_ionsph_top.count) };
					debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/_ionsph_top.csv", top_ionsph_downward.binnedData, top_ionsph_downward.bins.PA);
					cout << "complete\n";
				}

				debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/dn_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
					downbin.binnedData, downbin.bins.PA);
				debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/up_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
					upwdbin.binnedData, upwdbin.bins.PA);
				debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/scat_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
					scatbin.binnedData, scatbin.bins.PA);
				debug::outputVectorsToCSV("debug/_forjohn/dEvsPA/refl_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
					reflbin.binnedData, reflbin.bins.PA);
			}
			
			// End spit out csv files
			//
			//

			TESTVEC_ISZEROLASTHALF(ret.binnedData, "bsAtLevel:: " + to_string(level) + " ::backscatter_level");
			
			//ret += reflBinned;
			for (size_t ang = 0; ang < ret.binnedData.size(); ang++)
				for (size_t eny = 0; eny < ret.binnedData.at(0).size(); eny++)
					ret.binnedData.at(ang).at(eny) += reflBinned.binnedData.at(ang).at(eny);

			return ret;
			cout << "Warning: returing reflected without backscatter in bsAtLevel (near line 867)\n";
			return reflBinned;
			#undef count
		}

		DLLEXP percent scatterPct(percent sumCollideAbove, double Z, percm3 p, cm h, eV E, degrees PA)
		{
			if (sumCollideAbove >= 1.0) throw logic_error("ionosphere::multiLevelBS::scatterPct: sumCollideAbove is greater than / equal to 1.0.  "
				+ string("100% has scattered already.  Conditions in bsAtLevel should have prevented this from happening."));

			if (PA > 90.0) //for now there shouldn't be any scatter as particles move upward
				return 0.0;

			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	} //end namespace multiLevelBS

} //end namespace ionosphere
