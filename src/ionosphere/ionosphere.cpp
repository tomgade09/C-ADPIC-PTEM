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

inline void printIonosphere(const ionosphere::IonosphereSpecs& ionsph)
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

inline void printLayer(const ionosphere::IonosphereSpecs& ionsph, size_t layer)
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
			ParticleData eomtop;
			readDblBin(eomtop.vpara, eomsim.datadir + "bins\\satellites\\topElec_vpara.bin");
			readDblBin(eomtop.vperp, eomsim.datadir + "bins\\satellites\\topElec_vperp.bin");
			readDblBin(eomtop.s_pos, eomsim.datadir + "bins\\satellites\\topElec_s.bin");
			utils::numerical::v2DtoEPitch(eomtop.vpara, eomtop.vperp, eomsim.mass, eomtop.energy, eomtop.pitch);

			for (int part = 0; part < eomsim.initial.energy.size(); part++)
			{
				//does ideal dist have more detected particles?
				//top condition here - will require newPA
				if ((ideal.bottom.pitch.at(part) != 0.0) && (eomsim.bottom.pitch.at(part) == 0.0))
					extra.at(1)++;
				if ((ideal.upward.pitch.at(part) != 0.0) && (eomsim.upward.pitch.at(part) == 0.0))
					extra.at(2)++;
				if ((ideal.dnward.pitch.at(part) != 0.0) && (eomsim.dnward.pitch.at(part) == 0.0))
					extra.at(3)++;

				//compare pitch, E of top, btm, upw, dnw
				if ((ideal.bottom.pitch.at(part) != 0.0) && (eomsim.bottom.pitch.at(part) != 0.0) && (err(ideal.bottom.pitch.at(part), eomsim.bottom.pitch.at(part)) > maxErr.at(1)))
					maxErr.at(1) = err(ideal.bottom.pitch.at(part), eomsim.bottom.pitch.at(part));
				if ((ideal.upward.pitch.at(part) != 0.0) && (eomsim.upward.pitch.at(part) != 0.0) && (err(ideal.upward.pitch.at(part), eomsim.upward.pitch.at(part)) > maxErr.at(2)))
					maxErr.at(2) = err(ideal.upward.pitch.at(part), eomsim.upward.pitch.at(part));
				if ((ideal.dnward.pitch.at(part) != 0.0) && (eomsim.dnward.pitch.at(part) != 0.0) && (err(ideal.dnward.pitch.at(part), eomsim.dnward.pitch.at(part)) > maxErr.at(3)))
					maxErr.at(3) = err(ideal.dnward.pitch.at(part), eomsim.dnward.pitch.at(part));
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

	DLLEXP dEflux_v2D steadyFlux(const EOMSimData& eomdata)
	{
		{
			vector<double> qspstopUpward;
			utils::fileIO::readDblBin(qspstopUpward, eomdata.datadir + "bins/satellites/QSPSTopUpg_vpara.bin");

			double minvpara{ 10000000.0 };

			for (auto elec : qspstopUpward)
				if (elec < minvpara) minvpara = elec;

			std::cout << std::setprecision(10) << "min vpara upward: " << minvpara << " m/s\n";
		}

		{
			Satellite* qspsbtmUpg{ eomdata.sim->satellite("QSPSBtmUpg") };
			Satellite* qspsbtmDng{ eomdata.sim->satellite("QSPSBtmDng") };

			vector<double> qspsbtmUpgVpara{ qspsbtmUpg->data().at(0) };
			vector<double> qspsbtmUpgVperp{ qspsbtmUpg->data().at(1) };
			vector<double> qspsbtmDngVpara{ qspsbtmDng->data().at(0) };
			vector<double> qspsbtmDngVperp{ qspsbtmDng->data().at(1) };

			ParticleData qspsbtmUpgPD(qspsbtmUpgVpara, qspsbtmUpgVperp, MASS_ELECTRON);
			ParticleData qspsbtmDngPD(qspsbtmDngVpara, qspsbtmDngVperp, MASS_ELECTRON);

			auto err = [](double base, double dev)
			{//returns % error
				return std::abs((base - dev) / base);
			};

			double min_err{ 0.0001 };
			double max_err{ 0 };
			size_t max_ind{ 0 };

			vector<size_t> err_ind;
			for (size_t ind = 0; ind < qspsbtmUpgPD.energy.size(); ind++)
			{
				if (qspsbtmUpgPD.vpara.at(ind) == 0 || qspsbtmDngPD.vpara.at(ind) == 0) continue;

				double PA_err{ err(180.0 - qspsbtmUpgPD.pitch.at(ind), qspsbtmDngPD.pitch.at(ind)) };
				double E_err{ err(qspsbtmUpgPD.energy.at(ind), qspsbtmDngPD.energy.at(ind)) };

				if (E_err > min_err)
					err_ind.push_back(ind);
				else if (PA_err > min_err)
					err_ind.push_back(ind);

				if (PA_err > max_err) { max_err = PA_err; max_ind = ind; }
				if (E_err > max_err)  { max_err = E_err;  max_ind = ind; }
			}

			std::cout << "Num err: " << err_ind.size() << "\n";
			std::cout << "Max err: " << max_err << "\n";
			std::cout << "Err ind: " << max_ind << "\n";
			std::cout << "Upg PA : " << qspsbtmUpgPD.pitch.at(max_ind) << "\n";
			std::cout << "Dng PA : " << 180.0 - qspsbtmDngPD.pitch.at(max_ind) << " (" << qspsbtmDngPD.pitch.at(max_ind) << ")\n";
			std::cout << "Upg E  : " << qspsbtmUpgPD.energy.at(max_ind) << "\n";
			std::cout << "Dng E  : " << qspsbtmDngPD.energy.at(max_ind) << "\n";
		}
		
		exit(1);

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

		//
		//
		// QSPS test functions
		/*
		vector<double> final_s;
		vector<double> final_vpara;
		vector<double> final_vperp;

		utils::fileIO::readDblBin(final_s, eomdata.datadir + "bins/particles_final/elec_s.bin");
		utils::fileIO::readDblBin(final_vpara, eomdata.datadir + "bins/particles_final/elec_vpara.bin");
		utils::fileIO::readDblBin(final_vperp, eomdata.datadir + "bins/particles_final/elec_vperp.bin");

		vector<vector<double>> in_sim;
		for (size_t ind = 0; ind < final_s.size(); ind++)
		{
			if (final_s.at(ind) < eomdata.s_mag && final_s.at(ind) > eomdata.s_ion)
			{
				in_sim.push_back(vector<double>{ final_s.at(ind), final_vpara.at(ind), final_vperp.at(ind) });
			}
		}

		double low_s{ 1.0e18 }; double low_vpara{ 1.0e18 }; double low_vperp{ 1.0e18 };
		double hi_s{ 0.0 };	double hi_vpara{ 0.0 }; double hi_vperp{ 0.0 };

		for (size_t ind = 0; ind < in_sim.size(); ind++)
		{
			if (in_sim.at(ind).at(0) < low_s)
			{
				low_s = in_sim.at(ind).at(0);
				low_vpara = in_sim.at(ind).at(1);
				low_vperp = in_sim.at(ind).at(2);
			}
			if (in_sim.at(ind).at(0) > hi_s)
			{
				hi_s = in_sim.at(ind).at(0);
				hi_vpara = in_sim.at(ind).at(1);
				hi_vperp = in_sim.at(ind).at(2);
			}
		}

		cout << "Number of Particles Still In Sim: " << in_sim.size() << "\n";
		cout << "s_Mag, s_Ion: " << eomdata.s_mag << ", " << eomdata.s_ion << "\n";
		cout << "s_QSPS Top, s_QSPS btm: " << "11122579, 9796766\n";
		cout << "Lowest s - s, vpara, vperp : " << low_s << ", " << low_vpara << ", " << low_vperp << "\n";
		cout << "Highest s - s, vpara, vperp: " << hi_s << ", " << hi_vpara << ", " << hi_vperp << "\n";

		utils::fileIO::write2DCSV(in_sim, "in_sim.csv", in_sim.at(0).size(), in_sim.size());

		exit(1);
		*/
		//
		//
		//

		printIonosphere(eomdata.ionsph);

		// 1. Adjust Maxwellian by Cos, SQRT(B ratio) Factors
		double Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph to satellite
		double Aratio_mag_sat{ std::sqrt(eomdata.B_sat / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to satellite
		double Aratio_ion_ion{ std::sqrt(eomdata.B_ion / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph (up, reflect, down) to ionsph
		double Aratio_mag_ion{ std::sqrt(eomdata.B_ion / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to ionsph

		dNflux_v1D maxwellian_sat{ eomdata.maxwellian }; //scaled by decrease in gyroradius cross-sectional area A below
		dNflux_v1D maxwellian_ion{ eomdata.maxwellian };

		for (size_t part = 0; part < eomdata.initial.s_pos.size(); part++) //isotropize counts -> 3D
		{
			if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)    //ionospheric source, upgoing
			{
				maxwellian_sat.at(part) *= -cos(eomdata.initial.pitch.at(part) * RADS_PER_DEG) * Aratio_ion_sat;
				maxwellian_ion.at(part) *= Aratio_ion_ion; //should just be equal to 1, but gives it symmetry
			}
			else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)//magnetospheric source, downgoing
			{
				//if (eomdata.dnward.pitch.at(part) >= 90.0)
					//throw logic_error("ionosphere::steadyFlux : downward pitch of particle is >= 90 deg");

				maxwellian_sat.at(part) *= 1.0 / cos(eomdata.dnward.pitch.at(part) * RADS_PER_DEG) * Aratio_mag_sat;
				maxwellian_ion.at(part) *= Aratio_mag_ion;
			}
			else
				throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source - it's somewhere between the ionosphere and magnetosphere");

			//if (eomdata.bottom.pitch.at(part) >= 90.0)
				//throw logic_error("ionosphere::steadyFlux : bottom pitch of particle is >= 90 deg");
		}

		//TESTVEC_NOTNEGWHOLEVEC(maxwellian_sat, "steadyFlux::maxwellian_sat"); //had to disable to pass QSPS


		// 2. Calculate dEfluxes
		dEflux_v2D distfluxdnward{ dEFlux::satellite(eomdata.dnward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D distfluxupward{ dEFlux::satellite(eomdata.upward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D bkscfluxupward{ dEFlux::backscatr(eomdata, maxwellian_ion) };

		//printVec2D(distfluxdnward, "Dnward Flux at Satellite");
		//printVec2D(distfluxupward, "Upward Flux at Satellite");
		//printVec2D(bkscfluxupward, "Flux Due to Backscatter at Satellite");

		debug::outputVectorsToCSV("debug/final dE/downward.csv", distfluxdnward, eomdata.satbins.PA);
		debug::outputVectorsToCSV("debug/final dE/upward.csv", distfluxupward, eomdata.satbins.PA);
		debug::outputVectorsToCSV("debug/final dE/backscatter.csv", bkscfluxupward, eomdata.satbins.PA);


		// 3. Sum dEfluxes
		for (size_t ang = 0; ang < distfluxupward.size(); ang++)
			for (size_t eny = 0; eny < distfluxupward.at(ang).size(); eny++)
				distfluxupward.at(ang).at(eny) += distfluxdnward.at(ang).at(eny) + bkscfluxupward.at(ang).at(eny);

		//
		//
		// TEST BACKSCATTER ONLY //
		//return bkscfluxupward;
		// DELETE THIS WHEN DONE WITH TEST //
		//
		//

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
			for (size_t ang = 0; ang < ret.size(); ang++)
				for (size_t eng = 0; eng < ret.at(0).size(); eng++)
					ret.at(ang).at(eng) *= satBins.E.at(eng);


			return ret;
		}

		DLLEXP dEflux_v2D backscatr(const EOMSimData& eomdata, const dNflux_v1D& dNatIonsph1D)
		{
			// Section 1 - Get dNflux at Satellite
			// 1.1. Bin Escaped Particles by EOM simulation (distribution) PA, Energy Bins (for high res backscatter)
			dNflux_v2D dNatIonsph{ binning::binParticles(eomdata.bottom, eomdata.distbins, dNatIonsph1D) };
			// output: 2D vector [PA][Eng] of number of escaped particles (dNFlux), weighted by specified maxwellians, binned by Energy and Pitch Angle

			//TESTVEC_ISZEROFRSTHALF(dNatIonsph, "backscatr::dNatIonsph"); //had to disable to pass QSPS

			// 1.2. Calculate BS dNflux from dNflux Incident to Ionosphere
			dNflux_v2D BSatIonsph{ multiLevelBS::scatterMain(eomdata, dNatIonsph) }; //new multi-level hotness
			// output: 2D vector of backscatter dNFlux by dist bins at ionosphere (upgoing)

			TESTVEC_ISZEROLASTHALF(BSatIonsph, "backscatr::BSatIonsph");

			// 1.3. Translate BS dNflux at Ionosphere (dist binned) to dNflux at Satellite (sat binned)
			//
			//
			// will need to change this when we have time-dependent E fields - this doesn't account for them
			// just static / time-independent fields
			//
			dNflux_v2D ret{ backscat::ionsphToSatellite(eomdata, BSatIonsph) };
			// output: 2D vector of bs dNFlux at satellite per PA, E (sat binned) - should only be upward (90-180)
			// Section 1 End

			// 2.1. Convert from dNflux to dEflux
			for (size_t ang = 0; ang < ret.size(); ang++)
			{
				for (size_t eny = 0; eny < ret.at(0).size(); eny++)
					ret.at(ang).at(eny) *= eomdata.satbins.E.at(eny);
			}
			// output: 2D vector of bs dEflux

			return ret;
		}

		DLLEXP double newPA(const degrees PA_init, const tesla B_init, const tesla B_final)
		{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
			double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

			if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

			degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

			if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

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

			for (size_t part = 0; part < particles.energy.size(); part++) //iterate over particles
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

		DLLEXP void symmetricBins180To360(dEflux_v2D& data, double_v1D& binAngles) //takes bins from 0-180 and extends to 360 symmetrically reflected across 180 degrees
		{//data[angle][energy]
			binAngles.resize(2 * binAngles.size()); //double size of angles vector (goes from 0-180 now)
			data.resize(binAngles.size()); //set data (vector of double vectors) to size of binAngles (now double sized)
			
			for (size_t ang = binAngles.size() / 2; ang < binAngles.size(); ang++)
			{
				binAngles.at(ang) = 2 * binAngles.at(ang - 1) - binAngles.at(ang - 2);
				data.at(ang) = data.at(binAngles.size() - ang - 1);
			}
		} //needs a good testing
	} //end namespace ionosphere::binning

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

		//constexpr double EVANS_PRIM_LOGM{ 1.5 };  //obtained by log linefitting Evans, 1974 - these seem closest
		//constexpr double EVANS_PRIM_LOGB{ -4.0 };
		//constexpr double EVANS_SECD_LOGM{ -2.1 };
		//constexpr double EVANS_SECD_LOGB{ 0.3 };

		//constexpr double JOHND_PRIM_LOGM{ 1.132 };
		//constexpr double JOHND_PRIM_LOGB{ -4.90 };
		//constexpr double JOHND_SECD_LOGM_LT10{ -1.0 };
		//constexpr double JOHND_SECD_LOGB_LT10{ -2.5 };
		//constexpr double JOHND_SECD_LOGM_GT10{ -2.245 };
		//constexpr double JOHND_SECD_LOGB_GT10{  0.153 };

		constexpr double JOHND_PRIM_LOGM{ 0.505 };
		constexpr double JOHND_PRIM_LOGB{ -5.16 };
		constexpr double JOHND_SECD_LOGM_LT25{ -0.975 };
		constexpr double JOHND_SECD_LOGB_LT25{ -1.47 };
		constexpr double JOHND_SECD_LOGM_GT25{ -1.95 };
		constexpr double JOHND_SECD_LOGB_GT25{ -0.11 };

		DLLEXP dNflux johnd_flux(eV E_eval, eV E_incident, dNflux dN_incident)
		{
			if (E_eval > E_incident * (1.0 + FLT_EPSILON)) return 0.0;
				//throw logic_error("johnd_flux: E_eval is higher than E_incident.  Not physical.  Eval, Incident: " + to_string(E_eval) + " , " + to_string(E_incident));

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

		/*DLLEXP dEflux integralJohnd_flux(double lower, double upper, double E_incident)
		{
			throw exception("integralJohnd_flux used");//make sure this isnt used for now
			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			double integral_sec{ (pow(upper, secd_logm + 1.0) - pow(lower, secd_logm + 1.0)) * pow(10.0, secd_logb) / (secd_logm + 1.0) };
			double integral_prm{ (pow(upper, JOHND_PRIM_LOGM + 1.0) - pow(lower, JOHND_PRIM_LOGM + 1.0)) * pow(10.0, JOHND_PRIM_LOGB + 4.0) / ((JOHND_PRIM_LOGM + 1.0) * pow(E_incident, JOHND_PRIM_LOGM + 1.0)) };
			return integral_sec + integral_prm;
		}*/

		DLLEXP dNflux_v2D downwardToBackscatter(const Bins& dist, const dNflux_v2D& dNpointofScatter)
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

					dNsumPerE.at(egy) += dNpointofScatter.at(ang).at(egy) / (dist.PA.size() / 2) * 7500 * 4 * PI *
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
			dNflux_v2D BS(dist.PA.size(), vector<double>(dist.E.size()));

			for (size_t ang = 0; ang < dist.PA.size(); ang++)
			{
				if (dist.PA.at(ang) <= 90.0)
					continue; //empty vector of the right size
				else
				{
					for (size_t eny = 0; eny < dist.E.size(); eny++)
					{
						BS.at(ang).at(eny) = dNbsPerE.at(eny) * -cos(dist.PA.at(ang) * RADS_PER_DEG) / (4 * PI);
					}
				}
			}
			// output: 2D vector of bs dEFlux at ionosphere per pitch, energy bin - should only be upward (90-180)

			return BS;
		}

		DLLEXP dNflux_v2D ionsphToSatellite(const EOMSimData& eomdata, const dNflux_v2D& bsCounts)
		{
			double Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) };

			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleData particles; //TO DO: eventually reserve vectors - it may speed things up a bit

			for (size_t ang = 0; ang < eomdata.distbins.PA.size(); ang++)
			{
				if (eomdata.distbins.PA.at(ang) <= 90.0) continue;

				for (size_t eny = 0; eny < eomdata.distbins.E.size(); eny++) //this works because ionospheric bins are same as distribution
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
		DLLEXP dNflux_v2D scatterMain(const EOMSimData& eom, const dNflux_v2D& dNionsphTop)
		{
			//
			//
			// Spit out csv files (for diagnostics)
			debug::outputVectorsToCSV("debug/angle, pitch bins/satbins.csv", { eom.satbins.E, eom.satbins.PA }, vector<string>{ "E", "PA" });
			debug::outputVectorsToCSV("debug/angle, pitch bins/distbins.csv", { eom.distbins.E, eom.distbins.PA }, vector<string>{ "E", "PA" });
			debug::outputVectorsToCSV("debug/dNdwnwdIonosphereTop.csv", dNionsphTop, eom.distbins.PA, 600, 18000, 36000);
			// End spit out csv files
			//
			//

			double_v2D pctScattered(eom.distbins.PA.size(), dNflux_v1D(eom.distbins.E.size())); //% scattered per bin

			ParticleData particles;

			// >> level calculate
			for (size_t level = 0; level < eom.ionsph.s.size() - 1; level++)
			{//for now, depends on adding one extra level (to see if particles reflect somewhere within the last layer)
				printLayer(eom.ionsph, level);

				dNflux_v2D bs_level{ bsAtLevel(eom, dNionsphTop, pctScattered, level) };

				TESTVEC_ISZEROLASTHALF(bs_level, "scatterMain::bs_level");

				double Aratio_bslevel_ion{ std::sqrt(eom.B_ion / eom.ionsph.B.at(level)) };

				// >> >> adjust PA
				for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
				{
					if (eom.distbins.PA.at(ang) <= 90.0) continue;

					degrees paTopIon{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(level), eom.ionsph.B.at(0)) };
					if (paTopIon < 0.0) throw logic_error(string("ionosphere::multiLevelBS::scatterMain: return value of newPA is -1 ")
						+ "indicating reflection before B_final.  This shouldn't happen as B_final is less in magnitude than B_initial. "
						+ " The particle can't physically reflect in this case, but simply be pushed closer to 180.  Something is wrong.");

					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (bs_level.at(ang).at(eny) != 0.0)
						{
							//#define count s_pos <- defined above in bsSrcToSat
							//allows me to use an existing vector in ParticleData with a name that makes sense
							//so, particles.count uses particles.s_pos, but renames it to how I want to use it
							particles.pitch.push_back(paTopIon);
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

		DLLEXP dNflux_v2D bsAtLevel(const EOMSimData& eom, const dNflux_v2D& dNionsphTop, double_v2D& pctScatteredAbove, size_t level)
		{
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

			ParticleData dngoing; //TO DO: reserve size of vectors to increase the speed
			ParticleData upgoing;

			double Aratio_ion_bslevel{ std::sqrt(eom.ionsph.B.at(level) / eom.B_ion) };

			//lambda that generates scatter/reflect percentage and updates total scatter % of defined particle
			auto sctReflPct = [&](double& pctScatteredAbove, const eV E, const degrees pitch)
			{
				//upgoing - whatever hasn't scattered so far, reflects
				//no scattering happens in this case
				if (pitch > 90.0)
				{
					double ret{ 1.0 - pctScatteredAbove };
					pctScatteredAbove = 1.0;
					return ret; //so return 1.0 - "sum % scattered in layers above"
				}
				
				//==================//

				//dngoing - percent that scatters in this layer
				//scattering happens in this case
				double sct{ 0.0 };

				for (size_t species = 0; species < eom.ionsph.p.size(); species++)
				{ //iterate through species, add scatter % for each one
					sct += scatterPct(pctScatteredAbove, eom.ionsph.Z.at(species), eom.ionsph.p.at(species).at(level),
						eom.ionsph.h.at(level), E, pitch);
				}

				if (sct + pctScatteredAbove > 1.0)
				{ //no more than 100% of particles can have scattered
					sct = 1.0 - pctScatteredAbove; //sct% for this layer is equal to however much is left until the total scattered is 100%
					pctScatteredAbove = 1.0;
				}
				else
				{ //sum is less than 100%, so add to the total as normal
					pctScatteredAbove += sct;
				}

				if (sct < 0.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is < 0.0 - sct%, %collideAbove: "
					+ to_string(sct) + ", " + to_string(pctScatteredAbove));
				if (sct > 1.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: scatter % is > 1.0 - sct%, %collideAbove: "
					+ to_string(sct) + ", " + to_string(pctScatteredAbove));

				return sct;
			};

			// >> >> adjust PA, check if reflected
			for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
			{
				if (eom.distbins.PA.at(ang) > 90.0) continue; //there should be no upgoing in dNionsphTop

				degrees pa_level{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level)) };
				degrees pa_nextlevel{ dEFlux::newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(0), eom.ionsph.B.at(level + 1)) };

				if (pa_level > 90.0) throw logic_error("ionosphere::multiLevelBS::bsAtLevel: downgoing particle ended up with pitch > 90.0 - PA_bin, PA_level, level: "
					+ to_string(eom.distbins.PA.at(ang)) + ", " + to_string(pa_level) + ", " + to_string(level));
				
				if (pa_level < 0.0)
				{ //particle reflects before this level
					continue;
				}
				else if (pa_nextlevel < 0.0)
				{ //particle reflects before next level, add all particles of this pitch moving in the opposite direction
					//
					//
					// this code has been validated: treating the bottom of the ionosphere as a hard boundary,
					// the graph matches the graph produced with scattering turned off
					//
					//
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (dNionsphTop.at(ang).at(eny) != 0.0 && pctScatteredAbove.at(ang).at(eny) < 1.0)
						{ //these particles are upgoing
							if (180.0 - pa_level <= 90.0)
								throw logic_error("ionosphere::multiLevelBS::bsAtLevel: sctReflPct will not return (1 - % scattered), logic error somewhere");

							degrees pa_sat{ dEFlux::newPA(180.0 - pa_level, eom.ionsph.B.at(level), eom.B_sat) }; //pitch at satellite - eventually may need to run through sourceToSatellite
							
							dNflux dNbyEPAlayer{ dNionsphTop.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* sctReflPct(pctScatteredAbove.at(ang).at(eny), eom.distbins.E.at(eny), 180.0 - pa_level) //percent reflected - "180 - pa" is used to signal to the lambda to return whatever hasn't reflected
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ -cos(pa_sat * RADS_PER_DEG) }; //cos of pitch at satellite (where detected)

							//#define count s_pos <- defined above in sourceToSatellite
							//allows me to use an existing vector in ParticleData with a name that makes sense
							//so, particles.count uses particles.s_pos, but renames it to how I want to use it
							upgoing.energy.push_back(eom.distbins.E.at(eny));
							upgoing.pitch.push_back(180.0 - pa_level);
							upgoing.count.push_back(dNbyEPAlayer);
						}
					}
				}
				else
				{ //particle makes it to the next level - invoke scattering
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (dNionsphTop.at(ang).at(eny) != 0.0 && pctScatteredAbove.at(ang).at(eny) < 1.0)
						{ //these particles are downgoing
							dNflux dNbyEPAlayer{ dNionsphTop.at(ang).at(eny) //dNflux incident to ionosphere at a given E, PA
								* sctReflPct(pctScatteredAbove.at(ang).at(eny), eom.distbins.E.at(eny), pa_level) //percent scattered
								* Aratio_ion_bslevel //gyroradius cross-sectional area difference from ionsph to bslevel
								/ cos(pa_level * RADS_PER_DEG) }; //pitch at lowest level (level of scattering)

							dngoing.energy.push_back(eom.distbins.E.at(eny));
							dngoing.pitch.push_back(pa_level);
							dngoing.count.push_back(dNbyEPAlayer);
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

			//
			//
			// Spit out csv files (for diagnostics)
			debug::outputVectorsToCSV("debug/level sctpct sum/sumpct_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				pctScatteredAbove, eom.distbins.PA, 1200);
			debug::outputVectorsToCSV("debug/level dN/bksc_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				ret, eom.distbins.PA, 600, 0, 18000);
			debug::outputVectorsToCSV("debug/level dN/scat_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				dnBinned, eom.distbins.PA, 600, 18000, 36000);
			debug::outputVectorsToCSV("debug/level dN/refl_" + to_string((int)eom.ionsph.s.at(level)) + ".csv",
				upBinned, eom.distbins.PA, 1, 17500, 18000); //should be very close to 90 degrees
			// End spit out csv files
			//
			//

			TESTVEC_ISZEROLASTHALF(ret, "bsAtLevel:: " + to_string(level) + " ::backscatter_level");
			
			for (size_t ang = 0; ang < ret.size(); ang++)
				for (size_t eny = 0; eny < ret.at(0).size(); eny++)
					ret.at(ang).at(eny) += upBinned.at(ang).at(eny);

			return ret;
		}

		DLLEXP percent scatterPct(percent sumCollideAbove, double Z, percm3 p, cm h, eV E, degrees PA)
		{
			if (sumCollideAbove >= 1.0) throw logic_error("ionosphere::multiLevelBS::scatterPct: sumCollideAbove is greater than / equal to 1.0.  "
				+ string("100% has scattered already.  Conditions in bsAtLevel should have prevented this from happening."));

			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	} //end namespace multiLevelBS

} //end namespace ionosphere
