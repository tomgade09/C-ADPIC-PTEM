#include "utils/ionosphere.h"
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
using std::logic_error;
using utils::fileIO::readDblBin;
using utils::numerical::generateSpacedValues;

constexpr double EVANS_PRIM_LOGM{ 1.5 };  //obtained by log linefitting Evans, 1974 - these seem closest
constexpr double EVANS_PRIM_LOGB{ -4.0 };
constexpr double EVANS_SECD_LOGM{ -2.1 };
constexpr double EVANS_SECD_LOGB{ 0.3 };

inline void printIonosphere(const ionosphere::IonosphereSpecs& ionsph)
{
	cout << "==================== Backscatter Simulation ====================" << "\n";
	cout << "Min, Max s (m): " << *(ionsph.s.end() - 2) << ", " << ionsph.s.front() << "\n";
	cout << "Min, Max B (T): " << *(ionsph.B.end() - 2) << ", " << ionsph.B.front() << "\n";
	cout << "Num of Layers : " << ionsph.s.size() - 1 << "\n";
	cout << "Atomic Species: " << "Fix later.\n";
	cout << "================================================================" << "\n";
}

inline void printLayer(const ionosphere::IonosphereSpecs& ionsph, unsigned int layer)
{
	cout << "Layer: " << layer << " / " << ionsph.s.size() - 2 << ", s: ";
	cout << ionsph.s.at(layer) << ", B: " << ionsph.B.at(layer) << "\n";
}

namespace ionosphere
{
	namespace debug
	{
		DLLEXP EOMSimData generateIdealSatDists(const EOMSimData& eomdata)
		{
			auto newPA = [](const degrees PA_init, const tesla B_init, const tesla B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

				degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

				return ret;
			};

			ParticleData btm(eomdata.bottom.energy.size());
			ParticleData upw(eomdata.upward.energy.size());
			ParticleData dnw(eomdata.dnward.energy.size());

			for (int part = 0; part < eomdata.initial.energy.size(); part++)
			{
				if (eomdata.initial.s_pos.at(part) < eomdata.s_ion * 1.001)
				{
					double satPA{ newPA(eomdata.initial.pitch.at(part), eomdata.B_ion, eomdata.B_sat) };

					if (satPA > 0)
					{
						upw.pitch.at(part) = satPA;
						upw.energy.at(part) = eomdata.initial.energy.at(part);
					}
				}
				else if (eomdata.initial.s_pos.at(part) > eomdata.s_mag * 0.999)
				{
					double satPA{ newPA(eomdata.initial.pitch.at(part), eomdata.B_mag, eomdata.B_sat) };
					double btmPA{ newPA(eomdata.initial.pitch.at(part), eomdata.B_mag, eomdata.B_ion) };

					if (satPA > 0)
					{
						dnw.pitch.at(part) = satPA;
						dnw.energy.at(part) = eomdata.initial.energy.at(part);
						if (btmPA < 0)
						{
							upw.pitch.at(part) = 180.0 - satPA;
							upw.energy.at(part) = eomdata.initial.energy.at(part);
						}
					}
					if (btmPA > 0)
					{
						btm.pitch.at(part) = btmPA;
						btm.energy.at(part) = eomdata.initial.energy.at(part);
					}
				}
				else
				{
					throw logic_error("debug::generateIdealSatDists : particle is not ionospheric or magnetospheric source");
				}
			}

			EOMSimData ret{ eomdata };
			ret.bottom = btm;
			ret.upward = upw;
			ret.dnward = dnw;

			return ret;
		}
		
		DLLEXP void eomError(const EOMSimData& ideal, const EOMSimData& eomsim)
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

		void setMaxwellians(EOMSimData& eomdata)
		{
			std::vector<double> realMagMaxwellian = {
				511.7803363,
				455.366866,
				404.1659709,
				357.6960119,
				315.5198528,
				277.2407488,
				242.498614,
				210.9666344,
				182.3481933,
				156.3740811,
				132.7999634,
				111.4040819,
				91.7746781,
				73.55799355,
				57.02452115,
				42.01873297,
				34.36415945,
				28.31990953,
				22.99218389,
				18.56252397,
				14.73461663,
				11.99012355,
				9.875617326,
				8.350030687,
				7.817520858,
				7.352563096,
				7.057355075,
				6.732956547,
				5.444957563,
				4.258616122,
				3.148063762,
				2.667671894,
				2.31874741,
				2.413216721,
				2.510346733,
				2.16919973,
				1.822258224,
				1.538644415,
				1.454668358,
				1.448302276,
				1.421289335,
				1.392400846,
				1.356434703,
				1.3223143,
				1.340996193,
				1.24936111,
				1.082697097,
				1.027704468,
				1.022389203,
				0.954603057,
				0.853591162,
				0.787414014,
				0.712480444,
				0.618899297,
				0.613108903,
				0.676524001,
				0.741544197,
				0.809125578,
				0.826634801,
				0.844583081,
				0.909628356,
				0.96415381,
				1.006445782,
				1.033757784,
				1.034400169 * 1.480727,
				1.042600148 * 1.480727,
				1.055748627 * 1.428545,
				1.07152601  * 1.428545,
				1.078133553 * 1.470468,
				1.058734323 * 1.470468,
				1.039323547 * 1.438624,
				1.012305997 * 1.438624,
				0.994205528 * 1.442204,
				0.985836018 * 1.442204,
				0.910594196 * 1.312604,
				0.838235396 * 1.312604,
				0.759978758 * 1.234802,
				0.688727757 * 1.234802,
				0.582535504 * 1.223116,
				0.484989966 * 1.223116,
				0.393631204 * 1.319446,
				0.330308124 * 1.319446,
				0.27732655  * 1.287453,
				0.225898359 * 1.287453,
				0.178965883, //84
				0.142250867,
				0.110109027,
				0.082409802,
				0.060637842,
				0.042555514,
				0.027887484,
				0.0150541,
				0.008638964,
				0.004889727,
				0.002865042,
				0.010868959
			};

			/*
			std::vector<double> realIonMaxwellian = {
				628.8979407 / 2.4786570,
				565.8638197 / 2.4786570,
				508.6540186 / 2.1904623,
				456.7303733 / 2.1904623,
				409.6044458 / 1.9406344,
				366.8329293 / 1.9406344,
				328.0134787 / 1.7238276,
				292.780925  / 1.7238276,
				260.8038409 / 1.5356057,
				231.7814228 / 1.5356057,
				205.4406609 / 1.3725660,
				181.5337716 / 1.3725660,
				159.4832008 / 1.2364910,
				138.7981913 / 1.2364910,
				120.0244659 / 1.1444663,
				102.9854228 / 1.1444663,
				89.92642941 / 1.0706112,
				78.43829219 / 1.0706112,
				67.92556403,
				58.16316132,
				49.00339734,
				39.55476535,
				32.43565065,
				27.49226718,
				23.65374146,
				20.06519194,
				16.08473353,
				12.55060252,
				10.72488715,
				9.128684634,
				7.798530255,
				6.437362853,
				5.176560025,
				4.356319064,
				3.620845049,
				3.394283882,
				3.146182049,
				2.535331134,
				2.166546899,
				1.906599314,
				1.725538745,
				1.572503401,
				1.42385702,
				1.291013895,
				1.33179657,
				1.213164987,
				0.98581798,
				0.891060771,
				0.856747842,
				0.813801222,
				0.767419376,
				0.762158834,
				0.727575643,
				0.647540963,
				0.613472219,
				0.616017952,
				0.584045464,
				0.515672013,
				0.529429898,
				0.548721841,
				0.510828937,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0,
				0
			};
			*/

			std::vector<double> realIonMaxwellian(96, 0);

			for (size_t part = 0; part < eomdata.initial.pitch.size(); part++)
			{
				if (eomdata.initial.pitch.at(part) < 90.0)
				{
					eomdata.maxwellian.at(part) = realMagMaxwellian.at(part % 96);
				}
				else
				{
					eomdata.maxwellian.at(part) = realIonMaxwellian.at(part % 96);
				}
			}
		}

		void outputFunctionsToCSV(const EOMSimData& eom)
		{
			string rootfold{".\\debug\\"};
			cout << rootfold << "\n";

			dNflux_v2D dE_johnd_flux(eom.distbins.E.size(), dNflux_v1D(eom.distbins.E.size() + 1)); //outer: E incident, inner: E eval
			for (size_t Einc = 0; Einc < eom.distbins.E.size(); Einc++)
			{
				for (size_t Eeval = 0; Eeval < eom.distbins.E.size(); Eeval++)
				{ //test when a single particle hits the ionosphere and scatters - this should produce Evans' graph or something like it
					dE_johnd_flux.at(Einc).at(Eeval + 1) = 
						backscat::johnd_flux(eom.distbins.E.at(Eeval), eom.distbins.E.at(Einc), 1.0);
				}
			}
			
			{
				utils::fileIO::CSV bsdEOnePart(rootfold + "bsdEOnePart.csv");
				bsdEOnePart.add(dE_johnd_flux, vector<string>(eom.distbins.E.size()));
			}
		}
	}

	DLLEXP dEflux_v2D steadyFlux(const EOMSimData& eomdata)
	{
		/*{ //print maxwellian values for magsph and ionsph
			std::vector<double> print;
			for (int iii = 0; iii < 96; iii++)
				print.push_back(eomdata.maxwellian.at(iii));
			printVec2D({ print }, "Maxwellian ionsph");

			std::vector<double> print2;
			for (int iii = 1728000; iii < 1728096; iii++)
				print2.push_back(eomdata.maxwellian.at(iii));
			printVec2D({ print2 }, "Maxwellian magsph");
		}*/

		//std::cout << "B_ion: " << eom.B_ion << "\n";
		//std::cout << "B_sat: " << eom.B_sat << "\n";
		//std::cout << "B_mag: " << eom.B_mag << "\n";

		//debug::outputFunctionsToCSV(eom);

		//exit(1);

		//EOMSimData eomdata{ debug::generateIdealSatDists(eom) };
		//debug::setMaxwellians(eomdata);

		//
		//
		// Spit out csv files (for diagnostics)
		{ //outputs densities of ionosphere species by height to a csv
			utils::fileIO::CSV dens{ "debug/ionospheric densities/densities.csv" };

			dens.add(eomdata.ionsph.p.at(0), "N2 dens");
			dens.add(eomdata.ionsph.p.at(1), "He dens");
			dens.add(eomdata.ionsph.p.at(2), "O2 dens");
			dens.add(eomdata.ionsph.p.at(3), "O dens");
		}
		// End spit out csv files
		//
		//

		printIonosphere(eomdata.ionsph);

		// 1. Adjust Maxwellian by Cos, SQRT(B ratio) Factors
		double Aratio_ion_sat{ std::sqrt(eomdata.B_sat / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph to satellite
		double Aratio_mag_sat{ std::sqrt(eomdata.B_sat / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to satellite
		double Aratio_ion_ion{ std::sqrt(eomdata.B_ion / eomdata.B_ion) }; //gyroradius cross-sectional area ratio, ionsph (up, reflect, down) to ionsph
		double Aratio_mag_ion{ std::sqrt(eomdata.B_ion / eomdata.B_mag) }; //gyroradius cross-sectional area ratio, magsph to ionsph

		cout << "ion to sat: " << Aratio_ion_sat << "\n";
		cout << "mag to sat: " << Aratio_mag_sat << "\n";
		cout << "ion to ion: " << Aratio_ion_ion << "\n";
		cout << "mag to ion: " << Aratio_mag_ion << "\n";

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
				if (eomdata.dnward.pitch.at(part) >= 90.0)
					throw logic_error("ionosphere::steadyFlux : downward pitch of particle is >= 90 deg");

				maxwellian_sat.at(part) *= 1.0 / cos(eomdata.dnward.pitch.at(part) * RADS_PER_DEG) * Aratio_mag_sat;
				maxwellian_ion.at(part) *= Aratio_mag_ion;
			}
			else
				throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source - it's somewhere between the ionosphere and magnetosphere");

			if (eomdata.bottom.pitch.at(part) >= 90.0)
				throw logic_error("ionosphere::steadyFlux : bottom pitch of particle is >= 90 deg");
		}

		TESTVEC_NOTNEGWHOLEVEC(maxwellian_sat, "steadyFlux::maxwellian_sat");


		// 2. Calculate dEfluxes
		dEflux_v2D distfluxdnward{ dEFlux::satellite(eomdata.dnward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D distfluxupward{ dEFlux::satellite(eomdata.upward, eomdata.satbins, maxwellian_sat) };
		dEflux_v2D bkscfluxupward{ dEFlux::backscatr(eomdata, maxwellian_ion) };
		
		printVec2D(distfluxdnward, "Dnward Flux at Satellite");
		printVec2D(distfluxupward, "Upward Flux at Satellite");
		printVec2D(bkscfluxupward, "Flux Due to Backscatter at Satellite");


		// 3. Sum dEfluxes
		for (size_t ang = 0; ang < distfluxupward.size(); ang++)
			for (size_t eny = 0; eny < distfluxupward.at(ang).size(); eny++)
				distfluxupward.at(ang).at(eny) += distfluxdnward.at(ang).at(eny) + bkscfluxupward.at(ang).at(eny);

		//
		//
		// TEST BACKSCATTER ONLY //
		return bkscfluxupward;
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

			TESTVEC_ISZEROFRSTHALF(dNatIonsph, "backscatr::dNatIonsph");

			// 1.2. Calculate BS dNflux from dNflux Incident to Ionosphere
			dNflux_v2D BSatIonsph{ multiLevelBS::scatterMain(eomdata, dNatIonsph) }; //new multi-level hotness
			// output: 2D vector of backscatter dNFlux by dist bins at ionosphere (upgoing)

			TESTVEC_ISZEROLASTHALF(BSatIonsph, "backscatr::BSatIonsph");

			// 1.3. Translate BS dNflux at Ionosphere (dist binned) to dNflux at Satellite (sat binned)
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

		/*constexpr double JOHND_PRIM_LOGM{ 3.5 }; //was 3.5 (from John) //0.8, 0.4 matches high E bump
		constexpr double JOHND_PRIM_LOGB{ -5.8 }; //was -5.8 (from John) //-5.5 matches high E bump
		constexpr double JOHND_SECD_LOGM_LT10{ -1.0 }; //was -1.0 (from John)
		constexpr double JOHND_SECD_LOGB_LT10{ -2.5 }; //was -2.5 (from John)
		constexpr double JOHND_SECD_LOGM_GT10{ -1.8 }; //was -2.25 (from John, but he changed to -1.8) //was at -4.0
		constexpr double JOHND_SECD_LOGB_GT10{ -1.7 }; //was -1.7 (from John) //was at -1.0*/

		constexpr double JOHND_PRIM_LOGM{ 1.132 }; //was 3.5 (from John) //0.8, 0.4 matches high E bump
		constexpr double JOHND_PRIM_LOGB{ -4.90 }; //was -5.8 (from John) //-5.5 matches high E bump
		constexpr double JOHND_SECD_LOGM_LT10{ -1.0 }; //was -1.0 (from John)
		constexpr double JOHND_SECD_LOGB_LT10{ -2.5 }; //was -2.5 (from John)
		constexpr double JOHND_SECD_LOGM_GT10{ -2.245 }; //was -2.25 (from John, but he changed to -1.8) //was at -4.0
		constexpr double JOHND_SECD_LOGB_GT10{  0.153 }; //was -1.7 (from John) //was at -1.0

		DLLEXP dNflux johnd_flux(eV E_eval, eV E_incident, dNflux dN_incident)
		{
			if (E_eval > E_incident * (1 + FLT_EPSILON)) return 0.0;
				//throw logic_error("johnd_flux: E_eval is higher than E_incident.  Not physical.  Eval, Incident: " + to_string(E_eval) + " , " + to_string(E_incident));

			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			return dN_incident * pow(10.0, secd_logm * log10(E_eval) + secd_logb) + //secondary BS
				   dN_incident * (10000.0 / E_incident) * pow(10.0, JOHND_PRIM_LOGM * log10(E_eval / E_incident) + JOHND_PRIM_LOGB); //primary BS
		}

		/*DLLEXP dEflux integralJohnd_flux(double lower, double upper, double E_incident)
		{
			throw std::exception("integralJohnd_flux used");//make sure this isnt used for now
			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			double integral_sec{ (pow(upper, secd_logm + 1.0) - pow(lower, secd_logm + 1.0)) * pow(10.0, secd_logb) / (secd_logm + 1.0) };
			double integral_prm{ (pow(upper, JOHND_PRIM_LOGM + 1.0) - pow(lower, JOHND_PRIM_LOGM + 1.0)) * pow(10.0, JOHND_PRIM_LOGB + 4.0) / ((JOHND_PRIM_LOGM + 1.0) * pow(E_incident, JOHND_PRIM_LOGM + 1.0)) };
			return integral_sec + integral_prm;
		}*/

		DLLEXP dNflux_v2D downwardToBackscatter(const Bins& dist, const dNflux_v2D& dNpointofScatter)
		{ //converts downward dNflux at ionosphere (dist binned) to bs (upward) dNflux (also dist binned)
			// 1. Sum dNflux over PA Bins (dist bins), Per E Bin and Average
			dNflux_v1D dNsumPerE(dist.E.size());             //Sum of escaped particles at each energy, units of dNflux

			for (size_t egy = 0; egy < dist.E.size(); egy++) //iterate over energies
			{
				for (size_t ang = 0; ang < dist.PA.size(); ang++)
				{
					if (dist.PA.at(ang) > 90.0) continue;

					degrees dangle{ (ang == 0) ?
						(abs(dist.PA.at(ang + 1) - dist.PA.at(ang))) : //if we are at first index, ang - 1 doesn't exist
						(abs(dist.PA.at(ang) - dist.PA.at(ang - 1))) };//if we are at last index, ang + 1 doesn't exist

					dNsumPerE.at(egy) += dNpointofScatter.at(ang).at(egy) * 4 * PI * 27.85 *
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
			dNflux_v2D BS(dist.PA.size());

			for (size_t ang = 0; ang < dist.PA.size(); ang++)
			{
				BS.at(ang) = dNflux_v1D(dist.E.size());

				if (dist.PA.at(ang) <= 90.0)
					continue; //empty vector of the right size
				else
				{
					for (size_t eny = 0; eny < dist.E.size(); eny++)
						BS.at(ang).at(eny) = dNbsPerE.at(eny) * -cos(dist.PA.at(ang) * RADS_PER_DEG) / (4 * PI);
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
			{ //outputs pitch angle and energy bins to a CSV
				utils::fileIO::CSV bins{ "debug/angle, pitch bins - satellite/bins.csv" };

				bins.add(eom.satbins.E, "E");
				bins.add(eom.satbins.PA, "PA");
			}
			{
				utils::fileIO::CSV bins{ "debug/angle, pitch bins - simulation/bins.csv" };

				bins.add(eom.distbins.E, "E");
				bins.add(eom.distbins.PA, "PA");
			}
			// End spit out csv files
			//
			//

			double_v2D pctScattered(eom.distbins.PA.size(), dNflux_v1D(eom.distbins.E.size())); //% scattered per bin

			ParticleData particles;

			auto newPA = [](const degrees PA_init, const tesla B_init, const tesla B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle has reflects before B_final

				degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations

				return ret;
			};

			// >> level calculate
			for (size_t level = 0; level < eom.ionsph.s.size() - 1; level++)
			{//for now, depends on adding one extra level (to see if particles reflect somewhere within the last layer)
				printLayer(eom.ionsph, (unsigned int)level);

				dNflux_v2D bs_level{ bsAtLevel(eom, dNionsphTop, pctScattered, (unsigned int)level) };

				TESTVEC_ISZEROLASTHALF(bs_level, "scatterMain::bs_level");

				double Aratio_bslevel_ion{ std::sqrt(eom.B_ion / eom.ionsph.B.at(level)) };

				// >> >> adjust PA
				for (size_t ang = 0; ang < eom.distbins.PA.size(); ang++)
				{
					if (eom.distbins.PA.at(ang) <= 90.0) continue;

					degrees paTopIon{ newPA(eom.distbins.PA.at(ang), eom.ionsph.B.at(level), eom.ionsph.B.at(0)) };
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

		DLLEXP dNflux_v2D bsAtLevel(const EOMSimData& eom, const dNflux_v2D& dNionsphTop, double_v2D& pctScatteredAbove, unsigned int level)
		{
			//
			//
			// Spit out csv files (for diagnostics)
			{ //produces % scatter probability at each level, per PA bin, per E bin
				vector<string> labels(eom.satbins.PA.size() + 1);
				for (int label = 0; label < eom.satbins.PA.size(); label++)
					labels.at(label) = to_string(eom.satbins.PA.at(label));

				string levelattrs{ "" };
				for (size_t specs = 0; specs < eom.ionsph.names.size(); specs++) levelattrs += eom.ionsph.names.at(specs) + " ";
				levelattrs += "Z p h";
				labels.back() = levelattrs;

				vector<vector<double>> prob(eom.satbins.PA.size() + 1, vector<double>(eom.satbins.E.size()));
				prob.back().resize(eom.ionsph.p.size() * 3);
				for (size_t PAs = 0; PAs < eom.satbins.PA.size(); PAs++)
				{
					for (size_t Es = 0; Es < eom.satbins.E.size(); Es++)
					{
						for (size_t species = 0; species < eom.ionsph.p.size(); species++)
						{
							if (PAs == 0 && Es == 0)
							{
								prob.back().at(3 * species) = eom.ionsph.Z.at(species);
								prob.back().at(3 * species + 1) = eom.ionsph.p.at(species).at(level);
								prob.back().at(3 * species + 2) = eom.ionsph.h.at(level);
							}
							prob.at(PAs).at(Es) += scatterPct(0.0, eom.ionsph.Z.at(species),
								eom.ionsph.p.at(species).at(level),	eom.ionsph.h.at(species),
								eom.satbins.E.at(Es), eom.satbins.PA.at(PAs));
						}
					}
				}

				utils::fileIO::CSV scatprobs{ "debug/level scatter probs/scatter_" + to_string((int)eom.ionsph.s.at(level)) + ".csv" };
				scatprobs.add(eom.satbins.E, "Es + PAs");
				scatprobs.add(prob, labels);
			}
			// End spit out csv files
			//
			//

			ParticleData dngoing; //TO DO: reserve size of vectors to increase the speed
			ParticleData upgoing;

			double Aratio_ion_bslevel{ std::sqrt(eom.ionsph.B.at(level) / eom.B_ion) };

			//lambda that generates new PA from an initial PA based on B field strength change
			auto newPA = [](const degrees PA_init, const tesla B_init, const tesla B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

				degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

				return ret;
			};

			//lambda that generates scatter/reflect percentage and updates total scatter % of defined particle
			auto sctReflPct = [&](double& pctScatteredAbove, const eV E, const degrees pitch)
			{
				//upgoing - whatever hasn't scattered so far, reflects
				//no scattering happens in this case
				if (pitch > 90.0)
					return (1.0 - pctScatteredAbove); //so return 1.0 - "sum % scattered in layers above"
				
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
					for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
					{
						if (dNionsphTop.at(ang).at(eny) != 0.0 && pctScatteredAbove.at(ang).at(eny) < 1.0)
						{ //these particles are upgoing
							if (180.0 - pa_level <= 90.0)
								throw logic_error("ionosphere::multiLevelBS::bsAtLevel: sctReflPct will not return (1 - % scattered), logic error somewhere");

							degrees pa_sat{ newPA(180.0 - pa_level, eom.ionsph.B.at(level), eom.B_sat) }; //pitch at satellite - eventually may need to run through sourceToSatellite
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
			{ //produces dN refl and scat at each level, per PA bin, per E bin
				vector<string> scatlabels;
				vector<string> refllabels;
				vector<vector<double>> refl;
				vector<vector<double>> scat;
				vector<vector<double>> backscat;

				for (size_t ang = 0; ang < eom.distbins.PA.size(); ang += 5000)
				{
					if (eom.distbins.PA.at(ang) < 90.0)
					{
						scatlabels.push_back(to_string(eom.distbins.PA.at(ang)));

						vector<double> tmp(eom.distbins.E.size());
						vector<double> tmp2(eom.distbins.E.size());
						for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
						{
							tmp.at(eny) = dnBinned.at(ang).at(eny);
							tmp2.at(eny) = ret.at(ang).at(eny);
						}
						scat.push_back(tmp);
						backscat.push_back(tmp2);
					}
					else
					{
						refllabels.push_back(to_string(eom.distbins.PA.at(ang)));

						vector<double> tmp(eom.distbins.E.size());
						for (size_t eny = 0; eny < eom.distbins.E.size(); eny++)
						{
							tmp.at(eny) = upBinned.at(ang).at(eny);
						}
						refl.push_back(tmp);
					}
				}

				utils::fileIO::CSV bksccsv{ "debug/level dN upgoing/bksc_" + to_string((int)eom.ionsph.s.at(level)) + ".csv" };
				utils::fileIO::CSV scatcsv{ "debug/level dN upgoing/scat_" + to_string((int)eom.ionsph.s.at(level)) + ".csv" };
				utils::fileIO::CSV reflcsv{ "debug/level dN upgoing/refl_" + to_string((int)eom.ionsph.s.at(level)) + ".csv" };
				bksccsv.add(eom.distbins.E, "Es + PAs");
				bksccsv.add(backscat, refllabels);
				scatcsv.add(eom.distbins.E, "Es + PAs");
				scatcsv.add(scat, scatlabels);
				reflcsv.add(eom.distbins.E, "Es + PAs");
				reflcsv.add(refl, refllabels);
			}
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
