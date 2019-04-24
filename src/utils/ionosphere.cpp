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
	DLLEXP dEflux_v2D steadyFlux(const EOMSimData& eomdata2)
	{
		//generate "Ideal(tm)" distribution, compare with sim
		/*auto newPA = [](degrees PA_init, tesla B_init, tesla B_final)
		{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
			double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

			if (one < 0.0) return -1.0; //if this is the case, particle reflects before B_final

			degrees ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

			if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations, so invert ret to be above 90

			return ret;
		};

		ParticleData btm(3456000);
		ParticleData upw(3456000);
		ParticleData dnw(3456000);

		double btmPAerr{ 0.0 };
		double upwPAerr{ 0.0 };
		double dnwPAerr{ 0.0 };
		double btmEerr{ 0.0 };
		double upwEerr{ 0.0 };
		double dnwEerr{ 0.0 };
		int extrabtm{ 0 };
		int extraupw{ 0 };
		int extradnw{ 0 };

		for (int part = 0; part < eomdata2.initial.energy.size(); part++)
		{
			if (eomdata2.initial.s_pos.at(part) < eomdata2.s_ion * 1.001)
			{
				double satPA{ newPA(eomdata2.initial.pitch.at(part), eomdata2.B_ion, eomdata2.B_sat) };

				if (satPA > 0)
				{
					upw.pitch.at(part) = satPA;
					upw.energy.at(part) = eomdata2.initial.energy.at(part);
				}
			}
			else if (eomdata2.initial.s_pos.at(part) > eomdata2.s_mag * 0.999)
			{
				double satPA{ newPA(eomdata2.initial.pitch.at(part), eomdata2.B_mag, eomdata2.B_sat) };
				double btmPA{ newPA(eomdata2.initial.pitch.at(part), eomdata2.B_mag, eomdata2.B_ion) };

				if (satPA > 0)
				{
					dnw.pitch.at(part) = satPA;
					dnw.energy.at(part) = eomdata2.initial.energy.at(part);
					if (btmPA < 0)
					{
						upw.pitch.at(part) = 180.0 - satPA;
						upw.energy.at(part) = eomdata2.initial.energy.at(part);
					}
				}
				if (btmPA > 0)
				{
					btm.pitch.at(part) = btmPA;
					btm.energy.at(part) = eomdata2.initial.energy.at(part);
				}
			}
			else
			{
				throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source");
			}

			if ((btm.pitch.at(part) != 0.0) && (eomdata2.bottom.pitch.at(part) != 0.0) && (abs(btm.pitch.at(part) - eomdata2.bottom.pitch.at(part)) > btmPAerr))
				btmPAerr = abs(btm.pitch.at(part) - eomdata2.bottom.pitch.at(part));
			if ((upw.pitch.at(part) != 0.0) && (abs(upw.pitch.at(part) - eomdata2.upward.pitch.at(part)) > upwPAerr))
				upwPAerr = abs(upw.pitch.at(part) - eomdata2.upward.pitch.at(part));
			if ((dnw.pitch.at(part) != 0.0) && (eomdata2.dnward.pitch.at(part) != 0.0) && (abs(dnw.pitch.at(part) - eomdata2.dnward.pitch.at(part)) > dnwPAerr))
				dnwPAerr = abs(dnw.pitch.at(part) - eomdata2.dnward.pitch.at(part));

			if ((btm.pitch.at(part) != 0.0) && (eomdata2.bottom.pitch.at(part) == 0.0))
				extrabtm++;
			if ((upw.pitch.at(part) != 0.0) && (eomdata2.upward.pitch.at(part) == 0.0))
				extraupw++;
			if ((dnw.pitch.at(part) != 0.0) && (eomdata2.dnward.pitch.at(part) == 0.0))
				extradnw++;
		}

		std::cout << btmPAerr << "\n" << upwPAerr << "\n" << dnwPAerr << "\n\n";
		std::cout << extrabtm << "\n" << extraupw << "\n" << extradnw << "\n";

		//exit(1);

		//EOMSimData eomdata{ eomdata2 };
		//eomdata.bottom = btm;
		//eomdata.upward = upw;
		//eomdata.dnward = dnw;*/

		/*try //compare "ideal" array with EOM sim values, use with above "ideal" creation code
		{
			std::string simpath{ "..\\_dataout\\190202_11.37.17.620km.dtDivBy10\\" };
			//std::string simpath{ "..\\_dataout\\190112_13.18.14.QSPS800eV.500s\\"};

			ParticleData topData(3456000, false); //need to load data at top of sim
			ParticleData finalData(3456000, false);

			{
				Satellite top("topElec", { "vpara", "vperp", "s", "time", "index" }, eomdata2.s_mag, false, 3456000, nullptr);
				top.loadDataFromDisk(simpath + "bins\\satellites\\");

				topData.vpara = top.data().at(0).at(0);
				topData.vperp = top.data().at(0).at(1);
				utils::numerical::v2DtoEPitch(topData.vpara, topData.vperp, eomdata.mass, topData.energy, topData.pitch);

				utils::fileIO::readDblBin(finalData.vpara, simpath + "bins\\particles_final\\elec_vpara.bin");
				utils::fileIO::readDblBin(finalData.vperp, simpath + "bins\\particles_final\\elec_vperp.bin");
				utils::fileIO::readDblBin(finalData.s_pos, simpath + "bins\\particles_final\\elec_s.bin");
				utils::numerical::v2DtoEPitch(finalData.vpara, finalData.vperp, eomdata.mass, finalData.energy, finalData.pitch);
			}

			int btmIdealSize{ 0 };
			int topIdealSize{ 0 };
			int btmSimSize{ 0 };
			int topSimSize{ 0 };
			int simMissDn{ 0 };
			int simMissUp{ 0 };

			for (int part = 0; part < eomdata.bottom.energy.size(); part++)
			{
				if (eomdata2.bottom.energy.at(part) > 0.0)
					btmSimSize++;
				if (topData.energy.at(part) > 0.0)
					topSimSize++;

				if (eomdata.bottom.energy.at(part) > 0.0)
					btmIdealSize++;

				if (eomdata.dnward.energy.at(part) != 0.0 && eomdata2.dnward.energy.at(part) == 0.0)
					simMissDn++;
				if (eomdata.upward.energy.at(part) != 0.0 && eomdata2.upward.energy.at(part) == 0.0)
					simMissUp++;

				//if (eomdata.bottom.energy.at(part) == 0.0 && topEnergy.at(part) == 0.0)
				//{
					//if (part == 277949) continue;
					//std::cout << "Final: " << part << " ind, " << finalEnergy.at(part) << " eV, " << finalPitchA.at(part) << " deg, " << finals.at(part) << " m\n";
					//std::cout << "Init:  " << part << " ind, " << eomdata.initial.energy.at(part) << " eV, " << eomdata.initial.pitch.at(part) << " deg, " << eomdata.initial.s_pos.at(part) << " m\n";
					//std::cout << "Btm:   " << part << " ind, " << eomdata.bottom.energy.at(part)  << " eV, " << eomdata.bottom.pitch.at(part)  << " deg, " << eomdata.bottom.s_pos.at(part)  << " m\n";
					//std::cout << "       " << eomdata.bottom.vpara.at(part) << " vpara m/s, " << eomdata.bottom.vperp.at(part) << " vperp m/s\n";
					//std::cout << "Upw:   " << part << " ind, " << eomdata.upward.energy.at(part) << " eV, " << eomdata.upward.pitch.at(part) << " deg, " << eomdata.upward.s_pos.at(part) << " m\n";
					//std::cout << "       " << eomdata.upward.vpara.at(part) << " vpara m/s, " << eomdata.upward.vperp.at(part) << " vperp m/s\n";
					//std::cout << "Dnw:   " << part << " ind, " << eomdata.dnward.energy.at(part) << " eV, " << eomdata.dnward.pitch.at(part) << " deg, " << eomdata.dnward.s_pos.at(part) << " m\n";
					//std::cout << "       " << eomdata.dnward.vpara.at(part) << " vpara m/s, " << eomdata.dnward.vperp.at(part) << " vperp m/s\n";
					//exit(1);
				//}
			}

			std::cout << "\n\n\n" << btmSimSize << "\n" << topSimSize << "\n" << "\n";
			std::cout << btmIdealSize << "\n" << "\n";
			std::cout << simMissDn << "\n" << simMissUp << "\n";
		}
		catch (std::exception& e)
		{
			std::cout << e.what() << "\n";
		}

		exit(1);*/

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
1.034400169,
1.042600148,
1.055748627,
1.07152601,
1.078133553,
1.058734323,
1.039323547,
1.012305997,
0.994205528,
0.985836018,
0.910594196,
0.838235396,
0.759978758,
0.688727757,
0.582535504,
0.484989966,
0.393631204,
0.330308124,
0.27732655,
0.225898359,
0.178965883,
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

		std::vector<double> realIonMaxwellian = {
		628.8979407,
565.8638197,
508.6540186,
456.7303733,
409.6044458,
366.8329293,
328.0134787,
292.780925,
260.8038409,
231.7814228,
205.4406609,
181.5337716,
159.4832008,
138.7981913,
120.0244659,
102.9854228,
89.92642941,
78.43829219,
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

		EOMSimData eomdata{ eomdata2 };
		for (unsigned int part = 0; part < eomdata.initial.pitch.size(); part++)
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
				if (eomdata.dnward.pitch.at(part) >= 90.0)
				{
					std::cout << "Warning! Downgoing particle pitch >= 90.0. Taking 180-PA. ind, s, pitch, energy: "
						<< part << ", " << eomdata.dnward.s_pos.at(part) << ", " << eomdata.dnward.pitch.at(part)
						<< ", " << eomdata.dnward.energy.at(part) << "\n";
					eomdata.dnward.pitch.at(part) = 180.0 - eomdata.dnward.pitch.at(part);
				}

				maxwellian_sat.at(part) *= 1.0 / cos(eomdata.dnward.pitch.at(part) * RADS_PER_DEG) * Aratio_mag_sat;
				maxwellian_ion.at(part) *= Aratio_mag_ion;
			}
			else
				throw logic_error("ionosphere::steadyFlux : particle is not ionospheric or magnetospheric source");

			if (eomdata.bottom.pitch.at(part) >= 90.0)
			{
				std::cout << "Warning! Bottom particle pitch >= 90.0. Taking 180-PA. ind, s, pitch, energy: "
					<< part << ", " << eomdata.bottom.s_pos.at(part) << ", " << eomdata.bottom.pitch.at(part)
					<< ", " << eomdata.bottom.energy.at(part) << "\n";
				eomdata.bottom.pitch.at(part) = 180.0 - eomdata.bottom.pitch.at(part);
			}
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
		constexpr double JOHND_SECD_LOGM_GT10{ -1.8 }; //was -2.25 (from John, but he changed to -1.8) //was at -4.0
		constexpr double JOHND_SECD_LOGB_GT10{ -1.7 }; //was -1.7 (from John) //was at -1.0

		DLLEXP dNflux johnd_flux(double E_eval, double E_incident, dEflux dE_incident)
		{
			if (E_eval > E_incident * (1 + FLT_EPSILON))// return 0.0;
				throw logic_error("johnd_flux: E_eval is higher than E_incident.  Not physical.  Eval, Incident: " + to_string(E_eval) + " , " + to_string(E_incident));

			double secd_logm{ (E_incident <= 10.0) ? JOHND_SECD_LOGM_LT10 : JOHND_SECD_LOGM_GT10 };
			double secd_logb{ (E_incident <= 10.0) ? JOHND_SECD_LOGB_LT10 : JOHND_SECD_LOGB_GT10 };

			return dE_incident * pow(10.0, secd_logm * log10(E_eval) + secd_logb) + //secondary BS
				  (dE_incident / E_incident) * 10000.0 * pow(10.0, JOHND_PRIM_LOGM * log10(E_eval / E_incident) + JOHND_PRIM_LOGB); //primary BS
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
				for (unsigned int ang = 0; ang < dist.PA.size(); ang++)
				{
					if (dist.PA.at(ang) > 90.0) continue;

					degrees dangle{ (ang == 0) ?
						(abs(dist.PA.at(ang + 1) - dist.PA.at(ang))) : //if we are at first index, ang - 1 doesn't exist
						(abs(dist.PA.at(ang) - dist.PA.at(ang - 1))) };//if we are at last index, ang + 1 doesn't exist

					dNsumPerE.at(egy) += downwardAtIonsph.at(ang).at(egy) * 4 * PI * 27.85 *
						abs(sin(dist.PA.at(ang) * RADS_PER_DEG) * sin(dangle / 2.0 * RADS_PER_DEG) *
							cos(dist.PA.at(ang) * RADS_PER_DEG) * cos(dangle / 2.0 * RADS_PER_DEG));
				}
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

			for (unsigned int ang = 0; ang < dist.PA.size(); ang++)
			{
				BS.at(ang) = dNflux_v1D(dist.E.size());

				if (dist.PA.at(ang) <= 90.0)
					continue; //empty vector of the right size
				else
				{
					for (unsigned int eny = 0; eny < dist.E.size(); eny++)
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

		DLLEXP double scatterPct(double sumCollideAbove, double Z, double p, double h, eV E, degrees PA)
		{
			if (sumCollideAbove >= 1.0) throw logic_error("ionosphere::multiLevelBS::scatterPct: sumCollideAbove is greater than / equal to 1.0.  "
				+ string("100% has scattered already.  Conditions in bsAtLevel should have prevented this from happening."));

			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	} //end namespace multiLevelBS

} //end namespace ionosphere
