#include "utils/postprocess.h"
#include "utils/numerical.h"
#include <algorithm>
#include <cmath>

#include "utils/fileIO.h"

#include <iomanip>

using std::pow;
using std::string;
using utils::numerical::generateSpacedValues;

constexpr double EVANS_PRIM_LOGM{ 1.5 };  //obtained by log linefitting Evans, 1974 - these seem closest
constexpr double EVANS_PRIM_LOGB{ -4.0 };
constexpr double EVANS_SECD_LOGM{ -2.1 };
constexpr double EVANS_SECD_LOGB{ 0.3 };

#define TESTVEC_ISZEROFRSTHALF(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name, 0, (unsigned int)vec.size() / 2);
#define TESTVEC_ISZEROLASTHALF(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name, (unsigned int)vec.size() / 2);
#define TESTVEC_ISZEROWHOLEVEC(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name);
#define TESTVEC_NOTNEGWHOLEVEC(vec, name) vecTest(vec, [](double cnt) { return (cnt < 0.0); }, true, name);

inline dblVec serialize2DVec(const dblVec2D& in)
{
	dblVec out;

	for (auto& v : in)
		out.insert(std::end(out), std::begin(v), std::end(v));

	return out;
}

inline bool vecTest(const dblVec2D& vec, std::function<bool(double)> test, bool throwOnTrue=false,
	string label = "", unsigned int outStart=0, unsigned int outStop=0, unsigned int inStart=0, unsigned int inStop=0)
{
	if (outStop == 0) outStop = (unsigned int)vec.size();

	for (unsigned int out = outStart; out < outStop; out++)
	{
		if (inStop == 0) inStop = (unsigned int)vec.at(out).size();
		for (unsigned int in = 0; in < inStop; in++)
			if (test(vec.at(out).at(in)))
			{
				if (throwOnTrue)
					throw std::logic_error("dblVec2DTest: " + label + " condition met.  Throwing - out, in, data: " +
					std::to_string(out) + ", " + std::to_string(in) + ", " + std::to_string(vec.at(out).at(in)));
				
				return true;
			}
	}

	return false;
}

inline bool vecTest(const dblVec& vec, std::function<bool(double)> test, bool throwOnTrue = false,
	string label = "", unsigned int start = 0, unsigned int stop = 0)
{
	dblVec2D tmp;
	tmp.push_back(vec);

	return vecTest(tmp, test, throwOnTrue, label, 0, 0, start, stop);
}

void printIonosphere(const postprocess::Ionosphere& ionsph)
{
	std::cout << "==================== Backscatter Simulation ====================" << "\n";
	std::cout << "Min, Max s (m): " << *(ionsph.s.end() - 1) << ", " << ionsph.s.front() << "\n";
	std::cout << "Min, Max B (T): " << *(ionsph.B.end() - 1) << ", " << ionsph.B.front() << "\n";
	std::cout << "Num of Layers : " << ionsph.s.size() - 1 << "\n";
	std::cout << "Atomic Species: " << "Fix later.\n";
	std::cout << "================================================================" << "\n";
}

void printLayer(const postprocess::Ionosphere& ionsph, unsigned int layer)
{
	std::cout << "Layer: " << layer << " / " << ionsph.s.size() - 2 << ", s: " << ionsph.s.at(layer) << ", B: " << ionsph.B.at(layer) << "\n";
}

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

		printIonosphere(ppdata.ionsph);

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
				//weights_atIon.at(iii) *= 1.0 / cos(ppdata.bottom.pitch.at(iii) * RADS_PER_DEG) * Aratio_mag_bs * bsScale;
				weights_atIon.at(iii) *= Aratio_mag_bs * bsScale;
			}
			else
				throw std::logic_error("postprocess::steadyFlux : particle is not ionospheric or magnetospheric source");
		}


		////
		////
		utils::fileIO::writeDblBin(weights_atSat, "dat\\00weights_atSat.bin", (unsigned int)weights_atSat.size());
		utils::fileIO::writeDblBin(weights_atIon, "dat\\01weights_atIon.bin", (unsigned int)weights_atIon.size());
		////
		////


		TESTVEC_NOTNEGWHOLEVEC(weights_atSat, "steadyFlux::weights_atSat");
		TESTVEC_NOTNEGWHOLEVEC(weights_atIon, "steadyFlux::weights_atIon");

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

			utils::fileIO::writeDblBin(serialize2DVec(escapeCountBinned), "dat\\02escapeCntBin.bin", (unsigned int)(escapeCountBinned.size()*escapeCountBinned.front().size()));

			TESTVEC_ISZEROFRSTHALF(escapeCountBinned, "bksdEFlux::escapeCountBinned");

			// 1.2. Calculate BS dNflux from dNflux Incident to Ionosphere
			//dblVec2D dNflux_BS{ backscat::dNflux_bs_ion(ppdata.distbins, escapeCountBinned) }; //original, single level
			dblVec2D dNflux_BS{ multLevelBS::scatterMain(ppdata.ionsph, ppdata.distbins, escapeCountBinned, ppdata.B_sat) }; //new multi-level hotness
			// output: 2D vector of backscatter dNFlux by dist bins at ionosphere (upgoing)

			utils::fileIO::writeDblBin(serialize2DVec(dNflux_BS), "dat\\03dNflux_BS.bin", (unsigned int)(dNflux_BS.size()*dNflux_BS.front().size()));

			TESTVEC_ISZEROLASTHALF(dNflux_BS, "bksdEFlux::dNflux_BS");

			// 1.3. Translate BS dNflux at Ionosphere (dist binned) to dNflux at Satellite (sat binned)
			dblVec2D ret{ backscat::bsSrcToSat(ppdata.distbins, ppdata.satbins, dNflux_BS, ppdata.initial, ppdata.upward) };
			// output: 2D vector of bs dNFlux at satellite per PA, E (sat binned) - should only be upward (90-180)
			// Section 1 End
			
			utils::fileIO::writeDblBin(serialize2DVec(ret), "dat\\04dNflux_BS_atSat.bin", (unsigned int)(ret.size()*ret.front().size()));

			TESTVEC_ISZEROFRSTHALF(ret, "bksdEFlux::backscatter");

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
				if (!Eascending) engbin = (int)bins.E.size()  - 1 - engbin; //reverses the bin index if E is not ascending
				if (!Aascending) angbin = (int)bins.PA.size() - 1 - angbin; //ditto

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
			
			for (unsigned int iii = (unsigned int)binAngles.size() / 2; iii < (unsigned int)binAngles.size(); iii++)
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
			utils::fileIO::writeDblBin(dist.PA, "dat\\3.1distPA.bin", (unsigned int)dist.PA.size());
			utils::fileIO::writeDblBin(dist.E, "dat\\3.2distE.bin", (unsigned int)dist.E.size());
			utils::fileIO::writeDblBin(serialize2DVec(escapeCountBinned), "dat\\3.3binnedCountsIn.bin", (unsigned int)(escapeCountBinned.size()*escapeCountBinned.front().size()));

			// 1. Sum dNflux over PA Bins (dist bins), Per E Bin and Average
			dblVec escapeCountPerE(dist.E.size());                 //Sum of escaped particles at each energy, units of dNflux
			for (unsigned int egy = 0; egy < dist.E.size(); egy++) //iterate over energies
			{
				for (auto& angvec : escapeCountBinned)                //sum over angles, add to sum vector at energyBin
					escapeCountPerE.at(egy) += angvec.at(egy);
				escapeCountPerE.at(egy) /= (double)dist.PA.size() / 2; //isotropically space across pitch angle bins - divide by # ion PA bins, this value later put in each ionospheric angle bin
			}
			// output: 1D vector of total number of escaped particles (dNFlux) per energy, reduced by # of ionsph pitch bins

			utils::fileIO::writeDblBin(escapeCountPerE, "dat\\3aEscapeCntPerE.bin", (unsigned int)escapeCountPerE.size());
			
			// 2. Calculate upward dNflux (backscatter) per E bin
			double logEBinMin{ log10(dist.E.at(0)) };         //depends on an array where E is minimum at index 0, max at last index
			double dlogE{ log10(dist.E.at(1)) - logEBinMin }; //depends on a logarithmically spaced E, won't work otherwise

			dblVec upwardCountPerE{ escapeCountPerE };
			for (unsigned int ebin = 0; ebin < escapeCountPerE.size(); ebin++) //why is this the case!!!
				upwardCountPerE.at(ebin) *= dist.E.at(ebin);  //convert to dEflux escaping into the layer

			utils::fileIO::writeDblBin(upwardCountPerE, "dat\\3bEscapeEPerE.bin", (unsigned int)upwardCountPerE.size());

			dblVec dNfluxPerE_bs(dist.E.size());
			for (unsigned int dNFluxBin = 0; dNFluxBin < dNfluxPerE_bs.size(); dNFluxBin++)       //bins that contain the number flux of the backscatter in the energy bin of the same index
			{
				double engmin{ pow(10, (dNFluxBin - 0.5) * dlogE + logEBinMin) };                 //assumes evenly log spaced angle bins
				double engmax{ pow(10, (dNFluxBin + 0.5) * dlogE + logEBinMin) };

				for (unsigned int incEbin = dNFluxBin; incEbin < dNfluxPerE_bs.size(); incEbin++) //nFlux bins are contributed to by all incident particles with E higher than the E associated with the nFlux bin
				{//change to mid bin
					double incidentE{ dist.E.at(incEbin) };                                       //incident E is upper limit of bin
					double intF{ integralEvans_flux(engmin, engmax, incidentE) };
					
					//dNfluxPerE_bs.at(dNFluxBin) += intF * escapeCountPerE.at(incEbin) * dist.E.at(incEbin) / (engmax - engmin);
					dNfluxPerE_bs.at(dNFluxBin) += intF * upwardCountPerE.at(incEbin) / (engmax - engmin);
										
					//dNfluxPerE_bs.at(dNFluxBin) += intF * escapeCountPerE.at(incEbin) / (engmax - engmin);
					//dNfluxPerE_bs.at(dNFluxBin) += john_flux(binEnergies.at(dNFluxBin), binEnergies.at(incEbin), primary_logm, primary_logb, secondary_logm, secondary_logb) * binCntxE.at(incEbin);
				}
			}
			// output: 1D vector of the upgoing (backscatter) dNflux per E
			
			utils::fileIO::writeDblBin(dNfluxPerE_bs, "dat\\3cEscapeEPerE.bin", (unsigned int)dNfluxPerE_bs.size());

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
						dNfluxPerEPA_bs.at(ang).at(eny) *= -cos(dist.PA.at(ang) * RADS_PER_DEG);
				}
			}

			// normalize //
			/*for (unsigned int eny = 0; eny < dNfluxPerEPA_bs.front().size(); eny++)
			{
				double nonNormTotal{ 0.0 };
				for (unsigned int ang = 0; ang < dNfluxPerEPA_bs.size(); ang++)
					nonNormTotal += dNfluxPerEPA_bs.at(ang).at(eny);

				for (unsigned int ang = 0; ang < dNfluxPerEPA_bs.size(); ang++)
				{
					if (nonNormTotal != 0.0)
						dNfluxPerEPA_bs.at(ang).at(eny) *= escapeCountPerE.at(eny) / nonNormTotal;
				}
			}*/
			// remove if doesn't work //
			// output: 2D vector of bs dNFlux at ionosphere per pitch, energy bin - should only be upward (90-180)

			return dNfluxPerEPA_bs;
		}

		DLLEXP dblVec2D bsSrcToSat(const Bins& dist, const Bins& sat, const dblVec2D& bsCounts, const ParticleData& initialData, const ParticleData& satUpwardData)
		{
			const double ANGMAXERR{ 0.1 }; //max error in degrees
			const double ENYMAXERR{ 0.1 }; //max error in eV
			auto err = [](double base, double diff) { return std::abs((base - diff) / base); };

			const double logEMinBinMid{ log10(dist.E.at(0)) };
			const double dlogE{ log10(dist.E.at(1)) - logEMinBinMid };
			const double dangle{ dist.PA.at(1) - dist.PA.at(0) };

			#define count s_pos //allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleData particles;

			for (unsigned int ang = 0; ang < dist.PA.size(); ang++)
			{
				for (unsigned int eny = 0; eny < dist.E.size(); eny++) //this works because ionospheric bins are same as distribution
				{
					particles.pitch.push_back(satUpwardData.pitch.at(ang * dist.E.size() + eny));
					particles.energy.push_back(satUpwardData.energy.at(ang * dist.E.size() + eny));
					particles.count.push_back(bsCounts.at(ang).at(eny));
				}
			}

			return binning::binWeighted(particles, sat, particles.count);
		}
	} //end namespace postprocess::backscat

	namespace multLevelBS
	{
		DLLEXP dblVec2D scatterMain(const Ionosphere& ionsph, const Bins& distbins, const dblVec2D& topLevelCounts, double B_sat)
		{
			dblVec2D sumScatter(distbins.PA.size(), dblVec(distbins.E.size()));

			//#define count s_pos //defined above - allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleData particles;

			auto newPA = [](double PA_init, double B_init, double B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle has reflects before B_final

				double ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations

				return ret;
			};

			// >> level calculate
			for (unsigned int level = 0; level < ionsph.s.size() - 1; level++)
			{//for now, depends on adding one extra level (to see if particles reflect somewhere within the last layer)
				printLayer(ionsph, level);
				dblVec2D bs_level{ bsAtLevel(ionsph, distbins, topLevelCounts, sumScatter, B_sat, level) };
				
				TESTVEC_ISZEROLASTHALF(bs_level, "scatterMain::bs_level");

				// >> >> adjust PA
				for (unsigned int ang = 0; ang < distbins.PA.size(); ang++)
				{
					double newpa{ newPA(distbins.PA.at(ang), ionsph.B.at(level), ionsph.B.at(0)) };
					if (newpa < 0.0) throw std::logic_error("postprocess::multLevelBS::scatterMain: PA higher is out of range of newPA function.");

					for (unsigned int eny = 0; eny < distbins.E.size(); eny++)
					{
						if (bs_level.at(ang).at(eny) != 0.0)
						{
							particles.pitch.push_back(newpa);
							particles.energy.push_back(distbins.E.at(eny));
							particles.count.push_back(bs_level.at(ang).at(eny));
						}
					}
				}
			}
			
			return binning::binWeighted(particles, distbins, particles.count);
		}

		DLLEXP dblVec2D bsAtLevel(const Ionosphere& ionsph, const Bins& distbins, const dblVec2D& topLevelCounts, dblVec2D& sumCollideAbove, double B_sat, unsigned int level)
		{
			//#define count s_pos //defined above - allows me to use an existing vector in ParticleData with a name that makes sense
			ParticleData dngoing;
			ParticleData upgoing;
			vector<double*> sumCollideAboveDngoing1D;
			vector<double*> sumCollideAboveUpgoing1D;

			auto newPA = [](double PA_init, double B_init, double B_final)
			{ //relies on constant mu to calculate - if mu is not conserved, this function doesn't give accurate results
				double one{ B_init / B_final * (1.0 + 1.0 / pow(tan(PA_init * RADS_PER_DEG), 2.0)) - 1.0 };

				if (one < 0.0) return -1.0; //if this is the case, particle has reflects before B_final

				double ret{ atan(sqrt(1.0 / one)) / RADS_PER_DEG };

				if (PA_init > 90.0) ret = 180.0 - ret; //anything above 90 is returned as below 90 by these equations

				return ret;
			};

			// >> >> adjust PA, check if reflected
			for (unsigned int ang = 0; ang < distbins.PA.size(); ang++)
			{
				if (distbins.PA.at(ang) > 90.0) continue; //there should be no upgoing in topLevelCounts

				double pa_level{ newPA(distbins.PA.at(ang), ionsph.B.at(0), ionsph.B.at(level)) };

				if (pa_level < 0.0) continue; //particle reflects before this level
				if (pa_level > 90.0) throw std::logic_error("postprocess::multLevelBS::bsAtLevel: downgoing particle ended up with pitch > 90.0 - PA_bin, PA_level, level: "
					+ std::to_string(distbins.PA.at(ang)) + ", " + std::to_string(pa_level) + ", " + std::to_string(level));

				if (newPA(distbins.PA.at(ang), ionsph.B.at(0), ionsph.B.at(level + 1)) < 0.0)
				{ //particle reflects before next level, add all particles of this pitch moving in the opposite direction
					for (unsigned int eny = 0; eny < distbins.E.size(); eny++)
					{
						if (topLevelCounts.at(ang).at(eny) != 0.0)
						{
							upgoing.energy.push_back(distbins.E.at(eny));
							upgoing.pitch.push_back(180.0 - pa_level);
							upgoing.count.push_back(topLevelCounts.at(ang).at(eny));
							sumCollideAboveUpgoing1D.push_back(&sumCollideAbove.at(ang).at(eny));
						}
					}
				}
				else
				{ //particle makes it to the next level - invoke scattering
					for (unsigned int eny = 0; eny < distbins.E.size(); eny++)
					{
						if (topLevelCounts.at(ang).at(eny) != 0.0)
						{
							dngoing.energy.push_back(distbins.E.at(eny));
							dngoing.pitch.push_back(pa_level);
							dngoing.count.push_back(topLevelCounts.at(ang).at(eny));
							sumCollideAboveDngoing1D.push_back(&sumCollideAbove.at(ang).at(eny));
						}
					}
				}
			}
			
			// >> >> adjust by cos factor
			for (unsigned int up = 0; up < upgoing.count.size(); up++) //pitch at satellite
				upgoing.count.at(up) *= 1.0 / abs(cos(newPA(upgoing.pitch.at(up), ionsph.B.at(level), B_sat) * RADS_PER_DEG));

			for (unsigned int dn = 0; dn < dngoing.count.size(); dn++) //pitch at lowest level (level of scattering)
				dngoing.count.at(dn) *= 1.0 / cos(dngoing.pitch.at(dn) * RADS_PER_DEG);

			TESTVEC_NOTNEGWHOLEVEC(upgoing.count, "bsAtLevel:: " + std::to_string(level) + " :: upgoing.count");
			TESTVEC_NOTNEGWHOLEVEC(dngoing.count, "bsAtLevel:: " + std::to_string(level) + " :: dngoing.count");

			// >> >> calculate scattered
			for (unsigned int part = 0; part < dngoing.count.size(); part++)
			{
				double sct{ 0.0 };

				if (*sumCollideAboveDngoing1D.at(part) < 1.0)
				{
					for (size_t species = 0; species < ionsph.p.size(); species++)
					{
						sct += scatterPct(*sumCollideAboveDngoing1D.at(part), ionsph.Z.at(species), ionsph.p.at(species).at(level),
							ionsph.h.at(level), dngoing.energy.at(part), dngoing.pitch.at(part));
					}

					*sumCollideAboveDngoing1D.at(part) += sct;

					if (*sumCollideAboveDngoing1D.at(part) > 1.0)
					{//no more than 100% of particles can have scattered
						sct -= (*sumCollideAboveDngoing1D.at(part) - 1.0);
						*sumCollideAboveDngoing1D.at(part) = 1.0;
					}
				}

				if (sct < 0.0) throw std::logic_error("postprocess::multLevelBS::bsAtLevel: scatter % is < 0.0 - sct%, %collideAbove: "
					+ std::to_string(sct) + ", " + std::to_string(*sumCollideAboveDngoing1D.at(part)));

				//if (sct != 1.0) throw std::logic_error("postprocess::multLevelBS::bsAtLevel: scatter not 100%: sct, PA, E: " +
					//std::to_string(sct) + ", " + std::to_string(dngoing.pitch.at(part)) + ", " + std::to_string(dngoing.energy.at(part)));

				dngoing.count.at(part) *= sct;
			}

			for (unsigned int part = 0; part < upgoing.count.size(); part++)
			{
				upgoing.count.at(part) *= (1.0 - *sumCollideAboveUpgoing1D.at(part)); //whatever has not scattered is reflected
			}

			// >> >> bin particles at level
			dblVec2D dnBinned{ binning::binWeighted(dngoing, distbins, dngoing.count) };
			dblVec2D upBinned{ binning::binWeighted(upgoing, distbins, upgoing.count) };

			TESTVEC_ISZEROFRSTHALF(dnBinned, "bsAtLevel:: " + std::to_string(level) + " ::dnBinned");
			TESTVEC_ISZEROLASTHALF(upBinned, "bsAtLevel:: " + std::to_string(level) + " ::upBinned");

			// >> >> calculate backscatter
			dblVec2D ret{ backscat::dNflux_bs_ion(distbins, dnBinned) };

			TESTVEC_ISZEROLASTHALF(ret, "bsAtLevel:: " + std::to_string(level) + " ::backscatter_level");
			
			for (unsigned int ang = 0; ang < ret.size(); ang++)
				for (unsigned int eny = 0; eny < ret.at(0).size(); eny++)
					ret.at(ang).at(eny) += upBinned.at(ang).at(eny);

			return ret;
		}

		DLLEXP dblVec2D alt_reflect(const Bins& distbins, BField* B, double B_ion, double t)
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

			dblVec2D s_ref(distbins.PA.size(), dblVec(distbins.E.size()));

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

		DLLEXP double scatterPct(double sumCollideAbove, double Z, double p, double h, double E, double PA)
		{
			return (1.0 - sumCollideAbove) * 1.62e-14 * Z * p * h / (pow(E, 2.0) * cos(PA * RADS_PER_DEG));
		}
	} //end namespace multLevelBS

} //end namespace postprocess
