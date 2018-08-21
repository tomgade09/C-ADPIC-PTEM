#include <algorithm>
#include <iterator>


#include "utils/postprocess.h"
#include "utils/numerical.h"
#include "utils/silenceStreamMacros.h"
#include "ErrorHandling/simExceptionMacros.h"

using utils::fileIO::ParticleDistribution;

namespace postprocess
{
	//ParticleData
	ParticleData::ParticleData(const std::vector<double>& v_para, const std::vector<double>& v_perp, double mass) : vpara{ v_para }, vperp{ v_perp } //auto calculate E, Pitch
	{
		utils::numerical::v2DtoEPitch(v_para, v_perp, mass, energy, pitch);
	}

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
	//End ParticleData

	//Maxwellian
	Maxwellian::Maxwellian(double dlogEdist) : dlogE_dist{ dlogEdist } {}

	void Maxwellian::push_back_ion(double E_peak, double dE_magnitude, int partsAtE)
	{
		ionEPeak.push_back(E_peak);
		iondEMag.push_back(dE_magnitude / (double)partsAtE); //scales dE according to how many particles will be in the bin containing E
	}
	void Maxwellian::push_back_mag(double E_peak, double dE_magnitude, int partsAtE)
	{
		magEPeak.push_back(E_peak);
		magdEMag.push_back(dE_magnitude / (double)partsAtE);
	}

	std::vector<double> Maxwellian::counts(const ParticleData& init, const double s_ion, const double s_mag, const double B_ion, const double B_alt, const double B_mag)
	{
		std::vector<double> max(init.energy.size()); //array of maxwellian counts
		std::vector<double> maxtmp; //temporary holder allowing std::transform to be used twice instead of nested loops

		//multiply magnitudes by scaling factors
		//may need to be scaled for a number of reasons
		//ex: differing ranges of pitch angles for ion and mag sources
		std::transform(iondEMag.begin(), iondEMag.end(), iondEMag.begin(), [&](double dE) { return dE * ionModFactor * std::sqrt(B_ion/B_alt); });
		std::transform(magdEMag.begin(), magdEMag.end(), magdEMag.begin(), [&](double dE) { return dE * magModFactor * std::sqrt(B_alt/B_mag); });

		auto count_E = [&](double E_cnt, double E_peak, double dEflux_peak, bool zero)
		{ //calculates count for a given E (E_cnt) for a maxwellian peaked at E_peak with magnitude dEflux_peak
			if (E_peak <= 0.0)
				throw std::logic_error("Maxwellian::counts: invalid maxwellian peak E (le 0.0)");
			if (zero || E_cnt == 0.0)
				return 0.0;
			else
			{
				double binWidth_E { pow(10, log10(E_cnt)  + 0.5 * dlogE_dist) - pow(10, log10(E_cnt)  - 0.5 * dlogE_dist) };
				double binWidth_kT{ pow(10, log10(E_peak) + 0.5 * dlogE_dist) - pow(10, log10(E_peak) - 0.5 * dlogE_dist) };
				return exp(-E_cnt / E_peak) * binWidth_E / E_cnt * (dEflux_peak / (exp(-1.0) * binWidth_kT));
			}
		};

		auto genCounts = [&](const std::vector<double>& EPeak, const std::vector<double>& dEMag, std::function<bool(double)> zero)
		{ //iterates over particles and specified Peak/Magnitude values (ionospheric or magnetospheric) and adds values to "max" (after "maxtmp")
			for (unsigned int entr = 0; entr < EPeak.size(); entr++) //iterate over ionospheric maxwellian specifications
			{
				auto getCount = [&](double E, double s) { return count_E(E, EPeak.at(entr), dEMag.at(entr), zero(s)); };

				std::transform(init.energy.begin(), init.energy.end(), init.s_pos.begin(), std::back_inserter(maxtmp), getCount); //generate maxwellian count if ionospheric particle
				std::transform(maxtmp.begin(), maxtmp.end(), max.begin(), max.begin(), [](double x, double y) { return x + y; }); //add final vector and the tmp vector together
				maxtmp.clear();
			}
		};

		//run the above lambda for ionospheric and magnetospheric sources
		//lambda determines whether or not the particle in question is ion/mag source, returns boolean
		genCounts(ionEPeak, iondEMag, [&s_ion](double s) { return (s > s_ion * 1.001); });
		genCounts(magEPeak, magdEMag, [&s_mag](double s) { return (s < s_mag * 0.999); });

		return max;
	}
	//End Maxwellian

	//PPData
	PPData::PPData(Maxwellian maxwellian, std::vector<double> EBins, std::vector<double> PABins, std::string simDataDir, std::string particleName, std::string simBtmSat, std::string altUpgSat, std::string altDngSat) :
		energyBins{ EBins }, pitchBins{ PABins }
	{
		std::unique_ptr<Simulation> sim;
		SILENCE_COUT(SIM_API_EXCEP_CHECK(sim = std::make_unique<Simulation>(simDataDir)));

		SIM_API_EXCEP_CHECK(
			s_ion = sim->simMin();
			s_alt = sim->satellite(altUpgSat)->altitude();
			s_mag = sim->simMax();

			B_ion = sim->getBFieldAtS(s_ion, 0.0);
			B_alt = sim->getBFieldAtS(s_alt, 0.0);
			B_mag = sim->getBFieldAtS(s_mag, 0.0);

			mass = sim->particle(particleName)->mass();
			charge = sim->particle(particleName)->charge();

			int vparaind{ sim->particle(particleName)->getAttrIndByName("vpara") };
			int vperpind{ sim->particle(particleName)->getAttrIndByName("vperp") };
			int sind{ sim->particle(particleName)->getAttrIndByName("s") };

			initial = ParticleData(sim->particle(particleName)->data(true).at(vparaind), sim->particle(particleName)->data(true).at(vperpind), mass);
			initial.s_pos = sim->particle(particleName)->data(true).at(sind);
			bottom = ParticleData(sim->satellite(simBtmSat)->data().at(0).at(vparaind), sim->satellite(simBtmSat)->data().at(0).at(vperpind), mass);
			upward = ParticleData(sim->satellite(altUpgSat)->data().at(0).at(vparaind), sim->satellite(altUpgSat)->data().at(0).at(vperpind), mass);
			dnward = ParticleData(sim->satellite(altDngSat)->data().at(0).at(vparaind), sim->satellite(altDngSat)->data().at(0).at(vperpind), mass);
		); //end SIM_API_EXCEP_CHECK

		maxCounts = maxwellian.counts(initial, s_ion, s_mag, B_ion, B_alt, B_mag);
		std::cout << "PPData::PPData - multiplying counts by -cos(PA_src) for ion source, 1/cos(PA_sat) for mag source\n";
		for (unsigned int iii = 0; iii < maxCounts.size(); iii++) //isotropize counts -> 3D
		{
			if (initial.s_pos.at(iii) < s_ion * 1.001) //ionospheric source
				maxCounts.at(iii) *= -cos(initial.pitch.at(iii) * RADS_PER_DEG);
			else if (initial.s_pos.at(iii) > s_mag * 0.999)//magnetospheric source
				maxCounts.at(iii) *= 1.0 / cos(dnward.pitch.at(iii) * RADS_PER_DEG);
			else
				throw std::logic_error("PPData::PPData : particle is not ionospheric or magnetospheric source");
		}
	}
	//End PPData
}