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
	ParticleData::ParticleData()
	{
		//creates an empty ParticleData with vectors of size 0
	}

	ParticleData::ParticleData(vector<double>& v_para, vector<double>& v_perp, double mass) :
		vpara{ std::move(v_para) }, vperp{ std::move(v_perp) } //auto calculate E, Pitch
	{
		utils::numerical::v2DtoEPitch(vpara, vperp, mass, energy, pitch);
	}

	ParticleData::ParticleData(unsigned int size, bool EPA_only) //EPA_only defaults to true
	{
		energy.resize(size);
		pitch.resize(size);
		
		if (EPA_only) return;

		vpara.resize(size);
		vperp.resize(size);
		t_esc.resize(size);
		s_pos.resize(size);
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
	Maxwellian::Maxwellian(double dlogEdist) : dlogE_dist{ dlogEdist }
	{

	}

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

	vector<double> Maxwellian::counts(ParticleData& init, double s_ion, double s_mag)
	{
		vector<double> max(init.energy.size()); //array of maxwellian counts
		vector<double> maxtmp; //temporary holder allowing std::transform to be used twice instead of nested loops

		//multiply magnitudes by scaling factors
		//may need to be scaled for a number of reasons
		//ex: differing ranges of pitch angles for ion and mag sources
		std::transform(iondEMag.begin(), iondEMag.end(), iondEMag.begin(), [&](double dE) { return dE * ionModFactor; });
		std::transform(magdEMag.begin(), magdEMag.end(), magdEMag.begin(), [&](double dE) { return dE * magModFactor; });

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

		auto genCounts = [&](const vector<double>& EPeak, const vector<double>& dEMag, double modFactor, std::function<bool(double)> zero)
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
		genCounts(ionEPeak, iondEMag, ionModFactor, [&s_ion](double s) { return (s > s_ion * 1.001); });
		genCounts(magEPeak, magdEMag, magModFactor, [&s_mag](double s) { return (s < s_mag * 0.999); });

		return max;
	}
	//End Maxwellian


	//Bins
	Bins::Bins(vector<double>& E_bins, vector<double>& PA_bins) ://, dblVec2D& ind_1D) :
		E{ std::move(E_bins) }, PA{ std::move(PA_bins) }//, index_1D{ std::move(ind_1D) }
	{

	}

	Bins::Bins(const Bins& copy) //copy constructor
	{
		E = copy.E;
		PA = copy.PA;
		//index_1D.resize(copy.index_1D.size());
		//for (unsigned int cnt = 0; cnt < copy.index_1D.size(); cnt++)
			//index_1D.at(cnt) = copy.index_1D.at(cnt);
	}
	//End Bins
	

	//Ionosphere
	Ionosphere::Ionosphere(unsigned int numLayers, double s_max, double s_min)
	{
		h = vector<double>(numLayers);
		Z = vector<double>(numLayers);
		p = vector<double>(numLayers);
		s = vector<double>(numLayers);
		B = vector<double>(numLayers);

		for (unsigned int layer = 0; layer < numLayers; layer++) //endpoint inclusive
		{
			s.at(layer) = s_max - layer * (s_max - s_min) / (numLayers - 1); //in m
			h.at(layer) = 100 * (s_max - s_min) / (numLayers - 1); //in cm, hence * 100
		}
	}

	Ionosphere::Ionosphere(vector<double>& s_vec) : s{ std::move(s_vec) }
	{
		h = vector<double>(s_vec.size());
		Z = vector<double>(s_vec.size());
		p = vector<double>(s_vec.size());
		B = vector<double>(s_vec.size());
	}

	void Ionosphere::seth(double h_all)
	{
		if (h.at(0) != 0) std::cout << "Ionosphere::seth: Warning, array already has non-zero values.  Continuing\n";
		for (auto& h_layer : h)
			h_layer = h_all;
	}

	void Ionosphere::seth(vector<double>& h_vec)
	{
		if (h_vec.size() != s.size()) throw std::invalid_argument("Ionosphere::seth: h_vec.size does not match s.size");
		h = h_vec;
	}

	void Ionosphere::setZ(double Z_all)
	{
		for (auto& Z_layer : Z)
			Z_layer = Z_all;
	}

	void Ionosphere::setZ(vector<double>& Z_vec)
	{
		if (Z_vec.size() != s.size()) throw std::invalid_argument("Ionosphere::setZ: Z_vec.size does not match s.size");
		Z = Z_vec;
	}

	void Ionosphere::setp(double p_max, double p_min, bool log) //log defaults to true
	{
		if (log)
		{
			p_max = log10(p_max);
			p_min = log10(p_min);

			for (unsigned int layer = 0; layer < s.size(); layer++)
				p.at(layer) = pow(10.0, layer * (p_max - p_min) / ((double)s.size() - 1) + p_min); //start at low density (high alt) and go up in density (low alt)
		}
		else
		{
			for (unsigned int layer = 0; layer < s.size(); layer++)
				p.at(layer) = layer * (p_max - p_min) / ((double)s.size() - 1) + p_min;
		}
	}

	void Ionosphere::setp(vector<double>& p_vec)
	{
		if (p_vec.size() != s.size()) throw std::invalid_argument("Ionosphere::setp: p_vec.size does not match s.size");
		p = p_vec;
	}

	void Ionosphere::setB(BField* Bfield, double t)
	{
		for (unsigned int s_layer = 0; s_layer < s.size(); s_layer++)
		{
			B.at(s_layer) = Bfield->getBFieldAtS(s.at(s_layer), t);
		}
	}

	void Ionosphere::setB(vector<double>& B_vec)
	{
		if (B_vec.size() != s.size()) throw std::invalid_argument("Ionosphere::setp: B_vec.size does not match s.size");
		B = B_vec;
	}

	void Ionosphere::altToS(BField* B)
	{
		for (unsigned int s_ind = 0; s_ind < s.size(); s_ind++)
			s.at(s_ind) = B->getSAtAlt(s.at(s_ind));
	}
	//End Ionosphere


	//PPData
	PPData::PPData(Ionosphere& ionosphere, Maxwellian& maxwellian, const Bins& distBins, const Bins& satBins, string simDataDir, string particleName, string btmSatName, string upgSatName, string dngSatName) :
		ionsph{ std::move(ionosphere) }, distbins{ std::move(distBins) }, satbins{ std::move(satBins) }
	{
		std::unique_ptr<Simulation> sim;
		SILENCE_COUT(SIM_API_EXCEP_CHECK(sim = std::make_unique<Simulation>(simDataDir)));

		SIM_API_EXCEP_CHECK(
			s_ion = sim->simMin();
			s_sat = sim->satellite(upgSatName)->altitude();
			s_mag = sim->simMax();

			B_ion = sim->getBFieldAtS(s_ion, 0.0);
			B_sat = sim->getBFieldAtS(s_sat, 0.0);
			B_mag = sim->getBFieldAtS(s_mag, 0.0);

			mass = sim->particle(particleName)->mass();
			charge = sim->particle(particleName)->charge();

			int vparaind{ sim->particle(particleName)->getAttrIndByName("vpara") };
			int vperpind{ sim->particle(particleName)->getAttrIndByName("vperp") };
			int sind{ sim->particle(particleName)->getAttrIndByName("s") };

			initial = ParticleData(sim->particle(particleName)->__data(true).at(vparaind), sim->particle(particleName)->__data(true).at(vperpind), mass);
			initial.s_pos = sim->particle(particleName)->__data(true).at(sind);
			bottom = ParticleData(sim->satellite(btmSatName)->__data().at(0).at(vparaind), sim->satellite(btmSatName)->__data().at(0).at(vperpind), mass);
			upward = ParticleData(sim->satellite(upgSatName)->__data().at(0).at(vparaind), sim->satellite(upgSatName)->__data().at(0).at(vperpind), mass);
			dnward = ParticleData(sim->satellite(dngSatName)->__data().at(0).at(vparaind), sim->satellite(dngSatName)->__data().at(0).at(vperpind), mass);

			DipoleB dip(sim->Bmodel()->ILAT(), 1.0e-10, RADIUS_EARTH / 1000.0, false);
			ionsph.altToS(&dip);
			ionsph.setB(&dip, 0.0);
		); //end SIM_API_EXCEP_CHECK

		maxWeights = maxwellian.counts(initial, s_ion, s_mag);
	}
	//End PPData
}