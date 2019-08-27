#include <algorithm>
#include <iterator>


#include "utils/ionosphere.h"
#include "utils/numerical.h"
#include "utils/silenceStreamMacros.h"
#include "ErrorHandling/simExceptionMacros.h"

using utils::fileIO::ParticleDistribution;

namespace ionosphere
{
	//ParticleData
	ParticleData::ParticleData()
	{
		//creates an empty ParticleData with vectors of size 0
	}

	ParticleData::ParticleData(double_v1D& v_para, double_v1D& v_perp, double mass) :
		vpara{ std::move(v_para) }, vperp{ std::move(v_perp) } //auto calculate E, Pitch
	{
		utils::numerical::v2DtoEPitch(vpara, vperp, mass, energy, pitch);
	}

	ParticleData::ParticleData(size_t size, bool EPA_only) //EPA_only defaults to true
	{
		energy.resize(size);
		pitch.resize(size);
		
		if (EPA_only) return;

		vpara.resize(size);
		vperp.resize(size);
		t_esc.resize(size);
		s_pos.resize(size);
	}

	void ParticleData::clear()
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


	//MaxwellianSpecs
	MaxwellianSpecs::MaxwellianSpecs(double dlogEdist) : dlogE_dist{ dlogEdist }
	{

	}

	void MaxwellianSpecs::push_back_ion(eV E_peak, dEflux dE_magnitude, int partsAtE)
	{
		if (E_peak <= 0.0) throw std::logic_error("MaxwellianSpecs::push_back_ion: E_peak <= 0.0.  E_peak must be > 0.0 " + std::to_string(E_peak));
		if (dE_magnitude <= 0.0) throw std::logic_error("MaxwellianSpecs::push_back_ion: dE_mag <= 0.0.  dE_mag must be > 0.0 " + std::to_string(E_peak));
		if (partsAtE <= 0) throw std::logic_error("MaxwellianSpecs::push_back_ion: partsAtE <= 0.  Not physical. " + std::to_string(partsAtE));

		ionEPeak.push_back(E_peak);
		iondEMag.push_back(dE_magnitude / (double)partsAtE); //scales dE according to how many particles will be in the bin containing E
	}

	void MaxwellianSpecs::push_back_mag(eV E_peak, dEflux dE_magnitude, int partsAtE)
	{
		if (E_peak <= 0.0) throw std::logic_error("MaxwellianSpecs::push_back_mag: E_peak <= 0.0.  E_peak must be > 0.0 " + std::to_string(E_peak));
		if (dE_magnitude <= 0.0) throw std::logic_error("MaxwellianSpecs::push_back_mag: dE_mag <= 0.0.  dE_mag must be > 0.0 " + std::to_string(E_peak));
		if (partsAtE <= 0) throw std::logic_error("MaxwellianSpecs::push_back_mag: partsAtE <= 0.  Not physical. " + std::to_string(partsAtE));

		magEPeak.push_back(E_peak);
		magdEMag.push_back(dE_magnitude / (double)partsAtE);
	}

	dNflux_v1D MaxwellianSpecs::dNfluxAtE(ParticleData& init, meters s_ion, meters s_mag)
	{
		dNflux_v1D max(init.energy.size()); //array of maxwellian counts
		dNflux_v1D maxtmp; //temporary holder allowing std::transform to be used twice instead of nested loops

		//multiply magnitudes by scaling factors
		//may need to be scaled for a number of reasons
		//ex: differing ranges of pitch angles for ion and mag sources
		std::transform(iondEMag.begin(), iondEMag.end(), iondEMag.begin(), [&](double dE) { return dE * ionModFactor; });
		std::transform(magdEMag.begin(), magdEMag.end(), magdEMag.begin(), [&](double dE) { return dE * magModFactor; });

		auto count_E = [&](eV E_eval, eV E_peak, dEflux dEflux_peak, bool zero)
		{ //calculates count for a given E (E_eval) for a maxwellian peaked at E_peak with magnitude dEflux_peak
			if (E_peak <= 0.0)
				throw std::logic_error("MaxwellianSpecs::counts: invalid maxwellian peak E (le 0.0)");
			if (zero || E_eval == 0.0)
				return 0.0;
				
			//eV binWidth_E { pow(10, log10(E_eval) + 0.5 * dlogE_dist) - pow(10, log10(E_eval) - 0.5 * dlogE_dist) };
			//eV binWidth_kT{ pow(10, log10(E_peak) + 0.5 * dlogE_dist) - pow(10, log10(E_peak) - 0.5 * dlogE_dist) };
			//return exp(-E_eval / E_peak) * binWidth_E / E_eval * (dEflux_peak / (exp(-1.0) * binWidth_kT));
			return (dEflux_peak / (exp(-2.0) * E_peak * E_peak)) * (exp(-2.0 * E_eval / E_peak) * E_eval);
		};

		auto genCounts = [&](const double_v1D& E_peak, const dEflux_v1D& dEMag, double modFactor, std::function<bool(double)> zero)
		{ //iterates over particles and specified Peak/Magnitude values (ionospheric or magnetospheric) and adds values to "max" (after "maxtmp")
			for (size_t entr = 0; entr < E_peak.size(); entr++) //iterate over ionospheric maxwellian specifications
			{
				auto getCount = [&](double E, double s) { return count_E(E, E_peak.at(entr), dEMag.at(entr), zero(s)); };

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
	//End MaxwellianSpecs


	//Bins
	Bins::Bins(double_v1D& E_bins, degrees_v1D& PA_bins) :
		E{ std::move(E_bins) }, PA{ std::move(PA_bins) }
	{

	}

	Bins::Bins(const Bins& copy) //copy constructor
	{
		E = copy.E;
		PA = copy.PA;
	}
	//End Bins
	

	//IonosphereSpecs
	IonosphereSpecs::IonosphereSpecs(unsigned int numLayers, double s_max, double s_min)
	{
		s = double_v1D(numLayers + 1);
		h = double_v1D(numLayers + 1);
		B = double_v1D(numLayers + 1);

		for (unsigned int layer = 0; layer < numLayers + 1; layer++) //endpoint inclusive, adds one more at the bottom (sim needs)
		{
			s.at(layer) = s_max - layer * (s_max - s_min) / (numLayers - 1); //in m
			h.at(layer) = 100 * (s_max - s_min) / (numLayers - 1); //in cm, hence * 100
		}
	}

	void IonosphereSpecs::seth()
	{
		if (s.size() == 0)
			throw std::logic_error("IonosphereSpecs::seth: Error: s is not populated with values.");
		
		for (size_t layer = 0; layer < s.size() - 1; layer++)
			h.at(layer) = 100 * (s.at(layer) - s.at(layer + 1));
		
		h.at(s.size() - 1) = h.at(s.size() - 2);
	}

	void IonosphereSpecs::seth(double h_all)
	{
		bool warn{ false };
		for (auto& h_layer : h)
		{
			if (!warn && h.at(0) != 0)
			{
				std::cout << "IonosphereSpecs::seth: Warning: array already has non-zero values.  Continuing\n";
				warn = true;
			}
			h_layer = h_all;
		}
	}

	void IonosphereSpecs::setB(BField* Bfield, double t)
	{
		for (size_t s_layer = 0; s_layer < s.size(); s_layer++)
		{
			B.at(s_layer) = Bfield->getBFieldAtS(s.at(s_layer), t);
		}
	}

	void IonosphereSpecs::setB(double_v1D& B_vec)
	{
		if (B_vec.size() != s.size()) throw std::invalid_argument("IonosphereSpecs::setp: B_vec.size does not match s.size");
		B = B_vec;
	}

	void IonosphereSpecs::altToS(BField* B)
	{
		for (size_t s_ind = 0; s_ind < s.size(); s_ind++)
			s.at(s_ind) = B->getSAtAlt(s.at(s_ind));

		seth();
	}

	void IonosphereSpecs::addSpecies(string name, double Z_spec, function<double(double)> density_s)
	{
		if (names.size() != p.size()) throw std::logic_error("IonosphereSpecs::addSpecies: size of names and p vectors do not match - names, p: " +
			std::to_string(names.size()) + ", " + std::to_string(p.size()));

		names.push_back(name);
		Z.push_back(Z_spec);
		p.push_back(vector<double>(s.size()));

		for (size_t alt = 0; alt < p.back().size(); alt++)
			p.back().at(alt) = density_s(s.at(alt));
	}
	//End IonosphereSpecs


	//PPData
	EOMSimData::EOMSimData(IonosphereSpecs& ionosphere, MaxwellianSpecs& maxspecs, Bins& distribution, Bins& satellite, string dir_simdata, string name_particle, string name_btmsat, string name_upgsat, string name_dngsat) :
		ionsph{ std::move(ionosphere) }, distbins{ std::move(distribution) }, satbins{ std::move(satellite) }, datadir{ dir_simdata }
	{
		std::unique_ptr<Simulation> sim;
		SILENCE_COUT(SIM_API_EXCEP_CHECK(sim = std::make_unique<Simulation>(dir_simdata)));
		
		Particle* particle{ sim->particle(name_particle) };
		Satellite* sat_btm{ sim->satellite(name_btmsat) };
		Satellite* sat_dng{ sim->satellite(name_dngsat) };
		Satellite* sat_upg{ sim->satellite(name_upgsat) };

		SIM_API_EXCEP_CHECK(
			s_ion = sim->simMin();
			s_sat = sat_upg->altitude();
			s_mag = sim->simMax();

			B_ion = sim->getBFieldAtS(s_ion, 0.0);
			B_sat = sim->getBFieldAtS(s_sat, 0.0);
			B_mag = sim->getBFieldAtS(s_mag, 0.0);

			mass = particle->mass();
			charge = particle->charge();

			int vparaind{ particle->getAttrIndByName("vpara") };
			int vperpind{ particle->getAttrIndByName("vperp") };
			int sind{ particle->getAttrIndByName("s") };

			initial = ParticleData(particle->__data(true).at(vparaind), particle->__data(true).at(vperpind), mass);
			initial.s_pos = particle->__data(true).at(sind);
			bottom = ParticleData(sat_btm->__data().at(0).at(vparaind), sat_btm->__data().at(0).at(vperpind), mass);
			bottom.s_pos = std::move(sat_btm->__data().at(0).at(sind));
			upward = ParticleData(sat_upg->__data().at(0).at(vparaind), sat_upg->__data().at(0).at(vperpind), mass);
			upward.s_pos = std::move(sat_upg->__data().at(0).at(sind));
			dnward = ParticleData(sat_dng->__data().at(0).at(vparaind), sat_dng->__data().at(0).at(vperpind), mass);
			dnward.s_pos = std::move(sat_dng->__data().at(0).at(sind));

			DipoleB dip(sim->Bmodel()->ILAT(), 1.0e-10, RADIUS_EARTH / 1000.0, false);
			ionsph.altToS(&dip);
			ionsph.setB(&dip, 0.0);
		); //end SIM_API_EXCEP_CHECK

		maxwellian = maxspecs.dNfluxAtE(initial, s_ion, s_mag);

		//make sure EOMSimData has been properly formed
		if (s_ion <= 0.0) throw std::logic_error("EOMSimData::EOMSimData: s_ion <= 0.  Something is wrong. " + std::to_string(s_ion));
		if (s_sat <= 0.0) throw std::logic_error("EOMSimData::EOMSimData: s_sat <= 0.  Something is wrong. " + std::to_string(s_sat));
		if (s_mag <= 0.0) throw std::logic_error("EOMSimData::EOMSimData: s_mag <= 0.  Something is wrong. " + std::to_string(s_mag));

		if (B_ion >= 0.0) throw std::logic_error("EOMSimData::EOMSimData: B_ion >= 0.  Something is wrong. " + std::to_string(B_ion));
		if (B_sat >= 0.0) throw std::logic_error("EOMSimData::EOMSimData: B_sat >= 0.  Something is wrong. " + std::to_string(B_sat));
		if (B_mag >= 0.0) throw std::logic_error("EOMSimData::EOMSimData: B_mag >= 0.  Something is wrong. " + std::to_string(B_mag));

		if (mass <= 0.0)  throw std::logic_error("EOMSimData::EOMSimData: mass <= 0.  Something is wrong. " + std::to_string(mass));
		if (charge == 0.0) throw std::logic_error("EOMSimData::EOMSimData: charge = 0.  Something is wrong. " + std::to_string(charge));
	}
	//End PPData
}