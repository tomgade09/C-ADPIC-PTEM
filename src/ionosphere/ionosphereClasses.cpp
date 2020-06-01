#include <algorithm>
#include <iterator>


#include "ionosphere/ionosphere.h"
#include "utils/numerical.h"
#include "utils/silenceStreamMacros.h"
#include "ErrorHandling/simExceptionMacros.h"

using std::clog;
using std::move;
using std::make_unique;
using std::invalid_argument;
using utils::fileIO::ParticleDistribution;

namespace ionosphere
{
	//
	//ParticleList
	//
	ParticleList::ParticleList()
	{
		//creates an empty ParticleList with vectors of size 0
	}

	ParticleList::ParticleList(double_v1D& v_para, double_v1D& v_perp, double mass) :
		vpara{ move(v_para) }, vperp{ move(v_perp) } //auto calculate E, Pitch
	{
		utils::numerical::v2DtoEPitch(vpara, vperp, mass, energy, pitch);
	}

	ParticleList::ParticleList(size_t size, bool EPA_only) //EPA_only defaults to true
	{
		energy.resize(size);
		pitch.resize(size);

		if (EPA_only) return;

		vpara.resize(size);
		vperp.resize(size);
		t_esc.resize(size);
		s_pos.resize(size);
	}

	void ParticleList::clear()
	{
		vpara.clear(); //clear data in vectors
		vperp.clear();
		energy.clear();
		pitch.clear();
		t_esc.clear();
		s_pos.clear();
	}
	//
	//End ParticleList
	//


	//
	//ParticlesBinned<T>
	//
	template <typename T>
	ParticlesBinned<T>& ParticlesBinned<T>::operator+=(const ParticlesBinned<T>& rhs)
	{
		if (bins != rhs.bins)
			throw invalid_argument("ParticlesBinned::operator+=: Bins from lhs and rhs\
				of the operator are not equal and therefore these instances cannot be added.");

		for (size_t outer = 0; outer < binnedData.size(); outer++)
			for (size_t inner = 0; inner < binnedData.at(outer).size(); inner++)
				binnedData.at(outer).at(inner) += rhs.binnedData.at(outer).at(inner);

		return *this;
	}

	//template instantiation for operator+=() - covers dNflux, dEflux
	template
	ParticlesBinned<double>& ParticlesBinned<double>::operator+=(const ParticlesBinned<double>& rhs);
	//
	//End ParticlesBinned<T>
	//


	//
	//Bins
	//
	Bins::Bins() //protected
	{
		E = vector<eV>();
		PA = vector<degrees>();
	}

	Bins::Bins(vector<eV>& E_bins, vector<degrees>& PA_bins) :
		E{ E_bins }, PA{ PA_bins }
	{

	}

	Bins::Bins(const Bins& copy)
	{
		E = copy.E;
		PA = copy.PA;
	}

	template <typename T>
	ParticlesBinned<T> Bins::binParticleList(const ParticleList& particles, const vector<T>& weightPerParticle, const function<bool(T)>& conditionCheck) const
	{
		//guards to check sizes of vectors are equal?
		//guards to check bins are equal size and ascending/descending in order?
		
		ParticlesBinned<T> ret;
		ret.binnedData = vector<vector<T>>(PA.size(), vector<T>(E.size()));

		bool Eascending{ (E.back() > E.front()) };   //determines whether or not bin E's ascend from less E to more E as ind increases
		bool Aascending{ (PA.back() > PA.front()) }; //same for angle

		double Emax{ Eascending ? E.back() : E.front() };
		double Emin{ !Eascending ? E.back() : E.front() };
		double Amax{ Aascending ? PA.back() : PA.front() };
		double Amin{ !Aascending ? PA.back() : PA.front() };

		double dlogE_bin{ std::abs(log10(E.at(1)) - log10(E.at(0))) };
		double dangle_bin{ std::abs(PA.at(1) - PA.at(0)) };

		int outsideLimits{ 0 };

		for (size_t part = 0; part < particles.energy.size(); part++) //iterate over particles
		{
			eV      partEnerg{ particles.energy.at(part) };
			degrees partPitch{ particles.pitch.at(part) };

			if ((partEnerg == 0.0 && partPitch == 0.0) || weightPerParticle.at(part) == 0.0)
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
			size_t angbin{ static_cast<size_t>(std::floor(partPitch / dangle_bin)) }; //this should give the bin index
			size_t enybin{ static_cast<size_t>(std::floor((log10(partEnerg) - (log10(Emin) - 0.5 * dlogE_bin)) / dlogE_bin)) }; //ditto
			if (!Eascending) enybin = E.size() - 1 - enybin; //reverses the bin index if E is descending
			if (!Aascending) angbin = PA.size() - 1 - angbin; //ditto for PA

			if (partEnerg >= pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin) && //E bin min
				partEnerg <  pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin) && //E bin max
				partPitch >= PA.at(angbin) - 0.5 * dangle_bin && //A bin min
				partPitch <  PA.at(angbin) + 0.5 * dangle_bin) //A bin max
			{
				try
				{
					conditionCheck(ret.binnedData.at(angbin).at(enybin));
				}
				catch (std::exception& e)
				{
					eV      oldPartEnerg{ particles.energy.at(ret.binnedData.at(angbin).at(enybin)) };
					degrees oldPartPitch{ particles.pitch.at(ret.binnedData.at(angbin).at(enybin)) };
					size_t  oldangbin{ static_cast<size_t>(std::floor(oldPartPitch / dangle_bin)) };
					size_t  oldenybin{ static_cast<size_t>(std::floor((log10(oldPartEnerg) - (log10(Emin) - 0.5 * dlogE_bin)) / dlogE_bin)) };
					if (!Eascending) oldenybin = E.size() - 1 - oldenybin; //reverses the bin index if E is descending
					if (!Aascending) oldangbin = PA.size() - 1 - oldangbin; //ditto for PA

					cout << e.what() << "\n";
					cout << "part idx:   " << part << "\n";
					cout << "bin val:    " << ret.binnedData.at(angbin).at(enybin) << "\n";

					cout << "\n======== particle in bin  ========\n";
					cout << "part pitch: " << oldPartPitch << ", part eng: " << oldPartEnerg << "\n";
					cout << "ang bin:   " << oldangbin << "\n";
					cout << "eny bin:   " << oldenybin << "\n";
					cout << "ang min:   " << PA.at(oldangbin) - 0.5 * dangle_bin << ", ang max: "
						<< PA.at(oldangbin) + 0.5 * dangle_bin << "\n";
					cout << "eny min:   " << pow(10, log10(E.at(oldenybin)) - 0.5 * dlogE_bin) << ", eny max: " 
						<< pow(10, log10(E.at(oldenybin)) + 0.5 * dlogE_bin) << "\n";
					
					cout << "\n======== current particle ========\n";
					cout << "part pitch: " << partPitch << ", part eng: " << partEnerg << "\n";
					cout << "ang bin:    " << angbin << "\n";
					cout << "eny bin:    " << enybin << "\n";
					cout << "ang min:   " << PA.at(angbin) - 0.5 * dangle_bin << ", ang max: "
						<< PA.at(angbin) + 0.5 * dangle_bin << "\n";
					cout << "eny min:   " << pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin) << ", eny max: "
						<< pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin) << "\n";
					exit(1);
				}

				ret.binnedData.at(angbin).at(enybin) += weightPerParticle.at(part);
			}
			else //this shouldn't ever execute, guards should prevent zero and out of limits values
				throw logic_error("ionosphere::binning::binWeighted: Particle does not belong in bin identified for it.");
		}

		if (outsideLimits > 0) clog << "ionosphere::binning::binWeighted : Particles out of limits: " << outsideLimits << "\n";

		return ret;
	}

	//template instantiation for binParticleList - double covers dNflux, dEflux
	template
	ParticlesBinned<double> Bins::binParticleList(const ParticleList& particles, const vector<double>& weightPerParticle, const function<bool(double)>& conditionCheck) const;

	bool Bins::operator==(const Bins& other) const
	{
		return (E == other.E && PA == other.PA);
	}

	bool Bins::operator!=(const Bins& other) const
	{
		return !(*this == other);
	}
	//
	//End Bins
	//


	//
	//MaxwellianSpecs
	//
	MaxwellianSpecs::MaxwellianSpecs(double dlogEdist) : dlogE_dist{ dlogEdist }
	{

	}

	void MaxwellianSpecs::push_back_ion(eV E_peak, dEflux dE_magnitude, int partsAtE)
	{
		if (E_peak <= 0.0) throw logic_error("MaxwellianSpecs::push_back_ion: E_peak <= 0.0.  E_peak must be > 0.0 " + to_string(E_peak));
		if (dE_magnitude <= 0.0) throw logic_error("MaxwellianSpecs::push_back_ion: dE_mag <= 0.0.  dE_mag must be > 0.0 " + to_string(E_peak));
		if (partsAtE <= 0) throw logic_error("MaxwellianSpecs::push_back_ion: partsAtE <= 0.  Not physical. " + to_string(partsAtE));

		ionEPeak.push_back(E_peak);
		iondEMag.push_back(dE_magnitude / (double)partsAtE); //scales dE according to how many particles will be in the bin containing E
	}

	void MaxwellianSpecs::push_back_mag(eV E_peak, dEflux dE_magnitude, int partsAtE)
	{
		if (E_peak <= 0.0) throw logic_error("MaxwellianSpecs::push_back_mag: E_peak <= 0.0.  E_peak must be > 0.0 " + to_string(E_peak));
		if (dE_magnitude <= 0.0) throw logic_error("MaxwellianSpecs::push_back_mag: dE_mag <= 0.0.  dE_mag must be > 0.0 " + to_string(E_peak));
		if (partsAtE <= 0) throw logic_error("MaxwellianSpecs::push_back_mag: partsAtE <= 0.  Not physical. " + to_string(partsAtE));

		magEPeak.push_back(E_peak);
		magdEMag.push_back(dE_magnitude / (double)partsAtE);
	}

	dNflux_v1D MaxwellianSpecs::dNfluxAtE(ParticleList& init, meters s_ion, meters s_mag)
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
				throw logic_error("MaxwellianSpecs::counts: invalid maxwellian peak E (le 0.0)");
			if (zero || E_eval == 0.0)
				return 0.0;
				
			//eV binWidth_E { pow(10, log10(E_eval) + 0.5 * dlogE_dist) - pow(10, log10(E_eval) - 0.5 * dlogE_dist) };
			//eV binWidth_kT{ pow(10, log10(E_peak) + 0.5 * dlogE_dist) - pow(10, log10(E_peak) - 0.5 * dlogE_dist) };
			//return exp(-E_eval / E_peak) * binWidth_E / E_eval * (dEflux_peak / (exp(-1.0) * binWidth_kT));
			return (dEflux_peak / (exp(-2.0) * E_peak * E_peak)) * (exp(-2.0 * E_eval / E_peak) * E_eval);
		};

		auto genCounts = [&](const double_v1D& E_peak, const dEflux_v1D& dEMag, double modFactor, function<bool(double)> zero)
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
	//
	//End MaxwellianSpecs
	//

	
	//
	//IonosphereSpecs
	//
	IonosphereSpecs::IonosphereSpecs(int numLayers, double s_max, double s_min)
	{
		s = double_v1D(numLayers + 1);
		h = double_v1D(numLayers + 1);
		B = double_v1D(numLayers + 1);

		for (int layer = 0; layer < numLayers + 1; layer++) //endpoint inclusive, adds one more at the bottom (sim needs)
		{
			s.at(layer) = s_max - layer * (s_max - s_min) / (numLayers - 1); //in m
			h.at(layer) = 100 * (s_max - s_min) / (numLayers - 1); //in cm, hence * 100
		}
	}

	void IonosphereSpecs::seth()
	{
		if (s.size() == 0)
			throw logic_error("IonosphereSpecs::seth: Error: s is not populated with values.");
		
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
				cout << "IonosphereSpecs::seth: Warning: array already has non-zero values.  Continuing\n";
				warn = true;
			}
			h_layer = h_all;
		}
	}

	void IonosphereSpecs::setB(BModel* BModel, double t)
	{
		for (size_t s_layer = 0; s_layer < s.size(); s_layer++)
		{
			B.at(s_layer) = BModel->getBFieldAtS(s.at(s_layer), t);
		}
	}

	void IonosphereSpecs::setB(double_v1D& B_vec)
	{
		if (B_vec.size() != s.size()) throw invalid_argument("IonosphereSpecs::setp: B_vec.size does not match s.size");
		B = B_vec;
	}

	void IonosphereSpecs::altToS(BModel* B)
	{
		for (size_t s_ind = 0; s_ind < s.size(); s_ind++)
			s.at(s_ind) = B->getSAtAlt(s.at(s_ind));

		seth();
	}

	void IonosphereSpecs::addSpecies(string name, double Z_spec, function<double(double)> density_s)
	{
		if (names.size() != p.size()) throw logic_error("IonosphereSpecs::addSpecies: size of names and p vectors do not match - names, p: " +
			to_string(names.size()) + ", " + to_string(p.size()));

		names.push_back(name);
		Z.push_back(Z_spec);
		p.push_back(vector<double>(s.size()));

		for (size_t alt = 0; alt < p.back().size(); alt++)
			p.back().at(alt) = density_s(s.at(alt));
	}
	//
	//End IonosphereSpecs
	//


	//
	//PPData
	//
	EOMSimData::EOMSimData(IonosphereSpecs& ionosphere, MaxwellianSpecs& maxspecs, Bins& distribution, Bins& satellite, string dir_simdata, string name_particle, string name_btmsat, string name_upgsat, string name_dngsat) :
		ionsph{ move(ionosphere) }, distbins{ move(distribution) }, satbins{ move(satellite) }, datadir{ dir_simdata }
	{
		/*SILENCE_COUT(*/SIM_API_EXCEP_CHECK(sim = move(make_unique<Simulation>(dir_simdata)))/*)*/;

		Particles* particle{ sim->particles(name_particle) };
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

		size_t vparaind{ particle->getAttrIndByName("vpara") };
		size_t vperpind{ particle->getAttrIndByName("vperp") };
		size_t sind{ particle->getAttrIndByName("s") };
		
		initial = ParticleList(particle->__data(true).at(vparaind), particle->__data(true).at(vperpind), mass);
		initial.s_pos = particle->__data(true).at(sind);
		ionosph = ParticleList(sat_btm->__data().at(vparaind), sat_btm->__data().at(vperpind), mass);
		ionosph.s_pos = move(sat_btm->__data().at(sind));
		upgoing = ParticleList(sat_upg->__data().at(vparaind), sat_upg->__data().at(vperpind), mass);
		upgoing.s_pos = move(sat_upg->__data().at(sind));
		dngoing = ParticleList(sat_dng->__data().at(vparaind), sat_dng->__data().at(vperpind), mass);
		dngoing.s_pos = move(sat_dng->__data().at(sind));
		
		vector<size_t> idx;
		for (size_t part = 0; part < initial.energy.size(); part++)
			idx.push_back(part);

		function<bool(size_t)> cond = [](size_t a) { if (a != 0) throw logic_error("Bins::binParticleList: idx already assigned - " + to_string(a)); return true; };
		initial1DToDistbins2DMapping = distbins.binParticleList(initial, idx, cond);

		ionsph.altToS(sim->Bmodel());
		ionsph.setB(sim->Bmodel(), 0.0);
		qspsCount = sim->Efield()->qspsCount();
		); //end SIM_API_EXCEP_CHECK
		
		maxwellian = maxspecs.dNfluxAtE(initial, s_ion, s_mag);

		//make sure EOMSimData has been properly formed
		if (s_ion <= 0.0) throw logic_error("EOMSimData::EOMSimData: s_ion <= 0.  Something is wrong. " + to_string(s_ion));
		if (s_sat <= 0.0) throw logic_error("EOMSimData::EOMSimData: s_sat <= 0.  Something is wrong. " + to_string(s_sat));
		if (s_mag <= 0.0) throw logic_error("EOMSimData::EOMSimData: s_mag <= 0.  Something is wrong. " + to_string(s_mag));

		if (B_ion >= 0.0) throw logic_error("EOMSimData::EOMSimData: B_ion >= 0.  Something is wrong. " + to_string(B_ion));
		if (B_sat >= 0.0) throw logic_error("EOMSimData::EOMSimData: B_sat >= 0.  Something is wrong. " + to_string(B_sat));
		if (B_mag >= 0.0) throw logic_error("EOMSimData::EOMSimData: B_mag >= 0.  Something is wrong. " + to_string(B_mag));

		if (mass <= 0.0)  throw logic_error("EOMSimData::EOMSimData: mass <= 0.  Something is wrong. " + to_string(mass));
		if (charge == 0.0) throw logic_error("EOMSimData::EOMSimData: charge = 0.  Something is wrong. " + to_string(charge));
	}
	//
	//End PPData
	//
}