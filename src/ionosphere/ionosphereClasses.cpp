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

	ParticleList::ParticleList(const double_v1D& v_para, const double_v1D& v_perp, double mass) :
		vpara{ v_para }, vperp{ v_perp }
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
	Bins::Bins() : E{ vector<eV>() }, PA{ vector<degrees>() } //protected
	{
		E_bounds = vector<binBounds<eV>>();
		PA_bounds = vector<binBounds<degrees>>();
	}

	Bins::Bins(const vector<eV>& E_midbins, const vector<degrees>& PA_midbins) :
		E{ E_midbins }, PA{ PA_midbins }
	{
		dlogE_bin = std::abs(log10(E_midbins.at(1)) - log10(E_midbins.at(0)));
		dangle_bin = std::abs(PA_midbins.at(1) - PA_midbins.at(0));

		E_ascending = E_midbins.back() > E_midbins.front();
		A_ascending = PA_midbins.back() > PA_midbins.front();

		for (size_t eny = 0; eny < E_midbins.size(); eny++)
		{
			E_bounds.push_back(binBounds<eV>());
			E_bounds.at(eny).mid = E_midbins.at(eny);

			if (E_ascending)
			{ //calculate bin boundaries, set next lower bound to prev upper bound
				E_bounds.at(eny).upper = pow(10, log10(E_bounds.at(eny).mid) + (0.5 * dlogE_bin));
				
				if (eny == 0) E_bounds.at(eny).lower = pow(10, log10(E_bounds.at(eny).mid) - (0.5 * dlogE_bin));
				else          E_bounds.at(eny).lower = E_bounds.at(eny - 1).upper;
			}
			else
			{ //calculate bin boundaries, set next upper bound to prev lower bound
				E_bounds.at(eny).lower = pow(10, log10(E_bounds.at(eny).mid) - (0.5 * dlogE_bin));
				
				if (eny == 0) E_bounds.at(eny).upper = pow(10, log10(E_bounds.at(eny).mid) + (0.5 * dlogE_bin));
				else          E_bounds.at(eny).upper = E_bounds.at(eny - 1).lower;
			}
		}

		for (size_t ang = 0; ang < PA_midbins.size(); ang++)
		{
			PA_bounds.push_back(binBounds<degrees>());
			PA_bounds.at(ang).mid = PA_midbins.at(ang);

			if (A_ascending)
			{ //calculate bin boundaries, set next lower bound to prev upper bound
				PA_bounds.at(ang).upper = PA_bounds.at(ang).mid + (0.5 * dangle_bin);

				if (ang == 0) PA_bounds.at(ang).lower = PA_bounds.at(ang).mid - (0.5 * dangle_bin);
				else          PA_bounds.at(ang).lower = PA_bounds.at(ang - 1).upper;
			}
			else
			{ //calculate bin boundaries, set next upper bound to prev lower bound
				PA_bounds.at(ang).lower = PA_bounds.at(ang).mid - (0.5 * dangle_bin);

				if (ang == 0) PA_bounds.at(ang).upper = PA_bounds.at(ang).mid + (0.5 * dangle_bin);
				else          PA_bounds.at(ang).upper = PA_bounds.at(ang - 1).lower;
			}
		}
	}

	Bins::Bins(const Bins& copy) : E{ copy.E }, PA{ copy.PA }, E_bounds{ copy.E_bounds }, PA_bounds{ copy.PA_bounds },
		dlogE_bin{ copy.dlogE_bin }, dangle_bin{ copy.dangle_bin }, E_ascending{ copy.E_ascending }, A_ascending{ copy.A_ascending }
	{
		
	}

	template <typename T>
	ParticlesBinned<T> Bins::binParticleList(const ParticleList& particles, const vector<T>& weightPerParticle, const function<bool(T)>& conditionCheck) const
	{
		//guards to check sizes of vectors are equal?
		//guards to check bins are equal size and ascending/descending in order?
		
		ParticlesBinned<T> ret;
		ret.binnedData = vector<vector<T>>(PA.size(), vector<T>(E.size()));
		ret.bins = *this;

		//overall boundaries
		eV      E_max{ E_ascending ? E_bounds.back().upper : E_bounds.front().upper };
		eV      E_min{ E_ascending ? E_bounds.front().lower : E_bounds.back().lower };
		degrees A_max{ A_ascending ? PA_bounds.back().upper : PA_bounds.front().upper };
		degrees A_min{ A_ascending ? PA_bounds.front().lower : PA_bounds.back().lower };

		int outsideLimits{ 0 };
		int FPE{ 0 };

		for (size_t part = 0; part < particles.energy.size(); part++) //iterate over particles
		{
			const eV      partEnerg{ particles.energy.at(part) };
			const degrees partPitch{ particles.pitch.at(part) };
			
			if ((partEnerg == 0.0 && partPitch == 0.0) || weightPerParticle.at(part) == 0.0)
				continue;             //guards - if particle E, PA is zero, it wasnt detected - just skip it
			if (partEnerg >= E_max || //if particle is outside E, PA measured limits - skip it
				partEnerg <  E_min ||
				partPitch >= A_max || //PA shouldn't be an issue, there just in case
				partPitch <  A_min)   //also, code counts how many are outside limits for diagnostics
			{
				if (partPitch >= 180.0 ||
					partPitch <  0.0)
				{
					throw logic_error("Bins::binParticleList: particle outside 0-180 degrees - partPitch: " + to_string(partPitch));
				}
				
				outsideLimits++;
				continue;
			}
			
			//calculate bin index for E, PA of particle
			size_t angbin{ static_cast<size_t>(partPitch / dangle_bin) }; //this should give the bin index
			size_t enybin{ static_cast<size_t>((log10(partEnerg) - log10(E_ascending ? E_bounds.front().lower : E_bounds.back().lower)) / dlogE_bin) }; //ditto
			if (!E_ascending) enybin = E.size() - 1 - enybin; //reverses the bin index if E is descending
			if (!A_ascending) angbin = PA.size() - 1 - angbin; //ditto for PA
			
			const auto putInBin = [&](const size_t abin, const size_t ebin)
			{
				if (partEnerg >= E_bounds.at(ebin).lower && //E bin min
					partEnerg <  E_bounds.at(ebin).upper && //E bin max
					partPitch >= PA_bounds.at(abin).lower && //A bin min
					partPitch <  PA_bounds.at(abin).upper)   //A bin max
				{
					try
					{
						conditionCheck(ret.binnedData.at(abin).at(ebin));
					}
					catch (std::exception& e)
					{
						cout << "Condition Check failed.  " << e.what() << "\n";
						exit(1);
					}

					ret.binnedData.at(abin).at(ebin) += weightPerParticle.at(part);

					return true;
				}
				else
				{
					return false;
				}
			};
			
			if (!putInBin(angbin, enybin)) //has side effect of putting in bin / running this lambda
			{
				/*stringstream ss;

				ss << "================  PitchA  ================\n\n";

				ss << std::fixed;
				ss << std::setprecision(16);
				ss << "\npartPitch: " << partPitch << "\n";
				ss << "dangle_Bin: " << dangle_bin << "\n";
				ss << "partPitch / dangle_bin: " << partPitch / dangle_bin << "\n";
				ss << "static_cast<size_t> ^^: " << static_cast<size_t>(partPitch / dangle_bin) << "\n\n";

				ss << "A_ascending: " << A_ascending << "\n";
				ss << "pa bin:                 " << angbin << "\n";
				ss << "pa at bin:              " << PA.at(angbin) << "\n";
				ss << "0.5 * dangle_bin:       " << 0.5 * dangle_bin << "\n\n";



				ss << "================  Energy  ================\n\n";

				ss << "\npartEnerg: " << partEnerg << "\n";
				ss << "dlogE_Bin: " << dlogE_bin << "\n";
				ss << "E_ascending ? E_bounds.front().lower : E_bounds.back().lower" << (E_ascending ? E_bounds.front().lower : E_bounds.back().lower) << "\n";
				ss << "log10(partEnerg) - log10(low) / dlogE_bin: " << (log10(partEnerg) - log10(E_ascending ? E_bounds.front().lower : E_bounds.back().lower)) / dlogE_bin << "\n";
				ss << "static_cast<size_t> ^^: " << static_cast<size_t>((log10(partEnerg) - log10(E_ascending ? E_bounds.front().lower : E_bounds.back().lower)) / dlogE_bin) << "\n\n";

				ss << "E_ascending: " << E_ascending << "\n";
				ss << "E bin:                 " << enybin << "\n";
				ss << "E  at bin:             " << E.at(enybin) << "\n";
				ss << "0.5 * dlogE_bin:       " << 0.5 * dlogE_bin << "\n\n";

				ss << "gte energy: " << partEnerg << " >= " << pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin);
				ss << "  " << (partEnerg >= pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin)) << "\n";
				ss << "lt  energy: " << partEnerg << " < " << pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin);
				ss << "  " << (partEnerg < pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin)) << "\n";
				ss << "gte pitch: " << partPitch << " >= " << PA.at(angbin) - 0.5 * dangle_bin;
				ss << "  " << (partPitch >= PA.at(angbin) - 0.5 * dangle_bin) << "\n";
				ss << "gte pitch diff: " << partPitch - (PA.at(angbin) - 0.5 * dangle_bin) << "\n";
				ss << "lt  pitch: " << partPitch << " < " << PA.at(angbin) + 0.5 * dangle_bin;
				ss << "  " << (partPitch < PA.at(angbin) + 0.5 * dangle_bin) << "\n";
				ss << "lt  pitch diff: " << (PA.at(angbin) + 0.5 * dangle_bin) - partPitch << "\n";

				ss << "index: " << part << "\n";
				ss << "ang: " << partPitch << "\n";
				ss << "eny: " << partEnerg << "\n";

				ss << "\n================================================================\n\n";
				
				cout << ss.str();

				exit(1);*/

				if (partEnerg < E_bounds.at(enybin).lower) enybin -= E_ascending ? 1 : -1;
				else if (partEnerg >= E_bounds.at(enybin).upper) enybin += E_ascending ? 1 : -1;

				if (partPitch < PA_bounds.at(angbin).lower) angbin -= A_ascending ? 1 : -1;
				else if (partPitch >= PA_bounds.at(angbin).upper) angbin += A_ascending ? 1 : -1;
				
				try
				{
					if (enybin > E.size() - 1)
						throw logic_error("enybin");
					if (angbin > PA.size() - 1)
						throw logic_error("angbin");
				}
				catch (std::exception& e)
				{
					cout << e.what() << "\n";
					//cout << ss.str();
					exit(1);
				}

				FPE++;
				
				if (!putInBin(angbin, enybin))
				{
					cout << std::fixed;
					cout << std::setprecision(16);
					cout << "\npartPitch: " << partPitch << "\n";
					cout << "dangle_Bin: " << dangle_bin << "\n";
					cout << "partPitch / dangle_bin: " << partPitch / dangle_bin << "\n";
					cout << "static_cast<size_t> ^^: " << static_cast<size_t>(partPitch / dangle_bin) << "\n\n";

					cout << "pa at bin:              " << PA.at(angbin) << "\n";
					cout << "0.5 * dangle_bin:       " << 0.5 * dangle_bin << "\n\n";

					cout << "pa at bin precise:      " << PA.at(angbin) << "\n\n\n";

					cout << "gte energy: " << partEnerg << " >= " << pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin);
					cout << "  " << (partEnerg >= pow(10, log10(E.at(enybin)) - 0.5 * dlogE_bin)) << "\n";
					cout << "lt  energy: " << partEnerg << " < " << pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin);
					cout << "  " << (partEnerg < pow(10, log10(E.at(enybin)) + 0.5 * dlogE_bin)) << "\n";
					cout << "gte pitch: " << partPitch << " >= " << PA.at(angbin) - 0.5 * dangle_bin;
					cout << "  " << (partPitch >= PA.at(angbin) - 0.5 * dangle_bin) << "\n";
					cout << "gte pitch diff: " << partPitch - (PA.at(angbin) - 0.5 * dangle_bin) << "\n";
					cout << "lt  pitch: " << partPitch << " < " << PA.at(angbin) + 0.5 * dangle_bin;
					cout << "  " << (partPitch < PA.at(angbin) + 0.5 * dangle_bin) << "\n";
					cout << "lt  pitch diff: " << (PA.at(angbin) + 0.5 * dangle_bin) - partPitch << "\n";

					cout << "index: " << part << "\n";
					cout << "ang: " << partPitch << "\n";
					cout << "eny: " << partEnerg << "\n";
					exit(1);
					throw logic_error("ionosphere::binning::binWeighted: Particle does not belong in bin identified for it.");
				}
			}
		}

		if (outsideLimits > 0) cout << "ionosphere::binning::binWeighted : Particles out of limits: " << outsideLimits << "\n";
		if (FPE > 0) cout << "ionosphere::binning::binWeighted : Particles encountering FPE: " << FPE << "\n";
		
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

		initial = ParticleList(particle->data(true).at(vparaind), particle->data(true).at(vperpind), mass);
		initial.s_pos = particle->data(true).at(sind);
		ionosph = ParticleList(sat_btm->data().at(vparaind), sat_btm->data().at(vperpind), mass);
		ionosph.s_pos = move(sat_btm->data().at(sind));
		upgoing = ParticleList(sat_upg->data().at(vparaind), sat_upg->data().at(vperpind), mass);
		upgoing.s_pos = move(sat_upg->data().at(sind));
		dngoing = ParticleList(sat_dng->data().at(vparaind), sat_dng->data().at(vperpind), mass);
		dngoing.s_pos = move(sat_dng->data().at(sind));

		//
		//
		//
		vector<size_t> idx;
		for (size_t part = 0; part < initial.energy.size(); part++)
			idx.push_back(part);

		//function<bool(size_t)> cond = [](size_t a) { if (a != 0) throw logic_error("Bins::binParticleList: idx already assigned - " + to_string(a)); return true; };
		//function<bool(size_t)> cond = [](size_t a) { return true; };
		initial1DToDistbins2DMapping = distbins.binParticleList(initial, idx);// , cond);
		//
		//
		//

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