#ifndef IONOSPHERE_CLASSES_H
#define IONOSPHERE_CLASSES_H

#include <vector>
#include <string>
#include "Simulation/Simulation.h"
#include "utils/writeIOclasses.h"
#include "ionosphere/ionosphereUtils.h"

using std::string;
using std::vector;
using std::function;

namespace ionosphere
{
	struct DLLCLEXP Bins; //need this forward decls because I use it in ParticlesBinned

	struct DLLCLEXP ParticleList //1D particle data structure
	{
		vector<mpers>   vpara;
		vector<mpers>   vperp;
		vector<eV>      energy;
		vector<degrees> pitch;
		vector<seconds> t_esc;
		vector<meters>  s_pos;

		ParticleList(); //empty constructor for making an empty ParticleList
		ParticleList(size_t size, bool EPA_only = true);
		ParticleList(const double_v1D& v_para, const double_v1D& v_perp, double mass);

		void clear();
	};


	template <typename T>
	struct DLLCLEXP ParticlesBinned //2D particle data (binned) structure
	{
		Bins bins; //bins corresponding to the 2D data below
		vector<vector<T>> binnedData;

		ParticlesBinned<T>& operator+=(const ParticlesBinned<T>& rhs);
	};


	struct DLLCLEXP Bins
	{ //defines discrete bins (by the center value of that bin)
	protected:
		Bins(); //usually don't want to be able to create an empty bin, but it's necessary for an empty ParticleList

		//friend class ParticleList;
		template <typename T>
		friend class ParticlesBinned;

		template <typename T>
		struct binBounds
		{
			T lower{ 0.0 };
			T mid{ 0.0 };
			T upper{ 0.0 };
		};

	public:
		//bin boundaries
		//vector<eV>      E_min;
		//vector<degrees> PA_min;
		vector<eV>      E;
		vector<degrees> PA;
		//vector<eV>      E_max;
		//vector<degrees> PA_max;

		vector<binBounds<eV>>      E_bounds;
		vector<binBounds<degrees>> PA_bounds;

		//dist between bins
		eV      dlogE_bin{ 0.0 };
		degrees dangle_bin{ 0.0 };

		//bins ascending or descending?
		bool E_ascending{ false };
		bool A_ascending{ false };

		//maybe add machinery to (force) create the bins through this class instead feeding in preformed vectors
		//this would ensure distance between bins is even

		Bins(const vector<eV>& E_bins, const vector<degrees>& PA_bins);
		Bins(const Bins& copy);

		bool operator==(const Bins& other) const;
		bool operator!=(const Bins& other) const;

		template <typename T>
		ParticlesBinned<T> binParticleList(const ParticleList& particles, const vector<T>& weightPerParticle,
			const function<bool(T)>& conditionCheck = [](T a) { return true; }) const;
	};


	struct DLLCLEXP MaxwellianSpecs
	{ //at some point, needs to take in distribution data to generate weights/etc
		double dlogE_dist;          //dlogE of the particle distribution (not bins)
		double ionModFactor{ 1.0 }; //for now, the user needs to adjust these manually
		double magModFactor{ 1.0 }; //eventually, I'd like to calculate them automatically

		double_v1D ionEPeak;    //ionospheric source peak energy of maxwellian
		dEflux_v1D iondEMag;    //ionospheric source magnitude of peak

		double_v1D magEPeak;    //same for magnetosphere
		dEflux_v1D magdEMag;    //basically x, y axes (top, bottom) of a maxwellian graph dE vs E

		MaxwellianSpecs(double dlogEdist);

		void push_back_ion(eV E_peak, dEflux dE_magnitude, int partsAtE);
		void push_back_mag(eV E_peak, dEflux dE_magnitude, int partsAtE);
		
		dNflux_v1D dNfluxAtE(ParticleList& init, meters s_ion, meters s_mag);
	};


	struct DLLCLEXP IonosphereSpecs
	{
		double_v1D s;         //layer altitude (top of layer, in m, along field line)
		double_v1D h;         //layer height (in cm)
		double_v1D B;         //magnetic field strength at each layer

		vector<string> names; //atomic species name - i.e. H, He, O2, etc
		double_v1D Z;         //atomic number
		double_v2D p;         //density - outer = species, inner = level

		IonosphereSpecs(int numLayers, double s_max, double s_min);

		void seth(double h_all);
		void setB(BModel* B, double t);
		void setB(double_v1D& B_vec);
		void altToS(BModel* B);

		void addSpecies(string name, double Z, function<double(double)> density_s);

	private:
		void seth();
	};


	struct DLLCLEXP EOMSimData
	{
		unique_ptr<Simulation> sim{ nullptr };

		string datadir;          //directory where sim data is stored

		int qspsCount{ 0 };

		meters  s_ion{ 0.0 };    //distance **along the field line** (in m from Re) representing the top limit of the ionosphere (particle source)
		meters  s_sat{ 0.0 };    //distance **along the field line** (in m from Re) where the satellite resides
		meters  s_mag{ 0.0 };    //distance **along the field line** (in m from Re) representing outer limit of the sim

		tesla   B_ion{ 0.0 };    //B Field strength at ionospheric source
		tesla   B_sat{ 0.0 };    //B Field strength at satellite
		tesla   B_mag{ 0.0 };    //B Field strength at magnetospheric source

		kg      mass{ 0.0 };
		coulomb charge{ 0.0 };

		//eventually get rid of distbins and replace with the bins held in ParticleList
		Bins distbins;           //Original Distribution Bins - represents the binning of the equ of motion simulation
		Bins satbins;            //Satellite Bins - represents how the data is binned and output from the satellite
		IonosphereSpecs ionsph;  //Ionosphere Parameters (for multi-level scattering)

		dEflux_v1D maxwellian;   //maxwellian weights - "number of particles represented by the particle at same index"

		//Satellite Data
		ParticleList initial;    //initial particle data
		ParticleList ionosph;    //data of particles that escape to the ionosphere (bottom of PTEM)
		ParticleList upgoing;    //data of particles that are upgoing at satellite
		ParticleList dngoing;    //data of particles that are downgoing at satellite

		ParticlesBinned<size_t> initial1DToDistbins2DMapping;

		EOMSimData(IonosphereSpecs& ionosphere, MaxwellianSpecs& maxspecs, Bins& distribution, Bins& satellite, string dir_simdata, string name_particle, string name_btmsat, string name_upgsat, string name_dngsat);
	};

} //end namespace ionosphere

#endif /* !IONOSPHERE_CLASSES_H */