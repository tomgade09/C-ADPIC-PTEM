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
	struct DLLCLEXP ParticleData
	{
		double_v1D  vpara;
		double_v1D  vperp;
		double_v1D  energy;
		degrees_v1D pitch;
		double_v1D  t_esc;
		double_v1D  s_pos;

		ParticleData(); //empty constructor for making an empty ParticleData
		ParticleData(size_t size, bool EPA_only = true);
		ParticleData(double_v1D& v_para, double_v1D& v_perp, double mass);
		
		void clear();
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
		
		dNflux_v1D dNfluxAtE(ParticleData& init, meters s_ion, meters s_mag);
	};

	struct DLLCLEXP Bins
	{
		double_v1D  E;
		degrees_v1D PA;

		Bins(double_v1D& E_bins, degrees_v1D& PA_bins);
		Bins(const Bins& copy); //copy constructor
	};

	struct DLLCLEXP IonosphereSpecs
	{
		double_v1D s;         //layer altitude (top of layer, in m, along field line)
		double_v1D h;         //layer height (in cm)
		double_v1D B;         //magnetic field strength at each layer

		vector<string> names; //name - i.e. H, He, O2, etc
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
		meters s_ion;      //distance **along the field line** (in m from Re) representing the top limit of the ionosphere (particle source)
		meters s_sat;      //distance **along the field line** (in m from Re) where the satellite resides
		meters s_mag;      //distance **along the field line** (in m from Re) representing outer limit of the sim

		tesla  B_ion;      //B Field strength at ionospheric source
		tesla  B_sat;      //B Field strength at satellite
		tesla  B_mag;      //B Field strength at magnetospheric source

		kg      mass;
		coulomb charge;

		Bins distbins;     //Original Distribution Bins - represents the binning of the equ of motion simulation
		Bins satbins;      //Satellite Bins - represents how the data is binned and output from the satellite
		IonosphereSpecs ionsph; //Ionosphere Parameters (for multi-level scattering)

		//Satellite and Maxwellian Data
		dEflux_v1D   maxwellian; //maxwellian weights - "number of particles represented by the particle at same index"
		ParticleData initial;    //initial particle data
		ParticleData bottom;     //data on particles that escape out the bottom
		ParticleData upward;     //data on particles that are upgoing at satellite
		ParticleData dnward;     //data on particles that are downgoing at satellite

		string datadir;          //directory where sim data is stored

		EOMSimData(IonosphereSpecs& ionosphere, MaxwellianSpecs& maxspecs, Bins& distribution, Bins& satellite, string dir_simdata, string name_particle, string name_btmsat, string name_upgsat, string name_dngsat);
	};

} //end namespace ionosphere

#endif /* !IONOSPHERE_CLASSES_H */