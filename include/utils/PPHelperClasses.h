#ifndef POSTPROCESS_HELPER_H
#define POSTPROCESS_HELPER_H

#include <vector>
#include <string>
#include "Simulation/Simulation.h"
#include "utils/writeIOclasses.h"

using std::string;
using std::vector;

typedef vector<vector<double>> dblVec2D;
typedef vector<double> dblVec;


namespace postprocess
{
	struct DLLCLEXP ParticleData
	{
		vector<double> vpara;
		vector<double> vperp;
		vector<double> energy;
		vector<double> pitch;
		vector<double> t_esc;
		vector<double> s_pos;

		ParticleData(); //empty constructor for making an empty ParticleData
		ParticleData(unsigned int size, bool EPA_only = true);
		ParticleData(vector<double>& v_para, vector<double>& v_perp, double mass);
		void free();
	};

	struct DLLCLEXP Maxwellian
	{ //at some point, needs to take in distribution data to generate weights/etc
		double dlogE_dist;          //dlogE of the particle distribution (not bins)
		double ionModFactor{ 1.0 }; //for now, the user needs to adjust these manually
		double magModFactor{ 1.0 }; //eventually, I'd like to calculate them automatically

		vector<double> ionEPeak;    //ionospheric source peak energy of maxwellian
		vector<double> iondEMag;    //ionospheric source magnitude of peak

		vector<double> magEPeak;    //same for magnetosphere
		vector<double> magdEMag;    //basically x, y axes (top, bottom) of a maxwellian graph dE vs E

		Maxwellian(double dlogEdist);

		void push_back_ion(double E_peak, double dE_magnitude, int partsAtE = 1);
		void push_back_mag(double E_peak, double dE_magnitude, int partsAtE = 1);
		
		vector<double> counts(ParticleData& init, double s_ion, double s_mag);
	};

	struct DLLCLEXP Bins
	{
		vector<double> E;
		vector<double> PA;
		dblVec2D       counts;

		Bins(vector<double>& E_bins, vector<double>& PA_bins, dblVec2D& cnt_bins = vector<vector<double>>());
		Bins(const Bins& copy); //copy constructor
	};

	struct DLLCLEXP Ionosphere
	{
		vector<double> h;      //layer height (in cm)
		vector<double> Z;      //atomic number of dominant species/weighted average
		vector<double> p;      //density
		vector<double> s;      //layer altitude (top of layer, in m, along field line)
		vector<double> B;      //magnetic field strength at each layer

		Ionosphere(unsigned int numLayers, double s_max, double s_min);
		Ionosphere(vector<double>& s_vec);

		void seth(double h_all);
		void seth(vector<double>& h_vec);
		void setZ(double Z_all);
		void setZ(vector<double>& Z_vec);
		void setp(double p_max, double p_min, bool log = true);
		void setp(vector<double>& p_vec);
		void setB(BField* B, double t);
		void setB(vector<double>& B_vec);
	};

	struct DLLCLEXP PPData
	{
		double s_ion; //distance **along the field line** (in m from Re) representing the top limit of the ionosphere (particle source)
		double s_sat; //distance **along the field line** (in m from Re) where the satellite resides
		double s_mag; //distance **along the field line** (in m from Re) representing outer limit of the sim

		double B_ion; //B Field strength at ionospheric source
		double B_sat; //B Field strength at satellite
		double B_mag; //B Field strength at magnetospheric source

		double mass;
		double charge;

		Bins distbins; //Original Distribution Bins - represents the binning of the equ of motion simulation
		Bins satbins; //Satellite Bins - represents how the data is binned and output from the satellite
		Ionosphere ionsph; //Ionosphere Parameters (for multi-level scattering)

		//Satellite and Maxwellian Data
		vector<double> maxWeights; //maxwellian weights - "number of particles represented by the particle at same index"
		ParticleData   initial;    //initial particle data
		ParticleData   bottom;     //data on particles that escape out the bottom
		ParticleData   upward;     //data on particles that are upgoing at satellite
		ParticleData   dnward;     //data on particles that are downgoing at satellite

		PPData(const Ionosphere& ionosphere, Maxwellian& maxwellian, const Bins& distBins, const Bins& satBins, string simDataDir, string particleName, string btmSatName, string upgSatName, string dngSatName);
	};

} //end namespace postprocess

#endif /* !POSTPROCESS_HELPER_H */