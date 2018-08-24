#ifndef POSTPROCESS_HELPER_H
#define POSTPROCESS_HELPER_H

#include <vector>
#include <string>
#include "Simulation/Simulation.h"
#include "utils/writeIOclasses.h"

using std::string;
using std::vector;

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

		ParticleData() {} //empty constructor for making an empty ParticleData
		ParticleData(const vector<double>& v_para, const vector<double>& v_perp, double mass);
		void free();
	};

	struct DLLCLEXP Maxwellian
	{ //at some point, needs to take in distribution data to generate weights/etc
		double dlogE_dist;            //dlogE of the particle distribution (not bins)
		double ionModFactor{ 1.0 };   //for now, the user needs to adjust these manually
		double magModFactor{ 1.0 };   //eventually, I'd like to calculate them automatically

		vector<double> ionEPeak; //ionospheric source peak energy of maxwellian
		vector<double> iondEMag; //ionospheric source magnitude of peak

		vector<double> magEPeak; //same for magnetosphere
		vector<double> magdEMag; //basically x, y axes (top, bottom) of a maxwellian graph

		Maxwellian(double dlogEdist);

		void push_back_ion(double E_peak, double dE_magnitude, int partsAtE = 1);
		void push_back_mag(double E_peak, double dE_magnitude, int partsAtE = 1);
		
		vector<double> counts(const ParticleData& init, const double s_ion, const double s_mag);
	};

	struct DLLCLEXP PPData
	{
		double s_ion;
		double s_alt;
		double s_mag;

		double B_ion; //B Field strength at ionospheric source
		double B_alt; //B Field strength at satellite
		double B_mag; //B Field strength at magnetospheric source

		double mass;
		double charge;

		//Original Distribution Bins
		vector<double> distEBins;
		vector<double> distPABins;

		//Postprocessing Bins
		const vector<double> ppEBins;
		const vector<double> ppPABins;

		//Satellite and Maxwellian Data
		vector<double> maxCounts; //maxwellian weights - "number of particles represented by the particle at same index"
		ParticleData initial;     //initial particle data
		ParticleData bottom;      //data on particles that escape out the bottom
		ParticleData upward;      //data on particles that are upgoing at satellite
		ParticleData dnward;      //data on particles that are downgoing at satellite

		PPData(Maxwellian maxwellian, vector<double> ppEBins, vector<double> ppPABins, string simDataDir, string particleName, string simBtmSat, string upgoingAtAltSat, string dngoingAtAltSat);
	};
} //end namespace postprocess

#endif /* !POSTPROCESS_HELPER_H */