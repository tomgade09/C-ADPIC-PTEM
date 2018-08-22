#ifndef POSTPROCESS_HELPER_H
#define POSTPROCESS_HELPER_H

#include <vector>
#include <string>
#include "Simulation/Simulation.h"
#include "utils/writeIOclasses.h"

namespace postprocess
{
	struct DLLCLEXP ParticleData
	{
		std::vector<double> vpara;
		std::vector<double> vperp;
		std::vector<double> energy;
		std::vector<double> pitch;
		std::vector<double> t_esc;
		std::vector<double> s_pos;

		ParticleData() {} //empty constructor for making an empty ParticleData
		ParticleData(const std::vector<double>& v_para, const std::vector<double>& v_perp, double mass);
		void free();
	};

	struct DLLCLEXP Maxwellian
	{ //at some point, needs to take in distribution data to generate weights/etc
		double dlogE_dist;            //dlogE of the particle distribution (not bins)
		double ionModFactor{ 1.0 };   //for now, the user needs to adjust these manually
		double magModFactor{ 1.0 };   //eventually, I'd like to calculate them automatically

		std::vector<double> ionEPeak; //ionospheric source peak energy of maxwellian
		std::vector<double> iondEMag; //ionospheric source magnitude of peak

		std::vector<double> magEPeak; //same for magnetosphere
		std::vector<double> magdEMag; //basically x, y axes (top, bottom) of a maxwellian graph

		Maxwellian(double dlogEdist);

		void push_back_ion(double E_peak, double dE_magnitude, int partsAtE = 1);
		void push_back_mag(double E_peak, double dE_magnitude, int partsAtE = 1);
		
		std::vector<double> counts(const ParticleData& init, const double s_ion, const double s_mag);
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

		//postprocessing bins
		const std::vector<double> energyBins;
		const std::vector<double> pitchBins;

		//particle data
		std::vector<double> maxCounts; //maxwellian weights - "number of particles represented by the particle at same index"
		ParticleData initial;          //initial particle data
		ParticleData bottom;           //data on particles that escape out the bottom
		ParticleData upward;           //data on particles that are upgoing at satellite
		ParticleData dnward;           //data on particles that are downgoing at satellite

		PPData(Maxwellian maxwellian, std::vector<double> EBins, std::vector<double> PABins, std::string simDataDir, std::string particleName, std::string simBtmSat, std::string upgoingAtAltSat, std::string dngoingAtAltSat);
	};
} //end namespace postprocess

#endif /* !POSTPROCESS_HELPER_H */