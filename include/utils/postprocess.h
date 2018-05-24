#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include "dlldefines.h"

#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

typedef std::vector<std::vector<double>> dblVec2D;

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

	struct DLLCLEXP MaxwellianData
	{
		double s_ion;
		double s_mag;
		double dlogE_dist;

		std::vector<double> ionEPeak; //ionospheric source peak energy of maxwellian
		std::vector<double> iondEMag; //ionospheric source magnitude of peak

		std::vector<double> magEPeak; //same for magnetosphere
		std::vector<double> magdEMag; //basically x, y axes (top, bottom) of a maxwellian graph
		
		MaxwellianData(double sion, double smag, double dlogE);

		void push_back_ion(double peak, double magnitude);
		void push_back_mag(double peak, double magnitude);
	};

	struct DLLCLEXP PPData
	{
		double s_ion;
		double s_mag;

		double B_ion; //B Field strength at ionospheric source
		double B_sat; //B Field strength at satellite
		double B_mag; //B Field strength at magnetospheric source

		//postprocessing bins
		const std::vector<double> energyBins;
		const std::vector<double> pitchBins;

		//particle data
		std::vector<double> maxCounts; //maxwellian weights - "number of particles represented by the particle at same index"
		ParticleData initial;          //initial particle data
		ParticleData bottom;           //data on particles that escape out the bottom
		ParticleData upward;           //data on particles that are upgoing at satellite
		ParticleData dnward;           //data on particles that are downgoing at satellite

		PPData(double sion, double smag, double Bion, double Bsat, double Bmag,
			std::vector<double> EBins, std::vector<double> PABins, MaxwellianData maxData, ParticleData init, ParticleData btm, ParticleData up, ParticleData dn);
	};

	DLLEXP dblVec2D steadyFlux(PPData ppdata);

	namespace steady
	{
		DLLEXP dblVec2D matchIonBSToSatAndCount(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satDownData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies);
	}

	namespace EFlux
	{
		DLLEXP dblVec2D simEnergyFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWt);
		DLLEXP dblVec2D bsEnergyFlux(const ParticleData& initialData, const ParticleData& satData, const ParticleData& escapeData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge);
	}

	namespace maxwellian
	{
		DLLEXP double counts(double E, double binWidth_E, double sigma_kT, double dEflux_kT, double binWidth_kT);
		DLLEXP std::vector<double> formCountsVector(const ParticleData& init, const ParticleData& dnward, const MaxwellianData& maxData);
	}

	namespace binning
	{
		DLLEXP dblVec2D countInBinsWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts);
		DLLEXP void countsToEFlux(dblVec2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection);
		DLLEXP void divBinsByCosPitch(dblVec2D& data, std::vector<double> binAnglesDegrees);
	}

	namespace backscat
	{
		DLLEXP double F_flux(double evalE, double incidentE, double incidentCnt, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb);
		DLLEXP double integralF_flux(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb);
		DLLEXP std::vector<double> Nflux(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb);
	}
}

#endif /* !POSTPROCESS_H */