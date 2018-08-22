#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include "dlldefines.h"
#include "utils/PPHelperClasses.h"

//#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

typedef std::vector<std::vector<double>> dblVec2D;

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(const PPData& ppdata);

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satUpwardData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies);
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satdEFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWeight);
		DLLEXP dblVec2D bksdEFlux(const ParticleData& initial, const ParticleData& sat, const ParticleData& escape, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWeight);
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& counts);
		DLLEXP void symmetricBins0To360(dblVec2D& data, std::vector<double>& binAngles);
	}

	namespace backscat
	{
		DLLEXP double evans_flux(double E_eval, double E_incident, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb);
		DLLEXP double integralEvans_flux(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb);
		DLLEXP std::vector<double> dNflux_bs(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb);
	}
}

#endif /* !POSTPROCESS_H */