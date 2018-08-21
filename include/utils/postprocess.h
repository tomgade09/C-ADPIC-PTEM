#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include "dlldefines.h"
#include "utils/PPHelperClasses.h"

#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

typedef std::vector<std::vector<double>> dblVec2D;

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(PPData ppdata);

	namespace steady
	{
		DLLEXP dblVec2D bsSrcNFluxToSatNFlux(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satUpwardData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies);
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satEFlux(const ParticleData& sat, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& numWt);
		DLLEXP dblVec2D backEFlux(const ParticleData& initialData, const ParticleData& satData, const ParticleData& escapeData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge);
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts);
		DLLEXP void countsToEFlux(dblVec2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection);
		DLLEXP void divBinsByCosPitch(dblVec2D& data, std::vector<double> binAnglesDegrees);
		DLLEXP void symmetricBins0To360(dblVec2D& data, std::vector<double>& binAngles);
	}

	namespace backscat
	{
		DLLEXP double F_flux(double evalE, double incidentE, double incidentCnt, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb);
		DLLEXP double integralF_flux(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb);
		DLLEXP std::vector<double> Nflux(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb);
	}
}

#endif /* !POSTPROCESS_H */