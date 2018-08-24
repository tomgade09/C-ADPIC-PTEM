#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include "dlldefines.h"
#include "utils/PPHelperClasses.h"

//#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

using std::vector;
typedef vector<vector<double>> dblVec2D;

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(const PPData& ppdata);

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const dblVec2D& bsNumFluxBins, const ParticleData& initialData, const ParticleData& satUpwardData, const vector<double>& binAngles, const vector<double>& binEnergies);
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satdEFlux(const ParticleData& sat, const vector<double>& binAngles, const vector<double>& binEnergies, const vector<double>& numWeight);
		DLLEXP dblVec2D bksdEFlux(const ParticleData& initialData, const ParticleData& satData, const ParticleData& escapeData, const vector<double>& ppPABins, const vector<double>& ppEBins, const vector<double>& maxwCounts, const vector<double>& distPAbins, const vector<double>& distEbins);
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const vector<double>& particlePitches, const vector<double>& particleEnergies, const vector<double>& binAngles, const vector<double>& binEnergies, const vector<double>& counts);
		DLLEXP void symmetricBins0To360(dblVec2D& data, vector<double>& binAngles);
	}

	namespace backscat
	{
		DLLEXP double evans_flux(double E_eval, double E_incident, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb);
		DLLEXP double integralEvans_flux(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb);
		DLLEXP vector<double> dNflux_bs(const vector<double>& binCounts, const vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb);
	}
}

#endif /* !POSTPROCESS_H */