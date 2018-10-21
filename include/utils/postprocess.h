#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include "dlldefines.h"
#include "utils/PPHelperClasses.h"

#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

using std::vector;
typedef vector<vector<double>> dblVec2D;
typedef vector<double> dblVec;

namespace postprocess
{
	DLLEXP dblVec2D steadyFlux(const PPData& ppdata);

	namespace steady
	{
		DLLEXP dblVec2D bsSrcToSat(const Bins& dist, const Bins& sat, const dblVec2D& bsCounts, const ParticleData& initialData, const ParticleData& satUpwardData);
	}

	namespace EFlux
	{
		DLLEXP dblVec2D satdEFlux(const ParticleData& sat, const Bins& satBins, const dblVec& weights_atSat);
		DLLEXP dblVec2D bksdEFlux(const PPData& ppdata, const dblVec& numWeight);
	}

	namespace binning
	{
		DLLEXP dblVec2D binWeighted(const ParticleData& particles, const Bins& bins, const dblVec& weights_atIon);
		DLLEXP void     symmetricBins0To360(dblVec2D& data, dblVec& binAngles);
	}

	namespace backscat
	{
		DLLEXP double   evans_flux(double E_eval, double E_incident);
		DLLEXP double   integralEvans_flux(double lower, double upper, double incidentE);
		DLLEXP dblVec2D dNflux_bs_ion(const Bins& dist, const dblVec2D& escapeCountBinned);
	}

	namespace multLevelBS
	{
		//DLLEXP dblVec2D scatterMain(const Ionosphere& ionsph, const Bins& escape);
		//DLLEXP dblVec2D bsAtLevel(dblVec2D counts, const dblVec& PAbins, const dblVec& Ebins, double* sumCollideAbove, double Z, double p, double h, double E, double PA);
		//DLLEXP dblVec2D alt_reflect(const Bins& bins, BField* B, double B_ion, double t);
		//DLLEXP double   scatterPct(double sumCollideAbove, double Z, double p, double h, double E, double PA);
	}
}

#endif /* !POSTPROCESS_H */