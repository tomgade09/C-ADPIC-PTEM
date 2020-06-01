#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <vector>
#include "dlldefines.h"
#include "ionosphere/ionosphereClasses.h"

namespace ionosphere
{
	DLLEXP vector<vector<dEflux>> steadyFlux(const EOMSimData& eomdata);

	namespace dEFlux
	{
		DLLEXP ParticlesBinned<dEflux> satellite(const EOMSimData& eomdata);
		DLLEXP ParticlesBinned<dEflux> backscatr(const EOMSimData& eomdata);
		DLLEXP degrees newPA(const degrees PA_init, const tesla B_init, const tesla B_final);
	}

	namespace backscat
	{
		DLLEXP dNflux johnd_flux(eV E_eval, eV E_incident, dNflux dN_incident);
		DLLEXP ParticlesBinned<dNflux> downwardToBackscatter(const Bins& dist, const ParticlesBinned<dNflux>& dNpointofScatter);
		DLLEXP ParticlesBinned<dEflux> ionsphToSatellite(const EOMSimData& eomdata, const ParticlesBinned<dNflux>& bsCounts);
	}

	namespace multiLevelBS
	{
		DLLEXP ParticlesBinned<dNflux> scatterMain(const EOMSimData& eom, const ParticlesBinned<dNflux>& dNionsphTop);
		DLLEXP ParticlesBinned<dNflux> upgAtLevel(const EOMSimData& eom, const ParticlesBinned<dNflux>& dNionsphTop, double_v2D& pctScatteredAbove, size_t level);
		DLLEXP percent scatterPct(percent sumCollideAbove, double Z, percm3 p, cm h, eV E, degrees PA);
	}
}

#endif /* !IONOSPHERE_H */