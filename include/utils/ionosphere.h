#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <vector>
#include "dlldefines.h"
#include "utils/ionosphereClasses.h"

#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

namespace ionosphere
{
	DLLEXP dEflux_v2D steadyFlux(const EOMSimData& eomdata);

	namespace dEFlux
	{
		DLLEXP dEflux_v2D satellite(const ParticleData& particles, const Bins& satBins, const dNflux_v1D& dNatSat);
		DLLEXP dEflux_v2D backscatr(const EOMSimData& eomdata, const dNflux_v1D& dNatIonsph);
	}

	namespace binning
	{
		DLLEXP dNflux_v2D binParticles(const ParticleData& particles, const Bins& bins, const dNflux_v1D& countPerParticle);
		DLLEXP void       symmetricBins0To360(dEflux_v2D& data, double_v1D& binAngles);
	}

	namespace backscat
	{
		DLLEXP dNflux     johnd_flux(double E_eval, double E_incident, dEflux dE_incident);
		DLLEXP dEflux     integralJohnd_flux(double lower, double upper, double E_incident);
		DLLEXP dNflux_v2D downwardToBackscatter(const Bins& dist, const dNflux_v2D& escapeCountBinned);
		DLLEXP dNflux_v2D ionsphToSatellite(const EOMSimData& eomdata, const dNflux_v2D& bsCounts);
	}

	namespace multiLevelBS
	{
		DLLEXP dNflux_v2D scatterMain(const EOMSimData& eom, const dNflux_v2D& ionsphTopLvl);
		DLLEXP dNflux_v2D bsAtLevel(const EOMSimData& eom, const dNflux_v2D& ionsphTopLvl, double_v2D& sumCollideAbove, unsigned int level);
		DLLEXP double_v2D alt_reflect(const Bins& distbins, BField* B, double B_ion, double t);
		DLLEXP double     scatterPct (double sumCollideAbove, double Z, double p, double h, eV E, degrees PA);
	}
}

#endif /* !IONOSPHERE_H */