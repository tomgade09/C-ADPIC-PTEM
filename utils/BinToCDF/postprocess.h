#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#define PPCOUTDBG
//#define PPCOUTDBG_VERBOSE

#include <vector>
#include "FileIO\fileIO.h"

#define TRYCATCHSTDEXP(x) try{ x; }catch(std::exception& e){ std::cout << e.what() << " -> exiting." <<  std::endl; exit(1); }

typedef std::vector<std::vector<double>> vecDbl2D;

namespace postprocess
{
	vecDbl2D steadyFlux(std::string datarootfolder, std::vector<double> pitchBin_Min_Max, int pitchBinNum, std::vector<double> logEBin_Min_Max, int logEBinNum,
		std::vector<double> ionMaxwellian_kT_dEflux, std::vector<double> magMaxwellian_kT_dEflux, double mass, int particlecount, std::vector<double>& binAnglesOut, std::vector<double>& binEnergiesOut);
	//vecDbl2D timedepFlux();

	namespace steady
	{
		vecDbl2D simEnergyFlux(const vecDbl2D& particleData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge, double BatXSection);
		vecDbl2D bsEnergyFlux(const vecDbl2D& initialData, const vecDbl2D& satData, const vecDbl2D& escapeData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts, double mass, double charge, double BatXSection);
	}

	namespace timedep
	{
		vecDbl2D simEnergyFlux();
		vecDbl2D bsEnergyFlux();
	}

	namespace numerical
	{
		void vToEPitch(const std::vector<double>& vpara, const std::vector<double>& vperp, double mass, std::vector<double>& particlePitches, std::vector<double>& particleEnergies);
		std::vector<double> generatePitchBins(double pitchMin, double pitchMax, int numBins);
		std::vector<double> generateLogSpacedEnergyBins(double logEmidBinMin, double logEmidBinMax, int numBins);
		void splitIonMagEngs(const std::vector<double>& s_init, const std::vector<double>& E_init, std::vector<double>& ionsphE, std::vector<double>& magsphE);
		std::vector<double> maxwellianCounts(const std::vector<double>& initEnergies, double sigma_kT, double dEflux_kT);
		vecDbl2D countInBinsWeighted(const std::vector<double>& particlePitches, const std::vector<double>& particleEnergies, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, const std::vector<double>& maxwCounts);
		void countsToEFlux(vecDbl2D& energyData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies, double mass, double charge, double BatXSection);
		void divBinsByCosPitch(vecDbl2D& data, std::vector<double> binAnglesDegrees);
	}

	namespace backscat
	{
		double F_flux(double evalE, double incidentE, double incidentCnt, double prim_logm, double prim_logb, double scnd_logm, double scnd_logb);
		double integralF_flux(double lower, double upper, double incidentE, double prim_fact, double prim_logb, double scnd_fact, double scnd_logb);
		std::vector<double> sumIntegralsOfNumFluxFcnsPerBin(const std::vector<double>& binCounts, const std::vector<double>& binEnergies, double primary_logm, double primary_logb, double secondary_logm, double secondary_logb);
		vecDbl2D matchIonBSToSatAndCount(const std::vector<double>& bsEFluxBins, const vecDbl2D& initialData, const vecDbl2D& satDownData, const std::vector<double>& binAngles, const std::vector<double>& binEnergies);
	}

	namespace utils
	{
		void loadvFromDisk(std::string savefolder, std::string filePrefix, int count, std::vector<double>& vpara, std::vector<double>& vperp);
	}
}

#endif /* !POSTPROCESS_H */