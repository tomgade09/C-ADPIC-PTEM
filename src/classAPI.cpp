#include "classAPI.h"

namespace API
{
	DLLEXP_NOEXTC std::unique_ptr<Simulation> loadPreviousSimulation(std::string prevSimDir){
		return std::make_unique<Simulation>(prevSimDir); }

	DLLEXP_NOEXTC std::unique_ptr<DipoleB> newDipoleB(double ILATDegrees, double errorTolerance, double ds) {
		return std::make_unique<DipoleB>(ILATDegrees, errorTolerance, ds); }

	DLLEXP_NOEXTC std::unique_ptr<DipoleBLUT> newDipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_dipoleB, int numberOfMeasurements) {
		return std::make_unique<DipoleBLUT>(ILATDegrees, simMin, simMax, ds_dipoleB, numberOfMeasurements); }

	DLLEXP_NOEXTC std::unique_ptr<QSPS> newQSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude) {
		return std::make_unique<QSPS>(altMin, altMax, magnitude); }
}