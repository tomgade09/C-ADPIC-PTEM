#include "classAPI.h"

DLLEXP_NOEXTC std::unique_ptr<DipoleB> newDipoleBAPI(double ILATDegrees, double errorTolerance, double ds) {
	return std::make_unique<DipoleB>(ILATDegrees, errorTolerance, ds); }

DLLEXP_NOEXTC std::unique_ptr<DipoleBLUT> newDipoleBLUTAPI(double ILATDegrees, double simMin, double simMax, double ds_dipoleB, int numberOfMeasurements) {
	return std::make_unique<DipoleBLUT>(ILATDegrees, simMin, simMax, ds_dipoleB, numberOfMeasurements); }

DLLEXP_NOEXTC std::unique_ptr<QSPS> newQSPSAPI(std::vector<double> altMinMax, std::vector<double> magnitude) {
	return std::make_unique<QSPS>(altMinMax, magnitude); }