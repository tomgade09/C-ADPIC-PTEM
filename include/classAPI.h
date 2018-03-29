#ifndef CLASSAPI_H
#define CLASSAPI_H

#include <memory> //smart pointers
#include "dllexport.h"
#include "BField\allBModels.h"
#include "EField\allEModels.h"

DLLEXP_NOEXTC std::unique_ptr<DipoleB> newDipoleBAPI(double ILATDegrees, double errorTolerance, double ds);
DLLEXP_NOEXTC std::unique_ptr<DipoleBLUT> newDipoleBLUTAPI(double ILATDegrees, double simMin, double simMax, double ds_dipoleB, int numberOfMeasurements);
DLLEXP_NOEXTC std::unique_ptr<QSPS> newQSPSAPI(std::vector<double> altMinMax, std::vector<double> magnitude);

#endif /* !CLASSAPI_H */