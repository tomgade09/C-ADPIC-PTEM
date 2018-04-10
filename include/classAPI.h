#ifndef CLASSAPI_H
#define CLASSAPI_H

#include <memory> //smart pointers
#include "dllexport.h"
#include "SimulationClass\Simulation.h"
#include "BField\allBModels.h"
#include "EField\allEModels.h"

namespace API
{
	DLLEXP_NOEXTC std::unique_ptr<Simulation> loadPreviousSimulation(std::string prevSimDir);
	DLLEXP_NOEXTC std::unique_ptr<DipoleB> newDipoleB(double ILATDegrees, double errorTolerance, double ds);
	DLLEXP_NOEXTC std::unique_ptr<DipoleBLUT> newDipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_dipoleB, int numberOfMeasurements);
	DLLEXP_NOEXTC std::unique_ptr<QSPS> newQSPS(std::vector<double> altMinMax, std::vector<double> magnitude);
}

#endif /* !CLASSAPI_H */