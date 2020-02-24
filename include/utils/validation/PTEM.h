#ifndef VALIDATION_PTEM_H
#define VALIDATION_PTEM_H

#include <functional>

#include "BField/BModel.h"
#include "Simulation/Simulation.h"

using std::function;

namespace validation
{
	bool PTEM_dipB_noE(double maxError, bool printError = false);
	bool PTEM_dipB_QSPS(double maxError, bool printError = false);
	bool PTEM_dist_noE(double maxError, string simDataDir, bool printError = false); //take an existing distribution and validate it
	
	meters particleIdealMirrorAltitude(BModel* bmodel, meters sinit, degrees PAinit);
	//bool PTEM_rungeKutta(double maxError); //validate low error for the RK - "(old + new) / 2" factor I used
}

#endif