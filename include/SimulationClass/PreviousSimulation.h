#ifndef SIMULATION_PREVIOUS_H
#define SIMULATION_PREVIOUS_H

#include "SimulationClass\Simulation.h"

class PreviousSimulation : public Simulation
{
protected:
	std::string prevDataDir_m{ "" };

public:
	PreviousSimulation(std::string prevDataDir)
	{
		//load the previous data - depends on the previous sim saving data...maybe set a function returning a vector of characteristics
	}

};

#endif