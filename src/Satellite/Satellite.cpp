#include <omp.h>
#include <thread>

#include "Satellite/Satellite.h"
#include "ErrorHandling/simExceptionMacros.h"

void Satellite::satelliteDetectorCPU(const std::vector<std::vector<double>>& partdata, double simtime, double dt)
{
	std::vector<std::vector<double>>& detected{ data_m.at(0) };

	#pragma omp parallel for
	for (unsigned int partind = 0; partind < partdata.at(0).size(); partind++)
	{
		if (simtime == 0.0)
		{
			detected.at(3).at(partind) = -1.0; //t_esc
			detected.at(4).at(partind) = -1.0; //index
		}

		if (detected.at(3).at(partind) > -0.1)
			continue;

		double s_minus_vdt{ partdata.at(2).at(partind) - partdata.at(0).at(partind) * dt };

		if (
			//(detected.at(3).at(partind) < -0.1) &&
			((!upwardFacing_m) && (partdata.at(2).at(partind) >= altitude_m) && (s_minus_vdt < altitude_m))
			||
			((upwardFacing_m) && (partdata.at(2).at(partind) <= altitude_m) && (s_minus_vdt > altitude_m))
			)
		{
			detected.at(0).at(partind) = partdata.at(0).at(partind); //vpara
			detected.at(1).at(partind) = partdata.at(1).at(partind); //mu
			detected.at(2).at(partind) = partdata.at(2).at(partind); //s
			detected.at(3).at(partind) = simtime;                    //t_esc
			detected.at(4).at(partind) = (double)partind;            //index
		}
	}
}