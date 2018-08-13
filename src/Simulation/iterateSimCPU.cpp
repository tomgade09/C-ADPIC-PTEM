//Standard Library includes
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>

//Project specific includes
#include "physicalconstants.h"
#include "utils/loopmacros.h"
#include "Simulation/Simulation.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/SimFatalException.h"

//OpenMP Test
#include <omp.h>
#include <thread>
//OpenMP Test

//forward decls
namespace physics
{
	void vperpMuConvert(const double vpara, double* vperpOrMu, const double s, const double t_convert, BField* bfield, const double mass, const bool vperpToMu);
	void iterateParticle(double* vpara, double* mu, double* s, double* t_incident, double* t_escape, BField* bfield, EField* efield, const double simtime, const double dt, const double mass, const double charge, const double simmin, const double simmax);
}

void Simulation::__iterateSimCPU(int numberOfIterations, int checkDoneEvery)
{
	printSimAttributes(numberOfIterations, checkDoneEvery);

	using namespace physics;
	for (auto part = particles_m.begin(); part < particles_m.end(); part++)
	{
		(*part)->loadDataFromMem((*part)->data(true), false); //copies orig data to curr - normally this happens from CPU->GPU->CPU, but we aren't using the GPU here

		std::vector<std::vector<double>>& data{ (*part)->__data(false) };
		for (int ind = 0; ind < (*part)->getNumberOfParticles(); ind++)
		{//convert vperp to mu in Particle memory
			vperpMuConvert(data.at(0).at(ind), &data.at(1).at(ind), data.at(2).at(ind), data.at(4).at(ind), BFieldModel_m.get(), (*part)->mass(), true);
		}
	}

	bool done{ false };
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{
		//if (cudaloopind % checkDoneEvery == 0) { done = true; }
		for (auto part = particles_m.begin(); part < particles_m.end(); part++)
		{
			std::vector<std::vector<double>>& data{ (*part)->__data(false) }; //get a reference to the particle's curr data array

			//OpenMP Test
			omp_set_num_threads(std::thread::hardware_concurrency());
			
			#pragma omp parallel for
			for (int ind = 0; ind < (*part)->getNumberOfParticles(); ind++)
			{
				iterateParticle(&data.at(0).at(ind), &data.at(1).at(ind), &data.at(2).at(ind), &data.at(3).at(ind), &data.at(4).at(ind),
					BFieldModel_m.get(), EFieldModel_m.get(), simTime_m, dt_m, (*part)->mass(), (*part)->charge(), simMin_m, simMax_m);
				if ((cudaloopind % checkDoneEvery == 0) && done && (data.at(4).at(ind) < 0.0))
				{
					//#pragma omp atomic
					done = false; //maybe a problem - will have at least several simultaneous attempts to write...not thread-safe, but it's simply a flag, so here it doesn't need to be maybe?
					//OpenMP 2.0 (max supported version by VS, but very old) doesn't support done = false; as a legit expression following #pragma omp atomic
				}
			}
			//OpenMP Test
		}

		//iterate Satellites here (need a host-side detector)

		if (cudaloopind % checkDoneEvery == 0)
		{
			std::stringstream out;
			out << std::setw(std::to_string(numberOfIterations).size()) << cudaloopind;
			std::cout << out.str() << " / " << numberOfIterations << "  |  Sim Time (s): ";
			out.str(""); out.clear();
			out << std::setw(std::to_string((int)(numberOfIterations)* dt_m).size()) << std::fixed << simTime_m;
			std::cout << out.str() << "  |  Real Time Elapsed (s): ";
			logFile_m->printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << std::endl;

			if (done) { std::cout << "All particles finished early.  Breaking loop." << std::endl; break; }
		}

		incTime();
		cudaloopind++;
	}

	for (auto part = particles_m.begin(); part < particles_m.end(); part++)
	{
		std::vector<std::vector<double>> tmp{ (*part)->data(false) };
		for (int ind = 0; ind < (*part)->getNumberOfParticles(); ind++)
		{//convert mu to vperp in Particle memory
			vperpMuConvert(tmp.at(0).at(ind), &tmp.at(1).at(ind), tmp.at(2).at(ind), tmp.at(4).at(ind), BFieldModel_m.get(), (*part)->mass(), false);
		}
		(*part)->loadDataFromMem(tmp, false);
	}

	saveReady_m = true;
	saveDataToDisk();
	simTime_m = 0.0;

	std::cout << "Total sim time: "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;

	logFile_m->createTimeStruct("End Iterate " + std::to_string(numberOfIterations));
	logFile_m->writeLogFileEntry("Simulation::iterateSimulation: End Iteration of Sim:  " + std::to_string(numberOfIterations));
}