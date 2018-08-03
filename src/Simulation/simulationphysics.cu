//Standard Library includes
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"

//Project specific includes
#include "physicalconstants.h"
#include "utils/loopmacros.h"
#include "Simulation/Simulation.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/SimFatalException.h"

//CUDA Variables - if you change these, don't forget to change the associated curand code/blocks/etc
// For Geforce 960M (author's computer) - maximum 1024 threads per block - try this to see if it results in faster code execution sometime
constexpr int  BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance

//Commonly used values
extern const int SIMCHARSIZE{ 3 * sizeof(double) };

namespace physics
{
	__global__ void vperpMuConvert(double** dataToConvert, BField** bfield, double mass, bool vperpToMu, int timeInd = 4)
	{//dataToConvert[0] = vpara, [1] = vperp, [2] = s, [3] = t_incident, [4] = t_escape
		unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };

		if (dataToConvert[1][thdInd] != 0.0)
		{
			double B_s{ (*bfield)->getBFieldAtS(dataToConvert[2][thdInd], dataToConvert[timeInd][thdInd]) };
			if (vperpToMu)
				dataToConvert[1][thdInd] = 0.5 * mass * dataToConvert[1][thdInd] * dataToConvert[1][thdInd] / B_s;
			else
				dataToConvert[1][thdInd] = sqrt(2 * dataToConvert[1][thdInd] * B_s / mass);
		}
	}

	__host__ void vperpMuConvert(const double vpara, double* vperpOrMu, const double s, const double t_convert, BField* bfield, const double mass, const bool vperpToMu)
	{//dataToConvert[0] = vpara, [1] = vperp, [2] = s, [3] = t_incident, [4] = t_escape
		if (*vperpOrMu != 0.0)
		{
			double B_s{ bfield->getBFieldAtS(s, t_convert) };
			if (vperpToMu)
				*vperpOrMu = 0.5 * mass * (*vperpOrMu) * (*vperpOrMu) / B_s;
			else
				*vperpOrMu = sqrt(2 * (*vperpOrMu) * B_s / mass);
		}
	}

	__device__ __host__ double accel1dCUDA(const double vs_RK, const double t_RK, const double* args, BField** bfield, EField** efield) //made to pass into 1D Fourth Order Runge Kutta code
	{//args array: [s_0, mu, q, m, simtime]
		double F_lor, F_mir, stmp;
		stmp = args[0] + vs_RK * t_RK; //ps_0 + vs_RK * t_RK

		//Mirror force
		F_mir = -args[1] * (*bfield)->getGradBAtS(stmp, t_RK + args[4]); //-mu * gradB(pos, runge-kutta time + simtime)

		//Lorentz force - simply qE - v x B is taken care of by mu - results in kg.m/s^2 - to convert to Re equivalent - divide by Re
		F_lor = args[2] * (*efield)->getEFieldAtS(stmp, t_RK + args[4]); //q * EFieldatS

		return (F_lor + F_mir) / args[3];
	}//returns an acceleration in the parallel direction to the B Field

	__device__ __host__ double foRungeKuttaCUDA(const double y_0, const double h, const double* funcArg, BField** bfield, EField** efield)
	{
		// dy / dt = f(t, y), y(t_0) = y_0
		// funcArgs are whatever you need to pass to the equation
		// args array: [s_0, mu, q, m, simtime]
		double k1, k2, k3, k4; double y{ y_0 }; double t_RK{ 0.0 };
		
		k1 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k1 = f(t_n, y_n), returns units of dy / dt

		t_RK = h / 2;
		y = y_0 + k1 * t_RK;
		k2 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k2 = f(t_n + h/2, y_n + h/2 * k1)

		y = y_0 + k2 * t_RK;
		k3 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k3 = f(t_n + h/2, y_n + h/2 * k2)

		t_RK = h;
		y = y_0 + k3 * t_RK;
		k4 = accel1dCUDA(y, t_RK, funcArg, bfield, efield); //k4 = f(t_n + h, y_n + h k3)

		return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns delta y, not dy / dt, not total y
	}

	__global__ void simActiveCheck(double** currData_d, bool* simDone)
	{
		//Answers the question: Are there no particles left in the simulation?
		//stores the value to simDone, which is defaulted to true, and flipped to false
		//only if t_escape is less than zero for at least one particle
		//(in that case, the sim is not completely done iterating)
		if (*simDone)
		{
			const double* t_escape_d{ currData_d[4] }; //const double* t_incident_d{ currData_d[3] }; //to be implemented

			unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };

			if (t_escape_d[thdInd] >= 0.0) //particle has escaped the sim
				return;
			else
				(*simDone) = false;
		}
	}

	__global__ void iterateParticle(double** currData_d, BField** bfield, EField** efield,
		const double simtime, const double dt, const double mass, const double charge, const double simmin, const double simmax)
	{
		unsigned int thdInd{ blockIdx.x * blockDim.x + threadIdx.x };

		double* v_d{ currData_d[0] }; const double* mu_d{ currData_d[1] }; double* s_d{ currData_d[2] }; const double* t_incident_d{ currData_d[3] }; double* t_escape_d{ currData_d[4] };

		if (t_escape_d[thdInd] >= 0.0) //particle has escaped, t_escape is >= 0 iff it has both entered and is outside the sim boundaries
			return;
		else if (t_incident_d[thdInd] > simtime) //particle hasn't "entered the sim" yet
			return;
		else if (s_d[thdInd] < simmin * 0.999) //particle is out of sim to the bottom and t_escape not set yet
		{//eventually build in "fuzzy boundary" - maybe eventually create new particle with initial characteristics on escape
			t_escape_d[thdInd] = simtime;
			return;                      //fuzzyIonosphere(); if (t_escape_d[thdInd] >= 0.0 && t_escape_d[thdInd] < simtime) { return; }
		}
		else if (s_d[thdInd] > simmax * 1.001) //particle is out of sim to the top and t_escape not set yet
		{//maybe eventaully create new particle with initial characteristics on escape
			t_escape_d[thdInd] = simtime;
			return;
		}

		//args array: [ps_0, mu, q, m, simtime]
		const double args[]{ s_d[thdInd], mu_d[thdInd], charge, mass, simtime };

		//foRK (plus v0 in this case) gives v at the next time step (indicated vf in this note):
		//for downgoing (as an example), due to the mirror force, vf will be lower than v0 as the mirror force is acting in the opposite direction as v
		//along the path of the particle, ds, and so s will end up higher if we use ds = (vf * dt) than where it would realistically
		//if we use the ds = (v0 * dt), s will be lower down than where it would end up really (due to the fact that the mirror force acting along ds
		//will slow v down as the particle travels along ds), so I take the average of the two and it seems close enough s = (v0 + (v0 + dv)) / 2 * dt = v0 + dv/2 * dt
		//hence the /2 factor below - FYI, this was checked by the particle's energy (steady state, no E Field) remaining the same throughout the simulation
		double v_orig{ v_d[thdInd] };
		v_d[thdInd] += foRungeKuttaCUDA(v_d[thdInd], dt, args, bfield, efield);
		s_d[thdInd] += (v_d[thdInd] + v_orig) / 2 * dt;
	}

	__host__ void iterateParticle(double* vpara, double* mu, double* s, double* t_incident, double* t_escape, BField* bfield, EField* efield,
		const double simtime, const double dt, const double mass, const double charge, const double simmin, const double simmax)
	{
		if (simtime == 0.0) { *t_escape = -1.0; }
		if (*t_escape >= 0.0) //see above function for description of conditions
			return;
		else if (*t_incident > simtime)
			return;
		else if (*s < simmin * 0.999)
		{
			*t_escape = simtime;
			return;
		}
		else if (*s > simmax * 1.001)
		{
			*t_escape = simtime;
			return;
		}

		const double args[]{ *s, *mu, charge, mass, simtime };

		double v_orig{ *vpara };
		*vpara += foRungeKuttaCUDA(*vpara, dt, args, &bfield, &efield);
		*s += (*vpara + v_orig) / 2 * dt;
	}

	/*
	__device__ void fuzzyIonosphere(double& s_d, const double s_esc_absolute, double& v_d, double& t_escape_d, const double simtime)
	{
		if (v_d > 0.0) { return; } //or do we want upgoing to possibly collide???
		if (someRandomGenerator >/<(=)/== someCondition || s_d <= s_esc_absolute)
			t_escape_d = simtime;
		t_escape_d = simtime;
	}
	*/
}

//Simulation member functions
void Simulation::initializeSimulation()
{
	if (BFieldModel_m == nullptr)
		throw SimFatalException("Simulation::initializeSimulation: no Magnetic Field model specified", __FILE__, __LINE__);
	if (particles_m.size() == 0)
		throw SimFatalException("Simulation::initializeSimulation: no particles in simulation, sim cannot be initialized without particles", __FILE__, __LINE__);

	if (EFieldModel_m == nullptr) //make sure an EField (even if empty) exists
		EFieldModel_m = std::make_unique<EField>();
	
	EFieldModel_d = EFieldModel_m->getPtrGPU();

	if (tempSats_m.size() > 0)
	{ LOOP_OVER_1D_ARRAY(tempSats_m.size(), createSatellite(tempSats_m.at(iii).get())); } //create satellites
	else
		std::cerr << "Simulation::initializeSimulation: warning: no satellites created" << std::endl;

	initialized_m = true;
}

void Simulation::__iterateSimCPU(int numberOfIterations, int checkDoneEvery)
{
	using namespace physics;
	for (auto part = particles_m.begin(); part < particles_m.end(); part++)
	{
		std::vector<std::vector<double>> tmp{ (*part)->data(false) };
		for (int ind = 0; ind < (*part)->getNumberOfParticles(); ind++)
		{//convert vperp to mu in Particle memory
			vperpMuConvert(tmp.at(0).at(ind), &tmp.at(1).at(ind), tmp.at(2).at(ind), tmp.at(4).at(ind), BFieldModel_m.get(), (*part)->mass(), true);
		}
		(*part)->loadDataFromMem(tmp, false);
	}

	std::cout << "\tvpara\tvperp\ts\tE\n";

	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{
		bool done{ true };
		for (auto part = particles_m.begin(); part < particles_m.end(); part++)
		{
			std::vector<std::vector<double>> tmp{ (*part)->data(false) };
			for (int ind = 0; ind < (*part)->getNumberOfParticles(); ind++)
			{
				iterateParticle(&tmp.at(0).at(ind), &tmp.at(1).at(ind), &tmp.at(2).at(ind), &tmp.at(3).at(ind), &tmp.at(4).at(ind),
					BFieldModel_m.get(), EFieldModel_m.get(), simTime_m, dt_m, (*part)->mass(), (*part)->charge(), simMin_m, simMax_m);
				if ((cudaloopind % checkDoneEvery == 0) && done && (tmp.at(4).at(ind) < 0.0))
					done = false;
			}
			double vperp{ tmp.at(1).at(0) };
			vperpMuConvert(tmp.at(0).at(0), &vperp, tmp.at(2).at(0), tmp.at(4).at(0), BFieldModel_m.get(), (*part)->mass(), false);
			if (cudaloopind % checkDoneEvery == 0) { std::cout << tmp.at(0).at(0) << "  " << vperp << "  " << tmp.at(2).at(0) << "  " << tmp.at(1).at(0) << "  "; std::cout << std::setprecision(10) << 0.5 * (*part)->mass() * (tmp.at(0).at(0) * tmp.at(0).at(0) + vperp * vperp) / JOULE_PER_EV << std::setprecision(6) << "\n"; }
			(*part)->loadDataFromMem(tmp, false);
		}

		if (cudaloopind % checkDoneEvery == 0)
		{
			//what else here??
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
}

void Simulation::iterateSimulation(int numberOfIterations, int checkDoneEvery)
{//conducts iterative calculations of data previously copied to GPU - runs the data through the computeKernel
	using namespace physics;
	
	if (!initialized_m)
		throw SimFatalException("Simulation::iterateSimulation: sim not initialized with initializeSimulation()", __FILE__, __LINE__);
	
	printSimAttributes(numberOfIterations, checkDoneEvery);
	
	logFile_m->createTimeStruct("Start Iterate " + std::to_string(numberOfIterations));
	logFile_m->writeLogFileEntry("Simulation::iterateSimulation: Start Iteration of Sim:  " + std::to_string(numberOfIterations));

	//Copy data to device
	LOOP_OVER_1D_ARRAY(particles_m.size(), particles_m.at(iii)->copyDataToGPU());
	
	//convert particle vperp data to mu
	for (auto part = particles_m.begin(); part < particles_m.end(); part++)
		vperpMuConvert <<< (*part)->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> ((*part)->getCurrDataGPUPtr(), BFieldModel_d, (*part)->mass(), true);
	
	//Setup on GPU variable that checks to see if any threads still have a particle in sim and if not, end iterations early
	bool* simDone_d{ nullptr };
	CUDA_API_ERRCHK(cudaMalloc((void**)&simDone_d, sizeof(bool)));

	//Loop code
	long cudaloopind{ 0 };
	while (cudaloopind < numberOfIterations)
	{	
		if (cudaloopind % checkDoneEvery == 0) { CUDA_API_ERRCHK(cudaMemset(simDone_d, true, sizeof(bool))); } //if it's going to be checked in tnis iter (every checkDoneEvery iterations), set to true

		for (auto part = particles_m.begin(); part < particles_m.end(); part++)
		{
			iterateParticle <<< (*part)->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> ((*part)->getCurrDataGPUPtr(), BFieldModel_d, EFieldModel_d,
				simTime_m, dt_m, (*part)->mass(), (*part)->charge(), simMin_m, simMax_m);
			
			//kernel will set boolean to false if at least one particle is still in sim
			if (cudaloopind % checkDoneEvery == 0)
				simActiveCheck <<< (*part)->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> ((*part)->getCurrDataGPUPtr(), simDone_d);
		}

		CUDA_KERNEL_ERRCHK_WSYNC_WABORT(); //side effect: cudaDeviceSynchronize() needed for computeKernel to function properly, which this macro provides

		for (auto sat = satPartPairs_m.begin(); sat < satPartPairs_m.end(); sat++)
			(*sat)->satellite->iterateDetector(simTime_m, dt_m, BLOCKSIZE);
		
		if (cudaloopind % checkDoneEvery == 0)
		{
			std::stringstream out;
			out << std::setw(std::to_string(numberOfIterations).size()) << cudaloopind;
			std::cout << out.str() << " / " << numberOfIterations << "  |  Sim Time (s): ";
			out.str(""); out.clear();
			out << std::setw(std::to_string((int)(numberOfIterations) * dt_m).size()) << std::fixed << simTime_m;
			std::cout << out.str() << "  |  Real Time Elapsed (s): ";
			logFile_m->printTimeNowFromLastTS(); //need to add to log file as well?
			std::cout << std::endl;

			bool done{ false };
			CUDA_API_ERRCHK(cudaMemcpy(&done, simDone_d, sizeof(bool), cudaMemcpyDeviceToHost));
			if (done) { std::cout << "All particles finished early.  Breaking loop." << std::endl; break; }
		}

		incTime();
		cudaloopind++;
	}

	CUDA_API_ERRCHK(cudaFree(simDone_d));

	//Convert particle, satellite mu data to vperp
	for (auto part = particles_m.begin(); part < particles_m.end(); part++)
		vperpMuConvert <<< (*part)->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>> ((*part)->getCurrDataGPUPtr(), BFieldModel_d, (*part)->mass(), false); //nullptr will need to be changed if B ever becomes time dependent, would require loop to record when it stops tracking the particle

	for (auto sat = satPartPairs_m.begin(); sat < satPartPairs_m.end(); sat++)
		vperpMuConvert <<< (*sat)->particle->getNumberOfParticles() / BLOCKSIZE, BLOCKSIZE >>>  ((*sat)->satellite->get2DDataGPUPtr(), BFieldModel_d, (*sat)->particle->mass(), false, 3);

	//Copy data back to host
	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->copyDataToHost());
	LOOP_OVER_1D_ARRAY(getNumberOfSatellites(), satellite(iii)->copyDataToHost());

	saveReady_m = true;
	saveDataToDisk();
	simTime_m = 0.0;

	std::cout << "Total sim time: "; logFile_m->printTimeNowFromFirstTS(); std::cout << " s" << std::endl;

	logFile_m->createTimeStruct("End Iterate " + std::to_string(numberOfIterations));
	logFile_m->writeLogFileEntry("Simulation::iterateSimulation: End Iteration of Sim:  " + std::to_string(numberOfIterations));
}

void Simulation::freeGPUMemory()
{//used to free the memory on the GPU that's no longer needed
	if (!initialized_m)
		throw SimFatalException("Simulation::freeGPUMemory: simulation not initialized with initializeSimulation()", __FILE__, __LINE__);

	if (!dataOnGPU_m) { return; }

	logFile_m->writeLogFileEntry("Simulation::freeGPUMemory: Start free GPU Memory.");

	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->freeGPUMemory());
	LOOP_OVER_1D_ARRAY(getNumberOfSatellites(), satellite(iii)->freeGPUMemory());

	dataOnGPU_m = false;
	logFile_m->writeLogFileEntry("Simulation::freeGPUMemory: End free GPU Memory.");

	CUDA_API_ERRCHK(cudaProfilerStop()); //For profiling with the CUDA bundle
}
