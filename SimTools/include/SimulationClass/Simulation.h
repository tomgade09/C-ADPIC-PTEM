#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <vector>
#include <chrono>
#include "ParticleClass\Particle.h"
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "StandaloneTools\StandaloneTools.h"
#include "include\_simulationvariables.h" //remove this later, replace with passed in variables

struct SatAndPart
{
	Satellite* satellite;
	Particle*  particle;
};

class Simulation
{
protected:
	//Simulation Characteristics
	std::string  rootdir_m;
	double       simTime_m{ 0 };
	const double dt_m;
	const double simMin_m{ MIN_Z_SIM };
	const double simMax_m{ MAX_Z_SIM };
	const double tIon_m{ INITIAL_T_EV };
	const double tMag_m{ INITIAL_T_EV_MAG };
	const double vmean_m{ V_DIST_MEAN };
	const bool	 normalizedToRe_m{ true }; //change later if you want, make it an option
	
	std::vector<Particle*> particleTypes_m;

	//Satellites and data
	std::vector<Satellite*> satellites_m; //perhaps struct a satellite pointer and a particle pointer together??
	std::vector<std::vector<std::vector<std::vector<double>>>> satelliteData_m; //4D satelliteData[recorded measurement number][satellite number][attribute number][particle number]
	std::vector<std::vector<std::vector<double>>> preppedSatData_m;

	//GPU Memory Pointers
	std::vector<double*> gpuDblMemoryPointers_m { nullptr };
	std::vector<void*>	 gpuOtherMemoryPointers_m{ nullptr };

	//Flags
	bool initialized_m{ 0 };
	bool copied_m{ 0 };
	bool resultsPrepared_m{ 0 };

	//LogFile and Error Handling
	LogFile logFile_m{ "simulation.log", 20 };

	//Protected functions
	virtual void receiveSatelliteData();

	virtual void initializeFollowOn() { return; }
	virtual void copyDataToGPUFollowOn() { return; }
	virtual void iterateSimulationFollowOnPreLoop() { return; }
	virtual void iterateSimulationFollowOnInsideLoop() { return; }
	virtual void iterateSimulationFollowOnPostLoop() { return; }
	virtual void copyDataToHostFollowOn() { return; }
	virtual void freeGPUMemoryFollowOn() { return; }

public:
	Simulation(double dt, std::string rootdir):
		dt_m{ dt }, rootdir_m{ rootdir }
	{
		logFile_m.writeTimeDiffFromNow(0, "Simulation base class constructor");
	}

	virtual ~Simulation()
	{
		//
		//
		//Call from python
		writeSatelliteDataToCSV();
		//Call from python
		//
		//

		//Delete satellites and particles
		LOOP_OVER_1D_ARRAY(satellites_m.size(), delete satellites_m.at(iii);)
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), delete particleTypes_m.at(iii);)
		
		logFile_m.writeTimeDiffFromNow(0, "End Simulation Destructor");
	};
	//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	  getTime() { return simTime_m; }
	double	  getdt() { return dt_m; };
	void	  incTime() { simTime_m += dt_m; }

	size_t    getNumberOfParticleTypes() { return particleTypes_m.size(); }
	size_t    getNumberOfParticlesPerType(int partInd) { return particleTypes_m.at(partInd)->getNumberOfParticles(); }
	size_t    getNumberOfAttributesTracked(int partInd) { return particleTypes_m.at(partInd)->getNumberOfAttributes(); }

	bool	  areResultsPrepared() { return resultsPrepared_m; }
	bool	  getNormalized() { return normalizedToRe_m; }

	LogFile*  getLogFilePointer() { return &logFile_m; }

	virtual double getSimMin() { return simMin_m; }
	virtual double getSimMax() { return simMax_m; }

	double*   getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex, bool originalData);
	
	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double calculateBFieldAtZandTime(double z, double time) = 0;
	virtual double calculateEFieldAtZandTime(double z, double time) = 0;
	
	//Array tools
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& z, double mass);
	virtual void convertVPerpToMu(Particle* particle);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, double mass);
	virtual void convertMuToVPerp(Particle* particle);

	//Simulation management functions
	virtual void initializeSimulation() = 0;
	virtual void copyDataToGPU() = 0;
	virtual void iterateSimulation(int numberOfIterations) = 0;
	virtual void copyDataToHost() = 0;
	virtual void freeGPUMemory() = 0;
	virtual void prepareResults() = 0;

	//Satellite management functions
	virtual void	createSatellite(Particle* assignedPart, double altitude, bool upwardFacing, double** GPUdataPointers, bool elecTF, std::string name);
	virtual size_t	getNumberOfSatellites() { return satellites_m.size(); }
	virtual size_t	getNumberOfSatelliteMsmts() { return satelliteData_m.size(); }
	virtual double* getSatelliteDataPointers(int measurementInd, int satelliteInd, int attributeInd) { return satelliteData_m.at(measurementInd).at(satelliteInd).at(attributeInd).data(); }
	virtual void	writeSatelliteDataToCSV();

	virtual void    createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor);
	virtual Particle* getParticlePointer(int ind) { return particleTypes_m.at(ind); }
};//end class
#endif //end header guard