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

struct SatandPart
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
	
	bool useQSPS_m  { false };
	bool useAlfLUT_m{ false };
	bool useAlfCal_m{ false };

	std::vector<Particle*> particleTypes_m;

	//Satellites and data
	std::vector<Satellite*> satellites_m;
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
	double	  getTime() { return simTime_m; } //API
	double	  getdt() { return dt_m; }; //API
	double    getSimMin() { return simMin_m; } //API
	double    getSimMax() { return simMax_m; } //API
	void	  incTime() { simTime_m += dt_m; } //API
	void      setUseQSPS(bool qsps) { useQSPS_m = qsps; } //API

	size_t    getNumberOfParticleTypes() { return particleTypes_m.size(); } //API
	size_t    getNumberOfParticles(int partInd) { return particleTypes_m.at(partInd)->getNumberOfParticles(); } //API
	size_t    getNumberOfAttributes(int partInd) { return particleTypes_m.at(partInd)->getNumberOfAttributes(); } //API

	bool	  areResultsPrepared() { return resultsPrepared_m; } //API
	bool	  getNormalized() { return normalizedToRe_m; } //API, use to track whether or not values are normalized

	LogFile*  getLogFilePointer() { return &logFile_m; }
	double*   getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex, bool originalData); //API

	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double calculateBFieldAtZandTime(double z, double time) { return BFieldatZ(z, time); } //API
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(nullptr, z, time, 0.0, useQSPS_m, false); } //API
	
	//Array tools
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& z, double mass);
	virtual void convertVPerpToMu(Particle* particle);
	virtual void convertVPerpToMu(int partInd); //API
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, double mass);
	virtual void convertMuToVPerp(Particle* particle);
	virtual void convertMuToVPerp(int partInd); //API

	//Simulation management functions
	virtual void initializeSimulation(); //API
	virtual void copyDataToGPU(); //API
	virtual void iterateSimulation(int numberOfIterations); //API
	virtual void copyDataToHost(); //API
	virtual void freeGPUMemory(); //API
	virtual void prepareResults(); //API

	//Satellite management functions
	virtual void	createSatellite(int particleInd, double altitude, bool upwardFacing, double** GPUdataPointers, bool elecTF, std::string name); //API
	virtual size_t	getNumberOfSatellites() { return satellites_m.size(); } //API
	virtual size_t	getNumberOfSatelliteMsmts() { return satelliteData_m.size(); } //API
	virtual double* getSatelliteDataPointers(int measurementInd, int satelliteInd, int attributeInd) { return satelliteData_m.at(measurementInd).at(satelliteInd).at(attributeInd).data(); } //API
	virtual void	writeSatelliteDataToCSV(); //API

	virtual void    createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor); //API
	virtual Particle* getParticlePointer(int ind) { return particleTypes_m.at(ind); } //API - maybe??
};//end class
#endif //end header guard