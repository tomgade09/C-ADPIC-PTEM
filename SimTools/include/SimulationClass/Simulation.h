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

	~SatandPart() { delete satellite; delete particle; }
};

class Simulation
{
protected:
	//Simulation Characteristics
	std::string  rootdir_m;
	double       simTime_m{ 0.0 };
	const double dt_m;
	const double simMin_m;
	const double simMax_m;
	const double ionT_m;
	const double magT_m;
	const double vmean_m{ 0.0 };
	double constE_m{ 0.0 };
	
	bool useQSPS_m  { false };
	bool useAlfLUT_m{ false };
	bool useAlfCal_m{ false };

	std::vector<Particle*> particleTypes_m;

	//Satellites and data
	std::vector<SatandPart*> satellites_m;
	std::vector<std::vector<std::vector<double>>> satelliteData_m; //3D satelliteData[satellite number][attribute number][particle number]
	
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
	virtual void receiveSatelliteData(bool removeZeros);

	virtual void initializeFollowOn() { return; }
	virtual void copyDataToGPUFollowOn() { return; }
	virtual void iterateSimulationFollowOnPreLoop() { return; }
	virtual void iterateSimulationFollowOnInsideLoop() { return; }
	virtual void iterateSimulationFollowOnPostLoop() { return; }
	virtual void copyDataToHostFollowOn() { return; }
	virtual void freeGPUMemoryFollowOn() { return; }

public:
	Simulation(double dt, double simMin, double simMax, double ionT, double magT, std::string rootdir):
		dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, ionT_m{ ionT }, magT_m{ magT }, rootdir_m { rootdir }
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
		//LOOP_OVER_1D_ARRAY(particleTypes_m.size(), delete particleTypes_m.at(iii);)
		
		logFile_m.writeTimeDiffFromNow(0, "End Simulation Destructor");
	};
	//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	  getTime() { return simTime_m; }
	double	  getdt() { return dt_m; };
	double    getSimMin() { return simMin_m; }
	double    getSimMax() { return simMax_m; }
	void	  incTime() { simTime_m += dt_m; }
	void      setQSPS(double constE) { constE_m = constE; useQSPS_m = true; }

	size_t    getNumberOfParticleTypes() { return particleTypes_m.size(); }
	size_t    getNumberOfParticles(int partInd) { return particleTypes_m.at(partInd)->getNumberOfParticles(); }
	size_t    getNumberOfAttributes(int partInd) { return particleTypes_m.at(partInd)->getNumberOfAttributes(); }

	bool	  areResultsPrepared() { return resultsPrepared_m; }
	//bool	  getNormalized() { return normalizedToRe_m; }

	LogFile*  getLogFilePointer() { return &logFile_m; }
	double*   getPointerToParticleAttributeArray(int partIndex, int attrIndex, bool originalData);

	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double calculateBFieldAtZandTime(double z, double time) { return BFieldatZ(z, time); }
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(nullptr, z, time, 0.0, constE_m, useQSPS_m, false); }
	
	//Array tools
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& z, double mass);
	virtual void convertVPerpToMu(Particle* particle);
	virtual void convertVPerpToMu(int partInd);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, double mass);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& z, std::vector<double>& t, double mass);
	virtual void convertMuToVPerp(Particle* particle);
	virtual void convertMuToVPerp(int partInd);

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults(bool normalizeToRe);

	//Satellite management functions
	virtual void	createSatellite(int partInd, double altitude, bool upwardFacing, std::string name);
	virtual size_t	getNumberOfSatellites() { return satellites_m.size(); }
	//virtual size_t	getNumberOfSatelliteMsmts() { return satelliteData_m.size(); }
	virtual double* getSatelliteDataPointers(int satelliteInd, int attributeInd) { //some sort of check here to make sure you've received data
		return satelliteData_m.at(satelliteInd).at(attributeInd).data(); }
	virtual void	writeSatelliteDataToCSV();

	virtual void    createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir="");
	virtual Particle* getParticlePointer(int ind) { return particleTypes_m.at(ind); }
};//end class
#endif //end header guard