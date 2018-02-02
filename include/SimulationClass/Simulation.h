#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "BField\AllBModels.h"
#include "EField\AllEModels.h"
#include "ParticleClass\Particle.h"
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "StandaloneTools\StandaloneTools.h"

struct SatandPart
{
	Satellite* satellite;
	Particle*  particle;

	~SatandPart() { delete satellite; }
};

struct TempSat
{
	int particleInd;
	double altitude;
	bool upwardFacing;
	std::string name;
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

	//Classes tracked by Simulation
	std::vector<Particle*> particleTypes_m;
	BField* BFieldModel_m{ nullptr };
	std::vector<EField*> EFieldElements;

	//Satellites and data
	std::vector<TempSat*> tempSats_m; //holds data until the GPU data arrays are allocated, allows the user more flexibility of when to call createSatellitesAPI
	std::vector<SatandPart*> satellites_m;
	std::vector<std::vector<std::vector<double>>> satelliteData_m; //3D satelliteData[satellite number][attribute number][particle number]
	
	//GPU Memory Pointers
	double* simConstants_d;
	void*   curandRNGStates_d; //void* is used instead of curandStateMRG32k3a* so I don't have to include the cuda headers here (which are used in other sections of code not requiring CUDA)

	//Flags
	bool initialized_m{ false };
	bool copied_m{ false };
	bool freedGPUMem_m{ false };
	bool resultsPrepared_m{ false };
	bool errorEncountered{ false };

	//LogFile and Error Handling
	LogFile* logFile_m{ nullptr };

	//Protected functions
	virtual void receiveSatelliteData(bool removeZeros);
	virtual void createSatellite(int partInd, double altitude, bool upwardFacing, std::string name);

	std::streambuf* cerrStrBuf{ std::cerr.rdbuf() };
	std::ofstream   errLogOut{ "errors.log" };

public:
	Simulation(double dt, double simMin, double simMax, double ionT, double magT, std::string rootdir, std::string logName="simulation.log"):
		dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, ionT_m{ ionT }, magT_m{ magT }, rootdir_m { rootdir }
	{
		logFile_m = new LogFile(logName, 20); //use std::unique_ptr / std::make_unique here
		//write log entry??

		std::cerr.rdbuf(errLogOut.rdbuf());
	}

	virtual ~Simulation()
	{
		std::cerr.rdbuf(cerrStrBuf);

		if (initialized_m && !freedGPUMem_m) { freeGPUMemory(); }
		delete BFieldModel_m;

		//Delete satellites and particles
		LOOP_OVER_1D_ARRAY(satellites_m.size(), delete satellites_m.at(iii));
		LOOP_OVER_1D_ARRAY(particleTypes_m.size(), delete particleTypes_m.at(iii));

		logFile_m->writeTimeDiffFromNow(0, "End Simulation Destructor");
		delete logFile_m;
	}//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	  getTime()   { return simTime_m; }//not API worthy, not checking mid simulation
	double	  getdt()     { return dt_m; }	   //not API worthy, probably passed in
	double    getSimMin() { return simMin_m; } //not API worthy, probably passed in
	double    getSimMax() { return simMax_m; } //not API worthy, probably passed in
	void	  incTime() { simTime_m += dt_m; } //not API worthy, never need to increment time from outside cpp

	size_t    getNumberOfParticleTypes() { return particleTypes_m.size(); } //not API worthy, probably passed in
	size_t    getNumberOfParticles(int partInd) { return particleTypes_m.at(partInd)->getNumberOfParticles(); } //not API worthy, probably passed in
	size_t    getNumberOfAttributes(int partInd) { return particleTypes_m.at(partInd)->getNumberOfAttributes(); } //not API worthy, probably passed in

	bool	  areResultsPrepared() { return resultsPrepared_m; } //do I even use this??

	LogFile*  getLogFilePointer() { return logFile_m; }
	double*   getPointerToParticleAttributeArray(int partInd, int attrInd, bool originalData);

	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double calculateBFieldAtZandTime(double s, double time) { return BFieldModel_m->getBFieldAtS(s, time); }
	virtual double calculateEFieldAtZandTime(double s, double time) { return 0.0; }//EFieldatZ(nullptr, s, time, 0.0, constE_m, useQSPS_m, false); }
	
	//Array tools
	/* Are API functions needed for these?  Prob not... */
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, double mass);
	virtual void convertVPerpToMu(Particle* particle);
	virtual void convertVPerpToMu(int partInd);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, double mass);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, std::vector<double>& t, double mass);
	virtual void convertMuToVPerp(Particle* particle);
	virtual void convertMuToVPerp(int partInd);

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations, int itersBtwCouts);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults(bool normalizeToRe);

	//Satellite management functions
	virtual size_t	getNumberOfSatellites() { return satellites_m.size(); }
	virtual size_t  getSatelliteNumberOfDetectedParticles(int satInd) { return satelliteData_m.at(satInd).at(0).size(); } ///add index checking
	virtual double* getSatelliteDataPointers(int satInd, int attrInd) { //some sort of check here to make sure you've received data
		return satelliteData_m.at(satInd).at(attrInd).data(); } //also check indicies
	virtual void	writeSatelliteDataToCSV();
	virtual void    createTempSat(int partInd, double altitude, bool upwardFacing, std::string name)
	{
		if (initialized_m) { throw std::runtime_error ("createTempSat: initializeSimulation has already been called, no satellite will be created of name " + name); }
		if (partInd >= particleTypes_m.size()) { throw std::out_of_range ("createTempSat: no particle at the specifed index " + std::to_string(partInd)); }
		tempSats_m.push_back(new TempSat{ partInd, altitude, upwardFacing, name });
	}

	virtual void	setBFieldModel(std::string name, std::vector<double> args);
	virtual void	setBFieldModelOther(BField* bfieldptr) { BFieldModel_m = bfieldptr; } //add API function for this

	//virtual void	addEFieldModel(std::string name, std::vector<double> args);
	//virtual void	addEFieldModelOther(EField* efieldptr);

	//add function that saves simulation constants, data, etc to disk

	//going to restructure this
	virtual void    loadCompletedSimData(std::string fileDir, std::vector<std::string> partNames, std::vector<std::string> attrNames, std::vector<std::string> satNames, int numParts);

	virtual void    createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir="");
	virtual Particle* getParticlePointer(int ind) { return particleTypes_m.at(ind); }
};//end class
#endif //end header guard