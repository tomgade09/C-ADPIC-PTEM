#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <memory> //smart pointers
#include "BField\allBModels.h"
#include "EField\allEModels.h"
#include "ParticleClass\Particle.h"
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "StandaloneTools\StandaloneTools.h"
#include "ErrorHandling\simExceptionMacros.h"

struct SatandPart
{
	std::unique_ptr<Satellite> satellite;
	std::shared_ptr<Particle>  particle;

	SatandPart(std::unique_ptr<Satellite> sat, std::shared_ptr<Particle> part) :
		satellite{ std::move(sat) }, particle{ std::move(part) } {}
};

struct TempSat
{
	int particleInd;
	double altitude;
	bool upwardFacing;
	std::string name;

	TempSat(int partInd, double alt, bool upward, std::string nameStr) :
		particleInd{ partInd }, altitude{ alt }, upwardFacing{ upward }, name{ nameStr } {}
};

class Simulation
{
protected:
	//Simulation Characteristics
	std::string  saveRootDir_m;
	double       simTime_m{ 0.0 };
	double       dt_m;
	double       simMin_m;
	double       simMax_m;
	double       ionT_m;
	double       magT_m;
	double       vmean_m{ 0.0 };

	//Classes tracked by Simulation
	std::vector<std::shared_ptr<Particle>> particleTypes_m; //need shared_ptr because it will be assigned to SatandPart
	std::unique_ptr<BField> BFieldModel_m;
	std::unique_ptr<EField> EFieldModel_m;
	BField** BFieldModel_d{ nullptr };
	EField** EFieldModel_d{ nullptr };

	//Satellites and data
	std::vector<std::unique_ptr<TempSat>> tempSats_m; //holds data until the GPU data arrays are allocated, allows the user more flexibility of when to call createSatellitesAPI
	std::vector<std::unique_ptr<SatandPart>> satellites_m; //don't know if I'm crazy about this solution
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
	std::unique_ptr<LogFile> logFile_m;

	//Protected functions
	virtual void receiveSatelliteData(bool removeZeros);
	virtual void createSatellite(int partInd, double altitude, bool upwardFacing, std::string name);
	virtual void writeCharsToFiles(std::vector<double> chars, std::vector<std::string> charNames, std::string className, std::string folderFromSave="/_chars/");

	std::streambuf* cerrBufBak{ std::cerr.rdbuf() };
	std::ofstream   cerrLogOut;

	Simulation(std::string prevSaveDir) : saveRootDir_m{ prevSaveDir }
	{
		logFile_m = std::make_unique<LogFile>(saveRootDir_m + "/reloadSim.log", 20);
	}

public:
	Simulation(double dt, double simMin, double simMax, double ionT, double magT, std::string saveRootDir, std::string logName="simulation.log"):
		dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, ionT_m{ ionT }, magT_m{ magT }, saveRootDir_m { saveRootDir + "/" }
	{
		logFile_m = std::make_unique<LogFile>(saveRootDir_m + logName, 20);
		cerrLogOut.open(saveRootDir_m + "errors.log");
		std::cerr.rdbuf(cerrLogOut.rdbuf()); //set cerr output to "errors.log"

		writeCharsToFiles( {dt_m, simMin_m, simMax_m, ionT_m, magT_m, vmean_m}, {"dt", "simMin", "simMax", "T_ion", "T_mag", "v_mean"}, "Simulation");
	}

	virtual ~Simulation()
	{
		std::cerr.rdbuf(cerrBufBak); //restore cerr to normal

		if (initialized_m && !freedGPUMem_m) { freeGPUMemory(); }

		logFile_m->writeTimeDiffFromNow(0, "End Simulation Destructor");
	}//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	  getTime()   { return simTime_m; }
	double	  getdt()     { return dt_m; }
	double    getSimMin() { return simMin_m; }
	double    getSimMax() { return simMax_m; }
	void	  incTime()   { simTime_m += dt_m; }

	size_t    getNumberOfParticleTypes()         { return particleTypes_m.size(); }
	size_t    getNumberOfParticles(int partInd)  { return particleTypes_m.at(partInd)->getNumberOfParticles(); }
	size_t    getNumberOfAttributes(int partInd) { return particleTypes_m.at(partInd)->getNumberOfAttributes(); }
	const char* getParticleName(int partInd)     { return particleTypes_m.at(partInd)->getName().c_str(); }

	bool	  areResultsPrepared() { return resultsPrepared_m; } //do I even use this??

	LogFile*  getLogFilePointer() { return logFile_m.get(); }
	double*   getPointerToParticleAttributeArray(int partInd, int attrInd, bool originalData);

	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double getBFieldAtS(double s, double time) { return BFieldModel_m->getBFieldAtS(s, time); }
	virtual double getEFieldAtS(double s, double time) { if (EFieldModel_m == nullptr) { return 0.0; } return EFieldModel_m->getEFieldAtS(s, time); }
	
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
	const char* getSatelliteName(int satInd) { return satellites_m.at(satInd)->satellite->getName().c_str(); }
	virtual void	writeSatelliteDataToCSV();
	virtual void    createTempSat(int partInd, double altitude, bool upwardFacing, std::string name);

	virtual void	setBFieldModel(std::string name, std::vector<double> args);
	virtual void	setBFieldModelOther(std::unique_ptr<BField> bfieldptr) { BFieldModel_m = std::move(bfieldptr); } //add API function for this

	void Simulation::addEFieldModel(std::string name, std::vector<double> args, std::string altMinMaxStr, std::string magnitudeStr);
	//virtual void	addEFieldModelOther(std::unique_ptr<EField> efieldptr);

	//add function that saves simulation constants, data, etc to disk

	//going to restructure this
	virtual void    loadCompletedSimData(std::string fileDir, std::vector<std::string> partNames, std::vector<std::string> attrNames, std::vector<std::string> satNames, int numParts);

	virtual void    createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, int posDims, int velDims, double normFactor, std::string loadFilesDir="");
	virtual Particle* getParticlePointer(int ind) { return particleTypes_m.at(ind).get(); }
};//end class
#endif //end header guard