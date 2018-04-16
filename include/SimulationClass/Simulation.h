#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <vector>
#include <memory> //smart pointers
#include "BField\allBModels.h"
#include "EField\allEModels.h"
#include "ParticleClass\Particle.h"
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "FileIO\fileIO.h"
#include "utils\loopmacros.h"
#include "utils\numerical.h"
#include "SimAttributes\SimAttributes.h"
#include "ErrorHandling\simExceptionMacros.h"

class Simulation
{
protected:
	//Structs that fill various funcitons	
	struct TempSat
	{//Struct that holds data to create satellite - allows satellites to be created before particles through API
		int particleInd;
		double altitude;
		bool upwardFacing;
		std::string name;

		TempSat(int partInd, double alt, bool upward, std::string nameStr) :
			particleInd{ partInd }, altitude{ alt }, upwardFacing{ upward }, name{ nameStr } {}
	};

	struct SatandPart
	{//Satellite needs particle-specific data associated with it, so this struct holds a shared_ptr to the particle
		std::unique_ptr<Satellite> satellite;
		std::shared_ptr<Particle>  particle;

		SatandPart(std::unique_ptr<Satellite> sat, std::shared_ptr<Particle> part) :
			satellite{ std::move(sat) }, particle{ std::move(part) } {}
	};


	//Simulation Characteristics
	std::string  saveRootDir_m;
	double       simTime_m{ 0.0 };
	double       dt_m;
	double       simMin_m;
	double       simMax_m;

	//Classes tracked by Simulation
	std::vector<std::shared_ptr<Particle>> particles_m; //need shared_ptr because it will be assigned to SatandPart
	std::unique_ptr<BField> BFieldModel_m;
	std::unique_ptr<EField> EFieldModel_m;
	BField** BFieldModel_d{ nullptr };
	EField** EFieldModel_d{ nullptr };

	//Satellites and data
	std::vector<std::unique_ptr<TempSat>> tempSats_m; //holds data until the GPU data arrays are allocated, allows the user more flexibility of when to call createSatellitesAPI
	std::vector<std::unique_ptr<SatandPart>> satellites_m; //don't know if I'm crazy about this solution

	//Flags
	bool initialized_m{ false };
	bool dataOnGPU_m  { false };
	bool saveReady_m  { false };

	//Attribute saving, LogFile and Error Handling
	std::unique_ptr<SimAttributes> simAttr_m;
	std::unique_ptr<LogFile> logFile_m;
	std::streambuf* cerrBufBak{ std::cerr.rdbuf() };
	std::ofstream   cerrLogOut;

	//Protected functions
	virtual void createSatellite(TempSat* tmpsat, bool save = true);
	void incTime() { simTime_m += dt_m; }
	virtual void printSimAttributes(int numberOfIterations, int itersBtwCouts);


public:
	Simulation(double dt, double simMin, double simMax, std::string saveRootDir) :
		dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, saveRootDir_m{ saveRootDir + "/" },
		simAttr_m{ std::make_unique<SimAttributes>("Simulation.attr") },
		logFile_m{ std::make_unique<LogFile>(saveRootDir_m + "simulation.log", 20) }
	{
		size_t free, total;
		CUDA_API_ERRCHK(cudaMemGetInfo(&free, &total));
		std::cout << "Pre-Initialize cudaMemGetInfo: free: " << free << ", total: " << total << std::endl;

		cerrLogOut.open(saveRootDir_m + "errors.log");
		std::cerr.rdbuf(cerrLogOut.rdbuf()); //set cerr output to "errors.log"

		simAttr_m->addData("Simulation", "", {}, {}, { "dt", "simMin", "simMax" }, { dt_m, simMin_m, simMax_m });
	}

	Simulation(std::string prevSimDir); //for loading previous simulation data

	virtual ~Simulation()
	{
		std::cerr.rdbuf(cerrBufBak); //restore cerr to normal

		if (saveReady_m) { saveDataToDisk(); } //save data if it hasn't been done

		logFile_m->writeTimeDiffFromNow(0, "End Simulation Destructor");

		size_t free, total;
		CUDA_API_ERRCHK(cudaMemGetInfo(&free, &total));
		std::cout << "Post-Free      cudaMemGetInfo: free: " << free << ", total: " << total << std::endl;
	}//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	    simtime(){ return simTime_m; }
	double	    dt()     { return dt_m; }
	double      simMin() { return simMin_m; }
	double      simMax() { return simMax_m; }

	size_t      getNumberOfParticleTypes()         { return particles_m.size(); }
	size_t      getNumberOfSatellites()            { return satellites_m.size(); }
	size_t      getNumberOfParticles(int partInd)  { return particles_m.at(partInd)->getNumberOfParticles(); }
	size_t      getNumberOfAttributes(int partInd) { return particles_m.at(partInd)->getNumberOfAttributes(); }
	std::string getParticleName(int partInd)       { return particles_m.at(partInd)->name(); }
	std::string getSatelliteName(int satInd)       { return satellites_m.at(satInd)->satellite->name(); }

	LogFile*    log()                 { return logFile_m.get(); }
	Particle*   particle(int partInd) { return particles_m.at(partInd).get(); }
	Satellite*  satellite(int satInd) { return satellites_m.at(satInd)->satellite.get(); }
	BField*     Bmodel()              { return BFieldModel_m.get(); }
	EField*     Emodel()              { return EFieldModel_m.get(); }

	#define VEC(T) std::vector<T> //quick, lazy stand-in, easier on the eyes
	const virtual VEC(VEC(double))&       getParticleData(int partInd, bool originalData);
	const virtual VEC(VEC(VEC(double)))&  getSatelliteData(int satInd);
	#undef VEC

	///Forward decs for cpp file, or pure virtuals
	//Class creation functions
	virtual void   createParticleType(std::string name, std::vector<std::string> attrNames, double mass, double charge, long numParts, std::string loadFilesDir = "", bool save = true);
	virtual void   createTempSat(int partInd, double altitude, bool upwardFacing, std::string name);
	virtual void   setBFieldModel(std::string name, std::vector<double> args, bool save = true);
	virtual void   setBFieldModel(std::unique_ptr<BField> bfieldptr) { BFieldModel_m = std::move(bfieldptr); } //add API function for this
	virtual void   addEFieldModel(std::string name, std::vector<double> args, bool save = true);
	virtual void   addEFieldModel(std::unique_ptr<EElem> eelem) { EFieldModel_m->add(std::move(eelem)); }
	
	//vperp <-> mu
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, double mass);
	virtual void convertVPerpToMu(std::vector<double>& vperp, std::vector<double>& s, std::vector<double>& t, double mass);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, double mass);
	virtual void convertMuToVPerp(std::vector<double>& mu, std::vector<double>& s, std::vector<double>& t, double mass);
	
	//Other utilities
	virtual void saveDataToDisk();
	virtual void resetSimulation(bool fields=false);

	//Fields
	virtual double getBFieldAtS(double s, double time) { return BFieldModel_m->getBFieldAtS(s, time); } //need to throw if called too early?
	virtual double getEFieldAtS(double s, double time) { return EFieldModel_m->getEFieldAtS(s, time); }

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void iterateSimulation(int numberOfIterations, int itersBtwCouts);
	virtual void freeGPUMemory();
};//end class
#endif //end header guard