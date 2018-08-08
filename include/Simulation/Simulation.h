#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <memory> //smart pointers

#include "dlldefines.h"
#include "BField/allBModels.h"
#include "EField/allEModels.h"
#include "Particle/Particle.h"
#include "Satellite/Satellite.h"
#include "LogFile/LogFile.h"
#include "SimAttributes/SimAttributes.h"

class Simulation
{
protected:
	//Structs that fill various roles	
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

	//Simulation-specific classes tracked by Simulation
	std::vector<std::shared_ptr<Particle>>   particles_m;
	std::unique_ptr<BField>                  BFieldModel_m;
	std::unique_ptr<EField>                  EFieldModel_m;
	std::vector<std::unique_ptr<TempSat>>    tempSats_m; //holds data until the GPU data arrays are allocated, allows the user more flexibility of when to call createSatellitesAPI
	std::vector<std::unique_ptr<SatandPart>> satPartPairs_m;
	
	//GPU Pointers
	BField** BFieldModel_d{ nullptr };
	EField** EFieldModel_d{ nullptr };

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
	Simulation(double dt, double simMin, double simMax, std::string saveRootDir);
	Simulation(std::string prevSimDir); //for loading previous simulation data
	virtual ~Simulation();

	Simulation(const Simulation&) = delete;
	Simulation& operator=(const Simulation&) = delete;

	///One liner functions (usually access)
	double	    simtime() const { return simTime_m; }
	double	    dt()      const { return dt_m; }
	double      simMin()  const { return simMin_m; }
	double      simMax()  const { return simMax_m; }

	int         getNumberOfParticleTypes()         const { return (int)particles_m.size(); }
	int         getNumberOfSatellites()            const { return (int)satPartPairs_m.size(); }
	int         getNumberOfParticles(int partInd)  const { return (int)particles_m.at(partInd)->getNumberOfParticles(); }
	int         getNumberOfAttributes(int partInd) const { return (int)particles_m.at(partInd)->getNumberOfAttributes(); }
	std::string getParticleName(int partInd)       const { return particles_m.at(partInd)->name(); }
	std::string getSatelliteName(int satInd)       const { return satPartPairs_m.at(satInd)->satellite->name(); }
	int         getPartIndOfSat(int satInd)        const { return tempSats_m.at(satInd)->particleInd; }

	LogFile*    log()                       const { return logFile_m.get(); }
	Particle*   particle(int partInd)       const { return particles_m.at(partInd).get(); }
	Particle*   particle(std::string name)  const;  //search for name, return particle
	Satellite*  satellite(int satInd)       const { return satPartPairs_m.at(satInd)->satellite.get(); }
	Satellite*  satellite(std::string name) const; //search for name, return satellite
	BField*     Bmodel()                    const { return BFieldModel_m.get(); }
	EField*     Emodel()                    const { return EFieldModel_m.get(); }

	#define VEC(T) std::vector<T> //quick, lazy stand-in, easier on the eyes
	virtual const VEC(VEC(double))&      getParticleData(int partInd, bool originalData);
	virtual const VEC(VEC(VEC(double)))& getSatelliteData(int satInd);
	#undef VEC

	///Forward decs for cpp file, or pure virtuals
	//Class creation functions
	virtual void createParticleType(std::string name, double mass, double charge, long numParts, std::string loadFilesDir = "", bool save = true);
	virtual void createTempSat(std::string partName, double altitude, bool upwardFacing, std::string name);
	virtual void createTempSat(int partInd, double altitude, bool upwardFacing, std::string name);
	virtual void setBFieldModel(std::string name, std::vector<double> args, bool save = true);
	virtual void setBFieldModel(std::unique_ptr<BField> bfieldptr) { BFieldModel_m = std::move(bfieldptr); } //add API function for this
	virtual void addEFieldModel(std::string name, std::vector<double> args, bool save = true);
	//virtual void addEFieldModel(std::unique_ptr<EElem> eelemptr) { EFieldModel_m->add(std::move(eelemptr)); }
	//virtual void clearEFieldModels(){ EFieldModel_m = std::make_unique<EField>(); }
	

	//Fields
	virtual double getBFieldAtS(double s, double time) const;
	virtual double getEFieldAtS(double s, double time) const;


	//Simulation management functions
	virtual void initializeSimulation();
	virtual void __iterateSimCPU(int numberOfIterations, int checkDoneEvery);
	virtual void iterateSimulation(int numberOfIterations, int itersBtwCouts);
	virtual void saveDataToDisk();
	virtual void freeGPUMemory();
	virtual void resetSimulation(bool fields = false);
};//end class
#endif //end header guard