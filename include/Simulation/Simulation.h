#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <vector>
#include <memory> //smart pointers

#include "dlldefines.h"
#include "BField/allBModels.h"
#include "EField/allEModels.h"
#include "Particle/Particle.h"
#include "Satellite/Satellite.h"
#include "Log/Log.h"
#include "SimAttributes/SimAttributes.h"

using std::string;
using std::vector;
using std::unique_ptr;
using std::shared_ptr;

class Simulation
{
protected:
	//Structs and classes that fill various roles	
	enum class SimComp
	{
		BField,
		EField,
		Log,
		Particle,
		Satellite
	};

	struct TempSat
	{//Struct that holds data to create satellite - allows satellites to be added in any order, but ensure they are created before particles
		int particleInd;
		double altitude;
		bool upwardFacing;
		string name;

		TempSat(int partInd, double alt, bool upward, string nameStr) :
			particleInd{ partInd }, altitude{ alt }, upwardFacing{ upward }, name{ nameStr } {}
	};

	struct SatandPart
	{//Satellite needs particle-specific data associated with it, so this struct holds a shared_ptr to the particle
		unique_ptr<Satellite> satellite;
		shared_ptr<Particle>  particle;

		SatandPart(unique_ptr<Satellite> sat, shared_ptr<Particle> part) :
			satellite{ std::move(sat) }, particle{ std::move(part) } {}
	};

	//Simulation Characteristics
	double dt_m{ 0.0 };
	double simMin_m{ 0.0 };
	double simMax_m{ 0.0 };
	string saveRootDir_m{ "./" };
	double simTime_m{ 0.0 };

	//Simulation-specific classes tracked by Simulation
	vector<shared_ptr<Particle>>   particles_m;
	unique_ptr<BModel>             BFieldModel_m;
	unique_ptr<EField>             EFieldModel_m;
	vector<unique_ptr<TempSat>>    tempSats_m; //holds data until the GPU data arrays are allocated, allows the user more flexibility of when to call createSatellitesAPI
	vector<unique_ptr<SatandPart>> satPartPairs_m;

	//Flags
	bool initialized_m{ false };
	bool dataOnGPU_m  { false };
	bool saveReady_m  { false };

	//Attribute saving, Log, and Error Handling
	unique_ptr<SimAttributes> simAttr_m;
	unique_ptr<Log> log_m;

	//Protected functions
	virtual void createSatellite(TempSat* tmpsat, bool save = true);
	void incTime();
	virtual void printSimAttributes(int numberOfIterations, int itersBtwCouts, string GPUName);


public:
	Simulation(double dt, double simMin, double simMax, string saveRootDir);
	Simulation(string prevSimDir); //for loading previous simulation data
	virtual ~Simulation();

	Simulation(const Simulation&) = delete;
	Simulation& operator=(const Simulation&) = delete;

	class Environment;

	///One liner functions (usually access)
	double simtime() const;
	double dt()      const;
	double simMin()  const;
	double simMax()  const;

	int    getNumberOfParticleTypes()         const;
	int    getNumberOfSatellites()            const;
	int    getNumberOfParticles(int partInd)  const;
	int    getNumberOfAttributes(int partInd) const;
	string getParticleName(int partInd)       const;
	string getSatelliteName(int satInd)       const;
	int    getPartIndOfSat(int satInd)        const;

	Particle*  particle(int partInd)  const;
	Particle*  particle(string name)  const; //search for name, return particle
	Satellite* satellite(int satInd)  const;
	Satellite* satellite(string name) const; //search for name, return satellite
	BModel*    Bmodel()               const;
	EField*    Emodel()               const;

	virtual const vector<vector<double>>& getParticleData(int partInd, bool originalData);
	virtual const vector<vector<double>>& getSatelliteData(int satInd);

	///Forward decs for cpp file, or pure virtuals
	//Class creation functions
	virtual void createParticleType(string name, double mass, double charge, long numParts, string loadFilesDir = "", bool save = true);
	virtual void createTempSat(string partName, double altitude, bool upwardFacing, string name);
	virtual void createTempSat(int partInd, double altitude, bool upwardFacing, string name);
	virtual void setBFieldModel(string name, vector<double> args, bool save = true);
	virtual void setBFieldModel(unique_ptr<BModel> BModelptr);
	virtual void addEFieldModel(string name, vector<double> args, bool save = true);
	//virtual void addEFieldModel(unique_ptr<EModel> EModelptr) { EFieldModel_m->add(std::move(EModelptr)); }
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