#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <vector>
#include <chrono>
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "StandaloneTools\binaryfiletools.h"
#include "include\_simulationvariables.h" //remove this later, replace with passed in variables

//Cannot nest these - outer loops all use iii as a variable name
#define LOOP_OVER_3D_ARRAY(d1, d2, d3, x) for (int iii = 0; iii < d1; iii++)\
	{ for (int jjj = 0; jjj < d2; jjj++) \
		{ for (int kk = 0; kk < d3; kk++) {x} } }

#define LOOP_OVER_2D_ARRAY(d1, d2, x) for (int iii = 0; iii < d1; iii++)\
	{ for (int jjj = 0; jjj < d2; jjj++) {x} }

#define LOOP_OVER_1D_ARRAY(d1, x) for (int iii = 0; iii < d1; iii++) {x}

class Simulation
{
protected:
	//Simulation Characteristics
	std::string  rootdir_m;
	double       simTime_m{ 0 };
	const double dt_m;
	const int	 numberOfParticleTypes_m;
	const int	 numberOfParticlesPerType_m;
	const int	 numberOfAttributesTracked_m;
	const double simMin_m{ MIN_Z_SIM };
	const double simMax_m{ MAX_Z_SIM };
	const bool	 normalizedToRe_m{ true }; //change later if you want, make it an option

	//Particle data arrays and log file
	double*** particles_m{ form3Darray() };
	double*** particlesorig_m{ form3Darray() }; //initial data
	std::vector<double> mass_m;
	std::vector<void*> otherMemoryPointers_m;
	std::vector<Satellite*> satellites_m;
	std::vector<std::vector<std::vector<double*>>> satelliteData_m; //4D satelliteData[recorded measurement number][satellite number][attribute number][particle number]
	LogFile logFile_m{ "simulation.log", 20 };
	
	//GPU Memory Pointers
	std::vector<double*> gpuDblMemoryPointers_m { nullptr };
	std::vector<void*> gpuOtherMemoryPointers_m{ nullptr };

	//Flags
	bool initialized_m{ 0 };
	bool copied_m{ 0 };
	bool resultsPrepared_m{ 0 };
	bool mu_m{ 1 }; //comes off the gpu as mu, need to change if a dist of vperp is generated

	virtual void receiveSatelliteData();
	virtual double*** form3Darray();
	
public:
	Simulation(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir):
		numberOfParticleTypes_m{ numberOfParticleTypes }, numberOfParticlesPerType_m{ numberOfParticlesPerType },
		numberOfAttributesTracked_m{ numberOfAttributesTracked }, dt_m{ dt }, rootdir_m{ rootdir }
	{
		logFile_m.writeTimeDiffFromNow(0, "Simulation base class constructor");

		//Populate mass array
		mass_m.reserve(2);
		mass_m[0] = MASS_ELECTRON;
		mass_m[1] = MASS_PROTON;

		//Allocate room in vectors for GPU Memory Pointers
		gpuDblMemoryPointers_m.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m); //holds pointers to GPU memory for particle attributes
		gpuOtherMemoryPointers_m.reserve(2); //LUT data[0], random number generator[1]
	}

	virtual ~Simulation()
	{
		//Save init particle distributions to disk
		std::string fold{ "./bins/particles_init/" };
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, saveParticleAttributeToDisk(particlesorig_m[iii][jjj], numberOfParticlesPerType_m, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());)

		//Save final particle distributions to disk
		fold = "./bins/particles_final/";
		names = { "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, saveParticleAttributeToDisk(particles_m[iii][jjj], numberOfParticlesPerType_m, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());)
		
		//Delete arrays
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, delete[] particles_m[iii][jjj];)
		LOOP_OVER_1D_ARRAY(numberOfParticleTypes_m, delete[] particles_m[iii];)
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, delete[] particlesorig_m[iii][jjj];)
		LOOP_OVER_1D_ARRAY(numberOfParticleTypes_m, delete[] particlesorig_m[iii];)
		delete[] particles_m;
		delete[] particlesorig_m;

		LOOP_OVER_3D_ARRAY(satelliteData_m.size(), satelliteData_m[0].size(), satelliteData_m[0][0].size(), delete[] satelliteData_m[iii][jjj][kk];)
		
		//Delete satellites
		LOOP_OVER_1D_ARRAY(satellites_m.size(), delete satellites_m[iii];)

		logFile_m.writeTimeDiffFromNow(0, "End Simulation Destructor");
	};
	//Generally, when I'm done with this class, I'm done with the whole program, so the memory is returned anyway, but still good to get in the habit of returning memory

	///One liner functions (usually access)
	double	  getTime() { return simTime_m; }
	double	  getdt() { return dt_m; };
	void	  incTime() { simTime_m += dt_m; }

	int		  getNumberOfParticleTypes() { return numberOfParticleTypes_m; }
	int		  getNumberOfParticlesPerType() { return numberOfParticlesPerType_m; }
	int		  getNumberOfAttributesTracked() { return numberOfAttributesTracked_m; }

	bool	  areResultsPrepared() { return resultsPrepared_m; }
	bool	  getNormalized() { return normalizedToRe_m; }
	double*** getPointerTo3DParticleArray() { return particles_m; }
	double**  getPointerToSingleParticleTypeArray(int index) { if (index >= numberOfParticleTypes_m) { return nullptr; } return particles_m[index]; } //eventually print a warning so the user knows
	double*   getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex, bool originalData) { if ((partIndex > numberOfParticleTypes_m) || (attrIndex > numberOfAttributesTracked_m)) { 
		std::cout << numberOfParticleTypes_m << " particle types and " << numberOfAttributesTracked_m << " attributes per particle.  One index is out of bounds.\n"; return nullptr; } 
		if (originalData) { return particlesorig_m[partIndex][attrIndex]; } else { return particles_m[partIndex][attrIndex]; } }

	LogFile*  getLogFilePointer() { return &logFile_m; }

	virtual double getSimMin() { return simMin_m; }
	virtual double getSimMax() { return simMax_m; }

	///Forward decs for cpp file, or pure virtuals
	//Field tools
	virtual double calculateBFieldAtZandTime(double z, double time) = 0;
	virtual double calculateEFieldAtZandTime(double z, double time) = 0;
	
	//Array tools
	virtual void convertVPerpToMu(int vind = 1, int zind = 2);
	virtual void convertMuToVPerp(int vind = 1, int zind = 2);

	//Simulation management functions
	virtual void initializeSimulation() = 0;
	virtual void copyDataToGPU() = 0;
	virtual void iterateSimulation(int numberOfIterations) = 0;
	virtual void copyDataToHost() = 0;
	virtual void freeGPUMemory() = 0;
	virtual void prepareResults() = 0;

	//Satellite management functions
	virtual void	createSatellite(double altitude, bool upwardFacing, double** GPUdataPointers, bool elecTF, std::string name);
	virtual int		getNumberOfSatellites() { return satellites_m.size(); }
	virtual int		getNumberOfSatelliteMsmts() { return satelliteData_m.size(); }
	virtual double* getSatelliteDataPointers(int measurementInd, int satelliteInd, int attributeInd) { return satelliteData_m[measurementInd][satelliteInd][attributeInd]; }

};//end class
#endif //end header guard