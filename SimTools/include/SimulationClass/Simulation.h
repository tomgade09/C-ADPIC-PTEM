#ifndef SIMULATIONCLASS_H
#define SIMULATIONCLASS_H

#include <iostream>
#include <vector>
#include <chrono>
#include "SatelliteClass\Satellite.h"
#include "LogFile\LogFile.h"
#include "StandaloneTools\StandaloneTools.h"
#include "include\_simulationvariables.h" //remove this later, replace with passed in variables

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

	//Const variables
	const int LENGTHSATDATA{ (numberOfAttributesTracked_m + 1) * numberOfParticlesPerType_m };

	//Particle data arrays and log file
	std::vector<std::vector<std::vector<double>>> particles_m;
	std::vector<std::vector<std::vector<double>>> partInitData_m;
	std::vector<double> mass_m;

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
	bool mu_m{ 0 };

	//LogFile and Error Handling
	LogFile logFile_m{ "simulation.log", 20 };
	//ErrorHandler errors{}; //future feature to be added

	//Protected functions
	virtual void receiveSatelliteData();

public:
	Simulation(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir, bool loadDist=false):
		numberOfParticleTypes_m{ numberOfParticleTypes }, numberOfParticlesPerType_m{ numberOfParticlesPerType },
		numberOfAttributesTracked_m{ numberOfAttributesTracked }, dt_m{ dt }, rootdir_m{ rootdir }
	{
		logFile_m.writeTimeDiffFromNow(0, "Simulation base class constructor");

		particles_m = form3DvectorArray(numberOfParticleTypes_m, numberOfAttributesTracked_m, numberOfParticlesPerType_m);
		partInitData_m = form3DvectorArray(numberOfParticleTypes_m, numberOfAttributesTracked_m, numberOfParticlesPerType_m);

		//Populate mass array
		mass_m.resize(2);
		mass_m.at(0) = MASS_ELECTRON;
		mass_m.at(1) = MASS_PROTON;

		//Allocate room in vectors for GPU Memory Pointers
		gpuDblMemoryPointers_m.resize(2 * numberOfParticleTypes_m + 2);
		gpuOtherMemoryPointers_m.resize(2 * numberOfParticleTypes_m + 2);
		satelliteData_m.reserve(100); //not resize...Don't know the exact size here so need to use push_back

		std::cout << "Need to change flag in Simulation constructor.\n\n\n";
		if (false)//need to replace this condition with the passed in variable above
		{
			std::cout << "Loading initial particle data.\n";
			std::string fold{ "./../../in/epitchbin/" };
			std::vector<std::string> names{ "e_vpara.bin", "e_vperp.bin", "e_z.bin", "i_vpara.bin", "i_vperp.bin", "i_z.bin" };
			LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, fileIO::readDblBin(particles_m.at(iii).at(jjj), fold + names.at(iii * numberOfAttributesTracked_m + jjj), numberOfParticlesPerType_m);)
			LOOP_OVER_3D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, numberOfParticlesPerType_m, particles_m.at(iii).at(jjj).at(kk) *= RADIUS_EARTH;)
		}
	}

	virtual ~Simulation()
	{
		writeSatelliteDataToCSV();

		//Save init particle distributions to disk
		std::string fold{ "./bins/particles_init/" };
		std::vector<std::string> names{ "e_vpara.bin", "e_vperp.bin", "e_z.bin", "i_vpara.bin", "i_vperp.bin", "i_z.bin" };
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, fileIO::writeDblBin(partInitData_m.at(iii).at(jjj), fold + names.at(iii * numberOfAttributesTracked_m + jjj), numberOfParticlesPerType_m);)

		//Save final particle distributions to disk
		fold = "./bins/particles_final/";
		LOOP_OVER_2D_ARRAY(numberOfParticleTypes_m, numberOfAttributesTracked_m, fileIO::writeDblBin(particles_m.at(iii).at(jjj), fold + names.at(iii * numberOfAttributesTracked_m + jjj), numberOfParticlesPerType_m);)
		
		//Delete satellites
		LOOP_OVER_1D_ARRAY(satellites_m.size(), delete satellites_m.at(iii);)
		
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
	double*   getPointerToSingleParticleAttributeArray(int partIndex, int attrIndex, bool originalData) { if ((partIndex > numberOfParticleTypes_m) || (attrIndex > numberOfAttributesTracked_m)) { 
		std::cout << numberOfParticleTypes_m << " particle types and " << numberOfAttributesTracked_m << " attributes per particle.  One index is out of bounds.\n"; return nullptr; } 
		if (originalData) { return partInitData_m.at(partIndex).at(attrIndex).data(); } else { return particles_m.at(partIndex).at(attrIndex).data(); } }

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
	virtual size_t	getNumberOfSatellites() { return satellites_m.size(); }
	virtual size_t	getNumberOfSatelliteMsmts() { return satelliteData_m.size(); }
	virtual double* getSatelliteDataPointers(int measurementInd, int satelliteInd, int attributeInd) { return satelliteData_m.at(measurementInd).at(satelliteInd).at(attributeInd).data(); }
	virtual void	writeSatelliteDataToCSV();
};//end class
#endif //end header guard