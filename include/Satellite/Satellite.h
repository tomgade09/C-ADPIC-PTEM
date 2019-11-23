#ifndef SATELLITE_H
#define SATELLITE_H

#include <vector>
#include <string>

#include "dlldefines.h"

using std::vector;
using std::string;

#define STRVEC vector<string>
#define DBLVEC vector<double>
#define DBL2DV vector<vector<double>>
#define DBL3DV vector<vector<vector<double>>>

class Satellite
{
protected:
	string name_m;
	STRVEC attributeNames_m;
	
	double altitude_m{ 0.0 };
	bool   upwardFacing_m{ false };
	bool   initializedGPU_m{ false };

	long numberOfParticles_m{ -1 };

	DBL2DV   data_m; //[attribute][particle]
	double*  satCaptrData1D_d{ nullptr }; //flattened satellite capture data on GPU
	double** satCaptrData2D_d{ nullptr }; //2D satellite capture data on GPU
	double** particleData2D_d{ nullptr };

	virtual void   initializeGPU();
	virtual void   freeGPUMemory();
	virtual void   deserialize(string serialFolder, string name, double** particleData2D);
	virtual size_t getAttrIndByName(string name);

public:
	Satellite(string name, STRVEC attributeNames, double altitude, bool upwardFacing, long numberOfParticles, double** partDataGPUPtr);
	Satellite(string serialFolder, string name, double** particleData2D);
	virtual ~Satellite();
	Satellite(const Satellite&) = delete;
	Satellite& operator=(const Satellite&) = delete;

	//Access functions
	string        name() const;
	double        altitude() const;
	bool	      upward() const;
	DBL2DV&       __data();
	const DBL2DV& data() const;
	double**      get2DDataGPUPtr() const;
	double*       get1DDataGPUPtr() const;
	size_t  getNumberOfAttributes() const;
	long    getNumberOfParticles()  const;

	//Other functions
	virtual void iterateDetector(double simtime, double dt, int blockSize); //increment time, track overall sim time, or take an argument??
	virtual void iterateDetectorCPU(const DBL2DV& particleData, double simtime, double dt);
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??
	virtual void saveDataToDisk(string folder);
	virtual void loadDataFromDisk(string folder);

	virtual DBL2DV removeZerosData();
	virtual void serialize(string serialFolder);
};

#undef STRVEC
#undef DBLVEC
#undef DBL2DV
#undef DBL3DV

#endif