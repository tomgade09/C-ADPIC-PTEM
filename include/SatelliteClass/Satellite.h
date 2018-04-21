#ifndef SATELLITE_H
#define SATELLITE_H

#include <vector>
#include <string>

class Satellite
{
protected:
	double altitude_m;
	bool upwardFacing_m;
	bool dataReady_m{ false };
	bool dataOnGPU_m{ false };

	long numberOfParticles_m;

	std::vector<std::vector<std::vector<double>>> data_m; //[measurement][attribute][particle]
	double*  satCaptrData1D_d{ nullptr }; //flattened satellite capture data
	double** satCaptrData2D_d{ nullptr }; //2D satellite capture data
	double** particleData2D_d{ nullptr };
	
	std::vector<std::string> attrNames_m;
	std::string name_m;

	virtual void initializeGPU();
	virtual void dataAllocateNewMsmtVector() { data_m.push_back(std::vector<std::vector<double>>(attrNames_m.size(), std::vector<double>(numberOfParticles_m))); }

public:
	Satellite(std::string name, std::vector<std::string> attributeNames, double altitude, bool upwardFacing, long numberOfParticles, double** partDataGPUPtr):
		name_m{ name }, attrNames_m{ attributeNames }, altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfParticles_m{ numberOfParticles }, particleData2D_d{ partDataGPUPtr }
	{ initializeGPU(); }
	
	virtual ~Satellite() { freeGPUMemory(); }
	
	virtual void iterateDetector(double simtime, double dt, int blockSize); //increment time, track overall sim time, or take an argument??
	virtual void copyDataToHost(); //some sort of sim time check to verify I have iterated for the current sim time??
	virtual void freeGPUMemory();
	virtual void saveDataToDisk(std::string folder);
	virtual void loadDataFromDisk(std::string folder, bool expand);

	//Access functions
	std::string name() { return name_m; }
	double      altitude() { return altitude_m; }
	bool	    upward() { return upwardFacing_m; }
	std::vector<std::vector<double>> getConsolidatedData(bool removeZeros);
	const std::vector<std::vector<std::vector<double>>>& data() { return data_m; }
	double** get2DDataGPUPtr() { return satCaptrData2D_d; }
	double*  get1DDataGPUPtr() { return satCaptrData1D_d; }

	int     getNumberOfAttributes() { return (int)attrNames_m.size(); }
	long    getNumberOfParticles()  { return numberOfParticles_m; }
};
#endif