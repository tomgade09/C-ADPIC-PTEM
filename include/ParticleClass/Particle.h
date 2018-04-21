#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>

class Particle
{
protected:
	std::vector<std::vector<double>> origData_m;
	std::vector<std::vector<double>> currData_m;

	double*  origData1D_d{ nullptr };
	double*  currData1D_d{ nullptr };
	double** origData2D_d{ nullptr };
	double** currData2D_d{ nullptr };

	std::vector<std::string> attributeNames_m;

	long numberOfParticles_m;
	double mass_m;
	double charge_m;

	std::string name_m;

	void initializeGPU();

	bool initDataLoaded_m{ false };
	bool dataOnGPU_m{ false };

public:
	Particle(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts) :
		name_m{ name }, attributeNames_m{ attributeNames }, mass_m{ mass }, charge_m{ charge }, numberOfParticles_m{ numParts }
	{
		origData_m = std::vector<std::vector<double>>(attributeNames.size(), std::vector<double>(numParts));
		currData_m = std::vector<std::vector<double>>(attributeNames.size(), std::vector<double>(numParts));

		initializeGPU();
	}
	~Particle()
	{ freeGPUMemory(); }

	std::vector<std::vector<double>>& getOrigData() { return origData_m; }
	std::vector<std::vector<double>>& getCurrData() { return currData_m; }
	std::vector<std::string>& getAttrNames() { return attributeNames_m; }
	std::string name() { return name_m; }
	double   mass() { return mass_m; }
	double   charge() { return charge_m; }
	long     getNumberOfParticles() { return numberOfParticles_m; }
	int      getNumberOfAttributes() { return (int)attributeNames_m.size(); }
	bool     getInitDataLoaded() { return initDataLoaded_m; }
	double** getOrigDataGPUPtr() { return origData2D_d; }
	double** getCurrDataGPUPtr() { return currData2D_d; }

	int         getAttrIndByName(std::string searchName);
	std::string getAttrNameByInd(int searchIndx);

	void loadDataFromMem(std::vector<std::vector<double>> data, bool orig = true) { ((orig) ? origData_m = data : currData_m = data); numberOfParticles_m = ((orig) ? (int)origData_m.at(0).size() : (int)currData_m.at(0).size()); }
	void loadDataFromDisk(std::string folder, bool orig = true);
	void saveDataToDisk(std::string folder, bool orig);
	void generateRandomParticles(const std::vector<double>& s, int startInd, int length, double vmean, double kBT_eV, double mass);
	
	void copyDataToGPU();
	void copyDataToHost();
	void freeGPUMemory();
	void clearGPUMemory();
};

#endif