#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>

#include "dlldefines.h"

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
	Particle(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts);
	~Particle();

	const std::vector<std::vector<double>>& data(bool orig) const { return ((orig) ? origData_m : currData_m); }
	std::vector<double>& dataAttr(bool orig, int attr) { return ((orig) ? origData_m.at(attr) : currData_m.at(attr)); }
	const std::vector<std::string>& attrNames() const { return attributeNames_m; }
	std::string name()    const { return name_m; }
	double      mass()    const { return mass_m; }
	double      charge()  const { return charge_m; }
	long        getNumberOfParticles()  const { return numberOfParticles_m; }
	int         getNumberOfAttributes() const { return (int)attributeNames_m.size(); }
	bool        getInitDataLoaded() const { return initDataLoaded_m; }
	double**    getOrigDataGPUPtr() const { return origData2D_d; }
	double**    getCurrDataGPUPtr() const { return currData2D_d; }

	int         getAttrIndByName(std::string searchName) const;
	std::string getAttrNameByInd(int searchIndx) const;

	void loadDataFromMem(std::vector<std::vector<double>> data, bool orig = true) { ((orig) ? origData_m = data : currData_m = data); numberOfParticles_m = ((orig) ? (int)origData_m.at(0).size() : (int)currData_m.at(0).size()); }
	void loadDataFromDisk(std::string folder, bool orig = true);
	void saveDataToDisk(std::string folder, bool orig) const;
	void generateRandomParticles(const std::vector<double>& s, int startInd, int length, double vmean, double kBT_eV, double mass);
	
	void copyDataToGPU();
	void copyDataToHost();
	void freeGPUMemory();
	void clearGPUMemory();
};

#endif