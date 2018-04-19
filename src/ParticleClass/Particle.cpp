#include "ParticleClass\Particle.h"

void Particle::loadDataFromDisk(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));

	initDataLoaded_m = true;
}

void Particle::saveDataToDisk(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(fileIO::writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));
}

int Particle::getAttrIndByName(std::string searchName)
{
	for (int name = 0; name < attributeNames_m.size(); name++)
		if (searchName == attributeNames_m.at(name))
			return name;
	
	throw std::invalid_argument("Particle::getDimensionIndByName: specified name is not present in name array: " + searchName);
}

std::string Particle::getAttrNameByInd(int searchIndx)
{
	if (!(searchIndx <= (attributeNames_m.size() - 1) && (searchIndx >= 0)))
		throw std::invalid_argument("Particle::getDimensionNameByInd: specified index is invalid: " + std::to_string(searchIndx));

	return attributeNames_m.at(searchIndx);
}