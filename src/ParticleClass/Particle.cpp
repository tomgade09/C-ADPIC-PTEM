#include "ParticleClass\Particle.h"

void Particle::loadDataFromDisk(std::string folder, bool orig)
{
	FILE_RDWR_EXCEP_CHECK(
		for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
			fileIO::readDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m);
	);

	initDataLoaded_m = true;
}

void Particle::saveDataToDisk(std::string folder, bool orig)
{
	FILE_RDWR_EXCEP_CHECK(
		for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
			fileIO::writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m);
	);
}

int Particle::getDimensionIndByName(std::string searchName)
{
	for (int name = 0; name < attributeNames_m.size(); name++)
	{
		if (searchName == attributeNames_m.at(name))
			return name;
	}
	
	throw std::invalid_argument("Particle::getDimensionIndByName: specified name is not present in name array: " + searchName);
}

std::string Particle::getDimensionNameByInd(int searchIndx)
{
	if (!(searchIndx <= (attributeNames_m.size() - 1) && (searchIndx >= 0)))
		throw std::invalid_argument("Particle::getDimensionNameByInd: specified index is invalid: " + std::to_string(searchIndx));

	return attributeNames_m.at(searchIndx);
}