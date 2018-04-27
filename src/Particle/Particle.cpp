#include "Particle\Particle.h"

#include <iostream>
#include "utils\fileIO.h"

using utils::fileIO::readDblBin;
using utils::fileIO::writeDblBin;

//need a custom solution for this...
//file read/write exception checking (probably should mostly wrap fileIO functions)
#define FILE_RDWR_EXCEP_CHECK(x) \
	try{ x; } \
	catch(const std::invalid_argument& a) { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Invalid argument error: " << a.what() << ": continuing without loading file" << std::endl; std::cout << "FileIO exception: check log file for details" << std::endl; } \
	catch(...)                            { throw; }


void Particle::loadDataFromDisk(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(readDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));

	initDataLoaded_m = true;
}

void Particle::saveDataToDisk(std::string folder, bool orig) const
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));
}

int Particle::getAttrIndByName(std::string searchName) const
{
	for (int name = 0; name < attributeNames_m.size(); name++)
		if (searchName == attributeNames_m.at(name))
			return name;
	
	throw std::invalid_argument("Particle::getDimensionIndByName: specified name is not present in name array: " + searchName);
}

std::string Particle::getAttrNameByInd(int searchIndx) const
{
	if (!(searchIndx <= (attributeNames_m.size() - 1) && (searchIndx >= 0)))
		throw std::invalid_argument("Particle::getDimensionNameByInd: specified index is invalid: " + std::to_string(searchIndx));

	return attributeNames_m.at(searchIndx);
}