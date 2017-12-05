#include "StandaloneTools\binaryfiletools.h"

void saveParticleAttributeToDisk(double* arrayToSave, int length, const char* foldername, const char* name)
{
	std::string fn{ foldername };
	fn = fn + name + ".bin";
	fileIO::writeDblBin(fn.c_str(), arrayToSave, length);
}

void loadFileIntoParticleAttribute(double* arrayToLoadInto, int length, const char* foldername, const char* name)
{
	delete[] arrayToLoadInto;
	std::string fn{ foldername };
	fn = fn + name;
	arrayToLoadInto = fileIO::readDblBin(fn.c_str(), length);
}