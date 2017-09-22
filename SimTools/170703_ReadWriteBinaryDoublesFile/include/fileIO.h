#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <fstream>
#include <iostream>
#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

namespace dblBinIO
{
	DLLEXPORT double* readDblBin(const char* filename, long numOfDblsToRead);
	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements);
	DLLEXPORT void clrDataMemory(double* dataArray);
}

#endif