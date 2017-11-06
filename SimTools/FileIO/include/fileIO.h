#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

namespace fileIO
{
	DLLEXPORT double* readDblBin(const char* filename, long numOfDblsToRead);
	DLLEXPORT double** read2DCSV(const char* filename, int numofentries, int numofcols, const char delim);
	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements);
	DLLEXPORT void write2DCSV(const char* filename, double** dataarray, int numofentries, int numofcols, const char delim);
	DLLEXPORT void clrDataMemory(double* dataArray);
}

#endif