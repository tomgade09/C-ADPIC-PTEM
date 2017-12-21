#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

namespace fileIO
{
	DLLEXPORT void readDblBin(double* arrayToReadInto, const char* filename, long numOfDblsToRead);
	DLLEXPORT double** read2DCSV(const char* filename, int numofentries, int numofcols, const char delim);
	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements, bool overwrite=true);
	DLLEXPORT void write2DCSV(const char* filename, double** dataarray, int numofentries, int numofcols, const char delim, bool overwrite=true, int precision = 20);
	DLLEXPORT void writeTxtFile(const char* filename, const char* textToWrite, bool overwrite=false);
}

#endif