#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
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
	DLLEXPORT void readDblBin(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead);
	DLLEXPORT double** read2DCSV(std::string filename, int numofentries, int numofcols, const char delim);
	DLLEXPORT void writeDblBin(std::vector<double> dataarray, std::string filename, long numelements, bool overwrite=true);
	DLLEXPORT void write2DCSV(std::vector<std::vector<double>> dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite=true, int precision=20);
	DLLEXPORT void writeTxtFile(const char* filename, const char* textToWrite, bool overwrite=false);
}

#endif