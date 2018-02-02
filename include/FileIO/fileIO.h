#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace fileIO
{
	void readDblBin(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead);
	void read2DCSV(std::vector<std::vector<double>>& array2DToReadInto, std::string filename, int numofentries, int numofcols, const char delim);
	void writeDblBin(const std::vector<double>& dataarray, std::string filename, long numelements, bool overwrite=true);
	void write2DCSV(const std::vector<std::vector<double>>& dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite=true, int precision=20);
	void writeTxtFile(std::string textToWrite, std::string filename, bool overwrite=false);
}

#endif