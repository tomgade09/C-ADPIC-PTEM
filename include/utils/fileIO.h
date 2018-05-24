#ifndef UTILS_FILEIO_H
#define UTILS_FILEIO_H

#include <vector>
#include <string>

#include "dlldefines.h"
#include "utils/readIOclasses.h"
#include "utils/writeIOclasses.h"

namespace utils
{
	namespace fileIO
	{
		DLLEXP void readDblBin(std::vector<double>& arrayToReadInto, std::string filename);
		DLLEXP void readDblBin(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead);
		DLLEXP void read2DCSV(std::vector<std::vector<double>>& array2DToReadInto, std::string filename, int numofentries, int numofcols, const char delim);
		DLLEXP void readTxtFile(std::string& readInto, std::string filename);
		DLLEXP void writeDblBin(const std::vector<double>& dataarray, std::string filename, long numelements, bool overwrite = true);
		DLLEXP void write2DCSV(const std::vector<std::vector<double>>& dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite = true, int precision = 20);
		DLLEXP void writeTxtFile(std::string textToWrite, std::string filename, bool overwrite = false);
	}
}
#endif /* !UTILS_FILEIO_H */