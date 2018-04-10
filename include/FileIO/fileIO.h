#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

//file read/write exception checking (probably should mostly wrap fileIO functions)
#define FILE_RDWR_EXCEP_CHECK(x) \
	try{ x; } \
	catch(const std::invalid_argument& a) { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Invalid argument error: " << a.what() << ": continuing without loading file" << std::endl; std::cout << "Exception: check log for details" << std::endl; } \
	catch(...)                            { throw; }

namespace fileIO
{
	void readDblBin(std::vector<double>& arrayToReadInto, std::string filename);
	void readDblBin(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead);
	void read2DCSV(std::vector<std::vector<double>>& array2DToReadInto, std::string filename, int numofentries, int numofcols, const char delim);
	void readTxtFile(std::string& readInto, std::string filename);
	void writeDblBin(const std::vector<double>& dataarray, std::string filename, long numelements, bool overwrite=true);
	void write2DCSV(const std::vector<std::vector<double>>& dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite=true, int precision=20);
	void writeTxtFile(std::string textToWrite, std::string filename, bool overwrite=false);
	void writeAttrsToFiles(std::vector<double> chars, std::vector<std::string> charNames, std::string className, std::string saveFolder);
}

#endif