#ifndef FILEIOAPI_H
#define FILEIOAPI_H

#include <vector>
#include <string>

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif /* DLLFILE */

//write these at some point

DLLEXPORT void readDblBinAPI(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead);
DLLEXPORT void read2DCSVAPI(std::vector<std::vector<double>>& array2DToReadInto, std::string filename, int numofentries, int numofcols, const char delim);
DLLEXPORT void writeDblBinAPI(const std::vector<double>& dataarray, std::string filename, long numelements, bool overwrite = true);
DLLEXPORT void write2DCSVAPI(const std::vector<std::vector<double>>& dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite = true, int precision = 20);
DLLEXPORT void writeTxtFileAPI(std::string filename, std::string textToWrite, bool overwrite = false);

#endif /* !FILEIOAPI_H */