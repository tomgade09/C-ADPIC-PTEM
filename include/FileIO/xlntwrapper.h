#ifndef XLNTWRAPPER_FILEIO_H
#define XLNTWRAPPER_FILEIO_H
#include "xlnt\xlnt.hpp"
#include <vector>
#include <string>

//could never get working, so excluding from project...leaving around in case I come back in the future

namespace fileIO
{
	void writeXLSXDblArray(std::string filename, std::vector<std::vector<double>> dataarray, std::string sheetname = "Sheet1", std::vector<int> startColRow = { 1, 1 }, bool newfile = true);
}

#endif