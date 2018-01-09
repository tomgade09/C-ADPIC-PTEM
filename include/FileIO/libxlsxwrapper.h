#ifndef LIBXLSXWRAPPER_FILEIO_H
#define LIBXLSXWRAPPER_FILEIO_H
#include "xlsxwriter.h"
#include <vector>
#include <string>

namespace fileIO
{
	namespace xlsx
	{
		lxw_workbook* openNewWorkbook(std::string filename);
		lxw_workbook* writeDblArray(lxw_workbook* wb, std::string sheetname, std::vector<std::vector<double>> dataarray, std::vector<int> startColRow = { 0,0 });
		void closeAndSave(lxw_workbook* wb);
	}
}

#endif