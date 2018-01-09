#include "FileIO\libxlsxwrapper.h"

namespace fileIO
{
	namespace xlsx
	{
		lxw_workbook* openNewWorkbook(std::string filename)
		{
			return workbook_new(filename.c_str());
		}

		lxw_workbook* writeDblArray(lxw_workbook* wb, std::string sheetname, std::vector<std::vector<double>> dataarray, std::vector<int> startColRow)
		{
			lxw_worksheet* ws = workbook_get_worksheet_by_name(wb, sheetname.c_str());

			if (ws == NULL)
				workbook_add_worksheet(wb, sheetname.c_str());

			for (int col = 0; col < dataarray.size(); col++)
			{
				for (int row = 0; row < dataarray.at(col).size(); row++)
					worksheet_write_number(ws, startColRow.at(1) + row, startColRow.at(0) + col, dataarray.at(col).at(row), NULL);
			}

			return wb;
		}

		void closeAndSave(lxw_workbook* wb)
		{
			workbook_close(wb);
		}
	}//end fileIO::xlsx namespace
}//end fileIO namesapce