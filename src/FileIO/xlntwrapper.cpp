#include "FileIO\xlntwrapper.h"

//could never get working, so excluding from project...leaving around in case I come back in the future

namespace fileIO
{
	void writeXLSXDblArray(std::string filename, std::vector<std::vector<double>> dataarray, std::string sheetname, std::vector<int> startColRow, bool newfile)
	{
		xlnt::workbook wb; //workbook.cpp
		if (!newfile)
			wb.load(filename); //workbook.cpp
		
		xlnt::worksheet ws; //worksheet.cpp
		try
		{
			ws = wb.sheet_by_title(sheetname); //workbook.cpp
		}
		catch (...)
		{
			ws = wb.create_sheet(); //workbook.cpp
			ws.title(sheetname); //worksheet.cpp
		}

		for (int col = 0; col < dataarray.size(); col++)
		{
			for (int row = 0; row < dataarray.at(col).size(); row++)
				ws.cell(xlnt::cell_reference(col + startColRow.at(0), row + startColRow.at(1))).value(dataarray.at(col).at(row));
		}//worksheet.cpp
		
		wb.save(filename);
	}
}