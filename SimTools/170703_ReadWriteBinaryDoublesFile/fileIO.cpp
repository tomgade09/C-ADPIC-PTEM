#include "include\fileIO.h"
#define DEBUGSTDOUTPUT

namespace dblBinIO
{
	DLLEXPORT double* readDblBin(const char* filename, long numOfDblsToRead)
	{
		std::cout << "read: " << filename << "\n";
		std::ifstream binFile;
		binFile.open(filename, std::ios::in | std::ios::binary);

		if (!binFile.is_open())
		{
			std::cout << "Warning: Could not open file " << filename << " for reading!\n";
			return nullptr;
		}

		double* dataArray = new double[numOfDblsToRead];
		binFile.read(reinterpret_cast<char*>(dataArray), std::streamsize(numOfDblsToRead * sizeof(double)));

		binFile.close();

		return dataArray;
	}

	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements)
	{
		std::cout << "write: " << filename << "\n";
		std::ofstream binfile;
		binfile.open(filename, std::ios::binary | std::ios::out);

		if (!binfile.is_open())
		{
			std::cout << "Warning: Could not open (or create) file " << filename << " for writing!\n";
			return;
		}

		binfile.write(reinterpret_cast<const char*>(dataarray), std::streamsize(numelements * sizeof(double)));
		binfile.close();
	}

	DLLEXPORT void clrDataMemory(double* dataArray)
	{
		delete[] dataArray;
	}
}