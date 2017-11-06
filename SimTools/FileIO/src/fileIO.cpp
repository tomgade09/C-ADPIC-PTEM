#include "include\fileIO.h"

namespace fileIO
{
	DLLEXPORT double* readDblBin(const char* filename, long numOfDblsToRead)
	{
		//std::cout << "read: " << filename << "\n";
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

	DLLEXPORT double** read2DCSV(const char* filename, int numofentries, int numofcols, const char delim)
	{
		std::ifstream csv(filename);
		if (!csv.is_open())
		{
			std::cout << "Could not open file: " << filename << "\n";
			return nullptr;
		}

		double** ret = new double*[numofcols];

		for (int iii = 0; iii < numofcols; iii++)
		{
			ret[iii] = new double[numofentries];
		}

		for (int iii = 0; iii < numofentries; iii++)
		{
			std::string in;
			std::getline(csv, in);

			std::stringstream in_ss(in);

			for (int jjj = 0; jjj < numofcols; jjj++)
			{
				std::string val;
				std::getline(in_ss, val, delim);
				std::stringstream convert(val);
				convert >> ret[jjj][iii];
			}
		}

		return ret;
	}

	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements)
	{
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

	DLLEXPORT void write2DCSV(const char* filename, double** dataarray, int numofentries, int numofcols, const char delim)
	{
		std::ofstream csv(filename);
		if (!csv.is_open())
		{
			std::cout << "Could not open file: " << filename << "\n";
			return;
		}

		for (int iii = 0; iii < numofentries; iii++)
		{
			for (int jjj = 0; jjj < numofcols; jjj++)
				csv << dataarray[jjj][iii] << delim;
			
			csv << "\n";
		}
		
		return;
	}

	DLLEXPORT void clrDataMemory(double* dataArray)
	{
		delete[] dataArray;
	}
}