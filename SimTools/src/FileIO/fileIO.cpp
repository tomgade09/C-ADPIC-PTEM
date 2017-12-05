#include "FileIO\fileIO.h"

namespace fileIO
{
	DLLEXPORT double* readDblBin(const char* filename, long numOfDblsToRead)
	{
		std::ifstream binFile{ filename, std::ios::binary };
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
		std::ifstream csv{ filename };
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

	DLLEXPORT void writeDblBin(const char* filename, double* dataarray, long numelements, bool overwrite)//overwrite defaults to true
	{
		std::ofstream binfile{ filename, std::ios::binary | (overwrite ? (std::ios::trunc) : (std::ios::app)) };
		if (!binfile.is_open())
		{
			std::cout << "Warning: Could not open (or create) file " << filename << " for writing!\n";
			return;
		}

		binfile.write(reinterpret_cast<const char*>(dataarray), std::streamsize(numelements * sizeof(double)));
		binfile.close();
	}

	DLLEXPORT void write2DCSV(const char* filename, double** dataarray, int numofentries, int numofcols, const char delim, bool overwrite)//overwrite defaults to true
	{
		std::ofstream csv(filename, overwrite ? (std::ios::trunc) : (std::ios::app));
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

	DLLEXPORT void writeTxtFile(const char* filename, const char* textToWrite, bool overwrite)//overwrite defaults to false
	{//could return ofstream, leave file open, or provide a boolean option, but I don't know how to return a null ofstream
		std::ofstream txt(filename, overwrite ? (std::ios::trunc) : (std::ios::app));
		if (!txt.is_open())
		{
			std::cout << "Could not open file: " << filename << "\n";
			return;
		}

		txt << textToWrite;
		txt.close();
	}
}