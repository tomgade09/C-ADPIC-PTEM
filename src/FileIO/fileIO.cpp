#include "FileIO\fileIO.h"

namespace fileIO
{
	DLLEXPORT void readDblBin(std::vector<double>& arrayToReadInto, std::string filename, long numOfDblsToRead)
	{
		std::ifstream binFile{ filename, std::ios::binary };
		if (!binFile.is_open())
			throw std::invalid_argument ("fileIO::readDblBin: could not open file " + filename + " for reading");
		if (arrayToReadInto.size() < numOfDblsToRead)
		{
			binFile.close();
			throw std::invalid_argument ("fileIO::readDblBin: std::vector is not big enough to contain the data being read from file " + filename);
		}
		
		binFile.seekg(0, binFile.end);
		int length{ static_cast<int>(binFile.tellg()) };
		binFile.seekg(0, binFile.beg);

		if (length < numOfDblsToRead * 8)
		{
			binFile.close();
			throw std::invalid_argument ("fileIO::readDblBin: filesize of \"" + filename + "\" is smaller than specified number of doubles to read");
		}
		if (length > numOfDblsToRead * 8)
			std::cerr << "fileIO::readDblBin: warning: size of data read is less than the size of all data in file " << filename << ": continuing" << std::endl;

		binFile.read(reinterpret_cast<char*>(arrayToReadInto.data()), std::streamsize(numOfDblsToRead * sizeof(double)));
		binFile.close();
	}

	DLLEXPORT double** read2DCSV(std::string filename, int numofentries, int numofcols, const char delim)
	{
		std::ifstream csv{ filename };
		if (!csv.is_open())
			throw std::invalid_argument ("fileIO::read2DCSV: could not open file " + filename + " for reading");

		double** ret = new double*[numofcols];
		double* inner = new double[numofentries * numofcols];

		//1D array converted to 2D by using pointers to location of the start of the next column:
		// idx:  0				1			2		 "numofentries"  "numofe + 1"	"2*numofe"
		//[col 0 elem 0, col 0 elem 1, col 0 elem 2...col 1 elem 0, col 1 elem 1...col 2 elem 0...]
		//ptrs: [0]										[1]							[2]
		//So - 2nd dimension of array looks like this:
		//[pointer to 0, pointer to numofentries, pointer to 2 * numofentries...]
		for (int iii = 0; iii < numofcols; iii++)
			ret[iii] = &inner[iii * numofentries];

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
		csv.close();
		
		return ret;
	}

	DLLEXPORT void writeDblBin(std::vector<double> dataarray, std::string filename, long numelements, bool overwrite)//overwrite defaults to true
	{
		std::ofstream binfile{ filename, std::ios::binary | (overwrite ? (std::ios::trunc) : (std::ios::app)) };
		if (!binfile.is_open())
			throw std::invalid_argument ("fileIO::writeDblBin: could not open file " + filename + " for writing");
		if (dataarray.size() < numelements)
		{
			binfile.close();
			throw std::invalid_argument ("fileIO::writeDblBin: size of data vector is less than the number of doubles requested from it for filename " + filename);
		}

		binfile.write(reinterpret_cast<char*>(dataarray.data()), std::streamsize(numelements * sizeof(double)));
		binfile.close();
	}

	DLLEXPORT void write2DCSV(std::vector<std::vector<double>> dataarray, std::string filename, int numofentries, int numofcols, const char delim, bool overwrite, int precision)//overwrite defaults to true, precision to 20
	{
		std::ofstream csv(filename, overwrite ? (std::ios::trunc) : (std::ios::app));
		if (!csv.is_open())
			throw std::invalid_argument ("fileIO::write2DCSV: could not open file " + filename + " for writing");
		if (dataarray.size() < numofcols)
		{
			csv.close();
			throw std::invalid_argument ("fileIO::write2DCSV: size of data vector is less than the doubles requested from it for filename " + filename);
		}

		for (int iii = 0; iii < numofentries; iii++)
		{
			for (int jjj = 0; jjj < numofcols; jjj++)
				csv << std::setprecision(precision) << dataarray.at(jjj).at(iii) << delim;
			
			csv << "\n";
		}

		csv.close();
		
		return;
	}

	DLLEXPORT void writeTxtFile(std::string filename, std::string textToWrite, bool overwrite)//overwrite defaults to false
	{//could return ofstream, leave file open, or provide a boolean option, but I don't know how to return a null ofstream
		std::ofstream txt(filename, overwrite ? (std::ios::trunc) : (std::ios::app));
		if (!txt.is_open())
			throw std::invalid_argument ("fileIO::writeTxtFile: could not open file " + filename + " for writing");

		txt << textToWrite;
		txt.close();
	}
}