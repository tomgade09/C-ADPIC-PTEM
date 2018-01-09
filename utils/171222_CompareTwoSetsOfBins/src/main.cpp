#include "FileIO\fileIO.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

int main()
{
	std::vector<std::vector<std::vector<double>>> alldata;

	std::vector<std::string> filenames{ "e_vpara.bin", "e_vperp.bin", "e_z.bin", "i_vpara.bin", "i_vperp.bin", "i_z.bin" };
	std::string folder{ "./../../../in/" };

	long numberOfParts{ 0 };
	std::cout << "How many particles per file?? : ";
	std::cin >> numberOfParts;

	for (int dataind = 0; dataind < 2; dataind++)
	{
		std::vector<std::vector<double>> data;
		for (int fileind = 0; fileind < 6; fileind++)
		{
			std::vector<double> bins;
			bins.resize(numberOfParts);
			fileIO::readDblBin(bins, folder + ((dataind == 0) ? "data1/" : "data2/") + filenames[fileind], numberOfParts);
			data.push_back(bins);
		}
		alldata.push_back(data);
	}

	std::string tmp;
	long wrong{ 0 };
	for (int attrind = 0; attrind < 6; attrind++)
	{
		for (int partind = 0; partind < numberOfParts; partind++)
		{
			if (alldata.at(0).at(attrind).at(partind) != alldata.at(1).at(attrind).at(partind))
			{
				wrong++;
				if (abs((alldata.at(0).at(attrind).at(partind) - alldata.at(1).at(attrind).at(partind)) / alldata.at(0).at(attrind).at(partind)) > 1e-10) { std::cout << "Large Error: " << abs((alldata.at(0).at(attrind).at(partind) - alldata.at(1).at(attrind).at(partind)) / alldata.at(0).at(attrind).at(partind)) << "\n"; }
				if (wrong % 1000 == 0) { std::cout << std::setprecision(20) << "1000 wrongs: " << alldata.at(0).at(attrind).at(partind) << "\n             " << alldata.at(1).at(attrind).at(partind) << "\n"; }
			}
		}
		std::cout << "Wrong total: " << wrong << " Press enter to continue.\n";
		std::cin >> tmp;
		wrong = 0;
	}

	return 0;
}