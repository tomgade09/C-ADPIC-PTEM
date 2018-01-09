#include "FileIO\fileIO.h"
#include <iostream>
#include <string>
#include <vector>

int main()
{
	bool cond{ true };
	long doublecnt;

	std::cout << "Input num of doubles to read: ";
	std::cin >> doublecnt;

	std::vector<std::vector<double>> attrs;
	attrs.resize(3);
	std::vector<double> vparainbin;
	std::vector<double> vperpinbin;
	std::vector<double> zinbin;
	vparainbin.resize(doublecnt);
	vperpinbin.resize(doublecnt);
	zinbin.resize(doublecnt);

	std::cout << "Reading " << doublecnt << " Doubles\n";
	fileIO::readDblBin(vparainbin, "./../../../in/ions_vpara.bin", doublecnt);
	fileIO::readDblBin(vperpinbin, "./../../../in/ions_vperp.bin", doublecnt);
	fileIO::readDblBin(zinbin, "./../../../in/ions_z.bin", doublecnt);
	std::cout << "done\n";
	attrs.at(0) = vparainbin;
	attrs.at(1) = vperpinbin;
	attrs.at(2) = zinbin;

	while (cond)
	{
		int attrInd{ 0 };
		int partInd{ 0 };
		std::string instr;
		std::cout << "Input index of array to read, -1 = quit: ";
		std::cin >> partInd;
		if (partInd == -1)
			return 0;

		std::cout << "\nParticle[" << partInd << "] : ";
		for (int attribute = 0; attribute < 3; attribute++)
			std::cout << attrs[attribute][partInd] << "  ";
		std::cout << "\n\n";
	}

	return 0;
}