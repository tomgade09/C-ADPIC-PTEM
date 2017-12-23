#include "FileIO\fileIO.h"
#include <iostream>
#include <string>

int main()
{
	bool cond{ true };
	double* attrs[3];
	double* vparainbin = new double[1024 * 1024];
	double* vperpinbin = new double[1024 * 1024];
	double* zinbin = new double[1024 * 1024];
	std::cout << "Reading Doubles\n";
	fileIO::readDblBin(vparainbin, "./../../../in/e_vpara.bin", 1024 * 1024);
	fileIO::readDblBin(vperpinbin, "./../../../in/e_vperp.bin", 1024 * 1024);
	fileIO::readDblBin(zinbin, "./../../../in/e_z.bin", 1024 * 1024);
	std::cout << "done\n";
	attrs[0] = vparainbin;
	attrs[1] = vperpinbin;
	attrs[2] = zinbin;

	while (cond)
	{
		int attrInd{ 0 };
		int partInd{ 0 };
		std::string instr;
		std::cout << "Quit(q) or Continue (any other): ";
		std::cin >> instr;
		//if (instr.find('q') != -1)
			//cond = false;
		//else
		instr.clear();
		//std::cout << instr.find('q');
		std::cout << "Input array to read from.  0-v_para, 1-v_perp, 2-z: ";
		std::cin >> attrInd;
		std::cout << "Input index of array to read: ";
		std::cin >> partInd;

		std::cout << "Attributes[" << attrInd << "][" << partInd << "] : " << attrs[attrInd][partInd] << "\n\n";
	}
	delete[] vparainbin;
	delete[] vperpinbin;
	delete[] zinbin;

	return 0;
}