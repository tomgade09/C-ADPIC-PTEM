#include "StandaloneTools\binaryfiletools.h"

void saveParticleAttributeToDisk(double* arrayToSave, int length, const char* foldername, const char* name)
{
	std::string fn{ foldername };
	fn = fn + name;
	fileIO::writeDblBin(fn.c_str(), arrayToSave, length);
}

void loadFileIntoParticleAttribute(double* arrayToLoadInto, int length, const char* foldername, const char* name)
{
	std::string fn{ foldername };
	fn = fn + name;
	fileIO::readDblBin(arrayToLoadInto, fn.c_str(), length);
}

void stringPadder(std::string& in, int totalStrLen, int indEraseFrom)
{
	if (totalStrLen <= 0 || indEraseFrom <= 0)
		return;

	size_t txtlen = in.length();

	if ((totalStrLen - txtlen) > 0)
	{
		for (int iii = 0; iii < (totalStrLen - txtlen); iii++)
			in += ' ';
	}
	else
		in.erase(indEraseFrom, txtlen - totalStrLen);
}

double*** form3Darray(int d1, int d2, int d3)
{
	double*** array3D;
	array3D = new double**[d1];
	for (int iii = 0; iii < d1; iii++)
	{
		array3D[iii] = new double*[d2];
		array3D[iii][0] = new double[d2 * d3];
		for (int jjj = 0; jjj < d2; jjj++)
		{
			array3D[iii][jjj] = array3D[iii][0] + jjj * d3;
			for (int kk = 0; kk < d3; kk++)//values to initialize array
				array3D[iii][jjj][kk] = 0.0;
		}
	}
	return array3D;
}

void delete3Darray(double*** array3D, int d1)
{
	for (int iii = 0; iii < d1; iii++)
	{
		delete[] array3D[iii][0];
		delete[] array3D[iii];
	}
	delete[] array3D;
}