#include "include\Simulation170925.h"

void Simulation170925::setElecMagLUT(const char* filename, int rows, int cols)
{
	 elcFieldLUT_m = fileIO::read2DCSV("ez.out", rows, cols, ' ');
}

double Simulation170925::calculateBFieldAtZandTime(double z, double time)
{
	return BFieldatZ(z, time);
}

double Simulation170925::calculateEFieldAtZandTime(double z, double time)
{
	return EFieldatZ(getPointerToElectricFieldData(), z, time);
}

double* Simulation170925::returnResults()
{
	std::string fold;
	fold = "./particles_final/";
	std::vector<std::string> names;
	names.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m);
	names = { "v_e_para", "mu_e", "z_e", "v_i_para", "mu_i", "z_i" };
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
	}

	serializeParticleArray();
	return getPointerToSerializedParticleArray();
}