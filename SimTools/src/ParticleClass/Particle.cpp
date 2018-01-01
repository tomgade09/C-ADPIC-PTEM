#include "ParticleClass\Particle.h"

void Particle::loadFilesToArray(std::string folder)
{
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
	{
		fileIO::readDblBin(currData_m.at(attrs), folder + name_m + "_" + attributeNames_m.at(attrs) + ".bin", particleCount_m);
		//
		//
		//Need to remove this
		for (int parts = 0; parts < particleCount_m; parts++)
			currData_m.at(attrs).at(parts) *= 6.371e6;
		//Need to remove this
		//
		//
	}

	initDataLoaded_m = true;
}

void Particle::saveArrayToFiles(std::string folder, bool orig)
{
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
		fileIO::writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + name_m + "_" + attributeNames_m.at(attrs) + ".bin", particleCount_m);
}

void Particle::normalizeParticles(bool orig, bool curr, bool inverse)
{
	if (!orig && !curr)
		return;
	
	if (normalized_m == true)
	{
		std::cout << "Warning: At least one of the data sets is already normalized.  Make sure neither data set is being normalized twice!\n";
	}

	if (normFactor_m == 1)
	{
		std::cout << "Warning: Norm factor is 1.  Normalizing will have no effect.  Returning.\n";
		return;
	}

	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
	{
		for (int parts = 0; parts < particleCount_m; parts++)
		{//normalize -> divide by normalization factor
			if (orig) { origData_m.at(attrs).at(parts) /= (inverse ? (1 / normFactor_m) : (normFactor_m)); } //maybe a little less performant than separating /= and *=
			if (curr) { currData_m.at(attrs).at(parts) /= (inverse ? (1 / normFactor_m) : (normFactor_m)); }
		}
	}

	normalized_m = true;
}

int Particle::getDimensionIndByName(std::string searchName)
{
	for (int name = 0; name < attributeNames_m.size(); name++)
	{
		if (searchName == attributeNames_m.at(name))
			return name;
	}

	return -1;
}