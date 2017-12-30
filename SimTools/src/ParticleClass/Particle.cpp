#include "ParticleClass\Particle.h"

void Particle::loadFileToArray(std::vector<std::string> filenames, long numOfElements)
{
	//filenames check - make sure that length is equal to velDims + posDims
	if (filenames.size() > origData_m.size())
	{
		std::cout << "Warning: More filenames passed in than attributes of array.  Some files will not be loaded.";
		std::cout << "  Count of filenames: " << filenames.size() << " origData.size(): " << origData_m.size() << "\n";
	}
	else if (filenames.size() < origData_m.size())
	{
		std::cout << "Warning: Less filenames passed in than attributes of array.  Some dimensions will not have loaded data.";
		std::cout << "  Count of filenames: " << filenames.size() << " origData.size(): " << origData_m.size() << "\n";
	}

	//read particles into array
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
		fileIO::readDblBin(origData_m.at(attrs), filenames.at(attrs), particleCount_m);
}

void Particle::saveArrayToFile(std::vector<std::string> filenames, long numOfElements)
{//Should I ask for "labels" for the attributes: velocity and position, and base the savenames off of those??
	if (filenames.size() > origData_m.size())
	{
		std::cout << "Warning: More filenames passed in than attributes of array.  Some filenames will not be used.";
		std::cout << "  Count of filenames: " << filenames.size() << " origData.size(): " << origData_m.size() << "\n";
	}
	else if (filenames.size() < origData_m.size())
	{
		std::cout << "Warning: Less filenames passed in than attributes of array.  Creating generic names for the rest.";
		std::cout << "  Count of filenames: " << filenames.size() << " origData.size(): " << origData_m.size() << "\n";
		for (int iii = 1; iii < origData_m.size() - filenames.size() + 1; iii++)
			filenames.push_back(std::to_string(iii) + ".bin");
	}

	//write double binary files of the data
	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
		fileIO::writeDblBin(origData_m.at(attrs), filenames.at(attrs), particleCount_m);
}

void Particle::normalizeParticles(bool orig, bool curr, bool divide)
{
	if (!orig && !curr)
		return;
	
	if (normFactor_m == 1)
	{
		std::cout << "Warning: Norm factor is 1.  Normalizing will have no effect.  Returning.\n";
		return;
	}

	for (int attrs = 0; attrs < numberOfVelocityDims_m + numberOfPositionDims_m; attrs++)
	{
		for (int parts = 0; parts < particleCount_m; parts++)
		{
			if (orig) { origData_m.at(attrs).at(parts) /= (divide ? (normFactor_m) : (1 / normFactor_m)); } //maybe a little less performant than separating /= and *=
			if (curr) { currData_m.at(attrs).at(parts) /= (divide ? (normFactor_m) : (1 / normFactor_m)); }
		}
	}
}