//an encapsulant for all the vital elements of a simulation
//should cut down a bit on the complexity of writing new sims

#include "include\simulation.h"
#include "include\fileIO.h"

void Simulation::generateNormallyDistributedValues(int numberOfNormalAttributesPerPart, double* means, double* sigmas)
{
	if (numberOfNormalAttributesPerPart > numberOfAttributesTracked_m)
	{
		std::cout << "Error: More attributes specified for Simulation::generateNormallyDistributedValues than originally specified ";
		std::cout << "as 'total number of attributes tracked per particle'.  You WILL get an index error.  Returning without changes.";
		return;
	}

	if (initialized_m)
	{
		std::cout << "Warning: sim has been initialized.  Anything currently exisiting in particles_m will be overwritten!!\n";
	}
	
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfNormalAttributesPerPart; jjj++)
		{
			std::random_device randDev;  //means and sigmas need to be passed in as: [particle1-attr1 mean/sigma, attr2 mean/sigma...
			std::mt19937 mtgen(randDev());											//particle2-attr1 mean/sigma, attr2 mean/sigma...]

			std::normal_distribution<> data_nd(means[iii * numberOfNormalAttributesPerPart + jjj], sigmas[iii * numberOfNormalAttributesPerPart + jjj]);

			for (int kk = 0; kk < numberOfParticlesPerType_m; kk++)
				particles_m[iii][jjj][kk] = data_nd(mtgen);
		}//end for jjj
	}//end for iii
}//end function

double Simulation::calculateMeanOfParticleAttribute(int particleIndex, int attributeIndex, bool absValue)
{	
	double sum{ 0 };
	for (int iii = 0; iii < numberOfParticlesPerType_m; iii++)
	{
		if (absValue)
			sum += abs(particles_m[particleIndex][attributeIndex][iii]);
		else
			sum += particles_m[particleIndex][attributeIndex][iii];
	}
	return sum / numberOfParticlesPerType_m;
}

double Simulation::calculateStdDevOfParticleAttribute(int particleIndex, int attributeIndex)
{
	double stdDev{ 0 };
	double mean{ calculateMeanOfParticleAttribute(particleIndex, attributeIndex) };
	for (int iii = 0; iii < numberOfParticlesPerType_m; iii++)
	{
		stdDev += pow(particles_m[particleIndex][attributeIndex][iii] - mean, 2);
	}
	stdDev = sqrt(stdDev / numberOfParticlesPerType_m);
	return stdDev;
}

void Simulation::saveParticleAttributeToDisk(int particleIndex, int attributeIndex, const char* foldername, const char* name)
{
	std::string fn{ foldername };
	fn = fn + name + ".bin";
	fileIO::writeDblBin(fn.c_str(), particles_m[particleIndex][attributeIndex], numberOfParticlesPerType_m);
}

void Simulation::loadFileIntoParticleAttribute(int particleIndex, int attributeIndex, const char* foldername, const char* name)
{
	if (initialized_m)
	{
		std::cout << "Warning (loadFileToParticleAttribute): Simulation has been initialized.  Any existing data will be overwritten!!\n";
	}
	delete[] particles_m[particleIndex][attributeIndex];
	std::string fn{ foldername };
	fn = fn + name;
	particles_m[particleIndex][attributeIndex] = fileIO::readDblBin(fn.c_str(), numberOfParticlesPerType_m);
}

void Simulation::createSatellite(double altitude, bool upwardFacing, double** GPUdataPointers, std::string name)
{
	Satellite* newSat = new Satellite(altitude, upwardFacing, numberOfAttributesTracked_m, numberOfParticlesPerType_m, GPUdataPointers, name);
	satellites_m.push_back(newSat);
}

timeStruct* Simulation::createTimeStruct(std::string label)
{
	timeStruct* tS = new timeStruct;
	tS->label = label;
	tS->tp = std::chrono::steady_clock::now();
	return tS;
}

void Simulation::printTimeNowFromTimeStruct(timeStruct* tS, std::string label)
{
	std::cout << "Time measurement " << tS->label << " to " << label << ": ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tS->tp).count() << " ms\n";
}

void Simulation::printTimeDiffBtwTwoTimeStructs(timeStruct* startTS, timeStruct* endTS)
{
	std::cout << "Time measurement " << startTS->label << " to " << endTS->label << ": ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTS->tp - startTS->tp).count() << " ms\n";
}