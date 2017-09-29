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

void Simulation::serializeParticleArray(bool excludeOutOfSim)
{//Array is structured as: [numberOfParticleTypes_m, numberOfAttributesTracked_m, number of particles[type=0] in sim,
 //part[type=0][attr=0][num=0] data, part[type=0][attr=0][num=1] data...part[type=0][attr=1][num=0] data, part[type=0][attr=1][num=1] data...
 //number of particles[type=1] in sim,
 //part[type=1][attr=0][num=0] data, part[type=1][attr=0][num=1] data...part[type=1][attr=1][num=0] data, part[type=1][attr=1][num=1] data...]
 //So basically, the first element is the number of particle types in the simulation, the second is the number of attributes being tracked per particle,
 //Then the data of each particle is appended starting with the total number of particles of that type in the simulation (so a script can know the total size of the array)
 //Followed by the first attribute's data (particle by particle), followed by the second attribute (of that particle)'s data (particle by particle) and so on
 //Once the all of the particle's data is appended, the next particle's data is added in the same manner
 //Note: this is just one way of doing it.  The function is virtual so it can be overwritten in derived classes if desired.
	if (particlesSerialized_m != nullptr)
	{
		std::cout << "Warning: Serialized Particle array in class in NOT nullptr.  This means something has been assigned to it previously.  ";
		std::cout << "Assigning pointer to vector<void*> otherMemoryPointers.  It's your responsibility to go delete it (or use again).\n";
		otherMemoryPointers_m.push_back(particlesSerialized_m);
		std::cout << "Vector index at: " << otherMemoryPointers_m.size() - 1;
	}

	std::vector<int> partInSimCnt;
	partInSimCnt.reserve(numberOfParticleTypes_m);
	long arraySize{ 2 + numberOfParticleTypes_m };
	
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		partInSimCnt[iii] = 0;
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
		{
			if (excludeOutOfSim)
			{
				if (particlesInSim_m[iii][jjj])
					partInSimCnt[iii]++;
			}
			else
				partInSimCnt[iii]++;
		}
		arraySize += partInSimCnt[iii] * numberOfAttributesTracked_m;
	}
	particlesSerialized_m = new double[arraySize];

	particlesSerialized_m[0] = static_cast<double>(numberOfParticleTypes_m);
	particlesSerialized_m[1] = static_cast<double>(numberOfAttributesTracked_m);

	long prevPartInSim{ 0 };
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		particlesSerialized_m[2 + iii + iii * numberOfAttributesTracked_m * prevPartInSim] = static_cast<double>(partInSimCnt[iii]);
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
		{
			long particleIndex{ 0 };
			for (int kk = 0; kk < numberOfParticlesPerType_m; kk++)
			{
				if (excludeOutOfSim)
				{
					if (particlesInSim_m[iii][kk])
					{
						long eidx{ 3 + iii + iii * numberOfAttributesTracked_m * prevPartInSim + jjj * partInSimCnt[iii] + particleIndex };
						particleIndex++;
						particlesSerialized_m[eidx] = particles_m[iii][jjj][kk];
					}
				}
				else
				{
					long eidx{ 3 + iii * numberOfAttributesTracked_m * numberOfParticlesPerType_m + jjj * numberOfParticlesPerType_m + kk };
					particlesSerialized_m[eidx] = particles_m[iii][jjj][kk];
				}
			}//end kk
		}//end jjj
		prevPartInSim = partInSimCnt[iii];
	}//end iii
}//end function