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
	
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfNormalAttributesPerPart; jjj++)
		{
			std::random_device randDev;  //means and sigmas need to be passed in as: [particle1-attr1 mean/sigma, attr2 mean/sigma...
			std::mt19937 mtgen(randDev());											//particle2-attr1 mean/sigma, attr2 mean/sigma...]

			std::normal_distribution<> data_nd(means[iii * numberOfNormalAttributesPerPart + jjj], sigmas[iii * numberOfNormalAttributesPerPart + jjj]);

			for (int kkk = 0; kkk < numberOfParticlesPerType_m; kkk++)
			{
				particles_m[iii][jjj][kkk] = data_nd(mtgen);
			}//end for kkk
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

void Simulation::serializeParticleArray()
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
		otherMemoryPointers.push_back(particlesSerialized_m);
		std::cout << "Vector index at: " << otherMemoryPointers.size();
	}

	std::vector<int> part_in_sim;
	part_in_sim.reserve(numberOfParticleTypes_m);
	long arraySize{ 2 + numberOfParticleTypes_m };
	
	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		part_in_sim[iii] = 0;
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
		{
			if (particlesInSim_m[iii][jjj])
				part_in_sim[iii]++;
		}
		arraySize += part_in_sim[iii] * numberOfAttributesTracked_m;
	}
	particlesSerialized_m = new double[arraySize];

	particlesSerialized_m[0] = static_cast<double>(numberOfParticleTypes_m);
	particlesSerialized_m[1] = static_cast<double>(numberOfAttributesTracked_m);

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		particlesSerialized_m[2 + iii + iii * numberOfAttributesTracked_m * part_in_sim[iii]] = static_cast<double>(part_in_sim[iii]);
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
		{
			for (int kkk = 0; kkk < part_in_sim[iii]; kkk++)
			{
				if (particlesInSim_m[iii][kkk])
				{
					long eidx{ 3 + iii + iii * numberOfAttributesTracked_m * part_in_sim[iii] + jjj * part_in_sim[iii] + kkk };
					particlesSerialized_m[eidx] = particles_m[iii][jjj][kkk];
				}//end if
			}//end kkk
		}//end jjj
	}//end iii
}//end function