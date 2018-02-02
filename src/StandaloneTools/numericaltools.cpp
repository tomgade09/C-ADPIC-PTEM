#include "StandaloneTools\numericaltools.h"

void normalizeArray(std::vector<double>& arrayToNorm, double normFactor, bool inverse)
{//this is duplicated from particle class, I know, I know
	if (normFactor == 1.0)
		return;

	for (int elems = 0; elems < arrayToNorm.size(); elems++) //normalize -> divide by normalization factor
		arrayToNorm.at(elems) /= (inverse ? (1 / normFactor) : (normFactor));
}

void generateNormallyDistributedValues(double* arrayToPopulate, int length, double mean, double sigma)
{
	std::random_device randDev;
	std::mt19937 mtgen(randDev());

	std::normal_distribution<> data_nd(mean, sigma);

	for (int iii = 0; iii < length; iii++)
		arrayToPopulate[iii] = data_nd(mtgen);
}

double calculateMeanOfParticleAttribute(double* arrayToRead, int length, bool absValue)
{
	double sum{ 0 };
	for (int iii = 0; iii < length; iii++)
	{
		if (absValue)
			sum += abs(arrayToRead[iii]);
		else
			sum += arrayToRead[iii];
	}
	return sum / length;
}

double calculateStdDevOfParticleAttribute(double* arrayToRead, int length)
{
	double stdDev{ 0 };
	double mean{ calculateMeanOfParticleAttribute(arrayToRead, length, false) };
	for (int iii = 0; iii < length; iii++)
	{
		stdDev += pow(arrayToRead[iii] - mean, 2);
	}
	stdDev = sqrt(stdDev / length);
	return stdDev;
}