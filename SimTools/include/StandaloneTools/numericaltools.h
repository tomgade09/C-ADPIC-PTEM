#ifndef SA_NUMERICALTOOLS_H
#define SA_NUMERICALTOOLS_H

#include <random>
#include <cmath>
#include <iostream>

//Numerical tools
void generateNormallyDistributedValues(double* arrayToPopulate, int length, double mean, double sigma);
double calculateMeanOfParticleAttribute(double* arrayToRead, int length, bool absValue=false);
double calculateStdDevOfParticleAttribute(double* arrayToRead, int length);
void normalizeArray(std::vector<double>& arrayToNorm, double normFactor, bool inverse=false);

#endif