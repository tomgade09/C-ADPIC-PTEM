#ifndef SA_NUMERICALTOOLS_H
#define SA_NUMERICALTOOLS_H

#include <random>
#include <cmath>

//Numerical tools
void generateNormallyDistributedValues(double* arrayToPopulate, int length, double mean, double sigma);
double calculateMeanOfParticleAttribute(double* arrayToRead, int length, bool absValue=false);
double calculateStdDevOfParticleAttribute(double* arrayToRead, int length);

#endif