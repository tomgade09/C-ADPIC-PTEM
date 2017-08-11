#ifndef NUMERICALTOOLS_H
#define NUMERICALTOOLS_H

typedef double(*FncPnt1d_t)(double*, int);

double** normalDistribution_v_z(int numOfParticles, double vmean, double vsigma, double zmean, double zsigma);
double fourthOrderRungeKutta1D(FncPnt1d_t funcPointer, double* funcArg, int arrayLen, double h);

#endif