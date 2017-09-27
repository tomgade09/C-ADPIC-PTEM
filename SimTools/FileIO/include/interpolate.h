#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <Eigen\Dense>
#include <vector>
#include <iostream>

#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

DLLEXPORT unsigned long indexLambdaFcn(int xsteps, int ysteps, int zsteps, int xsize, int ysize, int zsize);
Eigen::Vector3d interpolateB3D(const Eigen::Vector3d& pos, const Eigen::Vector2d& xrng, const Eigen::Vector2d& yrng, const Eigen::Vector2d& zrng, const double& binSize, const double* dataArray);
DLLEXPORT double* interpolateB3DPyWrapper(double* pos, double* xrng, double* yrng, double* zrng, double binSize, double* dataArray);

#endif //end of header guard