#ifndef NUMERICAL_H
#define NUMERICAL_H

typedef Eigen::Vector3d (* FncPnt3dVector_t)(double*, int); //probably can move out to a general include file
//typedef Eigen::Vector2d (* someName)(double*, int);

Eigen::Vector3d calcLorentz(Eigen::Vector3d p, Eigen::Vector3d v, double q);
Eigen::Vector3d calcMirror(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m, double scaledLength);
Eigen::Vector3d calcBatP(Eigen::Vector3d p);
Eigen::Vector3d calcEatP(Eigen::Vector3d p);
Eigen::Vector3d fourthOrderRungeKutta3DVector(FncPnt3dVector_t funcPointer, double* funcArg, int arrayLen, double h);
Eigen::Vector3d totalAccelCB(double* args, int len);

#endif