#ifndef DIPOLEB_BFIELD_H
#define DIPOLEB_BFIELD_H

#include "BField/BField.h"
#include "physicalconstants.h"

constexpr double B0{ 3.12e-5 }; //won't change from sim to sim

class DipoleB : public BField
{
protected:
	//Field simulation constants
	double L_m{ 0.0 };
	double L_norm_m{ 0.0 };
	double s_max_m{ 0.0 };

	//specified variables
	double ILATDegrees_m{ 0.0 };
	double ds_m{ 0.0 };
	double errorTolerance_m{ 0.0 };

	//protected functions
	__host__            void   setupEnvironment() override;
	__host__            void   deleteEnvironment() override;
	__host__ __device__ double getSAtLambda(const double lambdaDegrees) const;
	__host__ __device__ double getLambdaAtS(const double s) const;

public:
	__host__ __device__ DipoleB(double ILATDegrees, double errorTolerance = 1e-4, double ds = RADIUS_EARTH / 1000.0);
	__host__ __device__ ~DipoleB();
	__host__ __device__ DipoleB(const DipoleB&) = delete;
	__host__ __device__ DipoleB& operator=(const DipoleB&) = delete;

	//for testing
	double ILAT()  const { return ILATDegrees_m; }
	double ds()    const { return ds_m; }
	double L()     const { return L_m; }
	double s_max() const { return s_max_m; }

	__host__            void setds(double ds) { ds_m = ds; }

	__host__ __device__ double getBFieldAtS(const double s, const double t) const override;
	__host__ __device__ double getGradBAtS (const double s, const double t) const override;
	//__host__ __device__ double getSAtBField(const double B, const double t) const override;

	__host__ double getErrTol() const { return errorTolerance_m; }
	__host__ double getds()     const { return ds_m; }
};

#endif