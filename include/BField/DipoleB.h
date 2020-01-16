#ifndef DIPOLEB_BFIELD_H
#define DIPOLEB_BFIELD_H

#include "BField/BField.h"
#include "physicalconstants.h"
#include "utils/unitsTypedefs.h"

class DipoleB : public BField
{
protected:
	//Field simulation constants
	meters  L_m{ 0.0 };
	meters  L_norm_m{ 0.0 };
	meters  s_max_m{ 0.0 };

	//specified variables
	degrees ILAT_m{ 0.0 };
	meters  ds_m{ 0.0 };
	double  lambdaErrorTolerance_m{ 0.0 };

	bool    useGPU_m{ true };

	//protected functions
	__host__            void    setupEnvironment() override;
	__host__            void    deleteEnvironment() override;
	__host__            void    deserialize(string serialFolder) override;

	__host__ __device__ meters  getSAtLambda(const degrees lambda) const;
	__host__ __device__ degrees getLambdaAtS(const meters s) const;

public:
	__host__ __device__ DipoleB(degrees ILAT, double lambdaErrorTolerance = 1e-4, meters ds = RADIUS_EARTH / 1000.0, bool useGPU = true);
	__host__            DipoleB(string serialFolder);
	__host__            ~DipoleB(); //device will have a default dtor created
	__host__ __device__ DipoleB(const DipoleB&) = delete;
	__host__ __device__ DipoleB& operator=(const DipoleB&) = delete;

	//for testing
	__host__            degrees ILAT()  const override;

	__host__ __device__ tesla   getBFieldAtS(const meters s, const seconds t) const override;
	__host__ __device__ double  getGradBAtS (const meters s, const seconds t) const override;
	__host__ __device__ meters  getSAtAlt   (const meters alt_fromRe) const override;

	__host__            double getErrTol() const;
	__host__            meters getds()     const;

	__host__            void serialize(string serialFolder) const override;
};

#endif