#include "API\classAPI.h"

#include "ErrorHandling\simExceptionMacros.h"

namespace API
{
	namespace Simulation_
	{
		DLLEXP std::unique_ptr<Simulation> create(double dt, double simMin, double simMax, std::string saveRootDir)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<Simulation>(dt, simMin, simMax, saveRootDir));
		}

		DLLEXP std::unique_ptr<Simulation> load(std::string prevSimDir)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<Simulation>(prevSimDir));
		}
	}

	namespace Particle_
	{
		DLLEXP std::unique_ptr<Particle> create(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<Particle>(name, attributeNames, mass, charge, numParts));
		}
	}

	namespace Satellite_
	{
		DLLEXP std::unique_ptr<Satellite> create(std::string name, std::vector<std::string> attributeNames, double altitude, bool upwardFacing, long numberOfParticles, double** partDataGPUPtr)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<Satellite>(name, attributeNames, altitude, upwardFacing, numberOfParticles, partDataGPUPtr));
		}
	}

	namespace BField_
	{
		DLLEXP std::unique_ptr<DipoleB> createDipoleB(double ILATDegrees, double errorTolerance, double ds)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<DipoleB>(ILATDegrees, errorTolerance, ds));
		}

		DLLEXP std::unique_ptr<DipoleB> createDipoleB(double ILATDegrees)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<DipoleB>(ILATDegrees));
		}

		DLLEXP std::unique_ptr<DipoleBLUT> createDipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_gradB, int numberOfMeasurements)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<DipoleBLUT>(ILATDegrees, simMin, simMax, ds_gradB, numberOfMeasurements));
		}
	}

	namespace EField_
	{
		DLLEXP std::unique_ptr<EField> create()
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<EField>());
		}

		DLLEXP std::unique_ptr<QSPS> createQSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<QSPS>(altMin, altMax, magnitude));
		}
	}

	namespace CSV_
	{
		DLLEXP std::unique_ptr<CSV> create(std::string filename)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<CSV>(filename));
		}
	}

	namespace ParticleDistribution_
	{
		DLLEXP std::unique_ptr<ParticleDistribution> create(std::string saveFolder, std::vector<std::string> attrNames, std::string particleName, double mass)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<ParticleDistribution>(saveFolder, attrNames, particleName, mass));
		}

		DLLEXP std::unique_ptr<ParticleDistribution> createdefault()
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<ParticleDistribution>());
		}
	}

	namespace DistributionFromDisk_
	{
		DLLEXP std::unique_ptr<DistributionFromDisk> create(std::string name, std::string folder, std::string partName, std::vector<std::string> attrNames)
		{
			SIM_API_EXCEP_CHECK(return std::make_unique<DistributionFromDisk>(name, folder, partName, attrNames));
		}
	}
}