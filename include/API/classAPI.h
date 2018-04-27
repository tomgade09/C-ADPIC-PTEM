#ifndef GPS_CLASSAPI_H
#define GPS_CLASSAPI_H

#include <memory>

#include "Simulation\Simulation.h"
#include "utils\fileIO.h"

using utils::fileIO::CSV;
using utils::fileIO::ParticleDistribution;
using utils::fileIO::DistributionFromDisk;

namespace API
{
	namespace Simulation_
	{
		DLLEXP std::unique_ptr<Simulation> create(double dt, double simMin, double simMax, std::string saveRootDir);
		DLLEXP std::unique_ptr<Simulation> load(std::string prevSimDir);
	}

	namespace Particle_
	{
		DLLEXP std::unique_ptr<Particle> create(std::string name, std::vector<std::string> attributeNames, double mass, double charge, long numParts);
	}

	namespace Satellite_
	{
		DLLEXP std::unique_ptr<Satellite> create(std::string name, std::vector<std::string> attributeNames, double altitude, bool upwardFacing, long numberOfParticles, double** partDataGPUPtr);
	}

	namespace BField_
	{
		DLLEXP std::unique_ptr<DipoleB> createDipoleB(double ILATDegrees, double errorTolerance, double ds);
		DLLEXP std::unique_ptr<DipoleB> createDipoleB(double ILATDegrees);
		DLLEXP std::unique_ptr<DipoleBLUT> createDipoleBLUT(double ILATDegrees, double simMin, double simMax, double ds_gradB, int numberOfMeasurements);
	}

	namespace EField_
	{
		DLLEXP std::unique_ptr<EField> create();
		DLLEXP std::unique_ptr<QSPS> createQSPS(std::vector<double> altMin, std::vector<double> altMax, std::vector<double> magnitude);
	}

	namespace CSV_
	{
		DLLEXP std::unique_ptr<CSV> create(std::string filename);
	}

	namespace ParticleDistribution_
	{
		DLLEXP std::unique_ptr<ParticleDistribution> create(std::string saveFolder, std::vector<std::string> attrNames, std::string particleName, double mass);
		DLLEXP std::unique_ptr<ParticleDistribution> createdefault();
	}

	namespace DistributionFromDisk_
	{
		DLLEXP std::unique_ptr<DistributionFromDisk> create(std::string name, std::string folder, std::string partName, std::vector<std::string> attrNames);
	}
}
#endif /* !GPS_CLASSAPI_H */