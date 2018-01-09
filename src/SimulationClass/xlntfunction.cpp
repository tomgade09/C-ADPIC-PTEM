#include "SimulationClass\Simulation.h"

void Simulation::writeDataToXLSX(std::string filename)
{
	//some sort of check that data is all where it should be...
	//header data - do here

	std::vector<std::vector<double>> writedata;

	//orig data
	for (int attrind = 0; attrind < numberOfAttributesTracked_m; attrind++)
	{
		std::vector<double> values;
		for (int partind = 0; partind < numberOfParticlesPerType_m; partind++)
			values.push_back(particlesorig_m[0][attrind][partind]);
		writedata.push_back(values);
	}

	fileIO::writeXLSXDblArray(filename, writedata, "EscapedElec", { 1,2 }, false);
	writedata.clear();

	//sat data - top escape
	for (int attrind = -1; attrind < numberOfAttributesTracked_m; attrind++)
	{
		std::vector<double> values;
		for (int partind = 0; partind < numberOfParticlesPerType_m; partind++)
		{
			if (attrind == -1)
				values.push_back(satelliteData_m[0][2][3][partind]);
			else
				values.push_back(satelliteData_m[0][2][attrind][partind]);
		}
		writedata.push_back(values);
	}

	fileIO::writeXLSXDblArray(filename, writedata, "EscapedElec", { 5,2 }, false);
	writedata.clear();

	//sat data - bottom escape
	for (int attrind = -1; attrind < numberOfAttributesTracked_m; attrind++)
	{
		std::vector<double> values;
		for (int partind = 0; partind < numberOfParticlesPerType_m; partind++)
		{
			if (attrind == -1)
				values.push_back(satelliteData_m[0][0][3][partind]);
			else
				values.push_back(satelliteData_m[0][0][attrind][partind]);
		}
		writedata.push_back(values);
	}

	fileIO::writeXLSXDblArray(filename, writedata, "EscapedElec", { 10,2 }, false);
	writedata.clear();

	//energy, pitches
	std::vector<double> energies;
	std::vector<double> pitches;
	for (int partind = 0; partind < numberOfParticlesPerType_m; partind++)
	{
		energies.push_back(0.5 * 9.109e-31 * (pow(particlesorig_m[0][0][partind], 2) + pow(particlesorig_m[0][1][partind], 2) / 1.60218e-19));
		pitches.push_back(atan2(-particlesorig_m[0][0][partind], abs(particlesorig_m[0][1][partind])) * 180 / PI);
	}
	writedata.push_back(energies);
	writedata.push_back(pitches);

	fileIO::writeXLSXDblArray(filename, writedata, "EscapedElec", { 15,2 }, false);
}