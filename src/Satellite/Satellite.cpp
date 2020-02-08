//STL includes
#include <iterator> //for back_inserter
#include <algorithm>
#include <fstream>
#include <filesystem>

//Project specific includes
#include "Satellite/Satellite.h"
#include "utils/fileIO.h"
#include "utils/arrayUtilsGPU.h"
#include "utils/serializationHelpers.h"
#include "ErrorHandling/simExceptionMacros.h"

using std::cerr;
using std::ofstream;
using std::streambuf;
using std::to_string;
using std::logic_error;
using std::runtime_error;
using std::invalid_argument;

using utils::fileIO::readDblBin;
using utils::fileIO::writeDblBin;
using namespace utils::fileIO::serialize;

Satellite::Satellite(std::string name, std::vector<std::string> attributeNames, meters altitude, bool upwardFacing, long numberOfParticles, double** partDataGPUPtr) :
	name_m{ name }, attributeNames_m{ attributeNames }, altitude_m{ altitude }, upwardFacing_m{ upwardFacing }, numberOfParticles_m{ numberOfParticles }, particleData2D_d{ partDataGPUPtr }
{
	data_m = vector<vector<double>>(attributeNames_m.size(), vector<double>(numberOfParticles_m));
	initializeGPU();
}

Satellite::Satellite(string serialFolder, string name, double** particleData2D)
{
	deserialize(serialFolder, name, particleData2D);
	
	data_m = vector<vector<double>>(attributeNames_m.size(), vector<double>(numberOfParticles_m));
	
	initializeGPU();
}

Satellite::~Satellite()
{
	freeGPUMemory();
}

void Satellite::initializeGPU()
{
	utils::GPU::setup2DArray(&satCaptrData1D_d, &satCaptrData2D_d, attributeNames_m.size(), numberOfParticles_m);
	initializedGPU_m = true;
}

size_t Satellite::getAttrIndByName(string name)
{
	for (size_t attr = 0; attr < attributeNames_m.size(); attr++)
	{
		if (name == attributeNames_m.at(attr))
			return attr;
	}
	
	throw invalid_argument("Satellite::getAttrIndByName: attribute name " + name + " doesn't exist.");
}

//void Satellite::iterateDetectorCPU is in Detector.cpp
//void Satellite::iterateDetector is in Detector.cu

void Satellite::copyDataToHost()
{// data_m array: [v_para, mu, s, time, partindex][particle number]

	utils::GPU::copy2DArray(data_m, &satCaptrData1D_d, false);
}

void Satellite::freeGPUMemory()
{
	if (!initializedGPU_m) { return; }

	utils::GPU::free2DArray(&satCaptrData1D_d, &satCaptrData2D_d);

	particleData2D_d = nullptr;
	initializedGPU_m = false;
}

vector<vector<double>> Satellite::removeZerosData()
{//GOING TO HAVE TO REMOVE TO SOMEWHERE - IMPLEMENTATION DEFINED, NOT GENERIC
	copyDataToHost();

	vector<vector<double>> dataCopy{ data_m }; //don't want to do this to the live data so create a copy
	vector<double> timeCopy{ dataCopy.at(getAttrIndByName("time")) }; //make a copy, because t_esc's zeroes are removed as well

	for (auto& attr : dataCopy)
	{//below searches time vector copy for -1.0 and removes the element if so (no negatives should exist except -1)
		auto checkIfNegOne = [&](double& x)
		{
			return (timeCopy.at(&x - &(*attr.begin())) < 0.0);
		};

		attr.erase(remove_if(attr.begin(), attr.end(), checkIfNegOne), attr.end());
	}

	return dataCopy;
}

void Satellite::saveDataToDisk(string folder) //move B and mass to getConsolidatedData and have it convert back (or in gpu?)
{
	vector<vector<double>> results{ removeZerosData() };

	for (int attr = 0; attr < results.size(); attr++)
		writeDblBin(results.at(attr), folder + name_m + "_" + attributeNames_m.at(attr) + ".bin", (int)results.at(attr).size());
}

void Satellite::loadDataFromDisk(string folder)
{
	data_m = vector<vector<double>>(attributeNames_m.size()); //this is done so readDblBin doesn't assume how many particles it's reading

	for (size_t attr = 0; attr < attributeNames_m.size(); attr++)
		readDblBin(data_m.at(attr), folder + name_m + "_" + attributeNames_m.at(attr) + ".bin");

	bool expand{ false };
	if (data_m.at(0).size() == numberOfParticles_m)
		expand = false;
	else if (data_m.at(0).size() < numberOfParticles_m)
		expand = true;
	else
		throw logic_error("Satellite::loadDataFromDisk: number of particles loaded from disk is greater than specified numberOfParticles_m.  "
			+ string("That means that the wrong data was loaded or wrong number of particles was specified.  Not loading data."));

	if (!expand)
	{
		return;
	}
	else //expand into a sparse array
	{
		size_t index{ getAttrIndByName("index") };
		size_t t_esc{ getAttrIndByName("time") };

		for (auto& attr : data_m) //add zeroes to the array where they are missing
		{
			attr.resize(numberOfParticles_m);

			for (size_t part = numberOfParticles_m - 1; part >= 0; part--) //particles, iterating backwards
			{
				int originalParticleIndex{ (int)data_m.at(index).at(part) }; //original index of the particle

				if (originalParticleIndex == 0 && (int)data_m.at(0).at(part) == 0 && (int)data_m.at(1).at(part) == 0)
				{
					data_m.at(t_esc).at(part) = -1;
					data_m.at(index).at(part) = -1;
				}
				else if ((originalParticleIndex != part) && (originalParticleIndex != -1))
				{
					if ((int)data_m.at(0).at(originalParticleIndex) != 0.0)
						throw runtime_error("Satellite::loadDataFromDisk: data is being overwritten in reconstructed array - something is wrong " + to_string(originalParticleIndex));
					
					for (size_t attr = 0; attr < data_m.size(); attr++)
					{
						data_m.at(attr).at(originalParticleIndex) = data_m.at(attr).at(part); //move attr at the current location in iteration - part - to the index where it should be - ind
						
						if (attr == index || attr == t_esc) data_m.at(attr).at(part) = -1.0;
						else data_m.at(attr).at(part) = 0.0; //overwrite the old data with 0s and -1s
					}
				}
			}
		}
	}
}

string Satellite::name() const
{
	return name_m;
}

double Satellite::altitude() const
{
	return altitude_m;
}

bool Satellite::upward() const
{
	return upwardFacing_m;
}

vector<vector<double>>& Satellite::__data()
{
	return data_m;
}

const vector<vector<double>>& Satellite::data() const
{
	return data_m;
}

double** Satellite::get2DDataGPUPtr() const
{
	return satCaptrData2D_d;
}

double* Satellite::get1DDataGPUPtr() const
{
	return satCaptrData1D_d;
}

size_t Satellite::getNumberOfAttributes() const
{
	return attributeNames_m.size();
}

long Satellite::getNumberOfParticles() const
{
	return numberOfParticles_m;
}

void Satellite::serialize(string serialFolder)
{
	string filename{ serialFolder + string("/Satellite_") + name_m + string(".ser") };

	if (std::filesystem::exists(filename))
		cerr << "Satellite::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("Satellite::serialize: unable to create file: " + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	out.write(reinterpret_cast<const char*>(this), sizeof(Satellite));
	writeStrBuf(serializeString(name_m));
	writeStrBuf(serializeStringVector(attributeNames_m));

	out.close();
}

void Satellite::deserialize(string serialFolder, string name, double** particleData2D)
{
	string filename{ serialFolder + string("/Satellite_") + name + string(".ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("Particle::deserialize: unable to open file: " + filename);

	Satellite* sat;
	vector<char> partchar(sizeof(Satellite));

	in.read(partchar.data(), sizeof(Satellite));
	sat = reinterpret_cast<Satellite*>(partchar.data());

	name_m = deserializeString(in);
	attributeNames_m = deserializeStringVector(in);

	altitude_m = sat->altitude();
	upwardFacing_m = sat->upward();

	numberOfParticles_m = sat->getNumberOfParticles();

	particleData2D_d = particleData2D;
}
