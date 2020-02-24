#include <fstream>
#include <filesystem>

//CUDA includes
#include "ErrorHandling/simExceptionMacros.h"

#include "Particle/Particle.h"
#include "utils/fileIO.h"
#include "utils/arrayUtilsGPU.h"
#include "utils/serializationHelpers.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::to_string;
using std::logic_error;
using std::runtime_error;
using std::invalid_argument;

using utils::fileIO::readDblBin;
using utils::fileIO::writeDblBin;
using namespace utils::fileIO::serialize;

Particle::Particle(string name, vector<string> attributeNames, double mass, double charge, long numParts) :
	name_m{ name }, attributeNames_m{ attributeNames }, mass_m{ mass }, charge_m{ charge }, numberOfParticles_m{ numParts }
{
	origData_m = vector<vector<double>>(attributeNames.size(), vector<double>(numParts));
	currData_m = vector<vector<double>>(attributeNames.size(), vector<double>(numParts));

	//multi-GPU: probably do like this: make this class aware of Environment, initialize on all GPUs where use_m = true with a loop or something
	initializeGPU();
}

Particle::Particle(string serialFolder, string name)
{ //for loading serialzed class
	deserialize(serialFolder, name);
	initializeGPU();
}

Particle::~Particle()
{
	freeGPUMemory(); //don't forget to free on multiple GPUs - probably do same as initializing
}

// ================ Particle - protected ================ //
void Particle::initializeGPU()
{
	utils::GPU::setup2DArray(&currData1D_d, &currData2D_d, attributeNames_m.size() + 1, numberOfParticles_m);
	initializedGPU_m = true;
}

void Particle::copyDataToGPU(bool origToGPU)
{
	if (!initializedGPU_m)
		throw logic_error("Particle::copyDataToGPU: GPU memory has not been initialized yet for particle " + name_m);
	if (!initDataLoaded_m)
		throw logic_error("Particle::copyDataToGPU: data not loaded from disk with Particle::loadDataFromDisk or generated with Particle::generateRandomParticles " + name_m);

	if (origToGPU) utils::GPU::copy2DArray(origData_m, &currData1D_d, true); //sending to orig(host) -> data(GPU)
	else           utils::GPU::copy2DArray(currData_m, &currData1D_d, true); //sending to curr(host) -> data(GPU)
}

void Particle::copyDataToHost()
{
	if (!initializedGPU_m)
		throw logic_error("Particle::copyDataToHost: GPU memory has not been initialized yet for particle " + name_m);

	utils::GPU::copy2DArray(currData_m, &currData1D_d, false); //coming back data(GPU) -> curr(host)
}

void Particle::freeGPUMemory()
{
	if (!initializedGPU_m) return;

	utils::GPU::free2DArray(&currData1D_d, &currData2D_d);
	
	initializedGPU_m = false;
}

//need a custom solution for this...
//file read/write exception checking (probably should mostly wrap fileIO functions)
#define FILE_RDWR_EXCEP_CHECK(x) \
	try{ x; } \
	catch(const invalid_argument& a) { cerr << __FILE__ << ":" << __LINE__ << " : " << "Invalid argument error: " << a.what() << ": continuing without loading file" << endl; cout << "FileIO exception: check log file for details" << endl; } \
	catch(...)                            { throw; }


// ================ Particle - public ================ //

vector<vector<double>>& Particle::__data(bool orig)
{ //returns a non-const version of data
	return ((orig) ? origData_m : currData_m);
}

const vector<vector<double>>& Particle::data(bool orig) const
{ //returns a const version of data
	return ((orig) ? origData_m : currData_m);
}

const vector<string>& Particle::attributeNames() const
{
	return attributeNames_m;
}

string Particle::name() const
{
	return name_m;
}

double Particle::mass() const
{
	return mass_m;
}

double Particle::charge() const
{
	return charge_m;
}

long Particle::getNumberOfParticles() const
{
	return numberOfParticles_m;
}

size_t Particle::getNumberOfAttributes() const
{
	return attributeNames_m.size();
}

bool Particle::getInitDataLoaded() const
{
	return initDataLoaded_m;
}

double** Particle::getCurrDataGPUPtr() const
{
	return currData2D_d;
}

size_t Particle::getAttrIndByName(string searchName) const
{
	for (size_t name = 0; name < attributeNames_m.size(); name++)
		if (searchName == attributeNames_m.at(name))
			return name;

	throw invalid_argument("Particle::getDimensionIndByName: specified name is not present in name array: " + searchName);
}

string Particle::getAttrNameByInd(size_t searchIndx) const
{
	if (!(searchIndx <= (attributeNames_m.size() - 1) && (searchIndx >= 0)))
		throw invalid_argument("Particle::getDimensionNameByInd: specified index is invalid: " + to_string(searchIndx));

	return attributeNames_m.at(searchIndx);
}

void Particle::setParticleSource_s(double s_ion, double s_mag)
{
	cout << "Particle::setParticleSource_s: This function is being depreciated and will be removed in a future release.\n";
	//need to create a ParticleDistribution and load into Particle::origData_m before I get rid of this
	size_t s_ind{ getAttrIndByName("s") };
	size_t v_ind{ getAttrIndByName("vpara") };

	for (size_t ind = 0; ind < origData_m.at(s_ind).size(); ind++)
	{
		if (origData_m.at(v_ind).at(ind) > 0.0)
			origData_m.at(s_ind).at(ind) = s_ion;
		else if (origData_m.at(v_ind).at(ind) < 0.0)
			origData_m.at(s_ind).at(ind) = s_mag;
		else
			throw std::logic_error("Particle::setParticleSource_s: vpara value is exactly 0.0 - load data to the origData_m array first.  Aborting.  index: " + std::to_string(ind));
	}

	copyDataToGPU();
}

void Particle::loadDataFromMem(vector<vector<double>> data, bool orig) //orig defaults to true
{
	((orig) ? origData_m = data : currData_m = data);
	numberOfParticles_m = ((orig) ? (int)origData_m.at(0).size() : (int)currData_m.at(0).size());

	if (orig) initDataLoaded_m = true; //copyDataToGPU uses this flag to ensure data is present in origData_m
	if (orig) copyDataToGPU();
}


void Particle::loadDataFromDisk(string folder, bool orig) //orig defaults to true
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(readDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));

	if (orig) initDataLoaded_m = true; //copyDataToGPU uses this flag to ensure data is present in origData_m
	if (orig) copyDataToGPU();
}

void Particle::saveDataToDisk(string folder, bool orig) const
{
	for (int attrs = 0; attrs < attributeNames_m.size(); attrs++)
		FILE_RDWR_EXCEP_CHECK(writeDblBin((orig ? origData_m.at(attrs) : currData_m.at(attrs)), folder + "/" + name_m + "_" + attributeNames_m.at(attrs) + ".bin", numberOfParticles_m));
}

void Particle::serialize(string serialFolder)
{ //saves necessary attributes about the particle to disk
	string filename{ serialFolder + string("/Particle_") + name_m + string(".ser") };

	if (std::filesystem::exists(filename))
		cerr << "Particle::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.\n";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("Particle::serialize: unable to create file: " + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	out.write(reinterpret_cast<const char*>(this), sizeof(Particle));
	writeStrBuf(serializeString(name_m));
	writeStrBuf(serializeStringVector(attributeNames_m));

	out.close();
}

void Particle::deserialize(string serialFolder, string name) //protected function
{ //recreates a saved "serialization"
	string filename{ serialFolder + string("/Particle_") + name + string(".ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("Particle::deserialize: unable to open file: " + filename);

	Particle* part{ nullptr };
	vector<char> partchar(sizeof(Particle), '\0');

	in.read(partchar.data(), sizeof(Particle));
	part = reinterpret_cast<Particle*>(partchar.data());

	name_m = deserializeString(in);
	attributeNames_m = deserializeStringVector(in);

	numberOfParticles_m = part->getNumberOfParticles();
	mass_m = part->mass();
	charge_m = part->charge();

	origData_m = vector<vector<double>>(attributeNames_m.size(), vector<double>(numberOfParticles_m));
	currData_m = vector<vector<double>>(attributeNames_m.size(), vector<double>(numberOfParticles_m));
}