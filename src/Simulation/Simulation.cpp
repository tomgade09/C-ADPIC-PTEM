#include "Simulation/Simulation.h"

#include <filesystem>

#include "utils/loopmacros.h"
#include "utils/serializationHelpers.h"
#include "ErrorHandling/simExceptionMacros.h"

using std::cout;
using std::cerr;
using std::endl;
using std::move;
using std::to_string;
using std::make_unique;
using std::make_shared;
using std::logic_error;
using std::out_of_range;
using std::runtime_error;
using std::invalid_argument;

using namespace utils::fileIO::serialize;

Simulation::Simulation(double dt, double simMin, double simMax, string saveRootDir) :
	dt_m{ dt }, simMin_m{ simMin }, simMax_m{ simMax }, saveRootDir_m{ saveRootDir + "/" },
	simAttr_m{ make_unique<SimAttributes>(saveRootDir + "/Simulation.attr") },
	Log_m{ make_unique<Log>(saveRootDir_m + "/simulation.log") }
{
	simAttr_m->addData("Simulation", "", {}, {}, { "dt", "simMin", "simMax" }, { dt_m, simMin_m, simMax_m });
}

Simulation::~Simulation()
{
	if (saveReady_m) saveDataToDisk(); //save data if it hasn't been done
	Log_m->createEntry("End simulation");
}


void Simulation::printSimAttributes(int numberOfIterations, int itersBtwCouts, string GPUName) //protected
{
	//Sim Header (folder) printed from Python - move here eventually
	cout << "GPU Name:       " << GPUName << "\n";
	cout << "Sim between:    " << simMin_m << "m - " << simMax_m << "m" << endl;
	cout << "dt:             " << dt_m << "s" << endl;
	cout << "BModel Model:   " << BFieldModel_m->name() << endl;
	cout << "EField Elems:   " << EFieldModel_m->getElementNames() << endl;
	cout << "Particles:      ";
	for (size_t iii = 0; iii < particles_m.size(); iii++) {
		cout << ((iii != 0) ? "                " : "") << particles_m.at(iii)->name() << ": #: " << particles_m.at(iii)->getNumberOfParticles() << ", loaded files?: " << (particles_m.at(iii)->getInitDataLoaded() ? "true" : "false") << std::endl;
	}
	cout << "Satellites:     ";
	for (int iii = 0; iii < getNumberOfSatellites(); iii++) {
		cout << ((iii != 0) ? "                " : "") << satellite(iii)->name() << ": alt: " << satellite(iii)->altitude() << " m, upward?: " << (satellite(iii)->upward() ? "true" : "false") << std::endl;
	}
	cout << "Iterations:     " << numberOfIterations << endl;
	cout << "Iters Btw Cout: " << itersBtwCouts << endl;
	cout << "Time to setup:  " << Log_m->timeElapsedTotal_s() << " s" << endl;
	cout << "===============================================================" << endl;
}

void Simulation::incTime()
{
	simTime_m += dt_m;
}

//Access functions
double Simulation::simtime() const
{
	return simTime_m;
}

double Simulation::dt() const
{
	return dt_m;
}

double Simulation::simMin() const
{
	return simMin_m;
}

double Simulation::simMax() const
{
	return simMax_m;
}

int Simulation::getNumberOfParticleTypes() const
{
	return (int)particles_m.size();
}

int Simulation::getNumberOfSatellites() const
{
	return (int)satPartPairs_m.size();
}

int Simulation::getNumberOfParticles(int partInd) const
{
	return (int)particles_m.at(partInd)->getNumberOfParticles();
}

int Simulation::getNumberOfAttributes(int partInd) const
{
	return (int)particles_m.at(partInd)->getNumberOfAttributes();
}

string Simulation::getParticlesName(int partInd) const
{
	return particles_m.at(partInd)->name();
}

string Simulation::getSatelliteName(int satInd) const
{
	return satPartPairs_m.at(satInd)->satellite->name();
}

int Simulation::getPartIndOfSat(int satInd) const
{
	return tempSats_m.at(satInd)->particleInd;
}

Particles* Simulation::particle(int partInd) const
{
	return particles_m.at(partInd).get();
}

const vector<vector<double>>& Simulation::getParticleData(int partInd, bool originalData)
{
	if (partInd > (particles_m.size() - 1))
		throw out_of_range("Simulation::getParticleData: no particle at the specifed index " + to_string(partInd));
	return ((originalData) ? (particles_m.at(partInd)->data(true)) : (particles_m.at(partInd)->data(false)));
}

const vector<vector<double>>& Simulation::getSatelliteData(int satInd)
{
	if (satInd > (satPartPairs_m.size() - 1))
		throw out_of_range("Simulation::getSatelliteData: no satellite at the specifed index " + to_string(satInd));
	return satellite(satInd)->data();
}

Particles* Simulation::particle(string name) const
{
	for (auto& part : particles_m)
		if (part->name() == name)
			return part.get();
	throw invalid_argument("Simulation::particle: no particle of name " + name);
}

Satellite* Simulation::satellite(string name) const
{
	for (auto& satPart : satPartPairs_m)
		if (satPart->satellite->name() == name)
			return satPart->satellite.get();
	throw invalid_argument("Simulation::satellite: no satellite of name " + name);
}

double Simulation::getBFieldAtS(double s, double time) const
{
	return BFieldModel_m->getBFieldAtS(s, time);
}

double Simulation::getEFieldAtS(double s, double time) const
{
	return EFieldModel_m->getEFieldAtS(s, time);
}

Satellite* Simulation::satellite(int satInd) const
{
	return satPartPairs_m.at(satInd)->satellite.get();
}

BModel* Simulation::Bmodel() const
{
	return BFieldModel_m.get();
}

EField* Simulation::Emodel() const
{
	return EFieldModel_m.get();
}

//Class creation functions
void Simulation::createParticlesType(string name, double mass, double charge, long numParts, string loadFilesDir, bool save)
{
	if (simAttr_m == nullptr)
		save = false;

	for (size_t part = 0; part < particles_m.size(); part++)
		if (name == particles_m.at(part).get()->name())
			throw invalid_argument("Simulation::createParticlesType: particle already exists with the name: " + name);

	Log_m->createEntry("Particles Type Created: " + name + ": Mass: " + to_string(mass) + ", Charge: " + to_string(charge) 
		+ ", Number of Parts: " + to_string(numParts) + ", Files Loaded?: " + ((loadFilesDir != "") ? "True" : "False"));
	
	vector<string> attrNames{ "vpara", "vperp", "s", "t_inc", "t_esc" };

	if (save)
	{
		vector<string> attrLabels;
		for (size_t atr = 0; atr < attrNames.size(); atr++)
			attrLabels.push_back("attrName");
		attrLabels.push_back("loadFilesDir");
		
		vector<string> namesCopy{ attrNames };
		namesCopy.push_back(loadFilesDir);
		simAttr_m->addData("Particles", name, attrLabels, namesCopy, { "mass", "charge", "numParts" }, { mass, charge, (double)numParts });
	}

	shared_ptr<Particles> newPart{ make_unique<Particles>(name, attrNames, mass, charge, numParts) };

	if (loadFilesDir != "")
		newPart->loadDataFromDisk(loadFilesDir);

	newPart->__data(true).at(4) = vector<double>(newPart->getNumberOfParticles(), -1.0); //sets t_esc to -1.0 - i.e. particles haven't escaped yet
	particles_m.push_back(move(newPart));
}

void Simulation::createTempSat(string partName, double altitude, bool upwardFacing, string name)
{
	for (size_t partInd = 0; partInd < particles_m.size(); partInd++)
	{
		if (particle((int)partInd)->name() == partName)
		{
			createTempSat(partInd, altitude, upwardFacing, name);
			return;
		}
	}
	throw invalid_argument("Simulation::createTempSat: no particle of name " + name);
}

void Simulation::createTempSat(int partInd, double altitude, bool upwardFacing, string name)
{//"Temp Sats" are necessary to ensure particles are created before their accompanying satellites
	if (initialized_m)
		throw runtime_error("Simulation::createTempSat: initializeSimulation has already been called, no satellite will be created of name " + name);
	if (partInd >= particles_m.size())
		throw out_of_range("Simulation::createTempSat: no particle at the specifed index " + to_string(partInd));

	tempSats_m.push_back(make_unique<TempSat>(partInd, altitude, upwardFacing, name));
}

void Simulation::createSatellite(TempSat* tmpsat, bool save) //protected
{
	int partInd{ tmpsat->particleInd };
	double altitude{ tmpsat->altitude };
	bool upwardFacing{ tmpsat->upwardFacing };
	string name{ tmpsat->name };

	if (particles_m.size() <= partInd)
		throw out_of_range("createSatellite: no particle at the specifed index " + to_string(partInd));
	if (particles_m.at(partInd)->getCurrDataGPUPtr() == nullptr)
		throw runtime_error("createSatellite: pointer to GPU data is a nullptr of particle " + particles_m.at(partInd)->name() + " - that's just asking for trouble");
	if (simAttr_m == nullptr)
		save = false;

	if (save) { simAttr_m->addData("Satellite", name, {}, {}, { "partInd", "altitude", "upwardFacing" }, { (double)partInd, altitude, (double)upwardFacing }); }

	Log_m->createEntry("Created Satellite: " + name + ", Particles tracked: " + particles_m.at(partInd)->name()
		+ ", Altitude: " + to_string(altitude) + ", " + ((upwardFacing) ? "Upward" : "Downward") + " Facing Detector");

	vector<string> attrNames{ "vpara", "vperp", "s", "time", "index" };
	shared_ptr<Particles>  part{ particles_m.at(partInd) };
	unique_ptr<Satellite> sat{ make_unique<Satellite>(name, attrNames, altitude, upwardFacing, part->getNumberOfParticles(), part->getCurrDataGPUPtr()) };
	satPartPairs_m.push_back(make_unique<SatandPart>(move(sat), move(part)));
}

void Simulation::setBFieldModel(string name, vector<double> args, bool save)
{//add log file messages
	if (BFieldModel_m)
		throw invalid_argument("Simulation::setBFieldModel: trying to assign B Field Model when one is already assigned - existing: " + BFieldModel_m->name() + ", attempted: " + name);
	if (args.empty())
		throw invalid_argument("Simulation::setBFieldModel: no arguments passed in");
	if (simAttr_m == nullptr)
		save = false;
	
	vector<string> attrNames;

	if (name == "DipoleB")
	{
		if (args.size() == 1)
		{ //for defaults in constructor of DipoleB
			BFieldModel_m = make_unique<DipoleB>(args.at(0));
			args.push_back(((DipoleB*)BFieldModel_m.get())->getErrTol());
			args.push_back(((DipoleB*)BFieldModel_m.get())->getds());
		}
		else if (args.size() == 3)
			BFieldModel_m = make_unique<DipoleB>(args.at(0), args.at(1), args.at(2));
		else
			throw invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleB: " + to_string(args.size()));

		attrNames = { "ILAT", "ds", "errTol" };
	}
	else if (name == "DipoleBLUT")
	{
		if (args.size() == 3)
			BFieldModel_m = make_unique<DipoleBLUT>(args.at(0), simMin_m, simMax_m, args.at(1), (int)args.at(2));
		else
			throw invalid_argument("setBFieldModel: wrong number of arguments specified for DipoleBLUT: " + to_string(args.size()));

		attrNames = { "ILAT", "ds", "numMsmts" };
	}
	else
	{
		cout << "Not sure what model is being referenced.  Using DipoleB instead of " << name << endl;
		BFieldModel_m = make_unique<DipoleB>(args.at(0));
		args.resize(3);
		args.at(1) = ((DipoleB*)BFieldModel_m.get())->getErrTol();
		args.at(1) = ((DipoleB*)BFieldModel_m.get())->getds();
		attrNames = { "ILAT", "ds", "errTol" };
	}

	if (save) { simAttr_m->addData("BField", name, {}, {}, attrNames, args); }
}

void Simulation::setBFieldModel(unique_ptr<BModel> BModelptr)
{
	BFieldModel_m = std::move(BModelptr);
}

void Simulation::addEFieldModel(string name, vector<double> args, bool save)
{
	if (EFieldModel_m == nullptr)
		EFieldModel_m = make_unique<EField>();
	if (simAttr_m == nullptr)
		save = false;
	
	vector<string> attrNames;

	if (name == "QSPS") //need to check to make sure args is formed properly, as well as save to disk
	{
		if (args.size() != 3) throw invalid_argument("Simulation::addEFieldModel: QSPS: Argument vector is improperly formed.  Proper format is: { altMin, altMax, mag }");
		
		vector<meters> altMin;
		vector<meters> altMax;
		vector<Vperm>  mag;

		for (size_t entry = 0; entry < args.size() / 3; entry++)
		{
			altMin.push_back(args.at(3 * entry));
			altMax.push_back(args.at(3 * entry + 1));
			mag.push_back(args.at(3 * entry + 2));
			attrNames.push_back("altMin");
			attrNames.push_back("altMax");
			attrNames.push_back("magnitude");
		}
		cout << "QSPS temporary fix - instantiates with [vector].at(0)\n";
		EFieldModel_m->add(make_unique<QSPS>(altMin.at(0), altMax.at(0), mag.at(0)));
	}
	else if (name == "AlfvenLUT")
	{
		cout << "AlfvenLUT not implemented quite yet.  Returning." << endl;
		return;
	}
	/*else if (name == "AlfvenCompute")
	{
		cout << "AlfvenCompute not implemented quite yet.  Returning." << endl;
		return;
	}*/

	if (save) { simAttr_m->addData("EField", name, {}, {}, attrNames, args); }
}


//Other utilities
void Simulation::saveDataToDisk()
{
	if (!initialized_m)
		throw logic_error("Simulation::saveDataToDisk: simulation not initialized with initializeSimulation()");
	if (!saveReady_m)
		throw logic_error("Simulation::saveDataToDisk: simulation not iterated and/or copied to host with iterateSmiulation()");

	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->saveDataToDisk(saveRootDir_m + "/bins/particles_init/", true));
	LOOP_OVER_1D_ARRAY(getNumberOfParticleTypes(), particles_m.at(iii)->saveDataToDisk(saveRootDir_m + "/bins/particles_final/", false));
	LOOP_OVER_1D_ARRAY(getNumberOfSatellites(), satellite(iii)->saveDataToDisk(saveRootDir_m + "/bins/satellites/"));

	simAttr_m = nullptr;
	saveReady_m = false;
}

void Simulation::resetSimulation(bool fields)
{
	while (satPartPairs_m.size() != 0)
		satPartPairs_m.pop_back();
	while (particles_m.size() != 0)
		particles_m.pop_back();

	if (fields)
	{
		BFieldModel_m.reset();
		EFieldModel_m.reset();
	}
}

void Simulation::save() const
{
	string filename{ saveRootDir_m + string("/simulation.ser") };

	if (std::filesystem::exists(filename))
		cerr << __func__ << ": Warning: filename exists: " << filename << " You are overwriting an existing file.";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument(__func__ + string(": unable to create file: ") + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	out.write(reinterpret_cast<const char*>(&dt_m), sizeof(double));
	out.write(reinterpret_cast<const char*>(&simMin_m), sizeof(double));
	out.write(reinterpret_cast<const char*>(&simMax_m), sizeof(double));
	writeStrBuf(serializeString(saveRootDir_m));
	out.write(reinterpret_cast<const char*>(&simTime_m), sizeof(double));
	
	//Write BField
	Component comp{ Component::BField };
	out.write(reinterpret_cast<const char*>(&comp), sizeof(Component));
	BModel::Type type{ BFieldModel_m->type() };
	out.write(reinterpret_cast<const char*>(&type), sizeof(BModel::Type));
	BFieldModel_m->serialize(out);

	//Write EField
	comp = Component::EField;
	out.write(reinterpret_cast<const char*>(&comp), sizeof(Component));
	EFieldModel_m->serialize(out);

	//Write Log
	//comp = Component::Log;
	//out.write(reinterpret_cast<const char*>(&comp), sizeof(Component));
	//Log_m->serialize(out); //Log::serialize is to be written

	//Write Particles
	comp = Component::Particles;
	for (const auto& part : particles_m)
	{
		out.write(reinterpret_cast<const char*>(&comp), sizeof(Component));
		part->serialize(out);
	}

	//Write Satellite
	comp = Component::Satellite;
	for (const auto& sat : satPartPairs_m)
	{
		out.write(reinterpret_cast<const char*>(&comp), sizeof(Component));
		for (size_t part = 0; part < particles_m.size(); part++) //write particle index
			if (sat->particle.get() == particles_m.at(part).get())
				out.write(reinterpret_cast<const char*>(&part), sizeof(size_t));
		sat->satellite->serialize(out);
	}
}

void Simulation::load(string saveRootDir)
{
	string filename{ saveRootDir + string("/simulation.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument(__func__ + string(": unable to open file: ") + filename);

	in.read(reinterpret_cast<char*>(&dt_m), sizeof(double));
	in.read(reinterpret_cast<char*>(&simMin_m), sizeof(double));
	in.read(reinterpret_cast<char*>(&simMax_m), sizeof(double));
	saveRootDir_m = deserializeString(in);
	in.read(reinterpret_cast<char*>(&simTime_m), sizeof(double));

	while (true)
	{
		Component comp{ 6 };
		in.read(reinterpret_cast<char*>(&comp), sizeof(Component));

		if (in.eof()) break;

		if (comp == Component::BField)
		{
			BModel::Type type{ BFieldModel_m->type() };
			in.read(reinterpret_cast<char*>(&type), sizeof(BModel::Type));

			if (type == BModel::Type::DipoleB)
				BFieldModel_m = make_unique<DipoleB>(in);
			else if (type == BModel::Type::DipoleBLUT)
				BFieldModel_m = make_unique<DipoleBLUT>(in);
			else throw std::runtime_error("Simulation::load: BModel type not recognized");
		}
		else if (comp == Component::EField)
		{
			EFieldModel_m = make_unique<EField>(in);
		}
		//else if (comp == Component::Log)
		//{
			//Log::deserialize to be written
			//Log_m = make_unique<Log>(in);
		//}
		else if (comp == Component::Particles)
		{
			particles_m.push_back(make_unique<Particles>(in));
		}
		else if (comp == Component::Satellite)
		{
			size_t part{ readSizetLength(in) };
			shared_ptr<Particles> particles{ particles_m.at(part) };

			satPartPairs_m.push_back(
				make_unique<SatandPart>(
					make_unique<Satellite>(in, particles_m.at(part)->getCurrDataGPUPtr()),
					move(particles)));
		}
		else throw std::runtime_error("Simulation::load: Simulation Component not recognized");
	}
}