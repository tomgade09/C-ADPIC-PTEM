#include <iostream>
#include <sstream>

#include "equalityOperators.h"
#include "testHelperUtils.h"
#include "memoryManagement.cuh"
#include "Simulation/Simulation.h" //includes for all other classes in this header
#include "utils/fileIO.h"
#include "utils/string.h"

using utils::fileIO::ParticleDistribution;
using utils::string::strVecToStr;

namespace test
{
	//Structs, templates for below function
	struct devMemMsmt
	{
		std::string label;
		size_t free;
		size_t total;

		devMemMsmt(std::string lbl) : label{ lbl }, free{ 0 }, total{ 0 } {}
	};

	template <typename T, typename...Args>
	void measure(std::function<void(std::string)> memchk, std::string label, Args... args)
	{//these create then destroy unique_ptr's to an instance of T, then check the GPU memory available to make sure there are no memory leaks
		{
			std::unique_ptr<T> tmp{ std::make_unique<T>(args...) };
		}
		memchk(label);
	}

	template <typename T>
	void measure(std::function<void(std::string)> memchk, std::string label, std::unique_ptr<T> ptr)
	{
		{
			ptr = nullptr;
		}
		memchk(label);
	}
	//End

	bool memLeakGPU()
	{
		std::vector<devMemMsmt> msmts;
		auto getMem = [&](std::string label)
		{
			msmts.push_back(devMemMsmt(label));
			checkGPUMemory(msmts.back().free, msmts.back().total);
		};

		measure<fieldshell<BField>, const char*, double>(getMem, "init", "Bfs", 1.0); //take an initial measurement
		measure<fieldshell<BField>, const char*, double>(getMem, "BField", "Bfs", 1.0); //run BField again, to make sure no leaks
		measure<fieldshell<EElem>, const char*, double>(getMem, "EElem", "Efs", 3.0);
		measure<DipoleB>(getMem, "DipoleB", std::make_unique<DipoleB>(72.0));
		measure<DipoleBLUT>(getMem, "DipoleBLUT", std::make_unique<DipoleBLUT>(72.0, 1.0e6, 2.0e7, 600.0, 1000000));
		measure<EField>(getMem, "EField", std::make_unique<EField>());
		measure<QSPS>(getMem, "QSPS", std::make_unique<QSPS>(std::vector<double>{ 1.0e6 }, std::vector<double>{ 1.0e7 }, std::vector<double>{ 0.02 }));
		measure<Particle>(getMem, "Particle", std::make_unique<Particle>("elec", std::vector<std::string>{ "vpara", "vperp", "s", "t_inc", "t_esc" }, 9.109e-31, -1.6e-19, 3456000));
		measure<Satellite>(getMem, "Satellite", std::make_unique<Satellite>("4e6ElecDn", std::vector<std::string>{"vpara", "vperp", "s", "time", "index"}, 4.0e6, false, 3456000, nullptr));

		size_t free_baseline{ msmts.front().free };
		auto compareSize = [&free_baseline](const devMemMsmt& compare)
		{ 
			bool eq{ compare.free == free_baseline };
			TEST_RESULTS("GPU Memory Leak Test: " + compare.label, eq);
			free_baseline = compare.free;
			return eq;
		};

		#ifdef TESTS_VERBOSE
		std::cout << "GPU Memory: Free/Total:\n";
		for (auto& msmt : msmts)
			std::cout << msmt.free << " / " << msmt.total << ": " << msmt.label << "\n";
		#endif

		bool pass{ true };
		for (auto& msmt : msmts)
			pass = compareSize(msmt) && pass; //I want this to execute everytime (and thus cout), so I put && pass at the end

		return pass;
	}

	//template for below function
	template <typename T>
	bool checkPrintEqual(std::string name, T x, T y)
	{
		TEST_RESULTS(name, (x == y));
		if (x != y)
		{
			std::cout << name << " : " << x << ", " << y << std::endl;
			return false;
		}
		return true;
	}

	bool simAttributeSaving(int runs)
	{
		///SETUP
		//grab cout, throw away output so it doesn't print to screen
		std::streambuf* coutbak{ std::cout.rdbuf() };
		std::stringstream throwaway; //temporary stringstream to capture anything written to cout
		std::cout.rdbuf(throwaway.rdbuf());

		//write particle distributions to "./out/" directory
		std::unique_ptr<ParticleDistribution> elec{ std::make_unique<ParticleDistribution>("./out/", std::vector<std::string>{ "vpara", "vperp", "s", "t_inc", "t_esc" }, "elec") };
		elec->addEnergyRange(256, 0.5, 4.5);
		elec->addPitchRange(256, 0.0, 180.0);
		elec->generate(101322.378940846, 19881647.2473464);

		elec = std::make_unique<ParticleDistribution>("./out/", std::vector<std::string>{ "vpara", "vperp", "s", "t_inc", "t_esc" }, "bselec"); //write previous distribution and create new distribution
		elec->addEnergyRange(256, 0.5, 4.5);
		elec->addPitchRange(256, 0.0, 180.0);
		elec->generate(101322.378940846, 19881647.2473464);
		elec = nullptr; //write previous distribution

		//setup and iterate new simulation
		std::unique_ptr<Simulation> sim{ std::make_unique<Simulation>(0.001, 101322.378940846, 19881647.2473464, "./out/") };

		int numParts{ 256 * 256 };
		TEST_EXCEP_CHECK(
			double simMin{ sim->simMin() };
			double simMax{ sim->simMax() };

		sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
		sim->addEFieldModel("QSPS", { 400000, 500000, 0.75261378926, 550000, 616161, 0.8456156454 });
		sim->createParticleType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, "./out/");
		sim->createParticleType("bselec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, "./out/");
		sim->createTempSat(0, simMin * 0.999, true, "btmElec");
		sim->createTempSat(0, simMax * 1.001, false, "topElec");
		sim->createTempSat(0, 4071307.04106411, false, "elecDn4e6"); //altitude = 4000km
		sim->createTempSat(0, 4071307.04106411, true, "elecUp4e6");
		); /* TEST_EXCEP_CHECK() */

		TEST_EXCEP_CHECK(
			sim->initializeSimulation();
			sim->iterateSimulation(runs, 500);
		);

		//load simulation characteristics into new sim
		std::unique_ptr<Simulation> simld{ std::make_unique<Simulation>("./out/") };

		std::cout.rdbuf(coutbak); //set cout buffer back to original (so it prints to terminal again)
		///END SETUP

		bool pass{ true };

		if (*sim == *simld)
		{
			TEST_RESULTS("Simulation create, save, load", true);
		}
		else
		{
			TEST_RESULTS("Simulation create, save, load", false);
			pass = false;
			#ifdef TESTS_VERBOSE
			checkPrintEqual<double>("dt", sim->dt(), simld->dt());
			checkPrintEqual<double>("s_min", sim->simMin(), simld->simMin());
			checkPrintEqual<double>("s_max", sim->simMax(), simld->simMax());
			checkPrintEqual<int>("# particle types", sim->getNumberOfParticleTypes(), simld->getNumberOfParticleTypes());
			checkPrintEqual<int>("# satellites", sim->getNumberOfSatellites(), simld->getNumberOfSatellites());
			#endif
		}

		for (int part = 0; part < sim->getNumberOfParticleTypes(); part++)
		{
			if (*sim->particle(part) == *simld->particle(part))
			{
				TEST_RESULTS("Particle " + sim->particle(part)->name() + " create, save, load", true);
			}
			else
			{
				TEST_RESULTS("Particle " + sim->particle(part)->name() + " create, save, load", false);
				pass = false;
				#ifdef TESTS_VERBOSE
				Particle* og{ sim->particle(part) };
				Particle* ld{ simld->particle(part) };

				checkPrintEqual<std::string>("name", og->name(), ld->name());
				checkPrintEqual<double>("mass", og->mass(), ld->mass());
				checkPrintEqual<double>("charge", og->charge(), ld->charge());
				checkPrintEqual<long>("# particles", og->getNumberOfParticles(), ld->getNumberOfParticles());
				
				TEST_EXCEP_CHECK(
				checkPrintEqual<std::string>("attributes", strVecToStr(og->attrNames()), strVecToStr(ld->attrNames()));
				);
				
				TEST_RESULTS("original data", (sim->particle(part)->data(true) == simld->particle(part)->data(true)));
				TEST_RESULTS("current data", (sim->particle(part)->data(false) == simld->particle(part)->data(false)));
				#endif
			}
		}

		for (int sat = 0; sat < sim->getNumberOfSatellites(); sat++)
		{
			if (*sim->satellite(sat) == *simld->satellite(sat))
			{
				TEST_RESULTS("Satellite " + sim->satellite(sat)->name() + " create, save, load", true);
			}
			else
			{
				TEST_RESULTS("Satellite " + sim->satellite(sat)->name() + " create, save, load", false);
				pass = false;
				#ifdef TESTS_VERBOSE
				Satellite* og{ sim->satellite(sat) };
				Satellite* ld{ simld->satellite(sat) };

				checkPrintEqual<std::string>("name", og->name(), ld->name());
				checkPrintEqual<double>("altitude", og->altitude(), ld->altitude());
				checkPrintEqual<bool>("upward", og->upward(), ld->upward());
				TEST_RESULTS("data", (og->data() == ld->data()));
				#endif
			}
		}

		if (*sim->Bmodel() == *simld->Bmodel())
		{
			TEST_RESULTS("BField " + sim->Bmodel()->name() + " create, save, load", true);
		}
		else
		{
			TEST_RESULTS("BField " + sim->Bmodel()->name() + " create, save, load", false);
			pass = false;
			#ifdef TESTS_VERBOSE
			checkPrintEqual("name", sim->Bmodel()->name(), simld->Bmodel()->name());
			#endif
		}

		if (*sim->Emodel() == *simld->Emodel())
		{
			TEST_RESULTS("EField [" + sim->Emodel()->getEElemsStr() + "] create, save, load", true);
		}
		else
		{
			TEST_RESULTS("EField [" + sim->Emodel()->getEElemsStr() + "] create, save, load", false);
			pass = false;
			#ifdef TESTS_VERBOSE
			checkPrintEqual("names", sim->Emodel()->getEElemsStr(), simld->Emodel()->getEElemsStr());
			#endif
		}

		return pass;
	}
}