#include <iostream>

#ifdef _WIN32
#include "Windows.h"
#endif

#include "memoryManagement.cuh"
#include "SimulationClass\Simulation.h" //includes for all other classes in this header

//Test Macros
#define TEST_RESULTS(name, pass) { std::string test{ name }; ((pass) ? coutColor("PASS", BACKGROUND_GREEN | BACKGROUND_INTENSITY) : coutColor("FAIL", 198/*Intense Red with white text*/)); std::cout << ": " << test << "\n"; }
#define EQUAL_EXCEP_CHECK(nm,x) try{std::string name{nm}; if(!x) throw std::logic_error(name + "Not equal"); else throw std::logic_error(name + " Equal");}catch(std::exception& e){std::cout << __FILE__ << ":" << __LINE__ << " : " << e.what() << std::endl;}
#define TEST_EXCEP_CHECK(x) try{x;}catch(std::exception& e){std::cout << __FILE__ << ":" << __LINE__ << " : " << e.what() << std::endl;}

namespace test
{
	//Pretty Console Colors
	#ifdef _WIN32
	void coutColor(std::string str, WORD color) //Windows way of changing console color
	{
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get handle to console
		
		CONSOLE_SCREEN_BUFFER_INFO consoleInfo;
		GetConsoleScreenBufferInfo(hConsole, &consoleInfo); //save original attributes

		SetConsoleTextAttribute(hConsole, color); //set color
		std::cout << str; //print string

		SetConsoleTextAttribute(hConsole, consoleInfo.wAttributes); //set back to normal
	}

	void printTextOptions()
	{
		for (int iii = 0; iii < 512; iii++)
			coutColor(std::to_string(iii) + " : Test text\n", iii);
	}
	#else
	void coutColor(std::string str, int color)
	{
		std::cout << str; //for right now, color not defined on other OSes...build *nix version eventually
	}
	void printTextOptions()
	{
		coutColor("All you got for now\n", 0);
	}
	#endif

	//Forward decls, helper functions, templates for below function
	void checkGPUMemory(size_t& free, size_t& total);

	struct devMemMsmt
	{
		std::string label;
		size_t free;
		size_t total;

		devMemMsmt(std::string lbl) : label{ lbl }, free{ 0 }, total{ 0 } {}
	};

	template <typename T, typename...Args>
	void measure(std::function<void(std::string)> memchk, std::string label, Args... args)
	{
		{
			std::unique_ptr<T> tmp{ std::make_unique<T>(args...) };
		}
		memchk(label);
	}
	//End

	void memLeakGPU()
	{
		std::vector<devMemMsmt> msmts;
		auto getMem = [&](std::string label) { msmts.push_back(devMemMsmt(label)); checkGPUMemory(msmts.back().free, msmts.back().total); };
		devMemMsmt init{ "init" };
		checkGPUMemory(init.free, init.total);

		//measure<fieldshell<BField>, double>(getMem, "BField", 1.0); //have to fix these
		measure<DipoleB, double>(getMem, "DipoleB", 72.0);
		measure<DipoleBLUT, double, double, double, double, int>(getMem, "DipoleBLUT", 72.0, 1.0e6, 2.0e7, 600.0, 1000000);
		measure<EField>(getMem, "EField");
		//measure<fieldshell<EElem>, double>(getMem, "EElem", 3.0);
		measure<QSPS, std::vector<double>, std::vector<double>, std::vector<double>>(getMem, "QSPS", { 1.0e6 }, { 1.0e7 }, { 0.02 });
		measure<Particle, std::string, std::vector<std::string>, double, double, int>(getMem, "Particle", "elec", { "vpara", "vperp", "s", "t_inc", "t_esc" }, 9.109e-31, -1.6e-19, 3456000);
		measure<Satellite, std::string, std::vector<std::string>, double, bool, int, double**>(getMem, "Satellite", "4e6ElecDn", std::vector<std::string>{"vpara", "vperp", "s", "time", "index"}, 4.0e6, false, 3456000, nullptr);
		
		size_t free_baseline{ init.free };

		auto compareSize = [&](devMemMsmt compare) { bool eq{ compare.free == free_baseline };
			TEST_RESULTS("GPU Memory Leak Test: " + compare.label, eq);
			free_baseline = compare.free;
			return eq; };

		#ifdef TESTS_VERBOSE
		std::cout << "GPU Memory: Free/Total:\n";
		std::cout << init.free << " / " << init.total << ": " << init.label << "\n";
		for (auto& msmt : msmts)
			std::cout << msmt.free << " / " << msmt.total << ": " << msmt.label << "\n";
		#endif

		bool pass{ true };
		for (auto& msmt : msmts)
			pass = compareSize(msmt) && pass;
			
		TEST_RESULTS("All GPU Memory Leak Tests", pass);
		std::cout << std::endl;
	}

	void simAttributeSaving(int runs)
	{
		std::unique_ptr<Simulation> sim{ std::make_unique<Simulation>(0.001, 101322.378940846, 19881647.2473464, "./out/") };

		int numParts{ 345600 };
		TEST_EXCEP_CHECK(
			double simMin{ sim->simMin() };
			double simMax{ sim->simMax() };

		sim->setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 });
		sim->addEFieldModel("QSPS", { 400000, 500000, 0.75261378926, 550000, 616161, 0.8456156454 });
		sim->createParticleType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, "./../../../../../_in/data/");
		sim->createParticleType("bselec", MASS_ELECTRON, -1 * CHARGE_ELEM, numParts, "./../../../../../_in/data/");
		sim->createTempSat(0, simMin * 0.999, true, "btmElec");
		sim->createTempSat(0, simMax * 1.001, false, "topElec");
		sim->createTempSat(0, 4071307.04106411, false, "elecDn4e6"); //altitude = 4000km
		sim->createTempSat(0, 4071307.04106411, true, "elecUp4e6");
		); /* TEST_EXCEP_CHECK() */

		TEST_EXCEP_CHECK(
			sim->initializeSimulation();
			sim->iterateSimulation(runs, 500);
		);

		std::unique_ptr<Simulation> simld{ std::make_unique<Simulation>("./out/") };

		auto compareStrs = [](std::string name, std::string x, std::string y) { if (x != y) { std::cout << name << " is not equal: " << x << ", " << y << std::endl; return false; } return true; };

		auto compareDbls = [](std::string name, double x, double y) { if (x != y) { std::cout << name << " is not equal: " << x << ", " << y << std::endl; return false; } return true; }; //T/F refers to equal?

		auto compareVecs = [&](std::string name, std::vector<double> x, std::vector<double> y) { if (x.size() != y.size()) std::cout << name << " sizes are not equal: " << x.size() << ", " << y.size() << std::endl;
		for (int elem = 0; elem < x.size(); elem++)	if (!compareDbls(name + "[" + std::to_string(elem) + "]", x.at(elem), y.at(elem))) /* if not equal: */return false;
		/* if equal: */return true; };

		auto compareVc2D = [&](std::string name, std::vector<std::vector<double>> x, std::vector<std::vector<double>> y) { if (x.size() != y.size()) std::cout << name << " sizes are not equal: " << x.size() << ", " << y.size() << std::endl;
		for (int vec = 0; vec < x.size(); vec++) if (!compareVecs(name + "[" + std::to_string(vec) + "]", x.at(vec), y.at(vec))) /* if not equal: */return false;
		/* if equal: */return true; };

		//Simulation Compares
		std::cout << "Simulation:\n";
		EQUAL_EXCEP_CHECK("dt", compareDbls(name, sim->dt(), simld->dt()));
		EQUAL_EXCEP_CHECK("simMin", compareDbls(name, sim->simMin(), simld->simMin()));
		EQUAL_EXCEP_CHECK("simMax", compareDbls(name, sim->simMax(), simld->simMax()));
		EQUAL_EXCEP_CHECK("Part Types", compareDbls(name, (double)sim->getNumberOfParticleTypes(), (double)simld->getNumberOfParticleTypes()));
		EQUAL_EXCEP_CHECK("Attributes0", compareDbls(name, (double)sim->getNumberOfAttributes(0), (double)simld->getNumberOfAttributes(0)));
		EQUAL_EXCEP_CHECK("Attributes1", compareDbls(name, (double)sim->getNumberOfAttributes(1), (double)simld->getNumberOfAttributes(1)));
		EQUAL_EXCEP_CHECK("#Satellites", compareDbls(name, (double)sim->getNumberOfSatellites(), (double)simld->getNumberOfSatellites()));
		std::cout << "\n";

		//Particle Compares
		std::cout << "Particle:\n";
		for (int part = 0; part < sim->getNumberOfParticleTypes(); part++)
		{
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " name", compareStrs(name, sim->particle(part)->name(), simld->particle(part)->name()));
			TEST_EXCEP_CHECK(if (sim->particle(part)->getAttrNames().size() != simld->particle(part)->getAttrNames().size()) std::cout << "Part " << part << " attrs name array not the same size: " << sim->particle(part)->getAttrNames().size() << ", " << simld->particle(part)->getAttrNames().size() << std::endl);
			TEST_EXCEP_CHECK(
				for (int str = 0; str < sim->particle(part)->getAttrNames().size(); str++)
					EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " attr", compareStrs(name, sim->particle(part)->getAttrNames().at(str), simld->particle(part)->getAttrNames().at(str)));
			);
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " mass", compareDbls(name, sim->particle(part)->mass(), simld->particle(part)->mass()));
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " charge", compareDbls(name, sim->particle(part)->charge(), simld->particle(part)->charge()));
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " #parts", compareDbls(name, (double)sim->getNumberOfParticles(part), (double)simld->getNumberOfParticles(part)));
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " orig", compareVc2D(name, sim->particle(part)->getOrigData(), simld->particle(part)->getOrigData()));
			EQUAL_EXCEP_CHECK("Part" + std::to_string(part) + " curr", compareVc2D(name, sim->particle(part)->getCurrData(), simld->particle(part)->getCurrData()));
		}
		std::cout << "\n";

		//Satellite Compares
		std::cout << "Satellite:\n";
		for (int sat = 0; sat < sim->getNumberOfSatellites(); sat++)
		{
			EQUAL_EXCEP_CHECK("Sat" + std::to_string(sat) + " name", compareStrs(name, sim->satellite(sat)->name(), simld->satellite(sat)->name()));
			EQUAL_EXCEP_CHECK("Sat" + std::to_string(sat) + " alt", compareDbls(name, sim->satellite(sat)->altitude(), simld->satellite(sat)->altitude()));
			EQUAL_EXCEP_CHECK("Sat" + std::to_string(sat) + " upw", compareDbls(name, (double)sim->satellite(sat)->upward(), (double)simld->satellite(sat)->upward()));
			TEST_EXCEP_CHECK(if (sim->satellite(sat)->data().size() != simld->satellite(sat)->data().size()) std::cout << "Sat" << sat << " data size (msmts) not equal: " << sim->satellite(sat)->data().size() << ", " << simld->satellite(sat)->data().size() << std::endl);
			TEST_EXCEP_CHECK(
				for (int msmt = 0; msmt < sim->satellite(sat)->data().size(); msmt++)
					EQUAL_EXCEP_CHECK("Sat" + std::to_string(sat) + " msmt" + std::to_string(msmt), compareVc2D(name, sim->satellite(sat)->data().at(msmt), simld->satellite(sat)->data().at(msmt)));
			);
			std::cout << "\n";
		}
		std::cout << "\n";

		//BField Compares
		std::cout << "BField:\n";
		EQUAL_EXCEP_CHECK("DipoleBLUT name", compareStrs(name, sim->Bmodel()->name(), simld->Bmodel()->name()));
		EQUAL_EXCEP_CHECK("DipoleBLUT errtol", compareDbls(name, ((DipoleBLUT*)sim->Bmodel())->getErrTol(), ((DipoleBLUT*)simld->Bmodel())->getErrTol()));
		EQUAL_EXCEP_CHECK("DipoleBLUT ILAT", compareDbls(name, ((DipoleBLUT*)sim->Bmodel())->ILAT(), ((DipoleBLUT*)simld->Bmodel())->ILAT()));
		EQUAL_EXCEP_CHECK("DipoleBLUT ds_gradB", compareDbls(name, ((DipoleBLUT*)sim->Bmodel())->ds_gradB(), ((DipoleBLUT*)simld->Bmodel())->ds_gradB()));
		EQUAL_EXCEP_CHECK("DipoleBLUT ds_msmt", compareDbls(name, ((DipoleBLUT*)sim->Bmodel())->ds_msmt(), ((DipoleBLUT*)simld->Bmodel())->ds_msmt()));
		std::cout << "\n";

		//EField Compares
		std::cout << "EField:\n";
		EQUAL_EXCEP_CHECK("EField QSPS altMin", compareVecs(name, ((QSPS*)sim->Emodel()->element(0))->altMin(), ((QSPS*)simld->Emodel()->element(0))->altMin()));
		EQUAL_EXCEP_CHECK("EField QSPS altMax", compareVecs(name, ((QSPS*)sim->Emodel()->element(0))->altMax(), ((QSPS*)simld->Emodel()->element(0))->altMax()));
		EQUAL_EXCEP_CHECK("EField QSPS magnitude", compareVecs(name, ((QSPS*)sim->Emodel()->element(0))->magnitude(), ((QSPS*)simld->Emodel()->element(0))->magnitude()));
		EQUAL_EXCEP_CHECK("EField QSPS EField Test", compareDbls(name, sim->Emodel()->getEFieldAtS(400001, 0.0), simld->Emodel()->getEFieldAtS(400001, 0.0)));
	}
}