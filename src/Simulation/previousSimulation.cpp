#include "Simulation/Simulation.h"

#include "utils/strings.h"

using utils::strings::findAttrInd;

Simulation::Simulation(std::string prevSimDir) : saveRootDir_m{ prevSimDir }, /*simAttr_m{ std::make_unique<SimAttributes>(prevSimDir + "/Simulation.attr", true) },*/
	logFile_m{ std::make_unique<LogFile>(saveRootDir_m + "/reload.log", 20) }
{
	std::cout << "Start Sim Load\n";

	if (prevSimDir == std::string("..\\_dataout\\190202_11.37.17.620km.dtDivBy10\\\\") ||
		prevSimDir == std::string("..\\_dataout\\190325_23.03.53.620km.dtDivBy10\\\\") ||
		prevSimDir == std::string("..\\_dataout\\191103_23.03.48.QSPS.100eV.dtDivBy10\\\\"))
	{ //hacky work around of a bug in the saved sim attributes file for the above sim results (0.0001s dt)
		std::cout << "special case\n";

		dt_m = 0.0001;
		simMin_m = 628565.8510817;
		simMax_m = 19881647.2473464;

		setBFieldModel("DipoleBLUT", { 72.0, 637.12, 1000000 }, false);

		createParticleType("elec", MASS_ELECTRON, -1 * CHARGE_ELEM, 3456000, prevSimDir + "/bins/particles_init/", false);
		particles_m.at(0)->loadDataFromDisk(prevSimDir + "/bins/particles_final/", false);

		std::unique_ptr<TempSat> btmElec{ std::make_unique<TempSat>(0, simMin_m * 0.999, true, "btmElec") };
		std::unique_ptr<TempSat> topElec{ std::make_unique<TempSat>(0, simMax_m * 1.001, false, "topElec") };
		std::unique_ptr<TempSat> upg{ std::make_unique<TempSat>(0, 4071307.04106411, false, "4e6ElecUpg") };
		std::unique_ptr<TempSat> dng{ std::make_unique<TempSat>(0, 4071307.04106411, true, "4e6ElecDng") };
		createSatellite(btmElec.get(), false);
		createSatellite(topElec.get(), false);
		createSatellite(upg.get(), false);
		createSatellite(dng.get(), false);
		satellite(0)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/");
		satellite(1)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/");
		satellite(2)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/");
		satellite(3)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/");
		tempSats_m.push_back(std::move(btmElec));
		tempSats_m.push_back(std::move(topElec));
		tempSats_m.push_back(std::move(upg));
		tempSats_m.push_back(std::move(dng));

		return;
	}

	simAttr_m = std::make_unique<SimAttributes>(prevSimDir + "/Simulation.attr", true);
	
	//Load Simulation attributes
	{//blocked so vectors lose scope afterward
		std::vector<std::string> simLbls{ simAttr_m->simAD.dblLabels_m.at(0) };
		std::vector<double>      simDbls{ simAttr_m->simAD.dblAttrs_m.at(0) };
		dt_m = simDbls.at(findAttrInd("dt", simLbls));
		simMin_m = simDbls.at(findAttrInd("simMin", simLbls));
		simMax_m = simDbls.at(findAttrInd("simMax", simLbls));
	}

	//Load BField Model
	setBFieldModel(simAttr_m->BAD.names_m.at(0), simAttr_m->BAD.dblAttrs_m.at(0), false);


	//Load EField Model
	for (size_t entry = 0; entry < simAttr_m->EAD.names_m.size(); entry++)
		addEFieldModel(simAttr_m->EAD.names_m.at(entry), simAttr_m->EAD.dblAttrs_m.at(entry), false);

	//Load Particles
	for (size_t part = 0; part < simAttr_m->partAD.names_m.size(); part++)
	{
		createParticleType(
			simAttr_m->partAD.names_m.at(part),
			simAttr_m->partAD.dblAttrs_m.at(part).at(findAttrInd("mass", simAttr_m->partAD.dblLabels_m.at(part))),
			simAttr_m->partAD.dblAttrs_m.at(part).at(findAttrInd("charge", simAttr_m->partAD.dblLabels_m.at(part))),
			(long)simAttr_m->partAD.dblAttrs_m.at(part).at(findAttrInd("numParts", simAttr_m->partAD.dblLabels_m.at(part))),
			prevSimDir + "/bins/particles_init/", false
		);
		particles_m.at(getNumberOfParticleTypes() - 1)->loadDataFromDisk(prevSimDir + "/bins/particles_final/", false);
	}
	
	//Load Satellites
	for (size_t sat = 0; sat < simAttr_m->satAD.names_m.size(); sat++)
	{
		std::unique_ptr<TempSat> tmpsat{ std::make_unique<TempSat>(
			(int)(simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("partInd", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("altitude", simAttr_m->satAD.dblLabels_m.at(sat))),
			(bool)(simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("upwardFacing", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.names_m.at(sat)
		)};
		createSatellite(tmpsat.get(), false);
		satellite(getNumberOfSatellites() - 1)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/");
		tempSats_m.push_back(std::move(tmpsat));
	}

	std::cout << "End Sim Load\n";
}