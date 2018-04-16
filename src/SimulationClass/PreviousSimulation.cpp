#include "SimulationClass\Simulation.h"

Simulation::Simulation(std::string prevSimDir) : saveRootDir_m{ prevSimDir }, simAttr_m{ std::make_unique<SimAttributes>(prevSimDir + '/' + "Simulation.attr", true) }
{
	//Load Simulation attributes
	{//blocked so vectors lose scope afterward
		std::vector<std::string> simLbls{ simAttr_m->simAD.dblLabels_m.at(0) };
		std::vector<double>      simDbls{ simAttr_m->simAD.dblAttrs_m.at(0) };
		dt_m = simDbls.at(utils::string::findAttrInd("dt", simLbls));
		simMin_m = simDbls.at(utils::string::findAttrInd("simMin", simLbls));
		simMax_m = simDbls.at(utils::string::findAttrInd("simMax", simLbls));
	}

	//Load BField Model
	setBFieldModel(simAttr_m->BAD.names_m.at(0), simAttr_m->BAD.dblAttrs_m.at(0), false);

	//Load EField Model
	for (int elem = 0; elem < simAttr_m->EAD.names_m.size(); elem++) //can be zero if there are no E Field elements
		addEFieldModel(simAttr_m->EAD.names_m.at(elem), simAttr_m->EAD.dblAttrs_m.at(elem));

	//Load Particles
	for (int part = 0; part < simAttr_m->partAD.names_m.size(); part++)
	{
		std::vector<std::string> attrNames;
		for (int strLbl = 0; strLbl < simAttr_m->partAD.strLabels_m.at(part).size(); strLbl++)
			if (simAttr_m->partAD.strLabels_m.at(part).at(strLbl) == "attrName") { attrNames.push_back(simAttr_m->partAD.strAttrs_m.at(part).at(strLbl)); }
		if (attrNames.size() == 0) { throw std::runtime_error("Simulation::Simulation (load simulation overload): Particle attrNames.size() is zero " + simAttr_m->partAD.names_m.at(part)); }
		
		createParticleType(
			simAttr_m->partAD.names_m.at(part), attrNames,
			simAttr_m->partAD.dblAttrs_m.at(part).at(utils::string::findAttrInd("mass", simAttr_m->partAD.dblLabels_m.at(part))),
			simAttr_m->partAD.dblAttrs_m.at(part).at(utils::string::findAttrInd("charge", simAttr_m->partAD.dblLabels_m.at(part))),
			simAttr_m->partAD.dblAttrs_m.at(part).at(utils::string::findAttrInd("numParts", simAttr_m->partAD.dblLabels_m.at(part))),
			simAttr_m->partAD.strAttrs_m.at(part).at(utils::string::findAttrInd("loadFilesDir", simAttr_m->partAD.strLabels_m.at(part))), false
		);
	}

	//Load Satellites
	for (int sat = 0; sat < simAttr_m->satAD.names_m.size(); sat++)
	{
		std::unique_ptr<TempSat> tmpsat{ std::make_unique<TempSat>(
			(int)(simAttr_m->satAD.dblAttrs_m.at(sat).at(utils::string::findAttrInd("partInd", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.dblAttrs_m.at(sat).at(utils::string::findAttrInd("altitude", simAttr_m->satAD.dblLabels_m.at(sat))),
			(bool)(simAttr_m->satAD.dblAttrs_m.at(sat).at(utils::string::findAttrInd("upwardFacing", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.names_m.at(sat)
		)};

		createSatellite(tmpsat.get(), false);
	}

	logFile_m = std::make_unique<LogFile>(saveRootDir_m + "/reloadSim.log", 20);
}