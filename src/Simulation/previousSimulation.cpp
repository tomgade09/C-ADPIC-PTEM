#include "Simulation/Simulation.h"

#include "utils/string.h"

using utils::string::findAttrInd;

Simulation::Simulation(std::string prevSimDir) : saveRootDir_m{ prevSimDir }, simAttr_m{ std::make_unique<SimAttributes>(prevSimDir + '/' + "Simulation.attr", true) },
	logFile_m{ std::make_unique<LogFile>(saveRootDir_m + "reload.log", 20) }
{
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
	for (int entry = 0; entry < (int)simAttr_m->EAD.names_m.size(); entry++)
		addEFieldModel(simAttr_m->EAD.names_m.at(entry), simAttr_m->EAD.dblAttrs_m.at(entry), false);

	//Load Particles
	for (int part = 0; part < simAttr_m->partAD.names_m.size(); part++)
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
	for (int sat = 0; sat < simAttr_m->satAD.names_m.size(); sat++)
	{
		std::unique_ptr<TempSat> tmpsat{ std::make_unique<TempSat>(
			(int)(simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("partInd", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("altitude", simAttr_m->satAD.dblLabels_m.at(sat))),
			(bool)(simAttr_m->satAD.dblAttrs_m.at(sat).at(findAttrInd("upwardFacing", simAttr_m->satAD.dblLabels_m.at(sat)))),
			simAttr_m->satAD.names_m.at(sat)
		)};
		createSatellite(tmpsat.get(), false);
		satellite(getNumberOfSatellites() - 1)->loadDataFromDisk(saveRootDir_m + "/bins/satellites/", true);
	}
}