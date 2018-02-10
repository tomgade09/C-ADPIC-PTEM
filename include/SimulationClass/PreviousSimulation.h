#ifndef SIMULATION_PREVIOUS_H
#define SIMULATION_PREVIOUS_H

#include "SimulationClass\Simulation.h"
#include "StandaloneTools\binaryFileTools.h"

constexpr int SIMATTRSTOREAD{ 6 }; //make sure these are up to date
constexpr int BDIPATTRTOREAD{ 1 };
//constexpr int BDIPLUTATRREAD{ 1? };
//constexpr int BRCUBATRTOREAD{ 1? };
//constexpr int BIGRFATRTOREAD{ 5? };
constexpr int EQSPSLUTATREAD{ 10 };
//constexpr int EALFLUTATRREAD{ 10 };
constexpr int PARTATTRTOREAD{ 6 };
constexpr int SATATTRSTOREAD{ 3 };

class PreviousSimulation : public Simulation
{
protected:
	int findAttrInd(std::string attr, std::vector<std::string> allAttrs);

public:
	PreviousSimulation(std::string prevDataDir, std::string bFieldModel, std::vector<std::string> eFieldElems, std::vector<std::string> particleNames, std::vector<std::string> satelliteNames) : Simulation(prevDataDir)
	{		
		/* Load Simulation Attributes from disk */
		std::vector<double> simAttrs(SIMATTRSTOREAD);
		std::string simAttrLblStr;
		fileIO::readDblBin(simAttrs, prevDataDir + "/_chars/Simulation.bin", SIMATTRSTOREAD);
		fileIO::readTxtFile(simAttrLblStr, prevDataDir + "/_chars/Simulation.txt");
		std::vector<std::string> simAttrLbls{ constCharToStrVec(simAttrLblStr.c_str()) };
		
		dt_m = simAttrs.at(findAttrInd("dt", simAttrLbls));
		simMin_m = simAttrs.at(findAttrInd("simMin", simAttrLbls));
		simMax_m = simAttrs.at(findAttrInd("simMax", simAttrLbls));
		ionT_m = simAttrs.at(findAttrInd("T_ion", simAttrLbls));
		magT_m = simAttrs.at(findAttrInd("T_mag", simAttrLbls));

		/* Load B Field Model Attributes from disk */
		if (bFieldModel == "DipoleB")
		{
			std::vector<double> dipoleBAttrs(BDIPATTRTOREAD);
			fileIO::readDblBin(dipoleBAttrs, prevDataDir + "/_chars/BField_DipoleB.bin", BDIPATTRTOREAD);
			setBFieldModel("DipoleB", dipoleBAttrs);
		}
		else if (bFieldModel == "DipoleBLUT")
			throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: invalid B Field Model specified");
		else if (bFieldModel == "InvRCubedB")
			throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: invalid B Field Model specified");
		else
			throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: invalid B Field Model specified");

		/* Load E Field Element Attributes from disk */
		for (int eelem = 0; eelem < eFieldElems.size(); eelem++)
		{
			if (eFieldElems.at(eelem) == "QSPS")
			{
				throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: QSPS not ready yet.");
				std::vector<double> qspsAttr(EQSPSLUTATREAD);
				std::string altMinMax;
				std::string magnitude;
				fileIO::readDblBin(qspsAttr, prevDataDir + "/_chars/EField_QSPS.bin", EQSPSLUTATREAD);
				fileIO::readTxtFile(altMinMax, prevDataDir + "/_chars/EField_QSPS_Alt.txt");
				fileIO::readTxtFile(magnitude, prevDataDir + "/_chars/EField_QSPS_Mag.txt");
				addEFieldModel("QSPS", qspsAttr, altMinMax, magnitude);
			}
			else if (eFieldElems.at(eelem) == "AlfvenLUT")
				throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: AlfvenLUT not implemented yet");
			else if (eFieldElems.at(eelem) == "AlfvenCompute")
				throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: AlfvenCompute not implemented yet");
			else
				throw std::invalid_argument ("PreviousSimulation::PreviousSimulation: invalid E Model specified");
		}

		for (int part = 0; part < particleNames.size(); part++)
		{
			std::vector<double> partAttr(PARTATTRTOREAD);
			std::string partAttrLblStr;

			fileIO::readDblBin(partAttr, prevDataDir + "/_chars/Particle_" + particleNames.at(part) + ".bin", PARTATTRTOREAD);
			fileIO::readTxtFile(partAttrLblStr, prevDataDir + "/_chars/Particle_" + particleNames.at(part) + ".txt");

			std::vector<std::string> partAttrLbls{ constCharToStrVec(partAttrLblStr.c_str()) };
			int attrNamesInd{ findAttrInd("attrNames", partAttrLbls) };
			int posDimsInd  { findAttrInd("posDims",   partAttrLbls) };
			int velDimsInd  { findAttrInd("velDims",   partAttrLbls) };

			std::vector<std::string> attrNames;
			for (int attr = 0; attr < partAttr.at(posDimsInd) + partAttr.at(velDimsInd); attr++)
				attrNames.push_back(partAttrLbls.at(attrNamesInd + 1 + attr));
			
			createParticleType(
				partAttrLbls.at(attrNamesInd + 1 + static_cast<int>(partAttr.at(posDimsInd)) + static_cast<int>(partAttr.at(velDimsInd))),
				attrNames,
				partAttr.at(findAttrInd("mass", partAttrLbls)),
				partAttr.at(findAttrInd("charge", partAttrLbls)),
				static_cast<long>(partAttr.at(findAttrInd("numParts", partAttrLbls))),
				static_cast<int>(partAttr.at(posDimsInd)),
				static_cast<int>(partAttr.at(velDimsInd)),
				partAttr.at(findAttrInd("normFactor", partAttrLbls)));
			particleTypes_m.at(part)->loadFilesToArray(prevDataDir + "/bins/particles_init/", true);
			particleTypes_m.at(part)->loadFilesToArray(prevDataDir + "/bins/particles_final/");
		}

		for (int sat = 0; sat < satelliteNames.size(); sat++)
		{
			std::vector<double> satAttr(SATATTRSTOREAD);
			std::string satStr;
			fileIO::readDblBin(satAttr, prevDataDir + "/_chars/Satellite_" + satelliteNames.at(sat) + ".bin", SATATTRSTOREAD);
			fileIO::readTxtFile(satStr, prevDataDir + "/_chars/Satellite_" + satelliteNames.at(sat) + ".txt");

			int partInd{ (int)(satAttr.at(findAttrInd("partInd", constCharToStrVec(satStr.c_str())))) };
			std::vector<std::string> attrNames{ particleTypes_m.at(partInd)->getAttrNames() };
			attrNames.push_back("time");

			std::vector<std::vector<double>> tmp2D;
			for (int attr = 0; attr < attrNames.size(); attr++)
			{
				int numPart{ particleTypes_m.at(partInd)->getNumberOfParticles() };
				std::vector<double> tmp1D(numPart);
				fileIO::readDblBin(tmp1D, prevDataDir + "/bins/satellites/" + satelliteNames.at(sat) + "_" + attrNames.at(attr) + ".bin", numPart);
				tmp2D.push_back(tmp1D);
			}
			satelliteData_m.push_back(tmp2D);
		}
	}
};

#endif