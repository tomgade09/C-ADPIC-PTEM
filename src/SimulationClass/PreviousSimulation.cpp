#include "SimulationClass\Simulation.h"

Simulation::Simulation(std::string prevSimDir) : saveRootDir_m{ prevSimDir }
{
	/* Load names from disk */
	std::string prevAttrDir{ saveRootDir_m + "/_chars/" };
	std::string prevDataDir{ saveRootDir_m + "/bins/" };
	
	std::string particlesStr;
	FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(particlesStr, prevAttrDir + "Particles.txt"));
	std::vector<std::string> particleNames{ utils::string::charToStrVec(particlesStr.c_str()) };

	std::string satellitesStr;
	FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(satellitesStr, prevAttrDir + "Satellites.txt"));
	std::vector<std::string> satelliteNames{ utils::string::charToStrVec(satellitesStr.c_str()) };

	std::string              bFieldModel{ utils::string::discoverBFieldType(prevAttrDir) };
	std::vector<std::string> eFieldElems{ utils::string::discoverEFieldTypes(prevAttrDir) };


	/* Load Simulation Attributes from disk */
	std::string simLblStr;
	FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(simLblStr, prevAttrDir + "Simulation.txt"));
	std::vector<std::string> simAttrLbls{ utils::string::charToStrVec(simLblStr.c_str()) };

	int simattrsize{ /*(int)simAttrLbls.size()*/ 3 }; //changle later, back to the commented out value
	std::vector<double> simAttrs(simattrsize);
	FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin(simAttrs, prevAttrDir + "Simulation.bin", simattrsize));

	dt_m =     simAttrs.at(utils::string::findAttrInd("dt", simAttrLbls));
	simMin_m = simAttrs.at(utils::string::findAttrInd("simMin", simAttrLbls));
	simMax_m = simAttrs.at(utils::string::findAttrInd("simMax", simAttrLbls));


	/* Load B Field Model Attributes from disk */
	if (bFieldModel == "DipoleB")
	{
		int battrsize{ utils::string::sizeofStrVecFromFile(prevAttrDir + "BField_DipoleB.txt") };
		std::vector<double> dipoleBAttrs(battrsize);
		FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin(dipoleBAttrs, prevAttrDir + "BField_DipoleB.bin", battrsize));
		setBFieldModel("DipoleB", dipoleBAttrs);
	}
	else if (bFieldModel == "DipoleBLUT")
	{
		int battrsize{ utils::string::sizeofStrVecFromFile(prevAttrDir + "BField_DipoleBLUT.txt") };
		std::vector<double> dipoleBLUTAttrs(battrsize);
		FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin(dipoleBLUTAttrs, prevAttrDir + "BField_DipoleBLUT.bin", battrsize));
		setBFieldModel("DipoleBLUT", dipoleBLUTAttrs);
	}
	else if (bFieldModel == "InvRCubedB")
		throw std::invalid_argument("PreviousSimulation::PreviousSimulation: invalid B Field Model specified");
	else
		throw std::invalid_argument("PreviousSimulation::PreviousSimulation: invalid B Field Model specified");


	/* Load E Field Element Attributes from disk */
	for (int eelem = 0; eelem < eFieldElems.size(); eelem++)
	{
		if (eFieldElems.at(eelem) == "QSPS")
		{
			throw std::invalid_argument("PreviousSimulation::PreviousSimulation: QSPS not ready yet.");
			std::string qspsLblStr;
			FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(qspsLblStr, prevAttrDir + "EField_QSPS.txt"));
			std::vector<std::string> qspsAttrLbls{ utils::string::charToStrVec(qspsLblStr.c_str()) };
			int regions{ stoi(qspsAttrLbls.at(utils::string::findAttrInd("numRegions", qspsAttrLbls) + 1)) };

			std::vector<double> qspsAttr(3 * regions); //change later
			FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin(qspsAttr, prevAttrDir + "EField_QSPS.bin", 3 * regions));

			std::vector<std::vector<double>> qsps(3); // altMin, altMax, magnitude
			for (int iii = 0; iii < regions; iii++)
			{
				qsps.at(0).push_back(qspsAttr.at(iii)); //altMin
				qsps.at(1).push_back(qspsAttr.at(iii + regions)); //altMax
				qsps.at(2).push_back(qspsAttr.at(iii + 2 * regions)); //magnitude
			}

			addEFieldModel("QSPS", qsps);
		}
		else if (eFieldElems.at(eelem) == "AlfvenLUT")
			throw std::invalid_argument("PreviousSimulation::PreviousSimulation: AlfvenLUT not implemented yet");
		else if (eFieldElems.at(eelem) == "AlfvenCompute")
			throw std::invalid_argument("PreviousSimulation::PreviousSimulation: AlfvenCompute not implemented yet");
		else
			throw std::invalid_argument("PreviousSimulation::PreviousSimulation: invalid E Model specified");
	}


	/* Create particles and load data from disk */
	for (int part = 0; part < particleNames.size(); part++)
	{
		std::string partAttrLblStr;
		FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(partAttrLblStr, prevAttrDir + "Particle_" + particleNames.at(part) + ".txt"));
		std::vector<std::string> partAttrLbls{ utils::string::charToStrVec(partAttrLblStr.c_str()) };

		int partattrsize{ (int)partAttrLbls.size() - 5 };
		std::vector<double> partAttr(partattrsize); //five string variables
		fileIO::readDblBin(partAttr, prevAttrDir + "Particle_" + particleNames.at(part) + ".bin", partattrsize);

		int attrNamesInd{ utils::string::findAttrInd("attrNames", partAttrLbls) };
		int posDimsInd{ utils::string::findAttrInd("posDims",   partAttrLbls) };
		int velDimsInd{ utils::string::findAttrInd("velDims",   partAttrLbls) };

		std::vector<std::string> attrNames;
		for (int attr = 0; attr < partAttr.at(posDimsInd) + partAttr.at(velDimsInd); attr++)
			attrNames.push_back(partAttrLbls.at(attrNamesInd + 1 + attr));

		createParticleType(
			partAttrLbls.at(attrNamesInd + 1 + (int)(partAttr.at(posDimsInd)) + (int)(partAttr.at(velDimsInd))),
			attrNames,
			partAttr.at(utils::string::findAttrInd("mass", partAttrLbls)),
			partAttr.at(utils::string::findAttrInd("charge", partAttrLbls)),
			(long)(partAttr.at(utils::string::findAttrInd("numParts", partAttrLbls))),
			(int)(partAttr.at(posDimsInd)),
			(int)(partAttr.at(velDimsInd)),
			partAttr.at(utils::string::findAttrInd("normFactor", partAttrLbls)));
		particles_m.at(part)->loadDataFromDisk(prevDataDir + "particles_init/", true);
		particles_m.at(part)->loadDataFromDisk(prevDataDir + "particles_final/", false);
	}


	/* Create satellites and load data from disk */
	for (int sat = 0; sat < satelliteNames.size(); sat++)
	{
		std::string satLblStr;
		FILE_RDWR_EXCEP_CHECK(fileIO::readTxtFile(satLblStr, prevAttrDir + "Satellite_" + satelliteNames.at(sat) + ".txt"));
		std::vector<std::string> satAttrLbls{ utils::string::charToStrVec(satLblStr.c_str()) };

		int satattrsize{ utils::string::sizeofStrVecFromFile(prevAttrDir + "Satellite_" + satelliteNames.at(sat)) - 1 };
		std::vector<double> satAttrs(satattrsize);
		FILE_RDWR_EXCEP_CHECK(fileIO::readDblBin(satAttrs, prevAttrDir + "Satellite_" + satelliteNames.at(sat) + ".bin", satattrsize));

		TempSat tmpsat{ (int)satAttrs.at(utils::string::findAttrInd("partInd", satAttrLbls)),
			satAttrs.at(utils::string::findAttrInd("altitude", satAttrLbls)),
			(int)satAttrs.at(utils::string::findAttrInd("upwardFacing", satAttrLbls)) != 0,
			satAttrLbls.at(3) };

		createSatellite(&tmpsat);

		std::vector<std::string> attrNames{ particles_m.at((int)satAttrs.at(utils::string::findAttrInd("partInd", satAttrLbls)))->getAttrNames() };
		attrNames.push_back("time");
		attrNames.push_back("index");
		satellites_m.at(sat)->satellite->loadDataFromDisk(prevDataDir + "satellite/", attrNames);
	}

	logFile_m = std::make_unique<LogFile>(saveRootDir_m + "/reloadSim.log", 20);
}