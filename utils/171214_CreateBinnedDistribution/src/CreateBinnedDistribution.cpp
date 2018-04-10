#include "FileIO\fileIO.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

constexpr double MASS_PROTON{ 1.6726219e-27 };   //kg
constexpr double MASS_ELECTRON{ 9.1093836e-31 }; //kg
constexpr int	 ENERGYBINS{ 96 };
constexpr int	 PITCHBINS{ 3600 }; //was 3600
constexpr int	 NUMPARTICLES{ ENERGYBINS * PITCHBINS };
constexpr double LOGEMAX{ 4.5 };
constexpr double LOGEMIN{ 0.5 };
constexpr double LOGEBINSIZE{ (LOGEMAX - LOGEMIN) / (ENERGYBINS - 1) };
//constexpr double PITCHBINSIZE{ 180.0 / (PITCHBINS) };
constexpr double ZTOP{ 19881647.2473464 };
constexpr double ZBOTTOM{ 101322.378940846 };
constexpr double PI_C{ 3.14159265358979323846 };

int main()
{
	if ((ENERGYBINS * PITCHBINS) != NUMPARTICLES)
	{
		std::cout << "Energy Bins * Pitch Bins != Number of Particles : data will be lost, index errors will occur.  Long story short, you don't want to go down this road.  Exiting.\n\n";
		return 0;
	}

	//std::cout << "NUMPARTICLES: " << NUMPARTICLES << "\nENERGYBINS: " << ENERGYBINS << "\nPITCHBINS: " << PITCHBINS << "\nLOGEMAX: " << LOGEMAX << "\nLOGEBINSIZE: " << LOGEBINSIZE;
	//std::cout << "\nPITCHBINSIZE: " << PITCHBINSIZE << "\nZTOP: " << ZTOP << "\nZTOPNORM: " << ZTOP / 6.3712e6 << "\nZBOTTOM: " << ZBOTTOM << "\nZBOTTOMNORM: " << ZBOTTOM / 6.3712e6 << std::endl;

	//Setup
	std::vector<double> energies; std::vector<double> pitches;
	std::vector<double> v_para; std::vector<double> v_perp; std::vector<double> s_top;

	energies.resize(ENERGYBINS); pitches.resize(PITCHBINS);
	v_para.resize(NUMPARTICLES); v_perp.resize(NUMPARTICLES); s_top.resize(NUMPARTICLES);

	std::vector<std::vector<double>> elec_vec; std::vector<std::vector<double>> ions_vec;

	for (int iii = 0; iii < 3; iii++)
	{
		elec_vec.push_back(std::vector<double>(NUMPARTICLES));
		ions_vec.push_back(std::vector<double>(NUMPARTICLES));
	}

	for (int iii = 0; iii < ENERGYBINS; iii++)
		energies.at(iii) = pow(10, iii * LOGEBINSIZE + LOGEMIN);
	for (int iii = 0; iii < PITCHBINS; iii++)
		pitches.at(iii) = (iii <  PITCHBINS / 2) ? (180.0 - (iii + 0.5) * 180.0 / (PITCHBINS)) : (16.0 - ((iii - PITCHBINS / 2) + 0.5) * 16.0 / (PITCHBINS / 2));
	
	//Populate Electron Data
	std::cout << "Populating and writing Electrons and Ions...\n";
	for (int iii = 0; iii < PITCHBINS; iii++)
	{
		for (int jjj = 0; jjj < ENERGYBINS; jjj++)
		{
			elec_vec.at(0).at(iii * ENERGYBINS + jjj) = -sqrt(2 * energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[iii] * PI_C / 180);
			elec_vec.at(1).at(iii * ENERGYBINS + jjj) =  sqrt(2 * energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[iii] * PI_C / 180);
			elec_vec.at(2).at(iii * ENERGYBINS + jjj) = (pitches[iii] <= 90) ? ZTOP : ZBOTTOM;
			ions_vec.at(0).at(iii * ENERGYBINS + jjj) = -sqrt(2 * energies[jjj] * 1.60218e-19 / MASS_PROTON) * cos(pitches[iii] * PI_C / 180);
			ions_vec.at(1).at(iii * ENERGYBINS + jjj) =  sqrt(2 * energies[jjj] * 1.60218e-19 / MASS_PROTON) * sin(pitches[iii] * PI_C / 180);
			ions_vec.at(2).at(iii * ENERGYBINS + jjj) = (pitches[iii] <= 90) ? ZTOP : ZBOTTOM;
		}
	}

	std::cout << "First two vpara/vperps: " << elec_vec.at(0).at(0) << "/" << elec_vec.at(1).at(0) << "  " << elec_vec.at(0).at(1) << "/" << elec_vec.at(1).at(1) << std::endl << std::endl;

	std::string folderout{ "./../../../out/" };
	std::vector<std::string> names{ "elec_vpara.bin", "elec_vperp.bin", "elec_s.bin" };
	for (int iii = 0; iii < 3; iii++)
	{
		std::string fullname{ folderout + names.at(iii) };
		fileIO::writeDblBin(elec_vec.at(iii), fullname.c_str(), NUMPARTICLES);//pointers[iii], NUMPARTICLES);
	}

	names = { "ions_vpara.bin", "ions_vperp.bin", "ions_s.bin" };
	for (int iii = 0; iii < 3; iii++)
	{
		std::string fullname{ folderout + names.at(iii) };
		fileIO::writeDblBin(ions_vec.at(iii), fullname.c_str(), NUMPARTICLES);//pointers[iii], NUMPARTICLES);
	}
	
	return 0;
}

/* OLD CODE */
	//double** elec = new double*[4];
	//double** ions = new double*[4];
	
	/*for (int iii = 0; iii < 4; iii++)
	{
		elec[iii] = new double[ENERGYBINS * PITCHBINS];
		ions[iii] = new double[ENERGYBINS * PITCHBINS];
	}*/

	//double* energies = new double[ENERGYBINS];
	//double* pitches = new double[PITCHBINS];

	//Setup energies and pitches
	/*for (int iii = 0; iii < ENERGYBINS; iii++)
		energies[iii] = pow(10, iii * LOGEBINSIZE);
	for (int iii = 0; iii < PITCHBINS; iii++)
		pitches[iii] = iii * 180.0 / (PITCHBINS - 1);*/ //PITCHBINS - 1: to get both 0 and 180 in
		//pitches[iii] = (iii < 256) ? (iii * 0.01 / 256) : (
			//(iii < 768) ? ((iii - 256) * 179.98 / 512 + 0.01) : ((iii - 768) * 0.01 / 256 + 179.98) );