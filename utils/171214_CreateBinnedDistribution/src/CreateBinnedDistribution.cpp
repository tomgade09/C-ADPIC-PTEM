#include "FileIO\fileIO.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

constexpr double MASS_PROTON{ 1.6726219e-27 }; //kg
constexpr double MASS_ELECTRON{ 9.1093836e-31 }; //kg
constexpr int	 NUMPARTICLES{ 1048576 };
constexpr int	 ENERGYBINS{ 1024 };
constexpr int	 PITCHBINS{ 1024 };
constexpr double LOGEMAX{ 5 };
constexpr double LOGEBINSIZE{ LOGEMAX / ENERGYBINS };
constexpr double PITCHBINSIZE{ 180 / PITCHBINS };
constexpr double ZTOP{ 4 * 6.371e6 };
constexpr double ZBOTTOM{ 8.371e6 };
constexpr double PI_C{ 3.1415927 };

int main()
{
	if ((ENERGYBINS * PITCHBINS) != NUMPARTICLES)
	{
		std::cout << "Energy Bins * Pitch Bins != Number of Particles : data will be lost, index errors will occur.  Long story short, you don't want to go down this road.  Exiting.\n\n";
		return 0;
	}

	std::cout << "NUMPARTICLES: " << NUMPARTICLES << "\nENERGYBINS: " << ENERGYBINS << "\nPITCHBINS: " << PITCHBINS << "\nLOGEMAX: " << LOGEMAX << "\nLOGEBINSIZE: " << LOGEBINSIZE;
	std::cout << "\nPITCHBINSIZE: " << PITCHBINSIZE << "\nZTOP: " << ZTOP << "\nZBOTTOM: " << ZBOTTOM << "\n";

	//Setup
	//std::vector<double> energies; std::vector<double> pitches;
	//std::vector<double> v_para; std::vector<double> v_perp; std::vector<double> z_top;

	//energies.reserve(ENERGYBINS); pitches.reserve(PITCHBINS);
	//v_para.reserve(NUMPARTICLES); v_perp.reserve(NUMPARTICLES); z_top.reserve(NUMPARTICLES);
	
	double** elec = new double*[4];
	double** ions = new double*[4];
	
	for (int iii = 0; iii < 4; iii++)
	{
		elec[iii] = new double[ENERGYBINS * PITCHBINS];
		ions[iii] = new double[ENERGYBINS * PITCHBINS];
	}

	double* energies = new double[ENERGYBINS];
	double* pitches = new double[PITCHBINS];

	//Setup energies and pitches
	for (int iii = 0; iii < ENERGYBINS; iii++)
		energies[iii] = pow(10, iii * LOGEBINSIZE);
	for (int iii = 0; iii < PITCHBINS; iii++)
		pitches[iii] = iii * 180.0 / (PITCHBINS - 1); //PITCHBINS - 1: to get both 0 and 180 in
		//pitches[iii] = (iii < 256) ? (iii * 0.01 / 256) : (
			//(iii < 768) ? ((iii - 256) * 179.98 / 512 + 0.01) : ((iii - 768) * 0.01 / 256 + 179.98) );

	/*std::cout << energies[0] << "  " << energies[1] << "  " << energies[2] << "  " << energies[3] << "\n";
	std::cout << pitches[0] << "  " << pitches[1] << "  " << pitches[2] << "  " << pitches[3] << "\n\n\n";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[0] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[1] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[2] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[3] * PI_C / 180) << "\n";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[0] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[1] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[2] * PI_C / 180) << "  ";
	std::cout << sqrt(energies[0] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[3] * PI_C / 180) << "\n\n";
	std::cout << cos(pitches[1] * PI_C / 180) << "  " << sin(pitches[1] * PI_C / 180) << "  " << pitches[1] * PI_C / 180 << "\n";*/
	
	//Populate Electron Data
	std::cout << "Populating and writing Electrons...\n";
	for (int iii = 0; iii < PITCHBINS; iii++)
	{
		for (int jjj = 0; jjj < ENERGYBINS; jjj++)
		{
			elec[0][iii * ENERGYBINS + jjj] = -sqrt(energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[iii] * PI_C / 180);
			elec[1][iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[iii] * PI_C / 180);
			elec[2][iii * ENERGYBINS + jjj] = (pitches[iii] <= 90) ? ZTOP : ZBOTTOM;
			//v_para[iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * cos(pitches[iii] * PI / 180);
			//v_perp[iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_ELECTRON) * sin(pitches[iii] * PI / 180);
			//z_top [iii * ENERGYBINS + jjj] = ((iii * ENERGYBINS + jjj) > (NUMPARTICLES / 2)) ? ZTOP : ZBOTTOM;
			//v_para[iii * ENERGYBINS + jjj] *= ((iii * ENERGYBINS + jjj) >(NUMPARTICLES / 2)) ? -1 : 1;
		}
	}

	//double* pointers[4];
	//pointers[0] = v_para.data();
	//pointers[1] = v_perp.data();
	//pointers[2] = z_top.data();

	std::string folderout{ "./../../../out/" };
	std::vector<std::string> names{ "e_vpara.bin", "e_vperp.bin", "e_z.bin" };
	for (int iii = 0; iii < 3; iii++)
	{
		std::string fullname{ folderout + names[iii] };
		fileIO::writeDblBin(fullname.c_str(), elec[iii], NUMPARTICLES);//pointers[iii], NUMPARTICLES);
	}

	//Populate Ion Data
	std::cout << "Populating and writing Ions...\n";
	for (int iii = 0; iii < PITCHBINS; iii++)
	{
		for (int jjj = 0; jjj < ENERGYBINS; jjj++)
		{
			elec[0][iii * ENERGYBINS + jjj] = -sqrt(energies[jjj] * 1.60218e-19 / MASS_PROTON) * cos(pitches[iii] * PI_C / 180);
			elec[1][iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_PROTON) * sin(pitches[iii] * PI_C / 180);
			elec[2][iii * ENERGYBINS + jjj] = (pitches[iii] <= 90) ? ZTOP : ZBOTTOM;
			//v_para[iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_PROTON) * cos(pitches[iii] * PI / 180);
			//v_perp[iii * ENERGYBINS + jjj] = sqrt(energies[jjj] * 1.60218e-19 / MASS_PROTON) * sin(pitches[iii] * PI / 180);
			//z_top [iii * ENERGYBINS + jjj] = ((iii * ENERGYBINS + jjj) >(NUMPARTICLES / 2)) ? ZTOP : ZBOTTOM;
			//v_para[iii * ENERGYBINS + jjj] *= ((iii * ENERGYBINS + jjj) > (NUMPARTICLES / 2)) ? -1 : 1;
		}
	}

	//pointers[0] = v_para.data();
	//pointers[1] = v_perp.data();
	//pointers[2] = z_top.data();

	names = { "i_vpara.bin", "i_vperp.bin", "i_z.bin" };
	for (int iii = 0; iii < 3; iii++)
	{
		std::string fullname{ folderout + names[iii] };
		fileIO::writeDblBin(fullname.c_str(), ions[iii], NUMPARTICLES);//pointers[iii], NUMPARTICLES);
	}
	
	return 0;
}