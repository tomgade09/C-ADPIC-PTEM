// 170810_SimpleDipoleB_ConstE.cpp
// Using Dipole B field, constant E that leads to a potential drop of 1kV over the range evaluated (~ 2 - 200 km above earths surface)
// Distances normalized to Re
// B = (B0 / r^3) * sqrt(1 + 3 * cos^2(theta)) - theta = 20 deg
// Tracking one position dimension for the particle: z, positive is up
// Tracking v_para, mu (const) for each particle -> mu provides v_perp as a function of B
// mu = 1/2 m v_perp^2 / B --- v_perp = sqrt(2 mu B / m)

// C std lib includes
#include <cmath>
#include <iostream>

//Other dependencies

//Project specific includes (mine)
#include "include\_simulationvariables.h" //as a bonus, also includes physicalconstants.h !
#include "include\numericaltools.h"

//Defines for making into a DLL File
#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

//Useful snippets
//pitch[iii] = atan2(vperp[iii], vpara[iii]) * 180 / C_PI;

double*** dllmain()
{
	bool loopBool{ true };
	unsigned long int loopIdx{ 0 };

	double** h_electrons;
	double** h_ions;
	bool* h_eInPlay = new bool[NUMPARTICLES];
	bool* h_iInPlay = new bool[NUMPARTICLES];

	//don't forget to deallocate memory later... returns array of pointers to arrays of [vpara, vperp, z, null] for particles
	h_electrons = normalDistribution_v_z(NUMPARTICLES, V_DIST_MEAN, V_SIGMA, Z_DIST_MEAN, Z_SIGMA);
	h_ions = normalDistribution_v_z(NUMPARTICLES, V_DIST_MEAN, V_SIGMA, Z_DIST_MEAN, Z_SIGMA);

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{//converting vperp (variable) to mu (constant) - only has to be done once
		h_electrons[1][iii] = 0.5 * MASS_ELECTRON * h_electrons[1][iii] * h_electrons[1][iii] / BFieldatZ(h_electrons[2][iii]);
		h_ions[1][iii] = 0.5 * MASS_PROTON * h_ions[1][iii] * h_ions[1][iii] / BFieldatZ(h_ions[2][iii]);
	}

	double accel{ 0.0 };
	double args[6];

	while (loopIdx < NUMITERATIONS)
	{
		for (int iii = 0; iii < NUMPARTICLES; iii++) //run every iteration
		{
			h_eInPlay[iii] = !((h_electrons[2][iii] > MAGSPH_MAX_Z) || (h_electrons[2][iii] < IONSPH_MIN_Z)); //Makes sure particles are within bounds
			h_iInPlay[iii] = !((h_ions[2][iii] > MAGSPH_MAX_Z) || (h_ions[2][iii] < IONSPH_MIN_Z));

			/*if (h_eInPlay[iii] == false)
				std::cout << "Electron outside of range: " << iii << "  " << h_electrons[2][iii] << " outside of " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << "\n";
			if (h_iInPlay[iii] == false)
				std::cout << "Ion outside of range: " << iii << "  " << h_ions[2][iii] << " outside of " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << "\n";*/
		}

		for (int iii = 0; iii < NUMPARTICLES; iii++) //updates position of every particle
		{//args array: [dt, vz, mu, q, m, pz_0]
			//electrons
			if (h_eInPlay[iii] == true)
			{//if electrons leave the simulation, stop doing calculations on them
				args[0] = 0.0;
				args[1] = h_electrons[0][iii];
				args[2] = h_electrons[1][iii];
				args[3] = CHARGE_ELEM * -1;
				args[4] = MASS_ELECTRON;
				args[5] = h_electrons[2][iii];
				accel = fourthOrderRungeKutta1D(accel1DCB, args, 6, DT);
				h_electrons[0][iii] += accel * DT;
			}

			//ions
			if (h_iInPlay[iii] == true)
			{
				args[1] = h_ions[0][iii];
				args[2] = h_ions[1][iii];
				args[3] = CHARGE_ELEM;
				args[4] = MASS_PROTON;
				args[5] = h_ions[2][iii];
				accel = fourthOrderRungeKutta1D(accel1DCB, args, 6, DT);
				h_ions[0][iii] += accel * DT;
			}
		}
		
		loopIdx++;
		if (loopIdx % 100 == 0)
			std::cout << loopIdx << "\n";
	}

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{//converting mu back to vperp
		h_electrons[1][iii] = sqrt(2 * h_electrons[1][iii] * BFieldatZ(h_electrons[2][iii]) / MASS_ELECTRON);
		h_ions[1][iii] = sqrt(2 * h_ions[1][iii] * BFieldatZ(h_ions[2][iii]) / MASS_PROTON);
	}

	double*** ret = new double**[2];
	ret[0] = h_electrons;
	ret[1] = h_ions;

	std::cout << "checking ret (pointer to pointer to pointer): " << ret[0][0][0] << ret[1][0][0] << "\n";
	std::cout << "checking ret (pointer to pointer to pointer): " << ret[0][0][NUMPARTICLES - 1] << ret[1][0][NUMPARTICLES - 1] << "\n";

	int ionCount{ 0 };
	int electronCount{ 0 };

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{
		if (h_eInPlay[iii] == true)
			electronCount++;
		if (h_iInPlay[iii] == true)
			ionCount++;
	}

	std::cout << "electrons left: " << electronCount << "  ions left: " << ionCount << "\n";

    return ret;
}

DLLEXPORT double* dllmainPyWrapper()
{
	double* yut = new double[NUMPARTICLES * 3 * 2];
	double*** yutyut{ nullptr };

	yutyut = dllmain();

	for (int iii = 0; iii < 2; iii++)
	{
		for (int jjj = 0; jjj < 3; jjj++)
		{
			for (int kkk = 0; kkk < NUMPARTICLES; kkk++)
			{
				yut[iii * 3 * NUMPARTICLES + jjj * NUMPARTICLES + kkk] = yutyut[iii][jjj][kkk];
			}
			delete[] yutyut[iii][jjj];
		}
		delete[] yutyut[iii];
	}
	delete[] yutyut;

	return yut;
}
	///Print particle characteristics
	/*for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
	{
		std::cout << "vpara, vperp, z, InPlay? : " << h_electrons[0][jjj] << ", " << h_electrons[1][jjj] << ", ";
		std::cout << h_electrons[2][jjj] << ", " << h_electrons[3][jjj] << " \n";
	}

	for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
	{
		std::cout << "vpara, vperp, z, InPlay? : " << h_ions[0][jjj] << ", " << h_ions[1][jjj] << ", ";
		std::cout << h_ions[2][jjj] << ", " << h_ions[3][jjj] << " \n";
	}*/

	///FOR CODE THAT REPLACES BOOL WITH DOUBLE AND APPENDS TO THE END OF H_ELEC/H_IONS - FROM MAIN
	//double* h_eInPlay = new double[NUMPARTICLES]; //positive - particle is in sim, negative - it isn't
	//double* h_iInPlay = new double[NUMPARTICLES]; //making this double instead of bool will allow me to "bundle" everything together without complicated structs, etc
	//h_electrons[3] = h_eInPlay;
	//h_ions[3] = h_iInPlay;

	/*for (int iii = 0; iii < NUMPARTICLES; iii++) //May want to do this differently later
	{
		if ((h_electrons[2][iii] > MAGSPH_MAX_Z) || (h_electrons[2][iii] < IONSPH_MIN_Z))
		{
			h_eInPlay[iii] = -10.0;
			std::cout << "Electrons outside of range: " << iii << "  " << h_electrons[2][iii] << " outside of " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << "\n";
		}
		else
			h_eInPlay[iii] = 10.0;
	}

	for (int iii = 0; iii < NUMPARTICLES; iii++) //May want to do this differently later
	{
		if ((h_ions[2][iii] > MAGSPH_MAX_Z) || (h_ions[2][iii] < IONSPH_MIN_Z))
		{
			h_iInPlay[iii] = -10.0;
			std::cout << "Ions outside of range: " << iii << "  " << h_ions[2][iii] << " outside of " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << "\n";
		}
		else
			h_iInPlay[iii] = 10.0;
	}*/

	/*for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
	{
		std::cout << "vpara, vperp, z, InPlay? : " << h_electrons[0][jjj] << ", " << h_electrons[1][jjj] << ", ";
		std::cout << h_electrons[2][jjj] << ", " << h_electrons[3][jjj] << " \n";
	}

	for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
	{
		std::cout << "vpara, vperp, z, InPlay? : " << h_ions[0][jjj] << ", " << h_ions[1][jjj] << ", ";
		std::cout << h_ions[2][jjj] << ", " << h_ions[3][jjj] << " \n";
	}*/