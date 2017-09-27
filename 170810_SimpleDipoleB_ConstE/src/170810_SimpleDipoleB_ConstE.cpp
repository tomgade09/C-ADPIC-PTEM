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
#include <chrono>
#include <string>

//Other dependencies

//Project specific includes (mine)
#include "include\_simulationvariables.h" //as a bonus, also includes physicalconstants.h !
#include "include\numericaltools.h"
#include "include\iowrapper.h"

//Defines for making into a DLL File
#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

//Useful snippets
//pitch[iii] = atan2(vperp[iii], vpara[iii]) * 180 / C_PI;

struct retStruct
{
	double*** dblRes{ nullptr };
	bool** inSim{ nullptr };
	double* B_z{ nullptr };
	double* E_z{ nullptr };
	double* B_E_z_dim{ nullptr };
};

retStruct dllmain()
{
	bool loopBool{ true };
	long loopIdx{ 0 };

	double** electrons;
	double** ions;
	bool* e_in_sim = new bool[NUMPARTICLES];
	bool* i_in_sim = new bool[NUMPARTICLES];
	double* B_z = new double[GRAPH_E_B_BINS];
	double* E_z = new double[GRAPH_E_B_BINS];
	double* B_E_z_dim = new double[GRAPH_E_B_BINS];

#ifdef NO_NORMALIZE_M
	std::string unitstring{ " m" };
#else
	std::string unitstring{ " Re" };
#endif

	std::cout << "Sim between:      " << IONSPH_MIN_Z << " - " << MAGSPH_MAX_Z << unitstring << "\n";
	std::cout << "E Field between:  " << (E_RNG_CENTER - E_RNG_DELTA) << " - " << (E_RNG_CENTER + E_RNG_DELTA) << unitstring << "\n";
	std::cout << "Const E:          " << CONSTEFIELD << " V/m\n\n";
	std::cout << "Particle Number:  " << NUMPARTICLES << "\n";
	std::cout << "Iteration Number: " << NUMITERATIONS << "\n";
	std::cout << "Replenish lost p: "; (REPLENISH_E_I) ? (std::cout << "True\n\n") : (std::cout << "False\n\n");

	//don't forget to deallocate memory later... returns array of pointers to arrays of [vpara, vperp, z, null] for particles
	electrons = normalDistribution_v_z(NUMPARTICLES, V_DIST_MEAN, V_SIGMA, Z_DIST_MEAN, Z_SIGMA);
	ions = normalDistribution_v_z(NUMPARTICLES, V_DIST_MEAN, V_SIGMA, Z_DIST_MEAN, Z_SIGMA);

	double*** dbls = new double**[2];
	dbls[0] = electrons;
	dbls[1] = ions;

	writeParticlesToBin(dbls, "./particles_init");

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{//converting vperp (variable) to mu (constant) - only has to be done once
		electrons[1][iii] = 0.5 * MASS_ELECTRON * electrons[1][iii] * electrons[1][iii] / BFieldatZ(electrons[2][iii], 0.0);
		ions[1][iii] = 0.5 * MASS_PROTON * ions[1][iii] * ions[1][iii] / BFieldatZ(ions[2][iii], 0.0);
		e_in_sim[iii] = true;
		i_in_sim[iii] = true;
	}

	std::chrono::steady_clock::time_point cudaBegin, cudaEnd;
	cudaBegin = std::chrono::steady_clock::now();
	
	//CUDA implementation
	mainCUDA(electrons, ions, e_in_sim, i_in_sim, B_z, E_z, B_E_z_dim);

	cudaEnd = std::chrono::steady_clock::now();

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{//converting mu back to vperp
		electrons[1][iii] = sqrt(2 * electrons[1][iii] * BFieldatZ(electrons[2][iii], DT * NUMITERATIONS) / MASS_ELECTRON);
		ions	 [1][iii] = sqrt(2 * ions	  [1][iii] * BFieldatZ(ions		[2][iii], DT * NUMITERATIONS) / MASS_PROTON);
	}

	std::cout << "Parallel Execution Time (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>(cudaEnd - cudaBegin).count() << std::endl;

	writeParticlesToBin(dbls, "./particles_final");

	retStruct ret;

	bool** bools = new bool*[2];
	bools[0] = e_in_sim;
	bools[1] = i_in_sim;

	ret.dblRes = dbls;
	ret.inSim = bools;
	ret.B_z = B_z;
	ret.E_z = E_z;
	ret.B_E_z_dim = B_E_z_dim;

    return ret;
}

DLLEXPORT double* dllmainPyWrapper(char* notusednow)
{
	retStruct yutyut;
	
	yutyut = dllmain();
	
	int e_in_sim{ 0 };
	int i_in_sim{ 0 };
	int eidx{ 1 };
	int iidx{ 2 };

	for (int iii = 0; iii < NUMPARTICLES; iii++)
	{
		if (yutyut.inSim[0][iii])
			e_in_sim++;
		if (yutyut.inSim[1][iii])
			i_in_sim++;
	}

	std::cout << "C++: " << e_in_sim << " " << i_in_sim << " " << ((e_in_sim + i_in_sim) * 3) + 2 << "\n";
	
	//Structure of returned array: [number of electrons remaining in simulation, v_para for electrons, v_perp for electrons, z for electrons,
	//								number of ions remaining in sim, v_para for ions, v_perp for ions, z for ions,
	//                              number of B and E bins, B(z) data, E(z) data, z bins]
	//This is a 1D array which has the number of elements appended before the elements are listed.  So, the number of electrons that have not
	//escaped is listed followed by the three electron properties measured: v_para, v_perp, and z.  This is followed by the ions:
	//number of ions, ion properties: v_para, v_perp, z.  Finally, B and E are measured along z in accordance with the number of bins in the
	//_simulationvariables.h header file.  This is appended to the array followed by the data.  This is done because Python doesn't handle
	//complex C++ structures and such well (at least, it's quite a bit more work).

	double* yut = new double[((e_in_sim + i_in_sim) * 3) + 2 + (GRAPH_E_B_BINS * 3) + 1];

	yut[0] = static_cast<double>(e_in_sim);
	yut[e_in_sim * 3 + 1] = static_cast<double>(i_in_sim);
	yut[((e_in_sim + i_in_sim) * 3) + 2] = static_cast<double>(GRAPH_E_B_BINS);

	for (int iii = 0; iii < 3; iii++)
	{
		for (int jjj = 0; jjj < NUMPARTICLES; jjj++)
		{
			if (yutyut.inSim[0][jjj])
			{
				yut[eidx] = yutyut.dblRes[0][iii][jjj];
				eidx++;
				/*if (eidx == e_in_sim * (iii + 1) + 1)
				{
					std::cout << "Done electron idx " << iii << "\n";
				}*/
				if (eidx > e_in_sim * (iii + 1) + 1)
				{
					std::cout << "Index error (electrons).  Your results are somewhat suspect. " << eidx << "\n";
				}
			}
			if (yutyut.inSim[1][jjj])
			{
				yut[3 * e_in_sim + iidx] = yutyut.dblRes[1][iii][jjj];
				iidx++;
				/*if (iidx == i_in_sim * (iii + 1) + 2)
				{
					std::cout << "Done ion idx " << iii << "\n";
				}*/
				if (iidx > i_in_sim * (iii + 1) + 2)
				{
					std::cout << "Index error (ions).  Your results are somewhat suspect." << iidx << "\n";
				}
			}
		}
	}
	
	for (int iii = 0; iii < GRAPH_E_B_BINS; iii++)
	{
		yut[((e_in_sim + i_in_sim) * 3) + 3 + iii] = yutyut.B_z[iii];
		yut[((e_in_sim + i_in_sim) * 3) + 3 + GRAPH_E_B_BINS + iii] = yutyut.E_z[iii];
		yut[((e_in_sim + i_in_sim) * 3) + 3 + 2 * GRAPH_E_B_BINS + iii] = yutyut.B_E_z_dim[iii];
	}

	return yut;
}

#ifndef DLLFILE
int main()
{//right now, doesn't do anything with the results - just included so the exe will compile and run
	dllmain();
	
	return 0;
}
#endif
