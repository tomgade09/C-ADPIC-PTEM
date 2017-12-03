#ifndef SIMULATIONCLASSEXTENSIONS_H
#define SIMULATIONCLASSEXTENSIONS_H

#include "include\_simulationvariables.h"
#include "include\Simulation.h"
#include "include\fileIO.h"

double BFieldatZ(double z, double simtime);
double EFieldatZ(double** LUT, double z, double simtime);

class Simulation170925 : public Simulation
{
protected:
	std::string LUTfilename_m;
	double**  elcFieldLUT_m{ nullptr };
	
	long totalElecEscaped_m{ 0 };
	long totalIonsEscaped_m{ 0 };
	bool mu_m{ 0 };

	std::vector<double> mass_m;

public:
	Simulation170925(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string rootdir, std::string LUTfilename) :
		Simulation(numberOfParticleTypes, numberOfParticlesPerType, numberOfAttributesTracked, dt, rootdir), LUTfilename_m{ LUTfilename }
	{
		//Generate particles first
		//2 particle types (electrons and ions) with v_para, mu, and z tracked, generate v_para and v_perp (eventually becomming mu) normally
		double means[4];
		double sigmas[4];
		means[0] = V_DIST_MEAN;
		means[1] = V_DIST_MEAN;
		means[2] = V_DIST_MEAN;
		means[3] = V_DIST_MEAN;
		sigmas[0] = V_SIGMA_ELEC;
		sigmas[1] = V_SIGMA_ELEC;
		sigmas[2] = V_SIGMA_IONS;
		sigmas[3] = V_SIGMA_IONS;
		generateNormallyDistributedValues(2, means, sigmas);
		
		//Generate z values and convert v_perp to mu here
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			{
				if (jjj % 2 == 0) //Generate z, every other particle starts at top/bottom of sim respectively
				{
					particles_m[iii][2][jjj] = MIN_Z_SIM + 0.01 * (RADIUS_EARTH / NORMFACTOR);
					particles_m[iii][0][jjj] = abs(particles_m[iii][0][jjj]);
				}
				else
				{
					particles_m[iii][2][jjj] = MAX_Z_SIM - 0.01 * (RADIUS_EARTH / NORMFACTOR);
					particles_m[iii][1][jjj] *= sqrt(T_RATIO); //higher energy magnetosphere particles - still velocity, so sqrt(tratio)
					particles_m[iii][0][jjj] = -abs(particles_m[iii][0][jjj]) * sqrt(T_RATIO);
				}
			}//end for jjj
		}//end for iii
		
		//loadDataFilesIntoParticleArray(); //To load other data (previously saved) into particle arrays

		//Populate mass array
		mass_m.reserve(2);
		mass_m[0] = MASS_ELECTRON;
		mass_m[1] = MASS_PROTON;

		///Test code:
		for (int iii = 0; iii < NUMPARTICLES; iii++)
		{
			if (iii % 2 == 0)
			{
				particles_m[0][0][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR);   //para
				particles_m[0][1][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR); //perp
				particles_m[0][2][iii] = 1.9 * (RADIUS_EARTH / NORMFACTOR); //z
				particles_m[1][0][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR);   //para
				particles_m[1][1][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR); //perp
				particles_m[1][2][iii] = 1.9 * (RADIUS_EARTH / NORMFACTOR); //z
			}
			else
			{
				particles_m[0][0][iii] = -0.15 * (RADIUS_EARTH / NORMFACTOR);
				particles_m[0][1][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR);
				particles_m[0][2][iii] = 2.1 * (RADIUS_EARTH / NORMFACTOR);
				particles_m[1][0][iii] = -0.15 * (RADIUS_EARTH / NORMFACTOR);
				particles_m[1][1][iii] = 0.15 * (RADIUS_EARTH / NORMFACTOR);
				particles_m[1][2][iii] = 2.1 * (RADIUS_EARTH / NORMFACTOR);
			}
		}

		//Populate E Field LUT
		std::string LUT{ rootdir_m + "\\in\\" + LUTfilename_m };
		setElecMagLUT(LUT.c_str(), 2951, 3);

		//Save particle distributions to disk
		std::string fold{ "./particles_init/" };
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}

		convertVPerpToMu();

		timeStructs_m.push_back(createTimeStruct("End Simulation170925 Constructor")); //index 1

		//test test test
		std::cout << "Mu:     " << particles_m[0][1][0] << " " << particles_m[0][1][1] << "\n";
		std::cout << "Mu:     " << particles_m[1][1][0] << " " << particles_m[1][1][1] << "\n";
		std::cout << "B at Z: " << BFieldatZ(1.9 * (RADIUS_EARTH / NORMFACTOR), 0.0) << " " << BFieldatZ(2.1 * (RADIUS_EARTH / NORMFACTOR), 0.0) << "\n";
	}//end constructor

	~Simulation170925() //Destructor
	{
		//Save final particle distributions to disk
		std::string fold{ "./particles_final/" };
		std::vector<std::string> names{ "e_vpara", "e_vperp", "e_z", "i_vpara", "i_vperp", "i_z" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}

		//Delete double arrays
		for (int iii = 0; iii < 3; iii++)
			delete[] elcFieldLUT_m[iii];
		delete[] elcFieldLUT_m;
		
		//Delete satellites
		for (int iii = 0; iii < satellites_m.size(); iii++)
			delete satellites_m[iii];
	}

	//One liners
	double*  getPointerToElectricFieldData(int column) { if (elcFieldLUT_m == nullptr)
		{ std::cout << "Array not initialized yet.  Run Simulation::setElecMagLUT.\n"; return nullptr; } return elcFieldLUT_m[column]; }
	virtual void resetParticlesEscapedCount() { totalElecEscaped_m = 0; totalIonsEscaped_m = 0; return; }

	//Array tools
	virtual void setElecMagLUT(const char* filename, int rows, int cols) { elcFieldLUT_m = fileIO::read2DCSV(filename, rows, cols, ' '); }
	virtual double calculateBFieldAtZandTime(double z, double time) { return BFieldatZ(z, time); }
	virtual double calculateEFieldAtZandTime(double z, double time) { return EFieldatZ(elcFieldLUT_m, z, time); }
	virtual void loadDataFilesIntoParticleArray();
	virtual void convertVPerpToMu();
	virtual void convertMuToVPerp();

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void freeGPUMemory();
	virtual void prepareResults();
};

#endif //end header guard