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
	//double**  magFieldLUT_m{ nullptr }; //not needed for now
public:
	Simulation170925(int numberOfParticleTypes, int numberOfParticlesPerType, int numberOfAttributesTracked, double dt, std::string LUTfilename = "") :
		Simulation(numberOfParticleTypes, numberOfParticlesPerType, numberOfAttributesTracked, dt), LUTfilename_m{ LUTfilename } 
	{
		//Generate particles first
		//2 particle types (electrons and ions) with v_para, mu, and z tracked, generate v_para and v_perp (eventually becomming mu) normally
		double means[4];
		double sigmas[4];
		means[0] = V_DIST_MEAN;
		means[1] = V_DIST_MEAN;
		means[2] = V_DIST_MEAN;
		means[3] = V_DIST_MEAN;
		sigmas[0] = V_SIGMA;
		sigmas[1] = V_SIGMA;
		sigmas[2] = V_SIGMA;
		sigmas[3] = V_SIGMA;
		generateNormallyDistributedValues(2, means, sigmas);

		//Generate z values and convert v_perp to mu here
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			{
				if (jjj % 2 == 0) //Generate z, every other particle starts at top/bottom of sim respectively
				{
					particles_m[iii][2][jjj] = IONSPH_MIN_Z + 0.01;
					particles_m[iii][0][jjj] = abs(particles_m[iii][0][jjj]);
				}
				else
				{
					particles_m[iii][2][jjj] = MAGSPH_MAX_Z + 0.01;
					particles_m[iii][0][jjj] = -abs(particles_m[iii][0][jjj]);
				}
				if (iii == 0) //Convert v_perp to mu
					particles_m[iii][1][jjj] = 0.5 * MASS_ELECTRON * particles_m[iii][1][jjj] * particles_m[iii][1][jjj] / BFieldatZ(particles_m[iii][2][jjj], 0.0);
				else
					particles_m[iii][1][jjj] = 0.5 * MASS_PROTON * particles_m[iii][1][jjj] * particles_m[iii][1][jjj] / BFieldatZ(particles_m[iii][2][jjj], 0.0);
			}//end for jjj
		}//end for iii

		std::string fold;
		fold = "./particles_init/";
		std::vector<std::string> names;
		names.reserve(numberOfParticleTypes_m * numberOfAttributesTracked_m);
		names = { "v_e_para", "mu_e", "z_e", "v_i_para", "mu_i", "z_i" };
		for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
		{
			for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
				saveParticleAttributeToDisk(iii, jjj, fold.c_str(), names[iii * numberOfAttributesTracked_m + jjj].c_str());
		}

		//setElecMagLUT(LUTfilename_m.c_str(), 2951, 3);
	}//end constructor
	~Simulation170925() {}

	//One liners
	double**  getPointerToElectricFieldData() { if (particlesSerialized_m == nullptr) { std::cout << "Array not serialized yet.  Run Simulation::setElecMagLUT.\n"; } return elcFieldLUT_m; }
	//double**  getPointerToMagneticFieldData() { if (particlesSerialized_m == nullptr) { std::cout << "Array not serialized yet.  Run Simulation::setElecMagLUT.\n"; } return magFieldLUT_m; }

	//Array tools
	virtual void setElecMagLUT(const char* filename, int rows, int cols);
	virtual double calculateBFieldAtZandTime(double z, double time);
	virtual double calculateEFieldAtZandTime(double z, double time);

	//Simulation management functions
	virtual void initializeSimulation();
	virtual void copyDataToGPU();
	virtual void iterateSimulation(int numberOfIterations);
	virtual void copyDataToHost();
	virtual void terminateSimulation(); //name is not great - maybe change - why would the returnResults function have to be called after the sim is terminated?  A function with a name like this needs to close everything down
	virtual double* returnResults();
};

#endif //end header guard