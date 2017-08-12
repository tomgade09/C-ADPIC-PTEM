#ifndef _SIMULATIONVARIABLES_H
#define _SIMULATIONVARIABLES_H
#include <cmath>
#include "physicalconstants.h"

constexpr double DT{ 0.1 }; //at approx v para of 0.15 and dt of 0.1, ~415 iterations is a bounce period
constexpr double TOTPOTDROP{ 1.0e3 }; //total electric potential drop in V across model
constexpr double DIPOLETHETA{ 20.0 }; // theta (in rad) - to calculate dipole electric field
constexpr double MIN_Z_FROM_RE{ 2.0e6 + RADIUS_EARTH }; //in m - how far up from earth's core is minimum height for sim
constexpr double MAX_Z_FROM_RE{ 2.0e8 + RADIUS_EARTH }; //in m - how far up is max height for sim
constexpr int NUMPARTICLES{ 100000 }; //number of particles in simulation
constexpr unsigned long int NUMITERATIONS{ 10000 };
constexpr double Z_DIST_MEAN{ 30.0 }; //in units of Re
constexpr double Z_SIGMA{ 0.5 }; //in units of Re
constexpr double V_DIST_MEAN{ 0.0 }; //in units of Re / s
constexpr double INITIAL_T_EV{ 2.5 }; //magical "2.5" is from Chiu/Schultz temperature (in eV) for the studied plasma, then convert to velocity - kT = 1/2 m v^2
constexpr double PARACONST{ 100 }; //for creating a more (or less) field aligned beam if desired - multiplied by V_SIGMA in distribution function

//Shouldn't need to change these, but they are still global const variables
constexpr double IONSPH_MIN_Z{ MIN_Z_FROM_RE / RADIUS_EARTH }; //normalized to Re
constexpr double MAGSPH_MAX_Z{ MAX_Z_FROM_RE / RADIUS_EARTH }; //normalized to Re
constexpr double CONSTEFIELD{ TOTPOTDROP / (MAX_Z_FROM_RE - MIN_Z_FROM_RE) }; //V / m
constexpr double DIPOLECONST{ BFIELD_EARTH *  1.910253 };//sqrt(1 + 3 * pow(cos(DIPOLETHETA * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))
const double V_SIGMA{ sqrt(INITIAL_T_EV * 1.60218e-19 * 2 / MASS_ELECTRON) / RADIUS_EARTH }; //in units of Re / s - for some reason, can't use constexpr with sqrt
																							 //for now I use the same values for ions as electrons
//Functions I can't bring myself to write a header for
double accel1DCB(double* args, int len);
double BFieldatZ(double z);
void mainCUDA(double** electrons, double** ions, bool* elec_in_sim_host, bool* ions_in_sim_host);

#endif