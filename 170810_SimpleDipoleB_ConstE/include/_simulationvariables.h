#ifndef _SIMULATIONVARIABLES_H
#define _SIMULATIONVARIABLES_H
#include <cmath>
#include "physicalconstants.h"
//#define NO_NORMALIZE_M

#ifdef NO_NORMALIZE_M
//not normalized
constexpr double NORMFACTOR{ 1.0 };
#else
//meters normalized to Re - i.e. instead of heights being measured in meters, they are measured in Re (from earth's center)
//also other values, such as velocity - in Re / s, not m / s
constexpr double NORMFACTOR{ RADIUS_EARTH };
#endif

//Simulation Variables
constexpr double DT{ 0.1 }; //in s, at v para of 0.15 (T = 2.5) and dt of 0.1, ~415 iterations is a bounce period
constexpr double MIN_Z_FROM_RE{ 2.0e6 + RADIUS_EARTH }; //in m - how far up from earth's core is minimum height for sim
constexpr double MAX_Z_FROM_RE{ 2.0e8 + RADIUS_EARTH }; //in m - how far up is max height for sim
constexpr int	 NUMPARTICLES{ 100000 }; //number of particles in simulation - needs to be an even number
constexpr unsigned long int NUMITERATIONS{ 1000000 };
constexpr double TOTPOTDROP{ 1.0e3 }; //total electric potential drop in V across model
constexpr double DIPOLETHETA{ 20.0 }; // theta (in rad) - to calculate dipole electric field
constexpr double INITIAL_T_EV{ 2.5 }; //magical "2.5" is from Chiu/Schultz temperature (in eV) for the studied plasma, then convert to velocity - kT = 1/2 m v^2
constexpr bool REPLENISH_E_I{ false }; //determines whether or not to replenish lost electrons/ions - same distribution is used that generates initial characteristics

//Distribution Variables - for now I use the same values for ions as electrons
constexpr double Z_DIST_MEAN{ 30.0 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m
constexpr double Z_SIGMA{ 0.5 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m
constexpr double V_DIST_MEAN{ 0.0 }; //in units of Re / s (if normalized), otherwise m / s
constexpr double VPARACONST{ 1 }; //for creating a more (or less) field aligned beam if desired - multiplied by V_SIGMA in distribution function for v para
const double V_SIGMA{ sqrt(INITIAL_T_EV * 1.60218e-19 * 2 / MASS_ELECTRON) / NORMFACTOR }; //in units of Re / s  (if normalized), otherwise m / s
constexpr double V_SIGMA_SQ{ INITIAL_T_EV * 1.60218e-19 * 2 / (MASS_ELECTRON * NORMFACTOR * NORMFACTOR) }; //need this for cuda code

//Shouldn't need to change these, but they are still global const variables
constexpr double IONSPH_MIN_Z{ MIN_Z_FROM_RE / NORMFACTOR }; //normalized to Re
constexpr double MAGSPH_MAX_Z{ MAX_Z_FROM_RE / NORMFACTOR }; //normalized to Re
constexpr double CONSTEFIELD{ TOTPOTDROP / (MAX_Z_FROM_RE - MIN_Z_FROM_RE) }; //V / m
constexpr double DIPOLECONST{ BFIELD_EARTH *  1.910253 };//sqrt(1 + 3 * pow(cos(DIPOLETHETA * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))

//Functions I can't bring myself to write a header for
double accel1DCB(double* args, int len);
double BFieldatZ(double z);
void mainCUDA(double** electrons, double** ions, bool* elec_in_sim_host, bool* ions_in_sim_host);

#endif