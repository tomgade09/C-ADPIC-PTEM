#ifndef _SIMULATIONVARIABLES_H
#define _SIMULATIONVARIABLES_H
#include <cmath>
#include "physicalconstants.h"

#define NO_NORMALIZE_M
#ifdef NO_NORMALIZE_M
//not normalized
constexpr double NORMFACTOR{ 1.0 };
#else
//meters normalized to Re - i.e. instead of heights being measured in meters, they are measured in Re (from earth's center)
//also other values, such as velocity - in Re / s, not m / s
constexpr double NORMFACTOR{ RADIUS_EARTH };
#endif

//Simulation Variables - Measured from center of earth
constexpr double DT{ 0.01 }; //in s, at v para of 0.15 (T = 2.5) and dt of 0.1, ~415 iterations is a bounce period
constexpr double MIN_Z_SIM{ (2.0e6 + RADIUS_EARTH) / NORMFACTOR }; //Min Z value of the simulation - particles below this are either ignored/removed from sim, or replaced with a new particle
constexpr double MAX_Z_SIM{ 10 * RADIUS_EARTH / NORMFACTOR }; //Max Z value of the simulation - same as above
constexpr int	 NUMPARTICLES{ 100352 }; //number of particles in simulation - best when it's a multiple of 64 (has to be a multiple of BLOCKSIZE)
constexpr long   NUMITERATIONS{ 10000 };
constexpr double INITIAL_T_EV{ 2.5 }; //"2.5" is from Chiu/Schultz temperature (in eV) for the studied plasma, then convert to velocity - kT = 1/2 m v^2
constexpr double INITIAL_T_EV_MAG{ 10 }; //Higher energy magnetospheric particles (Upper limit of Z)
constexpr double T_RATIO{ INITIAL_T_EV_MAG / INITIAL_T_EV }; //velocity proportional to sqrt(T(in eV))
//need one for magnetospheric plasmas
constexpr bool   REPLENISH_E_I{ false }; //determines whether or not to replenish lost electrons/ions - same distribution is used that generates initial characteristics

//E+M Variables
//constexpr double DIPOLETHETA{ 20.0 }; // theta (in deg) - to calculate dipole electric field
constexpr double DIPOLECONST{ BFIELD_EARTH *  1.9102530 };//sqrt(1 + 3 * pow(cos(20.0 * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))
constexpr double TOTPOTDROP{ 2.0e3 }; //total electric potential drop in V across model
constexpr int	 GRAPH_E_B_BINS{ 1000 }; //E, B are measured as a function of z at time 0 and passed out to graph

//CUDA Variables
constexpr int    BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance

//Distribution Variables - for now I use the same values for ions as electrons
constexpr double V_DIST_MEAN{ 0.0 }; //in units of Re / s (if normalized), otherwise m / s
constexpr double VPARACONST{ 1 }; //for creating a more (or less) field aligned beam if desired - multiplied by V_SIGMA in distribution function for v para
const double	 V_SIGMA_ELEC{ sqrt(INITIAL_T_EV * 1.60218e-19 * 2 / MASS_ELECTRON) / NORMFACTOR }; //in units of Re / s  (if normalized), otherwise m / s
constexpr double V_SIGMA_SQ_ELEC{ INITIAL_T_EV * 1.60218e-19 * 2 / (MASS_ELECTRON * NORMFACTOR * NORMFACTOR) }; //need this for cuda code
const double	 V_SIGMA_IONS{ sqrt(INITIAL_T_EV * 1.60218e-19 * 2 / MASS_PROTON) / NORMFACTOR };
constexpr double V_SIGMA_SQ_IONS{ INITIAL_T_EV * 1.60218e-19 * 2 / (MASS_PROTON * NORMFACTOR * NORMFACTOR) };
//Not used//constexpr double Z_DIST_MEAN{ 30.0 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m
//Not used//constexpr double Z_SIGMA{ 0.5 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m

//Functions I can't bring myself to write a header for
double BFieldatZ(double z, double simtime);
double EFieldatZ(double** LUT, double z, double simtime);

#endif