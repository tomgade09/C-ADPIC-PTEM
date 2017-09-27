#ifndef _SIMULATIONVARIABLES_H
#define _SIMULATIONVARIABLES_H
#include <cmath>
#include "physicalconstants.h"

//#define NO_NORMALIZE_M

#ifdef NO_NORMALIZE_M
//not normalized - not exactly working right now...need to look at, for now just leave the NO_NORMALIZE_M commented out
constexpr double NORMFACTOR{ 1.0 };
#else
//meters normalized to Re - i.e. instead of heights being measured in meters, they are measured in Re (from earth's center)
//also other values, such as velocity - in Re / s, not m / s
constexpr double NORMFACTOR{ RADIUS_EARTH };
#endif

//Simulation Variables - Measured from center of earth
constexpr double DT{ 0.01 }; //in s, at v para of 0.15 (T = 2.5) and dt of 0.1, ~415 iterations is a bounce period
constexpr double MIN_Z_FROM_RE{ 2.0e6 + RADIUS_EARTH }; //in m - how far up from earth's core is minimum height for sim
constexpr double MAX_Z_FROM_RE{ 10 * RADIUS_EARTH }; //in m - how far up is max height for sim
constexpr double IONSPH_MIN_Z{ MIN_Z_FROM_RE / NORMFACTOR }; //normalized to Re - don't change these
constexpr double MAGSPH_MAX_Z{ MAX_Z_FROM_RE / NORMFACTOR }; //normalized to Re - don't change these
constexpr int	 NUMPARTICLES{ 100352 }; //number of particles in simulation - best when it's a multiple of 64 (has to be a multiple of BLOCKSIZE)
constexpr long   NUMITERATIONS{ 10000 };
constexpr double INITIAL_T_EV{ 2.5 }; //magical "2.5" is from Chiu/Schultz temperature (in eV) for the studied plasma, then convert to velocity - kT = 1/2 m v^2
constexpr bool   REPLENISH_E_I{ false }; //determines whether or not to replenish lost electrons/ions - same distribution is used that generates initial characteristics

//E+M Variables
//constexpr double DIPOLETHETA{ 20.0 }; // theta (in deg) - to calculate dipole electric field
constexpr double DIPOLECONST{ BFIELD_EARTH *  1.9102530 };//sqrt(1 + 3 * pow(cos(20.0 * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))
constexpr double TOTPOTDROP{ 2.0e3 }; //total electric potential drop in V across model
constexpr int	 GRAPH_E_B_BINS{ 1000 }; //E, B are measured as a function of z at time 0 and passed out to graph
//E Field centered at 2Re, between about another 1000 km +/- the center
constexpr double E_RNG_CENTER{ (2 * RADIUS_EARTH) / NORMFACTOR }; //Where is the E Field centered?
constexpr double E_RNG_DELTA{ 1.0e6 / NORMFACTOR }; //in m, How far up and down from the center will the E field be "felt"?
constexpr double CONSTEFIELD{ TOTPOTDROP / (2 * E_RNG_DELTA * NORMFACTOR) }; //E Field centered at 2 Re, spread out 2000 km, V / m
//Const E Field across whole model (comment out three lines above, uncomment three lines below)
//constexpr double E_RNG_CENTER{ (MAGSPH_MAX_Z + IONSPH_MIN_Z) / 2 }; //Where is the E Field centered?
//constexpr double E_RNG_DELTA{ (MAGSPH_MAX_Z - IONSPH_MIN_Z) / 2 }; //How far up and down from the center will the E field be "felt"
//constexpr double CONSTEFIELD{ TOTPOTDROP / (MAX_Z_FROM_RE - MIN_Z_FROM_RE) }; //const E over whole sim range, V / m

//CUDA Variables
constexpr int    BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance

//Distribution Variables - for now I use the same values for ions as electrons
constexpr double Z_DIST_MEAN{ 30.0 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m
constexpr double Z_SIGMA{ 0.5 * (RADIUS_EARTH / NORMFACTOR) }; //in units of Re (if normalized), otherwise m
constexpr double V_DIST_MEAN{ 0.0 }; //in units of Re / s (if normalized), otherwise m / s
constexpr double VPARACONST{ 1 }; //for creating a more (or less) field aligned beam if desired - multiplied by V_SIGMA in distribution function for v para
const double	 V_SIGMA{ sqrt(INITIAL_T_EV * 1.60218e-19 * 2 / MASS_ELECTRON) / NORMFACTOR }; //in units of Re / s  (if normalized), otherwise m / s
constexpr double V_SIGMA_SQ{ INITIAL_T_EV * 1.60218e-19 * 2 / (MASS_ELECTRON * NORMFACTOR * NORMFACTOR) }; //need this for cuda code

//Functions I can't bring myself to write a header for
double BFieldatZ(double z, double simtime);
double EFieldatZ(double z, double simtime);

#endif