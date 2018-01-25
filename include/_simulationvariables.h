#ifndef _SIMULATIONVARIABLES_H
#define _SIMULATIONVARIABLES_H
#include <cmath>
#include "physicalconstants.h"

//Simulation Variables - Measured from center of earth in meters [/second] unless otherwise stated
//constexpr double DT{ 0.01 }; //in s, at v para of 0.15 (T = 2.5) and dt of 0.1, ~415 iterations is a bounce period
//constexpr double MIN_Z_SIM{ (2.0e6 + RADIUS_EARTH) }; //Min Z value of the simulation
//constexpr double MAX_Z_SIM{ 4 * RADIUS_EARTH }; //Max Z value of the simulation
//constexpr double MIN_Z_NORM{ MIN_Z_SIM / RADIUS_EARTH };
//constexpr double MAX_Z_NORM{ MAX_Z_SIM / RADIUS_EARTH };
//constexpr int	 NUMPARTICLES{ 100352 }; //number of particles in simulation - best when it's a multiple of 64 (has to be a multiple of BLOCKSIZE)
//constexpr double INITIAL_T_EV{ 2.5 }; //"2.5" is from Chiu/Schultz temperature (in eV) for the studied plasma, then convert to velocity - kT = 1/2 m v^2
//constexpr double INITIAL_T_EV_MAG{ 1000 }; //Higher energy magnetospheric particles (Upper limit of Z)
//constexpr double CONSTEFIELD{ 2e3 / (MAX_Z_SIM - MIN_Z_SIM) }; //strength of constant E Field (V / m)

//Distribution Variables - for now I use the same values for ions as electrons
//constexpr double V_DIST_MEAN{ 0.0 }; //mean of the velocity distribution
//constexpr double V_SIGMA_SQ_ELEC{ INITIAL_T_EV * JOULE_PER_EV / (MASS_ELECTRON) }; //sigma^2 of the velocity distribution - due to two maxwellians (para, perp), energy is 2 * what it should be, hence / 2
//constexpr double V_SIGMA_SQ_IONS{ INITIAL_T_EV * JOULE_PER_EV / (MASS_PROTON) };   //sigma^2 of the velocity distribution
//constexpr double T_RATIO{ INITIAL_T_EV_MAG / INITIAL_T_EV }; //velocity proportional to sqrt(T(in eV))
//divide sigma_sq by 2 - vpara^2 + vperp^2 = vtot^2, avg vpara^2 = avg vperp^2, avg vtot^2 = 2vpara^2 = 2vperp^2


//E+M Variables
//constexpr double DIPOLETHETA{ 20.0 }; // theta (in deg) - to calculate dipole electric field
//constexpr double B0ATTHETA{ BFIELD_EARTH *  1.9102530 };//sqrt(1 + 3 * pow(cos(20.0 * PI / 180),2)) }; //B0 * sqrt(1 + 3*cos^2(theta))

//CUDA Variables - if you change these, don't forget to change the associated curand code/blocks/etc
// For Geforce 960M (author's computer) - maximum 1024 threads per block - try this to see if it results in faster code execution sometime
//constexpr int partscount{ 100532 };
//constexpr int    BLOCKSIZE{ 256 }; //Number of threads per block - this is most efficient at a multiple of 128 (256 seems to work well), although 250 has been used with slightly less performance
//constexpr int	 NUMBLOCKS{ partscount / BLOCKSIZE }; //Number of blocks
//constexpr int	 NUMRNGSTATES{ 64 * BLOCKSIZE };

//constexpr int    LOOPS_BTW_PROGRESS_COUT{ 500 };

//Functions I can't bring myself to write a header for
//double BFieldatZ(double z, double simtime);
//double EFieldatZ(double** LUT, double z, double simtime, double omegaE, double constE, bool qsps, bool alfven);

#endif