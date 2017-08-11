#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

const auto pi_c{ 3.141592653589793238463 };	  //unitless
const auto c_c{ 299792458. };				  //[m/s]
const auto mu0_c{ 4e-7 * PI_C };			  //[H/m] or [m kg / (s^2 A^2)]
const auto eps0_c{ 1 / (c_c * c_c * mu0_c) }; //[F/m] or [A^2 s^4 / (m^3 kg)]

#endif