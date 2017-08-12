#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

//Source: Wikipedia
constexpr double MASS_PROTON  { 1.6726219e-27 }; //kg
constexpr double MASS_ELECTRON{ 9.1093836e-31 }; //kg
constexpr double CHARGE_ELEM  { 1.6021766e-19 }; //C
constexpr double RADIUS_EARTH { 6.371e6 };		 //m
constexpr double BFIELD_EARTH { 32.5e-9 };		 //T (at surface - 1 Re, Wiki mentioned a range from 25-65 nT, fit to dipole model, B0 would be about this)
constexpr double PI { 3.1415927 };

#endif
