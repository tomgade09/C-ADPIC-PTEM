# SIMULATION VARIABLES
# Set variables related to the simulation here

import os, inspect

#Directories
PYROOTDIR = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
ROOTDIR = os.path.abspath(PYROOTDIR + './../')
LIBDIR = os.path.abspath(ROOTDIR + './lib/')
DLLLOCATION = os.path.abspath(LIBDIR + './geoplasmasim.dll')

#Physical Constants
RADIUS_EARTH = 6.3712e6
J_PER_EV = 1.6021766209e-19
MASS_ELEC = 9.10938356e-31
MASS_PROT = 1.672621898e-27
CHARGE_ELEM = 1.6021766209e-19
PI = 3.14159265358979323846

#Simulation Specific Constants
DT = 0.01
MIN_S_SIM = 2030837.49610366 #2.0e6
MAX_S_SIM = 19881647.2473464 #3.0 * RADIUS_EARTH
MIN_S_NORM = MIN_S_SIM / RADIUS_EARTH
MAX_S_NORM = MAX_S_SIM / RADIUS_EARTH
NUMITER = 25000

#Distribution Quantities
RANDOM = True
LOADDIST = False
DISTINFOLDER = os.path.abspath(ROOTDIR + './_in/data/')

NUMPARTICLES = 100352
INITIAL_T_ION_EV = 2.5
INITIAL_T_MAG_EV = 1000

#B Field Model
DIPOLE = True
IGRF = False

ILAT = 72.0

#Satellites Data


#QSPS E


#Alfven E LUT
ALFVENLUTCSV = ""
#LUTDATAFILE = os.path.abspath(ROOTDIR + './_in/ez.out')