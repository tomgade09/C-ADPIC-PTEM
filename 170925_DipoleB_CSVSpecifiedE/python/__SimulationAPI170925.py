from __future__ import absolute_import, division, print_function
import os, sys, inspect
pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(os.path.normpath(pyfiledir + '/../../SimTools/SimulationClass/python/'))
from SimulationAPI import *

def getPointerToElectricFieldLUT(self, simulationptr): #This is not working.  Needs to be fixed in C++ and here.
    simDLL.getPointerToElectricFieldDataWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.getPointerToElectricFieldDataWrapper.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

    #convert to python-easy array
    #return simDLL.getPointerToElectricFieldDataWrapper(simulationptr)
    return
Simulation.getPointerToElectricFieldLUT = getPointerToElectricFieldLUT

if __name__ == '__main__':
    print("SimulationAPI.py is not meant to be called as main.  Run a simulation file and that will import this.")
