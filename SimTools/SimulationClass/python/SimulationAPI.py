from __future__ import absolute_import, division, print_function

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

import ctypes

dllLocation = './../vs/build/SimulationClass.dll'

simDLL = ctypes.CDLL(dllLocation)

#From Simulation
#One liner functions
def getSimulationTimeCPP(simulationptr):
    simDLL.getSimulationTimeWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.getSimulationTimeWrapper.restype = ctypes.c_double

    return simDLL.getSimulationTimeWrapper(simulationptr)

def incrementSimulationTimeByDtCPP(simulationptr):
    simDLL.incrementSimulationTimeByDtWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.incrementSimulationTimeByDtWrapper.restype = None

    simDLL.incrementSimulationTimeByDtWrapper(simulationptr)

def resetParticlesEscapedCountCPP(simulationptr):
    simDLL.resetParticlesEscapedCountWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.resetParticlesEscapedCountWrapper.restype = None

    simDLL.resetParticlesEscapedCountWrapper(simulationptr)


#Pointer one liners
def getPointerTo3DParticleArrayCPP(simulationptr): #Test to see if this is working - will python recognize a triple pointer?
    simDLL.getPointerTo3DParticleArrayWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.getPointerTo3DParticleArrayWrapper.restype = ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))

    #convert to python-easy array
    return simDLL.getPointerTo3DParticleArrayWrapper(simulationptr)

def getPointerToSingleParticleTypeArrayCPP(simulationptr, index): #Test to see if this is working - will python recognize a double pointer?
    simDLL.getPointerToSingleParticleTypeArrayWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int)
    simDLL.getPointerToSingleParticleTypeArrayWrapper.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

    #convert to python-easy array
    return simDLL.getPointerToSingleParticleTypeArrayWrapper(simulationptr, index)

def getPointerToSerializedParticleArrayCPP(simulationptr): #Test to see if this is working - will python recognize a double pointer?
    simDLL.getPointerToSerializedParticleArrayWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.getPointerToSerializedParticleArrayWrapper.restype = ctypes.POINTER(ctypes.c_double)

    #convert to python-easy array
    return simDLL.getPointerToSerializedParticleArrayWrapper(simulationptr)


#Numerical tools
def generateNormallyDistributedCPP(simulationptr, numberOfNormalAttributes, means, sigmas): #Test to see if this is working - will python recognize a double pointer?
    simDLL.generateNormallyDistributedValues.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
    simDLL.generateNormallyDistributedValues.restype = None

    C_DOUBLEA = ctypes.c_double * numberOfNormalAttributes
    meanscpp = C_DOUBLEA()
    sigmascpp = C_DOUBLEA()

    for iii in range(numberOfNormalAttributes):
        meanscpp[iii] = means[iii]
        sigmascpp[iii] = sigmas[iii]

    simDLL.generateNormallyDistributedValues(simulationptr, numberOfNormalAttributes, meanscpp, sigmascpp)

def calculateMeanOfParticleAttributeCPP(simulationptr, particleIndex, attributeIndex, absValueBool=False):
    simDLL.calculateMeanOfParticleAttributeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_bool)
    simDLL.calculateMeanOfParticleAttributeWrapper.restype = ctypes.c_double

    return simDLL.calculateMeanOfParticleAttributeWrapper(simulationptr, particleIndex, attributeIndex, absValueBool)

def calculateStdDevOfParticleAttributeCPP(simulationptr, particleIndex, attributeIndex):
    simDLL.calculateStdDevOfParticleAttributeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
    simDLL.calculateStdDevOfParticleAttributeWrapper.restype = ctypes.c_double

    return simDLL.calculateStdDevOfParticleAttributeWrapper(simulationptr, particleIndex, attributeIndex)


#Array tools
def serializeParticleArrayCPP(simulationptr):
    simDLL.serializeParticleArrayWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.serializeParticleArrayWrapper.restype = None

    simDLL.serializeParticleArrayWrapper(simulationptr)

def calculateBFieldAtZandTimeCPP(simulationptr, z, time):
    simDLL.calculateBFieldAtZandTimeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
    simDLL.calculateBFieldAtZandTimeWrapper.restype = ctypes.c_double

    return simDLL.calculateBFieldAtZandTimeWrapper(simulationptr, z, time)

def calculateEFieldAtZandTimeCPP(simulationptr, z, time):
    simDLL.calculateEFieldAtZandTimeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
    simDLL.calculateEFieldAtZandTimeWrapper.restype = ctypes.c_double

    return simDLL.calculateEFieldAtZandTimeWrapper(simulationptr, z, time)

#Simulation Management Function Wrappers
def initializeSimCPP(simulationptr):
    simDLL.initializeWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.initializeWrapper.restype = None

    simDLL.initializeWrapper(simulationptr)

def copyDataToGPUCPP(simulationptr):
    simDLL.copyDataToGPUWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.copyDataToGPUWrapper.restype = None

    simDLL.copyDataToGPUWrapper(simulationptr)

def iterateSimulationCPP(simulationptr, numberOfIterations):
    simDLL.iterateSimulationWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int)
    simDLL.iterateSimulationWrapper.restype = None

    simDLL.iterateSimulationWrapper(simulationptr, numberOfIterations)

def copyDataToHostCPP(simulationptr):
    simDLL.copyDataToHostWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.copyDataToHostWrapper.restype = None

    simDLL.copyDataToHostWrapper(simulationptr)

def terminateSimulationWrapper(simulationptr):
    simDLL.terminateSimulationWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.terminateSimulationWrapper.restype = None

    simDLL.terminateSimulationWrapper(simulationptr)

def returnResultsWrapper(simulationptr):
    simDLL.returnResultsWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.returnResultsWrapper.restype = ctypes.POINTER(ctypes.c_double)

    return simDLL.returnResultsWrapper(simulationptr)


#From Simulation170925
def getPointerToElectricFieldDataCPP(simulationptr): #Test to see if this is working - will python recognize a double pointer?
    simDLL.getPointerToElectricFieldDataWrapper.argtypes = (ctypes.c_void_p,)
    simDLL.getPointerToElectricFieldDataWrapper.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

    #convert to python-easy array
    return simDLL.getPointerToElectricFieldDataWrapper(simulationptr)

#def getPointerToMagneticFieldDataCPP(simulationptr): #Test to see if this is working - will python recognize a double pointer?
    #simDLL.getPointerToMagneticFieldDataWrapper.argtypes = (ctypes.c_void_p,)
    #simDLL.getPointerToMagneticFieldDataWrapper.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))

    #convert to python-easy array
    #return simDLL.getPointerToMagneticFieldDataWrapper(simulationptr)


if __name__ == '__main__':
    print("SimulationAPI.py is not meant to be called as main.  Run a simulation file and that will import this automatically.")
