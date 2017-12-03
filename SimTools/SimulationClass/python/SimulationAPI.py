from __future__ import absolute_import, division, print_function
import ctypes

class Simulation:
    def __init__(self, rootdir, DLLloc):
        self.dllLoc_m = DLLloc
        self.rootdir_m = rootdir
        self.simDLL_m = ctypes.CDLL(self.dllLoc_m)

        #One liner functions
        self.simDLL_m.getSimulationTimeAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimulationTimeAPI.restype = ctypes.c_double
        self.simDLL_m.getDtAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getDtAPI.restype = ctypes.c_double
        self.simDLL_m.incrementSimulationTimeByDtAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.incrementSimulationTimeByDtAPI.restype = None
        self.simDLL_m.getNumberOfParticleTypesAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfParticleTypesAPI.restype = ctypes.c_int
        self.simDLL_m.getNumberOfParticlesPerTypeAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfParticlesPerTypeAPI.restype = ctypes.c_int
        self.simDLL_m.getNumberOfAttributesTrackedAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfAttributesTrackedAPI.restype = ctypes.c_int
        self.simDLL_m.areResultsPreparedAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.areResultsPreparedAPI.restype = ctypes.c_bool
        self.simDLL_m.resetParticlesEscapedCountAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.resetParticlesEscapedCountAPI.restype = None
        self.simDLL_m.getNormalizedAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNormalizedAPI.restype = ctypes.c_bool
        self.simDLL_m.getSimMinAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimMinAPI.restype = ctypes.c_double
        self.simDLL_m.getSimMaxAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimMaxAPI.restype = ctypes.c_double

        #Pointer one liners
        self.simDLL_m.getPointerToParticlesInSimArrayAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getPointerToParticlesInSimArrayAPI.restype = ctypes.POINTER(ctypes.c_bool)
        self.simDLL_m.getPointerToSingleParticleAttributeArrayAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.getPointerToSingleParticleAttributeArrayAPI.restype = ctypes.POINTER(ctypes.c_double)

        #Numerical tools
        self.simDLL_m.generateNormallyDistributedValuesAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
        self.simDLL_m.generateNormallyDistributedValuesAPI.restype = None
        self.simDLL_m.calculateMeanOfParticleAttributeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_bool)
        self.simDLL_m.calculateMeanOfParticleAttributeAPI.restype = ctypes.c_double
        self.simDLL_m.calculateStdDevOfParticleAttributeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.calculateStdDevOfParticleAttributeAPI.restype = ctypes.c_double

        #Field tools
        self.simDLL_m.calculateBFieldAtZandTimeAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateBFieldAtZandTimeAPI.restype = ctypes.c_double
        self.simDLL_m.calculateEFieldAtZandTimeAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateEFieldAtZandTimeAPI.restype = ctypes.c_double

        #Simulation management
        self.simDLL_m.createSimulation170925.argtypes = (ctypes.c_char_p,)
        self.simDLL_m.createSimulation170925.restype = ctypes.c_void_p
        self.simDLL_m.initializeSimulationAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.initializeSimulationAPI.restype = None
        self.simDLL_m.copyDataToGPUAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToGPUAPI.restype = None
        self.simDLL_m.iterateSimulationAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.iterateSimulationAPI.restype = None
        self.simDLL_m.copyDataToHostAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToHostAPI.restype = None
        self.simDLL_m.freeGPUMemoryAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.freeGPUMemoryAPI.restype = None
        self.simDLL_m.prepareResultsAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.prepareResultsAPI.restype = ctypes.POINTER(ctypes.c_double)

        #Satellite functions
        self.simDLL_m.getNumberOfSatellitesAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfSatellitesAPI.restype = ctypes.c_int
        self.simDLL_m.getNumberOfSatelliteMsmtsAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfSatelliteMsmtsAPI.restype = ctypes.c_int
        self.simDLL_m.getSatelliteDataPointersAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.getSatelliteDataPointersAPI.restype = ctypes.POINTER(ctypes.c_double)


        #Now code for init
        crootdir = ctypes.create_string_buffer(bytes(self.rootdir_m, encoding='utf-8'))
        self.simulationptr = ctypes.c_void_p
        self.simulationptr = self.simDLL_m.createSimulation170925(crootdir)

        self.types_m = self.simDLL_m.getNumberOfParticleTypesAPI(self.simulationptr)
        self.attr_m = self.simDLL_m.getNumberOfAttributesTrackedAPI(self.simulationptr)
        self.numPart_m = self.simDLL_m.getNumberOfParticlesPerTypeAPI(self.simulationptr)
        self.dt_m = self.simDLL_m.getDtAPI(self.simulationptr)
        self.normalized_m = self.simDLL_m.getNormalizedAPI(self.simulationptr)
        self.simMin_m = self.simDLL_m.getSimMinAPI(self.simulationptr)
        self.simMax_m = self.simDLL_m.getSimMaxAPI(self.simulationptr)

        self.printSimCharacteristics()

        return
    
    

    #Run Simulation
    def runSim(self, iterations, inSimOnly=True):
        self.initializeSimulation()
        self.copyDataToGPU()
        self.iterateSimulation(iterations)
        self.copyDataToHost()
        self.freeGPUMemory()
        self.prepareResults()
        return self.getResultsfrom3D(inSimOnly)

    ###Member functions for Simulation class
    #One liner functions
    def getTime(self):
        return self.simDLL_m.getSimulationTimeAPI(self.simulationptr)

    def incTime(self):
        self.simDLL_m.incrementSimulationTimeByDtAPI(self.simulationptr)

    def resetEscapedCount(self):
        self.simDLL_m.resetParticlesEscapedCountAPI(self.simulationptr)
    
    #Pointer one liners (one liners in CPP obv, not Python)
    def getPartInSimBoolArray(self):
        ret = []
        partbool = []
        for iii in range(self.types_m):
            partbool_c = self.simDLL_m.getPointerToParticlesInSimArrayAPI(self.simulationptr, iii)
            for jjj in range(self.numPart_m):
                partbool.append(partbool_c[jjj])
            ret.append(partbool)
            partbool = []
        return ret

    def getResultsfrom3D(self, inSimOnly=True):
        if (inSimOnly):
            inSimBoolArray = self.getPartInSimBoolArray()
        
        ret = []
        partattr = []
        partdbl = []
        for iii in range(self.types_m):
            for jjj in range(self.attr_m):
                partdbl_c = self.simDLL_m.getPointerToSingleParticleAttributeArrayAPI(self.simulationptr, iii, jjj)
                for kk in range(self.numPart_m):
                   if(inSimOnly):
                       if(inSimBoolArray[iii][kk]):
                           partdbl.append(partdbl_c[kk])
                   else:
                       partdbl.append(partdbl_c[kk])
                partattr.append(partdbl)
                partdbl = []
            ret.append(partattr)
            partattr = []
        return ret

    #Numerical tools
    def __generateNormallyDistributedCPP(self, numberOfNormalAttributes, means, sigmas):
        #generally, a constructor in the derived c++ class will take care of this - that's why __
        C_DOUBLEA = ctypes.c_double * numberOfNormalAttributes
        meanscpp = C_DOUBLEA()
        sigmascpp = C_DOUBLEA()

        for iii in range(numberOfNormalAttributes):
            meanscpp[iii] = means[iii]
            sigmascpp[iii] = sigmas[iii]

        self.simDLL_m.generateNormallyDistributedValuesAPI(self.simulationptr, numberOfNormalAttributes, meanscpp, sigmascpp)

    def meanOfPartAttr(self, particleIndex, attributeIndex, absValueBool=False):
        return self.simDLL_m.calculateMeanOfParticleAttributeAPI(self.simulationptr, particleIndex, attributeIndex, absValueBool)

    def stddevOfPartAttr(self, particleIndex, attributeIndex):
        return self.simDLL_m.calculateStdDevOfParticleAttributeAPI(self.simulationptr, particleIndex, attributeIndex)

    #Field tools
    def BFieldatZandT(self, z, time):
        return self.simDLL_m.calculateBFieldAtZandTimeAPI(self.simulationptr, z, time)

    def EFieldatZandT(self, z, time):
        return self.simDLL_m.calculateEFieldAtZandTimeAPI(self.simulationptr, z, time)

    def fieldsAtAllZ(self, time, bins, binsize, z0):
        B_z = []
        E_z = []
        B_E_z_dim = []
        for iii in range(bins):
            B_z.append(self.BFieldatZandT(z0 + binsize * iii, time))
            E_z.append(self.EFieldatZandT(z0 + binsize * iii, time))
            B_E_z_dim.append(z0 + binsize * iii)
        return [B_z, E_z, B_E_z_dim]

    #Simulation management
    #Functions are prepended with __ because the intent is to simply runSim which will call them all
    #however if more refined control is needed, call them one by one and ignore runSim
    def initializeSimulation(self):
        self.simDLL_m.initializeSimulationAPI(self.simulationptr)

    def copyDataToGPU(self):
        self.simDLL_m.copyDataToGPUAPI(self.simulationptr)

    def iterateSimulation(self, numberOfIterations):
        print("Number Of Iterations: ", numberOfIterations)
        self.simDLL_m.iterateSimulationAPI(self.simulationptr, numberOfIterations)

    def copyDataToHost(self):
        self.simDLL_m.copyDataToHostAPI(self.simulationptr)

    def freeGPUMemory(self):
        self.simDLL_m.freeGPUMemoryAPI(self.simulationptr)

    def prepareResults(self):
        return self.simDLL_m.prepareResultsAPI(self.simulationptr)

    #Satellite functions
    def getNumberOfSatellites(self):
        return self.simDLL_m.getNumberOfSatellitesAPI(self.simulationptr)

    def getNumberOfSatelliteMsmts(self):
        return self.simDLL_m.getNumberOfSatelliteMsmtsAPI(self.simulationptr)

    def getSatelliteData(self):
        self.satMsmts_m = self.getNumberOfSatelliteMsmts()
        self.satNum_m = self.getNumberOfSatellites()
        
        data = []

        msmtptr = [] #constructs array of double pointers so z value can be checked before recording data
        for iii in range(self.satMsmts_m):
            satptr = []
            for jjj in range(self.satNum_m):
                attrptr = []
                for kk in range(self.attr_m):
                    attrptr.append(self.simDLL_m.getSatelliteDataPointersAPI(self.simulationptr, iii, jjj, kk))
                satptr.append(attrptr)
            msmtptr.append(satptr)

        for iii in range(self.satMsmts_m):
            sat = []
            for jjj in range(self.satNum_m):
                attr = []
                for kk in range(self.attr_m):
                    parts=[]
                    for lll in range(self.numPart_m):
                        #if (msmtptr[iii][jjj][2][lll] > self.simMin_m):
                        parts.append(msmtptr[iii][jjj][kk][lll])
                    attr.append(parts)
                sat.append(attr)
            data.append(sat)
        
        return data

    def printSimCharacteristics(self):
        print("Sim between:          ", self.simMin_m, " - ", self.simMax_m, " Re" if self.normalized_m else " m")
        print("Number of Particles:  ", self.numPart_m)
        print("Replenish lost part:  ", "True" if self.replenish_m else "False")
        print("Sim dt:               ", self.dt_m, " s")


if __name__ == '__main__':
    print("SimulationAPI.py is not meant to be called as main.  Run a simulation file and that will import this.")
