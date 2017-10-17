from __future__ import absolute_import, division, print_function
import ctypes

class Simulation:
    def __init__(self, rootdir, DLLloc):
        self.dllLoc_m = DLLloc
        self.rootdir_m = rootdir
        self.simDLL_m = ctypes.CDLL(self.dllLoc_m)

        #One liner functions
        self.simDLL_m.getSimulationTimeWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimulationTimeWrapper.restype = ctypes.c_double
        self.simDLL_m.getDtWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getDtWrapper.restype = ctypes.c_double
        self.simDLL_m.incrementSimulationTimeByDtWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.incrementSimulationTimeByDtWrapper.restype = None
        self.simDLL_m.getNumberOfParticleTypesWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfParticleTypesWrapper.restype = ctypes.c_int
        self.simDLL_m.getNumberOfParticlesPerTypeWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfParticlesPerTypeWrapper.restype = ctypes.c_int
        self.simDLL_m.getNumberOfAttributesTrackedWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfAttributesTrackedWrapper.restype = ctypes.c_int
        self.simDLL_m.areResultsPreparedWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.areResultsPreparedWrapper.restype = ctypes.c_bool
        self.simDLL_m.resetParticlesEscapedCountWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.resetParticlesEscapedCountWrapper.restype = None

        #Pointer one liners
        self.simDLL_m.getPointerToSerializedParticleArrayWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getPointerToSerializedParticleArrayWrapper.restype = ctypes.POINTER(ctypes.c_double)
        self.simDLL_m.getPointerToParticlesInSimArrayWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getPointerToParticlesInSimArrayWrapper.restype = ctypes.POINTER(ctypes.c_bool)
        self.simDLL_m.getPointerToSingleParticleAttributeArrayWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.getPointerToSingleParticleAttributeArrayWrapper.restype = ctypes.POINTER(ctypes.c_double)

        #Numerical tools
        self.simDLL_m.generateNormallyDistributedValuesWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
        self.simDLL_m.generateNormallyDistributedValuesWrapper.restype = None
        self.simDLL_m.calculateMeanOfParticleAttributeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_bool)
        self.simDLL_m.calculateMeanOfParticleAttributeWrapper.restype = ctypes.c_double
        self.simDLL_m.calculateStdDevOfParticleAttributeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.calculateStdDevOfParticleAttributeWrapper.restype = ctypes.c_double

        #Array tools
        self.simDLL_m.serializeParticleArrayWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.serializeParticleArrayWrapper.restype = None
        self.simDLL_m.calculateBFieldAtZandTimeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateBFieldAtZandTimeWrapper.restype = ctypes.c_double
        self.simDLL_m.calculateEFieldAtZandTimeWrapper.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateEFieldAtZandTimeWrapper.restype = ctypes.c_double

        #Simulation management
        self.simDLL_m.createSimulation170925.argtypes = (ctypes.c_char_p,)
        self.simDLL_m.createSimulation170925.restype = ctypes.c_void_p
        self.simDLL_m.initializeSimulationWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.initializeSimulationWrapper.restype = None
        self.simDLL_m.copyDataToGPUWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToGPUWrapper.restype = None
        self.simDLL_m.iterateSimulationWrapper.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.iterateSimulationWrapper.restype = None
        self.simDLL_m.copyDataToHostWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToHostWrapper.restype = None
        self.simDLL_m.freeGPUMemoryWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.freeGPUMemoryWrapper.restype = None
        self.simDLL_m.prepareResultsWrapper.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.prepareResultsWrapper.restype = ctypes.POINTER(ctypes.c_double)

        #Now code for init
        crootdir = ctypes.create_string_buffer(bytes(self.rootdir_m, encoding='utf-8'))
        self.simulationptr = ctypes.c_void_p
        self.simulationptr = self.simDLL_m.createSimulation170925(crootdir)

        self.types_m = self.simDLL_m.getNumberOfParticleTypesWrapper(self.simulationptr)
        self.attr_m = self.simDLL_m.getNumberOfAttributesTrackedWrapper(self.simulationptr)
        self.numPart_m = self.simDLL_m.getNumberOfParticlesPerTypeWrapper(self.simulationptr)
        self.dt_m = self.simDLL_m.getDtWrapper(self.simulationptr)

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
        return self.simDLL_m.getSimulationTimeWrapper(self.simulationptr)

    def incTime(self):
        self.simDLL_m.incrementSimulationTimeByDtWrapper(self.simulationptr)

    def resetEscapedCount(self):
        self.simDLL_m.resetParticlesEscapedCountWrapper(self.simulationptr)
    
    #Pointer one liners (one liners in CPP obv, not Python)
    def getResultsfromSerialized(self):
        serArray = self.simDLL_m.getPointerToSerializedParticleArrayWrapper(self.simulationptr)
        if (not(serArray)):
            return [[0,0,0],[0,0,0]]#needs to be tested
        
        attributesTracked = int(serArray[1])
        if not(attributesTracked == self.attr_m):
            print("Error in getResultsfromSerialized (python): Number of attributes reported in serialized array does not match what was reported from the C++ class.  Something is amiss.")
            print("Attributes in serialized array: ", attributesTracked, "  Attributes from C++ class: ", self.attr_m)
            print("Returning.")
            return [[0,0,0],[0,0,0]]

        electrons = int(serArray[2])
        ions = int(serArray[electrons * attributesTracked + 3])
        length = (electrons + ions) * attributesTracked + 4

        self.elecSerial = electrons
        self.ionsSerial = ions

        v_e_para = []#code is specific to electrons/ions and 3 attributes - could make more general with loops, arrays.append(size), etc
        v_e_perp = []
        z_e = []
        v_i_para = []
        v_i_perp = []
        z_i = []

        for iii in range(electrons):
            v_e_para.append(serArray[iii + 3]) 
            v_e_perp.append(serArray[iii + electrons + 3])
            z_e.append(serArray[iii + 2 * electrons + 3])
        
        for iii in range(ions):
            v_i_para.append(serArray[iii + 3 * electrons + 4])
            v_i_perp.append(serArray[iii + 3 * electrons + ions + 4])
            z_i.append(serArray[iii + 3 * electrons + 2 * ions + 4])

        return [[v_e_para, v_e_perp, z_e], [v_i_para, v_i_perp, z_i]]

    def getPartInSimBoolArray(self):
        ret = []
        partbool = []
        for iii in range(self.types_m):
            partbool_c = self.simDLL_m.getPointerToParticlesInSimArrayWrapper(self.simulationptr, iii)
            for jjj in range(self.numPart_m):
                partbool.append(partbool_c[jjj])
            ret.append(partbool)
            partbool = []
        return ret

    def getResultsfrom3D(self, inSimOnly=True):
        if (inSimOnly):
            inSimBoolArray = self.getPartInSimBoolArray()
        
        ret = [] #Generic form of serialized array above - gets pointers to 1D component arrays of the 3D c++ array and reconstructs the array in python
        partattr = []
        partdbl = []
        for iii in range(self.types_m):
            for jjj in range(self.attr_m):
                partdbl_c = self.simDLL_m.getPointerToSingleParticleAttributeArrayWrapper(self.simulationptr, iii, jjj)
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
        C_DOUBLEA = ctypes.c_double * numberOfNormalAttributes#generally, a constructor in the derived c++ class will take care of this - that's why __
        meanscpp = C_DOUBLEA()
        sigmascpp = C_DOUBLEA()

        for iii in range(numberOfNormalAttributes):
            meanscpp[iii] = means[iii]
            sigmascpp[iii] = sigmas[iii]

        self.simDLL_m.generateNormallyDistributedValuesWrapper(self.simulationptr, numberOfNormalAttributes, meanscpp, sigmascpp)

    def meanOfPartAttr(self, particleIndex, attributeIndex, absValueBool=False):
        return self.simDLL_m.calculateMeanOfParticleAttributeWrapper(self.simulationptr, particleIndex, attributeIndex, absValueBool)

    def stddevOfPartAttr(self, particleIndex, attributeIndex):
        return self.simDLL_m.calculateStdDevOfParticleAttributeWrapper(self.simulationptr, particleIndex, attributeIndex)

    #Array tools
    def serializeArray(self):
        self.simDLL_m.serializeParticleArrayWrapper(self.simulationptr)

    def BFieldatZandT(self, z, time):
        return self.simDLL_m.calculateBFieldAtZandTimeWrapper(self.simulationptr, z, time)

    def EFieldatZandT(self, z, time):
        return self.simDLL_m.calculateEFieldAtZandTimeWrapper(self.simulationptr, z, time)

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
        self.simDLL_m.initializeSimulationWrapper(self.simulationptr)

    def copyDataToGPU(self):
        self.simDLL_m.copyDataToGPUWrapper(self.simulationptr)

    def iterateSimulation(self, numberOfIterations):
        self.simDLL_m.iterateSimulationWrapper(self.simulationptr, numberOfIterations)

    def copyDataToHost(self):
        self.simDLL_m.copyDataToHostWrapper(self.simulationptr)

    def freeGPUMemory(self):
        self.simDLL_m.freeGPUMemoryWrapper(self.simulationptr)

    def prepareResults(self):
        return self.simDLL_m.prepareResultsWrapper(self.simulationptr)

    ###Tests
    def compareSerialWith3D(self):
        results3D = self.getResultsfrom3D()
        resultsSer = self.getResultsfromSerialized()
        
        e_wrong_para = 0
        e_wrong_perp = 0
        e_wrong_z = 0

        for iii in range(self.elecSerial):
            if (results3D[0][0][iii] != resultsSer[0][0][iii]):
                e_wrong_para += 1
            if (results3D[0][1][iii] != resultsSer[0][1][iii]):
                e_wrong_perp += 1
            if (results3D[0][2][iii] != resultsSer[0][2][iii]):
                e_wrong_z += 1

        i_wrong_para = 0
        i_wrong_perp = 0
        i_wrong_z = 0

        for iii in range(self.ionsSerial):
            if (results3D[1][0][iii] != resultsSer[1][0][iii]):
                i_wrong_para += 1
            if (results3D[1][1][iii] != resultsSer[1][1][iii]):
                i_wrong_perp += 1
            if (results3D[1][2][iii] != resultsSer[1][2][iii]):
                i_wrong_z += 1

        if (((e_wrong_para and e_wrong_perp and e_wrong_z) and (i_wrong_para and i_wrong_perp and i_wrong_z)) == 0):
            print("Arrays are identical.")
        if (not(e_wrong_para == 0) or not(e_wrong_perp == 0) or not(e_wrong_z == 0)):
            print("electrons wrong: ", e_wrong_para, e_wrong_perp, e_wrong_z)
        if (not(i_wrong_para == 0) or not(i_wrong_perp == 0) or not(i_wrong_z == 0)):
            print("ions wrong: ", i_wrong_para, i_wrong_perp, i_wrong_z)


if __name__ == '__main__':
    print("SimulationAPI.py is not meant to be called as main.  Run a simulation file and that will import this.")
