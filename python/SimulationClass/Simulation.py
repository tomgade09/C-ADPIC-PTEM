import ctypes, os, sys, inspect

sys.path.append(os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)) + './../'))
from __simulationvariables import *

class Simulation:
    def __init__(self, DLLloc, rootdir, dt, simMin, simMax, ionT, magT):
        self.dllLoc_m = DLLloc
        self.rootdir_m = rootdir
        self.dt_m = dt
        self.simMin_m = simMin
        self.simMax_m = simMax
        self.ionT_m = ionT
        self.magT_m = magT
        self.numTypes_m = 0
        self.numParts_m = []
        self.nameParts_m = []
        self.numAttrs_m = []
        self.satPartInd_m = []
        self.normalized_m = False

        self.simDLL_m = ctypes.CDLL(self.dllLoc_m)

        #One liner functions
        self.simDLL_m.getSimulationTimeAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimulationTimeAPI.restype = ctypes.c_double
        self.simDLL_m.getDtAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getDtAPI.restype = ctypes.c_double
        self.simDLL_m.getSimMinAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimMinAPI.restype = ctypes.c_double
        self.simDLL_m.getSimMaxAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getSimMaxAPI.restype = ctypes.c_double
        self.simDLL_m.incrementSimulationTimeByDtAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.incrementSimulationTimeByDtAPI.restype = None
        self.simDLL_m.getNumberOfParticleTypesAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfParticleTypesAPI.restype = ctypes.c_int
        self.simDLL_m.getNumberOfParticlesAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getNumberOfParticlesAPI.restype = ctypes.c_int
        self.simDLL_m.getNumberOfAttributesAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getNumberOfAttributesAPI.restype = ctypes.c_int
        self.simDLL_m.areResultsPreparedAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.areResultsPreparedAPI.restype = ctypes.c_bool
        
        #Pointer one liners
        self.simDLL_m.getPointerToParticleAttributeArrayAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_bool)
        self.simDLL_m.getPointerToParticleAttributeArrayAPI.restype = ctypes.POINTER(ctypes.c_double)

        #Mu<->VPerp functions
        self.simDLL_m.convertParticleVPerpToMuAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.convertParticleVPerpToMuAPI.restype = None
        self.simDLL_m.convertParticleMuToVPerpAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.convertParticleMuToVPerpAPI.restype = None

        #Field tools
        self.simDLL_m.calculateBFieldAtZandTimeAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateBFieldAtZandTimeAPI.restype = ctypes.c_double
        self.simDLL_m.calculateEFieldAtZandTimeAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.calculateEFieldAtZandTimeAPI.restype = ctypes.c_double

        #Simulation management
        self.simDLL_m.createSimulationAPI.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p)
        self.simDLL_m.createSimulationAPI.restype = ctypes.c_void_p
        self.simDLL_m.initializeSimulationAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.initializeSimulationAPI.restype = None
        self.simDLL_m.copyDataToGPUAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToGPUAPI.restype = None
        self.simDLL_m.iterateSimulationAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.iterateSimulationAPI.restype = None
        self.simDLL_m.copyDataToHostAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.copyDataToHostAPI.restype = None
        self.simDLL_m.freeGPUMemoryAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.freeGPUMemoryAPI.restype = None
        self.simDLL_m.prepareResultsAPI.argtypes = (ctypes.c_void_p, ctypes.c_bool)
        self.simDLL_m.prepareResultsAPI.restype = ctypes.POINTER(ctypes.c_double)
        self.simDLL_m.terminateSimulationAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.terminateSimulationAPI.restype = None
        self.simDLL_m.loadCompletedSimDataAPI.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int)
        self.simDLL_m.loadCompletedSimDataAPI.restype = None
        self.simDLL_m.runNormalSimulationAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_char_p)
        self.simDLL_m.runNormalSimulationAPI.restype = None
        self.simDLL_m.setBFieldModelAPI.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double)
        self.simDLL_m.setBFieldModelAPI.restype = None

        #Satellite functions
        self.simDLL_m.createSatelliteAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_bool, ctypes.c_char_p)
        self.simDLL_m.createSatelliteAPI.restype = None
        self.simDLL_m.getNumberOfSatellitesAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getNumberOfSatellitesAPI.restype = ctypes.c_int
        self.simDLL_m.getSatNumOfDetectedPartsAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getSatNumOfDetectedPartsAPI.restype = ctypes.c_int        
        self.simDLL_m.getSatelliteDataPointersAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.getSatelliteDataPointersAPI.restype = ctypes.POINTER(ctypes.c_double)
        self.simDLL_m.writeSatelliteDataToCSVAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.writeSatelliteDataToCSVAPI.restype = None

        #Particle Management
        self.simDLL_m.createParticleTypeAPI.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_char_p)
        self.simDLL_m.createParticleTypeAPI.restype = None

        #Log File
        self.simDLL_m.getLogFilePointerAPI.argtypes = (ctypes.c_void_p,)
        self.simDLL_m.getLogFilePointerAPI.restype = ctypes.c_void_p
        self.simDLL_m.writeLogFileEntryAPI.argtypes = (ctypes.c_void_p, ctypes.c_char_p)
        self.simDLL_m.writeLogFileEntryAPI.restype = None
        self.simDLL_m.writeTimeDiffFromNowAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p)
        self.simDLL_m.writeTimeDiffFromNowAPI.restype = None
        self.simDLL_m.writeTimeDiffAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        self.simDLL_m.writeTimeDiffAPI.restype = None

        #Now code for init
        rootdirBuf = ctypes.create_string_buffer(bytes(self.rootdir_m, encoding='utf-8'))
        self.simulationptr = ctypes.c_void_p
        self.simulationptr = self.simDLL_m.createSimulationAPI(dt, simMin, simMax, ionT, magT, rootdirBuf)
        self.logFileObj_m = self.simDLL_m.getLogFilePointerAPI(self.simulationptr)

        self.particles_m = []
        self.satellites_m = []

        return
    
    

    #Run Simulation
    def runSim(self, iterations, loadFile=False):
        if (loadFile):
            loadFileBuf = ctypes.create_string_buffer(bytes(DISTINFOLDER, encoding='utf-8'))
        else:
            loadFileBuf = ctypes.create_string_buffer(bytes("", encoding='utf-8'))

        #eventually replace with the functions to get the data from the CPP class, name these arrays better
        self.numTypes_m = 2
        self.numAttrs_m = [3, 3]
        self.numParts_m = [115200, 115200]
        self.nameParts_m = ["elec", "ions"]
        self.satellites_m = ["btmElec", "btmIons", "topElec", "topIons"]
        self.satPartInd_m = [0, 1, 0, 1]

        self.simDLL_m.runNormalSimulationAPI(self.simulationptr, iterations, 500, loadFileBuf)
        self.simDLL_m.writeSatelliteDataToCSVAPI(self.simulationptr)

        return self.getFinalDataAllParticles(), self.getOriginalDataAllParticles(), self.getSatelliteData()  #Returns final particle data, original particle data, satellite data

    ###Member functions for Simulation class
    #Particle, Satellite Mgmt
    #void loadCompletedSimDataAPI(Simulation* simulation, const char* fileDir, const char* partNames, const char* attrNames, const char* satNames, int numParts)
    def loadCompletedSimData(self, fileDir, partNames, attrNames, satNames, numParts, numPartTypes, numAttrs, numSats):
        fileDirBuf = ctypes.create_string_buffer(bytes(fileDir, encoding='utf-8'))
        partNamesBuf = ctypes.create_string_buffer(bytes(partNames, encoding='utf-8'))
        attrNamesBuf = ctypes.create_string_buffer(bytes(attrNames, encoding='utf-8'))
        satNamesBuf = ctypes.create_string_buffer(bytes(satNames, encoding='utf-8'))
        self.simDLL_m.loadCompletedSimDataAPI(self.simulationptr, fileDirBuf, partNamesBuf, attrNamesBuf, satNamesBuf, numParts)
        
        self.numTypes_m = numPartTypes
        
        for iii in range(numAttrs + 2):
            self.numAttrs_m.append(numAttrs)
        for iii in range(numPartTypes):
            self.numParts_m.append(numParts)
            self.nameParts_m.append("1")
        for iii in range(numSats):
            self.satPartInd_m.append(0)

        return self.getFinalDataAllParticles(), self.getOriginalDataAllParticles(), self.getSatelliteData()

    def createParticle(self, name, attrNames, mass, charge, numParts, posDims, velDims, normFactor, loadFileDir=""):
        nameBuf = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        attrNamesBuf = ctypes.create_string_buffer(bytes(attrNames, encoding='utf-8'))
        loadFileDirBuf = ctypes.create_string_buffer(bytes(loadFileDir, encoding='utf-8'))

        self.simDLL_m.createParticleTypeAPI(self.simulationptr, nameBuf, attrNamesBuf, mass, charge, numParts, posDims, velDims, normFactor, loadFileDirBuf)
        self.particles_m.append(name)

        self.numTypes_m += 1 #eventually check to see that C++ has created properly by calling Particle access functions or don't create a python one
        self.numAttrs_m.append(posDims + velDims)
        self.numParts_m.append(numParts)
        self.nameParts_m.append(name)

    def createSatellite(self, particleInd, altitude, upwardFacing, name):
        nameBuf = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        self.simDLL_m.createSatelliteAPI(self.simulationptr, particleInd, altitude, upwardFacing, nameBuf)
        self.satellites_m.append(name)
        self.satPartInd_m.append(particleInd)
    
    def convertParticleVPerpToMu(particleInd):
        self.simDLL_m.convertParticleVPerpToMuAPI(self.simulationptr, particleInd)

    def convertParticleMuToVPerp(particleInd):
        self.simDLL_m.convertParticleMuToVPerpAPI(self.simulationptr, particleInd)


    #One liner functions
    def getTime(self):
        return self.simDLL_m.getSimulationTimeAPI(self.simulationptr)

    def incTime(self):
        self.simDLL_m.incrementSimulationTimeByDtAPI(self.simulationptr)
    
    #Pointer one liners (one liners in CPP obv, not Python)
    def getParticleDataFromCPP(self, origData=False):
        ret = []
        partattr = []
        partdbl = []
        for iii in range(self.numTypes_m):
            for jjj in range(self.numAttrs_m[iii]):
                partdbl_c = self.simDLL_m.getPointerToParticleAttributeArrayAPI(self.simulationptr, iii, jjj, origData)
                for kk in range(self.numParts_m[iii]):
                   partdbl.append(partdbl_c[kk])
                partattr.append(partdbl)
                partdbl = []
            ret.append(partattr)
            partattr = []
        return ret

    def getOriginalDataAllParticles(self):
        return self.getParticleDataFromCPP(True)

    def getFinalDataAllParticles(self):
        return self.getParticleDataFromCPP()


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

    def iterateSimulation(self, numberOfIterations, itersBtwCouts):
        print("Number Of Iterations: ", numberOfIterations)
        self.simDLL_m.iterateSimulationAPI(self.simulationptr, numberOfIterations, itersBtwCouts)

    def copyDataToHost(self):
        self.simDLL_m.copyDataToHostAPI(self.simulationptr)

    def freeGPUMemory(self):
        self.simDLL_m.freeGPUMemoryAPI(self.simulationptr)

    def prepareResults(self, normalizeToRe=True):
        return self.simDLL_m.prepareResultsAPI(self.simulationptr, normalizeToRe)

    def terminateSimulation(self):
        self.simDLL_m.terminateSimulationAPI(self.simulationptr)

    #Satellite functions
    def getNumberOfSatellites(self):
        return self.simDLL_m.getNumberOfSatellitesAPI(self.simulationptr)

    def getSatNumberOfDetectedParticles(self, satInd):
        return self.simDLL_m.getSatNumOfDetectedPartsAPI(self.simulationptr, satInd)

    def getSatelliteData(self): #Need to add case where satellite captures more or less than the number of particles
        self.satNum_m = self.getNumberOfSatellites()

        satptr = [] #constructs array of double pointers so z value can be checked before recording data
        for jjj in range(self.satNum_m):
            attrptr = []
            for kk in range(self.numAttrs_m[self.satPartInd_m[jjj]] + 2):
                attrptr.append(self.simDLL_m.getSatelliteDataPointersAPI(self.simulationptr, jjj, kk))
            satptr.append(attrptr)
        
        satsdata = []
        for jjj in range(self.satNum_m):
            attr = []
            for kk in range(self.numAttrs_m[self.satPartInd_m[jjj]] + 2):
                parts=[]
                for lll in range(self.numParts_m[self.satPartInd_m[jjj]]):
                    parts.append(satptr[jjj][kk][lll])
                attr.append(parts)
            satsdata.append(attr)
        
        return satsdata

    def logWriteEntry(self, logMessage):
        logMessCbuf = ctypes.create_string_buffer(bytes(logMessage, encoding='utf-8'))
        self.simDLL_m.writeLogFileEntryAPI(self.logFileObj_m, logMessCbuf)

    def logWriteTimeDiffFromNow(self, startTSind, nowLabel):
        nowLabCbuf = ctypes.create_string_buffer(bytes(nowLabel, encoding='utf-8'))
        self.simDLL_m.writeTimeDiffFromNowAPI(self.logFileObj_m, startTSind, nowLabCbuf)

    def logWriteTimeDiff(self, startTSind, endTSind):
        self.simDLL_m.writeTimeDiffAPI(self.logFileObj_m, startTSind, endTSind)


if __name__ == '__main__':
    print("SimulationAPI.py is not meant to be called as main.  Run a simulation file and that will import this.")
