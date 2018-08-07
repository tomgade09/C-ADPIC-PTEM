import ctypes, os, sys, inspect

sys.path.append(os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)) + './../'))
from __simulationvariables import *
import _Simulation

class Simulation(_Simulation._SimulationCDLL):
    def __init__(self, DLLloc, savedir, dt=None, simMin=None, simMax=None):
        super().__init__(DLLloc)
        self.savedir_m = savedir
        self.dt_m = dt
        self.simMin_m = simMin
        self.simMax_m = simMax

        self.numPartTypes_m = 0
        self.numSats_m = 0
        self.numParts_m = []
        self.nameParts_m = []
        self.numAttrs_m = []
        self.satPartInd_m = []
        self.nameSat_m = []

        self.normalSim = False

        #Now code for init
        savedirBuf = ctypes.create_string_buffer(bytes(self.savedir_m, encoding='utf-8'))
        self.simulationptr = ctypes.c_void_p
        
        if dt is not None:
            self.simulationptr = self.simDLL_m.createSimulationAPI(dt, simMin, simMax, savedirBuf)
            self.logFileObj_m = self.simDLL_m.getLogFilePointerAPI(self.simulationptr)

        return


    #Run Simulation
    def setupExampleSim(self, numParts):
        if (LOADDIST):
            loadFileBuf = ctypes.create_string_buffer(bytes(DISTINFOLDER, encoding='utf-8'))
        else:
            loadFileBuf = ctypes.create_string_buffer(bytes("", encoding='utf-8'))

        self.simDLL_m.setupExampleSimulationAPI(self.simulationptr, numParts, loadFileBuf)
        if self.numAttrs_m == []:
            self.getSimChars()

        self.normalSim = True

    def runExampleSim(self, iterations, iterBtwCouts):
        if self.numAttrs_m == []:
            self.getSimChars()

        self.simDLL_m.runExampleSimulationAPI(self.simulationptr, iterations, iterBtwCouts)
        
        return self.getFinalDataAllParticles(), self.getOriginalDataAllParticles(), self.getSatelliteData()  #Returns final particle data, original particle data, satellite data

    def writeCommonCSV(self):
        self.simDLL_m.writeCommonCSVAPI(self.simulationptr)


    #Fields management
    def setBFieldModel(name, doublesString):
        name_c = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        doublesString_c = ctypes.create_string_buffer(bytes(doublesString, encoding='utf-8'))
        self.simDLL_m.setBFieldModelAPI(self.simulationptr, name_c, doublesString_c)

    def addEFieldModel(name, doublesString):
        name_c = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        doublesString_c = ctypes.create_string_buffer(bytes(doublesString, encoding='utf-8'))
        self.simDLL_m.addEFieldModelAPI(self.simulationptr, name_c, doublesString_c)

    ###Member functions for Simulation class
    #Particle, Satellite Mgmt
    def getSimChars(self): #call this after the simulation is set up
        if not self.numAttrs_m == []:
            return

        self.numPartTypes_m = self.simDLL_m.getNumberOfParticleTypesAPI(self.simulationptr)
        for part in range(self.numPartTypes_m):
            self.numAttrs_m.append(self.simDLL_m.getNumberOfAttributesAPI(self.simulationptr, part))
            self.numParts_m.append(self.simDLL_m.getNumberOfParticlesAPI(self.simulationptr, part))
            self.nameParts_m.append(self.simDLL_m.getParticleNameAPI(self.simulationptr, part))
        self.numSats_m = self.simDLL_m.getNumberOfSatellitesAPI(self.simulationptr)
        for sat in range(self.numSats_m):
            self.nameSat_m.append(self.simDLL_m.getSatelliteNameAPI(self.simulationptr, sat))
            if "lec" in self.nameSat_m[sat]:
                self.satPartInd_m.append(0)
            else:
                self.satPartInd_m.append(1)

    def loadCompletedSimData(self, fileDir):
        if self.dt_m is not None:
            print("Error, dt specified meaning that an instance of the Simulation class has already been created in c++ through Python. Returning")
            return

        fileDirBuf = ctypes.create_string_buffer(bytes(fileDir, encoding='utf-8'))
        
        self.simulationptr = self.simDLL_m.loadCompletedSimDataAPI(fileDirBuf)

        if self.numAttrs_m == []:
            self.getSimChars()

        self.dt_m = self.simDLL_m.getDtAPI(self.simulationptr)
        self.simMin_m = self.simDLL_m.getSimMinAPI(self.simulationptr)
        self.simMax_m = self.simDLL_m.getSimMaxAPI(self.simulationptr)

        return self.getFinalDataAllParticles(), self.getOriginalDataAllParticles(), self.getSatelliteData()

    def createParticle(self, name, attrNames, mass, charge, numParts, posDims, velDims, normFactor, loadFileDir=""):
        nameBuf = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        attrNamesBuf = ctypes.create_string_buffer(bytes(attrNames, encoding='utf-8'))
        loadFileDirBuf = ctypes.create_string_buffer(bytes(loadFileDir, encoding='utf-8'))

        self.simDLL_m.createParticleTypeAPI(self.simulationptr, nameBuf, attrNamesBuf, mass, charge, numParts, posDims, velDims, normFactor, loadFileDirBuf)

        self.numPartTypes_m += 1 #eventually check to see that C++ has created properly by calling Particle access functions or don't create a python one
        self.numAttrs_m.append(posDims + velDims)
        self.numParts_m.append(numParts)
        self.nameParts_m.append(name)

    def createSatellite(self, particleInd, altitude, upwardFacing, name):
        nameBuf = ctypes.create_string_buffer(bytes(name, encoding='utf-8'))
        self.simDLL_m.createSatelliteAPI(self.simulationptr, particleInd, altitude, upwardFacing, nameBuf)
        self.nameSat_m.append(name)
        self.satPartInd_m.append(particleInd)


    #One liner functions
    def getTime(self):
        return self.simDLL_m.getSimulationTimeAPI(self.simulationptr)
    
    
    #Pointer one liners (one liners in CPP obv, not Python)
    def getParticleDataFromCPP(self, origData=False):
        ret = []
        partattr = []
        partdbl = []
        for iii in range(self.numPartTypes_m):
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
    def getBFieldatS(self, s, time):
        return self.simDLL_m.getBFieldAtSAPI(self.simulationptr, s, time)

    def getEFieldatS(self, s, time):
        return self.simDLL_m.getEFieldAtSAPI(self.simulationptr, s, time)

    def fieldsAtAllZ(self, time, bins, binsize, z0):
        B_z = []
        E_z = []
        B_E_z_dim = []
        for iii in range(bins):
            B_z.append(self.getBFieldatS(z0 + binsize * iii, time))
            E_z.append(self.getEFieldatS(z0 + binsize * iii, time))
            B_E_z_dim.append(z0 + binsize * iii)
        return [B_z, E_z, B_E_z_dim]

    #Simulation management
    #Functions are prepended with __ because the intent is to simply runSim which will call them all
    #however if more refined control is needed, call them one by one and ignore runSim
    def initializeSimulation(self):
        self.simDLL_m.initializeSimulationAPI(self.simulationptr)

    def iterateSimulation(self, numberOfIterations, itersBtwCouts):
        print("Number Of Iterations: ", numberOfIterations)
        self.simDLL_m.iterateSimulationAPI(self.simulationptr, numberOfIterations, itersBtwCouts)

    def freeGPUMemory(self):
        self.simDLL_m.freeGPUMemoryAPI(self.simulationptr)

    def terminateSimulation(self):
        self.simDLL_m.terminateSimulationAPI(self.simulationptr)

    #Satellite functions
    def getNumberOfSatellites(self):
        return self.simDLL_m.getNumberOfSatellitesAPI(self.simulationptr)

    def getSatelliteData(self): #Need to add case where satellite captures more or less than the number of particles
        if self.numAttrs_m == []:
            self.getSimChars()

        satptr = [] #constructs array of double pointers so z value can be checked before recording data
        for jjj in range(self.numSats_m):
            attrptr = []
            for kk in range(self.numAttrs_m[self.satPartInd_m[jjj]] + 2):
                attrptr.append(self.simDLL_m.getSatelliteDataPointersAPI(self.simulationptr, jjj, kk))
            satptr.append(attrptr)
        
        satsdata = []
        for jjj in range(self.numSats_m):
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
    print("SimulationAPI.py is not meant to be called as main.  Run simulation.py and that will import this.")
