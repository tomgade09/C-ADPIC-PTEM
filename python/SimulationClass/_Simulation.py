import ctypes

class _SimulationCDLL: #acts as parent class for actual Simulation class - removes some of the messy argtypes/restype functions out to make it cleaner and easier to read
    def __init__(self, dllLoc):
        self.dllLoc_m = dllLoc
        
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
        self.simDLL_m.getParticleNameAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getParticleNameAPI.restype = ctypes.c_char_p
        self.simDLL_m.getSatelliteNameAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
        self.simDLL_m.getSatelliteNameAPI.restype = ctypes.c_char_p
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
        self.simDLL_m.getBFieldAtSAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.getBFieldAtSAPI.restype = ctypes.c_double
        self.simDLL_m.getEFieldAtSAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
        self.simDLL_m.getEFieldAtSAPI.restype = ctypes.c_double
        
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
        self.simDLL_m.loadCompletedSimDataAPI.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p)
        self.simDLL_m.loadCompletedSimDataAPI.restype = None
        self.simDLL_m.setupNormalSimulationAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p)
        self.simDLL_m.setupNormalSimulationAPI.restype = None
        self.simDLL_m.runNormalSimulationAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
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