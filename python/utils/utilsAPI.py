import ctypes
dllLoc = "./../../lib/geoplasmasim.dll"
simDLL = ctypes.CDLL(dllLoc)

#ParticleDistribution functions
simDLL.createParticleDistributionAPI.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double)
simDLL.createParticleDistributionAPI.restype = ctypes.c_void_p
simDLL.addPDEnergyRangeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_bool)
simDLL.addPDEnergyRangeAPI.restype = None
simDLL.addPDPitchRangeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_bool)
simDLL.addPDPitchRangeAPI.restype = None
simDLL.generatePDAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
simDLL.generatePDAPI.restype = None
simDLL.writePDAPI.argtypes = (ctypes.c_void_p,)
simDLL.writePDAPI.restype = None

#DistributionFromDisk functions
simDLL.loadDistributionFromDiskAPI.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p)
simDLL.loadDistributionFromDiskAPI.restype = ctypes.c_void_p
simDLL.DistFromDiskPrintAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
simDLL.DistFromDiskPrintAPI.restype = None
simDLL.DistFromDiskPrintDiffAPI.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int)
simDLL.DistFromDiskPrintDiffAPI.restype = None
simDLL.DistFromDiskZeroesAPI.argtypes = (ctypes.c_void_p,)
simDLL.DistFromDiskZeroesAPI.restype = None
simDLL.DistFromDiskCompareAPI.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
simDLL.DistFromDiskCompareAPI.restype = None
simDLL.deleteDistFromDisk.argtypes = (ctypes.c_void_p,)
simDLL.deleteDistFromDisk.restype = None