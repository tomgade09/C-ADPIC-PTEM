import ctypes
dllLoc = "./../../lib/geoplasmasim.dll"
simDLL = ctypes.CDLL(dllLoc)

#ParticleDistribution functions
simDLL.PDCreateAPI.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double)
simDLL.PDCreateAPI.restype = ctypes.c_void_p
simDLL.PDAddEnergyRangeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_bool)
simDLL.PDAddEnergyRangeAPI.restype = None
simDLL.PDAddPitchRangeAPI.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_bool)
simDLL.PDAddPitchRangeAPI.restype = None
simDLL.PDGenerateAPI.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
simDLL.PDGenerateAPI.restype = None
simDLL.PDFillExtraAttrsAPI.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p)
simDLL.PDFillExtraAttrsAPI.restype = None
simDLL.PDWriteAPI.argtypes = (ctypes.c_void_p,)
simDLL.PDWriteAPI.restype = None

#DistributionFromDisk functions
simDLL.DFDLoadAPI.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double)
simDLL.DFDLoadAPI.restype = ctypes.c_void_p
simDLL.DFDDataAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
simDLL.DFDDataAPI.restype = ctypes.POINTER(ctypes.c_double)
simDLL.DFDPrintAPI.argtypes = (ctypes.c_void_p, ctypes.c_int)
simDLL.DFDPrintAPI.restype = None
simDLL.DFDPrintDiffAPI.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int)
simDLL.DFDPrintDiffAPI.restype = None
simDLL.DFDZeroesAPI.argtypes = (ctypes.c_void_p,)
simDLL.DFDZeroesAPI.restype = None
simDLL.DFDCompareAPI.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
simDLL.DFDCompareAPI.restype = None
simDLL.DFDDeleteAPI.argtypes = (ctypes.c_void_p,)
simDLL.DFDDeleteAPI.restype = None