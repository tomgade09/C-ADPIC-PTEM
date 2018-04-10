import ctypes, os

print(os.getcwd())

MASS_ELEC = 9.10938356e-31
dllLoc = "./../../lib/geoplasmasim.dll"
simDLL = ctypes.CDLL(dllLoc)

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

PDptr = ctypes.c_void_p
saveBuf = ctypes.create_string_buffer(bytes("./out/", encoding='utf-8'))
attrBuf = ctypes.create_string_buffer(bytes("vpara,vperp,s", encoding='utf-8'))
partBuf = ctypes.create_string_buffer(bytes("elec", encoding='utf-8'))
PDptr = simDLL.createParticleDistributionAPI(saveBuf, attrBuf, partBuf, MASS_ELEC)

simDLL.addPDEnergyRangeAPI(PDptr, 96, 0.5, 4.5, True)
simDLL.addPDPitchRangeAPI(PDptr, 3600, 0.0, 180.0, True)
simDLL.generatePDAPI(PDptr, 101322.378940846, 19881647.2473464)
simDLL.writePDAPI(PDptr)

PDptr = None