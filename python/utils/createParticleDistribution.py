import utilsAPI

MASS_ELEC = 9.10938356e-31

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