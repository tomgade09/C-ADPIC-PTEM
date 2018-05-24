import ctypes
from utilsAPI import *

name_orig = ctypes.create_string_buffer(bytes('orig dist', encoding='utf-8'))
name_comp = ctypes.create_string_buffer(bytes('comp dist', encoding='utf-8'))
fold_orig = ctypes.create_string_buffer(bytes('./../../_dataout/180423_09.34.25/bins/particles_final', encoding='utf-8'))
fold_comp = ctypes.create_string_buffer(bytes('./../../_dataout/180524_13.45.53/bins/particles_final', encoding='utf-8'))
attrs = ctypes.create_string_buffer(bytes('vpara,vperp,s', encoding='utf-8'))
origname = ctypes.create_string_buffer(bytes('elec', encoding='utf-8'))
compname = ctypes.create_string_buffer(bytes('elec', encoding='utf-8'))

dist_orig = ctypes.c_void_p
dist_comp = ctypes.c_void_p
dist_orig = simDLL.loadDistributionFromDiskAPI(name_orig, fold_orig, attrs, origname)
dist_comp = simDLL.loadDistributionFromDiskAPI(name_comp, fold_comp, attrs, compname)

simDLL.DistFromDiskCompareAPI(dist_orig, dist_comp)

print("Distributions:")
print("dist_orig: original distribution; dist_comp: distribution to compare")
print("")
print("Available Functions:")
print("simDLL.DistFromDiskDataAPI(distptr, ind of data) - return double* of data attr")
print("simDLL.DistFromDiskPrintAPI(distptr, ind of data) - print index of data")
print("simDLL.DistFromDiskPrintDiffAPI(origdistptr, compdistptr, ind of data) - print difference at index")
print("simDLL.DistFromDiskZeroesAPI(distptr) - print how many zeroes")
print("simDLL.DistFromDiskCompareAPI(origdistptr, compdistptr) - print various aspects comparing two dists")
print("simDLL.deleteDistFromDiskAPI(distptr) - delete, free memory of dist")
print("")
