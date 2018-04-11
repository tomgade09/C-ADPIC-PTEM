import ctypes
from utilsAPI import *

name_this = ctypes.create_string_buffer(bytes('load dist', encoding='utf-8'))
name_other = ctypes.create_string_buffer(bytes('init dist', encoding='utf-8'))
fold_this = ctypes.create_string_buffer(bytes('./../../_in/data/Baseline Data', encoding='utf-8'))
fold_other = ctypes.create_string_buffer(bytes('./../../_dataout/180328_14.42.05/bins/particles_init', encoding='utf-8'))
attrs = ctypes.create_string_buffer(bytes('vpara,vperp,s', encoding='utf-8'))
partname = ctypes.create_string_buffer(bytes('elec', encoding='utf-8'))

this_d = ctypes.c_void_p
this_other = ctypes.c_void_p
this_d = simDLL.loadDistributionFromDiskAPI(name_this, fold_this, attrs, partname)
this_other = simDLL.loadDistributionFromDiskAPI(name_other, fold_other, attrs, partname)