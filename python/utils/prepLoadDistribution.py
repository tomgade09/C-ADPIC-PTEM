import ctypes
from utilsAPI import *

name_this = ctypes.create_string_buffer(bytes('init dist', encoding='utf-8'))
name_other = ctypes.create_string_buffer(bytes('final dist', encoding='utf-8'))
fold_this = ctypes.create_string_buffer(bytes('./../../_dataout/180416_15.35.35/bins/particles_init', encoding='utf-8'))
fold_other = ctypes.create_string_buffer(bytes('./../../_dataout/180416_15.35.35/bins/satellites', encoding='utf-8'))
attrs = ctypes.create_string_buffer(bytes('vpara,vperp,s', encoding='utf-8'))
partname = ctypes.create_string_buffer(bytes('elec', encoding='utf-8'))
satname = ctypes.create_string_buffer(bytes('4e6ElecUp', encoding='utf-8'))

this_d = ctypes.c_void_p
this_other = ctypes.c_void_p
this_d = simDLL.loadDistributionFromDiskAPI(name_this, fold_this, attrs, partname)
this_other = simDLL.loadDistributionFromDiskAPI(name_other, fold_other, attrs, satname)