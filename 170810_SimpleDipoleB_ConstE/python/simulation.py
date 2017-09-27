import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import math
import ctypes
import time
import particleFileTest

import os, sys, inspect, shutil
pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(pyfiledir)
dtg = '/' + time.strftime("%y%m%d") + "_" + time.strftime("%H.%M")
savedir = './../distgraphs' + dtg
if (not(os.path.isdir(savedir))):
    os.makedirs(savedir)
    os.makedirs(savedir + '/particles_init')
    os.makedirs(savedir + '/particles_final')
os.chdir(savedir)
srcfile = '../../include/_simulationvariables.h'
shutil.copy(srcfile, './')

cppDLL = ctypes.CDLL('./../../vs/x64/Release/170810_SimpleDipoleB_ConstE.dll')
cppDLL.dllmainPyWrapper.argtypes = [ctypes.c_char_p]
cppDLL.dllmainPyWrapper.restype = ctypes.POINTER(ctypes.c_double)

cppDLL.readDblBin.argtypes = [ctypes.c_char_p, ctypes.c_long]
cppDLL.readDblBin.restype = ctypes.POINTER(ctypes.c_double)

def cppDLLTest():
    v_e_para = []
    v_e_perp = []
    z_e = []
    inplay_e = []

    v_i_para = []
    v_i_perp = []
    z_i = []
    inplay_i = []

    B_z = []
    E_z = []
    B_E_z_dim = []

    dir_c = ctypes.create_string_buffer(bytes(savedir, encoding='utf-8')) #left in as an example - the dll doesn't need a save dir for now
    dllmainreturn = cppDLL.dllmainPyWrapper(dir_c) #Python changes directory to the newly created folder - making the above unnecessary

    electrons = int(dllmainreturn[0])
    ions = int(dllmainreturn[electrons * 3 + 1])
    length = (electrons + ions) * 3 + 2
    B_E_length = int(dllmainreturn[length])

    print("Py : "+str(electrons)+" "+str(ions)+" "+str(length))

    for iii in range(electrons):
        v_e_para.append(dllmainreturn[iii + 1])
        v_e_perp.append(dllmainreturn[iii + electrons + 1])
        z_e.append(dllmainreturn[iii + 2 * electrons + 1])
        
    for iii in range(ions):
        v_i_para.append(dllmainreturn[iii + 3 * electrons + 2])
        v_i_perp.append(dllmainreturn[iii + 3 * electrons + ions + 2])
        z_i.append(dllmainreturn[iii + 3 * electrons + 2 * ions + 2])

    for iii in range(B_E_length):
        B_z.append(dllmainreturn[length + 1 + iii])
        E_z.append(dllmainreturn[length + 1 + B_E_length + iii])
        B_E_z_dim.append(dllmainreturn[length + 1 + 2 * B_E_length + iii])

    return v_e_para, v_e_perp, z_e, v_i_para, v_i_perp, z_i, B_z, E_z, B_E_z_dim

def plotNormParticles(v_e_para, v_e_perp, v_i_para, v_i_perp, z_e, z_i, B_z, E_z, B_E_z_dim):
    plt.figure(1)
    plt.plot(v_e_para, v_e_perp, '.')
    #plt.axis([-5,5,0.0,5])
    plt.title('Electrons')
    plt.xlabel('Vpara (Re / s)')
    plt.ylabel('Vperp (Re / s)')
    plt.savefig('electrons.png')

    plt.figure(2)
    plt.plot(v_i_para, v_i_perp, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('Ions')
    plt.xlabel('Vpara')
    plt.ylabel('Vperp')
    plt.savefig('ions.png')

    plt.figure(3)
    plt.plot(v_e_perp, z_e, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('Electrons')
    plt.xlabel('Vperp (Re / s)')
    plt.ylabel('Z (Re)')
    plt.savefig('z_vprp_electrons.png')
    
    plt.figure(4)
    plt.plot(v_i_perp, z_i, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('Ions')
    plt.xlabel('Vperp (Re / s)')
    plt.ylabel('Z (Re)')
    plt.savefig('z_vprp_ions.png')

    plt.figure(5)
    plt.plot(v_e_para, z_e, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('Electrons')
    plt.xlabel('Vpara (Re / s)')
    plt.ylabel('Z (Re)')
    plt.savefig('z_vpra_electrons.png')
    
    plt.figure(6)
    plt.plot(v_i_para, z_i, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('Ions')
    plt.xlabel('Vpara (Re / s)')
    plt.ylabel('Z (Re)')
    plt.savefig('z_vpra_ions.png')

    plt.figure(7)
    plt.plot(B_E_z_dim, B_z, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('B Field')
    plt.xlabel('Z (Re)')
    plt.ylabel('B (T)')
    plt.savefig('B(z).png')

    plt.figure(8)
    plt.plot(B_E_z_dim, E_z, '.')
    #plt.axis([-0.3,0.7,0.0,0.4])
    plt.title('E Field')
    plt.xlabel('Z (Re)')
    plt.ylabel('E (V/m)')
    plt.savefig('E(z).png')

    plt.show()

if __name__ == '__main__':
    v_e_pr, v_e_pp, z_e, v_i_pr, v_i_pp, z_i, B_z, E_z, B_E_z_dim = cppDLLTest()
    #datadir = "./../170922_12.34/particles_final" #these three lines test distributions written to disk, ignore
    #v_e_pr, v_e_pp, z_e, v_i_pr, v_i_pp, z_i, B_z, E_z, B_E_z_dim = \
        #readFiles(datadir + "/e_vpara.bin", datadir + "/e_vperp.bin", datadir + "/e_z.bin", datadir + "/i_vpara.bin", datadir + "/i_vperp.bin", datadir + "/i_z.bin")
    plotNormParticles(v_e_pr, v_e_pp, v_i_pr, v_i_pp, z_e, z_i, B_z, E_z, B_E_z_dim)
