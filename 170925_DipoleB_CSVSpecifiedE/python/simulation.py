import time
from __plotParticles import *
from __SimulationAPI import *

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

def simulationRunMain():
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

    simpointer = initializeSimCPP()
    copyDataToGPUCPP(simpointer)
    iterateSimulationCPP(simpointer, 10000)
    copyDataToHostCPP(simpointer)
    terminateSimulationCPP(simpointer)
    dllmainreturn = returnResultsCPP(simpointer)

    attributesTracked = int(dllmainreturn[1])
    electrons = int(dllmainreturn[2])
    ions = int(dllmainreturn[electrons * attributesTracked + 3])
    length = (electrons + ions) * attributesTracked + 4

    print("Py : "+str(electrons)+" "+str(ions)+" "+str(length))

    for iii in range(electrons):
        v_e_para.append(dllmainreturn[iii + 3])
        v_e_perp.append(dllmainreturn[iii + electrons + 3])
        z_e.append(dllmainreturn[iii + 2 * electrons + 3])
        
    for iii in range(ions):
        v_i_para.append(dllmainreturn[iii + 3 * electrons + 4])
        v_i_perp.append(dllmainreturn[iii + 3 * electrons + ions + 4])
        z_i.append(dllmainreturn[iii + 3 * electrons + 2 * ions + 4])

    binsize = (10 - 8.371/6.371)/ 1000
    for iii in range(1000):
        B_z.append(getBatZCPP(8.371/6.371 + binsize * iii, 0.0))
        E_z.append(getEatZCPP(None, 8.371/6.371 + binsize * iii, 0.0))
        B_E_z_dim.append(binsize * iii)

    return v_e_para, v_e_perp, z_e, v_i_para, v_i_perp, z_i, B_z, E_z, B_E_z_dim

if __name__ == '__main__':
    v_e_pr, v_e_pp, z_e, v_i_pr, v_i_pp, z_i, B_z, E_z, B_E_z_dim = simulationRunMain()
    plotParticles(v_e_pr, v_e_pp, v_i_pr, v_i_pp, z_e, z_i, B_z, E_z, B_E_z_dim)
