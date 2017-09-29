import time
from __plotParticles import *
from __SimulationAPI170925 import *

#Setting up folders, changing directory
import os, sys, inspect, shutil
pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(pyfiledir)
rootdir = os.path.dirname(os.path.abspath('./')) #not sure why this works (as opposed to './../'), but it does.  It should be the earlier due to the chdir to one level above
dtg = '/' + time.strftime("%y%m%d") + "_" + time.strftime("%H.%M")
savedir = './../distgraphs' + dtg
if (not(os.path.isdir(savedir))):
    os.makedirs(savedir)
    os.makedirs(savedir + '/particles_init')
    os.makedirs(savedir + '/particles_final')
os.chdir(savedir)

srcfile = '../../include/_simulationvariables.h'
shutil.copy(srcfile, './')

dllLocation = './../../vs/x64/Release/170925_DipoleB_CSVSpecifiedE.dll'

def simulationRunMain():

    sim = Simulation(rootdir, dllLocation)

    results = sim.runSim(10000)

    electrons = len(results[0][0])
    ions = len(results[1][0])
    length = (electrons + ions) * sim.attr_m + 3
    print("Py : "+str(electrons)+" "+str(ions)+" "+str(length))

    fields = sim.fieldsAtAllZ(0.0, 1000, (10 - 8.371/6.371)/ 1000, 8.371/6.371)

    #wrong_para = 0
    #wrong_perp = 0
    #wrong_z = 0

    #for iii in range(electrons):
        #if (v_e_para[iii] != v_e_para_3D[iii]):
            #wrong_para += 1
        #if (v_e_perp[iii] != v_e_perp_3D[iii]):
            #wrong_perp += 1
        #if (z_e[iii] != z_e_3D[iii]):
            #wrong_z += 1

    #if (not(wrong_para == 0) or not(wrong_perp == 0) or not(wrong_z == 0)):
        #print("electrons wrong: ", wrong_para, wrong_perp, wrong_z)

    #wrong_para = 0
    #wrong_perp = 0
    #wrong_z = 0

    #for iii in range(ions):
        #if (v_i_para[iii] != v_i_para_3D[iii]):
            #wrong_para += 1
        #if (v_i_perp[iii] != v_i_perp_3D[iii]):
            #wrong_perp += 1
        #if (z_i[iii] != z_i_3D[iii]):
            #wrong_z += 1

    #if (not(wrong_para == 0) or not(wrong_perp == 0) or not(wrong_z == 0)):
        #print("ions wrong: ", wrong_para, wrong_perp, wrong_z)

    return results[0][0], results[0][1], results[0][2], results[1][0], results[1][1], \
        results[1][2], fields[0], fields[1], fields[2]

if __name__ == '__main__':
    v_e_pr, v_e_pp, z_e, v_i_pr, v_i_pp, z_i, B_z, E_z, B_E_z_dim = simulationRunMain()
    plotParticles(v_e_pr, v_e_pp, v_i_pr, v_i_pp, z_e, z_i, B_z, E_z, B_E_z_dim, False)
