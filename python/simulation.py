import time, csv
import os, sys, shutil

from __simulationvariables import *

sys.path.append(os.path.normpath(PYROOTDIR + '/SimulationClass/'))

from Simulation import *
from __plotParticles import *

#Setting up folders, changing directory
def setupFolders(): #moving to C++ eventually
    os.chdir(PYROOTDIR)
    dtg = time.strftime("%y%m%d") + "_" + time.strftime("%H.%M.%S")
    savedir = os.path.abspath(ROOTDIR + '/_dataout' + '/' + dtg)
    if (not(os.path.isdir(savedir))):
        os.makedirs(savedir)
        os.makedirs(savedir + '/bins/particles_init')
        os.makedirs(savedir + '/bins/particles_final')
        os.makedirs(savedir + '/bins/satellites')
        os.makedirs(savedir + '/graphs/allparticles')
        os.makedirs(savedir + '/graphs/EBfields')
        os.makedirs(savedir + '/graphs/satellites')
        os.makedirs(savedir + '/_chars')
    os.chdir(savedir)

    srcfile = PYROOTDIR + '/__simulationvariables.py' #change this
    shutil.copy(srcfile, './')

    return savedir, dtg

def simulationRunMain():
    savedir, dtg = setupFolders()
    print("================  SIMULATION ", dtg, " ================")

    sim = Simulation(DLLLOCATION, ROOTDIR, DT, MIN_S_SIM, MAX_S_SIM, INITIAL_T_ION_EV, INITIAL_T_MAG_EV)
    finalDat, origDat, satDat = sim.runSim(25, True)

    fields = sim.fieldsAtAllZ(0.0, 4000, (sim.simMax_m - sim.simMin_m) / (4000), sim.simMin_m)
    for iii in range(len(fields[2])):
        fields[2][iii] /= RADIUS_EARTH #Normalizes location where field is measured (makes graph reading easier)
    for iii in range(len(fields[0])):
        fields[0][iii] *= 1e9 #Displays B Field in nT

    sim.logWriteEntry('Python: Done getting data.  Plotting.')

    plotFields(fields[0], fields[1], fields[2], True)

    sim.logWriteEntry('Python: Done plotting data.  Terminating simulation.')

    sim.terminateSimulation()

    return


if __name__ == '__main__':
    simulationRunMain()

if __name__ != '__main__':
    print("Functions available in simulation.py:")
    print()
    print("setupFolders(): Sets up folders according to default structure.  Changes dir to the data root folder (/PROJECTROOT/distgraphs/DATE_TIME).  Returns dllLocation, savedir, rootdir, DATE_TIME string")
    print()
    print("simulationRunMain(): Calls setupFolders().  Runs the simulation with the variables specified in simulationvariables.py and ends the sim when done.  Returns nothing, and saves data to the default folder structure.")
    print()