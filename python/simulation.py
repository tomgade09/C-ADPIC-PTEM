import time, csv
import os, sys, inspect, shutil

pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(os.path.normpath(pyfiledir + '/SimulationClass/'))

from AlfvenLUT import *
from __plotParticles import *

#Setting up folders, changing directory
def setupFolders():
    os.chdir(pyfiledir)
    rootdir = os.path.dirname(os.path.abspath('./'))
    dtg = '/' + time.strftime("%y%m%d") + "_" + time.strftime("%H.%M.%S")
    savedir = os.path.abspath(rootdir + '/_dataout' + dtg)
    if (not(os.path.isdir(savedir))):
        os.makedirs(savedir)
        os.makedirs(savedir + '/bins/particles_init')
        os.makedirs(savedir + '/bins/particles_final')
        os.makedirs(savedir + '/bins/satellites')
        os.makedirs(savedir + '/graphs/allparticles')
        os.makedirs(savedir + '/graphs/EBfields')
        os.makedirs(savedir + '/graphs/satellites')
    os.chdir(savedir)

    #srcfile = '../../include/_simulationvariables.h' #change this
    #shutil.copy(srcfile, './')

    dllLocation = os.path.abspath('./../../lib/geoplasmasim.dll')

    return dllLocation, savedir, rootdir, dtg

def simulationRunMain():
    dllLocation, savedir, rootdir, dtg = setupFolders()
    print("================  SIMULATION ", dtg, " ================")

    sim = Simulation(dllLocation, rootdir, 0.01, 2030837.49610366, 19881647.2473464, 2.5, 1000, 0.0)#, "ez.out") #need to pass in either height from Re/geocentric or s - right now it's s
    finalDat, origDat, satDat = sim.runSim(25000)
    #                        time, bins,           binsize,                           z0
    fields = sim.fieldsAtAllZ(0.0, 4000, (sim.simMax_m - sim.simMin_m) / (4000), sim.simMin_m)
    for iii in range(len(fields[2])):
        fields[2][iii] /= 6.371e6

    sim.logWriteEntry('Python: Done getting data.  Plotting.')

    plotFields(fields[0], fields[1], fields[2], True)
    #plotAllParticles(results[0][0], results[0][1], results[0][2], results[1][0], results[1][1], \
        #results[1][2], fields[0], fields[1], fields[2], False)

    sim.logWriteEntry('Python: Plotting satellite data.')

    #Eventually, read names from satellites and construct into array
    #plotSatelliteData(satDat, 1, sim.satNum_m, sim.dt_m, ['downwardElectrons', 'downwardIons', 'upwardElectrons', 'upwardIons'])

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