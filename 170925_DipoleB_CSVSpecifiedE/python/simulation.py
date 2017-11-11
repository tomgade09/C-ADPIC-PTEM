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

    #sim.initializeSimulation() #if granular control over execution order is desired...
    #sim.copyDataToGPU()
    #sim.iterateSimulation(1)
    #sim.iterateSimulation(10000)
    #sim.copyDataToHost()
    #sim.freeGPUMemory()
    #sim.prepareResults()
    #results = sim.getResultsfrom3D()

    electrons = len(results[0][0])
    ions = len(results[1][0])
    length = (electrons + ions) * sim.attr_m + 3
    print("Py : "+str(electrons)+" "+str(ions)+" "+str(length))

    fields = sim.fieldsAtAllZ(0.0, 10000, (10 - 8.371/6.371)/ 10000, 8.371/6.371)

    plotAllParticles(results[0][0], results[0][1], results[0][2], results[1][0], results[1][1], \
        results[1][2], fields[0], fields[1], fields[2], False)

    #Eventually, read names from satellites and construct into array
    #plotSatelliteData(sim.getSatelliteData(), sim.satMsmts_m, sim.satNum_m, sim.dt_m, ['downwardElectrons', 'downwardIons', 'upwardElectrons', 'upwardIons'])

    sim.terminateSimulation170925()

    return

if __name__ == '__main__':
    simulationRunMain()