import time, csv
import os, sys, inspect, shutil

pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(os.path.normpath(pyfiledir + '/../../SimTools/python/SimulationClass/'))

from __SimulationAPI170925 import *
from __plotParticles import *

#Setting up folders, changing directory
os.chdir(pyfiledir)
rootdir = os.path.dirname(os.path.abspath('./'))
dtg = '/' + time.strftime("%y%m%d") + "_" + time.strftime("%H.%M.%S")
savedir = './../distgraphs' + dtg
if (not(os.path.isdir(savedir))):
    os.makedirs(savedir)
    os.makedirs(savedir + '/bins/particles_init')
    os.makedirs(savedir + '/bins/particles_final')
    os.makedirs(savedir + '/bins/satellites')
    os.makedirs(savedir + '/graphs/allparticles')
    os.makedirs(savedir + '/graphs/EBfields')
    os.makedirs(savedir + '/graphs/satellites')
os.chdir(savedir)

srcfile = '../../include/_simulationvariables.h'
shutil.copy(srcfile, './')

def simulationRunMain():
    dllLocation = './../../vs/x64/Release/170925_DipoleB_CSVSpecifiedE.dll'
    print("================  SIMULATION ", dtg, " ================")

    sim = Simulation(rootdir, dllLocation)
    results = sim.runSim(5000)
    orig = sim.getOriginalsfrom3D()
    satDat = sim.getSatelliteData()
    fields = sim.fieldsAtAllZ(0.0, 4000, (sim.simMax_m - sim.simMin_m) / (6.371e6 * 4000), sim.simMin_m)

    sim.logWriteEntry('Python', 'Done getting data.  Plotting.')

    plotAllParticles(results[0][0], results[0][1], results[0][2], results[1][0], results[1][1], \
        results[1][2], fields[0], fields[1], fields[2], False)

    sim.logWriteEntry('Python', 'Plotting satellite data.')

    #Eventually, read names from satellites and construct into array
    plotSatelliteData(satDat, sim.satMsmts_m, sim.satNum_m, sim.dt_m, ['downwardElectrons', 'downwardIons', 'upwardElectrons', 'upwardIons'])

    sim.logWriteEntry('Python', 'Done plotting data.  Terminating simulation.')

    #saveEscapedParticlesAndTimeToCSV(orig, satDat)

    sim.terminateSimulation170925()

    return

if __name__ == '__main__':
    simulationRunMain()