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

dllLocation = './../../vs/x64/Release/170925_DipoleB_CSVSpecifiedE.dll'

print("================  SIMULATION ", dtg, " ================")

#TO-DO
# DONE - Top at 4 Re
# DONE, NEEDS BETTER PHYSICS - Particles out top are lost completely (done), out bottom - scattering probability, distribution
# DONE, BUT FACTOR IN DENSITY? - Make distribution maxwellian
# DONE, CONSIDER OTHER INJECTION SCHEMES- More particles, e, p injected over simulation (bump up to a million?, 10000 particles injected every iteration/so many hundredths of second?)
# DONE - Stop generating a distribution on CPU at beginning - take out a whole bunch of copy code, maybe even normal distribution code, still 0 the arrays
# DONE - Fix satellite code
# DONE - Log file with success messages
# DONE - Spit run times to Log File, maybe print some on screen
# DONE - Either verify normalization or just normalize all the values at the end - Fix botched up normalization system
# MAYBE NOT, ALLOWS TIME BASED CAPTURE, COULD ALSO CAPTURE SIM TIME - Consolidate satellite data into one array, then pass to python
# DONE - Satellite fix - remove zero values
# DONE, NEEDS A GOOD TEST - Fix LUT code to be one function - make LUT 2D on proc - like 1D with pointers to start of next dimension
# - Write Satellite values to disk along with times
# - Enable qsps, alfven via switch
# ALMOST DONE - Enable load particle distribution via switch
# - Remove frivolous messages
# - Maybe, change the convertMuToVperp and vice versa to pass in array pointers and execute on whatever is passed in
# - Replace double* arrays with std::vectors, pass with std::vectors::data (I think)

# - XLSX writing with xlnt library?
# - Some sort of error handling system instead of cout messages
#       - Class of some sort, some common error codes (indexing errors, values out of range, etc)
#       - Encapsulate everything in try, except blocks?
#       - 
# - Change over current error messages to log file
# - Validation with known values (possibly magnetic bottle) - Calculations checked against excel and are good
# - Restructure - Combine most of two simulation classes - make the two particle, two source, 1D code the generic one
#     Or perhaps:
#       - Simulation container with common elements
#       |-- 1 position dimension, 2 velocity, two particle
#       | |--- 1708XX - Const E Field
#       | |--- 170925 - LUT E Field
#       | |--- 17XXXX - E Field generating code included
#       |-- Others if needed
# - Possibly more common code compiled to DLLs to make recompiling less complex - call already made DLLs from sims (ex: FileIO, Simulation)
# - Pass in variables instead of having to recompile every time (CSV or possibly XML if you want to be fancy)
# - Photoelectrons
# - Fix sideways CSV Writing

def simulationRunMain():
    sim = Simulation(rootdir, dllLocation)
    results = sim.runSim(25000)
    orig = sim.getOriginalsfrom3D()
    satDat = sim.getSatelliteData()
    fields = sim.fieldsAtAllZ(0.0, 4000, (sim.simMax_m - sim.simMin_m) / (6.371e6 * 4000), sim.simMin_m)

    sim.logWriteEntry('Python', 'Done getting data.  Plotting.')

    #save4DDataToCSV(satDat, './CSV')
    #with open("./BEfields.csv", "w", newline='\n') as f:
        #csvwriter = csv.writer(f)
        #csvwriter.writerows(fields)

    #sim.logWriteEntry('Python', 'Done writing CSVs')

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