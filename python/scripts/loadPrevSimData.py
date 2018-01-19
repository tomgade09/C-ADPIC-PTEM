import numpy as np
import math
import os, sys, inspect
import matplotlib.pyplot as plt

pyscriptdir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(os.path.normpath(pyscriptdir + '/../'))
from simulation import *

def printArrayNames():
	print("Returned from simulation                    : finalDat, origDat, satDat")
	print("Containing all calc quantities              : elecpitches, elecenergies, [ionspitches, ionsenergies]")
	print("Calc quantities subsets [elec|ions]Ionsph...: Energies, Pitches, Vpara, Vperp, Epara, Eperp, Z, Times, LogTimes, LogEnergies")
	print("Calc quantities subsets [elec|ions]Magsph...: Energies, Pitches, Vpara, Vperp, Epara, Eperp, Z, Times, LogTimes, LogEnergies")
	print("Example                                     : elecMagsphLogEnergies, elecIonsphZ, ,ionsMagsphLogTimes")


dllLocation, savedir, rootdir, dtg = setupFolders()
prevSimDir = "./../180108_10.19.23/bins/"
sim = Simulation(dllLocation, rootdir, 0.01, 8.371e6, 4 * 6.371e6, 2.5, 1000, 0.0)

finalDat, origDat, satDat = sim.loadCompletedSimData(prevSimDir, "elec,ions", "vpara,vperp,z", "bottomElectrons,bottomIons,topElectrons,topIons", 100352, 2, 3, 4)

MASS_ELEC = 9.10938356e-31
MASS_PROT = 1.67262190e-27
J_PER_EV = 1.60218e-19

atan2vec = np.vectorize(math.atan2)

elecpitches = atan2vec(abs(np.array(origDat[0][1])), -np.array(origDat[0][0])) * 180 / math.pi
elecenergies = 0.5 * MASS_ELEC * (np.square(np.array(origDat[0][0])) + np.square(np.array(origDat[0][1]))) / J_PER_EV
elecIonsphEnergies = elecenergies[0:50176]
elecIonsphPitches = elecpitches[0:50176]
elecIonsphVpara = np.array(origDat[0][0][0:50176])
elecIonsphVperp = np.array(origDat[0][1][0:50176])
elecIonsphZ = np.array(origDat[0][2][0:50176])
elecMagsphEnergies = elecenergies[50176:100352]
elecMagsphPitches = elecpitches[50176:100352]
elecMagsphVpara = np.array(origDat[0][0][50176:100352])
elecMagsphVperp = np.array(origDat[0][1][50176:100352])
elecMagsphZ = np.array(origDat[0][2][50176:100352])

ionspitches = atan2vec(abs(np.array(origDat[1][1])), -np.array(origDat[1][0])) * 180 / math.pi
ionsenergies = 0.5 * MASS_PROT * (np.square(np.array(origDat[1][0])) + np.square(np.array(origDat[1][1]))) / J_PER_EV
ionsIonsphEnergies = ionsenergies[0:50176]
ionsIonsphPitches = ionspitches[0:50176]
ionsIonsphVpara = np.array(origDat[1][0][0:50176])
ionsIonsphVperp = np.array(origDat[1][1][0:50176])
ionsIonsphZ = np.array(origDat[1][2])
ionsMagsphEnergies = ionsenergies[50176:100352]
ionsMagsphPitches = ionspitches[50176:100352]
ionsMagsphVpara = np.array(origDat[1][0][50176:100352])
ionsMagsphVperp = np.array(origDat[1][1][50176:100352])
ionsMagsphZ = np.array(origDat[1][2][50176:100352])

elecIonsphTimesPre = np.array([[]])
elecMagsphTimesPre = np.array([[]])
ionsIonsphTimesPre = np.array([[]])
ionsMagsphTimesPre = np.array([[]])

elecIonsphTimes = np.array([])
elecMagsphTimes = np.array([])
ionsIonsphTimes = np.array([])
ionsMagsphTimes = np.array([])

elecIonsphTimesPre = [np.array(satDat[0][3][0:50176]), np.array(satDat[2][3][0:50176])]
elecMagsphTimesPre = [np.array(satDat[0][3][50176:100352]), np.array(satDat[2][3][50176:100352])]

ionsIonsphTimesPre = [np.array(satDat[1][3][0:50176]), np.array(satDat[3][3][0:50176])]
ionsMagsphTimesPre = [np.array(satDat[1][3][50176:100352]), np.array(satDat[3][3][50176:100352])]

for iii in range(50176):
	if elecIonsphTimesPre[0][iii] == -1:
		if elecIonsphTimesPre[1][iii] == -1:
			elecIonsphTimes = np.append(elecIonsphTimes, -1)
		else:
			elecIonsphTimes = np.append(elecIonsphTimes, elecIonsphTimesPre[1][iii])
	else:
		elecIonsphTimes = np.append(elecIonsphTimes, elecIonsphTimesPre[0][iii])
	
	if elecMagsphTimesPre[0][iii] == -1:
		if elecMagsphTimesPre[1][iii] == -1:
			elecMagsphTimes = np.append(elecMagsphTimes, -1)
		else:
			elecMagsphTimes = np.append(elecMagsphTimes, elecMagsphTimesPre[1][iii])
	else:
		elecMagsphTimes = np.append(elecMagsphTimes, elecMagsphTimesPre[0][iii])

	if ionsIonsphTimesPre[0][iii] == -1:
		if ionsIonsphTimesPre[1][iii] == -1:
			ionsIonsphTimes = np.append(ionsIonsphTimes, -1)
		else:
			ionsIonsphTimes = np.append(ionsIonsphTimes, ionsIonsphTimesPre[1][iii])
	else:
		ionsIonsphTimes = np.append(ionsIonsphTimes, ionsIonsphTimesPre[0][iii])
	
	if ionsMagsphTimesPre[0][iii] == -1:
		if ionsMagsphTimesPre[1][iii] == -1:
			ionsMagsphTimes = np.append(ionsMagsphTimes, -1)
		else:
			ionsMagsphTimes = np.append(ionsMagsphTimes, ionsMagsphTimesPre[1][iii])
	else:
		ionsMagsphTimes = np.append(ionsMagsphTimes, ionsMagsphTimesPre[0][iii])


elecIonsphTimesPre = np.array([])
elecMagsphTimesPre = np.array([])
ionsIonsphTimesPre = np.array([])
ionsMagsphTimesPre = np.array([])

def log10Excep(x):
	if x == -1:
		return 10
	return math.log10(x)

logVec = np.vectorize(log10Excep)

elecIonsphLogEnergies = logVec(elecIonsphEnergies)
elecMagsphLogEnergies = logVec(elecMagsphEnergies)
ionsIonsphLogEnergies = logVec(ionsIonsphEnergies)
ionsMagsphLogEnergies = logVec(ionsMagsphEnergies)

elecIonsphLogTimes = logVec(elecIonsphTimes)
elecMagsphLogTimes = logVec(elecMagsphTimes)
ionsIonsphLogTimes = logVec(ionsIonsphTimes)
ionsMagsphLogTimes = logVec(ionsMagsphTimes)

def elecVtoE(x):
	multfact = 1
	#if x < 0:
		#multfact = -1
	return multfact * 0.5 * MASS_ELEC * x**2 / J_PER_EV

def ionsVtoE(x):
	multfact = 1
	#if x < 0:
		#multfact = -1
	return multfact * 0.5 * MASS_PROT * x**2 / J_PER_EV
	
vToEelecVec = np.vectorize(elecVtoE)
vToEionsVec = np.vectorize(ionsVtoE)

elecIonsphEpara = vToEelecVec(elecIonsphVpara)
elecIonsphEperp = vToEelecVec(elecIonsphVperp)
elecMagsphEpara = vToEelecVec(elecMagsphVpara)
elecMagsphEperp = vToEelecVec(elecMagsphVperp)

elecIonsphLogEpara = logVec(elecIonsphEpara)
elecIonsphLogEperp = logVec(elecIonsphEperp)
elecMagsphLogEpara = logVec(elecMagsphEpara)
elecMagsphLogEperp = logVec(elecMagsphEperp)

ionsIonsphEpara = vToEionsVec(ionsIonsphVpara)
ionsIonsphEperp = vToEionsVec(ionsIonsphVperp)
ionsMagsphEpara = vToEionsVec(ionsMagsphVpara)
ionsMagsphEperp = vToEionsVec(ionsMagsphVperp)

printArrayNames()

plt.scatter(elecMagsphPitches, elecMagsphLogEnergies, c=elecMagsphLogTimes, cmap=plt.cm.jet, s=5)
plt.title("Magnetospheric-Source Electrons Pitch Angle vs Log Energy")
plt.xlabel("Pitch Angles (degrees)")
plt.ylabel("Log Energies (log eV)")
plt.colorbar(label="Log Time in Simulation (log s)")
plt.savefig("MagsphPitchVsLogEnergy.jpg", dpi=300)
plt.clf()

plt.scatter(elecIonsphPitches, elecIonsphLogEnergies, c=elecIonsphLogTimes, cmap=plt.cm.jet, s=5)
plt.title("Ionospheric-Source Electrons Pitch Angle vs Log Energy")
plt.xlabel("Pitch Angles (degrees)")
plt.ylabel("Log Energies (log eV)")
plt.colorbar(label="Log Time in Simulation (log s)")
plt.savefig("IonsphPitchVsLogEnergy.jpg", dpi=300)
plt.clf()

plt.scatter(elecMagsphLogEpara, elecMagsphLogEperp, c=elecMagsphLogTimes, cmap=plt.cm.jet, s=5)
plt.title("Magnetospheric-Source Electrons E parallel vs E perpendicular")
plt.xlabel("Log E_parallel (log eV)")
plt.ylabel("Log E_perpendicular (log eV)")
plt.colorbar(label="Log Time in Simulation (log s)")
plt.savefig("MagsphEparaVsEperp.jpg", dpi=300)
plt.clf()

plt.scatter(elecIonsphLogEpara, elecIonsphLogEperp, c=elecIonsphLogTimes, cmap=plt.cm.jet, s=5)
plt.title("Ionospheric-Source Electrons E parallel vs E perpendicular")
plt.xlabel("Log E_parallel (log eV)")
plt.ylabel("Log E_perpendicular (log eV)")
plt.colorbar(label="Log Time in Simulation (log s)")
plt.savefig("IonsphEparavsEperp.jpg", dpi=300)

