import matplotlib.pyplot as plt
import numpy as np
import math
import ctypes

import os, sys, inspect
pyfiledir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(pyfiledir)
savedir = './distgraphs'
if (not(os.path.isdir(savedir))):
    os.makedirs(savedir)
os.chdir(savedir)

####################Multivariable normal distribution example
#mean = [0, 0]
#covfefe = [[1, 0], [0, 1]]  # diagonal covariance
#x, y = np.random.multivariate_normal(mean, covfefe, 5000).T
#plt.plot(x, y, 'x')
#plt.axis('equal')
#plt.show()

####################Two plot, separate windows example
#plt.plot(range(10))
#plt.figure()
#plt.plot(range(10), 'ro-')
#plt.show()

#DLLEXPORT double** normalDistribution_v_z(int numOfParticles, double vmean, double vsigma, double zmean, double zsigma)
def cppNormalDist_v_z(numOfParticles, vmean, vsigma, zmean, zsigma):
    cppND = ctypes.CDLL('./../cppnormdist.dll')
    cppND.normalDistribution_v_z.argtypes = (ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    cppND.normalDistribution_v_z.restype = ctypes.POINTER(ctypes.c_double)

    arraySize = numOfParticles * 4

    C_DOUBLENA = ctypes.c_double * arraySize

    print("cpp entered")
    params = C_DOUBLENA()
    params = cppND.normalDistribution_v_z(numOfParticles, vmean, vsigma, zmean, zsigma)

    print("cpp exited.")

    pyParams = []

    for iii in range(4):
        pyParamsData = []
        for jjj in range(numOfParticles):
            pyParamsData.append(params[jjj * 4 + iii])
            #print(iii, jjj, "   ", params[jjj * 4 + iii])
        pyParams.append(pyParamsData)
    
    print("pyParams written.")

    return pyParams

def plotNormParticles(vpara, vperp, z, pitch):
    plt.figure(1)
    plt.plot(vperp, vpara, '.')
    plt.title('Vpara vs Vperp')
    plt.xlabel('Vperp')
    plt.ylabel('Vpara')
    plt.savefig('vparaVsvperp.png')

    plt.figure(2)
    plt.plot(vperp, z, '.')
    plt.title('Z vs Vperp')
    plt.xlabel('Vperp')
    plt.ylabel('Z')
    plt.savefig('zVsvperp.png')

    plt.figure(3)
    plt.plot(vpara, z, '.')
    plt.title('Z vs Vpara')
    plt.xlabel('Vpara')
    plt.ylabel('Z')
    plt.savefig('zVsvpara.png')

    plt.figure(4)
    plt.hist(vperp, bins=100)
    plt.title('Vperp Histogram')
    plt.xlabel('Vperp')
    plt.ylabel('Number')
    plt.savefig('vperpHist.png')

    plt.figure(5)
    plt.hist(vpara, bins=100)
    plt.title('Vpara Histogram')
    plt.xlabel('Vpara')
    plt.ylabel('Number')
    plt.savefig('vparaHist.png')

    plt.figure(6)
    plt.hist(z, bins=100)
    plt.title('Z Histogram')
    plt.xlabel('Z')
    plt.ylabel('Number')
    plt.savefig('zHist.png')

    plt.figure(7)
    plt.hist2d(vpara, vperp, bins=100, range=[[-4,4],[-4,4]])
    plt.colorbar()
    plt.title('Vperp vs Vpara Histogram')
    plt.xlabel('Vpara')
    plt.ylabel('Vperp')
    plt.savefig('vperpVsvpara2DHist.png')

    plt.figure(8)
    plt.hist(pitch, bins=100, range=[-185,185])
    plt.title('Pitch Angle Histogram')
    plt.xlabel('Pitch Angle')
    plt.ylabel('Number')
    plt.savefig('pitchHist.png')

    #plt.show()

def createNormParticles(number, mean=[0,0,0], cov=[[1,0,0],[0,1,0],[0,0,1]]):
    vperp, vpara, z = np.random.multivariate_normal(mean, cov, number).T
    pitch = []

    for i in range(len(vperp)):
        pitch.append(math.atan2(vperp[i], vpara[i]) * 180 / math.pi)

    plotNormParticles(vperp, vpara, z, pitch)

    return [vpara, vperp, z, pitch]

if __name__ == '__main__':
    #createNormParticles(100000, mean=[0,0,2], cov=[[1,0,0],[0,1,0],[0,0,0.25]])
    pyParams = cppNormalDist_v_z(100000, 0.0, 1.0, 2.0, 0.5)
    plotNormParticles(pyParams[0], pyParams[1], pyParams[2], pyParams[3])
