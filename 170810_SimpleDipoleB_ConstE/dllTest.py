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
def cppDLLTest():
    cppDLL = ctypes.CDLL('./../170810_SimpleDipoleB_ConstE.dll')
    cppDLL.dllmainPyWrapper.argtypes = None
    cppDLL.dllmainPyWrapper.restype = ctypes.POINTER(ctypes.c_double)

    v_e_para = []
    v_e_perp = []
    z_e = []
    inplay_e = []

    v_i_para = []
    v_i_perp = []
    z_i = []
    inplay_i = []

    dllmainreturn = cppDLL.dllmainPyWrapper()

    for iii in range(100000):
        v_e_para.append(dllmainreturn[iii])
        v_e_perp.append(dllmainreturn[iii+100000])
        z_e.append(dllmainreturn[iii+200000])
        #inplay_e.append(dllmainreturn[iii+300000])
        v_i_para.append(dllmainreturn[iii+300000])
        v_i_perp.append(dllmainreturn[iii+400000])
        z_i.append(dllmainreturn[iii+500000])
        #inplay_i.append(dllmainreturn[iii+700000])

    return v_e_para, v_e_perp, z_e, v_i_para, v_i_perp, z_i

def plotNormParticles(v_e_para, v_e_perp, v_i_para, v_i_perp):
    plt.figure(1)
    #plt.scatter(v_e_para, v_e_perp, '.')
    plt.plot(v_e_para, v_e_perp, '.')
    plt.title('Electrons')
    plt.ylabel('Vperp')
    plt.xlabel('Vpara')
    plt.savefig('electrons.png')

    plt.figure(2)
    #plt.scatter(v_i_para, v_i_perp, '.')
    plt.plot(v_i_para, v_i_perp, '.')
    plt.title('Ions')
    plt.ylabel('Vperp')
    plt.xlabel('Vpara')
    plt.savefig('ions.png')

    #plt.show()

if __name__ == '__main__':
    v_e_pr, v_e_pp, z_e, v_i_pr, v_i_pp, z_i = cppDLLTest()
    plotNormParticles(v_e_pr, v_e_pp, v_i_pr, v_i_pp)
