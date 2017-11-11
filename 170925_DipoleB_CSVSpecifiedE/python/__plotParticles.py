import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os, sys

def plotXY(xdata, ydata, title, xlabel, ylabel, filename, showplot=False):
    plt.plot(xdata, ydata, '.')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename)

    if showplot:
        plt.show()


def plotAllParticles(v_e_para, v_e_perp, z_e, v_i_para, v_i_perp, z_i, B_z, E_z, B_E_z_dim, showplot=False):
    plt.figure(1)
    plotXY(v_e_para, v_e_perp, 'Electrons', 'Vpara (Re / s)', 'Vperp (Re / s)', 'electrons.png')

    plt.figure(2)
    plotXY(v_i_para, v_i_perp, 'Ions', 'Vpara', 'Vperp', 'ions.png')

    plt.figure(3)
    plotXY(v_e_perp, z_e, 'Electrons', 'Vperp (Re / s)', 'Z (Re)', 'z_vprp_electrons.png')
    
    plt.figure(4)
    plotXY(v_i_perp, z_i, 'Ions', 'Vperp (Re / s)', 'Z (Re)', 'z_vprp_ions.png')

    plt.figure(5)
    plotXY(v_e_para, z_e, 'Electrons', 'Vpara (Re / s)', 'Z (Re)', 'z_vpra_electrons.png')
    
    plt.figure(6)
    plotXY(v_i_para, z_i, 'Ions', 'Vpara (Re / s)', 'Z (Re)', 'z_vpra_ions.png')

    plt.figure(7)
    plotXY(B_E_z_dim, B_z, 'B Field', 'Z (Re)', 'B (T)', 'B(z).png')

    plt.figure(8)
    plotXY(B_E_z_dim, E_z, 'E Field', 'Z (Re)', 'E (V/m)', 'E(z).png')

    if showplot:
        plt.show()


def plotSatelliteData(dataArray4D, numberMsmts, numberSats, dt, satNamesTuple, showplot=False):
    #need to add satellite plotting here
    for iii in range(numberMsmts):
        for jjj in range(numberSats):
            if (not(os.path.isdir('./' + satNamesTuple[jjj]))):
                os.makedirs('./' + satNamesTuple[jjj])
            os.chdir('./' + satNamesTuple[jjj])
            if (not(os.path.isdir('./' + iii * 1000 * dt))):
                os.makedirs('./' + iii)
            os.chdir('./' + iii)
            
            plt.figure(1)
            plotXY(dataArray4D[iii][jjj][0], dataArray4D[iii][jjj][1], 'Particles', 'Vpara (Re / s)', 'Vperp (Re / s)', 'electrons.png')
            plt.figure(2)
            plotXY(dataArray4D[iii][jjj][1], dataArray4D[iii][jjj][2], 'Particles', 'Vperp (Re / s)', 'Z (Re)', 'z_vprp_electrons.png')
            plt.figure(3)
            plotXY(dataArray4D[iii][jjj][0], dataArray4D[iii][jjj][2], 'Particles', 'Vpara (Re / s)', 'Z (Re)', 'z_vpra_electrons.png')
            
        os.chdir('./../../')