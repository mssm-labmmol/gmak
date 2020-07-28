#!/usr/bin/python

import os 
import sys 
import numpy as np

# returns list with box lengths
def get_box (gro):
    gro = os.path.abspath(gro)
    fp = open(gro, "r")
    line = ""
    for line in fp:
        continue
    box = line.split()
    return box

def analyzeWrapper(inputData, prop, out):
    xtc = inputData['xtc']
    edr = inputData['edr']
    gro = inputData['gro']
    tpr = inputData['tpr']
    obtain_property(xtc, edr, gro, tpr, prop, out)
    return np.loadtxt(out, comments=['@','#'], usecols=(1,))

def obtain_property (xtc, edr, gro, tpr, name, output_file):
    xtc = os.path.abspath(xtc)
    edr = os.path.abspath(edr)
    gro = os.path.abspath(gro)
    tpr = os.path.abspath(tpr)
    output_file = os.path.abspath(output_file)

    path_of_preffix = '/'.join(output_file.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)
    
    # gamma is special
    if name == "gamma":
        # Obtain Lz.
        Lz = get_box(gro)[2]
        Lz_file = output_file + ".Lz"
        os.system("echo " + Lz + "> " + Lz_file)
        # Obtain Pxx, Pyy, Pzz series
        P = []
        for prop in ["Pres-XX", "Pres-YY", "Pres-ZZ"]: 
            os.system("echo " + prop + " | gmx energy -f " + edr + " -s " + tpr + " -o " + output_file + "." + prop + ".tmp")
            P.append(np.loadtxt(output_file + "." + prop + ".tmp.xvg", comments=['#','@'], usecols=(1,)))
        # Calculate gamma series.
        P = np.array(P)
        factor = 0.1 # conversion to mN/m
        gamma = factor * float(Lz) * 0.5 * (P[2,:] - 0.5 * (P[0,:] + P[1,:]))
        np.savetxt(output_file, np.column_stack((range(len(gamma)),gamma)))
        #for prop in ["Pres-XX", "Pres-YY", "Pres-ZZ"]: 
        #    os.system("rm " + prop + ".tmp.xvg")
        return

    if name == "density":
        prop = "Density"
    if name == "potential":
        prop = "Potential"
    if name == "pV":
        prop = "pV"
    if name == "volume":
        prop = "Volume"
    os.system("echo " + prop + " | gmx energy -f " + edr + " -s " + tpr + " -o " + output_file)

def obtain_gr (xtc, edr, tpr, ndx, g1, g2, output_file_preffix, cut=0, rmax=4.0, bin=0.002):
    xtc = os.path.abspath(xtc)
    edr = os.path.abspath(edr)
    tpr = os.path.abspath(tpr)
    ndx = os.path.abspath(ndx)
    output_file_preffix = os.path.abspath(output_file_preffix)

    path_of_preffix = '/'.join(output_file_preffix.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    output_list = []
    
    # dummy property, calculated only to extract the times
    os.system("echo \"Pot\" | gmx energy -f %s -s %s -o %s_dummy-times.xvg" % (edr,tpr,output_file_preffix))
    times = np.loadtxt("%s_dummy-times.xvg" % (output_file_preffix), usecols=(0,), comments=['@','#'])

    gr_data = []
    # calculate the rdf curves for all frames
    for i in range(len(times)):
        b = times[i]
        e = b
        os.system("gmx rdf -f %s -s %s -n %s -b %f -e %f -o %s_%d.xvg -ref %s -sel %s -cut %s -rmax %f -bin %f" % (xtc,tpr,ndx,b,e,output_file_preffix,i+1,g1,g2,cut,rmax,bin))
        idxs = np.loadtxt("%s_%d.xvg" % (output_file_preffix, i+1), usecols=(0,), comments=['@','#'])
        idxs_relevant = [j for j,x in enumerate(idxs) if ((x >= cut) and (x <= rmax))]
        gr_data.append(np.loadtxt("%s_%d.xvg" % (output_file_preffix, i+1),
            usecols=(1,), comments=['@','#'])[idxs_relevant])
        os.system("rm %s_%d.xvg" % (output_file_preffix, i+1))

    gr_data = np.array(gr_data)

    for i in range(len(gr_data[0,:])):
        np.savetxt("%s_dr_%d.xvg" % (output_file_preffix, i+1), np.column_stack((times,gr_data[:,i])))
        output_list.append("%s_dr_%d.xvg" % (output_file_preffix, i+1))

    return output_list


def obtain_polcorr (xtc, edr, gro, tpr, gas_dipole, gas_polarizability, output_file):
    xtc = os.path.abspath(xtc)
    edr = os.path.abspath(edr)
    gro = os.path.abspath(gro)
    tpr = os.path.abspath(tpr)
    output_file = os.path.abspath(output_file)

    path_of_preffix = '/'.join(output_file.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)
    Mtotfn = path_of_preffix + '/Mtot.xvg'

    AVOGADRO = 6.02e23
    DEBYE = 3.336e-30
    COUL = 8.987551e9

    # obtain dipoles
    os.system("echo 0 | gmx dipoles -f %s -s %s -o %s" % (xtc,tpr,Mtotfn))

    # get dummy times
    times = np.loadtxt(Mtotfn, comments=['@','#'], usecols=(0,))

    mu = np.loadtxt(Mtotfn, usecols=(4,), comments=['@','#'])
    corr = 0.5 * 1e-3 * AVOGADRO * DEBYE * DEBYE * (mu - gas_dipole) * (mu - gas_dipole) / (gas_polarizability * 1e-27 / COUL)

    outarray = np.column_stack ((times,corr))

    # write corrections to file
    np.savetxt(output_file, outarray)

    # remove all files at the end
    os.system("rm epsilon.xvg dipdist.xvg aver.xvg")

if __name__ == "__main__":
    xtc = "/home/yan/PROJECTS/mbar_par/systems/GEOMETRIES/MULTIPOLE-OPT/CAS-D/GRID_1/new-optimization-test-gas/544/sd/sd.xtc"
    gro = "/home/yan/PROJECTS/mbar_par/systems/GEOMETRIES/MULTIPOLE-OPT/CAS-D/GRID_1/new-optimization-test-gas/544/sd/sd.gro"
    edr = "/home/yan/PROJECTS/mbar_par/systems/GEOMETRIES/MULTIPOLE-OPT/CAS-D/GRID_1/new-optimization-test-gas/544/sd/sd.edr"
    tpr = "/home/yan/PROJECTS/mbar_par/systems/GEOMETRIES/MULTIPOLE-OPT/CAS-D/GRID_1/new-optimization-test-gas/544/sd/sd.tpr"
    ndx = "/home/yan/PROJECTS/mbar_par/systems/GEOMETRIES/MULTIPOLE-OPT/CAS-D/GRID_1/new-optimization-test-liq/544/md/oxygen.ndx"
    #obtain_property(xtc, edr, gro, tpr, "gamma", "density-test.xvg")
    #obtain_gr (xtc, edr, tpr, ndx, "OW", "OW", "./gr-files/gr_curve")
    gas_dipole = 1.85
    gas_polarizability = 0.00147
    obtain_polcorr (xtc, edr, gro, tpr, gas_dipole, gas_polarizability, "polcorr/polcorr.xvg")
