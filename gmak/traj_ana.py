import gmak.runcmd as runcmd
import os
import sys
import numpy as np
from gmak.config import ConfigVariables

# returns list with box lengths
def get_box (gro):
    gro = os.path.abspath(gro)
    fp = open(gro, "r")
    line = ""
    for line in fp:
        continue
    box = line.split()
    return [float(b) for b in box]

def analyzeWrapper(inputData, prop, out):
    xtc = inputData['xtc']
    edr = inputData['edr']
    gro = inputData['gro']
    tpr = inputData['tpr']
    if os.path.isfile(out) and (os.path.getmtime(out) > os.path.getmtime(tpr)):
        # i.e., the property file is newer than the trajectory
        pass
    else:
        obtain_property(xtc, edr, gro, tpr, prop, out)
    return np.loadtxt(out, comments=['@','#'], usecols=(-1,))

def obtain_property (xtc, edr, gro, tpr, name, output_file):
    xtc = os.path.abspath(xtc)
    edr = os.path.abspath(edr)
    gro = os.path.abspath(gro)
    tpr = os.path.abspath(tpr)
    output_file = os.path.abspath(output_file)

    path_of_preffix = '/'.join(output_file.split('/')[0:-1])
    # create path if it does not exist
    runcmd.run("mkdir -p " + path_of_preffix)
    runcmd.run("echo " + name + " | " + ConfigVariables.gmx + " energy -f " + edr + " -s " + tpr + " -o " + output_file)


def obtain_polcorr (xtc, edr, gro, tpr, gas_dipole, gas_polarizability, output_file):
    xtc = os.path.abspath(xtc)
    edr = os.path.abspath(edr)
    gro = os.path.abspath(gro)
    tpr = os.path.abspath(tpr)
    output_file = os.path.abspath(output_file)

    path_of_preffix = '/'.join(output_file.split('/')[0:-1])
    # create path if it does not exist
    runcmd.run("mkdir -p " + path_of_preffix)
    Mtotfn = path_of_preffix + '/Mtot.xvg'

    AVOGADRO = 6.02e23
    DEBYE = 3.336e-30
    COUL = 8.987551e9

    # obtain dipoles
    runcmd.run("echo 0 | %s dipoles -f %s -s %s -o %s" % (ConfigVariables.gmx, xtc,tpr,Mtotfn))

    # get dummy times
    times = np.loadtxt(Mtotfn, comments=['@','#'], usecols=(0,))

    mu = np.loadtxt(Mtotfn, usecols=(4,), comments=['@','#'])
    corr = 0.5 * 1e-3 * AVOGADRO * DEBYE * DEBYE * (mu - gas_dipole) * (mu - gas_dipole) / (gas_polarizability * 1e-27 / COUL)

    outarray = np.column_stack ((times,corr))

    # write corrections to file
    np.savetxt(output_file, outarray)

    # remove all files at the end
    runcmd.run("rm epsilon.xvg dipdist.xvg aver.xvg")

