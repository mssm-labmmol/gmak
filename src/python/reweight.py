import runcmd

import os
import sys
from config import ConfigVariables

def reweightWrapper(inputData, workdir):
    xtc     = inputData['xtc']
    gro     = inputData['gro']
    top     = inputData['top']
    mdp     = inputData['mdp']
    return reweight(xtc, gro, top, mdp, workdir)

def reweight (xtc, gro, top, mdp, workdir):
    xtc = os.path.abspath(xtc)
    gro = os.path.abspath(gro)
    top = os.path.abspath(top)
    mdp = os.path.abspath(mdp)
    workdir = os.path.abspath(workdir)
    check_log_loc = workdir + "/reweight.log"

    rw_tpr = workdir + "/reweight.tpr"
    rw_deffnm = workdir + "/reweight"

    runcmd.run("mkdir -p " + workdir)

    if (os.path.isfile(check_log_loc)):
        #print "log from reweight already exists, will not perform it"
        #print "if it did not exist, would issue commands:"
        #print("gmx grompp -f %s -c %s -o %s -p %s -maxwarn 5" % (mdp,gro,rw_tpr,top))
        #print("gmx mdrun -rerun %s -s %s -deffnm %s" % (xtc,rw_tpr,rw_deffnm))
        pass
    else:
        runcmd.run("%s grompp -f %s -c %s -o %s -p %s -maxwarn 5" % (ConfigVariables.gmx, mdp,gro,rw_tpr,top))
        runcmd.run("%s mdrun -rerun %s -s %s -deffnm %s" % (ConfigVariables.gmx, xtc,rw_tpr,rw_deffnm))

    output_files = {}
    output_files['xtc'] = xtc
    output_files['edr'] = rw_deffnm + ".edr"
    output_files['tpr'] = rw_deffnm + ".tpr"
    output_files['trr'] = rw_deffnm + ".trr"
    output_files['gro'] = gro
    return output_files
