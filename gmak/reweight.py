import gmak.runcmd as runcmd

import os
import sys
from gmak.config import ConfigVariables

def reweightWrapper(inputData, workdir):
    xtc     = inputData['xtc']
    gro     = inputData['gro']
    top     = inputData['top']
    mdp     = inputData['mdp']
    tpr     = inputData['tpr']
    return reweight(xtc, gro, top, mdp, tpr, workdir)

def reweight (xtc, gro, top, mdp, tpr, workdir):
    xtc = os.path.abspath(xtc)
    gro = os.path.abspath(gro)
    top = os.path.abspath(top)
    mdp = os.path.abspath(mdp)
    tpr = os.path.abspath(tpr)
    workdir = os.path.abspath(workdir)
    check_log_loc = workdir + "/reweight.log"

    rw_tpr = workdir + "/reweight.tpr"
    rw_deffnm = workdir + "/reweight"

    runcmd.run("mkdir -p " + workdir)

    if (os.path.isfile(check_log_loc)) and (os.path.getmtime(check_log_loc) > os.path.getmtime(tpr)):
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
    output_files['log'] = check_log_loc
    return output_files
