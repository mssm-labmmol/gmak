import sys
import gmak.runcmd as runcmd
import os
import subprocess
import re
from shutil import copyfile
from gmak.traj_ana import get_box
from gmak.mdputils import *
import gmak.logger as logger
from gmak.config import ConfigVariables

# gmx_nprocs is set via a cmd-line option '-np'
# It is the default number of processors in each gmx mdrun job.
mdrun_nprocs = -1

# Creates a call to 'gmx mdrun' with the appropriate number of processors.
def create_mdrun_call(nprocs):
    if (nprocs > 0):
        return ("%s mdrun -nt %d" % (ConfigVariables.gmx, nprocs))
    else:
        if (mdrun_nprocs > 0):
            return ("%s mdrun -nt %d" % (ConfigVariables.gmx, mdrun_nprocs))
        else:
            return ("%s mdrun" % ConfigVariables.gmx)

def read_nsteps_from_mdp (mdp):
    mu = mdpUtils()
    mu.parse_file(mdp)
    return mu.get_nsteps()

def check_simulation_state (workdir, label):
    check_log_loc = "%s/%s/%s.log" % (workdir, label, label)
    check_gro_loc = "%s/%s/%s.gro" % (workdir, label, label)
    check_cpt_loc = "%s/%s/%s.cpt" % (workdir, label, label)
    if (os.path.isfile(check_cpt_loc)):
        return 'FULL'
    elif (os.path.isfile(check_gro_loc)):
        return 'COMPLETE_MINIM'
    else:
        return 'NONE'

def extend_something (nsteps, workdir, label, nprocs=-1):
    logger.globalLogger.indent()
    logger.globalLogger.putMessage('Extending')
    logger.globalLogger.unindent()
    workdir = os.path.abspath(workdir)
    runcmd.run("%s convert-tpr -s %s/%s/%s.tpr -nsteps %d -o %s/%s/tmp.tpr" % (ConfigVariables.gmx, workdir, label, label, nsteps, workdir, label))
    runcmd.run("mv %s/%s/tmp.tpr %s/%s/%s.tpr" % (workdir, label, workdir, label, label))
    runcmd.run("%s -cpi %s/%s/%s.cpt -s %s/%s/%s.tpr -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label,
                                                                                    workdir, label, label,
                                                                                    workdir, label, label))

def simulate_something (conf, top, mdp, label, workdir, nprocs=-1):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top.path)
    mdp = os.path.abspath(mdp)
    workdir = os.path.abspath(workdir)
    runcmd.run("mkdir -p %s/%s" % (workdir, label))
    check_log_loc = "%s/%s/%s.log" % (workdir,label,label)
    simu_state = check_simulation_state(workdir, label)
    if (simu_state == 'COMPLETE_MINIM'):
        # This will only happen for minimization, I guess. We can just assume it finished alright with a warning.
        logger.globalLogger.putMessage('STEP {}: There is a log but no cpt; assumming a successful minimization.'.format(label))
    elif (simu_state == 'FULL'):
        logger.globalLogger.putMessage('STEP {}: Restarting from .cpt (possibly a complete simulation)'.format(label))
        runcmd.run("%s -s %s/%s/%s.tpr -cpi %s/%s/%s.cpt -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label, workdir, label, label, workdir, label, label))
    elif (simu_state == 'NONE'):
        logger.globalLogger.putMessage('STEP {}: Simulating from start'.format(label))
        command = "%s grompp -maxwarn 5 -f %s -c %s -p %s -o %s/%s/%s.tpr" % (ConfigVariables.gmx, mdp, conf, top, workdir, label, label)
        runcmd.run(command)
        command = "%s -s %s/%s/%s.tpr -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label, workdir, label, label)
        runcmd.run(command)


def fix_periodicity(in_conf, out_conf):
    runcmd.run("echo 0 | %s trjconv -f %s -s %s -o %s -pbc whole" %
               (ConfigVariables.gmx, in_conf, in_conf, out_conf))


def resize_box(in_conf, out_conf, box):
    runcmd.run("%s editconf -f %s -box %f %f %f -o %s" % (ConfigVariables.gmx,
                                                          in_conf,
                                                          float(box[0]),
                                                          float(box[1]),
                                                          float(box[2]),
                                                          out_conf))


def make_a_box(conf, nmols, box, outconf):
    conf = os.path.abspath(conf)
    outconf = os.path.abspath(outconf)
    # create box if it does not exist
    if not (os.path.isfile(outconf)):
        command = "%s insert-molecules -ci %s -nmol %d -box %f %f %f -o %s" % (
            ConfigVariables.gmx, conf,nmols,box[0],box[1],box[2],outconf)
        runcmd.run(command)


def simulate_protocol_solvation (conf, topo, mdps, nsteps, labels, workdir, simulate=True):
    topopath = topo
    confpath = conf
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    if (simulate):
        for i in range(len(mdps)):
            if i == 0:
                simulate_something (confpath, topopath, mdps[i], labels[i], workdir)
            elif i == len(mdps) - 1:
                if (read_nsteps_from_mdp(mdps[i]) == nsteps):
                    previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                    simulate_something (previous_conf, topopath, mdps[i], labels[i], workdir)
                else:
                    # Extend
                    extend_something(nsteps, workdir, labels[i])
            else:
                previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                simulate_something (previous_conf, topopath, mdps[i], labels[i], workdir)
    # create output dictionary
    output_dict['nsteps'] = nsteps
    output_dict['top']   =  topopath
    output_dict['dhdl']  =  workdir   +  "/"  +  labels[-1]  +  "/" + labels[-1] + ".xvg"
    for ext in ['xtc', 'tpr', 'trr', 'edr', 'gro', 'top']:
        output_dict[ext] = workdir + "/" + labels[-1] + "/" + labels[-1] + '.' + ext
    return output_dict


def simulate_protocol_general(conf, top, mdps, nsteps, labels, workdir):
    conf = os.path.abspath(conf)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (conf, top, mdps[i], labels[i], workdir)
        elif i == len(mdps) - 1:
            if (read_nsteps_from_mdp(mdps[i]) == nsteps):
                previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                simulate_something (previous_conf, top, mdps[i], labels[i], workdir)
            else:
                # Extend
                extend_something(nsteps, workdir, labels[i])
        else:
            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
            simulate_something (previous_conf, top, mdps[i], labels[i], workdir)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['trr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".trr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = top
    output_dict['nsteps'] = nsteps
    #
    return output_dict

