import sys
import runcmd
import os
import subprocess
import re
from shutil import copyfile
from traj_ana import get_box
from mdputils import *
from logger import *
from config import ConfigVariables

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
    globalLogger.indent()
    globalLogger.putMessage('Extending')
    globalLogger.unindent()
    workdir = os.path.abspath(workdir)
    runcmd.run("%s convert-tpr -s %s/%s/%s.tpr -nsteps %d -o %s/%s/tmp.tpr" % (ConfigVariables.gmx, workdir, label, label, nsteps, workdir, label))
    runcmd.run("mv %s/%s/tmp.tpr %s/%s/%s.tpr" % (workdir, label, workdir, label, label))
    runcmd.run("%s -cpi %s/%s/%s.cpt -s %s/%s/%s.tpr -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label,
                                                                                    workdir, label, label,
                                                                                    workdir, label, label))

def simulate_something (conf, top, mdp, label, workdir, nprocs=-1):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top)
    mdp = os.path.abspath(mdp)
    workdir = os.path.abspath(workdir)
    runcmd.run("mkdir -p %s/%s" % (workdir, label))
    check_log_loc = "%s/%s/%s.log" % (workdir,label,label)
    simu_state = check_simulation_state(workdir, label)
    if (simu_state == 'COMPLETE_MINIM'):
        # This will only happen for minimization, I guess. We can just assume it finished alright with a warning.
        globalLogger.putMessage('STEP {}: There is a log but no cpt; assumming a successful minimization.'.format(label))
    elif (simu_state == 'FULL'):
        globalLogger.putMessage('STEP {}: Restarting from .cpt (possibly a complete simulation)'.format(label))
        runcmd.run("%s -s %s/%s/%s.tpr -cpi %s/%s/%s.cpt -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label, workdir, label, label, workdir, label, label))
    elif (simu_state == 'NONE'):
        globalLogger.putMessage('STEP {}: Simulating from start'.format(label))
        command = "%s grompp -maxwarn 5 -f %s -c %s -p %s -o %s/%s/%s.tpr" % (ConfigVariables.gmx, mdp, conf, top, workdir, label, label)
        runcmd.run(command)
        command = "%s -s %s/%s/%s.tpr -deffnm %s/%s/%s" % (create_mdrun_call(nprocs), workdir, label, label, workdir, label, label)
        runcmd.run(command)


def get_molecule_name_from_itp (itp):
    itp = os.path.abspath(itp)
    fp = open(itp, "r")
    flag = False
    for line in fp:
        if (re.match(r".*\[ moleculetype \].*", line)):
            flag = True
            continue
        if (flag == True):
            if (re.match(r"^;", line) is None):
                return line.split()[0]
    fp.close()
    return "NOT-FOUND"

def make_topology (nmols, outtop, itp):
    outtop = os.path.abspath(outtop)
    itp = os.path.abspath(itp)
    # create box
    #command = "gmx insert-molecules -ci %s -nmol %d -box %f %f %f -o %s" % (conf,nmols,box[0],box[1],box[2],outconf)
    #runcmd.run(command)
    # make topology
    molecule_name = get_molecule_name_from_itp(itp)
    fp = open(outtop, 'w')
    fp.write("#include \"%s\"\n" % itp)
    fp.write("[ system ]\nSystem\n")
    fp.write("[ molecules ]\n")
    fp.write("%s %d\n" % (molecule_name, nmols))
    fp.close()
    # in itp file, substitute "include" line by appropriate thing
    #fi = open(itp, "r")
    #fp = open(itp + ".tmp", "w")
    #for line in fi:
    #    if (re.match(r".*include.*", line)):
    #        fp.write("#include \"%s/forcefield.itp\"\n" % ffdir)
    #    else:
    #        fp.write(line)
    #fp.close()
    #fi.close()
    #runcmd.run("mv %s.tmp %s" % (itp,itp))

def make_a_box_and_topology (conf, nmols, box, outconf, outtop, itp):
    conf = os.path.abspath(conf)
    outconf = os.path.abspath(outconf)
    outtop = os.path.abspath(outtop)
    itp = os.path.abspath(itp)
    # create box if it does not exist
    if not (os.path.isfile(outconf)):
        command = "%s insert-molecules -ci %s -nmol %d -box %f %f %f -o %s" % (ConfigVariables.gmx, conf,nmols,box[0],box[1],box[2],outconf)
        runcmd.run(command)
    # always make topology if not there
    if True:
    #if not (os.path.isfile(outtop)):
        molecule_name = get_molecule_name_from_itp(itp)
        fp = open(outtop, 'w')
        fp.write("#include \"%s\"\n" % itp)
        fp.write("[ system ]\nLiquid\n")
        fp.write("[ molecules ]\n")
        fp.write("%s %d\n" % (molecule_name, nmols))
        fp.close()

def make_solvation_box_and_topology (confs, nmols, outconf, outtop, itps, makeBox=True):
    # create box if it does not exist
    if not (os.path.isfile(outconf)) and makeBox:
        # first, create a box with one solute molecule
        command = "{} editconf -bt cubic -f {} -d 2.8 -o _box-1.gro".format(ConfigVariables.gmx, confs['solute'])
        runcmd.run(command)
        # complete solutes
        if (nmols['solute'] > 1):
            command = "{} insert-molecules -f _box-1.gro -ci {} -nmol {} -o _box-2.gro".format(ConfigVariables.gmx, confs['solute'], nmols['solute'] - 1)
            runcmd.run(command)
            os.rename("_box-2.gro", "_box-1.gro")
        # complete solvent
        command = "{} insert-molecules -f _box-1.gro -ci {} -nmol {} -o {}".format(ConfigVariables.gmx, confs['solvent'], nmols['solvent'], outconf)
        runcmd.run(command)
        os.remove("_box-1.gro")
    # always make topology if not there
    if True:
    #if not (os.path.isfile(outtop)):
        solute_name = get_molecule_name_from_itp(itps['solute'])
        solvent_name = get_molecule_name_from_itp(itps['solvent'])
        fp = open(outtop, 'w')
        fp.write("#include \"%s\"\n" % itps['solute'])
        fp.write("#include \"%s\"\n" % itps['solvent'])
        fp.write("[ system ]\nSolvation\n")
        fp.write("[ molecules ]\n")
        fp.write("%s %d\n" % (solute_name, nmols['solute']))
        fp.write("%s %d\n" % (solvent_name, nmols['solvent']))
        fp.close()

def simulate_protocol_liquid (conf, nmols, box, itp, mdps, nsteps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    make_a_box_and_topology (conf, nmols, box, workdir + "/liquid.gro", workdir + "/liquid.top", itp)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (workdir + "/liquid.gro", workdir + "/liquid.top", mdps[i], labels[i], workdir)
        elif i == len(mdps) - 1:
            if (read_nsteps_from_mdp(mdps[i]) == nsteps):
                previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                simulate_something (previous_conf, workdir + "/liquid.top", mdps[i], labels[i], workdir)
            else:
                # Extend
                extend_something(nsteps, workdir, labels[i])
        else:
            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
            simulate_something (previous_conf, workdir + "/liquid.top", mdps[i], labels[i], workdir)
    # create output dictionary
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['trr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".trr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = workdir + "/liquid.top"
    output_dict['nsteps'] = nsteps
    #
    return output_dict

def simulate_protocol_solvation (conf, nmols, itp, mdps, nsteps, labels, workdir, simulate=True, follows=False):
    conf = {k: os.path.abspath(v) for k,v in conf.items()}
    itp = {k: os.path.abspath(v) for k,v in itp.items()}
    workdir = os.path.abspath(workdir)
    topopath = workdir + "/solvation.top"
    confpath = workdir + "/solvation.gro"
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    if follows:
        make_solvation_box_and_topology (conf, nmols, confpath, topopath, itp, makeBox=False)
    else:
        make_solvation_box_and_topology (conf, nmols, confpath, topopath, itp, makeBox=simulate)
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

def dummy_protocol_solvation (conf, nmols, itp, mdps, nsteps, labels, workdir):
    simulate_protocol_solvation (conf, nmols, itp, mdps, nsteps, labels, workdir, simulate=False)

def dummy_protocol_liquid (conf, nmols, box, itp, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])

    output_dict = {}

    runcmd.run("mkdir -p " + workdir)
    make_topology (nmols, workdir + "/liquid.top", itp)
#    for i in range(len(mdps)):
#        if i == 0:
#            simulate_something (workdir + "/liquid.gro", workdir + "/liquid.top", mdps[i], labels[i], workdir)
#        else:
#            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
#            simulate_something (previous_conf, workdir + "/liquid.top", mdps[i], labels[i], workdir)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = workdir + "/liquid.top"
    #
    return output_dict

def simulate_protocol_gas (conf, itp, mdps, nsteps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    nmols = 1 # by default
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    make_topology (nmols, workdir + "/gas.top", itp)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
        elif i == len(mdps) - 1:
            if (read_nsteps_from_mdp(mdps[i]) == nsteps):
                previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                simulate_something (previous_conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
            else:
                extend_something (nsteps, workdir, labels[i])
        else:
            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
            simulate_something (previous_conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['trr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".trr"
    output_dict['top'] = workdir + "/gas.top"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['nsteps'] = nsteps
    #
    return output_dict

def dummy_protocol_gas (conf, itp, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    nmols = 1 # by default
    output_dict = {}
    runcmd.run("mkdir -p " + workdir)
    make_topology (nmols, workdir + "/gas.top", itp)
#    for i in range(len(mdps)):
#        if i == 0:
#            simulate_something (conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
#        else:
#            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
#            simulate_something (previous_conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['top'] = workdir + "/gas.top"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    #
    return output_dict

def simulate_protocol_slab (conf, top, liq_tpr, mdps, nsteps, labels, workdir, nprocs):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top)
    liq_tpr = os.path.abspath(liq_tpr)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    # extend initial configuration
    output_dict = {}
    pre_extended_conf = workdir + "/pre.slab.gro"
    extended_conf = workdir + "/slab.gro"
    box = get_box(conf)
    runcmd.run("mkdir -p " + workdir)
    # First, remove periodicity.
    runcmd.run("echo 0 | %s trjconv -f %s -s %s -o %s -pbc whole" % (ConfigVariables.gmx, conf, liq_tpr, pre_extended_conf))
    # Now, extend box.
    runcmd.run("%s editconf -f %s -box %f %f %f -o %s" % (ConfigVariables.gmx, pre_extended_conf, float(box[0]), float(box[1]), 5*float(box[2]), extended_conf))
    # Remove temporary conf.
    os.remove(pre_extended_conf)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (extended_conf, top, mdps[i], labels[i], workdir, nprocs)
        elif i == len(mdps) - 1:
            if (read_nsteps_from_mdp(mdps[i]) == nsteps):
                previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
                simulate_something (previous_conf, top, mdps[i], labels[i], workdir, nprocs)
            else:
                # Extend
                extend_something(nsteps, workdir, labels[i], nprocs)
        else:
            previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
            simulate_something (previous_conf, top, mdps[i], labels[i], workdir, nprocs)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['trr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".trr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = top
    output_dict['nsteps'] = nsteps
    return output_dict

def dummy_protocol_slab (conf, top, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    # extend initial configuration
    output_dict = {}
    #extended_conf = workdir + "/slab.gro"
    #box = get_box(conf)
    runcmd.run("mkdir -p " + workdir)
    #runcmd.run("gmx editconf -f %s -box %f %f %f -o %s" % (conf, float(box[0]), float(box[1]), 5*float(box[2]), extended_conf))
    #for i in range(len(mdps)):
        #if i == 0:
            #simulate_something (extended_conf, top, mdps[i], labels[i], workdir)
        #else:
            #previous_conf = workdir + "/" + labels[i-1] + "/" + labels[i-1] + ".gro"
            #simulate_something (previous_conf, top, mdps[i], labels[i], workdir)
    # create output dictionary 
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = top
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

def dummy_protocol_general(conf, top, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    # create output dictionary
    output_dict['xtc'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".xtc"
    output_dict['tpr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".tpr"
    output_dict['trr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".trr"
    output_dict['edr'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".edr"
    output_dict['gro'] = workdir + "/" + labels[-1] + "/" + labels[-1] + ".gro"
    output_dict['top'] = top
    #
    return output_dict
