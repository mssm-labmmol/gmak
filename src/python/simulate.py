#!/usr/bin/python

import sys
import os
import re
from traj_ana import get_box

# Functions.
#
def simulate_something (conf, top, mdp, label, workdir):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top)
    mdp = os.path.abspath(mdp)
    workdir = os.path.abspath(workdir)
    os.system("mkdir -p %s/%s" % (workdir, label))
    check_log_loc = "%s/%s/%s.log" % (workdir,label,label)
    if (os.path.isfile(check_log_loc)):
        print ("log from simulations " + check_log_loc + " already exists, will not perform it")
    else:
        print ("log from simulations " + check_log_loc + " does not exist, will perform it")
        command = "gmx grompp -maxwarn 5 -f %s -c %s -p %s -o %s/%s/%s.tpr" % (mdp, conf, top, workdir, label, label)
        print ("COMMAND: " + command)
        os.system(command)
        command = "gmx mdrun -s %s/%s/%s.tpr -deffnm %s/%s/%s" % (workdir, label, label, workdir, label, label)
        print ("COMMAND: " + command)
        os.system(command)

def get_molecule_name_from_itp (itp):
    itp = os.path.abspath(itp)
    fp = open(itp, "r")
    flag = False
    for line in fp:
        print ("itp line is " + line)
        if (re.match(r".*\[ moleculetype \].*", line)):
            print ("found moleculetype tag")
            flag = True
            continue
        if (flag == True):
            if (re.match(r"^;", line) is None):
                print ("Matched line is " + line)
                return line.split()[0]
    fp.close()
    return "NOT-FOUND"

def make_topology (nmols, outtop, itp):
    outtop = os.path.abspath(outtop)
    itp = os.path.abspath(itp)
    # create box
    #command = "gmx insert-molecules -ci %s -nmol %d -box %f %f %f -o %s" % (conf,nmols,box[0],box[1],box[2],outconf)
    #os.system(command)
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
    os.system("mv %s.tmp %s" % (itp,itp))

def make_a_box_and_topology (conf, nmols, box, outconf, outtop, itp):
    conf = os.path.abspath(conf)
    outconf = os.path.abspath(outconf)
    outtop = os.path.abspath(outtop)
    itp = os.path.abspath(itp)
    # create box if it does not exist
    if not (os.path.isfile(outconf)):
        command = "gmx insert-molecules -ci %s -nmol %d -box %f %f %f -o %s" % (conf,nmols,box[0],box[1],box[2],outconf)
        os.system(command)
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
        #os.system("mv %s.tmp %s" % (itp,itp))

def simulate_protocol_liquid (conf, nmols, box, itp, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    output_dict = {}
    os.system("mkdir -p " + workdir)
    make_a_box_and_topology (conf, nmols, box, workdir + "/liquid.gro", workdir + "/liquid.top", itp)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (workdir + "/liquid.gro", workdir + "/liquid.top", mdps[i], labels[i], workdir)
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
    #
    return output_dict

def dummy_protocol_liquid (conf, nmols, box, itp, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])

    output_dict = {}

    os.system("mkdir -p " + workdir)
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

def simulate_protocol_gas (conf, itp, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    itp = os.path.abspath(itp)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    nmols = 1 # by default
    output_dict = {}
    os.system("mkdir -p " + workdir)
    make_topology (nmols, workdir + "/gas.top", itp)
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (conf, workdir + "/gas.top", mdps[i], labels[i], workdir)
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
    os.system("mkdir -p " + workdir)
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

def simulate_protocol_slab (conf, top, mdps, labels, workdir):
    conf = os.path.abspath(conf)
    top = os.path.abspath(top)
    workdir = os.path.abspath(workdir)
    for i in range(len(mdps)):
        mdps[i] = os.path.abspath(mdps[i])
    # extend initial configuration
    output_dict = {}
    extended_conf = workdir + "/slab.gro"
    box = get_box(conf)
    os.system("mkdir -p " + workdir)
    os.system("gmx editconf -f %s -box %f %f %f -o %s" % (conf, float(box[0]), float(box[1]), 5*float(box[2]), extended_conf))
    for i in range(len(mdps)):
        if i == 0:
            simulate_something (extended_conf, top, mdps[i], labels[i], workdir)
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
    os.system("mkdir -p " + workdir)
    #os.system("gmx editconf -f %s -box %f %f %f -o %s" % (conf, float(box[0]), float(box[1]), 5*float(box[2]), extended_conf))
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

if __name__ == "__main__":
    labels = ["em", "nvt", "npt", "md"]
    mdps = [ "em_pme.mdp", "nvt_pme.mdp", "npt_pme.mdp", "md_pme.mdp" ]

    initial_configuration = "conf.gro"
    itp = "molecule.itp"
    nmols = 512
    initial_box = [ 3.40, 3.40, 3.40 ]

    ffdir = "./ff"

    # check if labels match mdps
    if (len(labels) != len(mdps)):
        print ("labels do not match mdps")
        exit(0)

    # simulate a full stage
    #simulate_protocol (initial_configuration, nmols, initial_box, ffdir, itp, mdps, labels, "./TEST-DIR-FINAL")
    simulate_protocol_slab ("./liquid.gro", "./liquid.top", mdps, labels, "./coisinha")
