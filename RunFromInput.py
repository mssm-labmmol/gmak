#!/usr/bin/python

import multiprocessing
from gridshifter import * 
from parameters import * 
from simulate import *
from traj_filter import * 
from traj_ana import *
from reweight import * 
from mbar_estimate import *
from grid_ana import * 
import numpy as np
import shlex
import os
import re
import sys

def reweightItpChanger (inputTopology, index, outputTopology):
    fp = open(inputTopology, 'r')
    fp_out = open(outputTopology, 'w')
    for line in fp:
        if (line[:8] == '#include'):
            pathToFile = line[line.find('\"')+1:line.rfind('\"')]
            newPathToFile = pathToFile[:pathToFile.find('.itp')] + ("_%d.itp" % (index))
            newLine = "#include \"%s\"\n" % newPathToFile
            fp_out.write(newLine)
        else:
            fp_out.write(line)
    fp_out.close()
    fp.close()

# This function must have global scope for Pool.map.
def core_parallel_clean_reweight_i_at_j (aux_input_list):
    i             =  aux_input_list[0]
    gp            =  aux_input_list[1]
    j             =  aux_input_list[2]
    gp_other      =  aux_input_list[3]
    protocol      =  aux_input_list[4]
    workdir       =  aux_input_list[5]
    workdir_mbar  =  aux_input_list[6]
    properties    =  aux_input_list[7]

    gp.reweight_with_protocol_at_other (protocol, gp_other,
            "%s/%d_at_%d/" % (workdir, i, j))
    # now get properties and clean
    gi = gp 
    gi_id = i
    gj_id = j 
    # get potential (and) pV
    out_file = workdir_mbar + \
            "/reweighted_properties/potential_%d_%d.xvg" % \
            (gi.id, gj_id)
    xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
    edr = gi.rw_outputs[gj_id][protocol.name]['edr']
    tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
    trr = gi.rw_outputs[gj_id][protocol.name]['trr']
    gro = gi.rw_outputs[gj_id][protocol.name]['gro']
    # If file already exists, no need to create it again.
    if not (os.path.isfile(out_file)):
        print "Getting potential from %s and storing in %s." % \
        (edr, out_file)
        # get
        obtain_property (xtc, edr, gro, tpr, 'potential', out_file)
    if (protocol.has_pv()):
        out_file = workdir_mbar + \
            "/reweighted_properties/pV_%d_%d.xvg" % \
            (gi.id, gj_id)
        # If file already exists, no need to create it again, just set the variable.
        if not (os.path.isfile(out_file)):
            print "Getting pV from %s and storing in %s." % \
            (edr, out_file)
            # get
            obtain_property (xtc, edr, gro, tpr, 'pV', out_file)

    if (j == i):
        # get the properties necessary for estimation
        for i_prop,prop in enumerate(properties):
            out_file = workdir_mbar + \
                    "/filtered_properties/%s_%d.xvg" % \
                    (prop, gi.id)
            if not (os.path.isfile(out_file)):
                print "Getting %s from %d and storing in %s." % \
                (prop, gi.id, out_file)
                # get
                if (prop == 'polcorr'):
                    obtain_polcorr (xtc, edr, gro, tpr, protocol.dipole,\
                            protocol.polar, out_file)
                else:
                    obtain_property (xtc, edr, gro, tpr, prop, out_file)

    # Delete trajectory files (the ones that are big).
    # We don't need to delete the xtc file because it corresponds to the original trajectory.
    for ff in [trr,edr]:
        if (os.path.isfile(ff)):
            os.remove (ff)

# Classes.

class gridOptimizer:
    
    def __init__ (self):
        # dictionary indexed by properties
        self.maxSteps = 5
        self.percentCutoff = 0.25
        self.referenceValues = {}
        self.referenceTolerances = {}
        self.referenceWeights = {}
        # this is a list of TUPLES where each tuple is
        # (gridpoint_id, score)
        # this will be sorted based on score
        self.stateScores = []

    def readFromStream (self, stream):
        for line in stream:
            if (line[0] == '#'):
                continue
            elif (line.rstrip() == '$end'):
                return
            elif (line.split()[0] == 'maxsteps'):
                self.maxSteps = int(line.split()[1])
            elif (line.split()[0] == 'ncut'):
                self.percentCutoff = float(line.split()[1])
            # In all other cases, the line should satisfy the syntax
            # PROPNAME    REFERENCEVALUE     WEIGHT      TOLERANCE
            else:
                splittedLine = line.split()
                propertyName = splittedLine[0]
                self.referenceValues[propertyName] = float(splittedLine[1])
                self.referenceWeights[propertyName] = float(splittedLine[2])
                self.referenceTolerances[propertyName] = float(splittedLine[3])

    def rankScores (self):
        self.stateScores.sort(cmp=lambda x,y: cmp(x[1],y[1]))  

    def fillWithScores (self, grid):
        # clean stateScores first
        self.stateScores = []
        for gP in grid.grid_points:
            gPScore = 0.0
            weightSum = 0.0
            for prop in self.referenceTolerances.keys():
                gPValue = gP.estimated_properties[prop].value
                weightSum += self.referenceWeights[prop]
                gPScore += self.referenceWeights[prop] * ((gPValue - self.referenceValues[prop]) / self.referenceValues[prop]) ** 2
            gPScore /= weightSum
            self.stateScores.append( (gP.id, gPScore) )
        self.rankScores()

    def printToFile (self, grid, filename):
        fp = open(filename, "w")
        properties = self.referenceTolerances.keys()
        fp.write("# %3s" % "id")
        for prop in properties:
            fp.write("%12s" % prop)
            fp.write("%12s" % "err")
        fp.write("%16s\n" % "score")
        for x in self.stateScores:
            fp.write("%5d" % x[0])
            for prop in properties:
                propValue = grid[x[0]].estimated_properties[prop].value
                propErr   = grid[x[0]].estimated_properties[prop].err
                fp.write("%12.4f%12.4f" % (propValue, propErr))
            fp.write("%16.6e\n" % x[1])
        fp.close()

    def determineNextSample (self, grid):
        nTested = 0
        properties = self.referenceTolerances.keys()
        for x in self.stateScores:
            # ignore sampled
            if (x[0] in grid.get_samples_id()):
                nTested += 1
                continue
            if (nTested > self.percentCutoff * grid.linear_size):
                return -1
            for prop in properties:
                propErr = grid[x[0]].estimated_properties[prop].err
                if (propErr > self.referenceTolerances[prop]):
                    return x[0]
            nTested += 1

class PropertyBase:

    def __init__ (self):
        self.name  = "property"
        self.units = "u.a."
        self.symbol = "$P$"
        self.value = 0.0
        self.err   = 0.0

    def set_textual_elements (self, name, units, symbol):
        self.name = name
        self.units = units
        self.symbol = symbol

    def set_value (self, value, err):
        self.value = value
        self.err   = err

    def get_label (self):
        return "%s [%s]" % (self.symbol, self.units)

    def get_label_err (self):
        return "$\\delta$%s [%s]" % (self.symbol, self.units)

class Density (PropertyBase):

    def __init__ (self, value, err):
        self.set_textual_elements ("density", "kg m$^{-3}$", "$\\rho_\\mathrm{liq}$")
        self.value = value
        self.err = err

class Gamma (PropertyBase):

    def __init__ (self, value, err):
        self.set_textual_elements ("gamma", "mN m$^{-1}$", "$\\gamma$")
        self.value = value
        self.err = err

class dHvap (PropertyBase):

    # corrs are taken into account as dHvap_corr = dHvap - value_pol + corrs
    def __init__ (self, value_liq, value_gas, value_pol, err_liq, err_gas, err_pol,\
            corrs, nmols, temperature):
        R = 8.3144598e-3
        self.set_textual_elements ("dhvap", "kJ mol$^{-1}$", "$\\Delta H_\\mathrm{vap}$")
        self.value = value_gas - value_liq/nmols - value_pol + corrs + R*temperature
        self.err   = np.sqrt(err_gas**2 + (err_liq/nmols)**2 + (err_pol)**2)
        print "************************************"
        print "Initialized dHvap = %.2f +/- %.2f " % (self.value, self.err)
        print "Components:"
        print "U_liq = %.2f (%.2f)" % (value_liq, value_liq/nmols)
        print "U_gas = %.2f " % (value_gas)
        print "Polcorr = %.2f " % (value_pol)
        print "Other corrections = %.2f " % (corrs)
        print "RT = %.2f " % (R*temperature)

# GridPoint
#
# Members
# -------
# id: 
#     A number identifying this gridpoint. This is necessary to use gridpoints
#     as keys in dictionaries.
# itp_path:
#     Path of itp file corresponding to this point.
# is_sample:
#     Flag which determines if this grid point is in samples list.
# protocol_outputs:
#     Dictionary containing final output files of the protocols executed.
# atomic_properties:
#     Dictionary containing output files for the properties extracted.
# rw_outputs:
#     A dictionary where keys are other gridpoints and values are
#     dictionaries indexing the output files of the reweighting. For instance,
#     rw_outputs[gp][pr]['xtc'] accesses the xtc file of the reweighting on 
#     grid point 'gp' simulated under protocol 'pr'.
class GridPoint:

    def __init__ (self, itp, idx):
        self.id = idx
        self.itp_path = itp
        self.is_sample = False
        self.protocol_outputs = {}
        self.atomic_properties = {}
        self.rw_outputs = {}
        # also initialize subdictionaries
        self.rw_outputs[self.id] = {}
        # dictionary of estimated_properties
        # e.g. self.estimated_propeties['gamma'] should be a PropertyBase::Gamma
        #     instance estimated from MBAR
        self.estimated_properties = {}

    def add_itp (self, itp):
        self.itp_path = itp

    def prepare_with_protocol_at_dir (self, protocol, workdir):
        protocol.prepare_gridpoint_at_dir (self, workdir)

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        protocol.run_gridpoint_at_dir (self, workdir)

    def reweight_with_protocol_at_other (self, protocol, gp_other, workdir):
        # from reweight.py
        if (protocol.type == 'slab'):
            out_rw = reweight (self.protocol_outputs[protocol.name]['trr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    gp_other.protocol_outputs[protocol.name]['top'],\
                    protocol.mdps[-1], workdir)
        else:
            out_rw = reweight (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    gp_other.protocol_outputs[protocol.name]['top'],\
                    protocol.mdps[-1], workdir)
        self.add_reweight_output (gp_other, protocol, out_rw)

    def set_as_sample (self):
        self.is_sample = True

    def is_sample (self):
        return self.is_sample

    # obj is a PropertyBase::X object (X=Density,Gamma,dHvap)
    def add_property_estimate (self, name, obj):
        if (name != obj.name):
            print "ERROR: expected %s but got %s.\n" % (name, obj.name)
            exit()
        self.estimated_properties[name] = obj

    # This is to ensure that any method that receives a grid point has access
    # to the file names of the simulations performed.
    #
    # protocol_outputs is a dictionary for the output (which are also
    # dictionaries) of a simulation, where keys are the protocol names.
    #
    def add_protocol_output (self, protocol, output):
        self.protocol_outputs[protocol.name] = output

    def add_reweight_output (self, gp_other, protocol, output):
        try:
            self.rw_outputs[gp_other.id][protocol.name] = output
        # key error means something has not been initialized
        except KeyError:
            self.rw_outputs[gp_other.id] = {}
            self.rw_outputs[gp_other.id][protocol.name] = output

    def get_atomic_property_from_protocol (self, name, protocol, output):
        # just to make sure absolute paths are used
        output = os.path.abspath(output)
        # from traj_ana.py
        if (name == 'polcorr'):
            obtain_polcorr (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['edr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    self.protocol_outputs[protocol.name]['tpr'],\
                    protocol.dipole, protocol.polar, output)
        else:
            obtain_property (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['edr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    self.protocol_outputs[protocol.name]['tpr'],\
                    name, output)
        self.add_atomic_property_output (name, output)

    def add_atomic_property_output (self, property_name, output):
        self.atomic_properties[property_name] = output

    def filter_xtc_in_protocol (self, protocol, properties, odir):
        for prop in properties:
            self.get_atomic_property_from_protocol (prop, protocol,\
                    odir + "/" + prop + ".xvg")
        # from traj_filter
        extract_uncorrelated_frames (self.protocol_outputs[protocol.name]['xtc'],\
                self.protocol_outputs[protocol.name]['tpr'], \
                [self.atomic_properties[x] for x in properties],\
                odir + '/filtered_trajectory.xtc')
        # update gridpoint trajectory
        self.protocol_outputs[protocol.name]['xtc'] = \
                os.path.abspath(odir + '/filtered_trajectory.xtc')

    def filter_trr_in_protocol (self, protocol, properties, odir):
        for prop in properties:
            self.get_atomic_property_from_protocol (prop, protocol,\
                    odir + "/" + prop + ".xvg")
        # from traj_filter
        extract_uncorrelated_frames (self.protocol_outputs[protocol.name]['trr'],\
                self.protocol_outputs[protocol.name]['tpr'], \
                [self.atomic_properties[x] for x in properties],\
                odir + '/filtered_trajectory.trr')
        # update gridpoint trajectory
        self.protocol_outputs[protocol.name]['trr'] = \
                os.path.abspath(odir + '/filtered_trajectory.trr')

class ParameterGrid:

    def __init__ (self):
        self.dim = 0
        self.size = []
        self.linear_size = 0;
        self.distmatrix_path = ""
        self.grid_points = []
        # indexed by protocol
        self.hashOfMBAR = {}
        self.xlabel = "$X$"
        self.ylabel = "$Y$"
        return

    def __getitem__ (self, i):
        return self.grid_points[i]

    # This assumes a stream is given.
    #
    # It will read from stream until line with terminating string '$end' is 
    # found.
    #
    def read_from_stream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if (re.match(r"^\$end.*",line)):
                return
            if (line.split()[0] == 'size'):
                # also sets dimension and linear_size
                self.set_size ([int(x) for x in line.split()[1:]])
            if (line.split()[0] == 'samples'):
                self.set_samples ([int(x) for x in line.split()[1:]])
            if (line.split()[0] == 'labels'):
                splitted = shlex.split(line)
                self.xlabel = splitted[1]
                self.ylabel = splitted[2]
            if (line.split()[0] == 'start'):
                # do loop
                prefix = line.split()[1][:-5]
                for i in range(self.linear_size):
                    itp_file = prefix + "%d.itp" % i
                    self.put_itp_at_position (itp_file, i)

    def get_number_of_samples (self):
        num = 0
        for x in self.grid_points:
            if x.is_sample:
                num += 1
        return num

    def get_samples (self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x)
        return sp

    def get_samples_id (self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x.id)
        return sp

    def set_dim (self, dim):
        self.dim = dim
        return

    def set_size (self, size):
        self.dim = len(size)
        self.size = size
        self.linear_size = 1
        for x in size:
            self.linear_size *= x
        for i in range(self.linear_size):
            new_gp = GridPoint("", i)
            self.grid_points.append(new_gp)

    def set_samples (self, samples_list):
        for x in samples_list:
            self.grid_points[x].set_as_sample()

    def add_sample (self, new_sample):
        self.grid_points[new_sample].set_as_sample()

    def put_itp_at_position (self, itp, pos):
        self.grid_points[pos].add_itp(itp)

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        for i,gp in enumerate(self.grid_points):
            # run 
            if gp.is_sample:
                # simulate
                gp.simulate_with_protocol_at_dir (protocol, \
                        workdir + "/" + str(i) + "/")
                # filter
                if (protocol.type == 'slab'):
                    gp.filter_trr_in_protocol (protocol, \
                            protocol.get_filtering_properties(), \
                            workdir + "/" + str(i) + "/")
                else:
                    gp.filter_xtc_in_protocol (protocol, \
                            protocol.get_filtering_properties(), \
                            workdir + "/" + str(i) + "/")
            # only prepare
            else:
                gp.prepare_with_protocol_at_dir (protocol, \
                        workdir + "/" + str(i) + "/")

    def reweight_with_protocol_at_dir (self, protocol, workdir):
        for i,gp in enumerate(self.grid_points):
            # run 
            if gp.is_sample:
                for j,gp_other in enumerate(self.grid_points):
                    gp.reweight_with_protocol_at_other (protocol, gp_other,
                            "%s/%d_at_%d/" % (workdir, i, j))

    def parallel_clean_reweight_with_protocol_at_dir (self, protocol, workdir, workdir_mbar):
        properties = protocol.get_properties()
        ncpus_per_process = 2 
        nprocesses = multiprocessing.cpu_count() / ncpus_per_process

        #arg_list = []
        #count = 0
        #for i,gp in enumerate(self.grid_points):
        #    if gp.is_sample:
        #        for j,gp_other in enumerate(self.grid_points):
        #            pinoffset = count * ncpus_per_process % multiprocessing.cpu_count()
        #            arg_list.append([i,gp,j,gp_other,protocol,workdir,workdir_mbar,properties,\
        #                    ncpus_per_process,pinoffset])
        #            count = count + 1

        #start_idx = 0
        #while ( start_idx < len(arg_list) ):
        #    processes = []
        #    for i in range(nprocesses):
        #        pi = multiprocessing.Process(target=core_parallel_clean_reweight_i_at_j,\
        #                args=arg_list[start_idx + i])
        #        pi.start()
        #    for p in processes:
        #        p.join()
        #    # end of all jobs
        #    start_idx = start_idx + nprocesses

        pool = multiprocessing.Pool(nprocesses)
        pool_arg_list = []
        for i,gp in enumerate(self.grid_points):
            if gp.is_sample:
                for j,gp_other in enumerate(self.grid_points):
                    pool_arg_list.append([i,gp,j,gp_other,protocol,workdir,workdir_mbar,properties])
        pool.map(core_parallel_clean_reweight_i_at_j, pool_arg_list)

    # 
    def fast_clean_reweight_with_protocol_at_dir (self, protocol, reweightHash, workdir, workdir_mbar):
        properties = protocol.get_properties()
        for i,gp in enumerate(self.grid_points):
            if (gp.is_sample):
                # Clean this list.
                preReweightOutputs = []
                listOfPotentialFiles = []
                # Before reweighting, rerun the trajectory progressively
                # deactivating the energy components associated to each parameter.
                numberOfRuns = 2 * int( reweightHash['npars']) + 1
                for j in range(numberOfRuns):
                    # Modify the topology of the original simulation to refer to
                    # the proper *.itp file. Only 'j' is needed, as the names of
                    # the files are hard-coded.
                    os.system("mkdir -p %s/decouple_%d_at_%d" % (workdir, i, j))
                    topologyName = "%s/decouple_%d_at_%d/topol.top" % (workdir, i, j)
                    reweightItpChanger (gp.protocol_outputs[protocol.name]['top'], j,\
                            topologyName)
                    # Run the reweighting.
                    if (protocol.type == 'slab'):
                        preReweightOutput = reweight (gp.protocol_outputs[protocol.name]['trr'],\
                                gp.protocol_outputs[protocol.name]['gro'],\
                                topologyName,\
                                protocol.mdps[-1],\
                                "%s/decouple_%d_at_%d" % (workdir, i, j))
                    else:
                        preReweightOutput = reweight (gp.protocol_outputs[protocol.name]['xtc'],\
                                gp.protocol_outputs[protocol.name]['gro'],\
                                topologyName,\
                                protocol.mdps[-1],\
                                "%s/decouple_%d_at_%d" % (workdir, i, j))
                    # Save.
                    xtc = preReweightOutput['xtc']
                    edr = preReweightOutput['edr']
                    gro = preReweightOutput['gro']
                    tpr = preReweightOutput['tpr']
                    outPotential =  "%s/decouple_%d_at_%d/potential.xvg" % (workdir, i, j)
                    preReweightOutputs.append(preReweightOutput)
                    listOfPotentialFiles.append(outPotential)
                    # Get the potential energy of the reweighted trajectory.
                    obtain_property (xtc, edr, gro, tpr, 'potential', outPotential)

                # Run the reweighing program to obtain the reweighted potential
                # energies.
                filesArg = " ".join(listOfPotentialFiles)
                fileGrid = reweightHash['parameters']
                anaRwBinary = reweightHash['program']
                command  = "%s -f %s -p %s -r %d -o %s/reweighted_properties/potential_%d_%%d.xvg -time" % (anaRwBinary, filesArg, fileGrid, i, workdir_mbar, i)
                print "ANA_RW: %s" % command 
                # Create containing directory.
                os.system("mkdir -p %s/reweighted_properties/" % workdir_mbar)
                os.system(command)

                # Redefine variables.
                # Note that rw = 0 corresponds to reweighting in the original state.
                gi  = gp
                xtc = preReweightOutputs[0]['xtc']
                edr = preReweightOutputs[0]['edr']
                tpr = preReweightOutputs[0]['tpr']
                gro = preReweightOutputs[0]['gro']
                # Calculate pV if necessary -- it is the same for all reweighted states.
                if (protocol.has_pv()):
                    for j,gp_other in enumerate(self.grid_points):
                        out_file = workdir_mbar + \
                                "/reweighted_properties/pV_%d_%d.xvg" % \
                                (gi.id, j)
                        # Do nothing if file exists.
                        if not (os.path.isfile(out_file)):
                            obtain_property (xtc, edr, gro, tpr, "pV", out_file)
                # Calculate filtered properties to use in MBAR estimation.
                for i_prop,prop in enumerate(properties):
                    out_file = workdir_mbar + \
                            "/filtered_properties/%s_%d.xvg" % \
                            (prop, gi.id)
                    if not (os.path.isfile(out_file)):
                        if (prop == 'polcorr'):
                            obtain_polcorr (xtc, edr, gro, tpr, protocol.dipole,\
                                    protocol.polar, out_file)
                        else:
                            obtain_property (xtc, edr, gro, tpr, prop, out_file)

    def clean_reweight_with_protocol_at_dir (self, protocol, workdir, workdir_mbar):
        properties = protocol.get_properties()
        for i,gp in enumerate(self.grid_points):
            # run 
            if gp.is_sample:
                for j,gp_other in enumerate(self.grid_points):
                    gp.reweight_with_protocol_at_other (protocol, gp_other,
                            "%s/%d_at_%d/" % (workdir, i, j))
                    # now get properties and clean
                    gi = gp 
                    gi_id = i
                    gj_id = j 
                    # get potential (and) pV
                    out_file = workdir_mbar + \
                            "/reweighted_properties/potential_%d_%d.xvg" % \
                            (gi.id, gj_id)
                    xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
                    edr = gi.rw_outputs[gj_id][protocol.name]['edr']
                    tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
                    trr = gi.rw_outputs[gj_id][protocol.name]['trr']
                    gro = gi.rw_outputs[gj_id][protocol.name]['gro']
                    # If file already exists, no need to create it again.
                    if not (os.path.isfile(out_file)):
                        print "Getting potential from %s and storing in %s." % \
                        (edr, out_file)
                        # get
                        obtain_property (xtc, edr, gro, tpr, 'potential', out_file)
                    if (protocol.has_pv()):
                        out_file = workdir_mbar + \
                            "/reweighted_properties/pV_%d_%d.xvg" % \
                            (gi.id, gj_id)
                        # If file already exists, no need to create it again, just set the variable.
                        if not (os.path.isfile(out_file)):
                            print "Getting pV from %s and storing in %s." % \
                            (edr, out_file)
                            # get
                            obtain_property (xtc, edr, gro, tpr, 'pV', out_file)

                    if (j == i):
                        # get the properties necessary for estimation
                        for i_prop,prop in enumerate(properties):
                            out_file = workdir_mbar + \
                                    "/filtered_properties/%s_%d.xvg" % \
                                    (prop, gi.id)
                            if not (os.path.isfile(out_file)):
                                print "Getting %s from %d and storing in %s." % \
                                (prop, gi.id, out_file)
                                # get
                                if (prop == 'polcorr'):
                                    obtain_polcorr (xtc, edr, gro, tpr, protocol.dipole,\
                                            protocol.polar, out_file)
                                else:
                                    obtain_property (xtc, edr, gro, tpr, prop, out_file)

                    # Delete trajectory files (the ones that are big).
                    # We don't need to delete the xtc file because it corresponds to the original trajectory.
                    for ff in [trr,edr]:
                        if (os.path.isfile(ff)):
                            os.remove (ff)


    def make_grid_for_protocol (self, protocol, workdir, reweightHash):
        simu_dir = workdir + "/simu"
        rw_dir   = workdir + "/rw"
        mbar_dir = workdir + "/mbar"
        mbar_out = []
        # simulate
        self.simulate_with_protocol_at_dir (protocol, simu_dir)
        # reweight 
        #self.reweight_with_protocol_at_dir (protocol, rw_dir)
        # reweight, obtain properties and clean
        # mbar_dir must be given to put the properties in the correct place

        #self.clean_reweight_with_protocol_at_dir (protocol, rw_dir, mbar_dir)

        if (reweightHash['type'] == 'standard'):
            self.parallel_clean_reweight_with_protocol_at_dir (protocol, rw_dir, mbar_dir)
        elif (reweightHash['type'] == 'fast'):
            # NOTE More arguments are probably needed in this function.
            self.fast_clean_reweight_with_protocol_at_dir (protocol, reweightHash, rw_dir, mbar_dir)
        else:
            print "|\n|\n|\n|\n|\n|ERROR: Reweight type must be 'standard' or 'fast'."
            return

        # estimate 
        mbar = MBARControl (self, protocol, reweightHash, mbar_dir)
        mbar.estimate ()
        self.hashOfMBAR[protocol.name] = mbar

    def compute_final_properties (self, prop, propProtocols):
        for gp in self.grid_points:
            if (prop == 'density'):
                propIdxInMBAR = self.hashOfMBAR[propProtocols].protocol.properties.index('density')
                estimateValue = self.hashOfMBAR[propProtocols].EA_k[propIdxInMBAR,gp.id]
                estimateErr   = self.hashOfMBAR[propProtocols].dEA_k[propIdxInMBAR,gp.id]
                estimateObj   = Density (estimateValue, estimateErr)
                gp.add_property_estimate (prop, estimateObj)
            if (prop == 'gamma'):
                propIdxInMBAR = self.hashOfMBAR[propProtocols].protocol.properties.index('gamma')
                estimateValue = self.hashOfMBAR[propProtocols].EA_k[propIdxInMBAR,gp.id]
                estimateErr   = self.hashOfMBAR[propProtocols].dEA_k[propIdxInMBAR,gp.id]
                estimateObj   = Gamma (estimateValue, estimateErr)
                gp.add_property_estimate (prop, estimateObj)
            if (prop == 'dhvap'):
                print "***************************\n\n\n\nhash of mbar is %s \n\n\n\n%s\n%s\n%s\n%s" % (self.hashOfMBAR,\
                     self.hashOfMBAR[propProtocols[0]].EA_k,  self.hashOfMBAR[propProtocols[1]].EA_k, self.hashOfMBAR[propProtocols[0]].protocol.properties,\
                      self.hashOfMBAR[propProtocols[1]].protocol.properties)
                propIdxInMBAR = self.hashOfMBAR[propProtocols[0]].protocol.properties.index('potential')
                estimateValueLiq = self.hashOfMBAR[propProtocols[0]].EA_k[propIdxInMBAR,gp.id]
                estimateErrLiq   = self.hashOfMBAR[propProtocols[0]].dEA_k[propIdxInMBAR,gp.id]
                nmols            = self.hashOfMBAR[propProtocols[0]].protocol.nmols
                #
                propIdxInMBAR = self.hashOfMBAR[propProtocols[1]].protocol.properties.index('potential')
                estimateValueGas = self.hashOfMBAR[propProtocols[1]].EA_k[propIdxInMBAR,gp.id]
                estimateErrGas   = self.hashOfMBAR[propProtocols[1]].dEA_k[propIdxInMBAR,gp.id]
                #
                propIdxInMBAR = self.hashOfMBAR[propProtocols[1]].protocol.properties.index('polcorr')
                estimateValuePol = self.hashOfMBAR[propProtocols[1]].EA_k[propIdxInMBAR,gp.id]
                estimateErrPol   = self.hashOfMBAR[propProtocols[1]].dEA_k[propIdxInMBAR,gp.id]
                corr             = self.hashOfMBAR[propProtocols[1]].protocol.other_corrections
                temp             = self.hashOfMBAR[propProtocols[1]].protocol.get_temperature()
                #
                estimateObj   = dHvap (estimateValueLiq, estimateValueGas, estimateValuePol,\
                        estimateErrLiq,estimateErrGas,estimateErrPol, corr, nmols, temp)
    #def __init__ (self, value_liq, value_gas, value_pol, err_liq, err_gas, err_pol,\
   #         corrs, nmols, temperature):
                gp.add_property_estimate (prop, estimateObj)

    def save_property_values_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def save_property_err_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def save_property_diff_to_file (self, prop, ref, filename):
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def plot_property_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.size[0], self.size[1])
        cbox_label = ""
        cbox_limits_colors = ()
        cbox_limits = ()
        title = self.grid_points[0].estimated_properties[prop].get_label()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())

    def plot_property_err_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.size[0], self.size[1])
        cbox_label = ""
        cbox_limits_colors = ('white', 'blue')
        cbox_limits = (0.0, np.max(data))
        title = self.grid_points[0].estimated_properties[prop].get_label_err()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())

    def plot_property_diff_to_file (self, prop, ref, filename):
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.size[0], self.size[1])
        cbox_label = ""
        # hard-coded preferences
        cbox_limits_colors = ('red', 'white', 'blue')
        # hard-coded preferences
        limits = {\
                'density': 20,\
                'dhvap': 4,\
                'gamma': 10}
        cbox_limits = (-limits[prop], limits[prop])
        title = "$\\Delta^\\mathrm{exp}$" + self.grid_points[0].estimated_properties[prop].get_label()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())

    def clean_reweight_residues (self, protocol):
        for gi in self.grid_points:
            if gi.is_sample:
                for gj_id in gi.rw_outputs:
                    gi_id = self.get_samples().index(gi)
                    xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
                    trr = gi.rw_outputs[gj_id][protocol.name]['trr']
                    edr = gi.rw_outputs[gj_id][protocol.name]['edr']
                    tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
                    gro = gi.rw_outputs[gj_id][protocol.name]['gro']
                    for res_file in [xtc,trr,edr,tpr]:
                        if (os.path.isfile(res_file)):
                            os.remove(res_file)

class LiquidProtocol:

    def __init__ (self):
        self.name = ""
        self.type = "liquid"
        self.mdps = []
        self.mdps_rw = []
        self.properties = []
        self.coords = ""
        self.box_size = 0.0

    def __init__ (self, name, nmols, coords, mdps, box, properties):
        self.name = name
        self.type = "liquid"
        self.mdps = mdps
        self.coords = coords
        self.box_size = box
        self.nmols = nmols
        self.properties = properties

    def read_from_stream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                return
            if line.split()[0] == 'name':
                self.name = line.split()[1]
            if line.split()[0] == 'mdps':
                self.mdps = line.split()[1:]
            if line.split()[0] == 'coords':
                self.coords = line.split()[1]
            if line.split()[0] == 'length':
                box_len = float(line.split()[1])
                self.box_size = [box_len, box_len, box_len]
            if line.split()[0] == 'nmols':
                self.nmols = int(line.split()[1])

    def has_pv (self):
        return True

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        out_liquid = simulate_protocol_liquid (self.coords, self.nmols,\
                self.box_size, gridpoint.itp_path, self.mdps, labels, workdir)
        gridpoint.add_protocol_output (self, out_liquid)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        out_liquid = dummy_protocol_liquid (self.coords, self.nmols, \
                self.box_size, gridpoint.itp_path, self.mdps, labels, workdir)
        gridpoint.add_protocol_output (self, out_liquid)
        return
    
    def get_temperature (self):
        temp = 0.0
        fp = open(self.mdps[-1], "r")
        for line in fp:
            if ((re.match(r"^ref_t.*", line)) or (re.match(r"^ref-t.*", line))):
                temp = float(line.split()[2])
        fp.close()
        return temp

    def get_filtering_properties (self):
        return ['potential']

    def set_properties (self, properties):
        self.properties = properties

    def get_properties (self):
        return self.properties

class GasProtocol:

    def __init__ (self):
        self.name = ""
        self.type = "gas"
        self.mdps = []
        self.coords = ""
        self.box_size = 0.0

    def __init__ (self, name, coords, dipole, polar, mdps, properties):
        self.name = name
        self.type = "gas"
        self.mdps = mdps
        self.coords = coords
        self.nmols = 1
        self.properties = properties
        self.dipole = dipole
        self.polar = polar

    def read_from_stream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                return
            if line.split()[0] == 'name':
                self.name = line.split()[1]
            if line.split()[0] == 'mdps':
                self.mdps = line.split()[1:]
            if line.split()[0] == 'coords':
                self.coords = line.split()[1]
            if line.split()[0] == 'gasdipole':
                self.dipole = float(line.split()[1])
            if line.split()[0] == 'polarizability':
                self.polar = float(line.split()[1])

    def has_pv (self):
        # NOTE: Formally, the gas simulations do have a pV term. However, they
        # are the same in all cases (pV = NkT), so they have no influence over
        # the distribution in phase space.
        return False

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        out_gas = simulate_protocol_gas (self.coords, gridpoint.itp_path,\
                self.mdps, labels, workdir)
        gridpoint.add_protocol_output (self, out_gas)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        out_gas = dummy_protocol_gas (self.coords, gridpoint.itp_path,\
                self.mdps, labels, workdir)
        gridpoint.add_protocol_output (self, out_gas)
        return

    def get_filtering_properties (self):
        return ['potential']

    def get_properties (self):
        return self.properties

    def get_temperature (self):
        temp = 0.0
        fp = open(self.mdps[-1], "r")
        for line in fp:
            if ((re.match(r"^ref_t.*", line)) or (re.match(r"^ref-t.*", line))):
                temp = float(line.split()[2])
        fp.close()
        return temp

    def set_other_corrections (self, corr):
        self.other_corrections = corr

class SlabProtocol:

    def __init__ (self):
        self.name = ""
        self.type = "slab"
        self.mdps = []
        self.box_size = 0.0

    def __init__ (self, name, mdps, factor, properties):
        self.name = name
        self.type = "slab"
        self.mdps = mdps
        self.properties = properties

    def read_from_stream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                return
            if line.split()[0] == 'name':
                self.name = line.split()[1]
            if line.split()[0] == 'mdps':
                self.mdps = line.split()[1:]
            if line.split()[0] == 'follow':
                self.follow = line.split()[1]

    def set_follow (self, pro_name):
        self.follow = pro_name
        
    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels   =  [str(x) for x in range(len(self.mdps))]
        conf     =  gridpoint.protocol_outputs[self.follow]['gro']
        top      =  gridpoint.protocol_outputs[self.follow]['top']
        out_slab = simulate_protocol_slab (conf, top, self.mdps, labels,\
                workdir)
        gridpoint.add_protocol_output (self, out_slab)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        labels   =  [str(x) for x in range(len(self.mdps))]
        conf     =  gridpoint.protocol_outputs[self.follow]['gro']
        top      =  gridpoint.protocol_outputs[self.follow]['top']
        out_slab = dummy_protocol_slab (conf, top, self.mdps, labels,\
                workdir)
        gridpoint.add_protocol_output (self, out_slab)
        return

    def has_pv (self):
        return False

    def get_properties (self):
        return self.properties

    def get_temperature (self):
        temp = 0.0
        fp = open(self.mdps[-1], "r")
        for line in fp:
            if ((re.match(r"^ref_t.*", line)) or (re.match(r"^ref-t.*", line))):
                temp = float(line.split()[2])
        fp.close()
        return temp

    def get_filtering_properties (self):
        return ['potential']


# **************************************************************************** #
# *                                   MBAR                                   * #
# **************************************************************************** #

class MBARControl:

    # Parameter grid is given to facilitate setting dimensions.
    # Protocol is given to retrieve temperature and properties.
    def __init__ (self, parameter_grid, protocol, reweightHash, workdir):
        self.n_samples = parameter_grid.get_number_of_samples()
        self.u_kn = [["" for j in range(self.n_samples)] for i in \
                range(parameter_grid.linear_size)]
        self.pv_kn = [["" for j in range(self.n_samples)] for i in \
                range(parameter_grid.linear_size)]
        self.temp = protocol.get_temperature()
        self.workdir = os.path.abspath(workdir)
        self.prop_matrix = [["" for i in range(self.n_samples)]\
                for j in range(len(protocol.get_properties()))]
        self.properties = protocol.get_properties()

        # Links to grid and protocol.
        self.parameter_grid = parameter_grid
        self.protocol = protocol

        # Estimated values and uncertainties.
        self.EA_k  = []
        self.dEA_k  = []

        # Head message.
        print \
"""
****************************************************************************
* MBAR control initialization                                              *
****************************************************************************
MBAR was initialized with %d samples reweighted over %d states.
The temperature read from protocol \"%s\" was %f.
The properties estimated for this protocol are %s.
""" % (self.n_samples, parameter_grid.linear_size, protocol.name, self.temp, \
        protocol.get_properties())

        # Fill u_kn and pv_kn with file names.
        for gi in parameter_grid.grid_points:
            if gi.is_sample:
                for gj_id,gj in enumerate(parameter_grid.grid_points):
                    gi_id = parameter_grid.get_samples().index(gi)
                    out_file = self.workdir + \
                            "/reweighted_properties/potential_%d_%d.xvg" % \
                            (gi.id, gj_id)
                    # If file already exists, no need to create it again.
                    if (os.path.isfile(out_file)):
                        self.u_kn[gj_id][gi_id] = out_file
                    else:
                        if reweightHash['type'] == 'fast':
                            print "UNEXPECTED BEHAVIOR: File %s does not exist." % out_file
                        xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
                        edr = gi.rw_outputs[gj_id][protocol.name]['edr']
                        tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
                        gro = gi.rw_outputs[gj_id][protocol.name]['gro']
                        print "Getting potential from %s and storing in %s." % \
                        (edr, out_file)
                        # get
                        obtain_property (xtc, edr, gro, tpr, 'potential', out_file)
                        # fill u_kn
                        self.u_kn[gj_id][gi_id] = out_file
                    if (protocol.has_pv()):
                        out_file = self.workdir + \
                            "/reweighted_properties/pV_%d_%d.xvg" % \
                            (gi.id, gj_id)
                        # If file already exists, no need to create it again, just set the variable.
                        if (os.path.isfile(out_file)):
                            self.pv_kn[gj_id][gi_id] = out_file
                        else:
                            print "Getting pV from %s and storing in %s." % \
                            (edr, out_file)
                            # get
                            obtain_property (xtc, edr, gro, tpr, 'pV', out_file)
                            # fill pv_kn
                            self.pv_kn[gj_id][gi_id] = out_file

    def estimate (self):
        # Create property files.
        for i_prop,prop in enumerate(self.properties):
            for gi in self.parameter_grid.grid_points:
                if gi.is_sample:
                    gi_id = self.parameter_grid.get_samples().index(gi)
                    out_file = self.workdir + \
                            "/filtered_properties/%s_%d.xvg" % \
                            (prop, gi.id)
                    if not (os.path.isfile(out_file)):
                        xtc = gi.rw_outputs[gi.id][self.protocol.name]['xtc']
                        edr = gi.rw_outputs[gi.id][self.protocol.name]['edr']
                        tpr = gi.rw_outputs[gi.id][self.protocol.name]['tpr']
                        gro = gi.rw_outputs[gi.id][self.protocol.name]['gro']
                        print "Getting %s from %d and storing in %s." % \
                        (prop, gi.id, out_file)
                        # get
                        if (prop == 'polcorr'):
                            obtain_polcorr (xtc, edr, gro, tpr, self.protocol.dipole,\
                                    self.protocol.polar, out_file)
                        else:
                            obtain_property (xtc, edr, gro, tpr, prop, out_file)
                        # fill matrix
                        print "Filling position %d,%d for array of shape %s" % (i_prop, gi_id, np.array(self.prop_matrix).shape)
                        self.prop_matrix[i_prop][gi_id] = out_file
                    else:
                        self.prop_matrix[i_prop][gi_id] = out_file

        # Actually estimate.
        # There are two cases: with and without pV.
        ukn_out = self.workdir + "/details/u_kn.dat"
        nk_out = self.workdir + "/details/n_k.dat"
        eff_out = self.workdir + "/details/eff_k.dat"
        eig_out = self.workdir + "/details/eig_k.dat"
        mat_out = self.workdir + "/details/mat_k.dat"
        out_preffixes = [self.workdir + "/estimated_properties/%s" % prop \
                for prop in self.properties]
        if (self.protocol.has_pv()):
            out = \
                    estimate_properties (self.u_kn, self.pv_kn, \
                    self.parameter_grid.get_samples_id(),\
                    self.temp, nk_out, ukn_out, eff_out, eig_out, mat_out,\
                    self.prop_matrix, out_preffixes)
            self.EA_k = out[:,0]
            self.dEA_k = out[:,1]
        else:
            out = \
            estimate_properties_no_pv (self.u_kn,\
                    self.parameter_grid.get_samples_id(),\
                    self.temp, nk_out, ukn_out, eff_out, eig_out, mat_out,\
                    self.prop_matrix, out_preffixes)
            self.EA_k = out[:,0]
            self.dEA_k = out[:,1]

def initialize_from_input (input_file):
    # will return a tuple containing, in order:
    # 0 = the grid
    # 1 = the list of protocols
    # 2 = specifications of grid optimization
    # 3 = specifications of the reweighting procedure (a hash)
    output_grid = ParameterGrid ()
    output_protocols = []
    output_properties = []
    output_protocolsHash = {}
    output_workdir = ""
    output_optimizer = gridOptimizer ()
    output_gridshifter = GridShifter ()
    output_paramLoop = ParameterLoop ()
    output_reweightHash = {}
    fp = open(input_file, "r")
    for line in fp:
        if line[0] == '#':
            continue
        # grid
        if len(line.split()) < 1:
            continue
        if (line.split()[0] == "workdir"):
            output_workdir = os.path.abspath(line.split()[1])
        if (line.rstrip() == "$grid"):
            output_grid.read_from_stream (fp)
        if (line.rstrip() == "$protocol"):
            line = next(fp)
            if (line.split()[0] == 'type'):
                typeRead = line.split()[1]
                if (typeRead == 'liquid'):
                    new_protocol = LiquidProtocol("",0,"",[""],[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'gas'):
                    new_protocol = GasProtocol("","",0.0,0.0,[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'slab'):
                    new_protocol = SlabProtocol("",[],5.0,[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                else:
                    print \
"""ERROR: Type \"%s\" is not supported.\n""" % typeRead
                    exit()
            else:
                print \
"""ERROR: First line after a $protocol flag MUST assign its type.\n"""
                exit()
        if (line.rstrip() == '$compute'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                propRead = line.split()[0]
                nameRead = line.split()[1]
                output_properties.append(propRead)
                if (propRead == 'density'):
                    output_protocolsHash[propRead] = nameRead
                    # find protocol with name given
                    protocols = filter (lambda x: x.name == nameRead, \
                            output_protocols)
                    for protocol in protocols:
                        protocol.properties.append('density') 
                        # potential ALWAYS
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                elif (propRead == 'dhvap'):
                    nameLiq = line.split()[1]
                    nameGas = line.split()[2]
                    corr    = float(line.split()[3])

                    output_protocolsHash[propRead] = [nameLiq,nameGas]
                    # find protocol with name given - Liq
                    protocols = filter (lambda x: x.name == nameLiq, \
                            output_protocols)
                    for protocol in protocols:
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')

                    # find protocol with name given - Gas
                    protocols = filter (lambda x: x.name == nameGas, \
                            output_protocols)
                    for protocol in protocols:
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                        protocol.properties.append('polcorr')
                        protocol.set_other_corrections(corr)

                elif (propRead == 'gamma'):
                    # find protocol with name given
                    output_protocolsHash[propRead] = nameRead
                    protocols = filter (lambda x: x.name == nameRead, \
                            output_protocols)
                    for protocol in protocols:
                        protocol.properties.append('gamma')
                        # potential ALWAYS
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                else:
                    print \
"""ERROR: Property \"%s\" is not supported.\n""" % typeRead
                    exit()
        if (line.rstrip() == '$optimize'):
            output_optimizer.readFromStream (fp)
        if (line.rstrip() == '$parameters'):
            output_paramLoop.readFromStream (fp)
        if (line.rstrip() == '$gridshift'):
            output_gridshifter.readFromStream (fp)
        if (line.rstrip() == '$reweight'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                option = line.split()[0]
                value  = line.split()[1]
                output_reweightHash[option] = value
    fp.close()
    return (output_workdir, output_grid, output_protocols, output_properties, \
            output_protocolsHash, output_optimizer, output_gridshifter,\
            output_paramLoop, output_reweightHash)

# **************************************************************************** #
# *                                   MAIN                                   * #
# **************************************************************************** #
#
if __name__ == "__main__":
    
    (base_workdir, grid, protocols, properties, protocolsHash, optimizer, gridShifter, paramLoop, reweightHash) = \
            initialize_from_input (sys.argv[1])

    for nshifts in range(gridShifter.maxshifts):
        nGrid = nshifts + 1
        nextSample = -1
        for nsteps in range(optimizer.maxSteps):
            if (nsteps != 0):
                grid.add_sample(nextSample)
            workdir = base_workdir + ("/grid_%d/" % nGrid)
            thisRunOutputs = workdir + "/step_%d" % (nsteps+1)
            os.system("mkdir -p %s" % thisRunOutputs)

            # create the grid file
            gridFile = workdir + "/grid.dat"
            paramLoop.loop.write_parameter_grid_to_file(gridFile)
            reweightHash['parameters'] = gridFile

            # create itp files for each grid point
            os.system("mkdir -p " + workdir + "/topo")
            for gridpoint in grid.grid_points:
                itpPathPreffix = workdir + "/topo/et_%d" % (gridpoint.id)
                paramLoop.create_full_itp_for_position (gridpoint.id, itpPathPreffix)
                gridpoint.itp_path = itpPathPreffix + ".itp"

            for protocol in protocols:
                grid.make_grid_for_protocol (protocol, workdir + "/" + protocol.name, reweightHash)
                # clean!
                #grid.clean_reweight_residues(protocol)
            #
            for prop in properties:
                referenceValue = optimizer.referenceValues[prop]
                # Parameters: property name and a list of protocols used for this
                # property.
                grid.compute_final_properties(prop, protocolsHash[prop])
                grid.save_property_values_to_file (prop, thisRunOutputs + '/' + prop + '_EA_k.dat')
                grid.save_property_err_to_file (prop, thisRunOutputs + '/' + prop + '_dEA_k.dat')
                grid.save_property_diff_to_file (prop, referenceValue, thisRunOutputs + '/' + prop + '_diff.dat')
                #
                grid.plot_property_to_file (prop, thisRunOutputs + "/" + prop + "_EA_k.pdf")
                grid.plot_property_err_to_file (prop, thisRunOutputs + "/" + prop + "_dEA_k.pdf")
                grid.plot_property_diff_to_file (prop, referenceValue, thisRunOutputs + "/" + prop + "_diff.pdf")
            #
            optimizer.fillWithScores (grid)
            optimizer.printToFile (grid, thisRunOutputs + "/optimizer_data.dat")
            nextSample = optimizer.determineNextSample (grid)
            print "Next sample is %d"  % nextSample
            if (nextSample == -1):
                break
        new_workdir = base_workdir + "/grid_%d" % (nGrid+1)
        if not (gridShifter.apply (optimizer, paramLoop, grid, workdir, new_workdir)):
            print ("End of grid shift!")
            break
                
    ## set things
    ##properties = ['density', 'potential']
    #mdps_liq = ["em_pme.mdp", "nvt_pme.mdp", "npt_pme.mdp", "md_pme.mdp"]
    ##properties = ['potential', 'polcorr']
    #properties_liq = ['density', 'potential']
    ##mdps = ['em_gas.mdp', 'sd_gas.mdp']
    ## also dipole and polarizability
    ##gas = GasProtocol ("my-gas", "conf.gro", 1.85, 0.00147, mdps, properties)
    #liq = LiquidProtocol ("my-liq", 512, "conf.gro",  mdps_liq, [3.4]*3, properties_liq)

    #grid = ParameterGrid()
    #grid.set_size([4])
    #grid.set_samples([1,2])
    #for i in range(4):
    #    grid[i].add_itp("water_itp_" + str(i) + ".itp")

    #grid.make_grid_for_protocol (liq, "compact_run/")

    ## do ~~EPIC~~ stuff
    #grid.simulate_with_protocol_at_dir (liq, "test_liq_RunFromInput/EPIC_simu/")
    #grid.reweight_with_protocol_at_dir (liq, "test_liq_RunFromInput/EPIC_rw")

    ## MBAR 
    #mbar_q = MBARControl (grid, liq, "test_liq_RunFromInput/EPIC_mbar")
    #mbar_q.estimate ()

    ## test slab
    #slab = SlabProtocol ("my-slab", ['em_slab.mdp','nvt_slab_eq.mdp','nvt_slab.mdp'],\
    #        5.0, ['gamma', 'potential'])
    #slab.set_follow("my-liq")
    #grid.simulate_with_protocol_at_dir (slab, "test_slab_RunFromInput/EPIC_simu/")
    #grid.reweight_with_protocol_at_dir (slab, "test_slab_RunFromInput/EPIC_rw")
    #mbar_s = MBARControl (grid, slab, "test_slab_RunFromInput/EPIC_mbar")
    #mbar_s.estimate ()
#
