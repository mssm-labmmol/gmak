#!/usr/bin/python

from gridshifter import * 
from parameters import * 
from grid_ana import * 
from coverage import *
from read_input import *
from mdputils import mdpUtils
from protocols import LiquidProtocol, GasProtocol, SlabProtocol
from gridbase import GridPoint, ParameterGrid
from gridoptimizer import gridOptimizer

import numpy as np
import os
import sys

if ('--legacy' in sys.argv):
    bool_legacy = True
else:
    bool_legacy = False

if __name__ == "__main__":
    
    (base_workdir,
     grid,
     protocols,
     properties,
     protocolsHash,
     optimizer,
     gridShifter,
     paramLoop,
     gaCoverInterface,
     surrogateModelHash,
     subgridHash,
    reweightHash) = initialize_from_input (sys.argv[1], bool_legacy)


    # *********************** Initial checks for run **************************

    # Consistency checks in input file.
    if (tuple(grid.size) != (paramLoop.loop.get_size())):
        raise ValueError("Sizes of grid and parameter loop are incompatible.")
    if (reweightHash['npars'] != grid.dim):
        raise ValueError("Expected dimension {} in fast reweighting, but got {}.".format(
            reweightHash['npars'], grid.dim))

    # In each protocol, check writing frequencies of energy and trajectory to avoid any problems.
    mdp_ut = mdpUtils()
    for protocol in protocols:
        # Check only production run.
        mdp_ut.parse_file(protocol.mdps[-1])
        efreq = mdp_ut.get_writefreq_energy()
        if (protocol.type == 'slab'):
            xfreq = mdp_ut.get_writefreq_pos()
        else:
            xfreq = mdp_ut.get_writefreq_compressedpos()
        if (efreq != xfreq):
            raise ValueError("For protocol {}, efreq ({}) and xfreq ({}) don't match.".format(
                protocol.name, efreq, xfreq))

    # *********************** End of checks for run  **************************        
    
    for nshifts in range(gridShifter.maxshifts):
        nGrid = nshifts + 1
        nextSample = -1

        # If necessary, fill corners before filling interior of grid.
        for protocol in protocols:
            if protocol.requires_corners():
                print("MESSAGE: This grid requires corners. They will be added automatically.")
                grid.add_corners()
                break

        """ 
        At the start of every (shifted) grid, fill it using the genetic
        algorithm approach.  If the $gacover block in the input file contains
        the 'nostart' directive, skip this for the initial grid.
        """
        if (nshifts > 0) or not (gaCoverInterface.noStart):
            gaCoverInterface.prepareForRun(grid.get_samples_id(), grid.size[0])
            gaNewSamples = gaCoverInterface.run()            
            if (len(gaNewSamples) > 0):
                for x in gaNewSamples:
                    grid.add_sample(x)

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

            # write current samples to file
            grid.save_samples_to_file("%s/samples.dat" % thisRunOutputs)
            
            for protocol in protocols:
                grid.make_grid_for_protocol (protocol, workdir + "/" + protocol.name, reweightHash)
                
            # remember: properties is a dict where keys are property id's and values are property types (e.g. properties['some-dens'] = 'density')
            # push nearest copies into optimizer
            optimizer.pushNearest(surrogateModelHash)

            for prop in properties:
                referenceValue = optimizer.referenceValues[prop]
                kind = surrogateModelHash[prop]
                grid.compute_final_properties(prop, properties[prop], protocolsHash[prop], protocols, kind)
                grid.save_property_values_to_file (prop, thisRunOutputs + '/' + prop + '_EA_k.dat')
                grid.save_property_err_to_file (prop, thisRunOutputs + '/' + prop + '_dEA_k.dat')
                grid.save_property_diff_to_file (prop, referenceValue, thisRunOutputs + '/' + prop + '_diff.dat')
                #
                if not ('--noplot' in sys.argv):
                    grid.plot_property_to_file (prop, thisRunOutputs + "/" + prop + "_EA_k.pdf")
                    grid.plot_property_err_to_file (prop, thisRunOutputs + "/" + prop + "_dEA_k.pdf")
                    grid.plot_property_diff_to_file (prop, properties[prop], referenceValue, thisRunOutputs + "/" + prop + "_diff.pdf")
            #
            optimizer.fillWithScores (grid)
            optimizer.printToFile (grid, thisRunOutputs + "/optimizer_data.dat")
            if not ('--noplot' in sys.argv):
                optimizer.plotToPdf (grid, thisRunOutputs + "/optimizer_score.pdf")
            nextSample = optimizer.determineNextSample (grid, surrogateModelHash)
            print ("Next sample is %d"  % nextSample)
            if (nextSample == -1):
                break

        # evaluate grid shifting
        new_workdir = base_workdir + "/grid_%d" % (nGrid+1)
        if (nshifts == gridShifter.maxshifts - 1) or not (gridShifter.apply (optimizer, paramLoop, grid, workdir, new_workdir)):
            print ("End of grid shift!")
            break

    # make refined
    subgrid = grid.create_refined_subgrid(subgridHash['factors'], subgridHash['method'], properties)
    subgridOutputs = "%s/subgrid" % workdir
    os.system("mkdir -p " + subgridOutputs)
    subgrid.save_samples_to_file("%s/subgrid/samples.dat" % workdir)

    refined_parameters = paramLoop.loop.create_refined(*subgridHash['factors'])
    refined_parameters.write_parameter_grid_to_file("%s/subgrid/subgrid.dat" % workdir)
    for prop in properties:
        referenceValue = optimizer.referenceValues[prop]
        subgrid.save_property_values_to_file (prop, subgridOutputs + '/' + prop + '_EA_k.dat')
        subgrid.save_property_err_to_file (prop, subgridOutputs + '/' + prop + '_dEA_k.dat')
        subgrid.save_property_diff_to_file (prop, referenceValue, subgridOutputs + '/' + prop + '_diff.dat')
        #
        if not ('--noplot' in sys.argv):
            subgrid.plot_property_to_file (prop, subgridOutputs + "/" + prop + "_EA_k.pdf")
            subgrid.plot_property_err_to_file (prop, subgridOutputs + "/" + prop + "_dEA_k.pdf")
            subgrid.plot_property_diff_to_file (prop, properties[prop], referenceValue, subgridOutputs + "/" + prop + "_diff.pdf")
            #
        optimizer.fillWithScores (subgrid)
        optimizer.printToFile (subgrid, subgridOutputs + "/optimizer_data.dat")
        optimizer.plotToPdf (subgrid, subgridOutputs + "/optimizer_score.pdf")
                
                
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
