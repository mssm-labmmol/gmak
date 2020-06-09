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

bool_legacy = True

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

    for nshifts in range(gridShifter.maxshifts):
        nGrid = nshifts + 1
        nextSample = -1
        for nsteps in range(optimizer.maxSteps):
            workdir = base_workdir + ("/grid_%d/" % nGrid)
            thisRunOutputs = workdir + "/step_%d" % (nsteps+1)

            # Read samples from file.
            grid.read_samples_from_file("%s/samples.dat" % thisRunOutputs)
            
            # remember: properties is a dict where keys are property id's and values are property types (e.g. properties['some-dens'] = 'density')
            # push nearest copies into optimizer
            optimizer.pushNearest(surrogateModelHash)

            for prop in properties:
                referenceValue = optimizer.referenceValues[prop]
                kind = surrogateModelHash[prop]
                grid.read_property_values_err_from_file (prop, properties[prop], thisRunOutputs + '/' + prop + '_EA_k.dat', thisRunOutputs + '/' + prop + '_dEA_k.dat')
                grid.plot_property_to_file (prop, thisRunOutputs + "/" + prop + "_EA_k.pdf")
                grid.plot_property_err_to_file (prop, thisRunOutputs + "/" + prop + "_dEA_k.pdf")
                grid.plot_property_diff_to_file (prop, properties[prop], referenceValue, thisRunOutputs + "/" + prop + "_diff.pdf")
            #
            optimizer.fillWithScores (grid)
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
