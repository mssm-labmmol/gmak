#!/usr/bin/env python3

from logger import *
from gridshifter import * 
from parameters import * 
from grid_ana import * 
from coverage import *
from read_input import *
from mdputils import mdpUtils
from protocols import LiquidProtocol, GasProtocol, SlabProtocol
from gridbase import GridPoint, ParameterGrid
from gridoptimizer import gridOptimizer
from subgrid import *
from state import *
from write_input import *
from results_assembler import *
from gridpoint_selector import *

import numpy as np
import os
import sys
import copy
import runcmd
import pickle
import simulate
import custom_atomic_properties
import custom_protocols

if ('--legacy' in sys.argv):
    bool_legacy = True
else:
    bool_legacy = False

#plotFlag = not ('--no-plot' in sys.argv)
plotFlag = False # never plot
validateFlag = ('--validate' in sys.argv)
writeNewInput = ('--write-input' in sys.argv)

try:
    cmdNprocs = sys.argv[sys.argv.index('-np') + 1]
    simulate.mdrun_nprocs = int(cmdNprocs)
except ValueError:
    pass

if __name__ == "__main__":

    if ('--restart' in sys.argv): # restart from a given binary file
        binFilename = sys.argv[sys.argv.index('--restart') + 1]
        binFile = open(binFilename, 'rb')
        globalState.setFromBinary(binFile) # load state
        binFile.close()

        (base_workdir, # read variables as if reading from an input file,
         grid,         # but actually reading from state
         protocols,    # 
         properties,   # this way the rest of this program can remain the same
         protocolsHash,
         optimizer,
         gaCoverInterface,
         surrogateModelHash,
         subgridHash) = globalState.getInitializationState()

        globalLogger.putMessage('Restarting from {}'.format(binFilename))

    else:
        globalState.setFromInput(initialize_from_input (sys.argv[1], bool_legacy, validateFlag))

        (base_workdir, # read variables from input file
         grid,
         protocols,
         properties,
         protocolsHash,
         optimizer,
         gaCoverInterface,
         surrogateModelHash,
         subgridHash) = globalState.getInitializationState()

        # In each protocol, check writing frequencies of energy and trajectory to avoid any problems
        mdp_ut = mdpUtils()
        for protocol in protocols:
            # Check only production run.
            mdp_ut.parse_file(protocol.mdps[-1])
            efreq = mdp_ut.get_writefreq_energy()
            xfreq = mdp_ut.get_writefreq_compressedpos()
            if (efreq != xfreq):
                raise ValueError("For protocol {}, efreq ({}) and xfreq ({}) don't match.".format(
                    protocol.name, efreq, xfreq))

        # *********************** End of checks for run  **************************        

    # Initialize results assembler - hard-coded for now, will go into input file
    resultsAssembler = ResultsAssembler(grid.getParameterNames())

    # Initialize selector - hard-coded for now, will go into input file
    selector = createSelector('best', gridOptimizer=optimizer, howMany=9)

    globalLogger.putMessage('BEGIN MAINLOOP', dated=True)
    globalLogger.indent()
    grid.run(protocols, optimizer, surrogateModelHash, properties, protocolsHash, resultsAssembler, plotFlag)
    globalLogger.unindent()
    globalLogger.putMessage('END MAINLOOP', dated=True)

    # DEBUG
    resultsAssembler.print(sys.stdout)

    # Write assembled results to a file
    resultsAssembler.writeToFile(grid.makeCurrentWorkdir() + "/assembled")


    # *********************** Write new input **********************************

    if (writeNewInput):
        mod = InputModifier(sys.argv[1])
        #
        mod.set_workdir(base_workdir + "_best_points")
        mod.set_main_variation(grid.parSpaceGen)
        mod.set_samples(selector.selectPoints())
        mod.write_to_file(sys.argv[1] + "_best_points")

    # *********************** Subgrid part           **************************        

    #subgridHash['parspacegen'] = copy.deepcopy(subgridHash['parspacegen'])
    #subgridHash['parspacegen'].rescale(subgridHash['factors'])
    #subgridObj = ParameterSubgrid(subgridHash['parspacegen'], grid, properties, subgridHash['method'], bool_legacy)
    #for prop in properties:
    #    refValue = optimizer.referenceValues[prop]
    #    subgridObj.save_property_values_to_file(prop)
    #    subgridObj.save_property_diff_to_file(prop, refValue)
    #    subgridObj.save_property_err_to_file(prop)
    #    if plotFlag:
    #        subgridObj.plot_property_to_file(prop)
    #        subgridObj.plot_property_diff_to_file(prop, properties[prop], refValue)
    #        subgridObj.plot_property_err_to_file(prop)
    #subgridObj.printAndPlotScores(optimizer, plotFlag)
    #subgridObj.writeParameters()
    #subgridObj.save_to_binary(optimizer)
