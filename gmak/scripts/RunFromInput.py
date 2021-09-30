#!/usr/bin/env python3

import numpy as np
import os
import sys
import copy
import pickle

import gmak.logger as logger
from gmak.gridshifter import *
from gmak.parameters import *
from gmak.grid_ana import *
from gmak.read_input import *
from gmak.mdputils import mdpUtils
from gmak.protocols import LiquidProtocol, GasProtocol, SlabProtocol
from gmak.gridbase import GridPoint, ParameterGrid
from gmak.gridoptimizer import gridOptimizer
from gmak.subgrid import *
from gmak.state import *
from gmak.write_input import *
from gmak.results_assembler import *
from gmak.gridpoint_selector import *

import gmak.runcmd as runcmd
import gmak.simulate as simulate
import gmak.custom_atomic_properties as custom_atomic_properties
import gmak.custom_protocols as custom_protocols
import gmak.custom_surrogate_models as custom_surrogate_models
import gmak.custom_scores as custom_scores
import gmak.custom_topologies as custom_topologies

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

def main():

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
         surrogateModelHash,
         subgridHash) = globalState.getInitializationState()

        logger.globalLogger.putMessage('Restarting from {}'.format(binFilename))

    else:
        globalState.setFromInput(initialize_from_input (sys.argv[1], bool_legacy, validateFlag))

        (base_workdir, # read variables from input file
         grid,
         protocols,
         properties,
         protocolsHash,
         optimizer,
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

    # Initialize logger.globalLogger with output file inside the workdir
    logger.globalLogger = logger.Logger(
        os.path.join(
            base_workdir,
            "gmak_{}.dat".format(os.getpid())))

    # Initialize results assembler - hard-coded for now, will go into input file
    resultsAssembler = ResultsAssembler(grid.getParameterNames())

    # Initialize selector - hard-coded for now, will go into input file
    selector = createSelector('best', gridOptimizer=optimizer, howMany=9)

    logger.globalLogger.putMessage('BEGIN MAINLOOP', dated=True)
    logger.globalLogger.indent()
    grid.run(protocols, optimizer, surrogateModelHash, properties, protocolsHash, resultsAssembler, plotFlag)
    logger.globalLogger.unindent()
    logger.globalLogger.putMessage('END MAINLOOP', dated=True)

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

if __name__ == "__main__":
    main()
